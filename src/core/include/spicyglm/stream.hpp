// The cell method image by image: the neighbour totals of every image without the whole data set in one Dataset,
// and the passes over the cell table's columns that lay the cells out by image and type. Plain C++17 over the
// caller's arrays (R vectors, NumPy arrays), shared by the R and Python packages; the readers (R data frames and
// Rarr, pandas and zarr) stay in each package.
// Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.
#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include "spicyglm/core.hpp"

namespace spicyglm {

// A column of the caller's memory: `size` elements of one of the common numeric types, element i at
// data + i * stride elements.
struct Column {
  enum class Type { F64, F32, I64, I32, I16, I8 };
  const void* data = nullptr;
  Type type = Type::F64;
  std::size_t size = 0;
  std::ptrdiff_t stride = 1;
  // A coordinate stored as I32 equal to the smallest int (R's NA_integer_) reads as this value.
  double int_na = std::numeric_limits<double>::quiet_NaN();

  double real(std::size_t i) const;
  long long code(std::size_t i) const;
};

// Integer codes of images or cell types: a code below `base` is missing (R factor codes: base 1, NA the smallest
// int; pandas / Zarr categorical codes: base 0, -1 missing).
struct Codes {
  Column column;
  int base = 0;
  // the 0-based code of row i, or -1 when missing
  long long at(std::size_t i) const {
    const long long c = column.code(i);
    return c < base ? -1 : c - base;
  }
};

// The cells of a table (or of some of its rows) in the caller's row order.
struct CellColumns {
  Column x, y;
  Codes image, type;
  std::size_t size() const { return image.column.size; }
};

// How codes map to the analysis: image_map[c] is the image (0-based, in the analysis' order) of image code c, or -1
// for a code without cells; type_map[c] the cell type of type code c, with one extra last entry for a missing type.
struct CodeMaps {
  const int* image_map = nullptr;
  std::size_t n_image_codes = 0;
  const int* type_map = nullptr;
  std::size_t n_type_codes = 0;  // type_map has n_type_codes + 1 entries
  int n_images = 0;
};

// What to sum. k > 0: the k nearest neighbours; else the radius graph at `radii` (one radius, or several from one
// neighbour search). allocation: per pair the number of `to` cells with a `from` neighbour and, per image and type,
// its expectation for a self-pair (Dataset::pair_neighbour_any_totals and self_any_expected); else the neighbour
// totals and their squared out-degree sums (pair_neighbour_totals and pair_neighbour_out_sq_totals).
struct StreamOptions {
  std::vector<double> radii;
  int k = 0;
  bool allocation = true;
  int n_types = 0;
  std::size_t n_results() const { return k > 0 ? 1 : radii.size(); }
  std::size_t block_a() const { return static_cast<std::size_t>(n_types) * n_types; }
  std::size_t block_b() const { return allocation ? static_cast<std::size_t>(n_types) : block_a(); }
};

// One pair of arrays per result (radius, or k): a[j] holds block_a() entries per image and b[j] block_b(), image
// after image, in the layout the whole-data Dataset gives (entry (img * n_types + from) * n_types + to). The
// functions below write the blocks of the images they see and leave the rest untouched.
struct StreamOutput {
  std::vector<double*> a, b;
};

// The totals of one batch of whole images already grouped image by image: image i of the batch is rows
// [offsets[i], offsets[i + 1]) of x, y and type (0-based cell types) and goes to output image image_index[i].
// Each image runs in a Dataset of its own, the images handed out in turn to up to n_threads threads, so the
// numbers are exactly those of the whole-data Dataset whatever the threads, and the memory is that of n_threads
// images' neighbour indices.
void batch_totals(const double* x, const double* y, const int* type, const int* offsets, int n_batch_images,
                  const int* image_index, const StreamOptions& opt, const StreamOutput& out, int n_threads = 1);

// The totals of the images in `rows`, cells in any row order (each image's cells in row order count as sorted
// stably by image); every image present must be in the rows whole. The cells are gathered image by image and run
// as batch_totals runs them.
void rows_totals(const CellColumns& rows, const CodeMaps& maps, const StreamOptions& opt, const StreamOutput& out,
                 int n_threads = 1);

// The totals of every image from the whole table in its own row order, without a sorted copy: batches of
// consecutive images (in the analysis' order) of about `budget` cells (one image at least), each gathered with one
// pass over the codes and run as batch_totals runs it.
void stream_totals(const CellColumns& cells, const CodeMaps& maps, const StreamOptions& opt, const StreamOutput& out,
                   int n_threads = 1, long long budget = 1LL << 21);

// One pass over a block of rows of the image and type codes (rows row_offset, row_offset + 1, ...): for every
// (image code, type code) the cells and the first row (-1 when none; a missing type is the extra last column), for
// every image code its last row and the first position within the block (-1 when absent), the number of changes
// of image between consecutive rows (the first row counts when its image differs from prev_image, the image code
// of the row before the block, -1 for none) and the block's last image code. Matrices are image-major
// (entry i * (n_type_codes + 1) + t). A missing or out-of-range image code, or a type code beyond the types, is
// an error.
struct CodeScan {
  std::vector<long long> counts, first, last, pos;
  long long runs = 0, last_image = -1;
};
CodeScan scan_codes(const Codes& image, std::size_t n_image_codes, const Codes& type, std::size_t n_type_codes,
                    long long row_offset = 0, long long prev_image = -1);

// Whether every row holds the same bytes (width bytes per element, contiguous) as its image's first row
// (first_row: per 0-based image code; codes of rows must be valid). Equal bytes are the same value; unequal bytes
// may still be equal values (strings, NaNs), which the caller checks.
bool constant_within(const void* values, std::size_t width, const Codes& image, const long long* first_row);

// The rows (0-based) whose bytes (width bytes per element) differ from those of ref at the row's image (ref
// indexed by 0-based image code, n_ref entries), with the same meaning of equal as constant_within.
std::vector<long long> rows_differing(const void* values, std::size_t width, const Codes& image, const void* ref,
                                      std::size_t n_ref);

// The rows (ascending) whose image code is selected by keep (indexed by 0-based image code, n_keep entries).
std::vector<long long> rows_of_images(const Codes& image, const bool* keep, std::size_t n_keep);

// Codes of keys (pointer-sized identities, such as R's cached strings) in order of first appearance, written to
// codes (n entries): 1 + the level of keys[i], or missing_code for the key `missing`; level_rows[l] is the first
// row of level l. Returns false (and leaves the result unspecified) when there are more than max_levels levels.
bool first_appearance_codes(const void* const* keys, std::size_t n, const void* missing, int missing_code,
                            std::size_t max_levels, int* codes, std::vector<std::size_t>& level_rows);

}  // namespace spicyglm

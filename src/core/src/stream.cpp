// The cell method image by image (spicyglm/stream.hpp). Every Dataset function used here works image by image on
// the image's cells in row order (rows relative to the image), so a Dataset holding one image gives that image
// exactly the numbers the whole-data Dataset gives it (cells sorted stably by image).
#include "spicyglm/stream.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>

namespace spicyglm {

double Column::real(std::size_t i) const {
  const std::ptrdiff_t at = static_cast<std::ptrdiff_t>(i) * stride;
  switch (type) {
    case Type::F64: return static_cast<const double*>(data)[at];
    case Type::F32: return static_cast<double>(static_cast<const float*>(data)[at]);
    case Type::I32: {
      const std::int32_t v = static_cast<const std::int32_t*>(data)[at];
      return v == std::numeric_limits<std::int32_t>::min() ? int_na : static_cast<double>(v);
    }
    case Type::I64: return static_cast<double>(static_cast<const std::int64_t*>(data)[at]);
    case Type::I16: return static_cast<double>(static_cast<const std::int16_t*>(data)[at]);
    case Type::I8: return static_cast<double>(static_cast<const std::int8_t*>(data)[at]);
  }
  return 0;
}

long long Column::code(std::size_t i) const {
  const std::ptrdiff_t at = static_cast<std::ptrdiff_t>(i) * stride;
  switch (type) {
    case Type::I32: return static_cast<const std::int32_t*>(data)[at];
    case Type::I64: return static_cast<const std::int64_t*>(data)[at];
    case Type::I16: return static_cast<const std::int16_t*>(data)[at];
    case Type::I8: return static_cast<const std::int8_t*>(data)[at];
    case Type::F64: case Type::F32: break;
  }
  throw std::invalid_argument("codes must be integers");
}

namespace {

// One image's cells, gathered in row order, and where its totals go.
struct ImageCells {
  std::vector<double> x, y;
  std::vector<int> type;
  int index = 0;
};

void check(const StreamOptions& opt, const StreamOutput& out) {
  if (opt.n_types < 1) throw std::invalid_argument("no cell types");
  if (opt.k <= 0 && opt.radii.empty()) throw std::invalid_argument("give a radius or k");
  if (out.a.size() != opt.n_results() || out.b.size() != opt.n_results() || out.c.size() != opt.n_results())
    throw std::invalid_argument("one set of output arrays (a, b, c) per radius");
}

// The totals of one image into its blocks of the outputs.
void image_totals(ImageCells& c, const StreamOptions& opt, const StreamOutput& out) {
  const int n = static_cast<int>(c.x.size());
  Dataset d(std::move(c.x), std::move(c.y), std::move(c.type), std::vector<int>{0, n}, opt.n_types);
  const std::size_t A = opt.block_a(), B = opt.block_b();
  const std::size_t at = static_cast<std::size_t>(c.index) * A, at_b = static_cast<std::size_t>(c.index) * B;
  const std::size_t at_c = static_cast<std::size_t>(c.index) * opt.block_c();
  auto put = [](const std::vector<double>& from, double* to) { std::copy(from.begin(), from.end(), to); };
  const bool knn = opt.k > 0;
  if (!knn && opt.radii.size() > 1) {
    Dataset::MultiRadiusTotals m = d.multi_radius_totals(opt.radii, opt.allocation);
    for (std::size_t j = 0; j < opt.radii.size(); ++j) {
      if (opt.allocation) { put(m.any[j], out.a[j] + at); put(m.self_expected[j], out.b[j] + at_b); }
      else { put(m.totals[j], out.a[j] + at); put(m.sq[j], out.b[j] + at_b); }
      put(m.self_var[j], out.c[j] + at_c);
    }
    return;
  }
  if (knn) d.build_knn(opt.k, 1); else d.build_radius_index(opt.radii[0]);
  if (opt.allocation) {
    put(d.pair_neighbour_any_totals(knn), out.a[0] + at);
    put(d.self_any_expected(knn), out.b[0] + at_b);
  } else {
    put(d.pair_neighbour_totals(knn), out.a[0] + at);
    put(d.pair_neighbour_out_sq_totals(knn), out.b[0] + at_b);
  }
  put(d.self_variances(knn, opt.allocation), out.c[0] + at_c);
}

// Each image in a Dataset of its own, the images handed out in turn to up to n_threads threads: make(i, cells)
// gives image i's cells when its thread takes it up, and they are released once it is done.
template <typename Make>
void run_images(int n, Make make, const StreamOptions& opt, const StreamOutput& out, int n_threads) {
  const int G = std::max(1, std::min(n_threads, n));
  std::atomic<int> next(0);
  std::vector<std::exception_ptr> err(G);
  auto work = [&](int g) {
    try {
      for (int i = next++; i < n; i = next++) {
        ImageCells c;
        make(i, c);
        image_totals(c, opt, out);
      }
    } catch (...) {
      err[g] = std::current_exception();
      next = n;
    }
  };
  std::vector<std::thread> pool;
  for (int g = 1; g < G; ++g) pool.emplace_back(work, g);
  work(0);
  for (std::thread& t : pool) t.join();
  for (const std::exception_ptr& e : err) if (e) std::rethrow_exception(e);
}

// run_images over images already gathered
void run_gathered(std::vector<ImageCells>& images, const StreamOptions& opt, const StreamOutput& out, int n_threads) {
  run_images(static_cast<int>(images.size()), [&](int i, ImageCells& c) { c = std::move(images[i]); }, opt, out,
             n_threads);
}

int image_of(const CellColumns& c, const CodeMaps& maps, std::size_t r) {
  const long long code = c.image.at(r);
  if (code < 0 || code >= static_cast<long long>(maps.n_image_codes) || maps.image_map[code] < 0 ||
      maps.image_map[code] >= maps.n_images)
    throw std::invalid_argument("missing or invalid image ID in row " + std::to_string(r + 1));
  return maps.image_map[code];
}

int type_of(const CellColumns& c, const CodeMaps& maps, int n_types, std::size_t r) {
  long long code = c.type.at(r);
  if (code < 0) code = static_cast<long long>(maps.n_type_codes);
  if (code > static_cast<long long>(maps.n_type_codes) || maps.type_map[code] < 0 || maps.type_map[code] >= n_types)
    throw std::invalid_argument("invalid cell type code in row " + std::to_string(r + 1));
  return maps.type_map[code];
}

void check_cells(const CellColumns& c) {
  const std::size_t n = c.size();
  if (c.x.size != n || c.y.size != n || c.type.column.size != n)
    throw std::invalid_argument("x, y and the codes differ in length");
}

}  // namespace

void batch_totals(const double* x, const double* y, const int* type, const int* offsets, int n_batch_images,
                  const int* image_index, const StreamOptions& opt, const StreamOutput& out, int n_threads) {
  check(opt, out);
  if (n_batch_images < 1) throw std::invalid_argument("no images in the batch");
  if (offsets[0] != 0) throw std::invalid_argument("batch: offsets must start at 0");
  for (int i = 0; i < n_batch_images; ++i)
    if (offsets[i + 1] < offsets[i]) throw std::invalid_argument("batch: offsets must be increasing");
  for (int r = 0; r < offsets[n_batch_images]; ++r)
    if (type[r] < 0 || type[r] >= opt.n_types) throw std::invalid_argument("batch: cell type out of range");
  // each image's cells copied when its thread takes it up
  run_images(n_batch_images, [&](int i, ImageCells& c) {
    const int s = offsets[i], e = offsets[i + 1];
    c.x.assign(x + s, x + e);
    c.y.assign(y + s, y + e);
    c.type.assign(type + s, type + e);
    c.index = image_index[i];
  }, opt, out, n_threads);
}

void rows_totals(const CellColumns& rows, const CodeMaps& maps, const StreamOptions& opt, const StreamOutput& out,
                 int n_threads) {
  check(opt, out);
  check_cells(rows);
  const std::size_t n = rows.size();
  // the images of the rows in order of first appearance, with their cells in row order
  std::vector<int> local(static_cast<std::size_t>(maps.n_images), -1);
  std::vector<ImageCells> images;
  std::vector<std::size_t> size;
  for (std::size_t r = 0; r < n; ++r) {
    int& l = local[image_of(rows, maps, r)];
    if (l < 0) {
      l = static_cast<int>(images.size());
      images.emplace_back();
      images.back().index = image_of(rows, maps, r);
      size.push_back(0);
    }
    ++size[l];
  }
  for (std::size_t i = 0; i < images.size(); ++i) {
    images[i].x.reserve(size[i]); images[i].y.reserve(size[i]); images[i].type.reserve(size[i]);
  }
  for (std::size_t r = 0; r < n; ++r) {
    ImageCells& c = images[local[image_of(rows, maps, r)]];
    c.x.push_back(rows.x.real(r)); c.y.push_back(rows.y.real(r));
    c.type.push_back(type_of(rows, maps, opt.n_types, r));
  }
  run_gathered(images, opt, out, n_threads);
}

void stream_totals(const CellColumns& cells, const CodeMaps& maps, const StreamOptions& opt, const StreamOutput& out,
                   int n_threads, long long budget) {
  check(opt, out);
  check_cells(cells);
  const std::size_t N = cells.size();
  const int n_images = maps.n_images;
  std::vector<long long> size(n_images, 0);
  for (std::size_t r = 0; r < N; ++r) ++size[image_of(cells, maps, r)];
  for (int lo = 0, hi; lo < n_images; lo = hi) {
    long long n = 0;
    for (hi = lo; hi < n_images && (hi == lo || n + size[hi] <= budget); ++hi) n += size[hi];
    // the batch's cells, image by image in row order (one pass over the codes)
    std::vector<ImageCells> images(hi - lo);
    for (int i = lo; i < hi; ++i) {
      ImageCells& c = images[i - lo];
      c.index = i;
      c.x.reserve(size[i]); c.y.reserve(size[i]); c.type.reserve(size[i]);
    }
    for (std::size_t r = 0; r < N; ++r) {
      const int m = maps.image_map[cells.image.at(r)];
      if (m < lo || m >= hi) continue;
      ImageCells& c = images[m - lo];
      c.x.push_back(cells.x.real(r)); c.y.push_back(cells.y.real(r));
      c.type.push_back(type_of(cells, maps, opt.n_types, r));
    }
    run_gathered(images, opt, out, n_threads);
  }
}

CodeScan scan_codes(const Codes& image, std::size_t n_image_codes, const Codes& type, std::size_t n_type_codes,
                    long long row_offset, long long prev_image) {
  const std::size_t N = image.column.size;
  if (type.column.size != N) throw std::invalid_argument("image and cell type codes differ in length");
  const std::size_t T = n_type_codes + 1;
  CodeScan s;
  s.counts.assign(n_image_codes * T, 0);
  s.first.assign(n_image_codes * T, -1);
  s.last.assign(n_image_codes, -1);
  s.pos.assign(n_image_codes, -1);
  long long last = prev_image;
  for (std::size_t r = 0; r < N; ++r) {
    const long long row = row_offset + static_cast<long long>(r);
    const long long i = image.at(r);
    if (i < 0 || i >= static_cast<long long>(n_image_codes))
      throw std::invalid_argument("missing or invalid image ID in row " + std::to_string(row + 1));
    long long t = type.at(r);
    if (t < 0) t = static_cast<long long>(n_type_codes);
    if (t > static_cast<long long>(n_type_codes))
      throw std::invalid_argument("invalid cell type code in row " + std::to_string(row + 1));
    const std::size_t at = static_cast<std::size_t>(i) * T + static_cast<std::size_t>(t);
    if (s.counts[at]++ == 0) s.first[at] = row;
    s.last[i] = row;
    if (s.pos[i] < 0) s.pos[i] = static_cast<long long>(r);
    if (last != i) ++s.runs;
    last = i;
  }
  s.last_image = last;
  return s;
}

bool constant_within(const void* values, std::size_t width, const Codes& image, const long long* first_row) {
  const char* v = static_cast<const char*>(values);
  for (std::size_t r = 0; r < image.column.size; ++r) {
    const long long i = image.at(r);
    if (i < 0) throw std::invalid_argument("missing image code");
    const std::size_t f = static_cast<std::size_t>(first_row[i]);
    if (std::memcmp(v + r * width, v + f * width, width) != 0) return false;
  }
  return true;
}

std::vector<long long> rows_differing(const void* values, std::size_t width, const Codes& image, const void* ref,
                                      std::size_t n_ref) {
  const char* v = static_cast<const char*>(values);
  const char* f = static_cast<const char*>(ref);
  std::vector<long long> out;
  for (std::size_t r = 0; r < image.column.size; ++r) {
    const long long i = image.at(r);
    if (i < 0 || i >= static_cast<long long>(n_ref)) throw std::invalid_argument("image code out of range");
    if (std::memcmp(v + r * width, f + static_cast<std::size_t>(i) * width, width) != 0)
      out.push_back(static_cast<long long>(r));
  }
  return out;
}

std::vector<long long> rows_of_images(const Codes& image, const bool* keep, std::size_t n_keep) {
  std::vector<long long> rows;
  for (std::size_t r = 0; r < image.column.size; ++r) {
    const long long i = image.at(r);
    if (i >= 0 && i < static_cast<long long>(n_keep) && keep[i]) rows.push_back(static_cast<long long>(r));
  }
  return rows;
}

bool first_appearance_codes(const void* const* keys, std::size_t n, const void* missing, int missing_code,
                            std::size_t max_levels, int* codes, std::vector<std::size_t>& level_rows) {
  level_rows.clear();
  std::unordered_map<const void*, int> code;
  const void* last = nullptr;
  int last_code = 0;
  bool have_last = false;
  for (std::size_t r = 0; r < n; ++r) {
    const void* s = keys[r];
    if (have_last && s == last) { codes[r] = last_code; continue; }
    int c = missing_code;
    if (s != missing) {
      auto it = code.find(s);
      if (it != code.end()) {
        c = it->second;
      } else {
        if (level_rows.size() >= max_levels) return false;
        level_rows.push_back(r);
        c = static_cast<int>(level_rows.size());
        code.emplace(s, c);
      }
    }
    codes[r] = last_code = c;
    last = s;
    have_last = true;
  }
  return true;
}

}  // namespace spicyglm

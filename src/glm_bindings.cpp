// R bindings for the spicyGLM numeric core, mirroring cpp/src/bindings.cpp
// (pybind11). The core itself is untouched and shared between the two.
#include <Rcpp.h>

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "spicyglm/core.hpp"
#include "spicyglm/stream.hpp"

using namespace Rcpp;
using namespace spicyglm;

namespace {

Window parse_window(const std::string& w) {
  if (w == "convex") return Window::Convex;
  if (w == "rectangle") return Window::Rectangle;
  stop("window must be 'convex' or 'rectangle'");
}

std::vector<int> ivec(const IntegerVector& v) { return std::vector<int>(v.begin(), v.end()); }
std::vector<double> dvec(const NumericVector& v) { return std::vector<double>(v.begin(), v.end()); }

template <typename T>
IntegerVector iv(const std::vector<T>& v) { return IntegerVector(v.begin(), v.end()); }
NumericVector nv(const std::vector<double>& v) { return NumericVector(v.begin(), v.end()); }

// PairFit -> list(beta, v_hat, df[, patient, image]), matching pair_fit_dict()
List pair_fit_list(const PairFit& p) {
  List out = List::create(
      _["beta"] = NumericVector::create(p.fit.beta[0], p.fit.beta[1]),
      _["v_hat"] = p.naive ? p.v_naive : p.cr2.v_hat,
      _["df"] = p.naive ? R_PosInf : p.cr2.df);
  if (!p.has_diagnostics) return out;
  const PairDiagnostics& D = p.diagnostics;
  out["patient"] = List::create(
      _["cluster_id"] = iv(p.cr2.cluster_id),
      _["group"] = iv(p.cr2.group),
      _["n_i"] = nv(D.n_i),
      _["T_i"] = nv(D.T),
      _["S_g"] = NumericVector::create(D.S_leverage[0], D.S_leverage[1]),
      _["l_i"] = nv(D.l),
      _["raw_residual_sum"] = nv(D.raw_sum),
      _["adjusted_residual_sum"] = nv(D.adjusted_sum),
      _["e_i"] = nv(p.cr2.e),
      _["influence_i"] = nv(D.influence),
      _["y_i"] = nv(D.y),
      _["d_i"] = nv(D.d),
      _["delta_i"] = nv(D.delta));
  out["image"] = List::create(
      _["cluster"] = iv(D.image_cluster),
      _["image_id"] = iv(D.image_id),
      _["n_ij"] = nv(D.n_ij),
      _["density_ij"] = nv(D.density_ij),
      _["l_ij"] = nv(D.l_ij),
      _["l_ij_group_share"] = nv(D.l_ij_group_share),
      _["raw_residual_sum_ij"] = nv(D.raw_sum_ij),
      _["adjusted_residual_sum_ij"] = nv(D.adjusted_sum_ij),
      _["e_ij"] = nv(D.e_ij),
      _["e_ij_share_within_patient"] = nv(D.e_share_ij),
      _["influence_ij"] = nv(D.influence_ij));
  return out;
}

}  // namespace

// [[Rcpp::export]]
SEXP dataset_create(NumericVector x, NumericVector y, IntegerVector cell_type,
                    IntegerVector image_offsets, int n_types) {
  XPtr<Dataset> p(new Dataset(dvec(x), dvec(y), ivec(cell_type), ivec(image_offsets), n_types), true);
  return p;
}

// [[Rcpp::export]]
NumericVector dataset_image_areas(SEXP ptr, std::string window) {
  XPtr<Dataset> d(ptr);
  std::vector<double> a = d->image_areas(parse_window(window));
  return NumericVector(a.begin(), a.end());
}

// [[Rcpp::export]]
void dataset_build_radius_index(SEXP ptr, double r) {
  XPtr<Dataset> d(ptr);
  d->build_radius_index(r);
}

// [[Rcpp::export]]
void dataset_build_knn(SEXP ptr, int k, int n_threads = 1) {
  XPtr<Dataset> d(ptr);
  d->build_knn(k, n_threads);
}

// Neighbour row indices (0-based, nearest first, -1 for images with at most k cells), N x k row-major.
// [[Rcpp::export]]
IntegerVector knn_rows(NumericVector x, NumericVector y, IntegerVector image_offsets, int k, int n_threads = 1) {
  return iv(knn_indices(dvec(x), dvec(y), ivec(image_offsets), k, n_threads));
}

// [[Rcpp::export]]
List dataset_poisson_model_data(SEXP ptr, NumericVector image_area, int from, int to) {
  XPtr<Dataset> d(ptr);
  ModelData md = d->poisson_model_data(dvec(image_area), from, to);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = IntegerVector(md.n.begin(), md.n.end()),
                      _["density"] = NumericVector(md.density.begin(), md.density.end()));
}

// [[Rcpp::export]]
void dataset_build_intensity(SEXP ptr, double sigma, double min_lambda, std::string window) {
  XPtr<Dataset> d(ptr);
  d->build_intensity(sigma, min_lambda, parse_window(window));
}

// [[Rcpp::export]]
List dataset_inhom_model_data(SEXP ptr, NumericVector image_area, int from, int to, bool edge_correct) {
  XPtr<Dataset> d(ptr);
  InhomModelData md = d->inhom_model_data(dvec(image_area), from, to, edge_correct);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = nv(md.n),
                      _["n_raw"] = IntegerVector(md.n_raw.begin(), md.n_raw.end()),
                      _["weight"] = nv(md.weight),
                      _["density"] = nv(md.density));
}

// [[Rcpp::export]]
void dataset_build_context(SEXP ptr, IntegerVector context_types, std::string window, bool edge_correct) {
  XPtr<Dataset> d(ptr);
  d->build_context(ivec(context_types), parse_window(window), edge_correct);
}

// [[Rcpp::export]]
List dataset_kontextual_model_data(SEXP ptr, int from, int to) {
  XPtr<Dataset> d(ptr);
  InhomModelData md = d->kontextual_model_data(from, to);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = nv(md.n),
                      _["n_raw"] = IntegerVector(md.n_raw.begin(), md.n_raw.end()),
                      _["weight"] = nv(md.weight),
                      _["density"] = nv(md.density));
}

// [[Rcpp::export]]
List dataset_kontextual_binomial_model_data(SEXP ptr, int from, int to) {
  XPtr<Dataset> d(ptr);
  BinomialTrialsModelData md = d->kontextual_binomial_model_data(from, to);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = IntegerVector(md.n.begin(), md.n.end()),
                      _["trials"] = IntegerVector(md.trials.begin(), md.trials.end()),
                      _["p0"] = nv(md.p0));
}

// [[Rcpp::export]]
List dataset_binomial_model_data(SEXP ptr, int from, int to) {
  XPtr<Dataset> d(ptr);
  BinomialModelData md = d->binomial_model_data(from, to);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = IntegerVector(md.n.begin(), md.n.end()),
                      _["p0"] = NumericVector(md.p0.begin(), md.p0.end()));
}

// n is numeric so the inhomogeneous design's pair-weighted counts pass through;
// integer counts convert exactly, leaving the homogeneous fit unchanged.
// [[Rcpp::export]]
List fit_pair_poisson_cpp(IntegerVector cluster, IntegerVector image, IntegerVector group,
                          NumericVector n, NumericVector density, std::string estimator,
                          std::string variance, bool diagnostics = false) {
  return pair_fit_list(fit_pair_poisson(ivec(cluster), ivec(image), ivec(group), dvec(n),
                                        dvec(density), estimator, variance, diagnostics));
}

// [[Rcpp::export]]
List fit_pair_binomial_cpp(IntegerVector cluster, IntegerVector image, IntegerVector group,
                           IntegerVector n, int k, NumericVector p0, std::string estimator,
                           std::string variance) {
  return pair_fit_list(fit_pair_binomial(ivec(cluster), ivec(image), ivec(group), ivec(n), k,
                                         dvec(p0), estimator, variance));
}

// [[Rcpp::export]]
List fit_pair_binomial_trials_cpp(IntegerVector cluster, IntegerVector image, IntegerVector group,
                                  IntegerVector n, IntegerVector trials, NumericVector p0,
                                  std::string estimator, std::string variance) {
  return pair_fit_list(fit_pair_binomial(ivec(cluster), ivec(image), ivec(group), ivec(n),
                                         ivec(trials), dvec(p0), estimator, variance));
}

// Every pair's per-image neighbour totals, as an n_types x n_types x n_images
// array [from, to, image] (see Dataset::pair_neighbour_totals).
// [[Rcpp::export]]
NumericVector dataset_pair_neighbour_totals(SEXP ptr, bool knn, int n_types) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->pair_neighbour_totals(knn);
  NumericVector out(v.begin(), v.end());
  // stored [img][from][to] row-major, i.e. column-major [to, from, img]
  out.attr("dim") = IntegerVector::create(n_types, n_types, static_cast<int>(v.size() / (n_types * n_types)));
  return out;
}

// Sum over each image's `to` cells of the squared number of `from` cells that
// count it as a neighbour, [to, from, image] like dataset_pair_neighbour_totals.
// [[Rcpp::export]]
NumericVector dataset_pair_neighbour_sq_totals(SEXP ptr, bool knn, int n_types) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->pair_neighbour_sq_totals(knn);
  NumericVector out(v.begin(), v.end());
  out.attr("dim") = IntegerVector::create(n_types, n_types, static_cast<int>(v.size() / (n_types * n_types)));
  return out;
}

// [[Rcpp::export]]
NumericVector dataset_pair_neighbour_out_sq_totals(SEXP ptr, bool knn, int n_types) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->pair_neighbour_out_sq_totals(knn);
  NumericVector out(v.begin(), v.end());
  out.attr("dim") = IntegerVector::create(n_types, n_types, static_cast<int>(v.size() / (n_types * n_types)));
  return out;
}

// The allocation design: the number of each image's `to` cells with at least one `from` cell among their
// neighbours, [to, from, image] like dataset_pair_neighbour_totals (see Dataset::pair_neighbour_any_totals).
// [[Rcpp::export]]
NumericVector dataset_pair_neighbour_any_totals(SEXP ptr, bool knn, int n_types) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->pair_neighbour_any_totals(knn);
  NumericVector out(v.begin(), v.end());
  out.attr("dim") = IntegerVector::create(n_types, n_types, static_cast<int>(v.size() / (n_types * n_types)));
  return out;
}

// Several radii from one neighbour search (Dataset::multi_radius_totals): a list with one element per radius, each
// with `any` ([to, from, image]) and `self_expected` (types x images) for the allocation design, or `totals` and `sq`
// ([to, from, image]) for the count, exactly as the single-radius bindings after dataset_build_radius_index.
// [[Rcpp::export]]
List dataset_multi_radius_totals(SEXP ptr, NumericVector radii, bool allocation, int n_types) {
  XPtr<Dataset> d(ptr);
  Dataset::MultiRadiusTotals m = d->multi_radius_totals(std::vector<double>(radii.begin(), radii.end()), allocation);
  const int K = radii.size();
  const int n_img = static_cast<int>((allocation ? m.any[0].size() : m.totals[0].size()) / (static_cast<std::size_t>(n_types) * n_types));
  auto cube = [&](const std::vector<double>& v) { NumericVector o(v.begin(), v.end());
    o.attr("dim") = IntegerVector::create(n_types, n_types, n_img); return o; };
  List out(K);
  for (int k = 0; k < K; ++k) {
    if (allocation) {
      NumericMatrix se(n_types, n_img); std::copy(m.self_expected[k].begin(), m.self_expected[k].end(), se.begin());
      out[k] = List::create(_["any"] = cube(m.any[k]), _["self_expected"] = se);
    } else {
      out[k] = List::create(_["totals"] = cube(m.totals[k]), _["sq"] = cube(m.sq[k]));
    }
  }
  return out;
}

// Exact null expectation of a self-pair's allocation count, types x images (Dataset::self_any_expected).
// [[Rcpp::export]]
NumericMatrix dataset_self_any_expected(SEXP ptr, bool knn, int n_types) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->self_any_expected(knn);
  NumericMatrix out(n_types, static_cast<int>(v.size() / n_types));
  std::copy(v.begin(), v.end(), out.begin());
  return out;
}

// [[Rcpp::export]]
List dataset_rl_model_data(SEXP ptr, int from, int to) {
  XPtr<Dataset> d(ptr);
  ModelData md = d->rl_model_data(from, to);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = IntegerVector(md.n.begin(), md.n.end()),
                      _["density"] = NumericVector(md.density.begin(), md.density.end()));
}

// [[Rcpp::export]]
NumericMatrix dataset_weighted_phi_sums(SEXP ptr, int from, int to, int design, bool edge_correct) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->weighted_phi_sums(from, to, design, edge_correct);
  NumericMatrix out(3, static_cast<int>(v.size() / 3));
  std::copy(v.begin(), v.end(), out.begin());
  return out;   // rows: sum c_b, sum c_b^2, number of candidates; one column per image
}

// [[Rcpp::export]]
NumericMatrix dataset_kontextual_sums(SEXP ptr, int from, int to) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->kontextual_sums(from, to);
  NumericMatrix out(7, static_cast<int>(v.size() / 7));
  std::copy(v.begin(), v.end(), out.begin());
  return out;   // rows: O, L, Q, M, D, unweighted pairs, context cells; one column per image
}

// [[Rcpp::export]]
NumericMatrix dataset_weighted_excess_sums(SEXP ptr, NumericVector a, IntegerVector mPtr, IntegerVector mLab,
                                           NumericVector mW, int n_labels, bool allocation) {
  XPtr<Dataset> d(ptr);
  std::vector<int> lab(mLab.size());
  for (R_xlen_t m = 0; m < mLab.size(); ++m) lab[m] = mLab[m] - 1;   // 1-based labels from R
  std::vector<double> v = d->weighted_excess_sums(std::vector<double>(a.begin(), a.end()),
                                                  std::vector<int>(mPtr.begin(), mPtr.end()), lab,
                                                  std::vector<double>(mW.begin(), mW.end()), n_labels, allocation);
  const int width = 3 + 4 * n_labels;
  NumericMatrix out(width, static_cast<int>(v.size() / width));
  std::copy(v.begin(), v.end(), out.begin());
  return out;   // rows: L, Q, M, then per label O, sum b, sum b^2 (candidates), sum b (image); one column per image
}

// [[Rcpp::export]]
NumericMatrix dataset_kontextual_sums_weighted(SEXP ptr, NumericVector a, NumericVector b, bool self) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->kontextual_sums_weighted(std::vector<double>(a.begin(), a.end()),
                                                      std::vector<double>(b.begin(), b.end()), self);
  NumericMatrix out(10, static_cast<int>(v.size() / 10));
  std::copy(v.begin(), v.end(), out.begin());
  return out;   // rows: O, L, Q, M, D, REF-TARGET weight, context cells, sum b, sum b^2 (candidates), sum b; per image
}

// [[Rcpp::export]]
NumericMatrix dataset_hac_phi_sums_ref(SEXP ptr, int from, int design, double h, int n_types) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->hac_phi_sums_ref(from, design, h);
  NumericMatrix out(n_types + 3, static_cast<int>(v.size() / (n_types + 3)));
  std::copy(v.begin(), v.end(), out.begin());
  return out;   // rows: the HAC variance for each TARGET type, then sum c_b, candidates and G; one column per image
}

// [[Rcpp::export]]
NumericMatrix dataset_hac_phi_sums(SEXP ptr, int from, int to, int design, double h) {
  XPtr<Dataset> d(ptr);
  std::vector<double> v = d->hac_phi_sums(from, to, design, h);
  NumericMatrix out(4, static_cast<int>(v.size() / 4));
  std::copy(v.begin(), v.end(), out.begin());
  return out;   // rows: HAC variance of O, sum c_b, number of candidates, and G with
                // E[HAC] = p (1 - p) G under random labelling; one column per image
}

// ---- the cell method image by image -------------------------------------------------------------------------------
// The loop over images and the passes over the codes live in the core (spicyglm/stream.hpp), shared with spicyr;
// these wrappers only describe R's vectors to it and convert the results.

namespace {

// An R numeric vector (double, or integer read as as.numeric() reads it) as a core column, without a copy.
spicyglm::Column r_column(SEXP v, R_xlen_t n, const char* what) {
  spicyglm::Column c;
  if (XLENGTH(v) != n) stop("%s: wrong length", what);
  if (TYPEOF(v) == REALSXP) { c.data = REAL(v); c.type = spicyglm::Column::Type::F64; }
  else if (TYPEOF(v) == INTSXP || TYPEOF(v) == LGLSXP) { c.data = INTEGER(v); c.type = spicyglm::Column::Type::I32; c.int_na = NA_REAL; }
  else stop("%s must be numeric", what);
  c.size = static_cast<std::size_t>(n);
  return c;
}

// Integer codes (a factor's, base 1; or 0-based Zarr codes, base 0) as core codes.
spicyglm::Codes r_codes(SEXP v, int base, const char* what) {
  if (TYPEOF(v) != INTSXP) stop("%s: expected integer codes", what);
  spicyglm::Codes c;
  c.column = r_column(v, XLENGTH(v), what);
  c.base = base;
  return c;
}

// The elements of an atomic vector as bytes: width and data (strings as their cached CHARSXP pointers).
std::size_t r_width(SEXP v) {
  switch (TYPEOF(v)) {
    case INTSXP: case LGLSXP: return sizeof(int);
    case REALSXP: return sizeof(double);
    case STRSXP: return sizeof(SEXP);
    default: return 0;
  }
}
const void* r_bytes(SEXP v) {
  switch (TYPEOF(v)) {
    case INTSXP: return INTEGER(v);
    case LGLSXP: return LOGICAL(v);
    case REALSXP: return REAL(v);
    default: return STRING_PTR_RO(v);
  }
}

spicyglm::StreamOptions stream_options(int n_types, const NumericVector& radii, int k, bool allocation) {
  spicyglm::StreamOptions o;
  o.radii.assign(radii.begin(), radii.end());
  o.k = k;
  o.allocation = allocation;
  o.n_types = n_types;
  return o;
}

// 0-based rows (-1 none) to 1-based (NA), as an integer vector
IntegerVector r_rows(const std::vector<long long>& v) {
  IntegerVector out(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) out[i] = v[i] < 0 ? NA_INTEGER : static_cast<int>(v[i] + 1);
  return out;
}

}  // namespace

// One pass over a block of rows of the image and type codes (base 1: factor codes, NA missing; base 0: Zarr codes,
// -1 missing), rows offset + 1, offset + 2, ...: for every (image level, type level) the cells and the first row
// (NA when none; a missing type is the extra last column), the number of changes of image between consecutive rows
// (counting the first row when its image differs from `prev`, the 1-based image code of the row before the block,
// NA for none), the last row's image code, and per image level the first position in the block (1-based; NA when
// absent). Codes are returned 1-based. With offset 0 and prev NA, the layout of the whole table.
// [[Rcpp::export]]
List cell_block_scan(SEXP image, int image_base, int n_image_levels, SEXP type, int type_base, int n_type_levels,
                     double offset = 0, int prev = NA_INTEGER) {
  const spicyglm::CodeScan s = spicyglm::scan_codes(r_codes(image, image_base, "image"), n_image_levels,
                                                    r_codes(type, type_base, "cell type"), n_type_levels,
                                                    static_cast<long long>(offset), prev == NA_INTEGER ? -1 : prev - 1);
  const int T = n_type_levels + 1;
  IntegerMatrix counts(n_image_levels, T), first(n_image_levels, T);
  for (int i = 0; i < n_image_levels; ++i)
    for (int t = 0; t < T; ++t) {
      const std::size_t at = static_cast<std::size_t>(i) * T + t;
      counts(i, t) = static_cast<int>(s.counts[at]);
      first(i, t) = s.first[at] < 0 ? NA_INTEGER : static_cast<int>(s.first[at] + 1);
    }
  return List::create(_["counts"] = counts, _["first"] = first, _["runs"] = static_cast<double>(s.runs),
                      _["last"] = s.last_image < 0 ? NA_INTEGER : static_cast<int>(s.last_image + 1),
                      _["pos"] = r_rows(s.pos));
}

// Whether v takes the same value (the same bits; for strings the same cached string) in every row of each image
// as in the image's first row (first_row: per image level, 1-based; factor codes). TRUE implies the values are
// constant within images as as.character() sees them; FALSE leaves the caller to check as.character() values.
// [[Rcpp::export]]
bool cell_constant_within(SEXP v, SEXP image_fac, IntegerVector first_row) {
  if (XLENGTH(v) != XLENGTH(image_fac) || r_width(v) == 0) return false;
  std::vector<long long> f(first_row.size());
  for (R_xlen_t i = 0; i < first_row.size(); ++i) f[i] = first_row[i] == NA_INTEGER ? 0 : first_row[i] - 1;
  return spicyglm::constant_within(r_bytes(v), r_width(v), r_codes(image_fac, 1, "image"), f.data());
}

// The rows (1-based) where v differs from ref at the row's image (ref indexed by the image code): the same bits
// (for strings the same cached string) count as equal, so a row not returned holds the same value as as.character()
// sees it; the caller compares the returned rows as.character(). Image codes are 1-based, or 0-based with
// image_base = 0.
// [[Rcpp::export]]
IntegerVector cell_rows_differing(SEXP v, SEXP image, SEXP ref, int image_base = 1) {
  if (XLENGTH(image) != XLENGTH(v) || TYPEOF(ref) != TYPEOF(v)) stop("cell_rows_differing: inconsistent input");
  if (r_width(v) == 0) stop("cell_rows_differing: unsupported type");
  return r_rows(spicyglm::rows_differing(r_bytes(v), r_width(v), r_codes(image, image_base, "image"), r_bytes(ref),
                                         XLENGTH(ref)));
}

// The neighbour totals of the cell method from the cells in input order (spicyglm::stream_totals): image_fac and
// type_fac are factor codes, mapped by image_map / type_map (indexed by code - 1; type_map's last entry is a missing
// type) to the 0-based image (in .cell_context()'s sorted order) and cell type. k > 0: the k nearest neighbours; else
// the radius graph at radii. One element per radius (one for kNN), each list(any, self_expected) (allocation) or
// list(totals, sq) (count), as the whole-data dataset gives them. Batches of about 2^21 cells per thread.
// [[Rcpp::export]]
List cell_stream_totals(SEXP x, SEXP y, SEXP image_fac, IntegerVector image_map, int n_images, SEXP type_fac,
                        IntegerVector type_map, int n_types, NumericVector radii, int k, bool allocation,
                        int n_threads = 1) {
  const R_xlen_t N = XLENGTH(image_fac);
  spicyglm::CellColumns cells{r_column(x, N, "x"), r_column(y, N, "y"), r_codes(image_fac, 1, "image"),
                              r_codes(type_fac, 1, "cell type")};
  if (XLENGTH(type_fac) != N) stop("cell type: wrong length");
  spicyglm::CodeMaps maps{image_map.begin(), static_cast<std::size_t>(image_map.size()), type_map.begin(),
                          static_cast<std::size_t>(type_map.size() - 1), n_images};
  const spicyglm::StreamOptions o = stream_options(n_types, radii, k, allocation);
  List res(o.n_results());
  spicyglm::StreamOutput out;
  for (std::size_t j = 0; j < o.n_results(); ++j) {
    NumericVector a(o.block_a() * n_images), b(o.block_b() * n_images);
    a.attr("dim") = IntegerVector::create(n_types, n_types, n_images);
    if (allocation) b.attr("dim") = IntegerVector::create(n_types, n_images);
    else b.attr("dim") = IntegerVector::create(n_types, n_types, n_images);
    out.a.push_back(a.begin()); out.b.push_back(b.begin());
    res[j] = allocation ? List::create(_["any"] = a, _["self_expected"] = b) : List::create(_["totals"] = a, _["sq"] = b);
  }
  spicyglm::stream_totals(cells, maps, o, out, n_threads, (1LL << 21) * std::max(1, n_threads));
  return res;
}

// The neighbour totals of one batch of whole images that the caller has already gathered (spicyglm::batch_totals):
// the cells image by image (each image's cells in row order), x, y (or x the N x 2 matrix of both and y empty) and
// 0-based cell types, offsets[i] the first cell of image i of the batch. One element per radius (one for kNN), each
// list(a, b): the any / totals block and the self_expected / sq block of every image of the batch.
// [[Rcpp::export]]
List cell_batch_totals(NumericVector x, NumericVector y, IntegerVector type, IntegerVector offsets, int n_types,
                       NumericVector radii, int k, bool allocation, int n_threads = 1) {
  const int nb = offsets.size() - 1;
  const R_xlen_t N = type.size();
  const bool matrix = y.size() == 0 && x.size() == 2 * N;   // y empty: x is the N x 2 matrix of the coordinates
  if (nb < 1 || (!matrix && (x.size() != N || y.size() != N)) || offsets[nb] != N) stop("batch: inconsistent lengths");
  const spicyglm::StreamOptions o = stream_options(n_types, radii, k, allocation);
  std::vector<int> index(nb);
  for (int i = 0; i < nb; ++i) index[i] = i;
  List res(o.n_results());
  spicyglm::StreamOutput out;
  for (std::size_t j = 0; j < o.n_results(); ++j) {
    NumericVector a(o.block_a() * nb), b(o.block_b() * nb);
    out.a.push_back(a.begin()); out.b.push_back(b.begin());
    res[j] = List::create(_["a"] = a, _["b"] = b);
  }
  spicyglm::batch_totals(x.begin(), matrix ? x.begin() + N : y.begin(), type.begin(), offsets.begin(), nb,
                         index.data(), o, out, n_threads);
  return res;
}

// factor(x, levels = unique(x)) for a character vector without attributes, in one pass with no hash of the whole
// vector (R's unique() hashes every element: about 3 GB transient at 100 million cells). Equal ASCII strings are
// the same cached CHARSXP, so strings are compared by pointer; NULL (the caller then uses factor()) when x is not
// such a vector, has a non-ASCII string (equal strings may then differ in their encoding mark), or has more than
// 2^20 distinct values.
// [[Rcpp::export]]
SEXP first_appearance_factor(SEXP x) {
  if (TYPEOF(x) != STRSXP) return R_NilValue;   // the caller passes only vectors without attributes
  const R_xlen_t N = XLENGTH(x);
  IntegerVector out(N);
  std::vector<std::size_t> rows;
  if (!spicyglm::first_appearance_codes(reinterpret_cast<const void* const*>(STRING_PTR_RO(x)), N, NA_STRING,
                                        NA_INTEGER, 1u << 20, out.begin(), rows))
    return R_NilValue;
  CharacterVector lev(rows.size());
  for (std::size_t l = 0; l < rows.size(); ++l) {
    SEXP s = STRING_ELT(x, rows[l]);
    for (const char* p = CHAR(s); *p; ++p) if (static_cast<unsigned char>(*p) >= 0x80) return R_NilValue;
    SET_STRING_ELT(lev, l, s);
  }
  out.attr("levels") = lev;
  out.attr("class") = "factor";
  return out;
}

// R bindings for the spicyGLM numeric core, mirroring cpp/src/bindings.cpp
// (pybind11). The core itself is untouched and shared between the two.
#include <Rcpp.h>

#include <algorithm>
#include <cstring>
#include <limits>
#include <exception>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "spicyglm/core.hpp"

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
// spicy()'s cell method needs only per-image summaries of the neighbour graph, so it never holds the whole data set
// in the core: the cells stay in R's columns, in input order, and are gathered a batch of whole images at a time
// into a Dataset of their own. Every Dataset function used here works image by image on the image's cells in row
// order (rows relative to the image), so a batch's Dataset gives each of its images exactly the numbers the
// whole-data Dataset gives it (cells sorted stably by image, as .cell_context() did).

namespace {

// A factor's codes (1-based, NA allowed) of length n, checked: an R integer vector.
const int* factor_codes(SEXP f, R_xlen_t n, const char* what) {
  if (TYPEOF(f) != INTSXP || XLENGTH(f) != n) stop("%s: expected integer codes of length %d", what, static_cast<int>(n));
  return INTEGER(f);
}

// Coordinates as R stores them (double or integer); an integer is converted as as.numeric() does.
struct Coord {
  const double* d = nullptr;
  const int* i = nullptr;
  Coord(SEXP v, R_xlen_t n, const char* what) {
    if (XLENGTH(v) != n) stop("%s: wrong length", what);
    if (TYPEOF(v) == REALSXP) d = REAL(v);
    else if (TYPEOF(v) == INTSXP) i = INTEGER(v);
    else stop("%s must be numeric", what);
  }
  double operator[](R_xlen_t r) const { return d ? d[r] : (i[r] == NA_INTEGER ? NA_REAL : static_cast<double>(i[r])); }
};

}  // namespace

// One pass over the cells: for every (image level, cell-type level) the number of cells and the first row (1-based;
// NA when there is none). Type NA is the extra last column. Images must not be NA.
// [[Rcpp::export]]
List cell_layout(SEXP image_fac, int n_image_levels, SEXP type_fac, int n_type_levels) {
  const R_xlen_t N = XLENGTH(image_fac);
  const int* im = factor_codes(image_fac, N, "image");
  const int* ty = factor_codes(type_fac, N, "cell type");
  const int T = n_type_levels + 1;
  IntegerMatrix counts(n_image_levels, T), first(n_image_levels, T);
  std::fill(first.begin(), first.end(), NA_INTEGER);
  for (R_xlen_t r = 0; r < N; ++r) {
    const int i = im[r], t = ty[r] == NA_INTEGER ? n_type_levels : ty[r] - 1;
    if (i == NA_INTEGER || i < 1 || i > n_image_levels) stop("missing or invalid image ID in row %d", static_cast<int>(r + 1));
    if (t < 0 || t > n_type_levels) stop("invalid cell type code in row %d", static_cast<int>(r + 1));
    const R_xlen_t at = static_cast<R_xlen_t>(t) * n_image_levels + (i - 1);
    if (counts[at]++ == 0) first[at] = static_cast<int>(r + 1);
  }
  return List::create(_["counts"] = counts, _["first"] = first);
}

// Whether v takes the same value (the same bits; for strings the same cached string) in every row of each image
// as in the image's first row (first_row: per image level, 1-based; codes 1-based). TRUE implies the values are
// constant within images as as.character() sees them; FALSE leaves the caller to check as.character() values.
// [[Rcpp::export]]
bool cell_constant_within(SEXP v, SEXP image_fac, IntegerVector first_row) {
  const R_xlen_t N = XLENGTH(image_fac);
  const int* im = factor_codes(image_fac, N, "image");
  if (XLENGTH(v) != N) return false;
  for (R_xlen_t r = 0; r < N; ++r) {
    const R_xlen_t f = first_row[im[r] - 1] - 1;
    switch (TYPEOF(v)) {
      case INTSXP: case LGLSXP: if (INTEGER(v)[r] != INTEGER(v)[f]) return false; break;
      case REALSXP: if (std::memcmp(REAL(v) + r, REAL(v) + f, sizeof(double)) != 0) return false; break;
      case STRSXP: if (STRING_ELT(v, r) != STRING_ELT(v, f)) return false; break;
      default: return false;
    }
  }
  return true;
}

// The neighbour totals of the cell method, image by image. Cells stay in input order: image_fac and type_fac are
// factor codes (1-based; a type NA is code n_type_levels + 1), mapped by image_map / type_map (indexed by code - 1)
// to the 0-based image (in .cell_context()'s sorted order) and cell type. k > 0: the k nearest neighbours; else the
// radius graph at radii (one radius, or several from one neighbour search). Returns one element per radius (one for
// kNN), each list(any, self_expected) (allocation) or list(totals, sq) (count), exactly as dataset_pair_neighbour_any_totals
// and dataset_self_any_expected, or dataset_pair_neighbour_totals and dataset_pair_neighbour_out_sq_totals, after
// dataset_build_radius_index / dataset_build_knn on the whole data set (or dataset_multi_radius_totals).
// Images are gathered in batches of about 2^21 cells per thread (one image at least), one pass over the codes each.
// [[Rcpp::export]]
List cell_stream_totals(SEXP x, SEXP y, SEXP image_fac, IntegerVector image_map, int n_images, SEXP type_fac,
                        IntegerVector type_map, int n_types, NumericVector radii, int k, bool allocation,
                        int n_threads = 1) {
  const R_xlen_t N = XLENGTH(image_fac);
  const int* im = factor_codes(image_fac, N, "image");
  const int* ty = factor_codes(type_fac, N, "cell type");
  const Coord X(x, N, "x"), Y(y, N, "y");
  const int n_type_levels = type_map.size() - 1;
  const bool knn = k > 0;
  const int K = knn ? 1 : radii.size();
  if (!knn && K < 1) stop("give a radius or k");
  auto image_of = [&](R_xlen_t r) {
    const int i = im[r];
    if (i == NA_INTEGER || i < 1 || i > image_map.size()) stop("missing or invalid image ID in row %d", static_cast<int>(r + 1));
    return image_map[i - 1];
  };
  std::vector<long long> size(n_images, 0);
  for (R_xlen_t r = 0; r < N; ++r) {
    const int m = image_of(r);
    if (m < 0 || m >= n_images) stop("image map out of range");
    ++size[m];
  }
  const std::size_t T = static_cast<std::size_t>(n_types), per = static_cast<std::size_t>(n_images) * T * T;
  std::vector<std::vector<double>> a(K, std::vector<double>(per, 0.0));   // any / totals
  std::vector<std::vector<double>> b(K, std::vector<double>(allocation ? static_cast<std::size_t>(n_images) * T : per, 0.0));
  const std::vector<double> rad(radii.begin(), radii.end());
  const long long budget = (1LL << 21) * std::max(1, n_threads);
  for (int lo = 0, hi; lo < n_images; lo = hi) {
    long long n = 0;
    for (hi = lo; hi < n_images && (hi == lo || n + size[hi] <= budget); ++hi) n += size[hi];
    // the batch's cells, image by image in row order
    std::vector<int> offsets(hi - lo + 1, 0);
    for (int i = lo; i < hi; ++i) offsets[i - lo + 1] = offsets[i - lo] + static_cast<int>(size[i]);
    std::vector<double> bx(n), by(n);
    std::vector<int> bt(n), fill(offsets.begin(), offsets.end() - 1);
    for (R_xlen_t r = 0; r < N; ++r) {
      const int m = image_map[im[r] - 1];
      if (m < lo || m >= hi) continue;
      const int t = ty[r] == NA_INTEGER ? n_type_levels : ty[r] - 1;
      if (t < 0 || t > n_type_levels) stop("invalid cell type code in row %d", static_cast<int>(r + 1));
      const int p = fill[m - lo]++;
      bx[p] = X[r]; by[p] = Y[r]; bt[p] = type_map[t];
    }
    Dataset d(std::move(bx), std::move(by), std::move(bt), std::move(offsets), n_types);
    const std::size_t at = static_cast<std::size_t>(lo) * T * T, at_se = static_cast<std::size_t>(lo) * T;
    auto put = [](const std::vector<double>& from, std::vector<double>& to, std::size_t at) {
      std::copy(from.begin(), from.end(), to.begin() + at); };
    if (K > 1) {
      Dataset::MultiRadiusTotals m = d.multi_radius_totals(rad, allocation);
      for (int j = 0; j < K; ++j) {
        if (allocation) { put(m.any[j], a[j], at); put(m.self_expected[j], b[j], at_se); }
        else { put(m.totals[j], a[j], at); put(m.sq[j], b[j], at); }
      }
      continue;
    }
    if (knn) d.build_knn(k, n_threads); else d.build_radius_index(rad[0]);
    if (allocation) { put(d.pair_neighbour_any_totals(knn), a[0], at); put(d.self_any_expected(knn), b[0], at_se); }
    else { put(d.pair_neighbour_totals(knn), a[0], at); put(d.pair_neighbour_out_sq_totals(knn), b[0], at); }
  }
  auto cube = [&](const std::vector<double>& v) { NumericVector o(v.begin(), v.end());
    o.attr("dim") = IntegerVector::create(n_types, n_types, n_images); return o; };
  List out(K);
  for (int j = 0; j < K; ++j) {
    if (allocation) {
      NumericMatrix se(n_types, n_images); std::copy(b[j].begin(), b[j].end(), se.begin());
      out[j] = List::create(_["any"] = cube(a[j]), _["self_expected"] = se);
    } else {
      out[j] = List::create(_["totals"] = cube(a[j]), _["sq"] = cube(b[j]));
    }
  }
  return out;
}

// One block of rows of a cell table read from a store: image and type codes as integers (base 0: Zarr categorical
// codes, -1 missing; base 1: factor codes, NA missing), rows offset + 1, offset + 2, ... For every (image level,
// type level) the cells and the first row (NA when none; a missing type is the extra last column), the number of
// changes of image between consecutive rows (counting the first row when its image differs from `prev`, the image
// code of the row before the block, NA for none), the last row's image code, and per image level the first position
// in the block (1-based; NA when absent). Codes are returned 1-based.
// [[Rcpp::export]]
List cell_block_scan(IntegerVector image, int image_base, int n_image_levels, IntegerVector type, int type_base,
                     int n_type_levels, int offset, int prev) {
  const R_xlen_t N = image.size();
  if (type.size() != N) stop("block: image and type codes differ in length");
  const int T = n_type_levels + 1;
  IntegerMatrix counts(n_image_levels, T), first(n_image_levels, T);
  std::fill(first.begin(), first.end(), NA_INTEGER);
  IntegerVector pos(n_image_levels, NA_INTEGER);
  double runs = 0;
  int last = prev;
  for (R_xlen_t r = 0; r < N; ++r) {
    const int iv = image[r], tv = type[r];
    const int i = (iv == NA_INTEGER || iv < image_base) ? -1 : iv - image_base;   // 0-based
    if (i < 0 || i >= n_image_levels) stop("missing or invalid image ID in row %d", static_cast<int>(offset + r + 1));
    const int t = (tv == NA_INTEGER || tv < type_base) ? n_type_levels : tv - type_base;
    if (t > n_type_levels) stop("invalid cell type code in row %d", static_cast<int>(offset + r + 1));
    const R_xlen_t at = static_cast<R_xlen_t>(t) * n_image_levels + i;
    if (counts[at]++ == 0) first[at] = static_cast<int>(offset + r + 1);
    if (pos[i] == NA_INTEGER) pos[i] = static_cast<int>(r + 1);
    if (last == NA_INTEGER || last != i + 1) ++runs;
    last = i + 1;
  }
  return List::create(_["counts"] = counts, _["first"] = first, _["runs"] = runs, _["last"] = last, _["pos"] = pos);
}

// The rows (1-based) where v differs from ref at the row's image (ref indexed by the image code, 1-based): the same
// bits (for strings the same cached string) count as equal, so a row not returned holds the same value as
// as.character() sees it; the caller compares the returned rows as.character(). Image codes are 1-based, or 0-based
// with image_base = 0.
// [[Rcpp::export]]
IntegerVector cell_rows_differing(SEXP v, IntegerVector image, SEXP ref, int image_base = 1) {
  const R_xlen_t N = XLENGTH(v);
  if (image.size() != N || TYPEOF(ref) != TYPEOF(v)) stop("cell_rows_differing: inconsistent input");
  const R_xlen_t M = XLENGTH(ref);
  std::vector<int> out;
  for (R_xlen_t r = 0; r < N; ++r) {
    const int i = image[r] == NA_INTEGER ? NA_INTEGER : image[r] - image_base + 1;
    if (i == NA_INTEGER || i < 1 || i > M) stop("cell_rows_differing: image code out of range");
    bool same;
    switch (TYPEOF(v)) {
      case INTSXP: case LGLSXP: same = INTEGER(v)[r] == INTEGER(ref)[i - 1]; break;
      case REALSXP: same = std::memcmp(REAL(v) + r, REAL(ref) + (i - 1), sizeof(double)) == 0; break;
      case STRSXP: same = STRING_ELT(v, r) == STRING_ELT(ref, i - 1); break;
      default: stop("cell_rows_differing: unsupported type");
    }
    if (!same) out.push_back(static_cast<int>(r + 1));
  }
  return IntegerVector(out.begin(), out.end());
}

// The neighbour totals of the cell method for one batch of whole images that the caller has already gathered: the
// cells image by image (each image's cells in row order), x, y (or x the N x 2 matrix of both and y empty) and
// 0-based cell types, with offsets[i] the first cell
// of image i of the batch. Returns one element per radius (one for kNN), each list(a, b): a the n_types x n_types
// block of every image of the batch (any / totals), b its self_expected (n_types per image) or sq block, in the
// layout of cell_stream_totals' output for these images. The images are shared among up to n_threads threads, each
// running a Dataset of its own whole images, so each image's numbers are those of the whole-data Dataset.
// [[Rcpp::export]]
List cell_batch_totals(NumericVector x, NumericVector y, IntegerVector type, IntegerVector offsets, int n_types,
                       NumericVector radii, int k, bool allocation, int n_threads = 1) {
  const int nb = offsets.size() - 1;
  if (nb < 1) stop("no images in the batch");
  const R_xlen_t N = type.size();
  // y empty: x is the N x 2 matrix of the coordinates
  const bool matrix = y.size() == 0 && x.size() == 2 * N;
  if ((!matrix && (x.size() != N || y.size() != N)) || offsets[0] != 0 || offsets[nb] != N)
    stop("batch: inconsistent lengths");
  for (int i = 0; i < nb; ++i) if (offsets[i + 1] < offsets[i]) stop("batch: offsets must be increasing");
  for (R_xlen_t r = 0; r < N; ++r) if (type[r] < 0 || type[r] >= n_types) stop("batch: cell type out of range");
  const bool knn = k > 0;
  const int K = knn ? 1 : radii.size();
  if (!knn && K < 1) stop("give a radius or k");
  const std::size_t T = static_cast<std::size_t>(n_types);
  const std::size_t per = static_cast<std::size_t>(nb) * T * T, per_b = allocation ? static_cast<std::size_t>(nb) * T : per;
  std::vector<std::vector<double>> a(K, std::vector<double>(per, 0.0)), b(K, std::vector<double>(per_b, 0.0));
  const std::vector<double> rad(radii.begin(), radii.end());
  const double* px = x.begin(); const double* py = matrix ? x.begin() + N : y.begin(); const int* pt = type.begin();
  const std::vector<int> off(offsets.begin(), offsets.end());
  // groups of consecutive images with about equal numbers of cells, one per thread
  const int G = std::max(1, std::min(n_threads, nb));
  std::vector<int> cut(1, 0);
  for (int g = 1; g < G; ++g) {
    const long long want = static_cast<long long>(N) * g / G;
    int i = cut.back() + 1;
    while (i < nb && off[i] < want) ++i;
    if (i < nb && i > cut.back()) cut.push_back(i);
  }
  cut.push_back(nb);
  std::vector<std::exception_ptr> err(cut.size() - 1);
  auto work = [&](std::size_t g) {
    try {
      const int lo = cut[g], hi = cut[g + 1];
      const std::size_t s = off[lo], e = off[hi];
      std::vector<int> o(hi - lo + 1);
      for (int i = lo; i <= hi; ++i) o[i - lo] = off[i] - off[lo];
      Dataset d(std::vector<double>(px + s, px + e), std::vector<double>(py + s, py + e), std::vector<int>(pt + s, pt + e),
                std::move(o), n_types);
      const std::size_t at = static_cast<std::size_t>(lo) * T * T, at_se = static_cast<std::size_t>(lo) * T;
      auto put = [](const std::vector<double>& from, std::vector<double>& to, std::size_t at) {
        std::copy(from.begin(), from.end(), to.begin() + at); };
      if (K > 1) {
        Dataset::MultiRadiusTotals m = d.multi_radius_totals(rad, allocation);
        for (int j = 0; j < K; ++j) {
          if (allocation) { put(m.any[j], a[j], at); put(m.self_expected[j], b[j], at_se); }
          else { put(m.totals[j], a[j], at); put(m.sq[j], b[j], at); }
        }
        return;
      }
      if (knn) d.build_knn(k, 1); else d.build_radius_index(rad[0]);
      if (allocation) { put(d.pair_neighbour_any_totals(knn), a[0], at); put(d.self_any_expected(knn), b[0], at_se); }
      else { put(d.pair_neighbour_totals(knn), a[0], at); put(d.pair_neighbour_out_sq_totals(knn), b[0], at); }
    } catch (...) { err[g] = std::current_exception(); }
  };
  if (err.size() == 1) work(0);
  else {
    std::vector<std::thread> pool;
    for (std::size_t g = 0; g < err.size(); ++g) pool.emplace_back(work, g);
    for (std::thread& t : pool) t.join();
  }
  for (const std::exception_ptr& e : err) if (e) std::rethrow_exception(e);
  List out(K);
  for (int j = 0; j < K; ++j)
    out[j] = List::create(_["a"] = NumericVector(a[j].begin(), a[j].end()), _["b"] = NumericVector(b[j].begin(), b[j].end()));
  return out;
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
  std::unordered_map<SEXP, int> code;
  std::vector<SEXP> levels;
  IntegerVector out(N);
  SEXP last = nullptr; int last_code = 0;
  for (R_xlen_t r = 0; r < N; ++r) {
    SEXP s = STRING_ELT(x, r);
    if (s == last) { out[r] = last_code; continue; }
    int c;
    if (s == NA_STRING) c = NA_INTEGER;
    else {
      auto it = code.find(s);
      if (it != code.end()) c = it->second;
      else {
        for (const char* p = CHAR(s); *p; ++p) if (static_cast<unsigned char>(*p) >= 0x80) return R_NilValue;
        if (levels.size() >= (1u << 20)) return R_NilValue;
        levels.push_back(s);
        c = static_cast<int>(levels.size());
        code.emplace(s, c);
      }
    }
    out[r] = last_code = c; last = s;
  }
  CharacterVector lev(levels.size());
  for (std::size_t i = 0; i < levels.size(); ++i) SET_STRING_ELT(lev, i, levels[i]);
  out.attr("levels") = lev;
  out.attr("class") = "factor";
  return out;
}

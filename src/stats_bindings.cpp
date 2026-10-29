// R bindings for the statistics module (core/include/spicyglm/stats.hpp). Thin: convert, call, return.
#include <Rcpp.h>

#include <string>
#include <vector>

#include "spicyglm/core.hpp"
#include "spicyglm/stats.hpp"

using namespace Rcpp;
using namespace spicyglm;

namespace {

std::vector<int> ivec(const IntegerVector& v) { return std::vector<int>(v.begin(), v.end()); }
std::vector<double> dvec(const NumericVector& v) { return std::vector<double>(v.begin(), v.end()); }
NumericVector nv(const std::vector<double>& v) { return NumericVector(v.begin(), v.end()); }

ImageRows rows_from(const List& d) {
  ImageRows r;
  IntegerVector img = d["img"];
  NumericVector O = d["O"], E = d["E"], n = d["n"], v = d["v"];
  r.image = ivec(img); r.O = dvec(O); r.E = dvec(E); r.n = dvec(n); r.v = dvec(v);
  return r;
}

// The image rows as the plain list of their columns.
List rows_list(const ImageRows& r) {
  return List::create(_["img"] = IntegerVector(r.image.begin(), r.image.end()), _["O"] = nv(r.O),
                      _["E"] = nv(r.E), _["n"] = nv(r.n), _["v"] = nv(r.v));
}

// ... and as a data.frame, made directly (the class and compact row names as.data.frame gives) rather than by
// DataFrame::create, which calls R's as.data.frame
List rows_to(const ImageRows& r) {
  List out = rows_list(r);
  const int n = static_cast<int>(r.O.size());
  out.attr("row.names") = n > 0 ? IntegerVector::create(NA_INTEGER, -n) : IntegerVector(0);
  out.attr("class") = "data.frame";
  return out;
}

// counts: R's images x types matrix, image-major
std::vector<double> image_major(const NumericMatrix& counts) {
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  return cnt;
}

List design_list(const DesignResult& d) {
  return List::create(_["ok"] = d.ok, _["reason"] = d.reason, _["theta"] = nv(d.theta), _["estimate"] = d.estimate,
                      _["se"] = d.se, _["df"] = d.df, _["p"] = d.p, _["tau2"] = d.tau2,
                      _["influence"] = nv(d.influence));
}

}  // namespace

// [[Rcpp::export]]
NumericVector stats_pt_two_sided(NumericVector t, double df) {
  NumericVector out(t.size());
  for (R_xlen_t i = 0; i < t.size(); ++i) out[i] = pt_two_sided(t[i], df);
  return out;
}

// [[Rcpp::export]]
NumericVector stats_norm_quantile(NumericVector p) {
  NumericVector out(p.size());
  for (R_xlen_t i = 0; i < p.size(); ++i) out[i] = norm_quantile(p[i]);
  return out;
}

// [[Rcpp::export]]
List stats_excess_image_rows(NumericVector totals, NumericVector out_sq_totals, NumericMatrix counts,
                                  int from, int to, bool knn, NumericMatrix psi, NumericVector self_var,
                                  Nullable<IntegerVector> image_group = R_NilValue) {
  // self_var: types x images (dataset_self_variances), needed for a self-pair; may be empty otherwise
  // counts: images x types (R layout); psi: types x images, or 0 x 0 for psi = 1; image_group: per image (NA: a
  // missing condition, no row), or NULL for every image
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  std::vector<double> ps;
  if (psi.nrow() > 0) {
    ps.resize(static_cast<std::size_t>(T) * n_images);
    for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) ps[static_cast<std::size_t>(t) * n_images + i] = psi(t, i);
  }
  std::vector<int> ig;
  if (image_group.isNotNull()) ig = ivec(IntegerVector(image_group));
  return rows_to(excess_image_rows(dvec(totals), dvec(out_sq_totals), cnt, T, n_images, from, to, knn, ps,
                                   dvec(self_var), ig));
}

// [[Rcpp::export]]
List stats_allocation_image_rows(NumericVector any_totals, NumericMatrix self_expected, NumericMatrix counts,
                                      int from, int to, NumericMatrix psi, NumericVector self_var,
                                      Nullable<IntegerVector> image_group = R_NilValue) {
  // self_var: types x images (dataset_self_variances), needed for a self-pair; may be empty otherwise
  // self_expected: types x images (dataset_self_any_expected); counts: images x types; psi: types x images, or 0 x 0;
  // image_group: per image (NA: a missing condition, no row and no part in the side), or NULL for every image
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  std::vector<double> ps;
  if (psi.nrow() > 0) {
    ps.resize(static_cast<std::size_t>(T) * n_images);
    for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) ps[static_cast<std::size_t>(t) * n_images + i] = psi(t, i);
  }
  std::vector<int> ig;
  if (image_group.isNotNull()) ig = ivec(IntegerVector(image_group));
  return rows_to(allocation_image_rows(dvec(any_totals), std::vector<double>(self_expected.begin(), self_expected.end()),
                                       cnt, T, n_images, from, to, ps, dvec(self_var), ig));
}

// [[Rcpp::export]]
List stats_weighted_image_rows(NumericMatrix sums, int n_labels, int to, bool allocation,
                               Nullable<IntegerVector> image_group = R_NilValue) {
  // sums: (3 + 4 n_labels) x images (dataset_weighted_excess_sums); to: 0-based label
  std::vector<int> ig;
  if (image_group.isNotNull()) ig = ivec(IntegerVector(image_group));
  return rows_to(weighted_image_rows(sums.begin(), sums.ncol(), n_labels, to, allocation, ig.empty() ? nullptr : ig.data()));
}

// [[Rcpp::export]]
List stats_kontextual_image_rows_weighted(NumericMatrix sums) {
  // sums: 10 x images (dataset_kontextual_sums_weighted)
  return rows_to(kontextual_image_rows_weighted(std::vector<double>(sums.begin(), sums.end()), sums.ncol()));
}

// [[Rcpp::export]]
List stats_kontextual_image_rows(NumericMatrix sums, NumericMatrix counts, int from, int to, NumericMatrix psi) {
  // sums: 7 x images (dataset_kontextual_sums); counts: images x types; psi: types x images, or 0 x 0
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  std::vector<double> ps;
  if (psi.nrow() > 0) {
    ps.resize(static_cast<std::size_t>(T) * n_images);
    for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) ps[static_cast<std::size_t>(t) * n_images + i] = psi(t, i);
  }
  return rows_to(kontextual_image_rows(std::vector<double>(sums.begin(), sums.end()), cnt, T, n_images, from, to, ps));
}

// [[Rcpp::export]]
List stats_image_heterogeneity_rows(List rows, IntegerVector unit) {
  double tau2 = 0;
  List out = rows_to(image_heterogeneity_rows(rows_from(rows), ivec(unit), &tau2));
  out.attr("tau2_image") = tau2;
  return out;
}

// [[Rcpp::export]]
IntegerVector stats_informative_units(List rows, IntegerVector unit, IntegerVector group, int n_groups) {
  // group: 0 .. n_groups - 1 per row (NA: left out)
  std::vector<int> r = informative_units(rows_from(rows), ivec(unit), ivec(group), n_groups);
  return IntegerVector(r.begin(), r.end());
}

// [[Rcpp::export]]
NumericVector stats_effective_units(List rows, IntegerVector unit, IntegerVector group, int n_groups, double tau2) {
  // group: 0 .. n_groups - 1 per row (NA: left out)
  std::vector<double> r = effective_units(rows_from(rows), ivec(unit), ivec(group), n_groups, tau2);
  return NumericVector(r.begin(), r.end());
}

// [[Rcpp::export]]
List stats_test_rules(bool allocation) {
  TestRules r = test_rules(allocation);
  return List::create(_["pooled_dispersion"] = r.pooled_dispersion, _["min_effective_units"] = r.min_effective_units,
                      _["min_n_eff"] = kMinEffectiveUnits);
}

// [[Rcpp::export]]
NumericMatrix stats_kontextual_clustering(SEXP ptr, IntegerVector from, IntegerVector to, NumericMatrix raw,
                                          NumericMatrix counts, double h) {
  // raw: pairs x images (unweighted REF-TARGET pairs within r); returns types x images
  XPtr<Dataset> d(ptr);
  const int n_images = counts.nrow(), T = counts.ncol(), K = raw.nrow();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T), rw(static_cast<std::size_t>(K) * n_images);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  for (int k = 0; k < K; ++k) for (int i = 0; i < n_images; ++i) rw[static_cast<std::size_t>(k) * n_images + i] = raw(k, i);
  std::vector<double> f = kontextual_clustering_factor(*d, ivec(from), ivec(to), rw, cnt, T, n_images, h);
  NumericMatrix out(T, n_images);
  for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) out(t, i) = f[static_cast<std::size_t>(t) * n_images + i];
  return out;
}

// [[Rcpp::export]]
NumericMatrix stats_label_clustering(SEXP ptr, IntegerVector from, IntegerVector to, NumericMatrix counts,
                                     bool knn, double h, bool allocation = false) {
  XPtr<Dataset> d(ptr);
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  std::vector<double> f = label_clustering_factor(*d, ivec(from), ivec(to), cnt, T, n_images, knn, h, allocation);
  NumericMatrix out(T, n_images);
  for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) out(t, i) = f[static_cast<std::size_t>(t) * n_images + i];
  return out;
}

// [[Rcpp::export]]
List stats_excess_test(List rows, IntegerVector unit, IntegerVector group, int n_units, bool frailty,
                       std::string variance, bool pooled_dispersion = false) {
  ExcessResult r = excess_test(rows_from(rows), ivec(unit), ivec(group), n_units, frailty,
                               variance == "hartung_knapp" ? Variance::HartungKnapp : Variance::CR2, pooled_dispersion);
  return List::create(_["ok"] = r.ok, _["reason"] = r.reason, _["coef_ref"] = r.coef_ref, _["coef_comp"] = r.coef_comp,
                      _["difference"] = r.difference, _["se"] = r.se, _["df"] = r.df, _["p"] = r.p, _["tau2"] = r.tau2,
                      _["influence"] = nv(r.influence), _["unit_summary"] = nv(r.unit_summary),
                      _["unit_info"] = nv(r.unit_info), _["image_weight"] = nv(r.image_weight),
                      _["n_eff"] = nv(r.n_eff), _["n_units"] = nv(r.n_units));
}

// [[Rcpp::export]]
List stats_design_test(List rows, IntegerVector unit, int n_units, NumericMatrix Z, NumericVector contrast,
                       double tau2) {
  const int N = Z.nrow(), q = Z.ncol();
  std::vector<double> z(static_cast<std::size_t>(N) * q);
  for (int i = 0; i < N; ++i) for (int j = 0; j < q; ++j) z[static_cast<std::size_t>(i) * q + j] = Z(i, j);
  return design_list(design_test(rows_from(rows), ivec(unit), n_units, z, q, dvec(contrast), tau2));
}

// [[Rcpp::export]]
List stats_design_tests(List rows, IntegerVector unit, int n_units, NumericMatrix Z, NumericMatrix contrasts,
                        double tau2, bool hartung_knapp, int n_groups = 0) {
  // contrasts: one row per contrast, q columns
  const int N = Z.nrow(), q = Z.ncol(), k = contrasts.nrow();
  std::vector<double> z(static_cast<std::size_t>(N) * q), c(static_cast<std::size_t>(k) * q);
  for (int i = 0; i < N; ++i) for (int j = 0; j < q; ++j) z[static_cast<std::size_t>(i) * q + j] = Z(i, j);
  for (int i = 0; i < k; ++i) for (int j = 0; j < q; ++j) c[static_cast<std::size_t>(i) * q + j] = contrasts(i, j);
  std::vector<DesignResult> r = design_tests(rows_from(rows), ivec(unit), n_units, z, q, c, k, tau2, hartung_knapp,
                                             n_groups);
  List out(k);
  for (int i = 0; i < k; ++i) out[i] = design_list(r[i]);
  return out;
}

// [[Rcpp::export]]
List stats_availability_test(List rows, IntegerVector unit, IntegerVector group, int n_units, NumericVector x,
                             double tau2) {
  return design_list(availability_test(rows_from(rows), ivec(unit), ivec(group), n_units, dvec(x), tau2));
}

// [[Rcpp::export]]
List stats_cox_fit(NumericVector time, IntegerVector event, NumericMatrix X) {
  const int n = X.nrow(), p = X.ncol();
  std::vector<double> x(static_cast<std::size_t>(n) * p);
  for (int i = 0; i < n; ++i) for (int j = 0; j < p; ++j) x[static_cast<std::size_t>(i) * p + j] = X(i, j);
  CoxResult r = cox_fit(dvec(time), ivec(event), x, p);
  return List::create(_["ok"] = r.ok, _["beta"] = nv(r.beta), _["se"] = nv(r.se), _["p"] = nv(r.p),
                      _["martingale"] = nv(r.martingale), _["loglik"] = r.loglik);
}

// [[Rcpp::export]]
List stats_survival_test(List rows, IntegerVector unit, int n_units, NumericVector martingale, NumericVector time,
                         IntegerVector event, NumericVector x) {
  SurvivalResult r = survival_test(rows_from(rows), ivec(unit), n_units, dvec(martingale), dvec(time), ivec(event), dvec(x));
  return List::create(_["ok"] = r.ok, _["reason"] = r.reason, _["score_coef"] = r.score_coef, _["score_se"] = r.score_se,
                      _["score_df"] = r.score_df, _["score_p"] = r.score_p, _["log_hr_sd"] = r.log_hr_sd,
                      _["hr_sd"] = r.hr_sd, _["hr_se"] = r.hr_se, _["hr_p"] = r.hr_p, _["log_hr_unit"] = r.log_hr_unit,
                      _["tau2"] = r.tau2);
}

// [[Rcpp::export]]
double stats_cauchy(NumericVector p) { return cauchy_combine(dvec(p)); }

// [[Rcpp::export]]
List stats_max_t(NumericMatrix influence, NumericVector t, NumericVector df) {
  // influence: units x radii
  std::vector<std::vector<double>> inf(influence.ncol(), std::vector<double>(influence.nrow()));
  for (int k = 0; k < influence.ncol(); ++k) for (int u = 0; u < influence.nrow(); ++u) inf[k][u] = influence(u, k);
  MaxTResult r = max_t(inf, dvec(t), dvec(df));
  return List::create(_["p"] = r.p, _["best"] = r.best + 1);
}

// [[Rcpp::export]]
List stats_davies(NumericMatrix influence, NumericVector t, NumericVector df) {
  // influence: units x radii, radii in increasing order
  std::vector<std::vector<double>> inf(influence.ncol(), std::vector<double>(influence.nrow()));
  for (int k = 0; k < influence.ncol(); ++k) for (int u = 0; u < influence.nrow(); ++u) inf[k][u] = influence(u, k);
  MaxTResult r = davies_max(inf, dvec(t), dvec(df));
  return List::create(_["p"] = r.p, _["best"] = r.best + 1);
}

// [[Rcpp::export]]
double stats_mvn_outside(NumericVector lower, NumericVector upper, NumericMatrix R) {
  const int K = R.nrow();
  std::vector<double> r(static_cast<std::size_t>(K) * K);
  for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) r[static_cast<std::size_t>(i) * K + j] = R(i, j);
  return mvn_outside(dvec(lower), dvec(upper), r, K);
}

// [[Rcpp::export]]
List stats_mvn_outside_qmc(NumericVector lower, NumericVector upper, NumericMatrix R, double abseps = 1e-4,
                           double max_points = 1e5) {
  const int K = R.nrow();
  std::vector<double> r(static_cast<std::size_t>(K) * K);
  for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) r[static_cast<std::size_t>(i) * K + j] = R(i, j);
  MvnResult m = mvn_outside_qmc(dvec(lower), dvec(upper), r, K, abseps, static_cast<long>(max_points));
  return List::create(_["p"] = m.p, _["error"] = m.error, _["points"] = static_cast<double>(m.points));
}

// [[Rcpp::export]]
double stats_mvn_outside_quadrature(NumericVector lower, NumericVector upper, NumericMatrix R) {
  const int K = R.nrow();
  std::vector<double> r(static_cast<std::size_t>(K) * K);
  for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) r[static_cast<std::size_t>(i) * K + j] = R(i, j);
  return mvn_outside_quadrature(dvec(lower), dvec(upper), r, K);
}

// [[Rcpp::export]]
List stats_pair_tests(List a, List b, List self_var, bool allocation, bool knn, NumericMatrix counts, NumericMatrix psi,
                      IntegerVector from, IntegerVector to, IntegerVector image_unit, IntegerVector image_group,
                      int n_units, bool image_heterogeneity, bool frailty, std::string variance, bool adjust,
                      NumericMatrix z, int keep_radius, std::string combine, int n_threads = 1) {
  // a, b: one vector per radius (allocation: any totals and self_expected; count: totals and sq); counts: images x
  // types; psi: types x images, or 0 x 0 for psi = 1; z: images x covariates (centred, NA where missing), or 0 x 0;
  // keep_radius and the radius chosen are 1-based. Per pair and radius: K x P matrices.
  const int K = a.size(), n_images = counts.nrow(), T = counts.ncol();
  if (b.size() != K || self_var.size() != K || K < 1) stop("stats_pair_tests: one a, b and self_var per radius");
  std::vector<NumericVector> av, bv, cv;
  PairTestInput in;
  for (int k = 0; k < K; ++k) {
    av.push_back(as<NumericVector>(a[k])); bv.push_back(as<NumericVector>(b[k])); cv.push_back(as<NumericVector>(self_var[k]));
    in.a.push_back(av.back().begin()); in.b.push_back(bv.back().begin());
    in.self_var.push_back(cv.back().size() ? cv.back().begin() : nullptr);
  }
  std::vector<double> cnt = image_major(counts), ps, zz;
  if (psi.nrow() > 0) {
    ps.resize(static_cast<std::size_t>(T) * n_images);
    for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) ps[static_cast<std::size_t>(t) * n_images + i] = psi(t, i);
    in.psi = ps.data();
  }
  if (z.ncol() > 0) {
    zz = image_major(z);
    in.z = zz.data(); in.n_z = z.ncol();
  }
  in.allocation = allocation; in.knn = knn; in.counts = cnt.data(); in.n_types = T; in.n_images = n_images;
  in.from = ivec(from); in.to = ivec(to);
  in.image_unit = image_unit.begin(); in.image_group = image_group.begin(); in.n_units = n_units;
  in.image_heterogeneity = image_heterogeneity; in.frailty = frailty;
  in.variance = variance == "hartung_knapp" ? Variance::HartungKnapp : Variance::CR2;
  in.adjust = adjust; in.keep_radius = keep_radius - 1; in.n_threads = n_threads;
  in.combine = combine == "maxT" ? Combine::MaxT : combine == "davies" ? Combine::Davies
             : combine == "cauchy" ? Combine::Cauchy : Combine::None;
  PairTests r = pair_tests(in);
  const int P = r.n_pairs, ne = 1 + r.n_z;
  auto km = [&](const std::vector<double>& v) { NumericVector x(v.begin(), v.end()); x.attr("dim") = Dimension(K, P); return x; };
  auto ki = [&](const std::vector<int>& v) { IntegerVector x(v.begin(), v.end()); x.attr("dim") = Dimension(K, P); return x; };
  auto ke = [&](const std::vector<double>& v) { NumericVector x(v.begin(), v.end()); x.attr("dim") = Dimension(ne, K, P); return x; };
  CharacterVector reason(r.design_reason.begin(), r.design_reason.end());
  reason.attr("dim") = Dimension(K, P);
  List rows(P), weight(P);
  for (int p = 0; p < P; ++p) { rows[p] = rows_list(r.rows[p]); weight[p] = nv(r.image_weight[p]); }
  IntegerVector best(P);
  for (int p = 0; p < P; ++p) best[p] = r.best[p] < 0 ? NA_INTEGER : r.best[p] + 1;
  List out;   // (more elements than List::create takes)
  out.push_back(ki(r.ok), "ok"); out.push_back(ki(r.side), "side"); out.push_back(ki(r.n_rows), "n_rows");
  out.push_back(km(r.tau2_image), "tau2_image"); out.push_back(km(r.coef_ref), "coef_ref");
  out.push_back(km(r.coef_comp), "coef_comp"); out.push_back(km(r.difference), "difference");
  out.push_back(km(r.se), "se"); out.push_back(km(r.df), "df"); out.push_back(km(r.p), "p");
  out.push_back(km(r.tau2), "tau2"); out.push_back(km(r.u_difference), "u_difference");
  out.push_back(km(r.u_se), "u_se"); out.push_back(km(r.u_df), "u_df"); out.push_back(km(r.u_p), "u_p");
  out.push_back(ki(r.design), "design"); out.push_back(ki(r.abundance), "abundance");
  CharacterVector not_tested(r.not_tested.begin(), r.not_tested.end());
  not_tested.attr("dim") = Dimension(K, P);
  out.push_back(reason, "design_reason"); out.push_back(not_tested, "not_tested"); out.push_back(ke(r.effect), "effect"); out.push_back(ke(r.effect_p), "effect_p");
  out.push_back(rows, "rows"); out.push_back(weight, "image_weight"); out.push_back(best, "best");
  out.push_back(nv(r.comb_p), "comb_p"); out.push_back(nv(r.comb_p_best), "comb_p_best");
  out.push_back(nv(r.comb_u_p), "comb_u_p");
  return out;
}

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

ImageRows rows_from(const DataFrame& d) {
  ImageRows r;
  IntegerVector img = d["img"];
  NumericVector O = d["O"], E = d["E"], n = d["n"], v = d["v"];
  r.image = ivec(img); r.O = dvec(O); r.E = dvec(E); r.n = dvec(n); r.v = dvec(v);
  return r;
}

DataFrame rows_to(const ImageRows& r) {
  return DataFrame::create(_["img"] = IntegerVector(r.image.begin(), r.image.end()), _["O"] = nv(r.O),
                           _["E"] = nv(r.E), _["n"] = nv(r.n), _["v"] = nv(r.v));
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
DataFrame stats_excess_image_rows(NumericVector totals, NumericVector out_sq_totals, NumericMatrix counts,
                                  int from, int to, bool knn, NumericMatrix psi) {
  // counts: images x types (R layout); psi: types x images, or 0 x 0 for psi = 1
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  std::vector<double> ps;
  if (psi.nrow() > 0) {
    ps.resize(static_cast<std::size_t>(T) * n_images);
    for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) ps[static_cast<std::size_t>(t) * n_images + i] = psi(t, i);
  }
  return rows_to(excess_image_rows(dvec(totals), dvec(out_sq_totals), cnt, T, n_images, from, to, knn, ps));
}

// [[Rcpp::export]]
DataFrame stats_allocation_image_rows(NumericVector any_totals, NumericMatrix self_expected, NumericMatrix counts,
                                      int from, int to, NumericMatrix psi) {
  // self_expected: types x images (dataset_self_any_expected); counts: images x types; psi: types x images, or 0 x 0
  const int n_images = counts.nrow(), T = counts.ncol();
  std::vector<double> cnt(static_cast<std::size_t>(n_images) * T);
  for (int i = 0; i < n_images; ++i) for (int t = 0; t < T; ++t) cnt[static_cast<std::size_t>(i) * T + t] = counts(i, t);
  std::vector<double> ps;
  if (psi.nrow() > 0) {
    ps.resize(static_cast<std::size_t>(T) * n_images);
    for (int t = 0; t < T; ++t) for (int i = 0; i < n_images; ++i) ps[static_cast<std::size_t>(t) * n_images + i] = psi(t, i);
  }
  return rows_to(allocation_image_rows(dvec(any_totals), std::vector<double>(self_expected.begin(), self_expected.end()),
                                       cnt, T, n_images, from, to, ps));
}

// [[Rcpp::export]]
DataFrame stats_kontextual_image_rows(NumericMatrix sums, NumericMatrix counts, int from, int to, NumericMatrix psi) {
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
List stats_excess_test(DataFrame rows, IntegerVector unit, IntegerVector group, int n_units, bool frailty,
                       std::string variance) {
  ExcessResult r = excess_test(rows_from(rows), ivec(unit), ivec(group), n_units, frailty,
                               variance == "hartung_knapp" ? Variance::HartungKnapp : Variance::CR2);
  return List::create(_["ok"] = r.ok, _["reason"] = r.reason, _["coef_ref"] = r.coef_ref, _["coef_comp"] = r.coef_comp,
                      _["difference"] = r.difference, _["se"] = r.se, _["df"] = r.df, _["p"] = r.p, _["tau2"] = r.tau2,
                      _["influence"] = nv(r.influence), _["unit_summary"] = nv(r.unit_summary),
                      _["unit_info"] = nv(r.unit_info), _["image_weight"] = nv(r.image_weight));
}

// [[Rcpp::export]]
List stats_design_test(DataFrame rows, IntegerVector unit, int n_units, NumericMatrix Z, NumericVector contrast,
                       double tau2) {
  const int N = Z.nrow(), q = Z.ncol();
  std::vector<double> z(static_cast<std::size_t>(N) * q);
  for (int i = 0; i < N; ++i) for (int j = 0; j < q; ++j) z[static_cast<std::size_t>(i) * q + j] = Z(i, j);
  return design_list(design_test(rows_from(rows), ivec(unit), n_units, z, q, dvec(contrast), tau2));
}

// [[Rcpp::export]]
List stats_design_tests(DataFrame rows, IntegerVector unit, int n_units, NumericMatrix Z, NumericMatrix contrasts,
                        double tau2, bool hartung_knapp) {
  // contrasts: one row per contrast, q columns
  const int N = Z.nrow(), q = Z.ncol(), k = contrasts.nrow();
  std::vector<double> z(static_cast<std::size_t>(N) * q), c(static_cast<std::size_t>(k) * q);
  for (int i = 0; i < N; ++i) for (int j = 0; j < q; ++j) z[static_cast<std::size_t>(i) * q + j] = Z(i, j);
  for (int i = 0; i < k; ++i) for (int j = 0; j < q; ++j) c[static_cast<std::size_t>(i) * q + j] = contrasts(i, j);
  std::vector<DesignResult> r = design_tests(rows_from(rows), ivec(unit), n_units, z, q, c, k, tau2, hartung_knapp);
  List out(k);
  for (int i = 0; i < k; ++i) out[i] = design_list(r[i]);
  return out;
}

// [[Rcpp::export]]
List stats_availability_test(DataFrame rows, IntegerVector unit, IntegerVector group, int n_units, NumericVector x,
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
List stats_survival_test(DataFrame rows, IntegerVector unit, int n_units, NumericVector martingale, NumericVector time,
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
double stats_mvn_outside(NumericVector lower, NumericVector upper, NumericMatrix R) {
  const int K = R.nrow();
  std::vector<double> r(static_cast<std::size_t>(K) * K);
  for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j) r[static_cast<std::size_t>(i) * K + j] = R(i, j);
  return mvn_outside(dvec(lower), dvec(upper), r, K);
}

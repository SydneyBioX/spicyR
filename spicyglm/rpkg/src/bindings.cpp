// R bindings for the spicyGLM numeric core, mirroring cpp/src/bindings.cpp
// (pybind11). The core itself is untouched and shared between the two.
#include <Rcpp.h>

#include <limits>
#include <string>
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

// PairFit -> list(beta, v_hat, df), matching pair_fit_dict() without diagnostics
List pair_fit_list(const PairFit& p) {
  return List::create(
      _["beta"] = NumericVector::create(p.fit.beta[0], p.fit.beta[1]),
      _["v_hat"] = p.naive ? p.v_naive : p.cr2.v_hat,
      _["df"] = p.naive ? R_PosInf : p.cr2.df);
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
List dataset_binomial_model_data(SEXP ptr, int from, int to) {
  XPtr<Dataset> d(ptr);
  BinomialModelData md = d->binomial_model_data(from, to);
  return List::create(_["row"] = IntegerVector(md.row.begin(), md.row.end()),
                      _["image"] = IntegerVector(md.image.begin(), md.image.end()),
                      _["n"] = IntegerVector(md.n.begin(), md.n.end()),
                      _["p0"] = NumericVector(md.p0.begin(), md.p0.end()));
}

// [[Rcpp::export]]
List fit_pair_poisson_cpp(IntegerVector cluster, IntegerVector image, IntegerVector group,
                          IntegerVector n, NumericVector density, std::string estimator,
                          std::string variance) {
  return pair_fit_list(fit_pair_poisson(ivec(cluster), ivec(image), ivec(group), ivec(n),
                                        dvec(density), estimator, variance, false));
}

// [[Rcpp::export]]
List fit_pair_binomial_cpp(IntegerVector cluster, IntegerVector image, IntegerVector group,
                           IntegerVector n, int k, NumericVector p0, std::string estimator,
                           std::string variance) {
  return pair_fit_list(fit_pair_binomial(ivec(cluster), ivec(image), ivec(group), ivec(n), k,
                                         dvec(p0), estimator, variance));
}

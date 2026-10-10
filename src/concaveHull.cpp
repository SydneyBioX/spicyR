// R binding of the C++ port of concaveman (concaveman.cpp, the same code as lisaClust's src/core).
#include <Rcpp.h>

#include <array>
#include <vector>

#include "concaveman.hpp"

// [[Rcpp::export(.concaveHull)]]
Rcpp::NumericMatrix concaveHullR(Rcpp::NumericVector x, Rcpp::NumericVector y, double concavity,
                                 double lengthThreshold) {
  std::vector<std::array<double, 2>> ring = concaveman::concaveman_r(
      std::vector<double>(x.begin(), x.end()), std::vector<double>(y.begin(), y.end()), concavity, lengthThreshold);
  Rcpp::NumericMatrix out(ring.size(), 2);
  for (std::size_t i = 0; i < ring.size(); ++i) {
    out(i, 0) = ring[i][0];
    out(i, 1) = ring[i][1];
  }
  return out;
}

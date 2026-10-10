// A C++ port of concaveman, giving the same polygons as the R package concaveman 1.2.0 (see concaveman.cpp).
#pragma once

#include <array>
#include <vector>

namespace concaveman {

// The concave hull of the points (x, y) as a closed ring (the first point repeated at the end), as the
// JavaScript concaveman(points, concavity, lengthThreshold) returns it.
std::vector<std::array<double, 2>> concave_hull(const std::vector<double>& x, const std::vector<double>& y,
                                                double concavity, double length_threshold);

// As concaveman::concaveman(cbind(x, y), concavity, length_threshold) in R: the coordinates rounded to
// 4 decimals and the parameters to 15 significant digits, as its R interface passes them to JavaScript.
std::vector<std::array<double, 2>> concaveman_r(const std::vector<double>& x, const std::vector<double>& y,
                                                double concavity, double length_threshold);

namespace detail {  // exposed for testing
double js_exp(double x);
double js_log(double x);
int orient3_sign(double ax, double ay, double bx, double by, double cx, double cy);
double json4(double v);
}  // namespace detail

}  // namespace concaveman

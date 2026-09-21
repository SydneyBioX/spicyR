// Plain C++ cores of the spatial statistic, free of R objects so that
// getPairwiseCpp() can run them on several threads.
#ifndef SPICYR_SPICYCORE_H
#define SPICYR_SPICYCORE_H

#include <cstddef>
#include <vector>

namespace spicy {

struct Pt {
  double x, y;
};
using Ring = std::vector<Pt>;

// Area of each disc of radius r (the npoly-gon spatstat.geom::disc() builds)
// intersected with the window `rings` (src/borderEdge.cpp).
void discAreas(const double* x, const double* y, std::size_t n, double r, int npoly,
               const std::vector<Ring>& rings, double* out);

// The pairwise L statistic of one image (src/inhomL.cpp). type is 0-based;
// edge is n x (Rs.size() - 1), column-major. Writes a K x K column-major
// matrix to out, NaN where a type combination has no value.
void inhomLCore(const double* x, const double* y, const int* type, int n, int K,
                const std::vector<double>& Rs, const std::vector<double>& labelVal,
                const std::vector<char>& isFrom, const std::vector<char>& isTo,
                const double* wt, const std::vector<double>& lam, double area,
                const double* edge, bool edgeCorrect, double* out);

}  // namespace spicy

#endif

#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <thread>
#include <vector>

#include "spicyCore.h"

// getPairwise() for square and convex windows without density weighting:
// every image's inhomLPair() in one call, with images spread over threads.
// The threads touch no R objects: inputs are copied into C++ containers first
// and the result is written into a plain buffer.
//
// Per image, as inhomLPair() does it:
// - the window is the bounding rectangle or the convex hull of all its cells;
// - the radii are capped at min(width, height) / 2.01, with a leading 0;
// - only cells whose type is in `from` or `to` are kept;
// - a kept cell is a border cell for radius r if it lies within r of the
//   window's boundary (spatstat's border(W, r) region, whose extra 2-unit
//   outer band holds no cells), and its edge weight is the fraction of its
//   disc inside the window. A cell at exactly distance r counts, as it does
//   in spatstat for rectangles and for all but a few cells of polygons (whose
//   border region spatstat builds with rounded polygon offsetting);
// - the statistic comes from spicy::inhomLCore(), and the row starts at
//   -sum(radii), with NA for type pairs absent from the image unless
//   includeZeroCells.

namespace {

using spicy::Pt;
using spicy::Ring;

// Convex hull, anticlockwise (Andrew's monotone chain).
Ring convexHull(std::vector<Pt> p) {
  std::sort(p.begin(), p.end(), [](const Pt& a, const Pt& b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
  auto cross = [](const Pt& o, const Pt& a, const Pt& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
  };
  Ring h(2 * p.size());
  std::size_t k = 0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    while (k >= 2 && cross(h[k - 2], h[k - 1], p[i]) <= 0) --k;
    h[k++] = p[i];
  }
  for (std::size_t i = p.size() - 1, t = k + 1; i-- > 0;) {
    while (k >= t && cross(h[k - 2], h[k - 1], p[i]) <= 0) --k;
    h[k++] = p[i];
  }
  h.resize(k > 1 ? k - 1 : k);
  return h;
}

double ringArea(const Ring& r) {
  double a = 0;
  for (std::size_t i = 0, j = r.size() - 1; i < r.size(); j = i++) a += (r[j].x + r[i].x) * (r[i].y - r[j].y);
  return a / 2;
}

// as.numeric(as.character(v)): R prints doubles with 15 significant digits.
double roundTrip15(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.15g", v);
  return std::strtod(buf, nullptr);
}

struct Input {
  const double *x, *y;
  const int* type;  // 0-based
  std::vector<int> offset;
  int K;
  std::vector<double> RsUser;  // sorted
  bool square, edgeCorrect, includeZeroCells;
  std::vector<char> isFrom, isTo;
  std::vector<int> labI, labJ;  // 0-based type of each output column
};

void image(const Input& in, int img, double* row, std::size_t stride, double na) {
  const int s = in.offset[img], e = in.offset[img + 1], nAll = e - s;
  const std::size_t nLab = in.labI.size();
  auto fillNA = [&]() { for (std::size_t l = 0; l < nLab; ++l) row[l * stride] = na; };
  if (nAll == 0) { fillNA(); return; }

  double x0 = in.x[s], x1 = in.x[s], y0 = in.y[s], y1 = in.y[s];
  for (int i = s; i < e; ++i) {
    x0 = std::min(x0, in.x[i]); x1 = std::max(x1, in.x[i]);
    y0 = std::min(y0, in.y[i]); y1 = std::max(y1, in.y[i]);
  }
  Ring W;
  if (in.square) {
    W = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  } else {
    std::vector<Pt> p(nAll);
    for (int i = 0; i < nAll; ++i) p[i] = {in.x[s + i], in.y[s + i]};
    W = convexHull(p);
  }
  if (W.size() < 3) { fillNA(); return; }
  const double area = ringArea(W);
  if (!(area > 0)) { fillNA(); return; }

  // Rs <- unique(pmin(c(0, sort(Rs)), maxR))
  const double maxR = std::min(x1 - x0, y1 - y0) / 2.01;
  std::vector<double> Rs;
  for (std::size_t k = 0; k <= in.RsUser.size(); ++k) {
    const double v = std::min(k == 0 ? 0.0 : in.RsUser[k - 1], maxR);
    if (std::find(Rs.begin(), Rs.end(), v) == Rs.end()) Rs.push_back(v);
  }
  const int nb = static_cast<int>(Rs.size()) - 1;
  if (nb < 1) { fillNA(); return; }
  std::vector<double> labelVal(nb);
  for (int k = 0; k < nb; ++k) labelVal[k] = roundTrip15(Rs[k + 1]);

  // cells whose type is in from or to
  std::vector<double> x, y;
  std::vector<int> t;
  std::vector<double> count(in.K, 0.0);
  for (int i = s; i < e; ++i) {
    const int ti = in.type[i];
    if (!in.isFrom[ti] && !in.isTo[ti]) continue;
    x.push_back(in.x[i]); y.push_back(in.y[i]); t.push_back(ti);
    count[ti] += 1;
  }
  const int n = static_cast<int>(x.size());
  if (n == 0) { fillNA(); return; }
  std::vector<double> lam(in.K);
  for (int k = 0; k < in.K; ++k) lam[k] = count[k] / area;

  std::vector<double> edge(static_cast<std::size_t>(n) * nb, 1.0);
  if (in.edgeCorrect) {
    // distance to the window boundary; for a convex window and a point in it,
    // the smallest distance to an edge's line
    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    const std::size_t m = W.size();
    for (std::size_t j = 0; j < m; ++j) {
      const Pt &p = W[j], &q = W[(j + 1) % m];
      const double len = std::hypot(q.x - p.x, q.y - p.y);
      if (len == 0) continue;
      for (int i = 0; i < n; ++i)
        dist[i] = std::min(dist[i], ((q.x - p.x) * (y[i] - p.y) - (q.y - p.y) * (x[i] - p.x)) / len);
    }
    const std::vector<Ring> rings(1, W);
    std::vector<double> bx, by, a;
    std::vector<int> idx;
    for (int r = 0; r < nb; ++r) {
      const double rad = Rs[r + 1];
      // distances within rounding of r are ties
      const double tol = 1e-9 * std::max(1.0, rad);
      bx.clear(); by.clear(); idx.clear();
      for (int i = 0; i < n; ++i)
        if (dist[i] <= rad + tol) { bx.push_back(x[i]); by.push_back(y[i]); idx.push_back(i); }
      a.resize(idx.size());
      spicy::discAreas(bx.data(), by.data(), idx.size(), rad, 128, rings, a.data());
      for (std::size_t k = 0; k < idx.size(); ++k)
        edge[idx[k] + static_cast<std::size_t>(r) * n] = a[k] / (M_PI * rad * rad);
    }
  }

  std::vector<double> L(static_cast<std::size_t>(in.K) * in.K);
  const std::vector<double> wt(n, 1.0);
  spicy::inhomLCore(x.data(), y.data(), t.data(), n, in.K, Rs, labelVal, in.isFrom, in.isTo,
                    wt.data(), lam, area, edge.data(), in.edgeCorrect, L.data());

  double sumRs = 0;
  for (double v : Rs) sumRs += v;
  for (std::size_t l = 0; l < nLab; ++l) {
    const int I = in.labI[l], J = in.labJ[l];
    double v = -sumRs;
    if (!in.includeZeroCells && (count[I] == 0 || count[J] == 0)) v = na;
    const double Lv = L[I + static_cast<std::size_t>(J) * in.K];
    if (!std::isnan(Lv)) v = Lv;
    row[l * stride] = v;
  }
}

}  // namespace

// offset: 0-based start of each image's cells (cells sorted by image), plus
// the total. labI/labJ: 1-based types of each output column. Returns the
// images x labels matrix getPairwise() returns.
// [[Rcpp::export]]
Rcpp::NumericMatrix getPairwiseCpp(Rcpp::NumericVector x, Rcpp::NumericVector y,
                                   Rcpp::IntegerVector type, Rcpp::IntegerVector offset,
                                   int nTypes, Rcpp::NumericVector Rs, bool square,
                                   Rcpp::LogicalVector isFrom, Rcpp::LogicalVector isTo,
                                   Rcpp::IntegerVector labI, Rcpp::IntegerVector labJ,
                                   bool edgeCorrect, bool includeZeroCells, int nThreads) {
  Input in;
  in.x = x.begin();
  in.y = y.begin();
  std::vector<int> t0(type.size());
  for (R_xlen_t i = 0; i < type.size(); ++i) t0[i] = type[i] - 1;
  in.type = t0.data();
  in.offset.assign(offset.begin(), offset.end());
  in.K = nTypes;
  in.RsUser.assign(Rs.begin(), Rs.end());
  std::sort(in.RsUser.begin(), in.RsUser.end());
  in.square = square;
  in.edgeCorrect = edgeCorrect;
  in.includeZeroCells = includeZeroCells;
  in.isFrom.assign(isFrom.begin(), isFrom.end());
  in.isTo.assign(isTo.begin(), isTo.end());
  for (R_xlen_t l = 0; l < labI.size(); ++l) { in.labI.push_back(labI[l] - 1); in.labJ.push_back(labJ[l] - 1); }

  const int nImg = static_cast<int>(in.offset.size()) - 1;
  const std::size_t nLab = in.labI.size();
  std::vector<double> out(static_cast<std::size_t>(nImg) * nLab);
  const double na = NA_REAL;
  std::atomic<int> next(0);
  auto work = [&]() {
    for (int img = next++; img < nImg; img = next++) image(in, img, out.data() + img, nImg, na);
  };
  const int nt = std::max(1, std::min(nThreads, nImg));
  std::vector<std::thread> pool;
  for (int k = 1; k < nt; ++k) pool.emplace_back(work);
  work();
  for (std::thread& th : pool) th.join();

  Rcpp::NumericMatrix res(nImg, static_cast<int>(nLab));
  std::copy(out.begin(), out.end(), res.begin());
  return res;
}

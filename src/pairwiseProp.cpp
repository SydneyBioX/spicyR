#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <thread>
#include <vector>

// getPairwiseProp(): the fixed-k observed/expected proportion ratio for every
// image in one call, with images spread over threads. The threads touch no R
// objects: inputs are copied into C++ containers first and the result is
// written into a plain buffer.
//
// Per image, as propPair() did it: each cell's k nearest neighbours (any type,
// itself excluded, found exactly as spatstat.geom::nnwhich() finds them); for
// the pair (A, B), piHat is the mean over A cells of the fraction of their
// neighbours that are B, p0 is B's share of all cells, and the value is
// piHat / p0. NA when the image has at most k cells, when A or B is absent,
// or when every cell is B; in those cases 0 instead when includeZeroCells and
// B is present.

namespace {

// Port of R's R_qsort_I (src/main/qsort-body.c, R 4.5.0): Singleton's CACM
// #347 quicksort with Peto's remark. Sorts v[i..j] ascending (1-indexed) and
// applies the same permutation to I. It is not stable, so it is ported exactly:
// spatstat's neighbour search depends on the order it leaves tied values in.
void r_qsort_I(double* v, int* I, int i, int j) {
  int il[40], iu[40];
  double vt, vtt;
  double R = 0.375;
  int ii, ij, k, l, m;
  int it, tt;

  --v;
  --I;
  ii = i;
  m = 1;

L10:
  if (i < j) {
    if (R < 0.5898437) R += 0.0390625; else R -= 0.21875;
  L20:
    k = i;
    ij = i + static_cast<int>((j - i) * R);
    it = I[ij];
    vt = v[ij];
    if (v[i] > vt) {
      I[ij] = I[i]; I[i] = it; it = I[ij];
      v[ij] = v[i]; v[i] = vt; vt = v[ij];
    }
    l = j;
    if (v[j] < vt) {
      I[ij] = I[j]; I[j] = it; it = I[ij];
      v[ij] = v[j]; v[j] = vt; vt = v[ij];
      if (v[i] > vt) {
        I[ij] = I[i]; I[i] = it; it = I[ij];
        v[ij] = v[i]; v[i] = vt; vt = v[ij];
      }
    }
    for (;;) {
      do l--; while (v[l] > vt);
      tt = I[l];
      vtt = v[l];
      do k++; while (v[k] < vt);
      if (k > l) break;
      I[l] = I[k]; I[k] = tt;
      v[l] = v[k]; v[k] = vtt;
    }
    m++;
    if (l - i <= j - k) {
      il[m] = k;
      iu[m] = j;
      j = l;
    } else {
      il[m] = i;
      iu[m] = l;
      i = k;
    }
  } else {
  L80:
    if (m == 1) return;
    i = il[m];
    j = iu[m];
    m--;
  }

  if (j - i > 10) goto L20;
  if (i == ii) goto L10;

  --i;
L100:
  do {
    ++i;
    if (i == j) goto L80;
    it = I[i + 1];
    vt = v[i + 1];
  } while (v[i] <= vt);
  k = i;
  do {
    I[k + 1] = I[k];
    v[k + 1] = v[k];
    --k;
  } while (vt < v[k]);
  I[k + 1] = it;
  v[k + 1] = vt;
  goto L100;
}

// k nearest neighbours of each point in one image, as spatstat.geom 3.5-0's
// nnwhich(): order points by y with sort.list(method = "quick"), then scan
// backward and forward from each point keeping a candidate only if it is
// strictly closer than the current k-th (src/knndist.h). Ties at the k-th
// distance are therefore resolved exactly as in R. Needs n > k.
void image_knn(const double* x, const double* y, int n, int k, int row_offset, int* out) {
  std::vector<double> ys(y, y + n);
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 1);
  r_qsort_I(ys.data(), order.data(), 1, n);
  for (int& o : order) --o;  // 0-based original index of the i-th point by y
  std::vector<double> xs(n);
  for (int i = 0; i < n; ++i) { xs[i] = x[order[i]]; ys[i] = y[order[i]]; }

  const double huge = std::sqrt(std::numeric_limits<double>::max());
  const double hu2 = huge * huge;
  std::vector<double> d2min(k);
  std::vector<int> which(k);
  const int last = k - 1;

  auto consider = [&](int j, double dy2, double xi, double& d2minK) {
    double dx = xs[j] - xi;
    double d2 = dx * dx + dy2;
    if (d2 < d2minK) {
      d2min[last] = d2;
      which[last] = j;
      for (int s = last; s > 0 && d2min[s] < d2min[s - 1]; --s) {
        std::swap(d2min[s], d2min[s - 1]);
        std::swap(which[s], which[s - 1]);
      }
      d2minK = d2min[last];
    }
  };

  for (int i = 0; i < n; ++i) {
    double d2minK = hu2;
    std::fill(d2min.begin(), d2min.end(), hu2);
    std::fill(which.begin(), which.end(), -1);
    double xi = xs[i], yi = ys[i];
    for (int left = i - 1; left >= 0; --left) {
      double dy = yi - ys[left], dy2 = dy * dy;
      if (dy2 > d2minK) break;
      consider(left, dy2, xi, d2minK);
    }
    for (int right = i + 1; right < n; ++right) {
      double dy = ys[right] - yi, dy2 = dy * dy;
      if (dy2 > d2minK) break;
      consider(right, dy2, xi, d2minK);
    }
    int* row = out + static_cast<std::size_t>(order[i]) * k;
    for (int s = 0; s < k; ++s) row[s] = row_offset + order[which[s]];
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix getPairwisePropCpp(Rcpp::NumericVector x, Rcpp::NumericVector y,
                                       Rcpp::IntegerVector type, Rcpp::IntegerVector offset,
                                       int nTypes, int k, Rcpp::IntegerVector labI,
                                       Rcpp::IntegerVector labJ, bool includeZeroCells,
                                       int nThreads) {
  // type, labI and labJ are 1-based; a label of NA is a type absent from the data.
  const std::vector<double> xs(x.begin(), x.end()), ys(y.begin(), y.end());
  std::vector<int> ty(type.size());
  for (R_xlen_t i = 0; i < type.size(); ++i) ty[i] = type[i] - 1;
  const std::vector<int> off(offset.begin(), offset.end());
  std::vector<int> li(labI.size()), lj(labJ.size());
  for (R_xlen_t l = 0; l < labI.size(); ++l) {
    li[l] = labI[l] == NA_INTEGER ? -1 : labI[l] - 1;
    lj[l] = labJ[l] == NA_INTEGER ? -1 : labJ[l] - 1;
  }

  const int nImg = static_cast<int>(off.size()) - 1;
  const std::size_t nLab = li.size();
  const double na = NA_REAL;
  std::vector<double> out(static_cast<std::size_t>(nImg) * nLab, na);

  auto image = [&](int img) {
    const int start = off[img], n = off[img + 1] - start;
    if (n <= k) return;
    std::vector<int> nn(static_cast<std::size_t>(n) * k);
    image_knn(xs.data() + start, ys.data() + start, n, k, 0, nn.data());
    const int* t = ty.data() + start;
    // hits[a * nTypes + b]: B neighbours summed over the A cells
    std::vector<double> hits(static_cast<std::size_t>(nTypes) * nTypes, 0.0);
    std::vector<int> count(nTypes, 0);
    for (int i = 0; i < n; ++i) {
      ++count[t[i]];
      double* row = hits.data() + static_cast<std::size_t>(t[i]) * nTypes;
      for (int s = 0; s < k; ++s) row[t[nn[static_cast<std::size_t>(i) * k + s]]] += 1.0;
    }
    for (std::size_t l = 0; l < nLab; ++l) {
      const int A = li[l], B = lj[l];
      const int nA = A < 0 ? 0 : count[A], nB = B < 0 ? 0 : count[B];
      double& o = out[l * nImg + img];
      if (nA == 0 || nB == 0 || nB == n) {
        if (includeZeroCells && nB > 0) o = 0.0;
        continue;
      }
      const double piHat = hits[static_cast<std::size_t>(A) * nTypes + B] / nA / k;
      o = piHat / (static_cast<double>(nB) / n);
    }
  };

  std::atomic<int> next(0);
  auto work = [&]() {
    for (int img = next++; img < nImg; img = next++) image(img);
  };
  const int nt = std::max(1, std::min(nThreads, nImg));
  std::vector<std::thread> pool;
  for (int t = 1; t < nt; ++t) pool.emplace_back(work);
  work();
  for (std::thread& th : pool) th.join();

  Rcpp::NumericMatrix res(nImg, static_cast<int>(nLab));
  std::copy(out.begin(), out.end(), res.begin());
  return res;
}

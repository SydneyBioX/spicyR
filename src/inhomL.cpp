#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "spicyCore.h"

// Pairwise L-function statistic for one image: the C++ core of inhomLPair().
//
// Reproduces what the R version did with spatstat's closepairs(), cut(),
// dplyr and data.table:
// - ordered pairs (i, j), i != j, with distance <= max(Rs), counted from type
//   I = type[i] (must be in `from`) to type J = type[j] (must be in `to`);
// - each pair adds wt[j] / e_i / (lam[J] * lam[I] * area), where e_i is the
//   edge-correction weight of cell i (1 without edge correction);
// - data.table's CJ() only crosses values that occur, so a type combination
//   or distance bin with no pairs is skipped rather than counted as zero.
//
// Rs: radii with a leading 0 (already capped and made unique in R).
// labelVal: as.numeric(as.character(Rs[-1])), the bin labels the R code
//   compared against each radius when edge correcting.
// edge: n x (length(Rs) - 1) edge-correction weights, one column per radius.
// Writes a K x K matrix of statistics to out, NaN where the combination has
// no value (the caller then fills in its defaults).
void spicy::inhomLCore(const double* x, const double* y, const int* type, int n, int K,
                       const std::vector<double>& Rs, const std::vector<double>& labelVal,
                       const std::vector<char>& isFrom, const std::vector<char>& isTo,
                       const double* wt, const std::vector<double>& lam, double area,
                       const double* edge, bool edgeCorrect, double* out) {
  const int nb = static_cast<int>(Rs.size()) - 1;  // number of distance bins / radii
  std::fill(out, out + K * K, std::numeric_limits<double>::quiet_NaN());
  if (n == 0 || nb < 1 || Rs[nb] <= 0) return;
  const double rmax = Rs[nb], r2max = rmax * rmax;

  // incl[k * nb + r]: pairs in distance bin k count towards radius r.
  std::vector<char> incl(nb * nb);
  for (int k = 0; k < nb; ++k)
    for (int r = 0; r < nb; ++r) incl[k * nb + r] = labelVal[k] <= Rs[r + 1];

  // Sums per (bin or radius, I, J) and which values occur, as CJ() sees them.
  std::vector<double> S(static_cast<std::size_t>(nb) * K * K, 0.0);
  std::vector<char> seenI(nb * K, 0), seenJ(nb * K, 0), seenBin(nb, 0);

  // Uniform grid with cells of side rmax / 2, so every neighbour of a cell lies
  // in the 5 x 5 block around it. Cells are copied into grid order so that the
  // neighbour scan reads memory sequentially.
  const double h = rmax / 2;
  const double x0 = *std::min_element(x, x + n);
  const double y0 = *std::min_element(y, y + n);
  const double x1 = *std::max_element(x, x + n);
  const double y1 = *std::max_element(y, y + n);
  const int gx = static_cast<int>((x1 - x0) / h) + 1;
  const int gy = static_cast<int>((y1 - y0) / h) + 1;
  std::vector<int> cellOf(n), start(static_cast<std::size_t>(gx) * gy + 1, 0), orig(n);
  for (int i = 0; i < n; ++i) {
    const int cx = std::min(static_cast<int>((x[i] - x0) / h), gx - 1);
    const int cy = std::min(static_cast<int>((y[i] - y0) / h), gy - 1);
    cellOf[i] = cy * gx + cx;
    ++start[cellOf[i] + 1];
  }
  for (std::size_t c = 1; c < start.size(); ++c) start[c] += start[c - 1];
  {
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (int i = 0; i < n; ++i) orig[fill[cellOf[i]]++] = i;
  }
  std::vector<double> sx(n), sy(n), sw(n);
  std::vector<int> st(n);
  for (int s = 0; s < n; ++s) {
    const int i = orig[s];
    sx[s] = x[i]; sy[s] = y[i]; sw[s] = wt[i]; st[s] = type[i];
  }

  const std::vector<char>& from = isFrom;
  const std::vector<char>& to = isTo;
  const std::vector<double>& R = Rs;
  // Radii run in increasing order, so a bin counts towards a suffix of them.
  std::vector<int> firstR(nb, nb);
  for (int k = 0; k < nb; ++k)
    for (int r = nb - 1; r >= 0 && incl[k * nb + r]; --r) firstR[k] = r;
  // Reciprocal edge weights in grid order.
  std::vector<double> ie(edgeCorrect ? static_cast<std::size_t>(n) * nb : 0);
  if (edgeCorrect)
    for (int s = 0; s < n; ++s)
      for (int r = 0; r < nb; ++r) ie[static_cast<std::size_t>(s) * nb + r] = 1.0 / edge[orig[s] + static_cast<std::size_t>(r) * n];

  // Count the pair from reference cell a (type I) to neighbour b (type J).
  auto add = [&](int a, int b, int k) {
    const int I = st[a], J = st[b];
    if (!from[I] || !to[J]) return;
    if (!edgeCorrect) {
      S[(static_cast<std::size_t>(k) * K + I) * K + J] += sw[b];
      seenBin[k] = 1;
      seenI[I] = 1;
      seenJ[J] = 1;
    } else {
      for (int r = firstR[k]; r < nb; ++r) {
        S[(static_cast<std::size_t>(r) * K + I) * K + J] +=
            sw[b] * ie[static_cast<std::size_t>(a) * nb + r];
        seenI[r * K + I] = 1;
        seenJ[r * K + J] = 1;
      }
    }
  };
  // Visit each unordered pair once and count it in both directions.
  auto visit = [&](int si, int from_, int to_) {
    const double xi = sx[si], yi = sy[si];
    for (int sj = from_; sj < to_; ++sj) {
      const double dx = sx[sj] - xi, dy = sy[sj] - yi;
      const double d2 = dx * dx + dy * dy;
      if (d2 > r2max) continue;
      const double d = std::sqrt(d2);
      // cut(d, Rs, include.lowest = TRUE): first bin k with d <= Rs[k + 1]
      int k = 0;
      while (k < nb && d > R[k + 1]) ++k;
      if (k == nb) continue;
      add(si, sj, k);
      add(sj, si, k);
    }
  };

  for (int cy = 0; cy < gy; ++cy) {
    for (int cx = 0; cx < gx; ++cx) {
      const int c = cy * gx + cx;
      for (int si = start[c]; si < start[c + 1]; ++si) {
        // the rest of this grid cell, then the forward half of the 5 x 5 block
        visit(si, si + 1, start[c + 1]);
        if (cx + 1 < gx) visit(si, start[c + 1], start[cy * gx + std::min(cx + 2, gx - 1) + 1]);
        for (int ny = cy + 1; ny <= std::min(cy + 2, gy - 1); ++ny)
          visit(si, start[ny * gx + std::max(cx - 2, 0)],
                start[ny * gx + std::min(cx + 2, gx - 1) + 1]);
      }
    }
  }

  // Density scaling, applied once per (I, J) instead of once per pair.
  for (int b = 0; b < nb; ++b)
    for (int I = 0; I < K; ++I)
      for (int J = 0; J < K; ++J)
        S[(static_cast<std::size_t>(b) * K + I) * K + J] /= lam[J] * lam[I] * area;

  for (int I = 0; I < K; ++I) {
    for (int J = 0; J < K; ++J) {
      if (!edgeCorrect) {
        // Cumulative K over the bins that occur, then sum of L minus sum(Rs).
        if (!seenI[I] || !seenJ[J]) continue;
        double cum = 0, total = 0, sumRs = 0;
        for (int k = 0; k <= nb; ++k) sumRs += Rs[k];
        for (int k = 0; k < nb; ++k) {
          if (!seenBin[k]) continue;
          cum += S[(static_cast<std::size_t>(k) * K + I) * K + J];
          total += std::sqrt(cum / M_PI);
        }
        out[I + J * K] = total - sumRs;
      } else {
        // L(r) - r at each radius where the combination occurs, averaged.
        double total = 0;
        int count = 0;
        for (int r = 0; r < nb; ++r) {
          if (!seenI[r * K + I] || !seenJ[r * K + J]) continue;
          total += std::sqrt(S[(static_cast<std::size_t>(r) * K + I) * K + J] / M_PI) - Rs[r + 1];
          ++count;
        }
        if (count > 0) out[I + J * K] = total / count;
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::NumericMatrix inhomLCpp(Rcpp::NumericVector x, Rcpp::NumericVector y,
                              Rcpp::IntegerVector type, int nTypes,
                              Rcpp::NumericVector Rs, Rcpp::NumericVector labelVal,
                              Rcpp::LogicalVector isFrom, Rcpp::LogicalVector isTo,
                              Rcpp::NumericVector wt, Rcpp::NumericVector lam,
                              double area, Rcpp::NumericMatrix edge,
                              bool edgeCorrect) {
  const int n = x.size();
  std::vector<int> t0(n);
  for (int i = 0; i < n; ++i) t0[i] = type[i] - 1;
  std::vector<char> from(isFrom.begin(), isFrom.end()), to(isTo.begin(), isTo.end());
  Rcpp::NumericMatrix out(nTypes, nTypes);
  spicy::inhomLCore(x.begin(), y.begin(), t0.data(), n, nTypes,
                    std::vector<double>(Rs.begin(), Rs.end()),
                    std::vector<double>(labelVal.begin(), labelVal.end()), from, to,
                    wt.begin(), std::vector<double>(lam.begin(), lam.end()), area,
                    edge.begin(), edgeCorrect, out.begin());
  for (double& v : out) if (std::isnan(v)) v = NA_REAL;
  return out;
}

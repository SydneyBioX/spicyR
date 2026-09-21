#include "spicyglm/core.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace spicyglm {

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

double expit(double eta) {
  return eta >= 0 ? 1.0 / (1.0 + std::exp(-eta)) : std::exp(eta) / (1.0 + std::exp(eta));
}

// Root of the (optionally Jeffreys-penalised) binomial score for one group.
// Cells sharing an offset are pooled: m[j] cells at offset o[j].
double solve_group(double Y, const std::vector<double>& m, const std::vector<double>& o, int k, bool firth) {
  auto score = [&](double beta, double& deriv) {
    double U = Y, I = 0.0, I1 = 0.0, I2 = 0.0;
    for (std::size_t j = 0; j < m.size(); ++j) {
      double p = expit(beta + o[j]);
      double v = k * m[j] * p * (1 - p);
      U -= k * m[j] * p;
      I += v;
      I1 += v * (1 - 2 * p);
      I2 += v * (1 - 6 * p + 6 * p * p);
    }
    if (!firth) { deriv = -I; return U; }
    // U* = U + 0.5 d/dbeta log I
    deriv = -I + 0.5 * (I2 / I - (I1 / I) * (I1 / I));
    return U + 0.5 * I1 / I;
  };

  // bracket the root (the score decreases from + to - in beta)
  double d, lo = -1.0, hi = 1.0;
  while (score(lo, d) < 0) { lo -= 2 * (hi - lo); if (lo < -1e3) throw std::runtime_error("binomial fit: no root"); }
  while (score(hi, d) > 0) { hi += 2 * (hi - lo); if (hi > 1e3) throw std::runtime_error("binomial fit: no root"); }

  double beta = 0.5 * (lo + hi);
  for (int iter = 0; iter < 200; ++iter) {
    double F = score(beta, d);
    if (F > 0) lo = beta; else hi = beta;
    double next = beta - F / d;
    if (!(d < 0) || !(next > lo && next < hi)) next = 0.5 * (lo + hi);  // bisect
    if (std::abs(next - beta) < 1e-13 * std::max(1.0, std::abs(beta)) || hi - lo < 1e-14) return next;
    beta = next;
  }
  throw std::runtime_error("binomial fit did not converge");
}

}  // namespace

std::vector<int> knn_indices(const std::vector<double>& x, const std::vector<double>& y,
                             const std::vector<int>& image_offsets, int k, int n_threads) {
  if (k < 1) throw std::invalid_argument("k must be a positive integer");
  if (n_threads < 1) throw std::invalid_argument("n_threads must be positive");
  std::vector<int> out(x.size() * static_cast<std::size_t>(k), -1);
  int n_images = static_cast<int>(image_offsets.size()) - 1;
  // images write to disjoint slices of `out`, so they can run concurrently
  std::atomic<int> next{0};
  auto worker = [&]() {
    for (int img = next++; img < n_images; img = next++) {
      int start = image_offsets[img], n = image_offsets[img + 1] - start;
      if (n <= k) continue;
      image_knn(x.data() + start, y.data() + start, n, k, start, out.data() + static_cast<std::size_t>(start) * k);
    }
  };
  std::vector<std::thread> threads;
  for (int t = 1; t < std::min(n_threads, std::max(n_images, 1)); ++t) threads.emplace_back(worker);
  worker();
  for (std::thread& t : threads) t.join();
  return out;
}

GlmFit fit_binomial(const std::vector<int>& n, int k, const std::vector<double>& p0,
                    const std::vector<int>& group, const std::string& estimator) {
  if (estimator != "mle" && estimator != "firth") throw std::invalid_argument("estimator must be 'mle' or 'firth'");
  std::array<double, 2> Y{0.0, 0.0};
  std::array<std::unordered_map<double, double>, 2> cells_at_offset;
  std::vector<double> offset(n.size());
  for (std::size_t i = 0; i < n.size(); ++i) {
    offset[i] = std::log(p0[i] / (1 - p0[i]));
    Y[group[i]] += n[i];
    cells_at_offset[group[i]][offset[i]] += 1.0;
  }
  GlmFit fit;
  for (int g = 0; g < 2; ++g) {
    std::vector<double> m, o;
    for (const auto& [off, count] : cells_at_offset[g]) { o.push_back(off); m.push_back(count); }
    fit.beta[g] = solve_group(Y[g], m, o, k, estimator == "firth");
  }
  fit.mu.resize(n.size());
  for (std::size_t i = 0; i < n.size(); ++i) fit.mu[i] = k * expit(fit.beta[group[i]] + offset[i]);
  return fit;
}

PairFit fit_pair_binomial(const std::vector<int>& cluster, const std::vector<int>& image,
                          const std::vector<int>& group, const std::vector<int>& n, int k,
                          const std::vector<double>& p0, const std::string& estimator,
                          const std::string& variance) {
  if (variance != "fast" && variance != "naive") throw std::invalid_argument("variance must be 'fast' or 'naive'");
  PairFit out;
  out.fit = fit_binomial(n, k, p0, group, estimator);
  std::vector<double> var(n.size()), resid(n.size());
  for (std::size_t c = 0; c < n.size(); ++c) {
    double p = out.fit.mu[c] / k;
    var[c] = k * p * (1 - p);
    resid[c] = n[c] - out.fit.mu[c];
  }
  if (variance == "naive") {
    out.naive = true;
    out.v_naive = naive_variance(group, var);
    return out;
  }
  out.cr2 = cr2_wald(cluster, image, group, var, resid);
  return out;
}

}  // namespace spicyglm

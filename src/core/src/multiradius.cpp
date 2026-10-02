// Several radii: the Cauchy combination and the max-T test with the sandwich correlation of the
// per-unit CR2 influences (Supplementary Methods, Section 4).
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "spicyglm/stats.hpp"

namespace spicyglm {

namespace {

const double kPi = 3.141592653589793238462643383279502884;

// 20-point Gauss-Legendre nodes and weights on [-1, 1].
const double kGLx[10] = {0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271,
                         0.6360536807265150, 0.7463319064601508, 0.8391169718222188, 0.9122344282513259,
                         0.9639719272779138, 0.9931285991850949};
const double kGLw[10] = {0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766,
                         0.1181945319615184, 0.1019301198172404, 0.0832767415767048, 0.0626720483341091,
                         0.0406014298003869, 0.0176140071391521};

double norm_cdf(double x) { return pnorm_upper(-x); }

// P(X outside [a, b]) for X ~ N(mean, cov) in K dimensions (cov row-major), by conditioning on the
// first coordinate: P(X1 outside) + int_{a1}^{b1} phi(x1) P(rest outside | x1) dx1.
// `nodes`: quadrature nodes per level, chosen so the whole recursion costs about 2e6 evaluations.
double outside(int K, const double* a, const double* b, const std::vector<double>& mean, const std::vector<double>& cov,
               int nodes) {
  double s1 = std::sqrt(std::max(cov[0], 0.0)), m1 = mean[0];
  if (s1 < 1e-10) {  // degenerate: X1 = m1
    if (m1 < a[0] || m1 > b[0]) return 1.0;
    if (K == 1) return 0.0;
    std::vector<double> mr(mean.begin() + 1, mean.end()), cr;
    for (int i = 1; i < K; ++i) for (int j = 1; j < K; ++j) cr.push_back(cov[i * K + j]);
    return outside(K - 1, a + 1, b + 1, mr, cr, nodes);
  }
  double pout = norm_cdf((a[0] - m1) / s1) + pnorm_upper((b[0] - m1) / s1);
  if (K == 1) return pout;
  // conditional covariance (does not depend on x1) and regression of the rest on x1
  std::vector<double> beta(K - 1), cr((K - 1) * (K - 1));
  for (int i = 1; i < K; ++i) beta[i - 1] = cov[i * K] / cov[0];
  for (int i = 1; i < K; ++i) for (int j = 1; j < K; ++j) cr[(i - 1) * (K - 1) + (j - 1)] = cov[i * K + j] - cov[i * K] * cov[j] / cov[0];
  // integrate over [max(a1, m1 - 9 s1), min(b1, m1 + 9 s1)] with panels of width <= s1 / 2
  double lo = std::max(a[0], m1 - 9 * s1), hi = std::min(b[0], m1 + 9 * s1);
  if (hi <= lo) return pout;
  int panels = std::min(std::max(1, nodes / 20), std::max(1, static_cast<int>(std::ceil((hi - lo) / (0.5 * s1)))));
  double h = (hi - lo) / panels, integral = 0;
  std::vector<double> mr(K - 1);
  for (int pnl = 0; pnl < panels; ++pnl) {
    double c = lo + (pnl + 0.5) * h, half = 0.5 * h;
    for (int k = 0; k < 20; ++k) {
      double x = c + (k < 10 ? -1 : 1) * half * kGLx[k % 10], w = half * kGLw[k % 10];
      double dens = std::exp(-0.5 * ((x - m1) / s1) * ((x - m1) / s1)) / (s1 * std::sqrt(2 * kPi));
      for (int i = 0; i < K - 1; ++i) mr[i] = mean[i + 1] + beta[i] * (x - m1);
      integral += w * dens * outside(K - 1, a + 1, b + 1, mr, cr, nodes);
    }
  }
  return std::min(1.0, pout + integral);
}

}  // namespace

double mvn_outside(const std::vector<double>& lower, const std::vector<double>& upper, const std::vector<double>& R, int K) {
  if (K < 1) return 0;
  if (K > 6) throw std::invalid_argument("max-T supports at most 6 radii");
  int nodes = K <= 2 ? 720 : std::max(20, static_cast<int>(std::pow(2e6, 1.0 / (K - 1))));
  return outside(K, lower.data(), upper.data(), std::vector<double>(K, 0.0), R, nodes);
}

double cauchy_combine(const std::vector<double>& p) {
  double T = 0; int K = 0;
  for (double x : p) {
    if (!std::isfinite(x)) continue;
    T += x < 1e-15 ? 1 / (x * kPi) : std::tan((0.5 - x) * kPi);
    ++K;
  }
  if (!K) return std::numeric_limits<double>::quiet_NaN();
  T /= K;
  return T > 1e15 ? 1 / (T * kPi) : 0.5 - std::atan(T) / kPi;
}

MaxTResult max_t(const std::vector<std::vector<double>>& infl, const std::vector<double>& t, const std::vector<double>& df) {
  MaxTResult res;
  const int K0 = static_cast<int>(t.size());
  std::vector<int> use;
  for (int k = 0; k < K0; ++k) if (std::isfinite(t[k]) && std::isfinite(df[k])) use.push_back(k);
  if (use.empty()) { res.p = std::numeric_limits<double>::quiet_NaN(); return res; }
  const std::size_t m = infl[use[0]].size();
  auto cov = [&](int k, int l) { double s = 0; for (std::size_t u = 0; u < m; ++u) s += infl[k][u] * infl[l][u]; return s; };
  // drop a radius whose correlation with an earlier kept one exceeds 0.995 (redundant)
  std::vector<int> keep;
  for (int k : use) {
    bool red = false;
    for (int l : keep) { double r = cov(k, l) / std::sqrt(cov(k, k) * cov(l, l)); if (std::fabs(r) > 0.995) red = true; }
    if (!red) keep.push_back(k);
  }
  const int K = static_cast<int>(keep.size());
  std::vector<double> R(K * K), z(K);
  for (int i = 0; i < K; ++i) for (int j = 0; j < K; ++j)
    R[i * K + j] = i == j ? 1.0 : cov(keep[i], keep[j]) / std::sqrt(cov(keep[i], keep[i]) * cov(keep[j], keep[j]));
  double c = 0;
  for (int i = 0; i < K; ++i) {
    double tail = pt_upper(std::fabs(t[keep[i]]), df[keep[i]]);
    z[i] = -norm_quantile(tail);
    if (z[i] > c) { c = z[i]; res.best = keep[i]; }
  }
  if (K == 1) { res.p = pt_two_sided(t[keep[0]], df[keep[0]]); return res; }
  res.p = mvn_outside(std::vector<double>(K, -c), std::vector<double>(K, c), R, K);
  return res;
}

}  // namespace spicyglm

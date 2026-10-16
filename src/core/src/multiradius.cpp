// Several radii: the Cauchy combination and the max-T test with the sandwich correlation of the
// per-unit CR2 influences (Supplementary Methods, Section 4). The multivariate normal probability of the max-T is
// computed by the quasi-Monte Carlo method of Genz (1992, J. Comput. Graph. Stat. 1, 141-149) and Genz & Bretz
// (2009, Computation of Multivariate Normal and t Probabilities, Springer), written from the published algorithm.
#include <algorithm>
#include <cstdint>
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

double mvn_outside_quadrature(const std::vector<double>& lower, const std::vector<double>& upper, const std::vector<double>& R,
                              int K) {
  if (K < 1) return 0;
  if (K > 6) throw std::invalid_argument("max-T supports at most 6 radii");
  int nodes = K <= 2 ? 720 : std::max(20, static_cast<int>(std::pow(2e6, 1.0 / (K - 1))));
  return outside(K, lower.data(), upper.data(), std::vector<double>(K, 0.0), R, nodes);
}

namespace {

// The standard normal distribution function and its upper tail, from erfc (accurate in both tails).
inline double Phi(double x) { return pnorm_upper(-x); }

// P(lo < Z < hi) for Z ~ N(0, 1), without cancellation in either tail.
double interval_prob(double lo, double hi) {
  if (!(hi > lo)) return 0;
  if (lo > 0) return pnorm_upper(lo) - pnorm_upper(hi);
  if (hi < 0) return Phi(hi) - Phi(lo);
  return 1 - pnorm_upper(-lo) - pnorm_upper(hi);
}

// Korobov lattice rules {n, a}: n points with generating vector (1, a, a^2, a^3, a^4) mod n, n prime and growing by
// about 1.5; a minimises the weighted P2 criterion, by an exhaustive search (spicyR inst/scripts/korobov_search.cpp).
const long long kLattice[][2] = {{67, 16}, {97, 18}, {149, 29}, {223, 38}, {331, 88}, {487, 190}, {733, 87},
                                 {1093, 120}, {1657, 393}, {2459, 414}, {3691, 781}, {5531, 2049}, {8297, 3294},
                                 {12451, 1375}, {18671, 1948}, {27997, 9012}, {41999, 1316}, {62987, 5820}};
const std::size_t kRules = sizeof(kLattice) / sizeof(kLattice[0]);

// SplitMix64 (Steele, Lea & Flood 2014): a fixed, portable stream for the random shifts.
struct SplitMix64 {
  std::uint64_t s;
  double next() {
    std::uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) * 0x1.0p-53;   // in [0, 1)
  }
};

const double kTol = 1e-10;   // conditional variance below this: the coordinate is deterministic given the others

// The integrand after separation of variables: for the cube point w (K - 1 coordinates in [0, 1]), the probability
// that the rest falls outside, 1 - prod_i P(coordinate i inside | the earlier ones), from the log of the product
// so that a small outside probability keeps its relative accuracy.
struct Separated {
  int K = 0;
  std::vector<double> L;     // lower-triangular Cholesky factor of the reordered R, row-major K x K
  std::vector<double> a, b;  // the reordered limits
  std::vector<char> det;     // deterministic coordinate (L[i][i] = 0)
  mutable std::vector<double> y;

  // the first coordinate's probabilities do not depend on the point
  double below0 = 0, above0 = 0, f0 = 0, log0 = 0;

  // P(lo < Z < hi) and the two tails, with two erfc calls and no cancellation in either tail
  static void probs(double lo, double hi, double& below, double& above, double& f) {
    if (lo > 0) {
      const double tl = pnorm_upper(lo);
      above = pnorm_upper(hi); below = 1 - tl; f = tl - above;
    } else if (hi < 0) {
      const double th = Phi(hi);
      below = Phi(lo); above = 1 - th; f = th - below;
    } else {
      below = Phi(lo); above = pnorm_upper(hi); f = 1 - below - above;
    }
  }
  static double log_inside(double below, double above, double f) {
    return f > 0.5 ? std::log1p(-(below + above)) : std::log(f);
  }
  // y_i = Phi^{-1}(below + w_i f), from the nearer tail
  static double draw(double below, double above, double f, double w) {
    const double u = below + w * f;
    if (u < 0.5) return norm_quantile(std::max(u, 1e-300));
    return -norm_quantile(std::max(above + (1 - w) * f, 1e-300));
  }

  void prepare() {
    if (det[0]) return;
    const double l = L[0];
    probs(a[0] / l, b[0] / l, below0, above0, f0);
    log0 = f0 > 0 ? log_inside(below0, above0, f0) : 0;
  }

  double outside(const double* w) const {
    double log_in = 0;
    for (int i = 0; i < K; ++i) {
      double s = 0;
      for (int j = 0; j < i; ++j) s += L[i * K + j] * y[j];
      if (det[i]) {
        if (s < a[i] || s > b[i]) return 1.0;
        y[i] = 0;
        continue;
      }
      double below, above, f;
      if (i == 0) {
        below = below0; above = above0; f = f0;
      } else {
        const double l = L[i * K + i];
        probs((a[i] - s) / l, (b[i] - s) / l, below, above, f);
      }
      if (!(f > 0)) return 1.0;
      log_in += i == 0 ? log0 : log_inside(below, above, f);
      if (i == K - 1) break;
      y[i] = draw(below, above, f, w[i]);
    }
    return -std::expm1(log_in);
  }
};

// The Cholesky factor with the variable reordering of Genz & Bretz (2009, Section 4.1.3): at each step the
// coordinate with the smallest expected interval probability, given the expected values of the earlier ones.
Separated separate(const std::vector<double>& lower, const std::vector<double>& upper, const std::vector<double>& R0,
                   int K) {
  Separated sp;
  sp.K = K; sp.a = lower; sp.b = upper; sp.L.assign(static_cast<std::size_t>(K) * K, 0.0); sp.det.assign(K, 0);
  sp.y.assign(K, 0.0);
  std::vector<double> R = R0, ey(K, 0.0);
  auto swap_vars = [&](int i, int j) {
    if (i == j) return;
    std::swap(sp.a[i], sp.a[j]); std::swap(sp.b[i], sp.b[j]);
    for (int k = 0; k < K; ++k) std::swap(R[i * K + k], R[j * K + k]);
    for (int k = 0; k < K; ++k) std::swap(R[k * K + i], R[k * K + j]);
    for (int k = 0; k < K; ++k) std::swap(sp.L[i * K + k], sp.L[j * K + k]);
  };
  for (int i = 0; i < K; ++i) {
    int best = -1; double best_p = 0;
    for (int j = i; j < K; ++j) {
      double v = R[j * K + j], s = 0;
      for (int k = 0; k < i; ++k) { v -= sp.L[j * K + k] * sp.L[j * K + k]; s += sp.L[j * K + k] * ey[k]; }
      // a deterministic coordinate goes last
      const double pj = v > kTol ? interval_prob((sp.a[j] - s) / std::sqrt(v), (sp.b[j] - s) / std::sqrt(v)) : 2.0;
      if (best < 0 || pj < best_p) { best = j; best_p = pj; }
    }
    swap_vars(i, best);
    double v = R[i * K + i], s = 0;
    for (int k = 0; k < i; ++k) { v -= sp.L[i * K + k] * sp.L[i * K + k]; s += sp.L[i * K + k] * ey[k]; }
    if (!(v > kTol)) { sp.det[i] = 1; ey[i] = 0; continue; }
    const double l = std::sqrt(v);
    sp.L[i * K + i] = l;
    for (int j = i + 1; j < K; ++j) {
      double r = R[j * K + i];
      for (int k = 0; k < i; ++k) r -= sp.L[j * K + k] * sp.L[i * K + k];
      sp.L[j * K + i] = r / l;
    }
    // the expected value of the coordinate truncated to its interval
    const double lo = (sp.a[i] - s) / l, hi = (sp.b[i] - s) / l, pr = interval_prob(lo, hi);
    auto phi = [](double x) { return std::isfinite(x) ? std::exp(-0.5 * x * x) / std::sqrt(2 * kPi) : 0.0; };
    if (pr > 1e-300) ey[i] = (phi(lo) - phi(hi)) / pr;
    else ey[i] = std::isfinite(lo) ? lo : (std::isfinite(hi) ? hi : 0.0);
  }
  return sp;
}

}  // namespace

MvnResult mvn_outside_qmc(const std::vector<double>& lower, const std::vector<double>& upper, const std::vector<double>& R,
                          int K, double abseps, long max_points) {
  MvnResult res;
  if (K < 1) return res;
  Separated sp = separate(lower, upper, R, K);
  sp.prepare();
  const int D = std::max(1, K - 1);      // cube dimension (the last coordinate needs no point)
  const int kShifts = 12;
  if (D > 5) throw std::invalid_argument("mvn_outside_qmc supports at most 6 dimensions");
  std::vector<double> w(D), wa(D), shift(D);
  std::vector<long long> z(D);
  SplitMix64 rng{0x5EED5EED2026ULL};
  double est = 0, var = -1;
  for (std::size_t it = 0;; ++it) {
    // one lattice rule of n points under each shift: the shift means and their spread
    const long long n = kLattice[it][0];
    z[0] = 1;
    for (int d = 1; d < D; ++d) z[d] = (z[d - 1] * kLattice[it][1]) % n;
    double means[kShifts], sum = 0;
    for (int m = 0; m < kShifts; ++m) {
      for (int d = 0; d < D; ++d) shift[d] = rng.next();
      double tot = 0;
      for (long long j = 0; j < n; ++j) {
        for (int d = 0; d < D; ++d) {
          double x = static_cast<double>((j * z[d]) % n) / static_cast<double>(n) + shift[d];
          if (x >= 1) x -= 1;
          w[d] = std::fabs(2 * x - 1);   // the baker's transformation
          wa[d] = 1 - w[d];              // its antithetic point
        }
        tot += 0.5 * (sp.outside(w.data()) + sp.outside(wa.data()));
      }
      means[m] = tot / static_cast<double>(n);
      sum += means[m];
    }
    res.points += n;
    const double mean = sum / kShifts;
    double ss = 0;
    for (double x : means) ss += (x - mean) * (x - mean);
    const double v = ss / (kShifts * (kShifts - 1.0));   // the variance of the mean over the shifts
    if (var < 0) { est = mean; var = v; }
    else if (var + v > 0) { est += (mean - est) * var / (var + v); var = var * v / (var + v); }
    else est = mean;
    res.error = 3.5 * std::sqrt(var);
    if (res.error < abseps || it + 1 == kRules || res.points + kLattice[it + 1][0] > max_points) break;
  }
  res.p = std::min(1.0, std::max(0.0, est));
  return res;
}

double mvn_outside(const std::vector<double>& lower, const std::vector<double>& upper, const std::vector<double>& R, int K) {
  if (K <= 2) return mvn_outside_quadrature(lower, upper, R, K);
  return mvn_outside_qmc(lower, upper, R, K).p;
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

MaxTResult davies_max(const std::vector<std::vector<double>>& infl, const std::vector<double>& t, const std::vector<double>& df) {
  MaxTResult res;
  const int K0 = static_cast<int>(t.size());
  std::vector<int> use;
  for (int k = 0; k < K0; ++k) if (std::isfinite(t[k]) && std::isfinite(df[k])) use.push_back(k);
  if (use.empty()) { res.p = std::numeric_limits<double>::quiet_NaN(); return res; }
  const std::size_t m = infl[use[0]].size();
  auto cov = [&](int k, int l) { double s = 0; for (std::size_t u = 0; u < m; ++u) s += infl[k][u] * infl[l][u]; return s; };
  double c = 0;
  for (int k : use) {
    double z = -norm_quantile(pt_upper(std::fabs(t[k]), df[k]));
    if (z > c) { c = z; res.best = k; }
  }
  if (use.size() == 1) { res.p = pt_two_sided(t[use[0]], df[use[0]]); return res; }
  double V = 0;
  for (std::size_t i = 1; i < use.size(); ++i) {
    double r = cov(use[i - 1], use[i]) / std::sqrt(cov(use[i - 1], use[i - 1]) * cov(use[i], use[i]));
    if (!std::isfinite(r)) r = 0;
    V += std::acos(std::min(1.0, std::max(-1.0, r)));
  }
  res.p = std::min(1.0, 2 * pnorm_upper(c) + V * std::exp(-c * c / 2) / kPi);
  return res;
}

}  // namespace spicyglm

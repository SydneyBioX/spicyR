#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

// REML fit of the linear mixed model spatialMEM() fits with lmerTest::lmer():
//   y = X beta + b_subject + e,  b ~ N(0, (sigma theta)^2),  e_i ~ N(0, sigma^2 / w_i)
// with one random intercept per subject. Per-subject sums reduce every
// quantity to p x p algebra: with phi = theta^2 and s_g = sum of weights in
// subject g, X'V^{-1}X = X'WX - sum_g phi / (1 + phi s_g) A_g A_g', where
// A_g = sum_{i in g} w_i x_i (Woodbury, per subject).
//
// Also returns the exact Hessian of lme4's REML deviance in (theta, sigma) and
// the derivatives of vcov(beta) = sigma^2 (X'V^{-1}X)^{-1}, which lmerTest
// approximates with numDeriv to compute Satterthwaite degrees of freedom.

namespace {

struct Sums {
  int n, p, G;
  std::vector<double> s, c, A;        // per subject: sum w, sum w y, sum w x (G x p)
  std::vector<double> XtWX, XtWy;     // p x p, p
  double ytWy;
};

// Cholesky-based inverse and log-determinant of a symmetric p x p matrix.
bool invSym(const std::vector<double>& M, int p, std::vector<double>& inv, double& logdet) {
  std::vector<double> L(p * p, 0.0);
  logdet = 0;
  for (int j = 0; j < p; ++j) {
    double d = M[j * p + j];
    for (int k = 0; k < j; ++k) d -= L[j * p + k] * L[j * p + k];
    if (!(d > 0)) return false;
    L[j * p + j] = std::sqrt(d);
    logdet += std::log(d);
    for (int i = j + 1; i < p; ++i) {
      double v = M[i * p + j];
      for (int k = 0; k < j; ++k) v -= L[i * p + k] * L[j * p + k];
      L[i * p + j] = v / L[j * p + j];
    }
  }
  // inverse of L, then inv = L^{-T} L^{-1}
  std::vector<double> Li(p * p, 0.0);
  for (int j = 0; j < p; ++j) {
    Li[j * p + j] = 1 / L[j * p + j];
    for (int i = j + 1; i < p; ++i) {
      double v = 0;
      for (int k = j; k < i; ++k) v -= L[i * p + k] * Li[k * p + j];
      Li[i * p + j] = v / L[i * p + i];
    }
  }
  inv.assign(p * p, 0.0);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j <= i; ++j) {
      double v = 0;
      for (int k = i; k < p; ++k) v += Li[k * p + i] * Li[k * p + j];
      inv[i * p + j] = inv[j * p + i] = v;
    }
  return true;
}

// X'V^{-1}X, X'V^{-1}y and y'V^{-1}y, or their derivatives: each subject's
// correction is weighted by kfun(s_g) (phi/(1+phi s), its first or second
// derivative in phi) and `base` says whether to include the X'WX terms.
template <class F>
void quad(const Sums& S, F kfun, bool base, std::vector<double>& M, std::vector<double>& m, double& yy) {
  const int p = S.p;
  M.assign(p * p, 0.0);
  m.assign(p, 0.0);
  yy = 0;
  if (base) { M = S.XtWX; m = S.XtWy; yy = S.ytWy; }
  for (int g = 0; g < S.G; ++g) {
    const double k = kfun(S.s[g]);
    const double* a = &S.A[g * p];
    for (int i = 0; i < p; ++i) {
      m[i] -= k * a[i] * S.c[g];
      for (int j = 0; j < p; ++j) M[i * p + j] -= k * a[i] * a[j];
    }
    yy -= k * S.c[g] * S.c[g];
  }
}

// Profiled REML criterion (up to a constant) at theta.
double profiled(const Sums& S, double theta, bool& ok) {
  const double phi = theta * theta;
  std::vector<double> M, m, Mi;
  double yy, logdetM;
  quad(S, [phi](double s) { return phi / (1 + phi * s); }, true, M, m, yy);
  ok = invSym(M, S.p, Mi, logdetM);
  if (!ok) return R_PosInf;
  double Q = yy;
  for (int i = 0; i < S.p; ++i)
    for (int j = 0; j < S.p; ++j) Q -= m[i] * Mi[i * S.p + j] * m[j];
  if (!(Q > 0)) { ok = false; return R_PosInf; }
  double ld = 0;
  for (int g = 0; g < S.G; ++g) ld += std::log1p(phi * S.s[g]);
  return ld + logdetM + (S.n - S.p) * std::log(Q);
}

// Brent's minimiser on [a, b] (as in R's optimize()).
template <class F>
double brent(F f, double a, double b, double tol) {
  const double c = 0.5 * (3 - std::sqrt(5.0));
  double v = a + c * (b - a), w = v, x = v, e = 0, d = 0;
  double fx = f(x), fv = fx, fw = fx;
  for (int it = 0; it < 500; ++it) {
    const double xm = 0.5 * (a + b), tol1 = 1e-12 * std::fabs(x) + tol / 3, tol2 = 2 * tol1;
    if (std::fabs(x - xm) <= tol2 - 0.5 * (b - a)) break;
    bool golden = true;
    if (std::fabs(e) > tol1) {
      double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw), pp = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if (q > 0) pp = -pp; else q = -q;
      r = e; e = d;
      if (std::fabs(pp) < std::fabs(0.5 * q * r) && pp > q * (a - x) && pp < q * (b - x)) {
        d = pp / q;
        const double u = x + d;
        if (u - a < tol2 || b - u < tol2) d = x < xm ? tol1 : -tol1;
        golden = false;
      }
    }
    if (golden) { e = (x < xm ? b : a) - x; d = c * e; }
    const double u = x + (std::fabs(d) >= tol1 ? d : (d > 0 ? tol1 : -tol1));
    const double fu = f(u);
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
      else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
    }
  }
  return x;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List lmerRandomIntercept(Rcpp::NumericMatrix X, Rcpp::NumericVector y,
                               Rcpp::NumericVector w, Rcpp::IntegerVector group, int G) {
  Sums S;
  S.n = X.nrow(); S.p = X.ncol(); S.G = G;
  const int n = S.n, p = S.p;
  S.s.assign(G, 0.0); S.c.assign(G, 0.0); S.A.assign(G * p, 0.0);
  S.XtWX.assign(p * p, 0.0); S.XtWy.assign(p, 0.0); S.ytWy = 0;
  for (int i = 0; i < n; ++i) {
    const int g = group[i] - 1;
    const double wi = w[i];
    S.s[g] += wi; S.c[g] += wi * y[i]; S.ytWy += wi * y[i] * y[i];
    for (int a = 0; a < p; ++a) {
      S.A[g * p + a] += wi * X(i, a);
      S.XtWy[a] += wi * X(i, a) * y[i];
      for (int b = 0; b < p; ++b) S.XtWX[a * p + b] += wi * X(i, a) * X(i, b);
    }
  }

  // Minimise over theta >= 0: coarse log-spaced scan (plus theta = 0), then
  // Brent between the neighbours of the best grid point.
  auto F = [&](double t) { bool ok; return profiled(S, t, ok); };
  std::vector<double> grid(1, 0.0);
  for (double e = -4; e <= 4 + 1e-9; e += 0.05) grid.push_back(std::pow(10.0, e));
  int best = 0;
  double fbest = F(grid[0]);
  for (std::size_t k = 1; k < grid.size(); ++k) {
    const double f = F(grid[k]);
    if (f < fbest) { fbest = f; best = static_cast<int>(k); }
  }
  if (!std::isfinite(fbest)) return Rcpp::List::create(Rcpp::Named("ok") = false);
  double theta = 0;
  if (best > 0 || F(grid[1] * 1e-3) < fbest) {
    const double lo = best > 0 ? grid[best - 1] : 0.0;
    const double hi = grid[std::min<std::size_t>(best + 1, grid.size() - 1)];
    theta = brent(F, lo, hi, 1e-10 * std::max(hi, 1e-8));
    if (F(0.0) <= F(theta)) theta = 0;
  }

  // Quantities and exact derivatives at the optimum (phi = theta^2).
  const double phi = theta * theta;
  std::vector<double> M, m, M1, m1, M2, m2, Mi;
  double yy, yy1, yy2, logdetM;
  quad(S, [phi](double s) { return phi / (1 + phi * s); }, true, M, m, yy);
  quad(S, [phi](double s) { const double q = 1 + phi * s; return 1 / (q * q); }, false, M1, m1, yy1);
  quad(S, [phi](double s) { const double q = 1 + phi * s; return -2 * s / (q * q * q); }, false, M2, m2, yy2);
  if (!invSym(M, p, Mi, logdetM)) return Rcpp::List::create(Rcpp::Named("ok") = false);

  auto mv = [p](const std::vector<double>& A, const std::vector<double>& v) {
    std::vector<double> r(p, 0.0);
    for (int i = 0; i < p; ++i) for (int j = 0; j < p; ++j) r[i] += A[i * p + j] * v[j];
    return r;
  };
  auto dot = [p](const std::vector<double>& a, const std::vector<double>& b) {
    double r = 0; for (int i = 0; i < p; ++i) r += a[i] * b[i]; return r;
  };
  auto mm = [p](const std::vector<double>& A, const std::vector<double>& B) {
    std::vector<double> r(p * p, 0.0);
    for (int i = 0; i < p; ++i) for (int k = 0; k < p; ++k) for (int j = 0; j < p; ++j)
      r[i * p + j] += A[i * p + k] * B[k * p + j];
    return r;
  };
  auto trace = [p](const std::vector<double>& A) { double r = 0; for (int i = 0; i < p; ++i) r += A[i * p + i]; return r; };

  const std::vector<double> beta = mv(Mi, m);
  const double Q = yy - dot(m, beta);
  std::vector<double> t1 = m1;  // m1 - M1 beta
  { const std::vector<double> M1b = mv(M1, beta); for (int i = 0; i < p; ++i) t1[i] -= M1b[i]; }
  const std::vector<double> betaPhi = mv(Mi, t1);
  const double Q1 = yy1 - 2 * dot(m1, beta) + dot(beta, mv(M1, beta));
  const double Q2 = yy2 - 2 * dot(m2, beta) - 2 * dot(m1, betaPhi) + dot(beta, mv(M2, beta)) +
                    2 * dot(beta, mv(M1, betaPhi));
  double Dl1 = 0, Dl2 = 0;
  for (int g = 0; g < G; ++g) {
    const double q = 1 + phi * S.s[g];
    Dl1 += S.s[g] / q;
    Dl2 -= S.s[g] * S.s[g] / (q * q);
  }
  const std::vector<double> MiM1 = mm(Mi, M1);
  const double D1 = Dl1 + trace(MiM1);
  const double D2 = Dl2 + trace(mm(Mi, M2)) - trace(mm(MiM1, MiM1));

  const double nu = n - p, sigma2 = Q / nu, sigma = std::sqrt(sigma2);
  // Hessian of lme4's REML deviance
  //   f = D(phi) + Q(phi) / sigma^2 + (n - p) log(2 pi sigma^2)
  // with respect to (theta, sigma), using d/dtheta = 2 theta d/dphi.
  Rcpp::NumericMatrix H(2, 2);
  H(0, 0) = 2 * (D1 + Q1 / sigma2) + 4 * phi * (D2 + Q2 / sigma2);
  H(0, 1) = H(1, 0) = -4 * theta * Q1 / (sigma2 * sigma);
  H(1, 1) = 6 * Q / (sigma2 * sigma2) - 2 * nu / sigma2;

  // vcov(beta) = sigma^2 Mi and its derivatives in theta and sigma.
  Rcpp::NumericMatrix V(p, p), dVtheta(p, p), dVsigma(p, p);
  const std::vector<double> MiM1Mi = mm(MiM1, Mi);
  for (int i = 0; i < p; ++i)
    for (int j = 0; j < p; ++j) {
      V(i, j) = sigma2 * Mi[i * p + j];
      dVtheta(i, j) = -sigma2 * MiM1Mi[i * p + j] * 2 * theta;
      dVsigma(i, j) = 2 * sigma * Mi[i * p + j];
    }
  return Rcpp::List::create(
      Rcpp::Named("ok") = true, Rcpp::Named("beta") = Rcpp::wrap(beta),
      Rcpp::Named("theta") = theta, Rcpp::Named("sigma") = sigma,
      Rcpp::Named("vcov") = V, Rcpp::Named("hessian") = H,
      Rcpp::Named("dvcov_theta") = dVtheta, Rcpp::Named("dvcov_sigma") = dVsigma);
}

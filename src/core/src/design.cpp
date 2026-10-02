// The excess with a general design (Supplementary Methods, Section 1), the availability adjustment
// (Section 2) and their CR2 tests. Ported from the R prototypes excess_design.R / adjust_cov.R.
//
// Per unit k (patient) with whitened design X_k = L_k^{-1} D_k (V_k = L_k L_k'), bread B = (sum X'X)^-1
// and b = B c, the Bell-McCaffrey adjustment A_k = (I - X_k B X_k')^{-1/2} gives the unit's CR2
// contribution s_k = a_k' e_k with a_k = A_k X_k b. The Satterthwaite Gram matrix of the loading
// vectors (Supplementary, "Satterthwaite Degrees of Freedom") is, with t_k = X_k' a_k,
//   Omega = diag(|a_k|^2) - T' B T,
// so tr Omega and tr Omega^2 cost O(m q^2) instead of O(images x m^2).
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

#include "spicyglm/stats.hpp"

namespace spicyglm {

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct Units {
  std::vector<std::vector<int>> rows;  // image rows of each unit present, in first-appearance order
};

Units group_units(const std::vector<int>& unit) {
  Units u;
  std::vector<int> seen;
  for (std::size_t i = 0; i < unit.size(); ++i) {
    auto it = std::find(seen.begin(), seen.end(), unit[i]);
    if (it == seen.end()) { seen.push_back(unit[i]); u.rows.push_back({static_cast<int>(i)}); }
    else u.rows[it - seen.begin()].push_back(static_cast<int>(i));
  }
  return u;
}

struct Fit {
  bool ok = false;
  VectorXd theta;
  MatrixXd B;                           // bread
  std::vector<MatrixXd> X;              // whitened design per unit
  std::vector<VectorXd> y;              // whitened response per unit
};

Fit gls(const ImageRows& d, const Units& U, const MatrixXd& D, double tau2) {
  Fit f;
  const int q = static_cast<int>(D.cols());
  MatrixXd F = MatrixXd::Zero(q, q);
  VectorXd g = VectorXd::Zero(q);
  for (const auto& r : U.rows) {
    const int J = static_cast<int>(r.size());
    MatrixXd V(J, J); VectorXd yy(J); MatrixXd Dk(J, q);
    for (int a = 0; a < J; ++a) {
      yy(a) = d.O[r[a]] - d.E[r[a]];
      Dk.row(a) = D.row(r[a]);
      for (int b = 0; b < J; ++b) V(a, b) = tau2 * d.n[r[a]] * d.n[r[b]] + (a == b ? d.v[r[a]] : 0.0);
    }
    Eigen::LLT<MatrixXd> llt(V);
    MatrixXd Xk = llt.matrixL().solve(Dk);
    VectorXd yk = llt.matrixL().solve(yy);
    F += Xk.transpose() * Xk; g += Xk.transpose() * yk;
    f.X.push_back(Xk); f.y.push_back(yk);
  }
  Eigen::ColPivHouseholderQR<MatrixXd> qr(F);
  if (qr.rank() < q) return f;
  f.B = F.inverse();
  f.theta = f.B * g;
  f.ok = true;
  return f;
}

// Pearson statistic of the rank-one frailty under the design (Supplementary Methods, eq. 2).
double pearson(const ImageRows& d, const Units& U, const MatrixXd& D, double tau2) {
  Fit f = gls(d, U, D, tau2);
  if (!f.ok) return std::numeric_limits<double>::quiet_NaN();
  double x = 0;
  for (const auto& r : U.rows) {
    double s = 0, J = 0;
    for (int i : r) {
      double e = d.O[i] - d.E[i] - D.row(i).dot(f.theta);
      s += d.n[i] * e / d.v[i]; J += d.n[i] * d.n[i] / d.v[i];
    }
    x += s * s / (J * (1 + tau2 * J));
  }
  return x;
}

double brent_root(const std::function<double(double)>& f, double a, double b) {
  // bisection safeguarded secant (Brent) to 1e-12 on [a, b]; f(a) > 0 > f(b)
  double fa = f(a), fb = f(b);
  for (int it = 0; it < 300 && b - a > 1e-12 * std::max(1.0, b); ++it) {
    double m = 0.5 * (a + b), s = b - fb * (b - a) / (fb - fa);
    double x = (s > a && s < b && std::fabs(s - m) < 0.5 * (b - a)) ? s : m;
    double fx = f(x);
    if (fx == 0) return x;
    if ((fx > 0) == (fa > 0)) { a = x; fa = fx; } else { b = x; fb = fx; }
    if (it % 3 == 2) {  // keep bisecting so the bracket shrinks
      double mm = 0.5 * (a + b), fm = f(mm);
      if ((fm > 0) == (fa > 0)) { a = mm; fa = fm; } else { b = mm; fb = fm; }
    }
  }
  return 0.5 * (a + b);
}

double design_tau2(const ImageRows& d, const Units& U, const MatrixXd& D) {
  double df = static_cast<double>(U.rows.size()) - D.cols();
  double q0 = pearson(d, U, D, 0);
  if (df < 1 || !std::isfinite(q0) || q0 <= df) return 0;
  double up = 1e-4;
  while (pearson(d, U, D, up) > df && up < 1e6) up *= 4;
  return brent_root([&](double t) { return pearson(d, U, D, t) - df; }, 0, up);
}

}  // namespace

std::vector<DesignResult> design_tests(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                                       const std::vector<double>& Z, int q, const std::vector<double>& contrasts,
                                       int k, double tau2, bool hartung_knapp) {
  std::vector<DesignResult> out(k);
  const int N = static_cast<int>(rows.O.size());
  MatrixXd D(N, q);
  for (int i = 0; i < N; ++i) for (int j = 0; j < q; ++j) D(i, j) = rows.n[i] * Z[static_cast<std::size_t>(i) * q + j];
  Units U = group_units(unit);
  if (static_cast<int>(U.rows.size()) <= q) { for (auto& r : out) r.reason = "too_few_units"; return out; }
  if (tau2 < 0) tau2 = design_tau2(rows, U, D);
  Fit f = gls(rows, U, D, tau2);
  if (!f.ok) { for (auto& r : out) r.reason = "design_not_full_rank"; return out; }
  // per-unit pieces shared by the contrasts: the adjustment A_k and the whitened residuals
  std::vector<MatrixXd> A(U.rows.size());
  std::vector<VectorXd> E(U.rows.size());
  for (std::size_t u = 0; u < U.rows.size(); ++u) {
    const MatrixXd& Xk = f.X[u];
    const int J = static_cast<int>(Xk.rows());
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(MatrixXd::Identity(J, J) - Xk * f.B * Xk.transpose());
    VectorXd ev = es.eigenvalues().cwiseMax(1e-12).cwiseSqrt().cwiseInverse();
    A[u] = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
    E[u] = f.y[u] - Xk * f.theta;
  }
  double X_hk = 0;
  for (int c = 0; c < k; ++c) {
    DesignResult& res = out[c];
    VectorXd cv = Eigen::Map<const VectorXd>(contrasts.data() + static_cast<std::size_t>(c) * q, q);
    VectorXd b = f.B * cv;
    double V = 0, sum_a2 = 0, sum_a4 = 0, sum_a2_tBt = 0, sum_tBt = 0;
    MatrixXd S = MatrixXd::Zero(q, q);
    res.influence.assign(n_units, 0.0);
    for (std::size_t u = 0; u < U.rows.size(); ++u) {
      VectorXd a = A[u] * (f.X[u] * b);
      double s = a.dot(E[u]);
      V += s * s;
      res.influence[unit[U.rows[u][0]]] = s;
      VectorXd t = f.X[u].transpose() * a;
      double a2 = a.squaredNorm(), tBt = t.dot(f.B * t);
      sum_a2 += a2; sum_a4 += a2 * a2; sum_a2_tBt += a2 * tBt; sum_tBt += tBt;
      S += t * t.transpose();
    }
    MatrixXd BS = f.B * S;
    double trO = sum_a2 - sum_tBt, trO2 = sum_a4 - 2 * sum_a2_tBt + (BS * BS).trace();
    res.ok = true; res.tau2 = tau2;
    res.theta.assign(f.theta.data(), f.theta.data() + q);
    res.estimate = cv.dot(f.theta); res.df = trO * trO / trO2;
    const double nu = static_cast<double>(U.rows.size()) - q;
    if (hartung_knapp && nu >= 1) {
      // the model-based variance, scaled by the Pearson statistic and floored at CR2, on m - q df
      if (c == 0) X_hk = pearson(rows, U, D, tau2);
      V = std::max(V, cv.dot(f.B * cv) * std::max(1.0, X_hk / nu)); res.df = nu;
    }
    res.se = std::sqrt(V);
    res.p = pt_two_sided(res.estimate / res.se, res.df);
  }
  return out;
}

DesignResult design_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                         const std::vector<double>& Z, int q, const std::vector<double>& contrast, double tau2) {
  return design_tests(rows, unit, n_units, Z, q, contrast, 1, tau2)[0];
}

DesignResult availability_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                               int n_units, const std::vector<double>& x, double tau2) {
  const std::size_t N = rows.O.size();
  double mx = 0; for (double v : x) mx += v; mx /= static_cast<double>(N);
  std::vector<double> Z(3 * N);
  for (std::size_t i = 0; i < N; ++i) { Z[3 * i] = group[i] == 0; Z[3 * i + 1] = group[i] == 1; Z[3 * i + 2] = x[i] - mx; }
  return design_test(rows, unit, n_units, Z, 3, {-1, 1, 0}, tau2);
}

}  // namespace spicyglm

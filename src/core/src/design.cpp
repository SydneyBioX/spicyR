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
  // unit codes are 0-based: look each one up in a table rather than searching the units seen so far
  int max_unit = -1;
  for (int c : unit) max_unit = std::max(max_unit, c);
  std::vector<int> pos_of(static_cast<std::size_t>(max_unit + 1), -1);
  for (std::size_t i = 0; i < unit.size(); ++i) {
    int& q = pos_of[unit[i]];
    if (q < 0) { q = static_cast<int>(u.rows.size()); u.rows.push_back({}); }
    u.rows[q].push_back(static_cast<int>(i));
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

// Per unit: the CR2 adjustment A_k = (I - X_k B X_k')^{-1/2} and the whitened residuals.
void unit_adjustments(const Fit& f, std::vector<MatrixXd>& A, std::vector<VectorXd>& E) {
  A.assign(f.X.size(), MatrixXd()); E.assign(f.X.size(), VectorXd());
  for (std::size_t u = 0; u < f.X.size(); ++u) {
    const MatrixXd& Xk = f.X[u];
    const int J = static_cast<int>(Xk.rows());
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(MatrixXd::Identity(J, J) - Xk * f.B * Xk.transpose());
    VectorXd ev = es.eigenvalues().cwiseMax(1e-12).cwiseSqrt().cwiseInverse();
    A[u] = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();
    E[u] = f.y[u] - Xk * f.theta;
  }
}

bool no_variation(const Fit& f, const std::vector<VectorXd>& E, int N) {
  // No variation: the design fits the response exactly (e.g. the same fraction in every image, as at a radius where
  // no `to` cell has a `from` neighbour, or every one has), so the residuals are rounding error and every test is 0 / 0
  // computed from noise. Fitting such a response in double precision leaves a residual of the order of n eps |y| (the
  // n rows are summed into F and g, then theta is subtracted back out of each), so a weighted residual sum of squares
  // at most (64 eps n)^2 times the weighted sum of squares of y about zero is taken as none; the factor 64 is room for
  // the q-dimensional solve. It is scale-free, and far below any response with real variation (for n = 1000 it asks
  // the residuals to be under 1.4e-11 of |y|).
  double rss = 0, yss = 0;
  for (std::size_t u = 0; u < E.size(); ++u) { rss += E[u].squaredNorm(); yss += f.y[u].squaredNorm(); }
  const double tol = 64 * std::numeric_limits<double>::epsilon() * N;
  return rss <= tol * tol * yss;
}

}  // namespace

std::vector<DesignResult> design_tests(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                                       const std::vector<double>& Z, int q, const std::vector<double>& contrasts,
                                       int k, double tau2, bool hartung_knapp, int n_groups) {
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
  std::vector<MatrixXd> A;
  std::vector<VectorXd> E;
  unit_adjustments(f, A, E);
  if (no_variation(f, E, N)) { for (auto& r : out) r.reason = "no_variation"; return out; }
  // the pooled dispersion: each condition's CR2 / model variance ratio of its coefficient, pooled by units
  const bool pooled = n_groups > 0 && !hartung_knapp;
  double phi = 1;
  if (pooled) {
    double num = 0, den = 0;
    for (int g = 0; g < n_groups; ++g) {
      const VectorXd bg = f.B.col(g);
      double Vg = 0, mg = 0;
      for (std::size_t u = 0; u < U.rows.size(); ++u) {
        const double s = (A[u] * (f.X[u] * bg)).dot(E[u]);
        Vg += s * s;
        bool in_g = false;
        for (int i : U.rows[u]) if (Z[static_cast<std::size_t>(i) * q + g] != 0) { in_g = true; break; }
        if (in_g) ++mg;
      }
      num += mg * Vg / f.B(g, g); den += mg;
    }
    phi = num / den;
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
    } else if (pooled) {
      V = phi * cv.dot(f.B * cv);
    }
    res.se = std::sqrt(V);
    res.p = pt_two_sided(res.estimate / res.se, res.df);
  }
  return out;
}

JointResult design_joint_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<double>& Z, int q,
                              const std::vector<double>& contrasts, int k, double tau2, bool hartung_knapp) {
  JointResult out;
  out.df1 = k;
  const int N = static_cast<int>(rows.O.size());
  MatrixXd D(N, q);
  for (int i = 0; i < N; ++i) for (int j = 0; j < q; ++j) D(i, j) = rows.n[i] * Z[static_cast<std::size_t>(i) * q + j];
  Units U = group_units(unit);
  const int m = static_cast<int>(U.rows.size());
  if (m <= q) { out.reason = "too_few_units"; return out; }
  if (tau2 < 0) tau2 = design_tau2(rows, U, D);
  Fit f = gls(rows, U, D, tau2);
  if (!f.ok) { out.reason = "design_not_full_rank"; return out; }
  std::vector<MatrixXd> A;
  std::vector<VectorXd> E;
  unit_adjustments(f, A, E);
  if (no_variation(f, E, N)) { out.reason = "no_variation"; return out; }
  MatrixXd C(k, q);
  for (int c = 0; c < k; ++c) for (int j = 0; j < q; ++j) C(c, j) = contrasts[static_cast<std::size_t>(c) * q + j];
  const VectorXd est = C * f.theta;
  const MatrixXd CBC = C * f.B * C.transpose();   // the model covariance of C theta: the expectation of CR2's
  out.tau2 = tau2;
  // CR2 loadings P_u = C B X_u' A_u (k x J_u), their projections T_u = P_u X_u, and the CR2 covariance
  std::vector<MatrixXd> P(m), T(m);
  MatrixXd V = MatrixXd::Zero(k, k);
  for (int u = 0; u < m; ++u) {
    P[u] = C * f.B * f.X[u].transpose() * A[u];
    T[u] = P[u] * f.X[u];
    const VectorXd s = P[u] * E[u];
    V += s * s.transpose();
  }
  Eigen::LDLT<MatrixXd> ldlt(V);
  if (ldlt.info() != Eigen::Success || ldlt.vectorD().minCoeff() <= 0) { out.reason = "singular_covariance"; return out; }
  const double wald_cr2 = est.dot(ldlt.solve(est));
  if (hartung_knapp) {
    // floored at CR2 as design_tests floors each contrast: the smaller of the two Wald statistics (with one
    // contrast, the t test of design_tests squared)
    const double nu = static_cast<double>(m - q);
    if (nu < 1) { out.reason = "too_few_units"; return out; }
    const double s = std::max(1.0, pearson(rows, U, D, tau2) / nu);
    out.wald = std::min(est.dot((s * CBC).ldlt().solve(est)), wald_cr2);
    out.statistic = out.wald / k; out.df2 = nu;
  } else {
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(CBC);
    if (es.eigenvalues().minCoeff() <= 0) { out.reason = "design_not_full_rank"; return out; }
    const MatrixXd W = es.eigenvectors() * es.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() *
                       es.eigenvectors().transpose();
    // normalised by W = (C B C')^-1/2, the covariance of the loadings of units u and v under the working model is
    // Omega_uv = [u == v] P_u P_u' - X_uv, X_uv = T_u B T_v' (normalised T); the entries of V have summed variance
    // sum_uv f(Omega_uv), f(O) = tr(O)^2 + tr(O O). Over all u, v, f(X_uv) sums to sum_u a_u' (G + H) a_u with
    // a_u = vec(T_u B) (index (i, a) = i q + a), G[(i,a),(j,b)] = sum_v T_v[i,a] T_v[j,b] (the tr^2 terms) and
    // H[(i,a),(j,b)] = sum_v T_v[j,a] T_v[i,b] (the tr(O O) terms), in O(m (k q)^2); then u == v is corrected to
    // f(P_u P_u' - X_uu).
    const int kq = k * q;
    std::vector<MatrixXd> TB(m), Tn(m);
    MatrixXd G = MatrixXd::Zero(kq, kq), H = MatrixXd::Zero(kq, kq);
    for (int u = 0; u < m; ++u) {
      Tn[u] = W * T[u]; TB[u] = Tn[u] * f.B;
      for (int i = 0; i < k; ++i)
        for (int a = 0; a < q; ++a)
          for (int j = 0; j < k; ++j)
            for (int b = 0; b < q; ++b) {
              G(i * q + a, j * q + b) += Tn[u](i, a) * Tn[u](j, b);
              H(i * q + a, j * q + b) += Tn[u](j, a) * Tn[u](i, b);
            }
    }
    const MatrixXd GH = G + H;
    auto fO = [](const MatrixXd& O) { const double tr = O.trace(); return tr * tr + (O * O).trace(); };
    double total = 0;
    for (int u = 0; u < m; ++u) {
      VectorXd au(kq);
      for (int i = 0; i < k; ++i) for (int b = 0; b < q; ++b) au(i * q + b) = TB[u](i, b);
      total += au.dot(GH * au);
      const MatrixXd X = TB[u] * Tn[u].transpose();
      const MatrixXd Pn = W * P[u];
      total += fO(Pn * Pn.transpose() - X) - fO(X);
    }
    const double eta = k * (k + 1.0) / total;
    if (!std::isfinite(eta)) { out.reason = "too_few_df"; out.p = std::numeric_limits<double>::quiet_NaN(); return out; }
    out.wald = wald_cr2;
    out.df2 = eta - k + 1;
    out.statistic = out.df2 / (eta * k) * out.wald;
  }
  if (!(out.df2 > 0)) { out.reason = "too_few_df"; out.p = std::numeric_limits<double>::quiet_NaN(); return out; }
  out.p = pf_upper(out.statistic, out.df1, out.df2);
  out.ok = true;
  return out;
}

DesignResult design_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                         const std::vector<double>& Z, int q, const std::vector<double>& contrast, double tau2) {
  return design_tests(rows, unit, n_units, Z, q, contrast, 1, tau2)[0];
}

DesignResult availability_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                               int n_units, const std::vector<double>& x, double tau2) {
  // rows with a missing condition (a negative group) are left out
  ImageRows r;
  std::vector<int> u;
  std::vector<double> xk;
  std::vector<int> gk;
  for (std::size_t i = 0; i < rows.O.size(); ++i) {
    if (group[i] < 0) continue;
    r.image.push_back(rows.image[i]); r.O.push_back(rows.O[i]); r.E.push_back(rows.E[i]); r.n.push_back(rows.n[i]);
    r.v.push_back(rows.v[i]); u.push_back(unit[i]); xk.push_back(x[i]); gk.push_back(group[i]);
  }
  const std::size_t N = r.O.size();
  double mx = 0; for (double v : xk) mx += v; mx /= static_cast<double>(N);
  std::vector<double> Z(3 * N);
  for (std::size_t i = 0; i < N; ++i) { Z[3 * i] = gk[i] == 0; Z[3 * i + 1] = gk[i] == 1; Z[3 * i + 2] = xk[i] - mx; }
  return design_test(r, u, n_units, Z, 3, {-1, 1, 0}, tau2);
}

}  // namespace spicyglm

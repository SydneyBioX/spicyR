// Cox proportional hazards (Efron ties) and the spicyR Cell survival tests (new_methods.pdf,
// Section 3): the score test (the excess regressed on the null-model martingale residuals) and the
// Cox model on the shrunken per-patient excess. Ported from the R prototype survival.R.
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

#include "spicyglm/stats.hpp"

namespace spicyglm {

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;

struct CoxPieces { double loglik = 0; VectorXd grad; MatrixXd info; bool ok = true; };

// Efron log partial likelihood, score and information at beta. `order` sorts subjects by time
// (ascending); tied times are handled as one block.
CoxPieces cox_pieces(const std::vector<double>& time, const std::vector<int>& event, const MatrixXd& X,
                     const VectorXd& beta, const std::vector<int>& order) {
  const int n = static_cast<int>(time.size()), p = static_cast<int>(X.cols());
  CoxPieces c; c.grad = VectorXd::Zero(p); c.info = MatrixXd::Zero(p, p);
  VectorXd eta = p ? VectorXd(X * beta) : VectorXd::Zero(n);
  // risk-set sums accumulated from the latest time backwards
  double S0 = 0; VectorXd S1 = VectorXd::Zero(p); MatrixXd S2 = MatrixXd::Zero(p, p);
  int i = n - 1;
  while (i >= 0) {
    int j = i; double t = time[order[i]];
    while (j >= 0 && time[order[j]] == t) --j;          // block (j, i] has time t
    double D0 = 0; VectorXd D1 = VectorXd::Zero(p); MatrixXd D2 = MatrixXd::Zero(p, p); int d = 0;
    VectorXd sx = VectorXd::Zero(p);
    for (int k = i; k > j; --k) {
      int s = order[k]; double w = std::exp(eta(s));
      S0 += w; if (p) { S1 += w * X.row(s).transpose(); S2 += w * X.row(s).transpose() * X.row(s); }
      if (event[s]) { ++d; D0 += w; if (p) { D1 += w * X.row(s).transpose(); D2 += w * X.row(s).transpose() * X.row(s); sx += X.row(s).transpose(); } c.loglik += eta(s); }
    }
    for (int l = 0; l < d; ++l) {
      double f = static_cast<double>(l) / d, den = S0 - f * D0;
      if (!(den > 0)) { c.ok = false; return c; }
      c.loglik -= std::log(den);
      if (p) {
        VectorXd m1 = (S1 - f * D1) / den;
        c.grad -= m1;
        c.info += (S2 - f * D2) / den - m1 * m1.transpose();
      }
    }
    if (p) c.grad += sx;
    i = j;
  }
  return c;
}

std::vector<double> martingale(const std::vector<double>& time, const std::vector<int>& event, const MatrixXd& X,
                               const VectorXd& beta, const std::vector<int>& order) {
  // M_i = d_i - exp(eta_i) H_i, with Efron's hazard increments: at a time with d events, a subject at
  // risk that does not fail gets sum_l 1 / den_l, one that fails gets sum_l (1 - l / d) / den_l.
  const int n = static_cast<int>(time.size()), p = static_cast<int>(X.cols());
  VectorXd eta = p ? VectorXd(X * beta) : VectorXd::Zero(n);
  // risk sums at each block, computed backwards, then hazards accumulated forwards
  struct Block { int lo, hi; double h_cens, h_event; };
  std::vector<Block> blocks;
  double S0 = 0; int i = n - 1;
  while (i >= 0) {
    int j = i; double t = time[order[i]];
    while (j >= 0 && time[order[j]] == t) --j;
    double D0 = 0; int d = 0;
    for (int k = i; k > j; --k) { int s = order[k]; double w = std::exp(eta(s)); S0 += w; if (event[s]) { ++d; D0 += w; } }
    double hc = 0, he = 0;
    for (int l = 0; l < d; ++l) { double f = static_cast<double>(l) / d, den = S0 - f * D0; hc += 1 / den; he += (1 - f) / den; }
    blocks.push_back({j + 1, i, hc, he});
    i = j;
  }
  std::reverse(blocks.begin(), blocks.end());
  std::vector<double> M(n);
  double H = 0;  // cumulative hazard of the earlier blocks (every subject of a later block was at risk)
  for (const auto& b : blocks) {
    for (int k = b.lo; k <= b.hi; ++k) {
      int s = order[k];
      double Hs = H + (event[s] ? b.h_event : b.h_cens);
      M[s] = event[s] - std::exp(eta(s)) * Hs;
    }
    H += b.h_cens;
  }
  return M;
}

}  // namespace

CoxResult cox_fit(const std::vector<double>& time, const std::vector<int>& event, const std::vector<double>& Xv, int p) {
  CoxResult res;
  const int n = static_cast<int>(time.size());
  MatrixXd X(n, p);
  for (int i = 0; i < n; ++i) for (int j = 0; j < p; ++j) X(i, j) = Xv[static_cast<std::size_t>(i) * p + j];
  std::vector<int> order(n); std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return time[a] < time[b]; });
  VectorXd beta = VectorXd::Zero(p);
  CoxPieces c = cox_pieces(time, event, X, beta, order);
  if (!c.ok) return res;
  bool converged = p == 0;
  for (int it = 0; it < 30 && p; ++it) {
    Eigen::LDLT<MatrixXd> ldlt(c.info);
    if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().minCoeff() > 0)) return res;
    VectorXd step = ldlt.solve(c.grad);
    VectorXd nb = beta + step;
    CoxPieces cn = cox_pieces(time, event, X, nb, order);
    int halve = 0;
    while ((!cn.ok || cn.loglik < c.loglik - 1e-12) && halve < 30) { step /= 2; nb = beta + step; cn = cox_pieces(time, event, X, nb, order); ++halve; }
    if (!cn.ok) return res;
    double change = std::fabs(cn.loglik - c.loglik);
    beta = nb; c = cn;
    if (change <= 1e-9 * (std::fabs(c.loglik) + 1)) { converged = true; break; }
  }
  if (!converged || !beta.allFinite() || (p > 0 && beta.cwiseAbs().maxCoeff() > 30)) return res;  // as coxph's infinite-beta warning
  res.ok = true; res.loglik = c.loglik;
  res.beta.assign(beta.data(), beta.data() + p);
  if (p) {
    MatrixXd Vb = c.info.inverse();
    for (int j = 0; j < p; ++j) {
      double se = std::sqrt(Vb(j, j));
      res.se.push_back(se); res.p.push_back(2 * pnorm_upper(std::fabs(beta(j) / se)));
    }
  }
  res.martingale = martingale(time, event, X, beta, order);
  return res;
}

SurvivalResult survival_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                             const std::vector<double>& M, const std::vector<double>& time,
                             const std::vector<int>& event) {
  SurvivalResult res;
  const std::size_t N = rows.O.size();
  // units with an outcome
  std::vector<int> keep;
  for (std::size_t i = 0; i < N; ++i) if (std::isfinite(M[unit[i]])) keep.push_back(static_cast<int>(i));
  ImageRows r;
  std::vector<int> u;
  for (int i : keep) { r.image.push_back(rows.image[i]); r.O.push_back(rows.O[i]); r.E.push_back(rows.E[i]); r.n.push_back(rows.n[i]); r.v.push_back(rows.v[i]); u.push_back(unit[i]); }
  std::vector<int> present(u); std::sort(present.begin(), present.end()); present.erase(std::unique(present.begin(), present.end()), present.end());
  if (present.size() < 10) { res.reason = "fewer_than_10_patients"; return res; }
  // score test: delta = theta0 + theta1 M_u, tau2 re-estimated under the design
  std::vector<double> Z(2 * r.O.size());
  for (std::size_t i = 0; i < r.O.size(); ++i) { Z[2 * i] = 1; Z[2 * i + 1] = M[u[i]]; }
  DesignResult sc = design_test(r, u, n_units, Z, 2, {0, 1}, -1);
  if (sc.ok) { res.ok = true; res.score_coef = sc.estimate; res.score_se = sc.se; res.score_df = sc.df; res.score_p = sc.p; }
  else { res.reason = sc.reason; return res; }
  // shrunken per-unit excess under the intercept-only frailty model, in a Cox model
  std::vector<double> Z1(r.O.size(), 1.0);
  DesignResult mu = design_test(r, u, n_units, Z1, 1, {1}, -1);
  res.tau2 = mu.tau2;
  std::vector<double> s(n_units, 0), J(n_units, 0);
  for (std::size_t i = 0; i < r.O.size(); ++i) { s[u[i]] += r.n[i] * (r.O[i] - r.E[i]) / r.v[i]; J[u[i]] += r.n[i] * r.n[i] / r.v[i]; }
  // With tau2 near 0 every patient is shrunk to the mean and the shrunken excess carries no information:
  // report no hazard ratio (new_methods.pdf, Section 3, check 4).
  double wmax = 0;
  for (int k : present) wmax = std::max(wmax, mu.ok ? mu.tau2 / (mu.tau2 + 1 / J[k]) : 1.0);
  if (wmax < 0.01) { res.hr_sd = res.log_hr_sd = res.hr_se = res.log_hr_unit = std::numeric_limits<double>::quiet_NaN();
                     res.hr_p = std::numeric_limits<double>::quiet_NaN(); return res; }
  std::vector<double> x, tt; std::vector<int> ee;
  for (int k : present) {
    double raw = s[k] / J[k];
    x.push_back(mu.ok ? mu.theta[0] + mu.tau2 / (mu.tau2 + 1 / J[k]) * (raw - mu.theta[0]) : raw);
    tt.push_back(time[k]); ee.push_back(event[k]);
  }
  double mean = std::accumulate(x.begin(), x.end(), 0.0) / x.size(), ss = 0;
  for (double v : x) ss += (v - mean) * (v - mean);
  double sd = std::sqrt(ss / (x.size() - 1));
  if (!(sd > 1e-12 * (std::fabs(mean) + 1))) { res.hr_sd = res.log_hr_sd = res.hr_se = res.log_hr_unit = res.hr_p = std::numeric_limits<double>::quiet_NaN(); return res; }
  for (double& v : x) v = (v - mean) / sd;
  CoxResult cx = cox_fit(tt, ee, x, 1);
  if (!cx.ok) { res.hr_sd = res.log_hr_sd = res.hr_se = res.log_hr_unit = res.hr_p = std::numeric_limits<double>::quiet_NaN(); return res; }
  res.log_hr_sd = cx.beta[0]; res.hr_sd = std::exp(cx.beta[0]); res.hr_se = cx.se[0]; res.hr_p = cx.p[0];
  res.log_hr_unit = cx.beta[0] / sd;
  return res;
}

}  // namespace spicyglm

// Cox proportional hazards (Efron ties) and the spicyR Cell survival tests (Supplementary Methods,
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

namespace {

// tau2 of a survival design (the rows' columns Z, row-major, q of them; patient_level[j]: column j is constant within
// units): Paule-Mandel under the design when every column is patient-level (then it is the between-patient PM of
// option (iv), DECISIONS 4 Nov 2026), else held at the Paule-Mandel tau2 of its patient-level columns (-1: estimate
// under the design).
// TODO(ext-AB merge): when a column varies within units, call adjusted_design() (ext-AB, design.cpp) instead of this
// hold: for the score test on (1, M_i, X) with q0 = 2, and for the hazard ratio's shrinkage model on (1, X) with
// q0 = 1 (no M_i), passing the rows before the image random effect so omega2 is B2 (within patients, after X); it adds
// the patient-mean column of each image-level column and gives tau2 by option (iv). Then run design_test on its rows,
// Z and q at its tau2 (the contrast padded with 0 for the mean columns), and use its q for the too_few_units rule
// (already counted that way below). Without image-level columns (patient-level covariates only, or none) this
// function is already exact and adjusted_design() does nothing different.
double survival_tau2(const ImageRows& r, const std::vector<int>& u, int n_units, const std::vector<double>& Z, int q,
                     const std::vector<char>& patient_level) {
  std::vector<int> cols;
  for (int j = 0; j < q; ++j) if (patient_level[j]) cols.push_back(j);
  if (static_cast<int>(cols.size()) == q) return -1;
  const int qp = static_cast<int>(cols.size());
  const std::size_t N = r.O.size();
  std::vector<double> Zp(N * qp), c(qp, 0.0);
  for (std::size_t i = 0; i < N; ++i) for (int k = 0; k < qp; ++k) Zp[i * qp + k] = Z[i * q + cols[k]];
  c[0] = 1;
  DesignResult d = design_test(r, u, n_units, Zp, qp, c, -1);
  return d.ok ? d.tau2 : -1;
}

}  // namespace

SurvivalResult survival_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                             const std::vector<double>& M, const std::vector<double>& time,
                             const std::vector<int>& event, const std::vector<double>& X, int p,
                             const std::vector<double>& W, int pw) {
  SurvivalResult res;
  const double kNaN = std::numeric_limits<double>::quiet_NaN();
  res.log_hr_sd = res.hr_sd = res.hr_se = res.log_hr_unit = kNaN;
  const std::size_t N = rows.O.size();
  // rows of units with an outcome, and with every covariate
  std::vector<int> keep;
  for (std::size_t i = 0; i < N; ++i) {
    if (!std::isfinite(M[unit[i]])) continue;
    bool ok = true;
    for (int j = 0; j < p && ok; ++j) ok = !std::isnan(X[i * p + j]);
    if (ok) keep.push_back(static_cast<int>(i));
  }
  ImageRows r;
  std::vector<int> u;
  std::vector<std::vector<double>> x(p);
  for (int i : keep) {
    r.image.push_back(rows.image[i]); r.O.push_back(rows.O[i]); r.E.push_back(rows.E[i]); r.n.push_back(rows.n[i]);
    r.v.push_back(rows.v[i]); u.push_back(unit[i]);
    for (int j = 0; j < p; ++j) x[j].push_back(X[static_cast<std::size_t>(i) * p + j]);
  }
  const std::size_t Nk = r.O.size();
  // the covariates centred over the rows kept; which of them vary within units
  std::vector<char> varies(p, 0);
  for (int j = 0; j < p; ++j) {
    double mx = 0; for (double v : x[j]) mx += v; if (Nk) mx /= static_cast<double>(Nk);
    for (double& v : x[j]) v -= mx;
    std::vector<double> at(n_units); std::vector<char> has(n_units, 0);
    for (std::size_t k = 0; k < Nk && !varies[j]; ++k) {
      if (!has[u[k]]) { has[u[k]] = 1; at[u[k]] = x[j][k]; } else if (x[j][k] != at[u[k]]) varies[j] = 1;
    }
  }
  // the informative units (an outcome and at least one row) and their deaths
  std::vector<int> present(u); std::sort(present.begin(), present.end()); present.erase(std::unique(present.begin(), present.end()), present.end());
  res.n_units = static_cast<int>(present.size());
  for (int k : present) res.n_deaths += event[k] != 0;
  res.q = 2 + p;
  for (int j = 0; j < p; ++j) res.q += varies[j];
  // enough information (DECISIONS 7 Nov 2026): a death and a survivor, and at least q + 2 informative units
  if (res.n_deaths == 0) { res.reason = "no_deaths"; return res; }
  if (res.n_deaths == res.n_units) { res.reason = "no_survivors"; return res; }
  if (res.n_units - res.q < 2) { res.reason = "too_few_units"; return res; }
  // score test: delta = theta0 + theta1 M_u + X beta, CR2 on Satterthwaite df
  const int q = 2 + p;
  std::vector<double> Z(Nk * q), c(q, 0.0);
  std::vector<char> pl(q, 1);
  for (std::size_t i = 0; i < Nk; ++i) {
    Z[i * q] = 1; Z[i * q + 1] = M[u[i]];
    for (int j = 0; j < p; ++j) Z[i * q + 2 + j] = x[j][i];
  }
  for (int j = 0; j < p; ++j) pl[2 + j] = !varies[j];
  c[1] = 1;
  DesignResult sc = design_test(r, u, n_units, Z, q, c, survival_tau2(r, u, n_units, Z, q, pl));
  if (!sc.ok) { res.reason = sc.reason; return res; }
  if (too_few_df(sc.df)) { res.reason = "too_few_df"; return res; }
  res.ok = true; res.score_coef = sc.estimate; res.score_se = sc.se; res.score_df = sc.df; res.score_p = sc.p;
  // the hazard ratio: the shrunken per-unit excess at the average covariates, from the frailty model on (1, X)
  // (without M, so the outcome does not enter its own predictor), in a Cox model with the null model's covariates W
  const int qh = 1 + p;
  std::vector<double> Zh(Nk * qh), ch(qh, 0.0);
  std::vector<char> plh(qh, 1);
  for (std::size_t i = 0; i < Nk; ++i) { Zh[i * qh] = 1; for (int j = 0; j < p; ++j) Zh[i * qh + 1 + j] = x[j][i]; }
  for (int j = 0; j < p; ++j) plh[1 + j] = !varies[j];
  ch[0] = 1;
  DesignResult mu = design_test(r, u, n_units, Zh, qh, ch, survival_tau2(r, u, n_units, Zh, qh, plh));
  std::vector<double> beta(p, 0.0);
  if (mu.ok) for (int j = 0; j < p; ++j) beta[j] = mu.theta[1 + j];
  res.tau2 = mu.ok ? mu.tau2 : kNaN;
  std::vector<double> s(n_units, 0), J(n_units, 0);
  for (std::size_t i = 0; i < Nk; ++i) {
    double fx = 0; for (int j = 0; j < p; ++j) fx += beta[j] * x[j][i];
    s[u[i]] += r.n[i] * (r.O[i] - r.E[i] - r.n[i] * fx) / r.v[i]; J[u[i]] += r.n[i] * r.n[i] / r.v[i]; }
  // With tau2 near 0 every patient is shrunk to the mean and the shrunken excess carries no information:
  // report no hazard ratio (Supplementary Methods, Section 3, check 4).
  double wmax = 0;
  for (int k : present) wmax = std::max(wmax, mu.ok ? mu.tau2 / (mu.tau2 + 1 / J[k]) : 1.0);
  if (wmax < 0.01) return res;
  std::vector<double> hx, tt; std::vector<int> ee;
  for (int k : present) {
    double raw = s[k] / J[k];
    hx.push_back(mu.ok ? mu.theta[0] + mu.tau2 / (mu.tau2 + 1 / J[k]) * (raw - mu.theta[0]) : raw);
    tt.push_back(time[k]); ee.push_back(event[k]);
  }
  double mean = std::accumulate(hx.begin(), hx.end(), 0.0) / hx.size(), ss = 0;
  for (double v : hx) ss += (v - mean) * (v - mean);
  double sd = std::sqrt(ss / (hx.size() - 1));
  if (!(sd > 1e-12 * (std::fabs(mean) + 1))) return res;
  // the Cox model: the standardised shrunken excess, then the covariates of the null model
  const int ph = 1 + pw;
  std::vector<double> Xc(present.size() * ph);
  for (std::size_t k = 0; k < present.size(); ++k) {
    Xc[k * ph] = (hx[k] - mean) / sd;
    for (int j = 0; j < pw; ++j) Xc[k * ph + 1 + j] = W[static_cast<std::size_t>(present[k]) * pw + j];
  }
  CoxResult cx = cox_fit(tt, ee, Xc, ph);
  if (!cx.ok) return res;
  res.log_hr_sd = cx.beta[0]; res.hr_sd = std::exp(cx.beta[0]); res.hr_se = cx.se[0];
  res.log_hr_unit = cx.beta[0] / sd;
  return res;
}

}  // namespace spicyglm

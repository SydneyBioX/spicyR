// The excess with a frailty GEE (Supplementary Part I): image rows, the label-clustering factor,
// Paule-Mandel, the closed-form CR2 variance with Satterthwaite df, and the Hartung-Knapp option.
// Ported from the R front end (R/excess.R, R/frailty.R of spicyR 1.99), which the plain-R reference
// checks to 1e-9.
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

#include "spicyglm/stats.hpp"

namespace spicyglm {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

double median_of(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Brent's root finder on [a, b] with f(a), f(b) of opposite sign (as R's uniroot / zeroin).
double brent(const std::function<double(double)>& f, double a, double b, double tol) {
  double fa = f(a), fb = f(b);
  if (fa == 0) return a;
  if (fb == 0) return b;
  double c = a, fc = fa, d = b - a, e = d;
  for (int it = 0; it < 1000; ++it) {
    if ((fb > 0) == (fc > 0)) { c = a; fc = fa; d = b - a; e = d; }
    if (std::fabs(fc) < std::fabs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    double tol1 = 2 * std::numeric_limits<double>::epsilon() * std::fabs(b) + 0.5 * tol, xm = 0.5 * (c - b);
    if (std::fabs(xm) <= tol1 || fb == 0) return b;
    if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
      double s = fb / fa, p, q;
      if (a == c) { p = 2 * xm * s; q = 1 - s; }
      else { double qq = fa / fc, r = fb / fc; p = s * (2 * xm * qq * (qq - r) - (b - a) * (r - 1)); q = (qq - 1) * (r - 1) * (s - 1); }
      if (p > 0) q = -q; else p = -p;
      if (2 * p < std::min(3 * xm * q - std::fabs(tol1 * q), std::fabs(e * q))) { e = d; d = p / q; }
      else { d = xm; e = d; }
    } else { d = xm; e = d; }
    a = b; fa = fb;
    b += std::fabs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
    fb = f(b);
  }
  return b;
}

// One group's fit at fixed tau2 (R: excess_fit_group). Units are the distinct codes among `idx`.
struct GroupFit {
  double beta = 0, B = 0;
  std::vector<int> units;              // unit codes of the group, in first-appearance order
  std::vector<double> r, J, a, s;      // per unit: score at beta, information, frailty weight, raw score
};

GroupFit fit_group(const ImageRows& d, const std::vector<int>& idx, const std::vector<int>& unit, double tau2) {
  GroupFit f;
  std::vector<int> pos(1, 0);
  std::vector<int> slot;  // per idx entry: its unit's position
  {
    std::vector<int> seen;
    for (int i : idx) {
      auto it = std::find(seen.begin(), seen.end(), unit[i]);
      if (it == seen.end()) { slot.push_back(static_cast<int>(seen.size())); seen.push_back(unit[i]); }
      else slot.push_back(static_cast<int>(it - seen.begin()));
    }
    f.units = seen;
  }
  std::size_t m = f.units.size();
  f.J.assign(m, 0); f.s.assign(m, 0); f.r.assign(m, 0); f.a.assign(m, 0);
  for (std::size_t k = 0; k < idx.size(); ++k) {
    int i = idx[k];
    f.J[slot[k]] += d.n[i] * d.n[i] / d.v[i];
    f.s[slot[k]] += d.n[i] * (d.O[i] - d.E[i]) / d.v[i];
  }
  double num = 0;
  for (std::size_t u = 0; u < m; ++u) { f.a[u] = 1 / (1 + tau2 * f.J[u]); num += f.a[u] * f.s[u]; f.B += f.a[u] * f.J[u]; }
  f.beta = num / f.B;
  for (std::size_t u = 0; u < m; ++u) f.r[u] = f.s[u] - f.beta * f.J[u];
  return f;
}

double pearson_group(const GroupFit& f) {
  double x = 0;
  for (std::size_t u = 0; u < f.units.size(); ++u) x += f.a[u] * f.r[u] * f.r[u] / f.J[u];
  return x;
}

// Closed-form CR2 of one group's beta and the two traces of its Satterthwaite df
// (Supplementary, Proposition "excess-cr2" and Corollary "linear-cost form").
struct GroupCR2 { double V = 0, EV = 0, trsq = 0; };
GroupCR2 cr2_group(const GroupFit& f) {
  GroupCR2 c;
  double sum_t1 = 0;
  for (std::size_t u = 0; u < f.units.size(); ++u) {
    double h = f.a[u] * f.J[u] / f.B;
    c.V += f.a[u] * f.a[u] / (f.B * f.B) * f.r[u] * f.r[u] / (1 - h);
    c.trsq += h * h * (1 - 2 * h) / ((1 - h) * (1 - h));
    sum_t1 += h * h / (1 - h);
  }
  c.EV = 1 / f.B;
  c.trsq = (c.trsq + sum_t1 * sum_t1) / (f.B * f.B);
  return c;
}

void split_groups(const std::vector<int>& group, std::vector<int>& g0, std::vector<int>& g1) {
  for (std::size_t i = 0; i < group.size(); ++i) (group[i] == 0 ? g0 : g1).push_back(static_cast<int>(i));
}

int distinct(const std::vector<int>& idx, const std::vector<int>& unit) {
  std::vector<int> u; for (int i : idx) u.push_back(unit[i]);
  std::sort(u.begin(), u.end()); return static_cast<int>(std::unique(u.begin(), u.end()) - u.begin());
}

}  // namespace

ImageRows excess_image_rows(const std::vector<double>& totals, const std::vector<double>& out_sq_totals,
                            const std::vector<double>& counts, int n_types, int n_images, int from, int to,
                            bool knn, const std::vector<double>& psi) {
  // totals[(img * T + b) * T + a] = sum over the image's b cells of their a neighbours (R: tot[a, b, img])
  const int T = n_types;
  const bool self = from == to;
  auto tot = [&](int a, int b, int img) { return totals[(static_cast<std::size_t>(img) * T + b) * T + a]; };
  auto sq = [&](int a, int b, int img) { return out_sq_totals[(static_cast<std::size_t>(img) * T + b) * T + a]; };
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    double nA = counts[static_cast<std::size_t>(img) * T + from], nB = counts[static_cast<std::size_t>(img) * T + to], N = 0;
    for (int t = 0; t < T; ++t) N += counts[static_cast<std::size_t>(img) * T + t];
    double O = tot(from, to, img), L = 0, Q = 0;
    for (int c = 0; c < T; ++c) {
      if (!self && c == from) continue;
      L += tot(from, c, img);
      Q += sq(c, from, img);
    }
    if (self && !knn) { O -= nA; L -= nA; }   // the radius totals count each cell as its own neighbour
    double M = self ? N : N - nA;
    double p = self ? std::max(nA - 1, 0.0) / std::max(M - 1, 1.0) : nB / std::max(M, 1.0);
    double v = p * (1 - p) * M / std::max(M - 1, 1.0) * std::max(Q - L * L / std::max(M, 1.0), 0.0);
    if (!psi.empty()) v *= psi[static_cast<std::size_t>(to) * n_images + img];
    if (nB > 0 && v > 0) {
      out.image.push_back(img); out.O.push_back(O); out.E.push_back(p * L); out.n.push_back(nB); out.v.push_back(v);
    }
  }
  return out;
}

std::vector<double> label_clustering_factor(const Dataset& data, const std::vector<int>& from,
                                            const std::vector<int>& to, const std::vector<double>& counts,
                                            int n_types, int n_images, bool knn, double h) {
  const int T = n_types, design = knn ? 4 : 3;
  // O of each pair from the REF's side, for the five-pair sparsity rule. For a self-pair on the
  // radius graph the totals count each cell as its own neighbour; subtract it (Supplementary rule;
  // spicyR <= 1.99.0 did not, LOG 2 Oct 2026).
  std::vector<double> totals = data.pair_neighbour_totals(knn);
  std::vector<std::vector<double>> by_ref(T);
  std::vector<std::vector<double>> ratio(from.size(), std::vector<double>(n_images, kNaN));
  for (std::size_t k = 0; k < from.size(); ++k) {
    int f = from[k], t = to[k];
    std::vector<double> V(n_images), Np(n_images), G(n_images);
    if (f == t) {
      std::vector<double> H = data.hac_phi_sums(f, t, design, h);
      for (int i = 0; i < n_images; ++i) { V[i] = H[4 * i]; Np[i] = H[4 * i + 2]; G[i] = H[4 * i + 3]; }
    } else {
      if (by_ref[f].empty()) by_ref[f] = data.hac_phi_sums_ref(f, design, h);
      const std::vector<double>& H = by_ref[f];
      for (int i = 0; i < n_images; ++i) {
        std::size_t o = static_cast<std::size_t>(i) * (T + 3);
        V[i] = H[o + t]; Np[i] = H[o + T + 1]; G[i] = H[o + T + 2];
      }
    }
    for (int i = 0; i < n_images; ++i) {
      double nB = counts[static_cast<std::size_t>(i) * T + t];
      double pr = f == t ? (nB - 1) / std::max(Np[i] - 1, 1.0) : nB / Np[i];
      double r = V[i] / (pr * (1 - pr) * G[i]);
      double O = totals[(static_cast<std::size_t>(i) * T + f) * T + t];
      if (f == t && !knn) O -= counts[static_cast<std::size_t>(i) * T + f];
      if (O < 5 || !std::isfinite(r) || r <= 0) r = kNaN;
      ratio[k][i] = r;
    }
  }
  std::vector<double> out(static_cast<std::size_t>(T) * n_images, 1.0);
  for (int t = 0; t < T; ++t) {
    bool any = false;
    for (std::size_t k = 0; k < to.size(); ++k) if (to[k] == t) any = true;
    if (!any) continue;
    for (int i = 0; i < n_images; ++i) {
      std::vector<double> z;
      for (std::size_t k = 0; k < to.size(); ++k) if (to[k] == t && std::isfinite(ratio[k][i])) z.push_back(ratio[k][i]);
      out[static_cast<std::size_t>(t) * n_images + i] = z.empty() ? 1.0 : std::max(1.0, median_of(z));
    }
  }
  return out;
}

ImageRows kontextual_image_rows(const std::vector<double>& sums, const std::vector<double>& counts, int n_types,
                                int n_images, int from, int to, const std::vector<double>& psi) {
  const int T = n_types;
  const bool self = from == to;
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    const double* s = sums.data() + static_cast<std::size_t>(img) * 7;
    const double nT = counts[static_cast<std::size_t>(img) * T + to];
    // Statial's weights divide by lambda_to = lambda_c n_to / n_context: scale O, L by n_context / n_to
    const double f = nT > 0 ? s[6] / nT : 0.0;
    const double O = s[0] * f, L = s[1] * f, Q = s[2] * f * f, M = s[3], D = s[4];
    const double p = self ? std::max(nT - 1, 0.0) / std::max(M - 1, 1.0) : nT / std::max(M, 1.0);
    double v = p * (1 - p) * M / std::max(M - 1, 1.0) * std::max(Q - L * L / std::max(M, 1.0), 0.0);
    if (!psi.empty()) v *= psi[static_cast<std::size_t>(to) * n_images + img];
    if (nT > 0 && D > 0 && v > 0) {
      out.image.push_back(img); out.O.push_back(O); out.E.push_back(p * L); out.n.push_back(D); out.v.push_back(v);
    }
  }
  return out;
}

std::vector<double> kontextual_clustering_factor(const Dataset& data, const std::vector<int>& from,
                                                 const std::vector<int>& to, const std::vector<double>& raw,
                                                 const std::vector<double>& counts, int n_types, int n_images,
                                                 double h) {
  const int T = n_types;
  std::vector<std::vector<double>> by_ref(T);
  std::vector<std::vector<double>> ratio(from.size(), std::vector<double>(n_images, kNaN));
  for (std::size_t k = 0; k < from.size(); ++k) {
    int f = from[k], t = to[k];
    std::vector<double> V(n_images), Np(n_images), G(n_images);
    if (f == t) {
      std::vector<double> H = data.hac_phi_sums(f, t, 0, h);
      for (int i = 0; i < n_images; ++i) { V[i] = H[4 * i]; Np[i] = H[4 * i + 2]; G[i] = H[4 * i + 3]; }
    } else {
      if (by_ref[f].empty()) by_ref[f] = data.hac_phi_sums_ref(f, 0, h);
      const std::vector<double>& H = by_ref[f];
      for (int i = 0; i < n_images; ++i) {
        std::size_t o = static_cast<std::size_t>(i) * (T + 3);
        V[i] = H[o + t]; Np[i] = H[o + T + 1]; G[i] = H[o + T + 2];
      }
    }
    for (int i = 0; i < n_images; ++i) {
      double nB = counts[static_cast<std::size_t>(i) * T + t];
      double pr = f == t ? (nB - 1) / std::max(Np[i] - 1, 1.0) : nB / Np[i];
      double r = V[i] / (pr * (1 - pr) * G[i]);
      if (raw[k * static_cast<std::size_t>(n_images) + i] < 5 || !std::isfinite(r) || r <= 0) r = kNaN;
      ratio[k][i] = r;
    }
  }
  std::vector<double> out(static_cast<std::size_t>(T) * n_images, 1.0);
  for (int t = 0; t < T; ++t) {
    bool any = false;
    for (std::size_t k = 0; k < to.size(); ++k) if (to[k] == t) any = true;
    if (!any) continue;
    for (int i = 0; i < n_images; ++i) {
      std::vector<double> z;
      for (std::size_t k = 0; k < to.size(); ++k) if (to[k] == t && std::isfinite(ratio[k][i])) z.push_back(ratio[k][i]);
      out[static_cast<std::size_t>(t) * n_images + i] = z.empty() ? 1.0 : std::max(1.0, median_of(z));
    }
  }
  return out;
}

double excess_tau2(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group, int) {
  std::vector<int> g0, g1; split_groups(group, g0, g1);
  int m = distinct(g0, unit) + distinct(g1, unit);
  double df = m - 2;
  auto X = [&](double tau2) { return pearson_group(fit_group(rows, g0, unit, tau2)) + pearson_group(fit_group(rows, g1, unit, tau2)); };
  if (df < 1 || X(0) <= df) return 0;
  double up = 1e-4;
  while (X(up) > df && up < 1e6) up *= 4;
  return brent([&](double t) { return X(t) - df; }, 0, up, 1e-12);
}

ExcessResult excess_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                         int n_units, bool frailty, Variance variance) {
  ExcessResult res;
  std::vector<int> g0, g1; split_groups(group, g0, g1);
  int m0 = distinct(g0, unit), m1 = distinct(g1, unit);
  if (m0 < 2 || m1 < 2) { res.reason = "one_patient_per_group"; return res; }
  double tau2 = frailty ? excess_tau2(rows, unit, group, n_units) : 0;
  GroupFit f0 = fit_group(rows, g0, unit, tau2), f1 = fit_group(rows, g1, unit, tau2);
  GroupCR2 c0 = cr2_group(f0), c1 = cr2_group(f1);
  double v = c0.V + c1.V, df = (c0.EV + c1.EV) * (c0.EV + c1.EV) / (c0.trsq + c1.trsq);
  double nu = m0 + m1 - 2;
  if (variance == Variance::HartungKnapp && nu >= 1) {
    double Vm = 1 / f0.B + 1 / f1.B, X = pearson_group(f0) + pearson_group(f1);
    v = std::max(v, Vm * std::max(1.0, X / nu)); df = nu;
  }
  res.ok = true; res.tau2 = tau2;
  res.coef_ref = f0.beta; res.coef_comp = f1.beta; res.difference = f1.beta - f0.beta;
  res.se = std::sqrt(v); res.df = df; res.p = pt_two_sided(res.difference / res.se, df);
  res.influence.assign(n_units, 0); res.unit_summary.assign(n_units, kNaN); res.unit_info.assign(n_units, 0);
  res.image_weight.assign(rows.O.size(), 0.0);
  for (int g = 0; g < 2; ++g) {
    const GroupFit& f = g == 0 ? f0 : f1;
    const std::vector<int>& idx = g == 0 ? g0 : g1;
    for (int i : idx) {
      auto it = std::find(f.units.begin(), f.units.end(), unit[i]);
      double a = f.a[it - f.units.begin()];
      res.image_weight[i] = a * rows.n[i] * rows.n[i] / rows.v[i] / f.B;   // a_u (n^2 / v) / S_g
    }
  }
  for (int g = 0; g < 2; ++g) {
    const GroupFit& f = g == 0 ? f0 : f1;
    for (std::size_t u = 0; u < f.units.size(); ++u) {
      double h = f.a[u] * f.J[u] / f.B;
      res.influence[f.units[u]] = (g == 0 ? -1 : 1) * f.a[u] / f.B * f.r[u] / std::sqrt(1 - h);
      res.unit_summary[f.units[u]] = f.s[u] / f.J[u];
      res.unit_info[f.units[u]] = f.J[u];
    }
  }
  return res;
}

}  // namespace spicyglm

// The excess with a frailty GEE (Supplementary Part I): image rows, the label-clustering factor,
// Paule-Mandel, the closed-form CR2 variance with Satterthwaite df, and the Hartung-Knapp option.
// Ported from the R front end (R/excess.R, R/frailty.R of spicyR 1.99), which the plain-R reference
// checks to 1e-9.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <stdexcept>
#include <thread>

#include <Eigen/Dense>

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
    // unit codes are 0-based: look each one up in a table rather than searching the units seen so far
    int max_unit = -1;
    for (int i : idx) max_unit = std::max(max_unit, unit[i]);
    std::vector<int> pos_of(static_cast<std::size_t>(max_unit + 1), -1), seen;
    for (int i : idx) {
      int& q = pos_of[unit[i]];
      if (q < 0) { q = static_cast<int>(seen.size()); seen.push_back(unit[i]); }
      slot.push_back(q);
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

// The rows of each condition; a row with a missing condition (a negative group: R's NA_integer_, Python's -1) is in
// neither, so it is left out of the test.
void split_groups(const std::vector<int>& group, std::vector<int>& g0, std::vector<int>& g1) {
  for (std::size_t i = 0; i < group.size(); ++i)
    if (group[i] >= 0) (group[i] == 0 ? g0 : g1).push_back(static_cast<int>(i));
}

int distinct(const std::vector<int>& idx, const std::vector<int>& unit) {
  std::vector<int> u; for (int i : idx) u.push_back(unit[i]);
  std::sort(u.begin(), u.end()); return static_cast<int>(std::unique(u.begin(), u.end()) - u.begin());
}

}  // namespace

ImageRows excess_image_rows(const std::vector<double>& totals, const std::vector<double>& out_sq_totals,
                            const std::vector<double>& counts, int n_types, int n_images, int from, int to,
                            bool knn, const std::vector<double>& psi, const std::vector<double>& self_var,
                            const std::vector<int>& image_group) {
  return excess_image_rows(totals.data(), out_sq_totals.data(), counts.data(), n_types, n_images, from, to, knn,
                           psi.empty() ? nullptr : psi.data(), self_var.empty() ? nullptr : self_var.data(),
                           image_group.empty() ? nullptr : image_group.data());
}

ImageRows excess_image_rows(const double* totals, const double* out_sq_totals, const double* counts, int n_types,
                            int n_images, int from, int to, bool knn, const double* psi, const double* self_var,
                            const int* image_group) {
  // totals[(img * T + b) * T + a] = sum over the image's b cells of their a neighbours (R: tot[a, b, img])
  const int T = n_types;
  const bool self = from == to;
  if (self && !self_var) throw std::invalid_argument("excess_image_rows: a self-pair needs self_var (Dataset::self_variances)");
  auto tot = [&](int a, int b, int img) { return totals[(static_cast<std::size_t>(img) * T + b) * T + a]; };
  auto sq = [&](int a, int b, int img) { return out_sq_totals[(static_cast<std::size_t>(img) * T + b) * T + a]; };
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    if (image_group && image_group[img] < 0) continue;   // a missing condition
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
    // cross pair: Proposition 1 (the scores x are fixed); self pair: the exact Var(O - pL) (Section 3.6)
    double v = self ? self_var[static_cast<std::size_t>(img) * T + from]
                    : p * (1 - p) * M / std::max(M - 1, 1.0) * std::max(Q - L * L / std::max(M, 1.0), 0.0);
    if (psi) v *= psi[static_cast<std::size_t>(to) * n_images + img];
    if (nB > 0 && v > 0) {
      out.image.push_back(img); out.O.push_back(O); out.E.push_back(p * L); out.n.push_back(nB); out.v.push_back(v);
    }
  }
  return out;
}

ImageRows weighted_image_rows(const double* sums, int n_images, int n_labels, int to, bool allocation,
                              const int* image_group) {
  const std::size_t width = 3 + 4 * static_cast<std::size_t>(n_labels);
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    if (image_group && image_group[img] < 0) continue;   // a missing condition
    const double* s = sums + static_cast<std::size_t>(img) * width;
    const double* l = s + 3 + 4 * static_cast<std::size_t>(to);
    const double L = s[0], Q = s[1], M = s[2], O = l[0], nB = l[1], B2 = l[2];
    // the variance of the TARGET weights over the candidates (p (1 - p) M / (M - 1) for 0/1 weights), times the
    // spread of the scores
    const double sb = std::max(B2 - nB * nB / std::max(M, 1.0), 0.0) / std::max(M - 1, 1.0);
    const double v = sb * std::max(Q - L * L / std::max(M, 1.0), 0.0);
    const double E = nB / std::max(M, 1.0) * L;
    if (nB > 0 && v > 0) {
      out.image.push_back(img); out.O.push_back(O); out.E.push_back(E); out.n.push_back(nB); out.v.push_back(v);
    }
  }
  if (!allocation) return out;
  double pooled = 0;
  for (std::size_t i = 0; i < out.O.size(); ++i) pooled += out.O[i] - out.E[i];
  const bool avoid = pooled < 0;
  ImageRows kept;
  for (std::size_t i = 0; i < out.O.size(); ++i) {
    double n = avoid ? out.E[i] : out.n[i] - out.E[i];
    if (n >= 1) {   // room guard: an image whose chosen room is below one cell is not used (Section 3.8)
      kept.image.push_back(out.image[i]); kept.O.push_back(out.O[i]); kept.E.push_back(out.E[i]);
      kept.n.push_back(n); kept.v.push_back(out.v[i]);
    }
  }
  return kept;
}

ImageRows allocation_image_rows(const std::vector<double>& any_totals, const std::vector<double>& self_expected,
                                const std::vector<double>& counts, int n_types, int n_images, int from, int to,
                                const std::vector<double>& psi, const std::vector<double>& self_var,
                                const std::vector<int>& image_group) {
  return allocation_image_rows(any_totals.data(), self_expected.data(), counts.data(), n_types, n_images, from, to,
                               psi.empty() ? nullptr : psi.data(), self_var.empty() ? nullptr : self_var.data(),
                               image_group.empty() ? nullptr : image_group.data());
}

ImageRows allocation_image_rows(const double* any_totals, const double* self_expected, const double* counts,
                                int n_types, int n_images, int from, int to, const double* psi,
                                const double* self_var, const int* image_group) {
  // any_totals[(img * T + a) * T + b] = the number of the image's b cells with an a cell among their neighbours
  const int T = n_types;
  const bool self = from == to;
  if (self && !self_var) throw std::invalid_argument("allocation_image_rows: a self-pair needs self_var (Dataset::self_variances)");
  auto any = [&](int a, int b, int img) { return any_totals[(static_cast<std::size_t>(img) * T + a) * T + b]; };
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    if (image_group && image_group[img] < 0) continue;   // a missing condition: no row, and no say in the side
    const double* cnt = counts + static_cast<std::size_t>(img) * T;
    double nA = cnt[from], nB = cnt[to], N = 0;
    for (int t = 0; t < T; ++t) N += cnt[t];
    double O = any(from, to, img), L = 0;
    for (int c = 0; c < T; ++c) if (self || c != from) L += any(from, c, img);
    double M = self ? N : N - nA;
    double p = self ? std::max(nA - 1, 0.0) / std::max(M - 1, 1.0) : nB / std::max(M, 1.0);
    double E = self ? self_expected[static_cast<std::size_t>(img) * T + from] : p * L;
    // cross pair: Proposition 1 with 0/1 scores; self pair: the exact Var(O) (Section 3.6)
    double v = self ? self_var[static_cast<std::size_t>(img) * T + from]
                    : p * (1 - p) * M / std::max(M - 1, 1.0) * std::max(L - L * L / std::max(M, 1.0), 0.0);
    if (psi) v *= psi[static_cast<std::size_t>(to) * n_images + img];
    if (nB > 0 && v > 0) {
      out.image.push_back(img); out.O.push_back(O); out.E.push_back(E); out.n.push_back(nB); out.v.push_back(v);
    }
  }
  // The side, once per pair from the pooled images with a condition (no group labels): attraction (sum O - E >= 0) scales by the
  // TARGET cells without a REF neighbour under the null, n = n_to - E; avoidance by those with one, n = E.
  double pooled = 0;
  for (std::size_t i = 0; i < out.O.size(); ++i) pooled += out.O[i] - out.E[i];
  const bool avoid = pooled < 0;
  ImageRows kept;
  for (std::size_t i = 0; i < out.O.size(); ++i) {
    double n = avoid ? out.E[i] : out.n[i] - out.E[i];
    if (n >= 1) {   // room guard: an image whose chosen room is below one cell is not used (Section 3.8)
      kept.image.push_back(out.image[i]); kept.O.push_back(out.O[i]); kept.E.push_back(out.E[i]);
      kept.n.push_back(n); kept.v.push_back(out.v[i]);
    }
  }
  return kept;
}

std::vector<double> label_clustering_factor(const Dataset& data, const std::vector<int>& from,
                                            const std::vector<int>& to, const std::vector<double>& counts,
                                            int n_types, int n_images, bool knn, double h, bool allocation) {
  const int T = n_types, design = allocation ? (knn ? 7 : 6) : (knn ? 4 : 3);
  // O of each pair from the REF's side, for the five-pair sparsity rule. For a self-pair on the
  // radius graph the totals count each cell as its own neighbour; subtract it (Supplementary rule;
  // spicyR <= 1.99.0 did not, LOG 2 Oct 2026). Allocation: the TARGET cells next to REF (self excluded).
  std::vector<double> totals = allocation ? data.pair_neighbour_any_totals(knn) : data.pair_neighbour_totals(knn);
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
      if (f == t && !knn && !allocation) O -= counts[static_cast<std::size_t>(i) * T + f];
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

ImageRows kontextual_image_rows_weighted(const std::vector<double>& sums, int n_images) {
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    const double* s = sums.data() + static_cast<std::size_t>(img) * 10;
    const double nT = s[9];
    const double f = nT > 0 ? s[6] / nT : 0.0;
    const double O = s[0] * f, L = s[1] * f, Q = s[2] * f * f, M = s[3], D = s[4], Sb = s[7], B2 = s[8];
    const double sb = std::max(B2 - Sb * Sb / std::max(M, 1.0), 0.0) / std::max(M - 1, 1.0);
    const double v = sb * std::max(Q - L * L / std::max(M, 1.0), 0.0);
    if (nT > 0 && D > 0 && v > 0) {
      out.image.push_back(img); out.O.push_back(O); out.E.push_back(Sb / std::max(M, 1.0) * L); out.n.push_back(D);
      out.v.push_back(v);
    }
  }
  return out;
}

double kontextual_self_variance(const double* ss, double M, double n) {
  // the chance pi_k that k given candidates are all TARGET cells (n of the M candidates, without replacement)
  auto pk = [&](int k) {
    if (k > n) return 0.0L;
    long double x = 1;
    for (int i = 0; i < k; ++i) x *= (n - i) / (M - i);
    return x;
  };
  const long double p1 = pk(1), p2 = pk(2), p3 = pk(3), p4 = pk(4);
  const long double U = ss[0], U2 = ss[1], A2 = ss[2], B2 = ss[3], C2 = ss[4];
  const long double p = std::max(n - 1, 0.0) / std::max(M - 1, 1.0);
  // E O^2 over pairs of links sharing both cells, one cell (S1) and no cell (S0)
  const long double S1 = A2 - 2 * U2, S0 = U * U - U2 - S1;
  const long double EO2 = p2 * U2 + p3 * S1 + p4 * S0;
  // R = O - p L with L = sum_l I_l r_l: E[O p L] and E[(p L)^2]; E R = 0
  const long double EOL = p3 * U * U + (p2 - p3) * C2;
  const long double EL2 = p1 * B2 + p2 * (U * U - B2);
  return static_cast<double>(EO2 - 2 * p * EOL + p * p * EL2);
}

ImageRows kontextual_image_rows(const std::vector<double>& sums, const std::vector<double>& counts, int n_types,
                                int n_images, int from, int to, const std::vector<double>& psi,
                                const std::vector<double>& self_sums) {
  const int T = n_types;
  const bool self = from == to;
  if (self && self_sums.size() != static_cast<std::size_t>(n_images) * 5)
    throw std::invalid_argument("kontextual_image_rows: a self triple needs Dataset::kontextual_self_sums()");
  ImageRows out;
  for (int img = 0; img < n_images; ++img) {
    const double* s = sums.data() + static_cast<std::size_t>(img) * 7;
    const double nT = counts[static_cast<std::size_t>(img) * T + to];
    // Statial's weights divide by lambda_to = lambda_c n_to / n_context: scale O, L by n_context / n_to
    const double f = nT > 0 ? s[6] / nT : 0.0;
    const double O = s[0] * f, L = s[1] * f, Q = s[2] * f * f, M = s[3], D = s[4];
    const double p = self ? std::max(nT - 1, 0.0) / std::max(M - 1, 1.0) : nT / std::max(M, 1.0);
    double v;
    bool informative;
    if (self) {
      // O is quadratic in the labels: the exact variance of O - p L (DECISIONS 8 Nov 2026), 0 below rounding; an
      // image where every relabelling gives the same O - p L (v <= 1e-10 (f U)^2) carries no information
      const double* ss = self_sums.data() + static_cast<std::size_t>(img) * 5;
      v = std::max(kontextual_self_variance(ss, M, nT), 0.0) * f * f;
      informative = v > 1e-10 * (f * ss[0]) * (f * ss[0]);
    } else {
      v = p * (1 - p) * M / std::max(M - 1, 1.0) * std::max(Q - L * L / std::max(M, 1.0), 0.0);
      informative = v > 0;
    }
    if (!psi.empty()) v *= psi[static_cast<std::size_t>(to) * n_images + img];
    if (nT > 0 && D > 0 && informative && v > 0) {
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

ImageRows image_heterogeneity_rows(const ImageRows& rows, const std::vector<int>& unit, double* tau2_image) {
  const std::size_t N = rows.O.size();
  if (unit.size() != N) throw std::invalid_argument("image_heterogeneity_rows: one unit per row");
  std::vector<double> y(N), s2(N);
  for (std::size_t i = 0; i < N; ++i) { y[i] = (rows.O[i] - rows.E[i]) / rows.n[i]; s2[i] = rows.v[i] / (rows.n[i] * rows.n[i]); }
  // Paule-Mandel within units: tau2 solves Q(tau2) = sum_u (k_u - 1), where
  // Q(t) = sum_u sum_j (y - ybar_u(t))^2 / (s2 + t) and ybar_u(t) is the unit's mean weighted by 1 / (s2 + t);
  // Q decreases in t, so tau2 = 0 when Q(0) <= df. Units with one row do not contribute.
  int n_units = 0;
  for (int u : unit) n_units = std::max(n_units, u + 1);
  std::vector<int> k(n_units, 0);
  for (int u : unit) ++k[u];
  double df = 0;
  for (int u = 0; u < n_units; ++u) if (k[u] > 1) df += k[u] - 1;
  std::vector<double> sw(n_units), swy(n_units);
  auto Q = [&](double t) {
    std::fill(sw.begin(), sw.end(), 0.0); std::fill(swy.begin(), swy.end(), 0.0);
    for (std::size_t i = 0; i < N; ++i) {
      const double w = 1 / (s2[i] + t); sw[unit[i]] += w; swy[unit[i]] += w * y[i];
    }
    double q = 0;
    for (std::size_t i = 0; i < N; ++i) {
      const int u = unit[i];
      if (k[u] > 1) { const double d = y[i] - swy[u] / sw[u]; q += d * d / (s2[i] + t); }
    }
    return q;
  };
  double tau2 = 0;
  if (df > 0 && Q(0) > df) {
    // widen the bracket until Q(up) <= df: Q -> 0 as tau2 grows, and the effect's scale is not fixed (Kontextual's
    // effects are in the thousands), so no cap below the overflow guard
    double up = 1e-4;
    while (Q(up) > df && up < 1e300) up *= 4;
    tau2 = brent([&](double t) { return Q(t) - df; }, 0, up, 0);
  }
  if (tau2_image) *tau2_image = tau2;
  ImageRows out;
  out.image = rows.image; out.O = y; out.E.assign(N, 0.0); out.n.assign(N, 1.0); out.v.resize(N);
  for (std::size_t i = 0; i < N; ++i) out.v[i] = s2[i] + tau2;
  return out;
}

std::vector<int> informative_units(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                                   int n_groups) {
  const std::size_t N = rows.O.size();
  if (unit.size() != N || group.size() != N) throw std::invalid_argument("informative_units: one unit and group per row");
  int n_units = 0;
  for (int u : unit) n_units = std::max(n_units, u + 1);
  std::vector<int> unit_group(n_units, -1);
  std::vector<char> informative(n_units, 0);
  for (std::size_t i = 0; i < N; ++i) {
    if (group[i] < 0 || group[i] >= n_groups) continue;
    unit_group[unit[i]] = group[i];
    const double y = (rows.O[i] - rows.E[i]) / rows.n[i];
    if (std::fabs(y) != 1) informative[unit[i]] = 1;
  }
  std::vector<int> out(n_groups, 0);
  for (int u = 0; u < n_units; ++u) if (unit_group[u] >= 0 && informative[u]) ++out[unit_group[u]];
  return out;
}

std::vector<double> effective_units(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                                    int n_groups, double tau2) {
  const std::size_t N = rows.O.size();
  if (unit.size() != N || group.size() != N) throw std::invalid_argument("effective_units: one unit and group per row");
  std::vector<double> out(n_groups, 0.0);
  for (int g = 0; g < n_groups; ++g) {
    std::vector<int> idx;
    for (std::size_t i = 0; i < N; ++i) if (group[i] == g) idx.push_back(static_cast<int>(i));
    if (idx.empty()) continue;
    GroupFit f = fit_group(rows, idx, unit, tau2);
    double sw = 0, sw2 = 0;
    for (std::size_t u = 0; u < f.units.size(); ++u) { const double w = f.a[u] * f.J[u]; sw += w; sw2 += w * w; }
    out[g] = sw * sw / sw2;
  }
  return out;
}

TestRules test_rules(bool allocation) {
  // Every effect (allocation, count, Kontextual) has the minimum effective units rule. The variance is each
  // condition's own CR2 (standard, Welch-type): the pooled dispersion (spicyR 1.99.27-1.99.47) is off
  // (DECISIONS 29 Oct 2026).
  constexpr bool kCountMinEffectiveUnits = true;
  TestRules r;
  r.pooled_dispersion = false;
  r.min_effective_units = allocation || kCountMinEffectiveUnits;
  return r;
}

double excess_tau2(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group, int) {
  std::vector<int> g0, g1; split_groups(group, g0, g1);
  int m = distinct(g0, unit) + distinct(g1, unit);
  double df = m - 2;
  auto X = [&](double tau2) { return pearson_group(fit_group(rows, g0, unit, tau2)) + pearson_group(fit_group(rows, g1, unit, tau2)); };
  if (df < 1 || X(0) <= df) return 0;
  double up = 1e-4;
  while (X(up) > df && up < 1e300) up *= 4;   // no cap on the effect's scale (Paule-Mandel above)
  return brent([&](double t) { return X(t) - df; }, 0, up, 1e-12);
}

ExcessResult excess_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                         int n_units, bool frailty, Variance variance, bool pooled_dispersion) {
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
  } else if (variance == Variance::CR2 && pooled_dispersion) {
    // each group's CR2 / model variance ratio, pooled by units, times the model variance (Student- rather than
    // Welch-type: the groups share one dispersion); the df stay those of CR2
    const double phi = (m0 * c0.V / c0.EV + m1 * c1.V / c1.EV) / (m0 + m1);
    v = phi * (c0.EV + c1.EV);
  }
  res.ok = true; res.tau2 = tau2;
  res.coef_ref = f0.beta; res.coef_comp = f1.beta; res.difference = f1.beta - f0.beta;
  res.se = std::sqrt(v); res.df = df; res.t = res.difference / res.se;
  if (variance == Variance::CR2 && !pooled_dispersion) {
    // the skewness correction of the per-condition CR2 t (excess_skewness)
    std::vector<int> gu; std::vector<double> Wu, eu;
    for (int g = 0; g < 2; ++g) {
      const GroupFit& f = g == 0 ? f0 : f1;
      for (std::size_t u = 0; u < f.units.size(); ++u) {
        gu.push_back(g); Wu.push_back(f.a[u] * f.J[u]); eu.push_back(f.r[u] / f.J[u]);
      }
    }
    const double S[2] = {f0.B, f1.B};
    res.skewness = excess_skewness(gu, Wu, eu, S);
    if (std::isfinite(res.t)) res.t = skew_corrected_t(res.t, res.skewness);
  }
  res.p = pt_two_sided(res.t, df);
  res.influence.assign(n_units, 0); res.unit_summary.assign(n_units, kNaN); res.unit_info.assign(n_units, 0);
  res.image_weight.assign(rows.O.size(), 0.0);
  res.n_eff.assign(2, 0.0); res.n_units.assign(2, 0.0);
  res.group_weight = {f0.B, f1.B};
  for (int g = 0; g < 2; ++g) {
    const GroupFit& f = g == 0 ? f0 : f1;
    double sw = 0, sw2 = 0;
    for (std::size_t u = 0; u < f.units.size(); ++u) { const double w = f.a[u] * f.J[u]; sw += w; sw2 += w * w; }
    res.n_eff[g] = sw * sw / sw2; res.n_units[g] = static_cast<double>(f.units.size());
  }
  std::vector<int> pos_of(static_cast<std::size_t>(std::max(n_units, 0)), -1);
  for (int g = 0; g < 2; ++g) {
    const GroupFit& f = g == 0 ? f0 : f1;
    const std::vector<int>& idx = g == 0 ? g0 : g1;
    for (std::size_t u = 0; u < f.units.size(); ++u) {
      if (f.units[u] >= static_cast<int>(pos_of.size())) pos_of.resize(f.units[u] + 1, -1);
      pos_of[f.units[u]] = static_cast<int>(u);
    }
    for (int i : idx) {
      double a = f.a[pos_of[unit[i]]];
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

double excess_skewness(const std::vector<int>& g, const std::vector<double>& W, const std::vector<double>& e,
                       const double S[2]) {
  if (!(S[0] > 0 && S[1] > 0)) return 0;
  double z2 = 0, z3 = 0, d2 = 0, d3 = 0, w15[2] = {0, 0};
  int M = 0, n_g[2] = {0, 0};
  std::vector<double> z(W.size(), kNaN);
  for (std::size_t u = 0; u < W.size(); ++u) {
    if ((g[u] != 0 && g[u] != 1) || !(W[u] > 0)) continue;
    const double h = W[u] / S[g[u]];
    z[u] = std::sqrt(W[u]) * e[u];
    z2 += z[u] * z[u]; z3 += z[u] * z[u] * z[u]; d2 += 1 - h; d3 += (1 - h) * (1 - 2 * h); ++M;
    w15[g[u]] += W[u] * std::sqrt(W[u]); ++n_g[g[u]];
  }
  // a condition with fewer than 3 units leaves its third moment without residual df: with two units, sum (1 - h)(1 - 2h)
  // is second order in their weight difference while sum z^3 is first order, so g would be unbounded
  if (n_g[0] < 3 || n_g[1] < 3) return 0;
  if (!(d2 > 0 && d3 > 0 && z2 > 0)) return 0;
  const double s2 = z2 / d2, s = std::sqrt(s2), g1 = (z3 / d3) / (s2 * s);
  // shrunk by its own uncertainty: g1 g1^2 / (g1^2 + V), V the delta-method (influence function) variance of the
  // skewness, mean_u IF_u^2 / n with IF_u = t^3 - 3t - g1 - (3/2) g1 (t^2 - 1), t = z / sigma, n = sum (1 - h) + 1
  // (6 / n for normal residuals; larger when a few units carry the third moment)
  double if2 = 0;
  for (std::size_t u = 0; u < W.size(); ++u) {
    if (std::isnan(z[u])) continue;
    const double t = z[u] / s, IF = t * t * t - 3 * t - g1 - 1.5 * g1 * (t * t - 1);
    if2 += IF * IF;
  }
  const double V = if2 / M / (d2 + 1);
  if (!(V > 0)) return 0;
  const double gamma = g1 * g1 * g1 / (g1 * g1 + V);
  const double k = gamma * (w15[1] / (S[1] * S[1] * S[1]) - w15[0] / (S[0] * S[0] * S[0])) /
                   std::pow(1 / S[0] + 1 / S[1], 1.5);
  return std::isfinite(k) ? k : 0;
}

TestCorrelation test_correlation(const double* y, const double* w, int n_tests, int n_images, const int* image_unit,
                                 const int* image_group, int n_units, const double* group_weight, const int* perms,
                                 int n_perms, int n_threads) {
  if (n_tests < 0 || n_images < 0 || n_units < 0 || n_perms < 0) throw std::invalid_argument("test_correlation: negative size");
  for (int i = 0; i < n_images; ++i)
    if (image_unit[i] < 0 || image_unit[i] >= n_units) throw std::invalid_argument("test_correlation: unit out of range");
  const std::size_t T = static_cast<std::size_t>(n_tests), N = static_cast<std::size_t>(n_images);
  TestCorrelation out;
  out.difference.assign(T, kNaN); out.se.assign(T, kNaN);
  // per test, its units: the unit, its condition, h, ybar, W = h S_g, and psi
  struct Unit { int unit, group; double h, ybar, W, psi; };
  std::vector<std::vector<Unit>> units(T);
  std::vector<double> h(n_units), sy(n_units);
  std::vector<int> ug(n_units);
  for (std::size_t t = 0; t < T; ++t) {
    std::fill(h.begin(), h.end(), 0.0); std::fill(sy.begin(), sy.end(), 0.0); std::fill(ug.begin(), ug.end(), -1);
    double beta[2] = {0, 0};
    bool any[2] = {false, false};
    for (std::size_t i = 0; i < N; ++i) {
      const double wi = w[t * N + i], yi = y[t * N + i];
      const int g = image_group[i];
      if (!(std::isfinite(wi) && std::isfinite(yi) && wi > 0) || g < 0 || g > 1) continue;
      h[image_unit[i]] += wi; sy[image_unit[i]] += wi * yi; ug[image_unit[i]] = g;
      beta[g] += wi * yi; any[g] = true;
    }
    if (!any[0] || !any[1]) continue;
    const double S[2] = {group_weight[2 * t], group_weight[2 * t + 1]};
    double ss = 0;
    for (int u = 0; u < n_units; ++u) {
      if (!(h[u] > 0)) continue;
      const int g = ug[u];
      const double psi = (g == 0 ? -1 : 1) * (sy[u] - h[u] * beta[g]) / std::sqrt(1 - h[u]);
      units[t].push_back({u, g, h[u], sy[u] / h[u], h[u] * S[g], psi});
      ss += psi * psi;
    }
    out.difference[t] = beta[1] - beta[0]; out.se[t] = std::sqrt(ss);
  }
  // R from a tests x columns matrix whose rows are each normalised to unit length (NaN rows: no test)
  auto gram = [&](std::vector<double>& M, std::size_t C) {
    std::vector<char> valid(T, 0);
    for (std::size_t t = 0; t < T; ++t) {
      double s = 0;
      for (std::size_t c = 0; c < C; ++c) s += M[t * C + c] * M[t * C + c];
      valid[t] = s > 0 && std::isfinite(s);
      for (std::size_t c = 0; c < C; ++c) M[t * C + c] = valid[t] ? M[t * C + c] / std::sqrt(s) : 0.0;
    }
    using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    Eigen::Map<const RowMat> X(M.data(), static_cast<Eigen::Index>(T), static_cast<Eigen::Index>(C));
    RowMat G = X * X.transpose();
    out.R.assign(T * T, kNaN);
    for (std::size_t a = 0; a < T; ++a) {
      if (!valid[a]) continue;
      out.R[a * T + a] = 1;
      for (std::size_t b = a + 1; b < T; ++b) if (valid[b]) out.R[a * T + b] = out.R[b * T + a] = G(a, b);
    }
  };
  if (n_perms == 0) {
    std::vector<double> M(T * static_cast<std::size_t>(n_units), 0.0);
    for (std::size_t t = 0; t < T; ++t) for (const Unit& x : units[t]) M[t * n_units + x.unit] = x.psi;
    gram(M, static_cast<std::size_t>(n_units));
    return out;
  }
  // relabellings: each test in closed form (fit_group and cr2_group on the units' summaries at fixed weights)
  const std::size_t B = static_cast<std::size_t>(n_perms);
  std::vector<int> generated;
  if (!perms) {
    // each relabelling permutes the units' conditions over the units with one (Fisher-Yates)
    std::vector<int> cond(n_units, -1), with;
    for (int i = 0; i < n_images; ++i) if (image_group[i] == 0 || image_group[i] == 1) cond[image_unit[i]] = image_group[i];
    for (int u = 0; u < n_units; ++u) if (cond[u] >= 0) with.push_back(u);
    std::vector<int> c(with.size());
    for (std::size_t k = 0; k < with.size(); ++k) c[k] = cond[with[k]];
    generated.assign(B * static_cast<std::size_t>(n_units), -1);
    SplitMix64 rng{0x5EED5EED2026ULL};
    for (std::size_t b = 0; b < B; ++b) {
      for (std::size_t k = c.size(); k > 1; --k) {
        const std::size_t j = std::min(k - 1, static_cast<std::size_t>(rng.next() * static_cast<double>(k)));
        std::swap(c[k - 1], c[j]);
      }
      for (std::size_t k = 0; k < with.size(); ++k) generated[b * n_units + with[k]] = c[k];
    }
    perms = generated.data();
  }
  out.z.assign(T * B, 0.0);
  auto one_test = [&](std::size_t t) {
    const std::vector<Unit>& us = units[t];
    if (us.empty()) return;
    for (std::size_t b = 0; b < B; ++b) {
      const int* pg = perms + b * static_cast<std::size_t>(n_units);
      double Wg[2] = {0, 0}, Wy[2] = {0, 0};
      int m[2] = {0, 0};
      for (const Unit& x : us) {
        const int g = pg[x.unit];
        if (g == 0 || g == 1) { Wg[g] += x.W; Wy[g] += x.W * x.ybar; ++m[g]; }
      }
      if (m[0] < 2 || m[1] < 2) continue;   // as excess_test, not tested ("one_patient_per_group"): z = 0
      const double beta[2] = {Wy[0] / Wg[0], Wy[1] / Wg[1]};
      double V[2] = {0, 0}, tr[2] = {0, 0}, t1[2] = {0, 0};
      for (const Unit& x : us) {
        const int g = pg[x.unit];
        if (g != 0 && g != 1) continue;
        const double hh = x.W / Wg[g], e = x.ybar - beta[g];
        V[g] += hh * hh * e * e / (1 - hh);
        tr[g] += hh * hh * (1 - 2 * hh) / ((1 - hh) * (1 - hh));
        t1[g] += hh * hh / (1 - hh);
      }
      double ev = 0, trsq = 0;
      for (int g = 0; g < 2; ++g) { ev += 1 / Wg[g]; trsq += (tr[g] + t1[g] * t1[g]) / (Wg[g] * Wg[g]); }
      const double D = beta[1] - beta[0], df = ev * ev / trsq;
      // the skewness correction, as excess_test
      std::vector<int> gu; std::vector<double> Wu, eu;
      for (const Unit& x : us) {
        const int g = pg[x.unit];
        if (g != 0 && g != 1) continue;
        gu.push_back(g); Wu.push_back(x.W); eu.push_back(x.ybar - beta[g]);
      }
      const double kappa = excess_skewness(gu, Wu, eu, Wg);
      const double tc = skew_corrected_t(D / std::sqrt(V[0] + V[1]), kappa);
      const double p = pt_two_sided(tc, df);
      const double z = (tc < 0 ? 1 : -1) * norm_quantile(std::max(p, 1e-300) / 2);
      out.z[t * B + b] = std::isfinite(z) ? z : 0.0;
    }
  };
  const int G = std::max(1, std::min(n_threads, n_tests));
  std::atomic<int> next(0);
  std::vector<std::exception_ptr> err(G);
  auto work = [&](int g) {
    try {
      for (int t = next++; t < n_tests; t = next++) one_test(static_cast<std::size_t>(t));
    } catch (...) {
      err[g] = std::current_exception();
      next = n_tests;
    }
  };
  std::vector<std::thread> pool;
  for (int g = 1; g < G; ++g) pool.emplace_back(work, g);
  work(0);
  for (std::thread& th : pool) th.join();
  for (const std::exception_ptr& e : err) if (e) std::rethrow_exception(e);
  // the correlation over the relabellings: each test's z centred
  std::vector<double> M(out.z);
  for (std::size_t t = 0; t < T; ++t) {
    double m = 0;
    for (std::size_t b = 0; b < B; ++b) m += M[t * B + b];
    m /= std::max<std::size_t>(B, 1);
    for (std::size_t b = 0; b < B; ++b) M[t * B + b] -= m;
  }
  gram(M, B);
  return out;
}

}  // namespace spicyglm

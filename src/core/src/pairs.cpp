// Every pair at every radius in one call (stats.hpp, pair_tests): the per-test steps of the front ends (R/cell.R
// .cell_rows, .cell_rows_test, .cell_design and .cell_combine of spicyR 1.99.23) on plain arrays, so that a test never
// goes back to R or Python. The arithmetic is theirs step for step: where R summed in long double (sum, colMeans),
// so does this file, and the results are identical.
// Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <thread>

#include "spicyglm/stats.hpp"

namespace spicyglm {

namespace {

const double kNaN = std::numeric_limits<double>::quiet_NaN();

// a missing condition: R's NA_integer_ or Python's -1
bool missing_group(int g) { return g < 0; }

ImageRows subset_rows(const ImageRows& r, const std::vector<std::size_t>& idx) {
  ImageRows out;
  for (std::size_t i : idx) {
    out.image.push_back(r.image[i]); out.O.push_back(r.O[i]); out.E.push_back(r.E[i]); out.n.push_back(r.n[i]);
    out.v.push_back(r.v[i]);
  }
  return out;
}

// One radius's main and unadjusted tests, as the combination over radii reads them.
struct RadiusTest {
  double t = 0, df = 0, p = 0;
  std::vector<double> influence;
};

struct Combined { int best = -1; double p = kNaN, p_best = kNaN; };

// R/spicy_main.R .cell_combine's comb(): the tests of the radii where the pair was tested (`at`: their radii).
Combined combine_radii(const std::vector<RadiusTest>& tests, const std::vector<int>& at, Combine how) {
  Combined c;
  if (tests.empty()) return c;
  std::vector<double> t, df, pv;
  for (const RadiusTest& x : tests) { t.push_back(x.t); df.push_back(x.df); pv.push_back(x.p); }
  if (how == Combine::MaxT || how == Combine::Davies) {
    std::vector<std::vector<double>> infl;
    for (const RadiusTest& x : tests) infl.push_back(x.influence);
    MaxTResult m = how == Combine::MaxT ? max_t(infl, t, df) : davies_max(infl, t, df);
    c.best = at[m.best]; c.p = m.p; c.p_best = pv[m.best];
  } else {
    int b = -1;   // which.min: the first smallest, NaN left out
    for (std::size_t i = 0; i < pv.size(); ++i) if (!std::isnan(pv[i]) && (b < 0 || pv[i] < pv[b])) b = static_cast<int>(i);
    c.p = cauchy_combine(pv);
    if (b >= 0) { c.best = at[b]; c.p_best = pv[b]; }
  }
  return c;
}

}  // namespace

PairDesign pair_design_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                            int n_units, const double* counts, int n_types, int from, bool adjust, const double* z,
                            int n_z, bool frailty, double tau2, bool hartung_knapp, bool fit) {
  PairDesign out;
  const std::size_t N = rows.O.size();
  std::vector<char> keep(N);
  for (std::size_t i = 0; i < N; ++i) keep[i] = !missing_group(group[i]);
  // the log share of the `from` type, when it is not constant over the rows with a condition (R: var(s[keep]) > 0)
  std::vector<double> share;
  if (adjust) {
    share.resize(N);
    for (std::size_t i = 0; i < N; ++i) {
      const double* c = counts + static_cast<std::size_t>(rows.image[i]) * n_types;
      double total = 0;
      for (int t = 0; t < n_types; ++t) total += c[t];
      share[i] = std::log(std::max(c[from], 0.5) / total);
    }
    bool seen = false; double first = 0;
    for (std::size_t i = 0; i < N; ++i) {
      if (!keep[i]) continue;
      if (!seen) { seen = true; first = share[i]; }
      else if (share[i] != first) { out.abundance = true; break; }
    }
  }
  const int qx = (out.abundance ? 1 : 0) + n_z;
  if (qx == 0) return out;
  if (!fit) {
    out.design = 2; out.tests.assign(1, DesignResult()); out.tests[0].reason = "no_variation";
    return out;
  }
  if (n_z > 0)   // rows with a missing covariate are left out
    for (std::size_t i = 0; i < N; ++i) {
      if (!keep[i]) continue;
      const double* zi = z + static_cast<std::size_t>(rows.image[i]) * n_z;
      for (int j = 0; j < n_z; ++j) if (std::isnan(zi[j])) { keep[i] = 0; break; }
    }
  std::vector<std::size_t> idx;
  for (std::size_t i = 0; i < N; ++i) if (keep[i]) idx.push_back(i);
  const std::size_t Nk = idx.size();
  // the added columns over the kept rows, centred by their means (R: colMeans, summed in long double)
  std::vector<std::vector<double>> X(qx, std::vector<double>(Nk));
  for (std::size_t k = 0; k < Nk; ++k) {
    const std::size_t i = idx[k];
    int j = 0;
    if (out.abundance) X[j++][k] = share[i];
    const double* zi = z ? z + static_cast<std::size_t>(rows.image[i]) * n_z : nullptr;
    for (int c = 0; c < n_z; ++c) X[j++][k] = zi[c];
  }
  for (int j = 0; j < qx; ++j) {
    long double s = 0;
    for (double x : X[j]) s += x;
    s /= static_cast<long double>(Nk);
    const double m = static_cast<double>(s);
    for (double& x : X[j]) x = x - m;
  }
  const int q = 2 + qx;
  std::vector<double> Z(Nk * q);
  std::vector<int> unit_k(Nk);
  for (std::size_t k = 0; k < Nk; ++k) {
    const std::size_t i = idx[k];
    Z[k * q] = group[i] == 0 ? 1.0 : 0.0; Z[k * q + 1] = group[i] == 1 ? 1.0 : 0.0;
    for (int j = 0; j < qx; ++j) Z[k * q + 2 + j] = X[j][k];
    unit_k[k] = unit[i];
  }
  // tau2: re-estimated when every added column is constant within units, else held
  double tau2_d = 0;
  if (frailty) {
    bool patient_level = true;
    std::vector<double> at(n_units);
    std::vector<char> has(n_units);
    for (int j = 0; j < qx && patient_level; ++j) {
      std::fill(has.begin(), has.end(), 0);
      for (std::size_t k = 0; k < Nk; ++k) {
        const int u = unit_k[k];
        if (!has[u]) { has[u] = 1; at[u] = X[j][k]; }
        else if (X[j][k] != at[u]) { patient_level = false; break; }
      }
    }
    tau2_d = patient_level ? -1 : tau2;
  }
  // the contrasts: the difference of the conditions, then each added column
  std::vector<double> C(static_cast<std::size_t>(q - 1) * q, 0.0);
  C[0] = -1; C[1] = 1;
  for (int j = 1; j <= qx; ++j) C[static_cast<std::size_t>(j) * q + 1 + j] = 1;
  out.tests = design_tests(subset_rows(rows, idx), unit_k, n_units, Z, q, C, q - 1, tau2_d, hartung_knapp);
  out.design = out.tests[0].ok ? 1 : 2;
  return out;
}

PairTests pair_tests(const PairTestInput& in) {
  PairTests res;
  const int P = static_cast<int>(in.from.size()), K = static_cast<int>(in.a.size()), nz = in.n_z, ne = 1 + nz;
  const std::size_t PK = static_cast<std::size_t>(P) * K;
  const bool adjusting = in.adjust || nz > 0, several = K > 1 && in.combine != Combine::None;
  res.n_pairs = P; res.n_radii = K; res.n_z = nz;
  res.ok.assign(PK, 0); res.side.assign(PK, 0); res.n_rows.assign(PK, 0); res.tau2_image.assign(PK, kNaN);
  for (auto* v : {&res.coef_ref, &res.coef_comp, &res.difference, &res.se, &res.df, &res.p, &res.tau2, &res.u_difference,
                  &res.u_se, &res.u_df, &res.u_p})
    v->assign(PK, kNaN);
  res.design.assign(PK, 0); res.abundance.assign(PK, 0); res.design_reason.assign(PK, std::string());
  res.effect.assign(PK * ne, kNaN); res.effect_p.assign(PK * ne, kNaN);
  res.rows.assign(P, ImageRows()); res.image_weight.assign(P, std::vector<double>());
  res.best.assign(P, -1); res.comb_p.assign(P, kNaN); res.comb_p_best.assign(P, kNaN); res.comb_u_p.assign(P, kNaN);

  auto one_pair = [&](int pr) {
    const int f = in.from[pr], t = in.to[pr];
    std::vector<RadiusTest> main, unadj;
    std::vector<int> at;
    for (int k = 0; k < K; ++k) {
      const std::size_t o = static_cast<std::size_t>(pr) * K + k;
      ImageRows r = in.allocation
        ? allocation_image_rows(in.a[k], in.b[k], in.counts, in.n_types, in.n_images, f, t, in.psi)
        : excess_image_rows(in.a[k], in.b[k], in.counts, in.n_types, in.n_images, f, t, in.knn, in.psi);
      const std::size_t N = r.O.size();
      std::vector<int> unit(N), group(N);
      for (std::size_t i = 0; i < N; ++i) { unit[i] = in.image_unit[r.image[i]]; group[i] = in.image_group[r.image[i]]; }
      // allocation: the pair's side from all images (R: sum(O - E) < 0, summed in long double)
      if (in.allocation) {
        long double s = 0;
        for (std::size_t i = 0; i < N; ++i) s += r.O[i] - r.E[i];
        res.side[o] = static_cast<double>(s) < 0 ? -1 : 1;
      }
      if (in.image_heterogeneity && N > 0) {
        double tau2_image = 0;
        r = image_heterogeneity_rows(r, unit, &tau2_image);
        res.tau2_image[o] = tau2_image;
      }
      res.n_rows[o] = static_cast<int>(N);
      ExcessResult u = excess_test(r, unit, group, in.n_units, in.frailty, in.variance);
      if (k == in.keep_radius) { res.rows[pr] = r; res.image_weight[pr] = u.image_weight; }
      if (!u.ok) continue;
      res.ok[o] = 1;
      res.u_difference[o] = u.difference; res.u_se[o] = u.se; res.u_df[o] = u.df; res.u_p[o] = u.p;
      res.coef_ref[o] = u.coef_ref; res.coef_comp[o] = u.coef_comp; res.difference[o] = u.difference;
      res.se[o] = u.se; res.df[o] = u.df; res.p[o] = u.p; res.tau2[o] = u.tau2;
      const std::vector<double>* influence = &u.influence;
      // allocation: a condition with no informative unit (every image of every one of its patients at -1 or +1, the
      // bounds of the fraction) has a CR2 variance of exactly 0, so the test's se would come from the other condition
      // alone. The radius is not tested: no design is fitted, p stays NaN (main and unadjusted), and nothing of it
      // enters the combination over radii. (A radius with no variation at all is "no_variation" below.)
      if (in.allocation && u.se > 0) {
        const std::vector<int> inf = informative_units(r, unit, group, 2);
        if (inf[0] == 0 || inf[1] == 0) {
          res.p[o] = kNaN; res.u_p[o] = kNaN;
          if (adjusting) { res.design[o] = 2; res.design_reason[o] = "no_information"; }
          if (several) {
            at.push_back(k);
            main.push_back({kNaN, res.df[o], kNaN, u.influence});
            if (adjusting) unadj.push_back({kNaN, u.df, kNaN, u.influence});
          }
          continue;
        }
      }
      PairDesign d;
      if (adjusting) {
        // no variation (e.g. no `to` cell has a `from` neighbour, or every one has): the unadjusted test is 0 / 0, and
        // an adjusted design would only turn its exact zeros into rounding noise. It is not fitted and the radius is
        // not tested: t and p stay NaN, so nothing of it enters the combination over radii.
        const bool fit = u.se > 0;
        d = pair_design_test(r, unit, group, in.n_units, in.counts, in.n_types, f, in.adjust, in.z, nz, in.frailty,
                             u.tau2, in.variance == Variance::HartungKnapp, fit);
        if (d.design == 2 && !fit) res.p[o] = kNaN;
        res.design[o] = d.design; res.abundance[o] = d.abundance;
        if (d.design == 2) res.design_reason[o] = d.tests[0].reason;
        if (d.design == 1) {
          const DesignResult& x = d.tests[0];
          res.coef_ref[o] = x.theta[0]; res.coef_comp[o] = x.theta[1]; res.difference[o] = x.estimate;
          res.se[o] = x.se; res.df[o] = x.df; res.p[o] = x.p; res.tau2[o] = x.tau2;
          influence = &x.influence;
          int j = 1;
          if (d.abundance) { res.effect[o * ne] = d.tests[j].estimate; res.effect_p[o * ne] = d.tests[j].p; ++j; }
          for (int c = 0; c < nz; ++c, ++j) {
            res.effect[o * ne + 1 + c] = d.tests[j].estimate; res.effect_p[o * ne + 1 + c] = d.tests[j].p;
          }
        }
      }
      if (several) {
        at.push_back(k);
        main.push_back({res.difference[o] / res.se[o], res.df[o], res.p[o], *influence});
        if (adjusting) unadj.push_back({u.difference / u.se, u.df, u.p, u.influence});
      }
    }
    if (several && !at.empty()) {
      Combined c = combine_radii(main, at, in.combine);
      res.best[pr] = c.best; res.comb_p[pr] = c.p; res.comb_p_best[pr] = c.p_best;
      if (adjusting) res.comb_u_p[pr] = combine_radii(unadj, at, in.combine).p;
    }
  };

  const int G = std::max(1, std::min(in.n_threads, P));
  std::atomic<int> next(0);
  std::vector<std::exception_ptr> err(G);
  auto work = [&](int g) {
    try {
      for (int pr = next++; pr < P; pr = next++) one_pair(pr);
    } catch (...) {
      err[g] = std::current_exception();
      next = P;
    }
  };
  std::vector<std::thread> pool;
  for (int g = 1; g < G; ++g) pool.emplace_back(work, g);
  work(0);
  for (std::thread& th : pool) th.join();
  for (const std::exception_ptr& e : err) if (e) std::rethrow_exception(e);
  return res;
}

}  // namespace spicyglm

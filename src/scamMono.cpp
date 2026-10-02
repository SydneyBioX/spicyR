#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

// The scam() fit calcWeights() uses: a Gaussian additive model with
// monotone-decreasing P-spline smooths ("mpd"), smoothing parameters chosen by
// GCV. The coefficients of the smooths enter as exp(beta) (the monotonicity
// reparameterisation), so the fit is non-linear in beta, but with a Gaussian
// identity link it only ever needs X'X, X'y and y'y: each Newton step of
// scam.fit() becomes q x q algebra.
//
// Inner fit (scam.fit(), Gaussian case): minimise
//   ||y - X b(beta)||^2 + beta' S beta,  b_j = exp(beta_j) for j in iv,
// with Newton steps when the Hessian is positive definite and scam's
// Fisher-type step otherwise, halving steps that increase the objective.
// GCV = n * deviance / (n - trA)^2 with trA as scam.fit() computes it.

namespace {

struct Problem {
  int q, n;
  std::vector<double> G, b;  // X'X (q x q), X'y
  double yty;
  std::vector<std::vector<double>> S;  // penalty matrices (q x q)
  std::vector<char> iv;                // exponentiated coefficients
};

bool chol(const std::vector<double>& A, int q, std::vector<double>& L) {
  L.assign(q * q, 0.0);
  for (int j = 0; j < q; ++j) {
    double d = A[j * q + j];
    for (int k = 0; k < j; ++k) d -= L[j * q + k] * L[j * q + k];
    if (!(d > 0)) return false;
    L[j * q + j] = std::sqrt(d);
    for (int i = j + 1; i < q; ++i) {
      double v = A[i * q + j];
      for (int k = 0; k < j; ++k) v -= L[i * q + k] * L[j * q + k];
      L[i * q + j] = v / L[j * q + j];
    }
  }
  return true;
}

// Solve A x = r given the Cholesky factor L of A.
std::vector<double> cholSolve(const std::vector<double>& L, int q, std::vector<double> r) {
  for (int i = 0; i < q; ++i) {
    for (int k = 0; k < i; ++k) r[i] -= L[i * q + k] * r[k];
    r[i] /= L[i * q + i];
  }
  for (int i = q - 1; i >= 0; --i) {
    for (int k = i + 1; k < q; ++k) r[i] -= L[k * q + i] * r[k];
    r[i] /= L[i * q + i];
  }
  return r;
}

struct Fit {
  std::vector<double> beta, bt, r;
  double dev, gcv, trA;
  bool ok;
};

class Model {
 public:
  explicit Model(const Problem& P) : P(P) {}

  // scam.fit() at smoothing parameters sp, warm-started from `start`.
  Fit fit(const std::vector<double>& sp, const std::vector<double>& start) const {
    const int q = P.q;
    std::vector<double> Sl(q * q, 0.0);
    for (std::size_t j = 0; j < P.S.size(); ++j)
      for (int k = 0; k < q * q; ++k) Sl[k] += sp[j] * P.S[j][k];
    Fit f;
    f.beta = start;
    f.ok = false;
    std::vector<double> bt(q), r(q), Cd(q), C1(q), A(q * q), H(q * q), L;
    auto transform = [&](const std::vector<double>& beta) {
      for (int j = 0; j < q; ++j) bt[j] = P.iv[j] ? std::exp(beta[j]) : beta[j];
    };
    // penalised deviance; also fills bt and r = X'(y - X bt)
    auto pdev = [&](const std::vector<double>& beta, double& dev) {
      transform(beta);
      double Gb, bb = 0, pen = 0;
      for (int i = 0; i < q; ++i) {
        Gb = 0;
        for (int k = 0; k < q; ++k) Gb += P.G[i * q + k] * bt[k];
        r[i] = P.b[i] - Gb;
        bb += bt[i] * (P.b[i] + r[i]);  // bt'b + bt'(b - G bt) = 2 bt'b - bt'G bt
        double Sb = 0;
        for (int k = 0; k < q; ++k) Sb += Sl[i * q + k] * beta[k];
        pen += beta[i] * Sb;
      }
      dev = P.yty - bb;
      return dev + pen;
    };
    double dev;
    double old = pdev(f.beta, dev);
    std::vector<double> betaOld = f.beta;
    for (int iter = 0; iter < 200; ++iter) {
      // A = C G C + S (C = diag(Cd)); Hessian H = A - diag(C1 * r)
      for (int j = 0; j < q; ++j) { Cd[j] = P.iv[j] ? bt[j] : 1.0; C1[j] = P.iv[j] ? bt[j] : 0.0; }
      for (int i = 0; i < q; ++i)
        for (int k = 0; k < q; ++k) {
          A[i * q + k] = Cd[i] * P.G[i * q + k] * Cd[k] + Sl[i * q + k];
          H[i * q + k] = A[i * q + k] - (i == k ? C1[i] * r[i] : 0.0);
        }
      std::vector<double> next(q);
      if (chol(H, q, L)) {
        // Newton: beta - H^{-1} grad, grad = -C r + S beta
        std::vector<double> g(q);
        for (int i = 0; i < q; ++i) {
          double Sb = 0;
          for (int k = 0; k < q; ++k) Sb += Sl[i * q + k] * f.beta[k];
          g[i] = -Cd[i] * r[i] + Sb;
        }
        const std::vector<double> step = cholSolve(L, q, g);
        for (int i = 0; i < q; ++i) next[i] = f.beta[i] - step[i];
      } else {
        // scam's step when the Hessian is not positive definite:
        // A^{-1} (C r + C G C beta)
        if (!chol(A, q, L)) return f;
        std::vector<double> rhs(q);
        for (int i = 0; i < q; ++i) {
          double v = Cd[i] * r[i];
          for (int k = 0; k < q; ++k) v += Cd[i] * P.G[i * q + k] * Cd[k] * f.beta[k];
          rhs[i] = v;
        }
        next = cholSolve(L, q, rhs);
      }
      for (double v : next) if (!std::isfinite(v)) return f;
      double pd = pdev(next, dev);
      const double thresh = 10 * (0.1 + std::fabs(old)) * std::sqrt(2.220446e-16);
      for (int h = 0; h < 100 && (!std::isfinite(pd) || pd - old > thresh); ++h) {
        for (int i = 0; i < q; ++i) next[i] = (next[i] + f.beta[i]) / 2;
        pd = pdev(next, dev);
      }
      betaOld = f.beta;
      f.beta = next;
      const bool small = std::fabs(pd - old) / (0.1 + std::fabs(pd)) < 1e-12;
      old = pd;
      if (small) {
        double gmax = 0, bmax = 0;
        for (int i = 0; i < q; ++i) {
          double Sb = 0;
          for (int k = 0; k < q; ++k) Sb += Sl[i * q + k] * f.beta[k];
          const double cd = P.iv[i] ? bt[i] : 1.0;
          gmax = std::max(gmax, std::fabs(-cd * r[i] + Sb));
          bmax = std::max(bmax, std::fabs(f.beta[i] + betaOld[i]) / 2);
        }
        if (gmax <= 1e-10 * std::max(bmax, 1.0)) { f.ok = true; break; }
      }
    }
    if (!f.ok) return f;
    pdev(f.beta, dev);
    f.bt = bt;
    f.dev = dev;
    // trA = tr(H^{-1} C G C) if H is positive definite, else tr(A^{-1} C G C)
    for (int j = 0; j < q; ++j) { Cd[j] = P.iv[j] ? bt[j] : 1.0; C1[j] = P.iv[j] ? bt[j] : 0.0; }
    std::vector<double> CGC(q * q);
    for (int i = 0; i < q; ++i)
      for (int k = 0; k < q; ++k) {
        CGC[i * q + k] = Cd[i] * P.G[i * q + k] * Cd[k];
        A[i * q + k] = CGC[i * q + k] + Sl[i * q + k];
        H[i * q + k] = A[i * q + k] - (i == k ? C1[i] * r[i] : 0.0);
      }
    if (!chol(H, q, L) && !chol(A, q, L)) { f.ok = false; return f; }
    double trA = 0;
    std::vector<double> col(q);
    for (int k = 0; k < q; ++k) {
      for (int i = 0; i < q; ++i) col[i] = CGC[i * q + k];
      trA += cholSolve(L, q, col)[k];
    }
    f.r = r;
    f.trA = trA;
    f.gcv = P.n * dev / ((P.n - trA) * (P.n - trA));
    return f;
  }

 private:
  const Problem& P;
};

}  // namespace

// Fit the model at scam's GCV-optimal smoothing parameters. `start` holds the
// starting coefficients (on the exp scale for iv, as pcls() returns them);
// blockStart/blockLen give each smooth's coefficients (0-based) for
// initial.sp.scam(). Returns the coefficients on the response scale (exp
// applied), the log smoothing parameters and the GCV.
// [[Rcpp::export]]
Rcpp::List scamMonoFit(Rcpp::NumericMatrix XtX, Rcpp::NumericVector Xty, double yty, int n,
                       Rcpp::List S, Rcpp::LogicalVector iv, Rcpp::NumericVector start,
                       Rcpp::IntegerVector blockStart, Rcpp::IntegerVector blockLen) {
  Problem P;
  P.q = XtX.nrow();
  P.n = n;
  const int q = P.q;
  P.G.resize(q * q);
  for (int i = 0; i < q; ++i) for (int k = 0; k < q; ++k) P.G[i * q + k] = XtX(i, k);
  P.b.assign(Xty.begin(), Xty.end());
  P.yty = yty;
  for (R_xlen_t j = 0; j < S.size(); ++j) {
    Rcpp::NumericMatrix Sj = S[j];
    std::vector<double> s(q * q);
    for (int i = 0; i < q; ++i) for (int k = 0; k < q; ++k) s[i * q + k] = Sj(i, k);
    P.S.push_back(s);
  }
  P.iv.resize(q);
  for (int j = 0; j < q; ++j) P.iv[j] = iv[j];
  const int m = static_cast<int>(P.S.size());
  Model model(P);

  std::vector<double> warm(start.begin(), start.end());
  for (int j = 0; j < q; ++j) if (P.iv[j]) warm[j] = std::log(warm[j]);

  // scam()'s starting point: the fit at sp = 0.05 from the penalised least
  // squares start (`start`, computed in R as pcls() does), then
  // initial.sp.scam()'s smoothing parameters from that fit's Hessian.
  Fit f0 = model.fit(std::vector<double>(m, 0.05), warm);
  if (!f0.ok) return Rcpp::List::create(Rcpp::Named("ok") = false);
  warm = f0.beta;
  std::vector<double> rho(m);
  for (int j = 0; j < m; ++j) {
    const int lo = blockStart[j], len = blockLen[j];
    double hn = 0, sn = 0;
    for (int a = lo; a < lo + len; ++a)
      for (int b = lo; b < lo + len; ++b) {
        const double cda = P.iv[a] ? f0.bt[a] : 1.0, cdb = P.iv[b] ? f0.bt[b] : 1.0;
        double h = cda * P.G[a * q + b] * cdb;
        if (a == b && P.iv[a]) h -= f0.bt[a] * f0.r[a];
        hn += h * h;
        sn += P.S[j][a * q + b] * P.S[j][a * q + b];
      }
    rho[j] = std::log(std::sqrt(hn / sn) + 1e-4);
  }

  // GCV and its gradient at rho. Like scam's fn(), each call warm-starts from
  // the previous fit and becomes the next start; the gradient is by central
  // differences from this fit (scam computes it analytically).
  struct Eval { double score, scaleEst; std::vector<double> grad; Fit fit; bool ok; };
  auto fn = [&](const std::vector<double>& r) {
    Eval e;
    std::vector<double> sp(m);
    for (int j = 0; j < m; ++j) sp[j] = std::exp(r[j]);
    e.fit = model.fit(sp, warm);
    e.ok = e.fit.ok;
    if (!e.ok) { e.score = R_PosInf; return e; }
    warm = e.fit.beta;
    e.score = e.fit.gcv;
    e.scaleEst = e.fit.dev / (P.n - e.fit.trA);
    e.grad.assign(m, 0.0);
    const double h = 1e-5;
    for (int j = 0; j < m; ++j) {
      std::vector<double> up = sp, dn = sp;
      up[j] = std::exp(r[j] + h);
      dn[j] = std::exp(r[j] - h);
      const Fit fu = model.fit(up, e.fit.beta), fd = model.fit(dn, e.fit.beta);
      if (!fu.ok || !fd.ok) { e.ok = false; e.score = R_PosInf; return e; }
      e.grad[j] = (fu.gcv - fd.gcv) / (2 * h);
    }
    return e;
  };

  // bfgs_gcv.ubre() with scam.control() defaults (typx = 1, typf = 1).
  const double maxNstep = 5, steptol = 1e-7, gradtol = 1e-6, c1 = 1e-4, c2 = 0.9;
  const int maxHalf = 30, maxStep = 200;
  Eval b = fn(rho);
  if (!b.ok) return Rcpp::List::create(Rcpp::Named("ok") = false);
  double score = b.score;
  std::vector<double> grad = b.grad;
  // initial inverse Hessian from finite differences of the gradient
  std::vector<double> B(m * m, 0.0);
  {
    const double feps = 1e-4;
    std::vector<double> H(m * m);
    for (int j = 0; j < m; ++j) {
      std::vector<double> r2 = rho;
      r2[j] += feps;
      const Eval b2 = fn(r2);
      if (!b2.ok) return Rcpp::List::create(Rcpp::Named("ok") = false);
      for (int i = 0; i < m; ++i) H[i * m + j] = (b2.grad[i] - grad[i]) / feps;
    }
    Rcpp::NumericMatrix Hs(m, m);
    for (int i = 0; i < m; ++i) for (int j = 0; j < m; ++j) Hs(i, j) = (H[i * m + j] + H[j * m + i]) / 2;
    Rcpp::Function eigen("eigen");
    Rcpp::List eh = eigen(Hs, Rcpp::Named("symmetric") = true);
    Rcpp::NumericVector ev = eh["values"];
    Rcpp::NumericMatrix vec = eh["vectors"];
    double mx = 0;
    for (int k = 0; k < m; ++k) { ev[k] = std::fabs(ev[k]); mx = std::max(mx, ev[k]); }
    for (int k = 0; k < m; ++k) if (ev[k] < mx * 1e-4) ev[k] = mx * 1e-4;
    for (int i = 0; i < m; ++i)
      for (int j = 0; j < m; ++j) {
        double v = 0;
        for (int k = 0; k < m; ++k) v += vec(i, k) * vec(j, k) / ev[k];
        B[i * m + j] = v;
      }
  }
  // fn() for B moved the warm start; resume from the fit at rho
  warm = b.fit.beta;
  std::vector<char> unconv(m);
  {
    const double sc = std::fabs(b.scaleEst) + std::fabs(score);
    bool any = false;
    for (int j = 0; j < m; ++j) { unconv[j] = std::fabs(grad[j]) > sc * gradtol; any = any || unconv[j]; }
    if (!any) std::fill(unconv.begin(), unconv.end(), 1);
  }
  int consecmax = 0;
  for (int it = 1; it <= maxStep; ++it) {
    std::vector<double> Nstep(m, 0.0);
    for (int i = 0; i < m; ++i) {
      if (!unconv[i]) continue;
      double v = 0;
      for (int j = 0; j < m; ++j) if (unconv[j]) v -= B[i * m + j] * grad[j];
      Nstep[i] = v;
    }
    double dot = 0;
    for (int i = 0; i < m; ++i) dot += Nstep[i] * grad[i];
    if (dot >= 0) for (int i = 0; i < m; ++i) Nstep[i] = -B[i * m + i] * grad[i];
    double Newtlen = 0;
    for (double v : Nstep) Newtlen += v * v;
    Newtlen = std::sqrt(Newtlen);
    if (Newtlen > maxNstep) {
      for (double& v : Nstep) v *= maxNstep / Newtlen;
      Newtlen = maxNstep;
    }
    bool maxtaken = false, curv = true;
    int retcode = 2;
    double initslope = 0, rhomax = 1, ms = 0;
    for (int i = 0; i < m; ++i) { initslope += Nstep[i] * grad[i]; rhomax = std::max(rhomax, std::fabs(rho[i])); ms = std::max(ms, std::fabs(Nstep[i])); }
    const double rellength = ms / rhomax;
    const double alphaMin = steptol / rellength;
    double alphaMax, alpha;
    if (ms - maxNstep > std::pow(2.220446e-16, 0.9)) { alpha = maxNstep / ms; alphaMax = alpha * 1.05; }
    else { alpha = 1; alphaMax = std::min(2.0, maxNstep / ms); }
    double oldAlpha = 0, oldScore1 = 0, score1 = 0, newslope = 0;
    std::vector<double> rho1(m), grad1(m);
    Eval bb;
    int ii = 0;
    auto at = [&](double a) { for (int i = 0; i < m; ++i) rho1[i] = rho[i] + a * Nstep[i]; };
    auto slope = [&](const std::vector<double>& g) { double v = 0; for (int i = 0; i < m; ++i) v += g[i] * Nstep[i]; return v; };
    auto near1 = [](double a) { return std::fabs(a - 1) < 1.5e-8; };
    for (;;) {
      at(alpha);
      bb = fn(rho1);
      score1 = bb.score;
      if (score1 <= score + c1 * alpha * initslope) {
        grad1 = bb.grad;
        newslope = slope(grad1);
        curv = true;
        if (newslope < c2 * initslope) {
          if (near1(alpha) && Newtlen < maxNstep) {
            for (int kk = 0; kk < 40; ++kk) {
              oldAlpha = alpha; oldScore1 = score1;
              alpha = std::min(2 * alpha, alphaMax);
              at(alpha);
              bb = fn(rho1);
              score1 = bb.score;
              if (score1 <= score + c1 * alpha * initslope) { grad1 = bb.grad; newslope = slope(grad1); }
              if (score1 > score + c1 * alpha * initslope) break;
              if (newslope >= c2 * initslope) break;
              if (alpha >= alphaMax) break;
            }
          }
          if ((!near1(alpha) && alpha < 1) || ((!near1(alpha) && alpha > 1) && score1 > score + c1 * alpha * initslope)) {
            double alphaLo = std::min(alpha, oldAlpha), alphaDiff = std::fabs(oldAlpha - alpha), scLo, scHi;
            if (alpha < oldAlpha) { scLo = score1; scHi = oldScore1; } else { scLo = oldScore1; scHi = score1; }
            for (int kk = 0; kk < 40; ++kk) {
              double incr = -newslope * alphaDiff * alphaDiff / (2 * (scHi - (scLo + newslope * alphaDiff)));
              if (incr < 0.2 * alphaDiff) incr = 0.2 * alphaDiff;
              alpha = alphaLo + incr;
              at(alpha);
              bb = fn(rho1);
              score1 = bb.score;
              if (score1 > score + c1 * alpha * initslope) { alphaDiff = incr; scHi = score1; }
              else {
                grad1 = bb.grad;
                newslope = slope(grad1);
                if (newslope < c2 * initslope) { alphaLo = alpha; alphaDiff -= incr; scLo = score1; }
              }
              if (std::fabs(newslope) <= -c2 * initslope) break;
              if (alphaDiff < alphaMin) break;
            }
            if (newslope < c2 * initslope) {
              curv = false;
              score1 = scLo;
              at(alphaLo);
              bb = fn(rho1);
            }
          }
        }
        retcode = 0;
        if (newslope < c2 * initslope) curv = false;
        if (alpha * Newtlen > 0.99 * maxNstep) maxtaken = true;
      } else if (alpha < alphaMin) {
        retcode = 1;
        rho1 = rho;
        bb = fn(rho1);
      } else {
        ++ii;
        double tmp;
        if (alpha == 1) {
          tmp = -initslope / (score1 - score - initslope) / 2;
        } else {
          const double a11 = 1 / (alpha * alpha), a12 = -1 / (oldAlpha * oldAlpha);
          const double a21 = -oldAlpha / (alpha * alpha), a22 = alpha / (oldAlpha * oldAlpha);
          const double b1 = score1 - score - alpha * initslope, b2 = oldScore1 - score - oldAlpha * initslope;
          const double ab1 = (a11 * b1 + a12 * b2) / (alpha - oldAlpha), ab2 = (a21 * b1 + a22 * b2) / (alpha - oldAlpha);
          const double disc = ab2 * ab2 - 3 * ab1 * initslope;
          tmp = ab1 == 0 ? -initslope / ab2 / 2 : (-ab2 + std::sqrt(disc)) / (3 * ab1);
          if (tmp > 0.5 * alpha) tmp = 0.5 * alpha;
        }
        oldAlpha = alpha;
        oldScore1 = score1;
        alpha = tmp <= 0.1 * alpha ? 0.1 * alpha : tmp;
      }
      if (ii == maxHalf) break;
      if (retcode < 2) break;
    }
    if (!bb.ok) return Rcpp::List::create(Rcpp::Named("ok") = false);
    std::vector<double> step(m), oldRho = rho, oldGrad = grad;
    for (int i = 0; i < m; ++i) step[i] = alpha * Nstep[i];
    const double oldScore = score;
    rho = rho1;
    score = score1;
    grad = bb.grad;
    b = bb;
    std::vector<double> yg(m);
    for (int i = 0; i < m; ++i) yg[i] = grad[i] - oldGrad[i];
    bool skip = true;
    for (int j = 0; j < m; ++j) {
      double By = 0;
      for (int k = 0; k < m; ++k) By += B[j * m + k] * yg[k];
      if (std::fabs(step[j] - By) >= gradtol * std::max(std::fabs(grad[j]), std::fabs(oldGrad[j]))) skip = false;
    }
    if (!curv) skip = true;
    if (!skip) {
      if (it == 1) for (double& v : B) v *= alpha;
      double ys = 0;
      for (int i = 0; i < m; ++i) ys += yg[i] * step[i];
      const double rr = 1 / ys;
      // B <- B - rr * step (yg' B); B <- B - rr * (B yg) step' + rr * step step'
      std::vector<double> ygB(m, 0.0), Byg(m, 0.0);
      for (int j = 0; j < m; ++j) for (int k = 0; k < m; ++k) ygB[j] += yg[k] * B[k * m + j];
      for (int i = 0; i < m; ++i) for (int j = 0; j < m; ++j) B[i * m + j] -= rr * step[i] * ygB[j];
      for (int i = 0; i < m; ++i) for (int k = 0; k < m; ++k) Byg[i] += B[i * m + k] * yg[k];
      for (int i = 0; i < m; ++i) for (int j = 0; j < m; ++j) B[i * m + j] += -rr * Byg[i] * step[j] + rr * step[i] * step[j];
    }
    double rhoScale = 1, gmax = 0, dmax = 0;
    for (int i = 0; i < m; ++i) rhoScale = std::max(rhoScale, std::fabs(rho[i]));
    for (int i = 0; i < m; ++i) {
      gmax = std::max(gmax, std::fabs(grad[i]) * rhoScale / std::max(std::fabs(score), 1.0));
      dmax = std::max(dmax, std::fabs(rho[i] - oldRho[i]) / rhoScale);
    }
    int term = 0;
    if (retcode == 1) term = gmax <= gradtol * 6.0554 ? 1 : 3;
    else if (gmax <= gradtol * 6.0554) term = 1;
    else if (dmax <= steptol) term = 2;
    else if (it == maxStep) term = 4;
    else if (maxtaken) { if (++consecmax == 5) term = 5; }
    else consecmax = 0;
    if (term > 0) break;
    const double sc = std::fabs(b.scaleEst) + std::fabs(score);
    bool conv = true;
    for (int j = 0; j < m; ++j) { unconv[j] = std::fabs(grad[j]) > sc * gradtol; if (unconv[j]) conv = false; }
    if (std::fabs(oldScore - score) > sc * gradtol) {
      if (conv) std::fill(unconv.begin(), unconv.end(), 1);
    }
  }
  return Rcpp::List::create(Rcpp::Named("ok") = true, Rcpp::Named("coef") = Rcpp::wrap(b.fit.bt),
                            Rcpp::Named("rho") = Rcpp::wrap(rho), Rcpp::Named("gcv") = score);
}

#include "spicyglm/core.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace spicyglm {

namespace {

template <typename Count>
GlmFit fit_poisson_counts(const std::vector<Count>& n, const std::vector<double>& density,
                          const std::vector<int>& group, const std::string& estimator) {
  double pseudo;
  if (estimator == "mle") pseudo = 0.0;
  else if (estimator == "firth") pseudo = 0.5;
  else throw std::invalid_argument("estimator must be 'mle' or 'firth'");

  std::array<double, 2> Y{0.0, 0.0}, D{0.0, 0.0};
  for (std::size_t i = 0; i < n.size(); ++i) {
    Y[group[i]] += n[i];
    D[group[i]] += density[i];
  }
  GlmFit fit;
  for (int g = 0; g < 2; ++g) fit.beta[g] = std::log((Y[g] + pseudo) / D[g]);
  fit.mu.resize(n.size());
  for (std::size_t i = 0; i < n.size(); ++i) fit.mu[i] = std::exp(fit.beta[group[i]]) * density[i];
  return fit;
}

}  // namespace

GlmFit fit_poisson(const std::vector<int>& n, const std::vector<double>& density,
                   const std::vector<int>& group, const std::string& estimator) {
  return fit_poisson_counts(n, density, group, estimator);
}

GlmFit fit_poisson(const std::vector<double>& n, const std::vector<double>& density,
                   const std::vector<int>& group, const std::string& estimator) {
  return fit_poisson_counts(n, density, group, estimator);
}

double naive_variance(const std::vector<int>& group, const std::vector<double>& var) {
  std::array<double, 2> S{0.0, 0.0};
  for (std::size_t c = 0; c < group.size(); ++c) S[group[c]] += var[c];
  return 1.0 / S[0] + 1.0 / S[1];
}

namespace {

// One cluster's cells, grouped by image (Section 6.1) and into CR2 blocks:
// cells sharing a working variance, on which A_i reduces exactly. With an
// image-level offset the blocks are the images; a per-cell offset (the
// inhomogeneous design) makes every cell its own block (Section 14.7).
struct Patient {
  int group = 0;
  std::vector<int> image_id;       // image code of each local image
  std::vector<double> n_ij;        // cells per image
  std::vector<double> var_ij;      // working variance of each image's first cell
  std::vector<int> local_image;    // local image index of each cell
  std::vector<double> var;         // working variance of each cell
  std::vector<double> resid;       // residual of each cell
  std::vector<double> n_b, var_b;  // cells and working variance per block
  std::vector<int> local_block;    // block of each cell
  Eigen::MatrixXd Abar;            // reduced correction, A_i = I + P_i Abar P_i' (eq. 11)
};

// Abar = Dbar Vbar Lambda^{-1/2} Vbar' Dbar - I, from the eigenpairs of
// Gbar = diag(v^2) - f f' / S_g with f_b = sqrt(n_b) v_b^{3/2} (Prop. 6).
Eigen::MatrixXd reduced_correction(const Patient& p, double S_g) {
  Eigen::Map<const Eigen::VectorXd> n(p.n_b.data(), p.n_b.size());
  Eigen::Map<const Eigen::VectorXd> v(p.var_b.data(), p.var_b.size());
  Eigen::VectorXd f = (n.array().sqrt() * v.array().pow(1.5)).matrix();
  Eigen::MatrixXd Gbar = Eigen::MatrixXd(v.array().square().matrix().asDiagonal()) - f * f.transpose() / S_g;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Gbar);
  if (eig.info() != Eigen::Success || eig.eigenvalues().minCoeff() <= 0.0)
    throw std::runtime_error("CR2 inner matrix is not positive definite (a group has a single cluster?)");
  Eigen::MatrixXd DV = v.array().sqrt().matrix().asDiagonal() * eig.eigenvectors();
  Eigen::VectorXd inv_sqrt = eig.eigenvalues().array().rsqrt();
  Eigen::MatrixXd Abar = DV * inv_sqrt.asDiagonal() * DV.transpose();
  Abar.diagonal().array() -= 1.0;
  return Abar;
}

// A_i u for a per-cell vector u: u + P_i Abar P_i' u.
Eigen::VectorXd apply_A(const Patient& p, const Eigen::VectorXd& u) {
  std::size_t n_blocks = p.n_b.size();
  Eigen::VectorXd d = Eigen::VectorXd::Zero(n_blocks);
  for (Eigen::Index c = 0; c < u.size(); ++c) d[p.local_block[c]] += u[c];
  for (std::size_t j = 0; j < n_blocks; ++j) d[j] /= std::sqrt(p.n_b[j]);
  Eigen::VectorXd corr = p.Abar * d;
  Eigen::VectorXd out = u;
  for (Eigen::Index c = 0; c < u.size(); ++c) {
    int j = p.local_block[c];
    out[c] += corr[j] / std::sqrt(p.n_b[j]);
  }
  return out;
}

// A_i U for the columns of U without forming Abar: Abar = Dbar C Dbar, where C
// is the resolvent integral (2 rho / pi) int_0^Inf E_t f f' E_t / (1 - rho f' E_t f) dt,
// E_t = (diag(v^2) + t^2)^{-1}, rho = 1 / S_g, evaluated by the trapezoid rule
// in u = log t (Section 14.8). The columns share each node's E_t f; O(blocks)
// per node and column, and no eigendecomposition.
Eigen::MatrixXd apply_A_integral(const Patient& p, const Eigen::MatrixXd& U, double S_g) {
  const double pi = 3.141592653589793238462643383279502884;
  const double du = 0.2;
  const double rho = 1.0 / S_g;
  std::size_t n_blocks = p.n_b.size();
  Eigen::Index q = U.cols();
  std::vector<double> Lambda(n_blocks), f(n_blocks), dbar(n_blocks), root_n(n_blocks), ef(n_blocks);
  double lev = 0.0;  // the patient's leverage share
  for (std::size_t b = 0; b < n_blocks; ++b) {
    double v = p.var_b[b];
    Lambda[b] = v * v;
    root_n[b] = std::sqrt(p.n_b[b]);
    f[b] = root_n[b] * std::pow(v, 1.5);
    dbar[b] = std::sqrt(v);
    lev += f[b] * f[b] / Lambda[b];
  }
  lev *= rho;
  if (!(lev < 1.0))
    throw std::runtime_error("CR2 inner matrix is not positive definite (a group has a single cluster?)");

  // Raw arrays in the loops below: this runs ~300 times per patient, and
  // unoptimised builds (devtools::load_all) make every checked element access
  // a function call. Z = Dbar P_i' U and acc are n_blocks x q, column-major.
  const Eigen::Index rows = U.rows();
  const double* u = U.data();
  const int* block = p.local_block.data();
  std::vector<double> Zv(n_blocks * q, 0.0), accv(n_blocks * q, 0.0), dotv(q);
  double* Z = Zv.data();
  double* acc = accv.data();
  double* dot = dotv.data();
  const double* L = Lambda.data();
  const double* fb = f.data();
  double* e = ef.data();
  for (Eigen::Index j = 0; j < q; ++j)
    for (Eigen::Index c = 0; c < rows; ++c) Z[j * n_blocks + block[c]] += u[j * rows + c];
  for (Eigen::Index j = 0; j < q; ++j)
    for (std::size_t b = 0; b < n_blocks; ++b) Z[j * n_blocks + b] *= dbar[b] / root_n[b];

  auto [lmin, lmax] = std::minmax_element(Lambda.begin(), Lambda.end());
  double lo = 0.5 * std::log(*lmin * (1.0 - lev)) - 36.0;
  double hi = 0.5 * std::log(*lmax) + 13.0;
  long n_nodes = static_cast<long>((hi - lo) / du + 1e-10);  // as R's seq(lo, hi, by = du)
  for (long k = 0; k <= n_nodes; ++k) {
    double uk = lo + k * du, t2 = std::exp(2.0 * uk);
    double wt = (2.0 / pi) * du * std::exp(uk) * rho;
    double ff = 0.0;
    for (std::size_t b = 0; b < n_blocks; ++b) {
      e[b] = fb[b] / (L[b] + t2);
      ff += e[b] * fb[b];
    }
    for (Eigen::Index j = 0; j < q; ++j) {
      const double* Zj = Z + j * n_blocks;
      double s = 0.0;
      for (std::size_t b = 0; b < n_blocks; ++b) s += e[b] * Zj[b];
      dot[j] = s;
    }
    for (Eigen::Index j = 0; j < q; ++j) {
      double coef = wt * dot[j] / (1.0 - rho * ff);
      double* accj = acc + j * n_blocks;
      for (std::size_t b = 0; b < n_blocks; ++b) accj[b] += coef * e[b];
    }
  }
  Eigen::MatrixXd out = U;
  double* o = out.data();
  for (Eigen::Index j = 0; j < q; ++j)
    for (Eigen::Index c = 0; c < rows; ++c) {
      int b = block[c];
      o[j * rows + c] += dbar[b] * acc[j * n_blocks + b] / root_n[b];
    }
  return out;
}

}  // namespace

CR2Result cr2_wald(const std::vector<int>& cluster, const std::vector<int>& image,
                   const std::vector<int>& group, const std::vector<double>& var,
                   const std::vector<double>& resid) {
  std::size_t N = cluster.size();
  if (image.size() != N || group.size() != N || var.size() != N || resid.size() != N)
    throw std::invalid_argument("cr2_wald: input lengths differ");

  // group cells into clusters and images, in first-appearance order
  std::vector<Patient> patients;
  std::vector<int> cluster_codes;
  std::unordered_map<int, std::size_t> cluster_pos;
  std::vector<std::unordered_map<int, int>> image_pos;
  for (std::size_t c = 0; c < N; ++c) {
    auto [it, fresh] = cluster_pos.try_emplace(cluster[c], patients.size());
    if (fresh) {
      patients.emplace_back();
      patients.back().group = group[c];
      cluster_codes.push_back(cluster[c]);
      image_pos.emplace_back();
    }
    Patient& p = patients[it->second];
    if (p.group != group[c]) throw std::invalid_argument("cr2_wald: a cluster spans both groups");
    auto [jt, new_image] = image_pos[it->second].try_emplace(image[c], static_cast<int>(p.n_ij.size()));
    if (new_image) {
      p.image_id.push_back(image[c]);
      p.n_ij.push_back(0.0);
      p.var_ij.push_back(var[c]);
    }
    p.n_ij[jt->second] += 1.0;
    p.local_image.push_back(jt->second);
    p.var.push_back(var[c]);
    p.resid.push_back(resid[c]);
  }

  // CR2 blocks (Section 14.7). Cell-level blocks have no image-level reduction,
  // so every patient then takes the integral route, as spicyR's buildGLM() does.
  bool integral = false;
  for (Patient& p : patients) {
    bool image_blocks = true;
    for (std::size_t c = 0; c < p.var.size(); ++c)
      image_blocks = image_blocks && p.var[c] == p.var_ij[p.local_image[c]];
    if (image_blocks) {
      p.n_b = p.n_ij;
      p.var_b = p.var_ij;
      p.local_block = p.local_image;
    } else {
      p.n_b.assign(p.var.size(), 1.0);
      p.var_b = p.var;
      p.local_block.resize(p.var.size());
      std::iota(p.local_block.begin(), p.local_block.end(), 0);
      integral = true;
    }
  }

  CR2Result out;
  out.S = {0.0, 0.0};
  for (const Patient& p : patients)
    for (std::size_t j = 0; j < p.n_b.size(); ++j) out.S[p.group] += p.n_b[j] * p.var_b[j];

  const std::array<double, 2> w{-1.0 / out.S[0], 1.0 / out.S[1]};  // L B with L = (-1, 1)
  std::size_t m = patients.size();
  Eigen::VectorXd P_diag(m), h(m);
  out.v_hat = 0.0;

  for (std::size_t i = 0; i < m; ++i) {
    Patient& p = patients[i];
    double S_g = out.S[p.group];
    std::size_t Ni = p.resid.size();

    Eigen::Map<const Eigen::VectorXd> r(p.resid.data(), Ni);
    Eigen::VectorXd v_cell(Ni);
    for (std::size_t c = 0; c < Ni; ++c) v_cell[c] = p.var_b[p.local_block[c]];

    Eigen::VectorXd Ar, A1, Av;
    if (integral) {
      Eigen::MatrixXd U(Ni, 3);
      U << r, Eigen::VectorXd::Ones(Ni), v_cell;
      Eigen::MatrixXd AU = apply_A_integral(p, U, S_g);
      Ar = AU.col(0);
      A1 = AU.col(1);
      Av = AU.col(2);
    } else {
      p.Abar = reduced_correction(p, S_g);
      Ar = apply_A(p, r);
      A1 = apply_A(p, Eigen::VectorXd::Ones(Ni));
      Av = apply_A(p, v_cell);
    }

    double e = w[p.group] * Ar.sum();
    out.e.push_back(e);
    out.v_hat += e * e;

    // Satterthwaite pieces (Section 9.5, Remark 15: tau uses the variance)
    P_diag[i] = w[p.group] * w[p.group] * (v_cell.array() * A1.array().square()).sum();
    h[i] = Av.sum() / std::pow(S_g, 1.5);

    std::vector<double> adj(p.n_ij.size(), 0.0), raw(p.n_ij.size(), 0.0);
    for (std::size_t c = 0; c < Ni; ++c) {
      adj[p.local_image[c]] += Ar[c];
      raw[p.local_image[c]] += r[c];
    }
    out.adjusted_image_sum.push_back(std::move(adj));
    out.raw_image_sum.push_back(std::move(raw));
    out.image_id.push_back(p.image_id);
    out.group.push_back(p.group);
  }
  out.cluster_id = cluster_codes;

  // nu = tr(P)^2 / sum(P^2), P_jk = -h_j h_k within a group, P_jj = P_diag - h_j^2
  Eigen::MatrixXd P = Eigen::MatrixXd::Zero(m, m);
  for (std::size_t j = 0; j < m; ++j)
    for (std::size_t k = 0; k < m; ++k)
      if (patients[j].group == patients[k].group) P(j, k) = -h[j] * h[k];
  P.diagonal() = P_diag - h.array().square().matrix();
  double tr = P.trace();
  out.df = tr * tr / P.squaredNorm();
  return out;
}

}  // namespace spicyglm

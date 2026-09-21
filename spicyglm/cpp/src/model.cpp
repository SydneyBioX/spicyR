#include "spicyglm/core.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace spicyglm {

GlmFit fit_poisson(const std::vector<int>& n, const std::vector<double>& density,
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

double naive_variance(const std::vector<int>& group, const std::vector<double>& var) {
  std::array<double, 2> S{0.0, 0.0};
  for (std::size_t c = 0; c < group.size(); ++c) S[group[c]] += var[c];
  return 1.0 / S[0] + 1.0 / S[1];
}

namespace {

// One cluster's cells, grouped by image (Section 6.1).
struct Patient {
  int group = 0;
  std::vector<int> image_id;       // image code of each local image
  std::vector<double> n_ij;        // cells per image
  std::vector<double> var_ij;      // working variance per image
  std::vector<int> local_image;    // local image index of each cell
  std::vector<double> resid;       // residual of each cell
  Eigen::MatrixXd Abar;            // reduced correction, A_i = I + P_i Abar P_i' (eq. 11)
};

// Abar = Dbar Vbar Lambda^{-1/2} Vbar' Dbar - I, from the eigenpairs of
// Gbar = diag(v^2) - f f' / S_g with f_j = sqrt(n_ij) v_ij^{3/2} (Prop. 6).
Eigen::MatrixXd reduced_correction(const Patient& p, double S_g) {
  Eigen::Map<const Eigen::VectorXd> n(p.n_ij.data(), p.n_ij.size());
  Eigen::Map<const Eigen::VectorXd> v(p.var_ij.data(), p.var_ij.size());
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
  std::size_t n_img = p.n_ij.size();
  Eigen::VectorXd d = Eigen::VectorXd::Zero(n_img);
  for (Eigen::Index c = 0; c < u.size(); ++c) d[p.local_image[c]] += u[c];
  for (std::size_t j = 0; j < n_img; ++j) d[j] /= std::sqrt(p.n_ij[j]);
  Eigen::VectorXd corr = p.Abar * d;
  Eigen::VectorXd out = u;
  for (Eigen::Index c = 0; c < u.size(); ++c) {
    int j = p.local_image[c];
    out[c] += corr[j] / std::sqrt(p.n_ij[j]);
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
    p.resid.push_back(resid[c]);
  }

  CR2Result out;
  out.S = {0.0, 0.0};
  for (const Patient& p : patients)
    for (std::size_t j = 0; j < p.n_ij.size(); ++j) out.S[p.group] += p.n_ij[j] * p.var_ij[j];

  const std::array<double, 2> w{-1.0 / out.S[0], 1.0 / out.S[1]};  // L B with L = (-1, 1)
  std::size_t m = patients.size();
  Eigen::VectorXd P_diag(m), h(m);
  out.v_hat = 0.0;

  for (std::size_t i = 0; i < m; ++i) {
    Patient& p = patients[i];
    double S_g = out.S[p.group];
    p.Abar = reduced_correction(p, S_g);
    std::size_t Ni = p.resid.size();

    Eigen::Map<const Eigen::VectorXd> r(p.resid.data(), Ni);
    Eigen::VectorXd v_cell(Ni);
    for (std::size_t c = 0; c < Ni; ++c) v_cell[c] = p.var_ij[p.local_image[c]];

    Eigen::VectorXd Ar = apply_A(p, r);
    Eigen::VectorXd A1 = apply_A(p, Eigen::VectorXd::Ones(Ni));
    Eigen::VectorXd Av = apply_A(p, v_cell);

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

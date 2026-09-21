#include "spicyglm/core.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace spicyglm {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

PairDiagnostics compute_diagnostics(const std::vector<int>& cluster, const std::vector<int>& image,
                                    const std::vector<int>& group, const std::vector<int>& n,
                                    const std::vector<double>& density, const CR2Result& cr2) {
  std::size_t m = cr2.cluster_id.size();
  std::unordered_map<int, std::size_t> cluster_pos;
  std::vector<std::unordered_map<int, std::size_t>> image_pos(m);
  std::vector<std::size_t> image_offset(m + 1, 0);
  for (std::size_t i = 0; i < m; ++i) {
    cluster_pos[cr2.cluster_id[i]] = i;
    for (std::size_t j = 0; j < cr2.image_id[i].size(); ++j) image_pos[i][cr2.image_id[i][j]] = j;
    image_offset[i + 1] = image_offset[i] + cr2.image_id[i].size();
  }

  PairDiagnostics D;
  std::size_t n_images = image_offset[m];
  D.n_i.assign(m, 0.0); D.T.assign(m, 0.0); D.y.assign(m, 0.0); D.d.assign(m, 0.0);
  D.n_ij.assign(n_images, 0.0); D.density_ij.assign(n_images, 0.0);
  for (std::size_t c = 0; c < cluster.size(); ++c) {
    std::size_t i = cluster_pos.at(cluster[c]);
    std::size_t k = image_offset[i] + image_pos[i].at(image[c]);
    D.n_i[i] += 1.0;
    D.y[i] += n[c];
    D.d[i] += density[c];
    D.n_ij[k] += 1.0;
    D.density_ij[k] = density[c];  // constant within an image
  }

  // leverage (pre-fit; exp(beta) cancels within a group)
  D.S_leverage = {0.0, 0.0};
  std::array<double, 2> Y{0.0, 0.0}, Dg{0.0, 0.0};
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t k = image_offset[i]; k < image_offset[i + 1]; ++k) D.T[i] += D.n_ij[k] * D.density_ij[k];
    D.S_leverage[cr2.group[i]] += D.T[i];
    Y[cr2.group[i]] += D.y[i];
    Dg[cr2.group[i]] += D.d[i];
  }
  const std::array<double, 2> w{-1.0 / cr2.S[0], 1.0 / cr2.S[1]};

  for (std::size_t i = 0; i < m; ++i) {
    int g = cr2.group[i];
    double e_i = cr2.e[i];
    D.l.push_back(D.T[i] / D.S_leverage[g]);
    D.influence.push_back(e_i * e_i / cr2.v_hat);

    bool sole = D.y[i] == Y[g] && D.d[i] == Dg[g];
    double beta_g = std::log((Y[g] + 0.5) / Dg[g]);
    D.delta.push_back(sole ? NaN : beta_g - std::log((Y[g] - D.y[i] + 0.5) / (Dg[g] - D.d[i])));

    double raw = 0.0, adjusted = 0.0;
    for (std::size_t j = 0; j < cr2.image_id[i].size(); ++j) {
      std::size_t k = image_offset[i] + j;
      double e_ij = w[g] * cr2.adjusted_image_sum[i][j];
      D.image_cluster.push_back(static_cast<int>(i));
      D.image_id.push_back(cr2.image_id[i][j]);
      D.l_ij.push_back(D.n_ij[k] * D.density_ij[k] / D.T[i]);
      D.l_ij_group_share.push_back(D.l[i] * D.l_ij.back());
      D.raw_sum_ij.push_back(cr2.raw_image_sum[i][j]);
      D.adjusted_sum_ij.push_back(cr2.adjusted_image_sum[i][j]);
      D.e_ij.push_back(e_ij);
      D.e_share_ij.push_back(e_i == 0.0 ? NaN : e_ij / e_i);
      D.influence_ij.push_back(e_i == 0.0 ? NaN : D.influence[i] * e_ij / e_i);
      raw += cr2.raw_image_sum[i][j];
      adjusted += cr2.adjusted_image_sum[i][j];
    }
    D.raw_sum.push_back(raw);
    D.adjusted_sum.push_back(adjusted);
  }
  return D;
}

}  // namespace

PairFit fit_pair_poisson(const std::vector<int>& cluster, const std::vector<int>& image,
                         const std::vector<int>& group, const std::vector<int>& n,
                         const std::vector<double>& density, const std::string& estimator,
                         const std::string& variance, bool diagnostics) {
  if (variance != "fast" && variance != "naive") throw std::invalid_argument("variance must be 'fast' or 'naive'");
  if (diagnostics && (estimator != "firth" || variance != "fast"))
    throw std::invalid_argument("diagnostics require estimator = 'firth' and variance = 'fast'");
  PairFit out;
  out.fit = fit_poisson(n, density, group, estimator);
  if (variance == "naive") {
    out.naive = true;
    out.v_naive = naive_variance(group, out.fit.mu);  // Poisson working variance = mean
    return out;
  }
  std::vector<double> resid(n.size());
  for (std::size_t c = 0; c < n.size(); ++c) resid[c] = n[c] - out.fit.mu[c];
  out.cr2 = cr2_wald(cluster, image, group, out.fit.mu, resid);
  if (diagnostics) {
    out.diagnostics = compute_diagnostics(cluster, image, group, n, density, out.cr2);
    out.has_diagnostics = true;
  }
  return out;
}

}  // namespace spicyglm

// spicyGLM numeric core. Equation and proposition numbers refer to
// spicyClub_Supplementary_Math.pdf (spicyR, gee branch).
#pragma once

#include <array>
#include <string>
#include <vector>

namespace spicyglm {

// ---------------------------------------------------------------- spatial ---

enum class Window { Rectangle, Convex };

// Area of the observation window of one image (all of its cells).
double window_area(const double* x, const double* y, std::size_t n, Window window);

// The k nearest neighbours of every cell within its own image (Euclidean, no
// edge correction, the cell itself excluded), as row indices in a flat
// N x k array, nearest first. Ties at equal distance are broken exactly as
// spatstat.geom's nnwhich(). Rows of images with at most k cells are -1.
// Images are processed on up to n_threads threads.
std::vector<int> knn_indices(const std::vector<double>& x, const std::vector<double>& y,
                             const std::vector<int>& image_offsets, int k, int n_threads = 1);

// Poisson model data for one (from, to) pair. Images missing either cell type
// contribute no rows.
struct ModelData {
  std::vector<int> row;         // row of each reference cell
  std::vector<int> image;       // image index of each reference cell
  std::vector<int> n;           // TARGET neighbours within r (inclusive, as spatstat::crosspairs)
  std::vector<double> density;  // expected neighbours under CSR: (n_to / area) * pi r^2
};

// Binomial model data for one pair (Section 12.1): for each reference cell,
// how many of its k nearest neighbours are TARGET. Images missing either type,
// with at most k cells, or made up only of TARGET cells contribute no rows.
struct BinomialModelData {
  std::vector<int> row, image, n;
  std::vector<double> p0;  // image background TARGET proportion
};

// Every cell of an analysis, indexed once and shared by all cell-type pairs.
// Cells must be sorted by image, image i occupying rows
// [image_offsets[i], image_offsets[i + 1]). Queries are const and thread-safe.
class Dataset {
 public:
  Dataset(std::vector<double> x, std::vector<double> y, std::vector<int> cell_type,
          std::vector<int> image_offsets, int n_types);

  std::vector<double> image_areas(Window window) const;

  // Per-image grids for radius counts (family poisson).
  void build_radius_index(double r);
  ModelData poisson_model_data(const std::vector<double>& image_area, int from, int to) const;

  // k-nearest-neighbour lists (family binomial).
  void build_knn(int k, int n_threads = 1);
  BinomialModelData binomial_model_data(int from, int to) const;

 private:
  struct Grid {
    double xmin = 0, ymin = 0, side = 1;
    long long nbx = 0, nby = 0;
    std::size_t bin_offset = 0;  // into bin_start_
  };
  int n_images() const { return static_cast<int>(image_offsets_.size()) - 1; }
  int type_count(int img, int type) const {
    return type_start_[img * n_types_ + type + 1] - type_start_[img * n_types_ + type];
  }

  std::vector<double> x_, y_;
  std::vector<int> type_, image_offsets_;
  int n_types_;
  std::vector<int> type_start_, type_rows_;  // rows of each (image, type), in row order

  double r_ = 0;
  std::vector<Grid> grids_;
  std::vector<int> bin_start_;     // per image, nbins + 1 positions into the sorted arrays
  std::vector<double> gx_, gy_;    // coordinates sorted by (image, bin, type)
  std::vector<int> gtype_;

  int k_ = 0;
  std::vector<int> knn_;
};

// ---------------------------------------------------------------- fitting ---

struct GlmFit {
  std::array<double, 2> beta;
  std::vector<double> mu;  // fitted mean on the count scale
};

// Closed-form fit of n ~ 0 + group, offset log(density); group is 0 or 1.
// estimator "mle": log(Y_g / D_g); "firth": log((Y_g + 0.5) / D_g) (eq. 16).
GlmFit fit_poisson(const std::vector<int>& n, const std::vector<double>& density,
                   const std::vector<int>& group, const std::string& estimator);

// Fit of n ~ Binomial(k, pi), logit(pi) = beta_g + logit(p0). The design is
// one coefficient per group, so each group is a one-dimensional root-find.
// estimator "firth" solves the Jeffreys-penalised score (Section 12.4).
GlmFit fit_binomial(const std::vector<int>& n, int k, const std::vector<double>& p0,
                    const std::vector<int>& group, const std::string& estimator);

// -------------------------------------------------------------------- CR2 ---

// CR2 variance of beta_2 - beta_1 and its Satterthwaite degrees of freedom,
// clustered by `cluster`. var is the working variance per cell (the fitted
// mean for Poisson) and must be constant within an image.
struct CR2Result {
  std::array<double, 2> S;                            // group totals S_g
  std::vector<int> cluster_id;                        // cluster code, in first-appearance order
  std::vector<int> group;                             // group of each cluster
  std::vector<double> e;                              // e_i = w_g 1' A_i r_i
  std::vector<std::vector<int>> image_id;             // images of each cluster
  std::vector<std::vector<double>> adjusted_image_sum;  // sum of (A_i r_i) within each image
  std::vector<std::vector<double>> raw_image_sum;       // sum of r_i within each image
  double v_hat;                                       // L V_CR2 L' (eq. 17)
  double df;                                          // Satterthwaite nu (Prop. 11)
};

CR2Result cr2_wald(const std::vector<int>& cluster, const std::vector<int>& image,
                   const std::vector<int>& group, const std::vector<double>& var,
                   const std::vector<double>& resid);

// ------------------------------------------------------------ diagnostics ---

// Per-pair QC diagnostics (Section 11). Cluster rows follow CR2Result's
// cluster order; image rows are flattened in cluster order, then each
// cluster's image order. NaN marks an undefined value.
struct PairDiagnostics {
  std::array<double, 2> S_leverage;  // sum of T_i per group
  std::vector<double> n_i, T, l;     // cells, T_i = sum_j n_ij den_ij, l_i = T_i / S_g (Def. 3)
  std::vector<double> raw_sum, adjusted_sum, influence;  // Infl_i = e_i^2 / v_hat (Def. 5)
  std::vector<double> y, d, delta;   // leave-one-out Firth shift (Def. 7)

  std::vector<int> image_cluster, image_id;
  std::vector<double> n_ij, density_ij, l_ij, l_ij_group_share;  // Def. 4
  std::vector<double> raw_sum_ij, adjusted_sum_ij, e_ij, e_share_ij, influence_ij;  // Def. 6
};

// Model-based ("naive") variance of beta_2 - beta_1, ignoring clustering: the
// inverse Fisher information is diag(1 / S_1, 1 / S_2), S_g the summed working
// variance of group g, so the contrast variance is 1 / S_1 + 1 / S_2.
double naive_variance(const std::vector<int>& group, const std::vector<double>& var);

struct PairFit {
  GlmFit fit;
  bool naive = false;       // variance "naive": v_naive set, cr2 not computed
  double v_naive = 0.0;
  CR2Result cr2;
  bool has_diagnostics = false;
  PairDiagnostics diagnostics;
};

// Fit one pair and compute the variance of the contrast: variance "fast" runs
// CR2 with Satterthwaite degrees of freedom, "naive" the model-based variance.
// Diagnostics require estimator "firth" (the point-estimate shift uses the
// closed-form Firth estimate) and variance "fast".
PairFit fit_pair_poisson(const std::vector<int>& cluster, const std::vector<int>& image,
                         const std::vector<int>& group, const std::vector<int>& n,
                         const std::vector<double>& density, const std::string& estimator,
                         const std::string& variance, bool diagnostics);

// Binomial fit with the variance substitution mu -> k pi (1 - pi) (Prop. 17).
// Diagnostics are not available for this family.
PairFit fit_pair_binomial(const std::vector<int>& cluster, const std::vector<int>& image,
                          const std::vector<int>& group, const std::vector<int>& n, int k,
                          const std::vector<double>& p0, const std::string& estimator,
                          const std::string& variance);

}  // namespace spicyglm

// spicyGLM numeric core. Equation and proposition numbers refer to
// the Supplementary Methods of the paper describing spicyR Cell (in preparation).
// Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.
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

// Inhomogeneous (intensity-reweighted) Poisson model data for one pair
// (Section 14.4). Every REF-TARGET pair within r counts c_ab * w_a * w_b, where
// c_ab is the translation edge correction and w the mean-one inverse disc-kernel
// intensity of a cell's type at the cell; the offset keeps the homogeneous
// image total, spread over REF cells in proportion to w_a, so that
// sum(n) / sum(density) is the image's inhomogeneous cross-K over pi r^2. A
// self-pair excludes the cell itself and has n_to - 1 target cells; images
// missing either type, or with fewer than two cells of a self-pair's type,
// contribute no rows.
struct InhomModelData {
  std::vector<int> row, image;
  std::vector<double> n;        // pair-weighted TARGET neighbours within r
  std::vector<int> n_raw;       // unweighted TARGET neighbours within r
  std::vector<double> weight;   // w_a
  std::vector<double> density;  // w_a * (n_to - [from == to]) / area * pi r^2
};

// Kontextual binomial model data for one pair (Section 15): for each REF cell,
// how many of its k nearest neighbours are context cells (the trials) and how
// many of those are TARGET (the successes), against the image's TARGET share
// of the context. REF cells with no context cell among their neighbours, and
// images missing either type, with at most k cells, or whose context is all
// TARGET, contribute no rows.
struct BinomialTrialsModelData {
  std::vector<int> row, image, n, trials;
  std::vector<double> p0;  // n_to / n_context of the image
};

// An image's window as a polygon (Section 14.10): its vertices anticlockwise
// and, for a convex hull, the centred vertices, edge vectors, outward unit
// normals and the determinants joining consecutive edge lines, which give the
// translation overlap |W intersect (W + v)| in O(#vertices).
struct WindowGeometry {
  bool rectangle = true;
  double area = 0;
  std::vector<double> vx, vy;
  std::vector<double> px, py, ex, ey, nx, ny, det;
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

  // Sufficient statistics of every pair at once, for image-level tests: entry
  // (img * n_types + from) * n_types + to is the sum over the image's `from`
  // cells of their `to` neighbours, within r (radius index, the cell itself
  // included as in poisson_model_data) or among the k nearest (kNN). One pass
  // over the neighbour graph.
  std::vector<double> pair_neighbour_totals(bool knn) const;
  // The same layout, entry (img, from, to) = sum over the image's `to` cells b
  // of c_b^2, where c_b is the number of `from` cells that count b as a
  // neighbour (other than b itself). Under independent random labelling of the
  // non-`from` cells, Var(O) = p (1 - p) sum_b c_b^2 for the pair's total O.
  std::vector<double> pair_neighbour_sq_totals(bool knn) const;
  // The same layout with c_b counted from b's side: the number of `from` cells
  // among b's own neighbours (its k nearest, or those within r). For the radius
  // graph this equals pair_neighbour_sq_totals; for k-NN it is the out-degree
  // version the TARGET-centred excess design (effect = "excess") needs.
  std::vector<double> pair_neighbour_out_sq_totals(bool knn) const;
  // Poisson design under the random-labelling null (Section 17): each REF cell's
  // offset is its number of candidate neighbours within r times the TARGET share
  // of the candidates. Candidates are the non-REF cells, or for a self-pair all
  // other cells, with share (n_A - 1) / (N - 1). Needs build_radius_index.
  ModelData rl_model_data(int from, int to) const;
  // Within-image variance of the weighted designs (Section 17): for every image,
  // sum_b c_b and sum_b c_b^2 over the candidate cells b, where c_b is the total
  // weight the REF cells give b as a neighbour, plus the number of candidates.
  // design 0: Kontextual Poisson (candidates: context cells other than REF-type
  // ones, or all other context cells for a self-pair; weight e_a lambda_c(a) /
  // lambda_c(b)). design 1: Kontextual binomial (candidates as for 0; weight 1 per
  // k-NN slot). design 2: inhomogeneous (candidates: the TARGET cells; weight
  // e_ab w_a w_b, as in inhom_model_data). Output rows: image-major, 3 per image.
  std::vector<double> weighted_phi_sums(int from, int to, int design, bool edge_correct) const;
  // Spatial HAC (Conley) estimate of Var(O_i) for every image, allowing the labels of
  // nearby candidates to be correlated: sum_b sum_b' k(d_bb') c_b c_b' e_b e_b' over
  // candidates within bandwidth h, Bartlett k(d) = 1 - d / h, with e_b the residual of
  // the TARGET indicator on c_b within the image (so the co-localisation itself is not
  // counted as variance). design 0/1: Kontextual Poisson / binomial (as in
  // weighted_phi_sums); 3: radius design without context; 4: k-NN design without
  // context. Output: image-major, 3 per image (HAC variance, sum_b c_b, candidates).
  std::vector<double> hac_phi_sums(int from, int to, int design, double h) const;
  // hac_phi_sums for one REF and every non-self TARGET in one pass: per image, the T HAC sums,
  // then sum c, n and G (layout T + 3 per image).
  std::vector<double> hac_phi_sums_ref(int from, int design, double h) const;

  // Inhomogeneous design (Section 14): the window of each image as a polygon,
  // and every cell's weight 1 / lambda, scaled to mean one within its (image,
  // type). lambda is the number of other cells of the type within sigma over
  // the area of that disc inside the window, floored at min_lambda times the
  // type's average intensity. Needs build_radius_index for the pairs.
  void build_intensity(double sigma, double min_lambda, Window window);
  InhomModelData inhom_model_data(const std::vector<double>& image_area, int from, int to,
                                  bool edge_correct) const;

  // Kontextual design (Section 15): the context (parent) population, which
  // must contain the TARGET type. With a radius index built, also every cell's
  // count of context cells within r (itself included) and the area of that
  // disc inside its image's window (pi r^2 without edge correction).
  void build_context(const std::vector<int>& context_types, Window window, bool edge_correct);
  // Poisson: every TARGET neighbour b of REF cell a within r counts
  // e_a lambda_c(x_a) / lambda_c(x_b), with lambda_c the context count over the
  // disc area and e_a = pi r^2 / area at a; the offset is the TARGET share of
  // the context times lambda_c(x_a) pi r^2. A self-pair leaves the cell itself
  // out of every count. REF cells with no context cell within r contribute no rows.
  InhomModelData kontextual_model_data(int from, int to) const;
  // Binomial, over the k nearest neighbours (needs build_knn).
  BinomialTrialsModelData kontextual_binomial_model_data(int from, int to) const;

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
  // Per-image grids of bin side >= radius over all cells, sorted by (image, bin,
  // type), so every neighbour within radius lies in the surrounding 3 x 3 bins.
  void build_grid(double radius, std::vector<Grid>& grids, std::vector<int>& bin_start,
                  std::vector<double>& gx, std::vector<double>& gy, std::vector<int>& gtype,
                  std::vector<int>& grow) const;

  std::vector<double> x_, y_;
  std::vector<int> type_, image_offsets_;
  int n_types_;
  std::vector<int> type_start_, type_rows_;  // rows of each (image, type), in row order

  double r_ = 0;
  std::vector<Grid> grids_;
  std::vector<int> bin_start_;     // per image, nbins + 1 positions into the sorted arrays
  std::vector<double> gx_, gy_;    // coordinates sorted by (image, bin, type)
  std::vector<int> gtype_, grow_;  // and their cell types and rows

  int k_ = 0;
  std::vector<int> knn_;

  std::vector<WindowGeometry> windows_;
  std::vector<double> weight_;     // per row, once build_intensity has run
  void build_windows(Window window);

  std::vector<char> is_context_;   // per type, once build_context has run
  std::vector<int> context_count_; // per row: context cells within r, itself included
  std::vector<double> context_area_;  // per row: |b(x, r) intersect W|, or pi r^2
};

// ---------------------------------------------------------------- fitting ---

struct GlmFit {
  std::array<double, 2> beta;
  std::vector<double> mu;  // fitted mean on the count scale
};

// Closed-form fit of n ~ 0 + group, offset log(density); group is 0 or 1.
// estimator "mle": log(Y_g / D_g); "firth": log((Y_g + 0.5) / D_g) (eq. 16).
// The double overload takes the inhomogeneous design's pair-weighted counts.
GlmFit fit_poisson(const std::vector<int>& n, const std::vector<double>& density,
                   const std::vector<int>& group, const std::string& estimator);
GlmFit fit_poisson(const std::vector<double>& n, const std::vector<double>& density,
                   const std::vector<int>& group, const std::string& estimator);

// Fit of n ~ Binomial(k, pi), logit(pi) = beta_g + logit(p0). The design is
// one coefficient per group, so each group is a one-dimensional root-find.
// estimator "firth" solves the Jeffreys-penalised score (Section 12.4).
GlmFit fit_binomial(const std::vector<int>& n, int k, const std::vector<double>& p0,
                    const std::vector<int>& group, const std::string& estimator);
// The same with a number of trials per cell (the Kontextual binomial design).
GlmFit fit_binomial(const std::vector<int>& n, const std::vector<int>& trials,
                    const std::vector<double>& p0, const std::vector<int>& group,
                    const std::string& estimator);

// -------------------------------------------------------------------- CR2 ---

// CR2 variance of beta_2 - beta_1 and its Satterthwaite degrees of freedom,
// clustered by `cluster`. var is the working variance per cell (the fitted
// mean for Poisson). When it is constant within every image, A_i reduces to
// one dimension per image and is formed by eigendecomposition; otherwise
// (the inhomogeneous design) every cell is its own block and A_i is applied by
// the resolvent-integral quadrature, linear in the number of cells (Section 14.8).
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
PairFit fit_pair_poisson(const std::vector<int>& cluster, const std::vector<int>& image,
                         const std::vector<int>& group, const std::vector<double>& n,
                         const std::vector<double>& density, const std::string& estimator,
                         const std::string& variance, bool diagnostics);

// Binomial fit with the variance substitution mu -> k pi (1 - pi) (Prop. 17).
// Diagnostics are not available for this family.
PairFit fit_pair_binomial(const std::vector<int>& cluster, const std::vector<int>& image,
                          const std::vector<int>& group, const std::vector<int>& n, int k,
                          const std::vector<double>& p0, const std::string& estimator,
                          const std::string& variance);
// Per-cell trials: the working variance t_c pi (1 - pi) then varies by cell,
// and CR2 takes the integral route (Section 14.8).
PairFit fit_pair_binomial(const std::vector<int>& cluster, const std::vector<int>& image,
                          const std::vector<int>& group, const std::vector<int>& n,
                          const std::vector<int>& trials, const std::vector<double>& p0,
                          const std::string& estimator, const std::string& variance);

}  // namespace spicyglm

// spicyR Cell statistics: the excess with a frailty GEE (Supplementary Part I), the general design
// (covariates, availability adjustment, survival score test), Cox models, and the combination of
// several radii. Plain C++17 with Eigen, shared by the R and Python packages.
//
// Every routine has a plain-R reference written from the equations (tests/testthat/helper-reference.R in
// spicyR); the package tests check the two agree.
// Written with AI assistance (Claude, Anthropic), directed by the authors; see NEWS.
#pragma once

#include <string>
#include <vector>

#include "spicyglm/core.hpp"

namespace spicyglm {

// ----------------------------------------------------------- distributions ---

double pnorm_upper(double z);                 // P(Z > z), accurate in the far tail
double norm_quantile(double p);                       // Wichura AS241
double pt_upper(double t, double df);         // P(T > t), Student t (df = inf: normal)
double pt_two_sided(double t, double df);     // P(|T| > |t|)

// ------------------------------------------------------------ image rows ---

// One row per image that carries information for a pair: O, E = p L, n and the working variance
// v = sigma^2 (L2 - L^2 / M), times the label-clustering factor psi (Supplementary, Sections 2-3).
struct ImageRows {
  std::vector<int> image;           // image index
  std::vector<double> O, E, n, v;
};

// `totals` and `out_sq_totals` are Dataset::pair_neighbour_totals / pair_neighbour_out_sq_totals,
// `counts` is image-major (img * n_types + type), psi is type-major (type * n_images + img), or
// empty for psi = 1. knn: the k-nearest-neighbour graph (no self-loops), else the radius graph.
ImageRows excess_image_rows(const std::vector<double>& totals, const std::vector<double>& out_sq_totals,
                            const std::vector<double>& counts, int n_types, int n_images, int from, int to,
                            bool knn, const std::vector<double>& psi);

// The allocation effect: the fraction of TARGET cells moved next to (or away from) REF cells. `any_totals` is
// Dataset::pair_neighbour_any_totals, `self_expected` Dataset::self_any_expected. Every candidate (a non-REF cell;
// for a self-pair every cell) has the binary score g_b = 1{b has a REF cell among its neighbours}; O = sum of g over
// the TARGET cells, E = p L with L = sum of g over the candidates (self-pair: the exact expectation `self_expected`),
// q = E / n_to the null share of TARGET cells next to REF, and v = p (1 - p) M / (M - 1) (L - L^2 / M) psi (g^2 = g).
// The side is chosen once per pair from the sum of O - E over all images (no group labels):
//   attraction (sum >= 0): n = n_to - E, so (O - E) / n = (share - q) / (1 - q); a fraction f of TARGET cells moved
//     next to REF cells gives f;
//   avoidance (sum < 0): n = E, so (O - E) / n = (share - q) / q; a fraction f of the TARGET cells with a REF
//     neighbour moved away gives -f;
// in both cases however many REF cells there are and however densely they are packed.
ImageRows allocation_image_rows(const std::vector<double>& any_totals, const std::vector<double>& self_expected,
                                const std::vector<double>& counts, int n_types, int n_images, int from, int to,
                                const std::vector<double>& psi);

// The two above on the caller's arrays, without copies (psi: nullptr for psi = 1).
ImageRows excess_image_rows(const double* totals, const double* out_sq_totals, const double* counts, int n_types,
                            int n_images, int from, int to, bool knn, const double* psi);
ImageRows allocation_image_rows(const double* any_totals, const double* self_expected, const double* counts,
                                int n_types, int n_images, int from, int to, const double* psi);

// Label-clustering factor (Supplementary, "Label Clustering"): per (TARGET type, image) the median
// over the requested REF types of HAC / (sigma^2 R), floored at 1. Returns type-major n_types x
// n_images. Needs the radius index built at the HAC bandwidth's scale (radius graph: r).
// allocation: the HAC of the binary scores of allocation_image_rows (designs 6, 7), with the sparsity
// rule on the number of TARGET cells next to REF.
std::vector<double> label_clustering_factor(const Dataset& data, const std::vector<int>& from,
                                            const std::vector<int>& to, const std::vector<double>& counts,
                                            int n_types, int n_images, bool knn, double h,
                                            bool allocation = false);

// The Kontextual excess (Dataset::kontextual_sums, 7 values per image), on Statial's scale (O, L times
// n_context / n_to, Q times its square): image rows O, E = p L,
// v = p (1 - p) M / (M - 1) (Q - L^2 / M) psi, n = D (the REF cells' context exposure), p = n_to / M
// (self-pair: (n_to - 1) / (M - 1)). Images without TARGET cells, exposure or variance give no row.
ImageRows kontextual_image_rows(const std::vector<double>& sums, const std::vector<double>& counts, int n_types,
                                int n_images, int from, int to, const std::vector<double>& psi);

// Label-clustering factor of the Kontextual excess: as label_clustering_factor, with the HAC of the
// Kontextual scores over the context candidates (Dataset::hac_phi_sums design 0; build_context first) and
// the sparsity rule on the unweighted pair counts `raw` (pair-major: k * n_images + img).
std::vector<double> kontextual_clustering_factor(const Dataset& data, const std::vector<int>& from,
                                                 const std::vector<int>& to, const std::vector<double>& raw,
                                                 const std::vector<double>& counts, int n_types, int n_images,
                                                 double h);

// Images within a unit (Supplementary Methods, Section 1): a second random effect, for the spread of the effect
// between a unit's images beyond their sampling variance. Each row becomes its effect y = (O - E) / n with
// O = y, E = 0, n = 1 and v = v / n^2 + tau2_image, so an image's weight 1 / (v / n^2 + tau2_image) is at most
// 1 / tau2_image however many cells it has. tau2_image is the DerSimonian-Laird estimate from the spread within
// units (weights n^2 / v; units with one row do not contribute; 0 when no unit has two rows). `unit` is per row.
// With one row per unit the rows are an exact rescaling of the input and every test gives the same result.
ImageRows image_heterogeneity_rows(const ImageRows& rows, const std::vector<int>& unit, double* tau2_image = nullptr);

// ------------------------------------------------------ the two-group test ---

enum class Variance { CR2, HartungKnapp };

struct ExcessResult {
  bool ok = false;
  std::string reason;               // why a pair was not tested
  double coef_ref = 0, coef_comp = 0, difference = 0, se = 0, df = 0, p = 1, tau2 = 0;
  // per unit (in unit-code order, units absent from the pair have a 0 influence):
  std::vector<double> influence;    // CR2-adjusted influence on the difference (sum of squares = se^2)
  std::vector<double> unit_summary; // the patient summary delta~_i
  std::vector<double> unit_info;    // I_i
  std::vector<double> image_weight; // per image row: its share of its group's information (sums to 1 in each group)
};

// unit: unit code of every image (0 .. n_units - 1); group: 0 or 1 per image.
ExcessResult excess_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                         int n_units, bool frailty, Variance variance);

// Paule-Mandel tau^2 of the two-group model (exposed for the tests).
double excess_tau2(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                   int n_units);

// ------------------------------------------------------- the general design ---

// delta_ij = z_ij' theta + b_i (Supplementary Methods, Section 1). Z is row-major, one row per image row,
// q columns. tau2 < 0: estimate it by Paule-Mandel under the design (exact for patient-level
// designs); tau2 >= 0: hold it fixed. Tests the contrast c' theta with CR2 on Satterthwaite df.
struct DesignResult {
  bool ok = false;
  std::string reason;
  std::vector<double> theta;
  double estimate = 0, se = 0, df = 0, p = 1, tau2 = 0;
  std::vector<double> influence;  // per unit (unit-code order): CR2 contribution to the contrast; sum of squares = se^2
};

DesignResult design_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                         const std::vector<double>& Z, int q, const std::vector<double>& contrast, double tau2);
// Several contrasts of one fit: `contrasts` is row-major k x q. All results share theta and tau2.
std::vector<DesignResult> design_tests(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                                       const std::vector<double>& Z, int q, const std::vector<double>& contrasts,
                                       int k, double tau2, bool hartung_knapp = false);

// Option 2 (Supplementary Methods, Section 2): the two-group test with a slope on the centred covariate x
// (the log share of the REF type), tau2 held at the unadjusted value.
DesignResult availability_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                               int n_units, const std::vector<double>& x, double tau2);

// ------------------------------------------------------------------ Cox ---

// Cox proportional hazards with Efron ties, Newton-Raphson with step halving. X is row-major n x p
// (p may be 0: the null model). Martingale residuals are returned for every subject.
struct CoxResult {
  bool ok = false;
  std::vector<double> beta, se, p;  // Wald
  std::vector<double> martingale;
  double loglik = 0;
};

CoxResult cox_fit(const std::vector<double>& time, const std::vector<int>& event, const std::vector<double>& X,
                  int p);

// Survival (Supplementary Methods, Section 3) for one pair. time and event per unit (unit-code order);
// covariates (row-major, per unit) enter the null Cox model.
struct SurvivalResult {
  bool ok = false;
  std::string reason;
  double score_coef = 0, score_se = 0, score_df = 0, score_p = 1;  // theta_1 of delta ~ 1 + M
  double log_hr_sd = 0, hr_sd = 1, hr_se = 0, hr_p = 1;             // Cox on the standardised BLUP
  double log_hr_unit = 0;   // per unit of the effect (count: per extra cell per target cell; allocation: per 100% of target cells)
  double tau2 = 0;
};

// x: optional per-row covariate (the log share of the `to` type); when given, the score test and the hazard
// ratio are adjusted for it (empty: unadjusted).
SurvivalResult survival_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                             const std::vector<double>& martingale, const std::vector<double>& time,
                             const std::vector<int>& event, const std::vector<double>& x = {});

// ------------------------------------------------------ several radii ---

// Cauchy combination of p-values (equal weights).
double cauchy_combine(const std::vector<double>& p);

// max-T over radii with the sandwich correlation of the CR2 influences (PLAN_multiradius_maxT.md).
// influence: K vectors (one per radius) over the same units; t and df per radius. Returns the p
// value and the index of the radius with the largest |z|.
struct MaxTResult {
  double p = 1;
  int best = 0;
};
MaxTResult max_t(const std::vector<std::vector<double>>& influence, const std::vector<double>& t,
                 const std::vector<double>& df);

// Davies (1987) upper bound for the maximum of a smooth Gaussian process over the radii, used as the p value of the
// max-T for many radii: p = 2 Phi(-c) + V exp(-c^2 / 2) / pi, with c the largest |z| and V = sum_k arccos(rho_{k,k+1})
// the length of the standardised process (correlations of adjacent radii, in the given order, from the influences).
// Same inputs and result as max_t; any number of radii.
MaxTResult davies_max(const std::vector<std::vector<double>>& influence, const std::vector<double>& t,
                      const std::vector<double>& df);

// P(Z outside the rectangle [lower, upper]) for Z ~ N(0, R), R a correlation matrix (row-major K x K), as max_t
// uses it: for K <= 2 the quadrature below (accurate to about 1e-10); for K >= 3 Genz-Bretz quasi-Monte Carlo
// (mvn_outside_qmc) with an absolute error target of 1e-4.
double mvn_outside(const std::vector<double>& lower, const std::vector<double>& upper,
                   const std::vector<double>& R, int K);

// Genz (1992) and Genz & Bretz (2009, Computation of Multivariate Normal and t Probabilities, Springer, Chapter 4)
// quasi-Monte Carlo for the same probability (K <= 6): separation of variables on the Cholesky factor, with the
// variables reordered so the narrowest expected interval comes first, integrated over the unit cube by Korobov
// lattice rules of growing size with the baker's transformation and antithetic points, each rule under 12 random
// shifts. The shifts come from a fixed-seed generator, so the result is the same on every call. The rules are
// combined by inverse-variance weights until 3.5 standard errors are below `abseps`, or at most `max_points` points
// (per shift, over all rules). Coordinates with conditional variance below 1e-10 (a singular R) are deterministic
// given the others and enter as indicators.
struct MvnResult {
  double p = 0;
  double error = 0;   // 3.5 standard errors of p over the shifts
  long points = 0;    // lattice points used, per shift
};
MvnResult mvn_outside_qmc(const std::vector<double>& lower, const std::vector<double>& upper,
                          const std::vector<double>& R, int K, double abseps = 1e-4, long max_points = 100000);

// The quadrature max_t used up to spicyR 1.99.24, kept for comparison: recursive conditioning with Gauss-Legendre
// quadrature (about 2e6 evaluations; K <= 6), relative accuracy kept in the tail.
double mvn_outside_quadrature(const std::vector<double>& lower, const std::vector<double>& upper,
                              const std::vector<double>& R, int K);

// ------------------------------------------- every pair at every radius ---

// spicy()'s two-condition test of every pair at every radius (or at the one k) in one call, with no per-test trip
// back to the front end: per pair and radius the image rows (allocation_image_rows or excess_image_rows), the side,
// the image random effect (image_heterogeneity_rows), excess_test and, when there is something to adjust for, the
// design of the adjusted test (pair_design_test); then, with several radii, their combination (max_t, davies_max
// or cauchy_combine) for the main and the unadjusted test. The numbers are exactly those of the front end calling
// these functions one test at a time (R/cell.R of spicyR 1.99.23), whatever n_threads.
enum class Combine { None, MaxT, Davies, Cauchy };

struct PairTestInput {
  // per radius (or the one k): allocation: a = any_totals, b = self_expected; count: a = totals, b = out_sq_totals
  // (the layouts of allocation_image_rows / excess_image_rows)
  std::vector<const double*> a, b;
  bool allocation = true, knn = false;
  const double* counts = nullptr;     // image-major n_images x n_types
  int n_types = 0, n_images = 0;
  const double* psi = nullptr;        // type-major n_types x n_images, or nullptr for psi = 1 (every radius)
  std::vector<int> from, to;          // the pairs (0-based types)
  const int* image_unit = nullptr;    // per image: its unit (patient) code
  const int* image_group = nullptr;   // per image: its condition, 0 or 1, or negative when missing
  int n_units = 0;
  bool image_heterogeneity = true, frailty = true;
  Variance variance = Variance::CR2;
  // the adjusted test: the abundance (log share of the `from` type) and the covariates z (image-major
  // n_images x n_z, already centred over the images; NaN where missing)
  bool adjust = false;
  const double* z = nullptr;
  int n_z = 0;
  int keep_radius = 0;                // the radius whose image rows and image weights are returned
  Combine combine = Combine::None;    // with several radii
  int n_threads = 1;                  // pairs handed out in turn to up to n_threads threads
};

// The adjusted test of one pair (R/cell.R .cell_design, .cell_design_tau2 and .cell_contrasts of spicyR 1.99.23, two
// conditions): rows (after the image random effect) with a group of 0 or 1 (negative: missing; left out) and
// complete covariates (NaN: missing), the design of the two group indicators, the centred log share of the `from`
// type (adjust; when not constant) and the centred covariates; tau2 re-estimated (-1) when every added column is
// constant within units, else held at `tau2` (frailty; 0 without). design = 0: nothing to adjust for (no column
// added); 1: adjusted, `tests` = the condition contrast, then one per added column (abundance first); 2: the design
// failed (`tests[0].reason`). fit = false: the design is not fitted (the response has no variation); design 2 with
// reason "no_variation" when there is something to adjust for, else 0.
struct PairDesign {
  int design = 0;
  bool abundance = false;             // the abundance column entered the design
  std::vector<DesignResult> tests;
};
PairDesign pair_design_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                            int n_units, const double* counts, int n_types, int from, bool adjust, const double* z,
                            int n_z, bool frailty, double tau2, bool hartung_knapp, bool fit = true);

// Results, per pair and radius at index p * K + k (K radii); the main test is the adjusted one when design == 1,
// else the unadjusted one. Where the unadjusted test has no variation (se 0: 0 / 0) no design is fitted: design 2,
// reason "no_variation", p NaN.
struct PairTests {
  int n_pairs = 0, n_radii = 0, n_z = 0;
  std::vector<int> ok;                // the unadjusted test was done
  std::vector<int> side;              // allocation: 1 attraction, -1 avoidance; count: 0
  std::vector<int> n_rows;            // image rows tested
  std::vector<double> tau2_image;     // NaN without the image random effect or rows
  std::vector<double> coef_ref, coef_comp, difference, se, df, p, tau2;   // the main test
  std::vector<double> u_difference, u_se, u_df, u_p;   // the unadjusted test
  std::vector<int> design;            // PairDesign::design
  std::vector<int> abundance;         // PairDesign::abundance
  std::vector<std::string> design_reason;   // design == 2: why
  // the added columns (abundance, then the n_z covariates) at ((p * K + k) * (1 + n_z) + j): estimate and p, NaN
  // where not in the design
  std::vector<double> effect, effect_p;
  // at keep_radius, per pair: the image rows tested and the unadjusted test's image weights (empty when not done)
  std::vector<ImageRows> rows;
  std::vector<std::vector<double>> image_weight;
  // several radii, per pair with at least one radius tested: the radius chosen (0-based; -1 for none), the
  // combined p of the main and of the unadjusted test, and the main test's p at the chosen radius
  // (comb_u_p only when there is something to adjust for, else NaN)
  std::vector<int> best;
  std::vector<double> comb_p, comb_p_best, comb_u_p;
};
PairTests pair_tests(const PairTestInput& in);

}  // namespace spicyglm

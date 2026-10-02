// spicyR Cell statistics: the excess with a frailty GEE (Supplementary Part I), the general design
// (covariates, availability adjustment, survival score test), Cox models, and the combination of
// several radii. Plain C++17 with Eigen, shared by the R and Python packages.
//
// Every routine has a plain-R reference written from the equations (R/reference/ in the project
// repository ProjectSpicyRCell); the package tests check the two agree.
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
double pchisq_upper(double x, double df);     // P(X > x)

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

// Label-clustering factor (Supplementary, "Label Clustering"): per (TARGET type, image) the median
// over the requested REF types of HAC / (sigma^2 R), floored at 1. Returns type-major n_types x
// n_images. Needs the radius index built at the HAC bandwidth's scale (radius graph: r).
std::vector<double> label_clustering_factor(const Dataset& data, const std::vector<int>& from,
                                            const std::vector<int>& to, const std::vector<double>& counts,
                                            int n_types, int n_images, bool knn, double h);

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
};

// unit: unit code of every image (0 .. n_units - 1); group: 0 or 1 per image.
ExcessResult excess_test(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                         int n_units, bool frailty, Variance variance);

// Paule-Mandel tau^2 of the two-group model (exposed for the tests).
double excess_tau2(const ImageRows& rows, const std::vector<int>& unit, const std::vector<int>& group,
                   int n_units);

// ------------------------------------------------------- the general design ---

// delta_ij = z_ij' theta + b_i (new_methods.pdf, Section 1). Z is row-major, one row per image row,
// q columns. tau2 < 0: estimate it by Paule-Mandel under the design (exact for patient-level
// designs); tau2 >= 0: hold it fixed. Tests the contrast c' theta with CR2 on Satterthwaite df.
struct DesignResult {
  bool ok = false;
  std::string reason;
  std::vector<double> theta;
  double estimate = 0, se = 0, df = 0, p = 1, tau2 = 0;
};

DesignResult design_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                         const std::vector<double>& Z, int q, const std::vector<double>& contrast, double tau2);

// Option 2 (new_methods.pdf, Section 2): the two-group test with a slope on the centred covariate x
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

// Survival (new_methods.pdf, Section 3) for one pair. time and event per unit (unit-code order);
// covariates (row-major, per unit) enter the null Cox model.
struct SurvivalResult {
  bool ok = false;
  std::string reason;
  double score_coef = 0, score_se = 0, score_df = 0, score_p = 1;  // theta_1 of delta ~ 1 + M
  double log_hr_sd = 0, hr_sd = 1, hr_se = 0, hr_p = 1;             // Cox on the standardised BLUP
  double log_hr_unit = 0;                                          // per extra cell per target cell
  double tau2 = 0;
};

SurvivalResult survival_test(const ImageRows& rows, const std::vector<int>& unit, int n_units,
                             const std::vector<double>& martingale, const std::vector<double>& time,
                             const std::vector<int>& event);

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

// P(Z outside the rectangle [lower, upper]) for Z ~ N(0, R), R a correlation matrix (row-major K x K):
// recursive conditioning with Gauss-Legendre quadrature, relative accuracy kept in the tail.
double mvn_outside(const std::vector<double>& lower, const std::vector<double>& upper,
                   const std::vector<double>& R, int K);

}  // namespace spicyglm

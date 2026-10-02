// Distribution functions for the statistics module (no dependency beyond the standard library). Written from the
// formulas of the NIST Digital Library of Mathematical Functions (DLMF, https://dlmf.nist.gov) and, for the normal
// quantile, Wichura (1988), Applied Statistics 37, 477-484 (Algorithm AS 241).
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "spicyglm/stats.hpp"

namespace spicyglm {

namespace {

const double kInf = std::numeric_limits<double>::infinity();

// The continued fraction of the regularised incomplete beta function (DLMF 8.17.22),
//   I_x(a, b) = x^a (1 - x)^b / (a B(a, b)) / (1 + d_1 / (1 + d_2 / (1 + ...))),
//   d_{2m+1} = -(a + m)(a + b + m) x / ((a + 2m)(a + 2m + 1)),   d_{2m} = m (b - m) x / ((a + 2m - 1)(a + 2m)),
// which converges fast for x < (a + 1) / (a + b + 2). Returns the value of the fraction 1 + d_1 / (1 + ...),
// computed from its convergents A_n / B_n by the fundamental recurrences A_n = A_{n-1} + d_n A_{n-2} (and the
// same for B_n; DLMF 1.12.5), rescaled to stay within floating-point range.
double incomplete_beta_fraction(double a, double b, double x) {
  double A_before = 1, B_before = 0;   // A_{-1}, B_{-1}
  double A_now = 1, B_now = 1;         // A_0, B_0
  double convergent = 1;
  for (int n = 1; n <= 20000; ++n) {
    const double m = static_cast<double>(n / 2);
    const double d = (n % 2 == 1) ? -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
                                  : m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
    const double A_next = A_now + d * A_before, B_next = B_now + d * B_before;
    A_before = A_now; B_before = B_now; A_now = A_next; B_now = B_next;
    const double scale = std::max(std::fabs(A_now), std::fabs(B_now));
    if (scale > 1e150 || (scale > 0 && scale < 1e-150)) {
      A_now /= scale; B_now /= scale; A_before /= scale; B_before /= scale;
    }
    const double next = A_now / B_now;
    if (n > 1 && std::fabs(next - convergent) <= 1e-15 * std::fabs(next)) return next;
    convergent = next;
  }
  return convergent;
}

// I_x(a, b), the regularised incomplete beta function, given x and y = 1 - x (passed separately so that a y too
// small to be represented as 1 - x keeps its precision); for x beyond (a + 1) / (a + b + 2) through the symmetry
// I_x(a, b) = 1 - I_y(b, a) (DLMF 8.17.4).
double ibeta(double a, double b, double x, double y) {
  if (x <= 0) return 0;
  if (y <= 0) return 1;
  const double log_front = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log(y);
  if (x < (a + 1) / (a + b + 2)) return std::exp(log_front) / (a * incomplete_beta_fraction(a, b, x));
  return 1 - std::exp(log_front) / (b * incomplete_beta_fraction(b, a, y));
}

}  // namespace

double pnorm_upper(double z) { return 0.5 * std::erfc(z / std::sqrt(2.0)); }

double norm_quantile(double p) {
  if (p <= 0) return -kInf;
  if (p >= 1) return kInf;
  // Wichura (1988), Algorithm AS 241, PPND16.
  double q = p - 0.5, r, val;
  if (std::fabs(q) <= 0.425) {
    r = 0.180625 - q * q;
    val = q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r +
                  45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
                133.14166789178437745) * r + 3.387132872796366608) /
          (((((((r * 5226.495278852545925 + 28729.085735721942674) * r + 39307.89580009271061) * r +
               21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
            42.313330701600911252) * r + 1.0);
    return val;
  }
  r = q < 0 ? p : 1 - p;
  r = std::sqrt(-std::log(r));
  if (r <= 5) {
    r -= 1.6;
    val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r + 0.24178072517745061177) * r +
               1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r +
             4.6303378461565452959) * r + 1.42343711074968357734) /
          (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r +
               0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r +
            2.05319162663775882187) * r + 1.0);
  } else {
    r -= 5;
    val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r +
               0.026532189526576123093) * r + 0.29656057182850489123) * r + 1.7848265399172913358) * r +
             5.4637849111641143699) * r + 6.6579046435011037772) /
          (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
               7.868691311456132591e-4) * r + 0.0148753612908506148525) * r + 0.13692988092273580531) * r +
            0.59983220655588793769) * r + 1.0);
  }
  return q < 0 ? -val : val;
}

double pt_upper(double t, double df) {
  if (std::isnan(t) || std::isnan(df) || df <= 0) return std::numeric_limits<double>::quiet_NaN();
  if (std::isinf(df) || df > 1e10) return pnorm_upper(t);
  const double x = df / (df + t * t), y = t * t / (df + t * t);
  double tail = 0.5 * ibeta(df / 2, 0.5, x, y);   // P(T > |t|)
  return t >= 0 ? tail : 1 - tail;
}

double pt_two_sided(double t, double df) {
  if (std::isnan(t) || std::isnan(df) || df <= 0) return std::numeric_limits<double>::quiet_NaN();
  if (std::isinf(df) || df > 1e10) return std::erfc(std::fabs(t) / std::sqrt(2.0));
  return ibeta(df / 2, 0.5, df / (df + t * t), t * t / (df + t * t));
}


}  // namespace spicyglm

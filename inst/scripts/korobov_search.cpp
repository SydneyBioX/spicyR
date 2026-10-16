// The Korobov lattice rules of the max-T's quasi-Monte Carlo (src/core/src/multiradius.cpp, kLattice).
//
// For each prime n in the sequence 64, 96, 144, ... (each the next prime at or above 1.5 times the previous
// target), an exhaustive search over a in [2, n/2] for the rank-1 lattice with generating vector
// (1, a, a^2, a^3, a^4) mod n that minimises the P2 criterion (the worst-case error in the weighted Korobov space
// of smoothness 1) with product weights gamma_j = 0.5^j, j = 0..4:
//   P2(a) = -1 + (1/n) sum_k prod_j (1 + gamma_j 2 pi^2 B2({k z_j / n})),   B2(x) = x^2 - x + 1/6.
// The weights favour the first coordinates, which carry most of the variance after the variable reordering.
//   g++ -O2 korobov_search.cpp -o korobov_search && ./korobov_search 0.5      (about a minute)
// prints "n a P2" per rule.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

static bool is_prime(long n) {
  if (n < 2) return false;
  for (long d = 2; d * d <= n; ++d) if (n % d == 0) return false;
  return true;
}

int main(int argc, char** argv) {
  const int D = 5;
  const double gamma = argc > 1 ? std::atof(argv[1]) : 0.5, c = 2 * M_PI * M_PI;
  double g[D];
  for (int j = 0; j < D; ++j) g[j] = std::pow(gamma, j);
  std::vector<long> ns;
  long target = 64, total = 0;
  while (total + target <= 200000) {
    long n = target;
    while (!is_prime(n)) ++n;
    ns.push_back(n); total += n; target = target * 3 / 2;
  }
  for (long n : ns) {
    std::vector<double> B(n);
    for (long k = 0; k < n; ++k) { double x = static_cast<double>(k) / n; B[k] = c * (x * x - x + 1.0 / 6); }
    double best = 1e300; long best_a = 1;
    for (long a = 2; a <= n / 2; ++a) {
      long z[D]; z[0] = 1;
      for (int j = 1; j < D; ++j) z[j] = (z[j - 1] * a) % n;
      double s = 0;
      for (long k = 0; k < n; ++k) {
        double p = 1;
        for (int j = 0; j < D; ++j) p *= 1 + g[j] * B[(k * z[j]) % n];
        s += p;
        if (s / n - 1 > best) break;   // already worse
      }
      const double P = s / n - 1;
      if (P < best) { best = P; best_a = a; }
    }
    std::printf("%ld %ld %.6g\n", n, best_a, best);
    std::fflush(stdout);
  }
}

#include "spicyglm/core.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace spicyglm {

namespace {

double cross(double ox, double oy, double ax, double ay, double bx, double by) {
  return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox);
}

// Andrew's monotone chain: indices of the hull's vertices, anticlockwise,
// collinear points dropped.
std::vector<std::size_t> convex_hull(const double* x, const double* y, std::size_t n) {
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  if (n < 3) return idx;
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
    return x[a] < x[b] || (x[a] == x[b] && y[a] < y[b]);
  });
  std::vector<std::size_t> hull(2 * n);
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t p = idx[i];
    while (k >= 2 && cross(x[hull[k - 2]], y[hull[k - 2]], x[hull[k - 1]], y[hull[k - 1]], x[p], y[p]) <= 0) --k;
    hull[k++] = p;
  }
  for (std::size_t i = n - 1, lower = k + 1; i-- > 0;) {
    std::size_t p = idx[i];
    while (k >= lower && cross(x[hull[k - 2]], y[hull[k - 2]], x[hull[k - 1]], y[hull[k - 1]], x[p], y[p]) <= 0) --k;
    hull[k++] = p;
  }
  hull.resize(k - 1);
  return hull;
}

// The shoelace formula over the convex hull.
double convex_hull_area(const double* x, const double* y, std::size_t n) {
  if (n < 3) return 0.0;
  std::vector<std::size_t> hull = convex_hull(x, y, n);
  double twice = 0.0;
  for (std::size_t i = 0; i < hull.size(); ++i) {
    std::size_t a = hull[i], b = hull[(i + 1) % hull.size()];
    twice += x[a] * y[b] - x[b] * y[a];
  }
  return std::abs(twice) / 2.0;
}

// Signed area of disc(0, R) intersected with the triangle (0, a, b). The segment
// a -> b is split where it crosses the circle; a piece inside the disc adds the
// triangle it spans with the origin, a piece outside adds the circular sector.
double disc_triangle_area(double ax, double ay, double bx, double by, double R) {
  double dx = bx - ax, dy = by - ay;
  double A = dx * dx + dy * dy;
  if (!(A > 0)) return 0.0;
  double B = ax * dx + ay * dy, C = ax * ax + ay * ay - R * R;
  double root = std::sqrt(std::max(B * B - A * C, 0.0));
  double t1 = std::min(std::max((-B - root) / A, 0.0), 1.0);
  double t2 = std::min(std::max((-B + root) / A, 0.0), 1.0);
  auto piece = [&](double s, double e) {
    double px = ax + s * dx, py = ay + s * dy, qx = ax + e * dx, qy = ay + e * dy;
    double cr = px * qy - py * qx;
    double mx = ax + (s + e) / 2 * dx, my = ay + (s + e) / 2 * dy;
    return mx * mx + my * my <= R * R ? cr / 2 : R * R * std::atan2(cr, px * qx + py * qy) / 2;
  };
  return piece(0.0, t1) + piece(t1, t2) + piece(t2, 1.0);
}

// Area of the disc of radius R about (x, y) inside a polygon, exactly: the
// signed areas of disc-triangle intersections summed over the polygon's edges.
double disc_window_area(const std::vector<double>& vx, const std::vector<double>& vy,
                        double x, double y, double R) {
  double area = 0.0;
  for (std::size_t k = 0, m = vx.size(); k < m; ++k) {
    std::size_t k2 = (k + 1) % m;
    area += disc_triangle_area(vx[k] - x, vy[k] - y, vx[k2] - x, vy[k2] - y, R);
  }
  return std::abs(area);
}

// |W intersect (W + (dx, dy))| for a convex W, by clipping W + v to W's edge
// half-planes (Sutherland-Hodgman). The fallback of translation_overlap.
double clipped_overlap(const WindowGeometry& W, double dx, double dy) {
  std::size_t m = W.px.size();
  std::vector<double> sx(m), sy(m), ox, oy;
  for (std::size_t k = 0; k < m; ++k) {
    sx[k] = W.px[k] + dx;
    sy[k] = W.py[k] + dy;
  }
  for (std::size_t k = 0; k < m && !sx.empty(); ++k) {
    double c = W.nx[k] * W.px[k] + W.ny[k] * W.py[k];
    ox.clear();
    oy.clear();
    for (std::size_t i = 0, s = sx.size(); i < s; ++i) {
      std::size_t j = (i + 1) % s;
      double fi = W.nx[k] * sx[i] + W.ny[k] * sy[i] - c, fj = W.nx[k] * sx[j] + W.ny[k] * sy[j] - c;
      if (fi <= 0) {
        ox.push_back(sx[i]);
        oy.push_back(sy[i]);
      }
      if ((fi <= 0) != (fj <= 0)) {
        double t = fi / (fi - fj);
        ox.push_back(sx[i] + t * (sx[j] - sx[i]));
        oy.push_back(sy[i] + t * (sy[j] - sy[i]));
      }
    }
    sx.swap(ox);
    sy.swap(oy);
  }
  double twice = 0.0;
  for (std::size_t i = 0, s = sx.size(); i < s; ++i) {
    std::size_t j = (i + 1) % s;
    twice += sx[i] * sy[j] - sy[i] * sx[j];
  }
  return std::abs(twice) / 2.0;
}

// |W intersect (W + (dx, dy))|, exactly (Section 14.10). W + v has the same edge
// normals as W, so for a convex W the intersection is W with every edge whose
// outward normal n_k has n_k . v < 0 moved inward by -n_k . v, and its vertices
// are the intersections of consecutive moved edge lines. If a move makes an
// edge vanish (its end vertices cross over), the intersection is clipped
// directly instead. h, qx and qy are scratch space.
double translation_overlap(const WindowGeometry& W, double dx, double dy, std::vector<double>& h,
                           std::vector<double>& qx, std::vector<double>& qy) {
  if (W.rectangle)
    return std::max(W.vx[1] - W.vx[0] - std::abs(dx), 0.0) * std::max(W.vy[2] - W.vy[1] - std::abs(dy), 0.0);
  std::size_t m = W.px.size();
  h.resize(m);
  qx.resize(m);
  qy.resize(m);
  for (std::size_t k = 0; k < m; ++k)
    h[k] = W.nx[k] * W.px[k] + W.ny[k] * W.py[k] + std::min(0.0, dx * W.nx[k] + dy * W.ny[k]);
  for (std::size_t k = 0; k < m; ++k) {  // vertex k joins edge k and edge k + 1
    std::size_t k2 = (k + 1) % m;
    qx[k] = (h[k] * W.ny[k2] - h[k2] * W.ny[k]) / W.det[k];
    qy[k] = (h[k2] * W.nx[k] - h[k] * W.nx[k2]) / W.det[k];
  }
  bool kept = true;
  double twice = 0.0;
  for (std::size_t k = 0; k < m; ++k) {  // edge k runs from vertex k - 1 to vertex k
    std::size_t k0 = (k + m - 1) % m;
    kept = kept && (qx[k] - qx[k0]) * W.ex[k] + (qy[k] - qy[k0]) * W.ey[k] >= 0;
    twice += qx[k0] * qy[k] - qy[k0] * qx[k];
  }
  return kept ? twice / 2 : clipped_overlap(W, dx, dy);
}

// Translation edge correction |W| / |W intersect (W + v)|, capped at spatstat's
// default maxedgewt of 100.
double translation_weight(const WindowGeometry& W, double dx, double dy, std::vector<double>& h,
                          std::vector<double>& qx, std::vector<double>& qy) {
  return std::min(W.area / translation_overlap(W, dx, dy, h, qx, qy), 100.0);
}

// The window of one image's cells: its bounding rectangle or convex hull.
WindowGeometry window_geometry(const double* x, const double* y, std::size_t n, Window window) {
  WindowGeometry W;
  W.area = window_area(x, y, n, window);
  if (n == 0) return W;
  if (window == Window::Rectangle) {
    auto [xmin, xmax] = std::minmax_element(x, x + n);
    auto [ymin, ymax] = std::minmax_element(y, y + n);
    W.vx = {*xmin, *xmax, *xmax, *xmin};
    W.vy = {*ymin, *ymin, *ymax, *ymax};
    return W;
  }
  W.rectangle = false;
  for (std::size_t i : convex_hull(x, y, n)) {
    W.vx.push_back(x[i]);
    W.vy.push_back(y[i]);
  }
  std::size_t m = W.vx.size();
  double cx = std::accumulate(W.vx.begin(), W.vx.end(), 0.0) / m;
  double cy = std::accumulate(W.vy.begin(), W.vy.end(), 0.0) / m;
  for (std::size_t k = 0; k < m; ++k) {  // centred, for accuracy
    W.px.push_back(W.vx[k] - cx);
    W.py.push_back(W.vy[k] - cy);
  }
  for (std::size_t k = 0; k < m; ++k) {
    std::size_t k2 = (k + 1) % m;
    double ex = W.px[k2] - W.px[k], ey = W.py[k2] - W.py[k], len = std::sqrt(ex * ex + ey * ey);
    W.ex.push_back(ex);
    W.ey.push_back(ey);
    W.nx.push_back(ey / len);  // outward normal of edge k, anticlockwise polygon
    W.ny.push_back(-ex / len);
  }
  for (std::size_t k = 0; k < m; ++k) {
    std::size_t k2 = (k + 1) % m;
    W.det.push_back(W.nx[k] * W.ny[k2] - W.ny[k] * W.nx[k2]);
  }
  return W;
}


// Exact random-labelling variance of a self-pair a -> a (Supplementary Section 3.6, Appendix B), for any neighbour
// relation: N(l) = the out-neighbours of cell l (the cells within r, or its k nearest; l itself excluded).
// Built once per image from the neighbour lists; evaluated for every type count n with self_graph_variance().
struct SelfGraph {
  int N = 0;
  // count: m links, b of them with their reverse a link, w in-degree, t = out + in degree
  long double m = 0, b = 0, S = 0, Swt = 0, Sw2 = 0;
  // fraction: number of cells with each out-degree d, and the coefficient of h2(x) in
  //   sum over ordered pairs l != k of T_lk (Proposition "self-frac"), apart from the constant N (N - 1)
  std::vector<double> deg_count, h2_coef;
};

// ptr / nb: out-neighbour lists of the image's cells, local indices 0..N-1 (CSR).
SelfGraph self_graph(int N, const std::vector<int>& ptr, const std::vector<int>& nb, bool allocation) {
  SelfGraph g; g.N = N;
  std::vector<int> d(N), w(N, 0);
  int dmax = 0;
  for (int l = 0; l < N; ++l) { d[l] = ptr[l + 1] - ptr[l]; dmax = std::max(dmax, d[l]); }
  for (int e = 0; e < ptr[N]; ++e) ++w[nb[e]];
  // in-neighbour lists
  std::vector<int> iptr(N + 1, 0), inb(ptr[N]);
  for (int e = 0; e < ptr[N]; ++e) ++iptr[nb[e] + 1];
  for (int l = 0; l < N; ++l) iptr[l + 1] += iptr[l];
  { std::vector<int> pos(iptr.begin(), iptr.end() - 1);
    for (int l = 0; l < N; ++l) for (int e = ptr[l]; e < ptr[l + 1]; ++e) inb[pos[nb[e]]++] = l; }
  std::vector<char> out_mark(N, 0), in_mark(N, 0);
  if (!allocation) {
    for (int l = 0; l < N; ++l) {
      for (int e = ptr[l]; e < ptr[l + 1]; ++e) out_mark[nb[e]] = 1;
      for (int e = iptr[l]; e < iptr[l + 1]; ++e) if (out_mark[inb[e]]) g.b += 1;   // l -> k and k -> l
      for (int e = ptr[l]; e < ptr[l + 1]; ++e) out_mark[nb[e]] = 0;
      const long double t = d[l] + w[l];
      g.m += d[l]; g.S += t * (t - 1); g.Swt += w[l] * t; g.Sw2 += static_cast<long double>(w[l]) * w[l];
    }
    return g;
  }
  g.deg_count.assign(dmax + 1, 0.0);
  for (int l = 0; l < N; ++l) g.deg_count[d[l]] += 1;
  std::vector<double>& A = g.h2_coef;
  A.assign(2 * dmax + 1, 0.0);
  // every ordered pair as if unlinked and sharing no neighbour: 1 - h2(d_l) - h2(d_k) + h2(d_l + d_k)
  for (int x = 0; x <= dmax; ++x) {
    if (!g.deg_count[x]) continue;
    A[x] -= 2.0 * (N - 1) * g.deg_count[x];
    A[2 * x] -= g.deg_count[x];                                   // the l = k terms of the double sum below
    for (int y = 0; y <= dmax; ++y) if (g.deg_count[y]) A[x + y] += g.deg_count[x] * g.deg_count[y];
  }
  // corrections for the pairs that are linked or share a neighbour: T_lk minus the term above. T_lk = T_kl, so each
  // unordered pair is visited once (k > l) and counted twice; each in-neighbour list is in increasing order (built by
  // l above), so the k > l part of it starts at upper_bound(l).
  std::vector<int> common(N, 0), touched;
  for (int l = 0; l < N; ++l) {
    touched.clear();
    auto touch = [&](int k) { if (!out_mark[k] && !in_mark[k] && !common[k]) touched.push_back(k); };
    for (int e = ptr[l]; e < ptr[l + 1]; ++e) if (nb[e] > l) { touch(nb[e]); out_mark[nb[e]] = 1; }
    for (int e = iptr[l]; e < iptr[l + 1]; ++e) if (inb[e] > l) { touch(inb[e]); in_mark[inb[e]] = 1; }
    for (int e = ptr[l]; e < ptr[l + 1]; ++e) {
      const int v = nb[e];
      const int* first = std::upper_bound(inb.data() + iptr[v], inb.data() + iptr[v + 1], l);
      for (const int* f = first; f < inb.data() + iptr[v + 1]; ++f) { touch(*f); ++common[*f]; }
    }
    for (int k : touched) {
      const int s = d[l] + d[k];
      if (out_mark[k] && in_mark[k]) { A[d[l]] += 2; A[d[k]] += 2; A[s] -= 2; }      // linked both ways: T = 1
      else if (out_mark[k]) { A[d[l]] += 2; A[s] -= 2; }                             // k in N(l) only
      else if (in_mark[k]) { A[d[k]] += 2; A[s] -= 2; }                              // l in N(k) only
      else { A[s - common[k]] += 2; A[s] -= 2; }                                     // c common neighbours
    }
    for (int k : touched) { out_mark[k] = 0; in_mark[k] = 0; common[k] = 0; }
  }
  return g;
}

// Var(O) (allocation) or Var(O - p L) (count) of the self-pair of a type with n of the image's N cells.
double self_graph_variance(const SelfGraph& g, int n, bool allocation) {
  const int N = g.N;
  if (n < 2 || n >= N) return 0;
  auto lchoose = [](double a, double c) { return std::lgamma(a + 1) - std::lgamma(c + 1) - std::lgamma(a - c + 1); };
  long double pk[5] = {1, 0, 0, 0, 0};
  for (int k = 1; k <= 4; ++k) pk[k] = n >= k ? pk[k - 1] * (n - k + 1) / static_cast<long double>(N - k + 1) : 0;
  if (!allocation) {
    const long double p = static_cast<long double>(n - 1) / (N - 1), m = g.m, b = g.b;
    const long double EO2 = (m + b) * pk[2] + (g.S - 2 * b) * pk[3] + (m * m - m + b - g.S) * pk[4];
    const long double EOL = pk[2] * g.Swt + pk[3] * (m * m - g.Swt);
    const long double EL2 = pk[1] * g.Sw2 + pk[2] * (m * m - g.Sw2);
    const long double V = EO2 - 2 * p * EOL + p * p * EL2;
    return V > 1e-11L * (1 + EO2) ? static_cast<double>(V) : 0.0;
  }
  // h1(s) = C(N-1-s, n-1) / C(N-1, n-1); h2(s) = C(N-2-s, n-2) / C(N-2, n-2)
  const double l1 = lchoose(N - 1, n - 1), l2 = lchoose(N - 2, n - 2);
  long double E = 0;
  for (std::size_t x = 0; x < g.deg_count.size(); ++x) if (g.deg_count[x])
    E += g.deg_count[x] * (N - 1 - static_cast<double>(x) >= n - 1 ? -std::expm1(lchoose(N - 1 - x, n - 1) - l1) : 1.0);
  E *= static_cast<long double>(n) / N;
  long double sum = static_cast<long double>(N) * (N - 1);
  for (std::size_t x = 0; x < g.h2_coef.size(); ++x) if (g.h2_coef[x] != 0) {
    const double h2 = N - 2 - static_cast<double>(x) >= n - 2 ? std::exp(lchoose(N - 2 - x, n - 2) - l2) : 0.0;
    sum += g.h2_coef[x] * static_cast<long double>(h2);
  }
  const long double V = E - E * E + pk[2] * sum;
  return V > 1e-11L * (1 + E * E) ? static_cast<double>(V) : 0.0;
}

}  // namespace

double window_area(const double* x, const double* y, std::size_t n, Window window) {
  if (n == 0) return 0.0;
  if (window == Window::Convex) return convex_hull_area(x, y, n);
  auto [xmin, xmax] = std::minmax_element(x, x + n);
  auto [ymin, ymax] = std::minmax_element(y, y + n);
  return (*xmax - *xmin) * (*ymax - *ymin);
}

Dataset::Dataset(std::vector<double> x, std::vector<double> y, std::vector<int> cell_type,
                 std::vector<int> image_offsets, int n_types)
    : x_(std::move(x)), y_(std::move(y)), type_(std::move(cell_type)),
      image_offsets_(std::move(image_offsets)), n_types_(n_types) {
  std::size_t N = x_.size();
  if (y_.size() != N || type_.size() != N || image_offsets_.empty() ||
      image_offsets_.front() != 0 || static_cast<std::size_t>(image_offsets_.back()) != N)
    throw std::invalid_argument("Dataset: inconsistent input lengths or image offsets");
  for (int t : type_)
    if (t < 0 || t >= n_types_) throw std::invalid_argument("Dataset: cell type code out of range");

  // counting sort of rows by (image, type), stable in row order
  type_start_.assign(static_cast<std::size_t>(n_images()) * n_types_ + 1, 0);
  for (int img = 0; img < n_images(); ++img)
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row)
      ++type_start_[img * n_types_ + type_[row] + 1];
  std::partial_sum(type_start_.begin(), type_start_.end(), type_start_.begin());
  type_rows_.resize(N);
  std::vector<int> fill(type_start_.begin(), type_start_.end() - 1);
  for (int img = 0; img < n_images(); ++img)
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row)
      type_rows_[fill[img * n_types_ + type_[row]]++] = row;
}

std::vector<double> Dataset::image_areas(Window window) const {
  std::vector<double> areas(n_images());
  for (int img = 0; img < n_images(); ++img) {
    int start = image_offsets_[img];
    areas[img] = window_area(x_.data() + start, y_.data() + start, image_offsets_[img + 1] - start, window);
  }
  return areas;
}

void Dataset::build_grid(double radius, std::vector<Grid>& grids, std::vector<int>& bin_start,
                         std::vector<double>& gx, std::vector<double>& gy, std::vector<int>& gtype,
                         std::vector<int>& grow) const {
  std::size_t N = x_.size();
  grids.assign(n_images(), Grid{});
  bin_start.clear();
  gx.resize(N);
  gy.resize(N);
  gtype.resize(N);
  grow.resize(N);

  std::vector<long long> bin_of;
  for (int img = 0; img < n_images(); ++img) {
    int start = image_offsets_[img], n = image_offsets_[img + 1] - start;
    Grid& g = grids[img];
    g.bin_offset = bin_start.size();
    if (n == 0) {
      bin_start.push_back(start);
      continue;
    }
    auto [xmin, xmax] = std::minmax_element(x_.begin() + start, x_.begin() + start + n);
    auto [ymin, ymax] = std::minmax_element(y_.begin() + start, y_.begin() + start + n);
    double ex = *xmax - *xmin, ey = *ymax - *ymin;
    // side >= radius keeps every neighbour in the 3x3 block; the other terms
    // cap the number of bins at about 9n
    g.side = std::max({radius, std::sqrt(ex * ey / n), std::max(ex, ey) / (4.0 * n)});
    g.xmin = *xmin;
    g.ymin = *ymin;
    g.nbx = static_cast<long long>(ex / g.side) + 1;
    g.nby = static_cast<long long>(ey / g.side) + 1;
    std::size_t nbins = static_cast<std::size_t>(g.nbx * g.nby);

    // stable counting sort by bin over rows already ordered by type -> (bin, type)
    std::vector<int> counts(nbins + 1, 0);
    bin_of.resize(n);
    int first = img * n_types_;
    for (int i = type_start_[first]; i < type_start_[first + n_types_]; ++i) {
      int row = type_rows_[i];
      long long bx = static_cast<long long>((x_[row] - g.xmin) / g.side);
      long long by = static_cast<long long>((y_[row] - g.ymin) / g.side);
      bin_of[row - start] = by * g.nbx + bx;
      ++counts[bin_of[row - start] + 1];
    }
    std::partial_sum(counts.begin(), counts.end(), counts.begin());
    for (int c : counts) bin_start.push_back(start + c);
    for (int i = type_start_[first]; i < type_start_[first + n_types_]; ++i) {
      int row = type_rows_[i];
      int pos = start + counts[bin_of[row - start]]++;
      gx[pos] = x_[row];
      gy[pos] = y_[row];
      gtype[pos] = type_[row];
      grow[pos] = row;
    }
  }
}

void Dataset::build_radius_index(double r) {
  if (!(r > 0)) throw std::invalid_argument("r must be positive");
  r_ = r;
  build_grid(r, grids_, bin_start_, gx_, gy_, gtype_, grow_);
}

ModelData Dataset::poisson_model_data(const std::vector<double>& image_area, int from, int to) const {
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_;
  ModelData out;
  for (int img = 0; img < n_images(); ++img) {
    int n_from = type_count(img, from), n_to = type_count(img, to);
    if (n_from == 0 || n_to == 0) continue;
    const Grid& g = grids_[img];
    double density = (static_cast<double>(n_to) / image_area[img]) * (pi * r_ * r_);
    int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      int row = type_rows_[i];
      double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      int count = 0;
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          auto lo = gtype_.begin() + bin_start_[b], hi = gtype_.begin() + bin_start_[b + 1];
          auto [t0, t1] = std::equal_range(lo, hi, to);
          for (auto p = t0 - gtype_.begin(); p < t1 - gtype_.begin(); ++p) {
            double dx = gx_[p] - px, dy = gy_[p] - py;
            count += dx * dx + dy * dy <= r2;
          }
        }
      }
      out.row.push_back(row);
      out.image.push_back(img);
      out.n.push_back(count);
      out.density.push_back(density);
    }
  }
  return out;
}

std::vector<double> Dataset::pair_neighbour_totals(bool knn) const {
  const std::size_t T = static_cast<std::size_t>(n_types_);
  std::vector<double> out(static_cast<std::size_t>(n_images()) * T * T, 0.0);
  if (knn) {
    if (k_ == 0) throw std::logic_error("call build_knn first");
    for (int img = 0; img < n_images(); ++img) {
      double* M = out.data() + static_cast<std::size_t>(img) * T * T;
      for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row) {
        double* m = M + static_cast<std::size_t>(type_[row]) * T;
        const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
        for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0) m[type_[nb[slot]]] += 1.0;
      }
    }
    return out;
  }
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const double r2 = r_ * r_;
  for (int img = 0; img < n_images(); ++img) {
    const Grid& g = grids_[img];
    double* M = out.data() + static_cast<std::size_t>(img) * T * T;
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row) {
      double* m = M + static_cast<std::size_t>(type_[row]) * T;
      double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy <= r2) m[gtype_[p]] += 1.0;
          }
        }
      }
    }
  }
  return out;
}

ModelData Dataset::rl_model_data(int from, int to) const {
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const double r2 = r_ * r_;
  const bool self = from == to;
  ModelData out;
  for (int img = 0; img < n_images(); ++img) {
    const int n_all = image_offsets_[img + 1] - image_offsets_[img];
    const int n_from = type_count(img, from), n_to = type_count(img, to);
    if (n_from == 0 || n_to == 0) continue;
    const int cand = self ? n_all - 1 : n_all - n_from;
    if (cand <= 0 || (self && n_to < 2)) continue;
    const double share = self ? static_cast<double>(n_to - 1) / cand : static_cast<double>(n_to) / cand;
    const Grid& g = grids_[img];
    const int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      const int row = type_rows_[i];
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      int count = 0, neighbours = 0;
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            if (grow_[p] == row) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2) continue;
            if (self || gtype_[p] != from) ++neighbours;
            count += gtype_[p] == to;
          }
        }
      }
      out.row.push_back(row);
      out.image.push_back(img);
      out.n.push_back(count);
      out.density.push_back(share * neighbours);
    }
  }
  return out;
}

std::vector<double> Dataset::weighted_phi_sums(int from, int to, int design, bool edge_correct) const {
  const bool self = from == to;
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_, disc = pi * r_ * r_;
  if (design == 1 && k_ == 0) throw std::logic_error("call build_knn first");
  if (design != 1 && grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  if (design <= 1 && is_context_.empty()) throw std::logic_error("call build_context first");
  if (design == 0 && context_count_.empty() && n_images() > 0) throw std::logic_error("call build_context after build_radius_index");
  if (design == 2 && windows_.empty() && n_images() > 0) throw std::logic_error("call build_intensity first");
  std::vector<double> out(static_cast<std::size_t>(n_images()) * 3, 0.0), c;
  std::vector<double> h, qx, qy;
  auto candidate = [&](int row, int a_row) {
    const int t = type_[row];
    if (design == 2) return t == to && row != a_row;
    if (!is_context_[t]) return false;
    return self ? row != a_row : t != from;
  };
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    const int n_from = type_count(img, from);
    if (n_from == 0) continue;
    c.assign(static_cast<std::size_t>(end - start), 0.0);
    const int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      const int row = type_rows_[i];
      if (design == 1) {
        const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
        for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && candidate(nb[slot], row)) c[nb[slot] - start] += 1.0;
        continue;
      }
      double wa = 1.0, la = 0.0, ea = 1.0;
      if (design == 0) {
        la = (context_count_[row] - self) / context_area_[row];
        if (!(la > 0)) continue;
        ea = disc / context_area_[row];
      } else {
        wa = weight_[row];
      }
      const Grid& g = grids_[img];
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            const int j = grow_[p];
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2 || !candidate(j, row)) continue;
            double w;
            if (design == 0) {
              double lb = (context_count_[j] - self) / context_area_[j];
              if (!(lb > 0)) continue;
              w = ea * la / lb;
            } else {
              double e = edge_correct ? translation_weight(windows_[img], dx, dy, h, qx, qy) : 1.0;
              w = e * wa * weight_[j];
            }
            c[j - start] += w;
          }
        }
      }
    }
    double s1 = 0.0, s2 = 0.0, n = 0.0;
    for (int row = start; row < end; ++row) {
      // a candidate for some REF cell; for a self-pair every other context (or TARGET) cell
      const int t = type_[row];
      const bool cand = design == 2 ? t == to : (is_context_[t] && (self || t != from));
      if (!cand) continue;
      s1 += c[row - start]; s2 += c[row - start] * c[row - start]; n += 1.0;
    }
    out[static_cast<std::size_t>(img) * 3] = s1;
    out[static_cast<std::size_t>(img) * 3 + 1] = s2;
    out[static_cast<std::size_t>(img) * 3 + 2] = n;
  }
  return out;
}

std::vector<double> Dataset::kontextual_sums(int from, int to) const {
  const bool self = from == to;
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_, disc = pi * r_ * r_;
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  if (is_context_.empty() || (context_count_.empty() && n_images() > 0))
    throw std::logic_error("call build_context after build_radius_index");
  if (!is_context_[to]) throw std::invalid_argument("Kontextual: the TARGET cell type must belong to the context");
  std::vector<double> out(static_cast<std::size_t>(n_images()) * 7, 0.0), c;
  auto is_cand = [&](int t) { return is_context_[t] && (self || t != from); };
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    const int n_from = type_count(img, from);
    double* o = out.data() + static_cast<std::size_t>(img) * 7;
    c.assign(static_cast<std::size_t>(end - start), 0.0);
    double D = 0.0, raw = 0.0;
    const Grid& g = grids_[img];
    const int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      const int row = type_rows_[i];
      const double la = (context_count_[row] - self) / context_area_[row];
      if (!(la > 0)) continue;
      D += la;
      const double ea = disc / context_area_[row];
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            const int j = grow_[p];
            if (j == row || !is_cand(type_[j])) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2) continue;
            const double lb = (context_count_[j] - self) / context_area_[j];
            if (!(lb > 0)) continue;
            c[j - start] += ea * la / lb;
            if (type_[j] == to) raw += 1.0;
          }
        }
    }
    double O = 0.0, L = 0.0, Q = 0.0, M = 0.0, nP = 0.0;
    for (int row = start; row < end; ++row) {
      if (is_context_[type_[row]]) nP += 1.0;
      if (!is_cand(type_[row])) continue;
      const double s = c[row - start];
      L += s; Q += s * s; M += 1.0;
      if (type_[row] == to) O += s;
    }
    o[0] = O; o[1] = L; o[2] = Q; o[3] = M; o[4] = D; o[5] = raw; o[6] = nP;
  }
  return out;
}

std::vector<double> Dataset::kontextual_sums_weighted(const std::vector<double>& a, const std::vector<double>& b,
                                                      bool self) const {
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_, disc = pi * r_ * r_;
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  if (is_context_.empty() || (context_count_.empty() && n_images() > 0))
    throw std::logic_error("call build_context after build_radius_index");
  if (a.size() != x_.size() || b.size() != x_.size())
    throw std::invalid_argument("Kontextual: one REF and one TARGET weight per cell");
  std::vector<double> out(static_cast<std::size_t>(n_images()) * 10, 0.0), c;
  const double s1 = self ? 1.0 : 0.0;
  auto is_cand = [&](int row) { return is_context_[type_[row]] && (self || a[row] < 1.0); };
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    double* o = out.data() + static_cast<std::size_t>(img) * 10;
    c.assign(static_cast<std::size_t>(end - start), 0.0);
    double D = 0.0, raw = 0.0;
    const Grid& g = grids_[img];
    // REF cells in the order kontextual_sums() visits them (by type, then row), so 0/1 weights add up alike
    for (int t = 0; t < n_types_; ++t) {
      for (int i = type_start_[img * n_types_ + t]; i < type_start_[img * n_types_ + t + 1]; ++i) {
        const int row = type_rows_[i];
        const double wa = a[row];
        if (!(wa > 0)) continue;
        const double la = (context_count_[row] - s1) / context_area_[row];
        if (!(la > 0)) continue;
        D += wa * la;
        const double ea = disc / context_area_[row];
        const double px = x_[row], py = y_[row];
        long long bx = static_cast<long long>((px - g.xmin) / g.side);
        long long by = static_cast<long long>((py - g.ymin) / g.side);
        for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
          for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
            std::size_t bin = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
            for (std::size_t p = bin_start_[bin]; p < static_cast<std::size_t>(bin_start_[bin + 1]); ++p) {
              const int j = grow_[p];
              if (j == row || !is_cand(j)) continue;
              double dx = gx_[p] - px, dy = gy_[p] - py;
              if (dx * dx + dy * dy > r2) continue;
              const double lb = (context_count_[j] - s1) / context_area_[j];
              if (!(lb > 0)) continue;
              c[j - start] += wa * (ea * la / lb);
              raw += wa * b[j];
            }
          }
      }
    }
    double O = 0.0, L = 0.0, Q = 0.0, M = 0.0, nP = 0.0, Sb = 0.0, Sb2 = 0.0, St = 0.0;
    for (int row = start; row < end; ++row) {
      St += b[row];
      if (is_context_[type_[row]]) nP += 1.0;
      if (!is_cand(row)) continue;
      const double s = c[row - start];
      L += s; Q += s * s; M += 1.0;
      O += b[row] * s;
      Sb += b[row]; Sb2 += b[row] * b[row];
    }
    o[0] = O; o[1] = L; o[2] = Q; o[3] = M; o[4] = D; o[5] = raw; o[6] = nP; o[7] = Sb; o[8] = Sb2; o[9] = St;
  }
  return out;
}

std::vector<double> Dataset::weighted_excess_sums(const std::vector<double>& a, const std::vector<int>& ptr,
                                                  const std::vector<int>& lab, const std::vector<double>& w,
                                                  int n_labels, bool allocation) const {
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const std::size_t N = x_.size();
  if (a.size() != N || ptr.size() != N + 1 || lab.size() != w.size() || static_cast<std::size_t>(ptr[N]) != lab.size())
    throw std::invalid_argument("weighted_excess_sums: one REF weight and one membership row per cell");
  const double r2 = r_ * r_;
  const std::size_t K = static_cast<std::size_t>(n_labels), width = 3 + 4 * K;
  std::vector<double> out(static_cast<std::size_t>(n_images()) * width, 0.0);
  for (int img = 0; img < n_images(); ++img) {
    double* o = out.data() + static_cast<std::size_t>(img) * width;
    const Grid& g = grids_[img];
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row) {
      // the cell's score: over its neighbours a (itself excluded), prod(1 - a) or sum a
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      double prod = 1.0, sum = 0.0;
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t bin = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[bin]; p < static_cast<std::size_t>(bin_start_[bin + 1]); ++p) {
            const int j = grow_[p];
            if (j == row) continue;
            const double wa = a[j];
            if (!(wa > 0)) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2) continue;
            prod *= 1.0 - wa;
            sum += wa;
          }
        }
      const double h = allocation ? 1.0 - prod : sum;
      const bool cand = a[row] < 1.0;
      if (cand) { o[0] += h; o[1] += h * h; o[2] += 1.0; }
      for (int m = ptr[row]; m < ptr[row + 1]; ++m) {
        double* ol = o + 3 + 4 * static_cast<std::size_t>(lab[m]);
        const double wb = w[m];
        ol[3] += wb;
        if (!cand) continue;
        ol[0] += wb * h; ol[1] += wb; ol[2] += wb * wb;
      }
    }
  }
  return out;
}

std::vector<double> Dataset::hac_phi_sums(int from, int to, int design, double h) const {
  const bool self = from == to;
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_, disc = pi * r_ * r_;
  const bool knn = design == 1 || design == 4 || design == 7, context = design <= 1 || design == 5;
  const bool any = design == 6 || design == 7;  // allocation: binary scores
  if (knn && k_ == 0) throw std::logic_error("call build_knn first");
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first (the HAC bandwidth uses it)");
  if (context && is_context_.empty()) throw std::logic_error("call build_context first");
  if (design == 0 && context_count_.empty() && n_images() > 0) throw std::logic_error("call build_context after build_radius_index");
  std::vector<double> out(static_cast<std::size_t>(n_images()) * 4, 0.0), c, e;
  std::vector<char> cand;
  auto is_cand = [&](int t) { return context ? (is_context_[t] && (self || t != from)) : (self || t != from); };
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    const int n_from = type_count(img, from);
    if (n_from == 0) continue;
    const int nc = end - start;
    c.assign(nc, 0.0); cand.assign(nc, 0);
    for (int row = start; row < end; ++row) cand[row - start] = is_cand(type_[row]);
    const Grid& g = grids_[img];
    const int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      const int row = type_rows_[i];
      if (knn) {
        const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
        for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && nb[slot] != row && cand[nb[slot] - start]) c[nb[slot] - start] += 1.0;
        continue;
      }
      double la = 0.0, ea = 1.0;
      if (design == 0) { la = (context_count_[row] - self) / context_area_[row]; if (!(la > 0)) continue; ea = disc / context_area_[row]; }
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            const int j = grow_[p];
            if (j == row || !cand[j - start]) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2) continue;
            double w = 1.0;
            if (design == 0) { double lb = (context_count_[j] - self) / context_area_[j]; if (!(lb > 0)) continue; w = ea * la / lb; }
            c[j - start] += w;
          }
        }
    }
    // allocation: c_b = 1{b has a REF cell among its own neighbours} (k-NN: b's k nearest, not the REF cells')
    if (any) {
      if (knn) {
        std::fill(c.begin(), c.end(), 0.0);
        for (int row = start; row < end; ++row) {
          if (!cand[row - start]) continue;
          const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
          for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && nb[slot] != row && type_[nb[slot]] == from) { c[row - start] = 1.0; break; }
        }
      } else {
        for (double& v : c) v = v > 0 ? 1.0 : 0.0;
      }
    }
    // residual of the TARGET indicator on (1, c_b) over the candidates
    double n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int row = start; row < end; ++row) {
      if (!cand[row - start]) continue;
      const double xv = c[row - start], yv = type_[row] == to ? 1.0 : 0.0;
      n += 1; sx += xv; sy += yv; sxx += xv * xv; sxy += xv * yv;
    }
    if (n < 3) continue;
    const double vx = sxx - sx * sx / n;
    const double slope = vx > 0 ? (sxy - sx * sy / n) / vx : 0.0, icpt = (sy - slope * sx) / n;
    e.assign(nc, 0.0);
    for (int row = start; row < end; ++row)
      if (cand[row - start]) e[row - start] = (type_[row] == to ? 1.0 : 0.0) - icpt - slope * c[row - start];
    // HAC sum over candidate pairs within h (the grid side is at least r; look far enough)
    const long long reach = static_cast<long long>(std::ceil(h / g.side));
    const double h2 = h * h, cbar = sx / n;
    // V and, for its expectation under random labelling, S1 = sum K c_a c_j and
    // S2 = sum K c_a (c_a - cbar) c_j (c_j - cbar) over the same pairs (a = j included)
    double V = 0.0, S1 = 0.0, S2 = 0.0;
    for (int row = start; row < end; ++row) {
      const int a = row - start;
      if (!cand[a] || c[a] == 0.0) continue;
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      for (long long yy = std::max(0LL, by - reach); yy <= std::min(g.nby - 1, by + reach); ++yy)
        for (long long xx = std::max(0LL, bx - reach); xx <= std::min(g.nbx - 1, bx + reach); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            const int j = grow_[p] - start;
            if (!cand[j] || c[j] == 0.0) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py, d2 = dx * dx + dy * dy;
            if (d2 > h2) continue;
            const double K = 1.0 - std::sqrt(d2) / h, cc = K * c[a] * c[j];
            V += cc * e[a] * e[j]; S1 += cc; S2 += cc * (c[a] - cbar) * (c[j] - cbar);
          }
        }
    }
    // Under random labelling of the candidates, Cov(y) = p(1-p) n/(n-1) (I - 11'/n) and
    // e = (I - H) y with H the projection on (1, c), so E[e e'] = p(1-p) n/(n-1) (I - H)
    // and E[V] = p(1-p) G with G = n/(n-1) sum K c_a c_j (I - H)_aj.
    const double G = n / (n - 1) * (sxx - S1 / n - (vx > 0 ? S2 / vx : 0.0));
    out[static_cast<std::size_t>(img) * 4] = V;
    out[static_cast<std::size_t>(img) * 4 + 1] = sx;
    out[static_cast<std::size_t>(img) * 4 + 2] = n;
    out[static_cast<std::size_t>(img) * 4 + 3] = G;
  }
  return out;
}

// hac_phi_sums for one REF type and every non-self TARGET at once: one neighbour pass per image.
// With mu = icpt_t + slope_t c and e = y_t - mu, V_t = sum K c_a c_j e_a e_j expands into
//   V1_t - 2 (icpt_t A1_t + slope_t A2_t) + icpt_t^2 S1 + 2 icpt_t slope_t Sc + slope_t^2 Scc,
// where V1_t sums K c_a c_j over pairs of TARGET cells, A1_t = sum_{a in t} c_a W1_a and
// A2_t = sum_{a in t} c_a W2_a with W1_a = sum_j K c_j, W2_a = sum_j K c_j^2, S1 = sum_a c_a W1_a,
// Sc = sum_a c_a W2_a and Scc = sum_a c_a^2 W2_a. Output per image: V_t for t = 0..T-1 (NaN for
// the REF type and non-candidates), then sum c, n and G, as hac_phi_sums.
std::vector<double> Dataset::hac_phi_sums_ref(int from, int design, double h) const {
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_, disc = pi * r_ * r_;
  const bool knn = design == 1 || design == 4 || design == 7, context = design <= 1 || design == 5;
  const bool any = design == 6 || design == 7;  // allocation: binary scores
  if (knn && k_ == 0) throw std::logic_error("call build_knn first");
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first (the HAC bandwidth uses it)");
  if (context && is_context_.empty()) throw std::logic_error("call build_context first");
  if (design == 0 && context_count_.empty() && n_images() > 0) throw std::logic_error("call build_context after build_radius_index");
  const int T = n_types_; const std::size_t W = static_cast<std::size_t>(T) + 3;
  std::vector<double> out(static_cast<std::size_t>(n_images()) * W, std::nan("")), c;
  std::vector<char> cand;
  std::vector<double> sy(T), sxy(T), V1(T), A1(T), A2(T);
  auto is_cand = [&](int t) { return context ? (is_context_[t] && t != from) : (t != from); };
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    const int n_from = type_count(img, from);
    double* o = out.data() + static_cast<std::size_t>(img) * W;
    for (std::size_t k = static_cast<std::size_t>(T); k < W; ++k) o[k] = 0.0;
    if (n_from == 0) continue;
    const int nc = end - start;
    c.assign(nc, 0.0); cand.assign(nc, 0);
    for (int row = start; row < end; ++row) cand[row - start] = is_cand(type_[row]);
    const Grid& g = grids_[img];
    const int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      const int row = type_rows_[i];
      if (knn) {
        const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
        for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && nb[slot] != row && cand[nb[slot] - start]) c[nb[slot] - start] += 1.0;
        continue;
      }
      double la = 0.0, ea = 1.0;
      if (design == 0) { la = context_count_[row] / context_area_[row]; if (!(la > 0)) continue; ea = disc / context_area_[row]; }
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            const int j = grow_[p];
            if (j == row || !cand[j - start]) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2) continue;
            double w = 1.0;
            if (design == 0) { double lb = context_count_[j] / context_area_[j]; if (!(lb > 0)) continue; w = ea * la / lb; }
            c[j - start] += w;
          }
        }
    }
    // allocation: c_b = 1{b has a REF cell among its own neighbours} (k-NN: b's k nearest, not the REF cells')
    if (any) {
      if (knn) {
        std::fill(c.begin(), c.end(), 0.0);
        for (int row = start; row < end; ++row) {
          if (!cand[row - start]) continue;
          const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
          for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && nb[slot] != row && type_[nb[slot]] == from) { c[row - start] = 1.0; break; }
        }
      } else {
        for (double& v : c) v = v > 0 ? 1.0 : 0.0;
      }
    }
    double n = 0, sx = 0, sxx = 0;
    std::fill(sy.begin(), sy.end(), 0.0); std::fill(sxy.begin(), sxy.end(), 0.0);
    for (int row = start; row < end; ++row) {
      if (!cand[row - start]) continue;
      const double xv = c[row - start]; const int t = type_[row];
      n += 1; sx += xv; sxx += xv * xv; sy[t] += 1; sxy[t] += xv;
    }
    if (n < 3) continue;
    const long long reach = static_cast<long long>(std::ceil(h / g.side));
    const double h2 = h * h;
    std::fill(V1.begin(), V1.end(), 0.0); std::fill(A1.begin(), A1.end(), 0.0); std::fill(A2.begin(), A2.end(), 0.0);
    double S1 = 0.0, Sc = 0.0, Scc = 0.0;
    for (int row = start; row < end; ++row) {
      const int a = row - start;
      if (!cand[a] || c[a] == 0.0) continue;
      const int ta = type_[row];
      const double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      double W1 = 0.0, W2 = 0.0, same = 0.0;
      for (long long yy = std::max(0LL, by - reach); yy <= std::min(g.nby - 1, by + reach); ++yy)
        for (long long xx = std::max(0LL, bx - reach); xx <= std::min(g.nbx - 1, bx + reach); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
            const int j = grow_[p] - start;
            if (!cand[j] || c[j] == 0.0) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py, d2 = dx * dx + dy * dy;
            if (d2 > h2) continue;
            const double Kc = (1.0 - std::sqrt(d2) / h) * c[j];
            W1 += Kc; W2 += Kc * c[j];
            if (gtype_[p] == ta) same += Kc;
          }
        }
      V1[ta] += c[a] * same; A1[ta] += c[a] * W1; A2[ta] += c[a] * W2;
      S1 += c[a] * W1; Sc += c[a] * W2; Scc += c[a] * c[a] * W2;
    }
    const double vx = sxx - sx * sx / n, cbar = sx / n;
    for (int t = 0; t < T; ++t) {
      if (t == from || !is_cand(t)) continue;
      const double slope = vx > 0 ? (sxy[t] - sx * sy[t] / n) / vx : 0.0, icpt = (sy[t] - slope * sx) / n;
      o[t] = V1[t] - 2.0 * (icpt * A1[t] + slope * A2[t]) + icpt * icpt * S1 + 2.0 * icpt * slope * Sc + slope * slope * Scc;
    }
    const double S2 = Scc - 2.0 * cbar * Sc + cbar * cbar * S1;
    o[T] = sx; o[T + 1] = n; o[T + 2] = n / (n - 1) * (sxx - S1 / n - (vx > 0 ? S2 / vx : 0.0));
  }
  return out;
}

std::vector<double> Dataset::pair_neighbour_sq_totals(bool knn) const {
  const std::size_t T = static_cast<std::size_t>(n_types_);
  std::vector<double> out(static_cast<std::size_t>(n_images()) * T * T, 0.0);
  if (knn && k_ == 0) throw std::logic_error("call build_knn first");
  if (!knn && grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const double r2 = r_ * r_;
  std::vector<double> c;  // c[(row - start) * T + from]
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    c.assign(static_cast<std::size_t>(end - start) * T, 0.0);
    if (knn) {
      // c_b[from] = number of `from` cells with b among their k nearest
      for (int row = start; row < end; ++row) {
        const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
        for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0)
          c[static_cast<std::size_t>(nb[slot] - start) * T + type_[row]] += 1.0;
      }
    } else {
      const Grid& g = grids_[img];
      for (int row = start; row < end; ++row) {
        double* cb = c.data() + static_cast<std::size_t>(row - start) * T;
        double px = x_[row], py = y_[row];
        long long bx = static_cast<long long>((px - g.xmin) / g.side);
        long long by = static_cast<long long>((py - g.ymin) / g.side);
        for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
          for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
            std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
            for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
              if (grow_[p] == row) continue;
              double dx = gx_[p] - px, dy = gy_[p] - py;
              if (dx * dx + dy * dy <= r2) cb[gtype_[p]] += 1.0;
            }
          }
        }
      }
    }
    double* M = out.data() + static_cast<std::size_t>(img) * T * T;
    for (int row = start; row < end; ++row) {
      const double* cb = c.data() + static_cast<std::size_t>(row - start) * T;
      const std::size_t to = static_cast<std::size_t>(type_[row]);
      for (std::size_t from = 0; from < T; ++from) M[from * T + to] += cb[from] * cb[from];
    }
  }
  return out;
}

std::vector<double> Dataset::pair_neighbour_out_sq_totals(bool knn) const {
  if (!knn) return pair_neighbour_sq_totals(false);
  if (k_ == 0) throw std::logic_error("call build_knn first");
  const std::size_t T = static_cast<std::size_t>(n_types_);
  std::vector<double> out(static_cast<std::size_t>(n_images()) * T * T, 0.0), cnt(T);
  for (int img = 0; img < n_images(); ++img) {
    double* M = out.data() + static_cast<std::size_t>(img) * T * T;
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row) {
      std::fill(cnt.begin(), cnt.end(), 0.0);
      const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
      for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && nb[slot] != row) cnt[type_[nb[slot]]] += 1.0;
      const std::size_t to = static_cast<std::size_t>(type_[row]);
      for (std::size_t from = 0; from < T; ++from) M[from * T + to] += cnt[from] * cnt[from];
    }
  }
  return out;
}

void Dataset::own_neighbour_counts(int img, bool knn, std::vector<double>& c) const {
  const std::size_t T = static_cast<std::size_t>(n_types_);
  const int start = image_offsets_[img], end = image_offsets_[img + 1];
  c.assign(static_cast<std::size_t>(end - start) * T, 0.0);
  if (knn) {
    for (int row = start; row < end; ++row) {
      double* cb = c.data() + static_cast<std::size_t>(row - start) * T;
      const int* nb = knn_.data() + static_cast<std::size_t>(row) * k_;
      for (int slot = 0; slot < k_; ++slot) if (nb[slot] >= 0 && nb[slot] != row) cb[type_[nb[slot]]] += 1.0;
    }
    return;
  }
  const Grid& g = grids_[img];
  const double r2 = r_ * r_;
  for (int row = start; row < end; ++row) {
    double* cb = c.data() + static_cast<std::size_t>(row - start) * T;
    const double px = x_[row], py = y_[row];
    long long bx = static_cast<long long>((px - g.xmin) / g.side);
    long long by = static_cast<long long>((py - g.ymin) / g.side);
    for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
      for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
        std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
        for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
          if (grow_[p] == row) continue;
          double dx = gx_[p] - px, dy = gy_[p] - py;
          if (dx * dx + dy * dy <= r2) cb[gtype_[p]] += 1.0;
        }
      }
  }
}

std::vector<double> Dataset::pair_neighbour_any_totals(bool knn) const {
  if (knn && k_ == 0) throw std::logic_error("call build_knn first");
  if (!knn && grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const std::size_t T = static_cast<std::size_t>(n_types_);
  std::vector<double> out(static_cast<std::size_t>(n_images()) * T * T, 0.0), c;
  for (int img = 0; img < n_images(); ++img) {
    own_neighbour_counts(img, knn, c);
    double* M = out.data() + static_cast<std::size_t>(img) * T * T;
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row) {
      const double* cb = c.data() + static_cast<std::size_t>(row - image_offsets_[img]) * T;
      const std::size_t to = static_cast<std::size_t>(type_[row]);
      for (std::size_t from = 0; from < T; ++from) if (cb[from] > 0) M[from * T + to] += 1.0;
    }
  }
  return out;
}

std::vector<double> Dataset::self_any_expected(bool knn) const {
  if (knn && k_ == 0) throw std::logic_error("call build_knn first");
  if (!knn && grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const std::size_t T = static_cast<std::size_t>(n_types_);
  std::vector<double> out(static_cast<std::size_t>(n_images()) * T, 0.0), c;
  for (int img = 0; img < n_images(); ++img) {
    if (image_offsets_[img + 1] - image_offsets_[img] < 2) continue;
    own_neighbour_counts(img, knn, c);
    self_any_from_counts(img, c, out.data() + static_cast<std::size_t>(img) * T);
  }
  return out;
}

void Dataset::self_any_from_counts(int img, const std::vector<double>& c, double* out_img) const {
  const std::size_t T = static_cast<std::size_t>(n_types_);
  const int start = image_offsets_[img], N = image_offsets_[img + 1] - start;
  if (N < 2) return;
  auto lchoose = [](double n, double k) { return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1); };
  std::vector<double> deg_count(static_cast<std::size_t>(N), 0.0);  // number of cells with each neighbour count
  for (int b = 0; b < N; ++b) {
    double d = 0;
    for (std::size_t t = 0; t < T; ++t) d += c[static_cast<std::size_t>(b) * T + t];
    deg_count[static_cast<std::size_t>(d)] += 1.0;
  }
  for (std::size_t a = 0; a < T; ++a) {
    const int na = type_count(img, static_cast<int>(a));
    if (na < 2) continue;
    // P(none of the other na - 1 a-labels among d of the N - 1 other cells) = C(N-1-d, na-1) / C(N-1, na-1)
    const double m = na - 1, K = N - 1, lden = lchoose(K, m);
    double S = 0;
    for (int d = 1; d < N; ++d) {
      if (deg_count[d] == 0) continue;
      // 1 - C(K - d, m) / C(K, m), with expm1 to keep accuracy when the ratio is close to 1
      S += deg_count[d] * (K - d >= m ? -std::expm1(lchoose(K - d, m) - lden) : 1.0);
    }
    out_img[a] = static_cast<double>(na) / N * S;
  }
}

void Dataset::own_neighbour_lists(int img, bool knn, std::vector<int>& ptr, std::vector<int>& nb) const {
  const int start = image_offsets_[img], end = image_offsets_[img + 1], N = end - start;
  ptr.assign(N + 1, 0); nb.clear();
  if (knn) {
    for (int row = start; row < end; ++row) {
      const int* q = knn_.data() + static_cast<std::size_t>(row) * k_;
      for (int slot = 0; slot < k_; ++slot) if (q[slot] >= 0 && q[slot] != row) nb.push_back(q[slot] - start);
      ptr[row - start + 1] = static_cast<int>(nb.size());
    }
    return;
  }
  const Grid& g = grids_[img];
  const double r2 = r_ * r_;
  for (int row = start; row < end; ++row) {
    const double px = x_[row], py = y_[row];
    long long bx = static_cast<long long>((px - g.xmin) / g.side);
    long long by = static_cast<long long>((py - g.ymin) / g.side);
    for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
      for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
        std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
        for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
          if (grow_[p] == row) continue;
          double dx = gx_[p] - px, dy = gy_[p] - py;
          if (dx * dx + dy * dy <= r2) nb.push_back(grow_[p] - start);
        }
      }
    ptr[row - start + 1] = static_cast<int>(nb.size());
  }
}

void Dataset::self_variance_from_lists(int img, const std::vector<int>& ptr, const std::vector<int>& nb,
                                       bool allocation, double* out_img) const {
  const std::size_t T = static_cast<std::size_t>(n_types_);
  const int N = image_offsets_[img + 1] - image_offsets_[img];
  bool any = false;
  for (std::size_t a = 0; a < T; ++a) { const int na = type_count(img, static_cast<int>(a)); if (na >= 2 && na < N) any = true; }
  if (!any) return;
  SelfGraph g = self_graph(N, ptr, nb, allocation);
  for (std::size_t a = 0; a < T; ++a) out_img[a] = self_graph_variance(g, type_count(img, static_cast<int>(a)), allocation);
}

void Dataset::self_variance_approx(int img, const std::vector<int>& ptr, const std::vector<int>& nb, double r,
                                   int method, int cells, double* out_img) const {
  const std::size_t T = static_cast<std::size_t>(n_types_);
  const int start = image_offsets_[img], N = image_offsets_[img + 1] - start;
  bool any = false;
  for (std::size_t a = 0; a < T; ++a) { const int na = type_count(img, static_cast<int>(a)); if (na >= 2 && na < N) any = true; }
  if (!any) return;
  SelfGraph g; g.N = N;
  std::vector<int> d(N);
  int dmax = 0;
  for (int l = 0; l < N; ++l) { d[l] = ptr[l + 1] - ptr[l]; dmax = std::max(dmax, d[l]); }
  g.deg_count.assign(dmax + 1, 0.0);
  for (int l = 0; l < N; ++l) g.deg_count[d[l]] += 1;
  std::vector<double>& A = g.h2_coef;
  A.assign(2 * dmax + 2, 0.0);
  for (int x = 0; x <= dmax; ++x) {
    if (!g.deg_count[x]) continue;
    A[x] -= 2.0 * (N - 1) * g.deg_count[x];
    A[2 * x] -= g.deg_count[x];
    for (int y = 0; y <= dmax; ++y) if (g.deg_count[y]) A[x + y] += g.deg_count[x] * g.deg_count[y];
  }
  // weight w at a fractional index x, split between its two neighbouring integers
  auto add = [&](double x, double w) {
    if (x < 0) x = 0;
    const std::size_t f = static_cast<std::size_t>(std::floor(x));
    const double frac = x - static_cast<double>(f);
    A[f] += w * (1 - frac);
    if (frac > 0) A[f + 1] += w * frac;
  };
  // linked pairs (a radius graph links both ways): T = 1, each unordered pair once, counted twice
  for (int l = 0; l < N; ++l)
    for (int e = ptr[l]; e < ptr[l + 1]; ++e) if (nb[e] > l) { const int k = nb[e]; A[d[l]] += 2; A[d[k]] += 2; A[d[l] + d[k]] -= 2; }
  // unlinked pairs at distance (r, 2r]: the geometric guess of their shared neighbours, from a grid of side 2r
  const double r2 = r * r, side = 2 * r, disc = 3.14159265358979323846 * r2;   // pi (MSVC has no M_PI)
  const double* X = x_.data() + start; const double* Y = y_.data() + start;
  double xmin = X[0], ymin = Y[0], xmax = X[0], ymax = Y[0];
  for (int i = 1; i < N; ++i) { xmin = std::min(xmin, X[i]); xmax = std::max(xmax, X[i]); ymin = std::min(ymin, Y[i]); ymax = std::max(ymax, Y[i]); }
  const long long nbx = static_cast<long long>((xmax - xmin) / side) + 1, nby = static_cast<long long>((ymax - ymin) / side) + 1;
  std::vector<int> bstart(static_cast<std::size_t>(nbx * nby) + 1, 0), bcell(N), bin(N);
  for (int i = 0; i < N; ++i) { bin[i] = static_cast<int>(static_cast<long long>((Y[i] - ymin) / side) * nbx + static_cast<long long>((X[i] - xmin) / side)); ++bstart[bin[i] + 1]; }
  for (std::size_t b = 0; b + 1 < bstart.size(); ++b) bstart[b + 1] += bstart[b];
  { std::vector<int> pos(bstart.begin(), bstart.end() - 1); for (int i = 0; i < N; ++i) bcell[pos[bin[i]]++] = i; }
  auto geo_pairs = [&](int l, bool upper, double w) {   // l's unlinked partners at distance (r, 2r], k > l if upper
    const long long bx = static_cast<long long>((X[l] - xmin) / side), by = static_cast<long long>((Y[l] - ymin) / side);
    for (long long yy = std::max(0LL, by - 1); yy <= std::min(nby - 1, by + 1); ++yy)
      for (long long xx = std::max(0LL, bx - 1); xx <= std::min(nbx - 1, bx + 1); ++xx) {
        const long long b = yy * nbx + xx;
        for (int q = bstart[b]; q < bstart[b + 1]; ++q) {
          const int k = bcell[q];
          if (k == l || (upper && k < l)) continue;
          const double dx = X[k] - X[l], dy = Y[k] - Y[l], dd = dx * dx + dy * dy;
          if (dd <= r2 || dd > 4 * r2) continue;
          const double t = std::sqrt(dd);
          const double lens = 2 * r2 * std::acos(t / (2 * r)) - (t / 2) * std::sqrt(4 * r2 - dd);
          const double c = (d[l] + d[k]) / (2 * disc) * lens;
          const int s = d[l] + d[k];
          add(s - c, w); add(s, -w);
        }
      }
  };
  if (method == 1 || method == 3)
    for (int l = 0; l < N; ++l) geo_pairs(l, true, 2.0);
  if (method == 2 || method == 3) {
    // a pseudo-random sample of cells: the m with the smallest hash of their coordinates, so the sample does not depend
    // on the order of the rows (cells at the same place contribute the same) and is the same in R and Python
    auto mix = [](std::uint64_t z) {
      z += 0x9E3779B97F4A7C15ULL; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
      return z ^ (z >> 31); };
    std::vector<std::pair<std::uint64_t, int>> key(N);
    for (int i = 0; i < N; ++i) {
      std::uint64_t bx, by; const double xi = X[i] + 0.0, yi = Y[i] + 0.0;   // + 0.0: -0 and 0 hash alike
      std::memcpy(&bx, &xi, 8); std::memcpy(&by, &yi, 8);
      key[i] = {mix(bx ^ mix(by)), i};
    }
    const int m = std::min(N, std::max(1, cells));
    std::nth_element(key.begin(), key.begin() + (m - 1), key.end());
    std::vector<int> sample(m);
    for (int i = 0; i < m; ++i) sample[i] = key[i].second;
    const double scale = static_cast<double>(N) / m;
    std::vector<int> common(N, 0), touched;
    std::vector<char> linked(N, 0);
    std::vector<long long> hits(A.size(), 0);   // whole counts, scaled once at the end: the same sums in any row order
    for (int i = 0; i < m; ++i) {
      const int l = sample[i];
      touched.clear();
      for (int e = ptr[l]; e < ptr[l + 1]; ++e) linked[nb[e]] = 1;
      for (int e = ptr[l]; e < ptr[l + 1]; ++e) {
        const int v = nb[e];
        for (int f = ptr[v]; f < ptr[v + 1]; ++f) { const int k = nb[f]; if (k == l || linked[k]) continue; if (!common[k]) touched.push_back(k); ++common[k]; }
      }
      for (int k : touched) { const int s = d[l] + d[k]; ++hits[s - common[k]]; --hits[s]; common[k] = 0; }
      for (int e = ptr[l]; e < ptr[l + 1]; ++e) linked[nb[e]] = 0;
      if (method == 3) geo_pairs(l, false, -scale);
    }
    for (std::size_t x = 0; x < hits.size(); ++x) if (hits[x]) A[x] += scale * static_cast<double>(hits[x]);
  }
  for (std::size_t a = 0; a < T; ++a) out_img[a] = self_graph_variance(g, type_count(img, static_cast<int>(a)), true);
}

std::vector<double> Dataset::self_variances(bool knn, bool allocation, int self_var_method, int self_var_cells) const {
  if (knn && k_ == 0) throw std::logic_error("call build_knn first");
  if (!knn && grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  const std::size_t T = static_cast<std::size_t>(n_types_);
  std::vector<double> out(static_cast<std::size_t>(n_images()) * T, 0.0);
  std::vector<int> ptr, nb;
  for (int img = 0; img < n_images(); ++img) {
    if (image_offsets_[img + 1] - image_offsets_[img] < 3) continue;
    own_neighbour_lists(img, knn, ptr, nb);
    if (!knn && allocation && self_var_method != 0)
      self_variance_approx(img, ptr, nb, r_, self_var_method, self_var_cells, out.data() + static_cast<std::size_t>(img) * T);
    else self_variance_from_lists(img, ptr, nb, allocation, out.data() + static_cast<std::size_t>(img) * T);
  }
  return out;
}

Dataset::MultiRadiusTotals Dataset::multi_radius_totals(const std::vector<double>& radii, bool allocation,
                                                       int self_var_method, int self_var_cells) {
  const int K = static_cast<int>(radii.size());
  if (K < 1) throw std::invalid_argument("radii must not be empty");
  for (int k = 0; k < K; ++k)
    if (!(radii[k] > 0) || (k > 0 && !(radii[k] > radii[k - 1]))) throw std::invalid_argument("radii must be positive and increasing");
  build_radius_index(radii.back());
  const std::size_t T = static_cast<std::size_t>(n_types_);
  const std::size_t per = static_cast<std::size_t>(n_images()) * T * T;
  MultiRadiusTotals out;
  out.self_var.assign(K, std::vector<double>(static_cast<std::size_t>(n_images()) * T, 0.0));
  if (allocation) {
    out.any.assign(K, std::vector<double>(per, 0.0));
    out.self_expected.assign(K, std::vector<double>(static_cast<std::size_t>(n_images()) * T, 0.0));
  } else {
    out.totals.assign(K, std::vector<double>(per, 0.0));
    out.sq.assign(K, std::vector<double>(per, 0.0));
  }
  std::vector<double> r2(K);
  for (int k = 0; k < K; ++k) r2[k] = radii[k] * radii[k];
  std::vector<int> e_row, e_col, e_type, e_k, order, k_start, lptr, lnb;
  std::vector<double> c;
  for (int img = 0; img < n_images(); ++img) {
    const int start = image_offsets_[img], end = image_offsets_[img + 1];
    const Grid& g = grids_[img];
    const int N = end - start;
    c.assign(static_cast<std::size_t>(N) * T, 0.0);
    if (K == 1) {
      // one radius: the neighbour lists from one search, and each cell's own counts c_b[t] from them
      own_neighbour_lists(img, false, lptr, lnb);
      for (int l = 0; l < N; ++l)
        for (int e = lptr[l]; e < lptr[l + 1]; ++e) c[static_cast<std::size_t>(l) * T + static_cast<std::size_t>(type_[start + lnb[e]])] += 1.0;
    } else {
      e_row.clear(); e_col.clear(); e_type.clear(); e_k.clear();
      // every neighbour pair (b, p), p != b, within the largest radius, with the smallest radius that contains it
      for (int row = start; row < end; ++row) {
        const double px = x_[row], py = y_[row];
        long long bx = static_cast<long long>((px - g.xmin) / g.side);
        long long by = static_cast<long long>((py - g.ymin) / g.side);
        for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy)
          for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
            std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
            for (std::size_t p = bin_start_[b]; p < static_cast<std::size_t>(bin_start_[b + 1]); ++p) {
              if (grow_[p] == row) continue;
              double dx = gx_[p] - px, dy = gy_[p] - py, d2 = dx * dx + dy * dy;
              if (!(d2 <= r2[K - 1])) continue;
              e_row.push_back(row - start); e_col.push_back(grow_[p] - start); e_type.push_back(gtype_[p]);
              e_k.push_back(static_cast<int>(std::lower_bound(r2.begin(), r2.end(), d2) - r2.begin()));
            }
          }
      }
      // counting sort of the pairs by radius
      k_start.assign(K + 1, 0);
      for (int k : e_k) ++k_start[k + 1];
      for (int k = 0; k < K; ++k) k_start[k + 1] += k_start[k];
      order.assign(e_k.size(), 0);
      { std::vector<int> pos(k_start.begin(), k_start.end() - 1);
        for (std::size_t e = 0; e < e_k.size(); ++e) order[pos[e_k[e]]++] = static_cast<int>(e); }
    }
    // cumulative own counts c_b[t] over the radii, and each radius's totals from them
    for (int k = 0; k < K; ++k) {
      if (K > 1) for (int q = k_start[k]; q < k_start[k + 1]; ++q) { const int e = order[q];
        c[static_cast<std::size_t>(e_row[e]) * T + static_cast<std::size_t>(e_type[e])] += 1.0; }
      // the self-pair variances at this radius: the neighbour lists of the pairs within radii[k]
      const bool approx = allocation && self_var_method != 0;
      if (N >= 3 && K == 1) {
        if (approx) self_variance_approx(img, lptr, lnb, radii[k], self_var_method, self_var_cells, out.self_var[k].data() + static_cast<std::size_t>(img) * T);
        else self_variance_from_lists(img, lptr, lnb, allocation, out.self_var[k].data() + static_cast<std::size_t>(img) * T);
      } else if (N >= 3) {
        lptr.assign(N + 1, 0);
        for (int q = 0; q < k_start[k + 1]; ++q) ++lptr[e_row[order[q]] + 1];
        for (int l = 0; l < N; ++l) lptr[l + 1] += lptr[l];
        lnb.assign(lptr[N], 0);
        { std::vector<int> pos(lptr.begin(), lptr.end() - 1);
          for (int q = 0; q < k_start[k + 1]; ++q) { const int e = order[q]; lnb[pos[e_row[e]]++] = e_col[e]; } }
        if (approx) self_variance_approx(img, lptr, lnb, radii[k], self_var_method, self_var_cells, out.self_var[k].data() + static_cast<std::size_t>(img) * T);
        else self_variance_from_lists(img, lptr, lnb, allocation, out.self_var[k].data() + static_cast<std::size_t>(img) * T);
      }
      if (allocation) {
        double* M = out.any[k].data() + static_cast<std::size_t>(img) * T * T;
        for (int row = start; row < end; ++row) {
          const double* cb = c.data() + static_cast<std::size_t>(row - start) * T;
          const std::size_t to = static_cast<std::size_t>(type_[row]);
          for (std::size_t from = 0; from < T; ++from) if (cb[from] > 0) M[from * T + to] += 1.0;
        }
        self_any_from_counts(img, c, out.self_expected[k].data() + static_cast<std::size_t>(img) * T);
      } else {
        double* M = out.sq[k].data() + static_cast<std::size_t>(img) * T * T;
        double* Mt = out.totals[k].data() + static_cast<std::size_t>(img) * T * T;
        for (int row = start; row < end; ++row) {
          const double* cb = c.data() + static_cast<std::size_t>(row - start) * T;
          const std::size_t ty = static_cast<std::size_t>(type_[row]);
          for (std::size_t from = 0; from < T; ++from) M[from * T + ty] += cb[from] * cb[from];
          // pair_neighbour_totals: the row as `from`, its neighbours of each type, the cell itself included
          double* m = Mt + ty * T;
          for (std::size_t t = 0; t < T; ++t) m[t] += cb[t];
          m[ty] += 1.0;
        }
      }
    }
  }
  return out;
}

void Dataset::build_windows(Window window) {
  windows_.assign(n_images(), WindowGeometry{});
  for (int img = 0; img < n_images(); ++img) {
    int start = image_offsets_[img];
    windows_[img] = window_geometry(x_.data() + start, y_.data() + start, image_offsets_[img + 1] - start, window);
  }
}

void Dataset::build_intensity(double sigma, double min_lambda, Window window) {
  if (!(sigma > 0)) throw std::invalid_argument("sigma must be positive");
  if (!(min_lambda > 0)) throw std::invalid_argument("min_lambda must be positive");
  build_windows(window);

  std::vector<Grid> grids;
  std::vector<int> bin_start, gtype, grow;
  std::vector<double> gx, gy;
  build_grid(sigma, grids, bin_start, gx, gy, gtype, grow);
  const double s2 = sigma * sigma;
  weight_.assign(x_.size(), 0.0);
  for (int img = 0; img < n_images(); ++img) {
    const Grid& g = grids[img];
    const WindowGeometry& W = windows_[img];
    for (int type = 0; type < n_types_; ++type) {
      int n_type = type_count(img, type);
      if (n_type == 0) continue;
      double floor = min_lambda * (static_cast<double>(n_type) / W.area);
      int base = type_start_[img * n_types_ + type];
      double inverse_sum = 0.0;
      for (int i = base; i < base + n_type; ++i) {
        int row = type_rows_[i];
        double px = x_[row], py = y_[row];
        long long bx = static_cast<long long>((px - g.xmin) / g.side);
        long long by = static_cast<long long>((py - g.ymin) / g.side);
        int count = -1;  // the cell itself is found at distance 0
        for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
          for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
            std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
            auto lo = gtype.begin() + bin_start[b], hi = gtype.begin() + bin_start[b + 1];
            auto [t0, t1] = std::equal_range(lo, hi, type);
            for (auto p = t0 - gtype.begin(); p < t1 - gtype.begin(); ++p) {
              double dx = gx[p] - px, dy = gy[p] - py;
              count += dx * dx + dy * dy <= s2;
            }
          }
        }
        double lambda = std::max(count / disc_window_area(W.vx, W.vy, px, py, sigma), floor);
        weight_[row] = 1.0 / lambda;
        inverse_sum += weight_[row];
      }
      double mean = inverse_sum / n_type;
      for (int i = base; i < base + n_type; ++i) weight_[type_rows_[i]] /= mean;
    }
  }
}

InhomModelData Dataset::inhom_model_data(const std::vector<double>& image_area, int from, int to,
                                         bool edge_correct) const {
  if (grids_.empty() && n_images() > 0) throw std::logic_error("call build_radius_index first");
  if (windows_.empty() && n_images() > 0) throw std::logic_error("call build_intensity first");
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_;
  const bool self = from == to;
  std::vector<double> h, qx, qy;  // scratch for translation_weight
  InhomModelData out;
  for (int img = 0; img < n_images(); ++img) {
    int n_from = type_count(img, from), n_to = type_count(img, to);
    if (n_from == 0 || n_to == 0 || (self && n_from < 2)) continue;
    const Grid& g = grids_[img];
    const WindowGeometry& W = windows_[img];
    double unit = (static_cast<double>(n_to - self) / image_area[img]) * (pi * r_ * r_);
    int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      int row = type_rows_[i];
      double px = x_[row], py = y_[row], wa = weight_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      double n = 0.0;
      int raw = 0;
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          auto lo = gtype_.begin() + bin_start_[b], hi = gtype_.begin() + bin_start_[b + 1];
          auto [t0, t1] = std::equal_range(lo, hi, to);
          for (auto p = t0 - gtype_.begin(); p < t1 - gtype_.begin(); ++p) {
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2 || (self && grow_[p] == row)) continue;
            double e = edge_correct ? translation_weight(W, dx, dy, h, qx, qy) : 1.0;
            n += e * wa * weight_[grow_[p]];
            ++raw;
          }
        }
      }
      out.row.push_back(row);
      out.image.push_back(img);
      out.n.push_back(n);
      out.n_raw.push_back(raw);
      out.weight.push_back(wa);
      out.density.push_back(wa * unit);
    }
  }
  return out;
}

void Dataset::build_context(const std::vector<int>& context_types, Window window, bool edge_correct) {
  is_context_.assign(n_types_, 0);
  for (int t : context_types) {
    if (t < 0 || t >= n_types_) throw std::invalid_argument("context cell type code out of range");
    is_context_[t] = 1;
  }
  context_count_.clear();
  context_area_.clear();
  if (grids_.empty()) return;  // the binomial design needs only the context types

  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_;
  if (edge_correct) build_windows(window);
  context_count_.assign(x_.size(), 0);
  context_area_.assign(x_.size(), pi * r_ * r_);
  for (int img = 0; img < n_images(); ++img) {
    const Grid& g = grids_[img];
    for (int row = image_offsets_[img]; row < image_offsets_[img + 1]; ++row) {
      double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      int count = 0;
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          for (int p = bin_start_[b]; p < bin_start_[b + 1]; ++p) {
            if (!is_context_[gtype_[p]]) continue;
            double dx = gx_[p] - px, dy = gy_[p] - py;
            count += dx * dx + dy * dy <= r2;
          }
        }
      }
      context_count_[row] = count;
      if (edge_correct) context_area_[row] = disc_window_area(windows_[img].vx, windows_[img].vy, px, py, r_);
    }
  }
}

InhomModelData Dataset::kontextual_model_data(int from, int to) const {
  if (context_count_.empty() && n_images() > 0)
    throw std::logic_error("call build_radius_index, then build_context");
  if (!is_context_[to]) throw std::invalid_argument("Kontextual: the TARGET cell type must belong to the context");
  const double pi = 3.141592653589793238462643383279502884;
  const double r2 = r_ * r_, disc = pi * r_ * r_;
  const bool self = from == to;
  InhomModelData out;
  for (int img = 0; img < n_images(); ++img) {
    int n_from = type_count(img, from), n_to = type_count(img, to), n_context = 0;
    for (int t = 0; t < n_types_; ++t)
      if (is_context_[t]) n_context += type_count(img, t);
    if (n_from == 0 || n_to == 0 || (self && n_from < 2)) continue;
    // the TARGET share of the context; a self-pair leaves the cell itself out
    double share = self ? static_cast<double>(n_to - 1) / (n_context - 1)
                        : static_cast<double>(n_to) / n_context;
    const Grid& g = grids_[img];
    int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      int row = type_rows_[i];
      double la = (context_count_[row] - self) / context_area_[row];  // lambda_c(x_a)
      if (!(la > 0)) continue;  // no context cell within r: weight zero
      double ea = disc / context_area_[row];
      double px = x_[row], py = y_[row];
      long long bx = static_cast<long long>((px - g.xmin) / g.side);
      long long by = static_cast<long long>((py - g.ymin) / g.side);
      double n = 0.0;
      int raw = 0;
      for (long long yy = std::max(0LL, by - 1); yy <= std::min(g.nby - 1, by + 1); ++yy) {
        for (long long xx = std::max(0LL, bx - 1); xx <= std::min(g.nbx - 1, bx + 1); ++xx) {
          std::size_t b = g.bin_offset + static_cast<std::size_t>(yy * g.nbx + xx);
          auto lo = gtype_.begin() + bin_start_[b], hi = gtype_.begin() + bin_start_[b + 1];
          auto [t0, t1] = std::equal_range(lo, hi, to);
          for (auto p = t0 - gtype_.begin(); p < t1 - gtype_.begin(); ++p) {
            double dx = gx_[p] - px, dy = gy_[p] - py;
            if (dx * dx + dy * dy > r2 || (self && grow_[p] == row)) continue;
            int j = grow_[p];
            // lambda_c(x_b) > 0: b is a context cell, and for a self-pair a lies within r of b
            n += ea * la / ((context_count_[j] - self) / context_area_[j]);
            ++raw;
          }
        }
      }
      out.row.push_back(row);
      out.image.push_back(img);
      out.n.push_back(n);
      out.n_raw.push_back(raw);
      out.weight.push_back(la * disc);
      out.density.push_back(share * la * disc);
    }
  }
  return out;
}

void Dataset::build_knn(int k, int n_threads) {
  knn_ = knn_indices(x_, y_, image_offsets_, k, n_threads);
  k_ = k;
}

BinomialModelData Dataset::binomial_model_data(int from, int to) const {
  if (k_ == 0) throw std::logic_error("call build_knn first");
  BinomialModelData out;
  for (int img = 0; img < n_images(); ++img) {
    int n_all = image_offsets_[img + 1] - image_offsets_[img];
    int n_from = type_count(img, from), n_to = type_count(img, to);
    if (n_from == 0 || n_to == 0 || n_all <= k_ || n_to == n_all) continue;
    double p0 = static_cast<double>(n_to) / n_all;
    int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      int row = type_rows_[i];
      int count = 0;
      for (int slot = 0; slot < k_; ++slot) count += type_[knn_[static_cast<std::size_t>(row) * k_ + slot]] == to;
      out.row.push_back(row);
      out.image.push_back(img);
      out.n.push_back(count);
      out.p0.push_back(p0);
    }
  }
  return out;
}

BinomialTrialsModelData Dataset::kontextual_binomial_model_data(int from, int to) const {
  if (k_ == 0) throw std::logic_error("call build_knn first");
  if (is_context_.empty()) throw std::logic_error("call build_context first");
  if (!is_context_[to]) throw std::invalid_argument("Kontextual: the TARGET cell type must belong to the context");
  BinomialTrialsModelData out;
  for (int img = 0; img < n_images(); ++img) {
    int n_all = image_offsets_[img + 1] - image_offsets_[img];
    int n_from = type_count(img, from), n_to = type_count(img, to), n_context = 0;
    for (int t = 0; t < n_types_; ++t)
      if (is_context_[t]) n_context += type_count(img, t);
    if (n_from == 0 || n_to == 0 || n_all <= k_ || n_to == n_context) continue;
    double p0 = static_cast<double>(n_to) / n_context;
    int base = type_start_[img * n_types_ + from];
    for (int i = base; i < base + n_from; ++i) {
      int row = type_rows_[i];
      int trials = 0, count = 0;
      for (int slot = 0; slot < k_; ++slot) {
        int t = type_[knn_[static_cast<std::size_t>(row) * k_ + slot]];
        trials += is_context_[t];
        count += t == to;
      }
      if (trials == 0) continue;  // no context cell among the neighbours
      out.row.push_back(row);
      out.image.push_back(img);
      out.n.push_back(count);
      out.trials.push_back(trials);
      out.p0.push_back(p0);
    }
  }
  return out;
}

}  // namespace spicyglm

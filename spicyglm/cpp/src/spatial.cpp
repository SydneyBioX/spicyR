#include "spicyglm/core.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace spicyglm {

namespace {

double cross(double ox, double oy, double ax, double ay, double bx, double by) {
  return (ax - ox) * (by - oy) - (ay - oy) * (bx - ox);
}

// Andrew's monotone chain, then the shoelace formula.
double convex_hull_area(const double* x, const double* y, std::size_t n) {
  if (n < 3) return 0.0;
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
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
  double twice = 0.0;
  for (std::size_t i = 0; i < hull.size(); ++i) {
    std::size_t a = hull[i], b = hull[(i + 1) % hull.size()];
    twice += x[a] * y[b] - x[b] * y[a];
  }
  return std::abs(twice) / 2.0;
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

void Dataset::build_radius_index(double r) {
  if (!(r > 0)) throw std::invalid_argument("r must be positive");
  r_ = r;
  std::size_t N = x_.size();
  grids_.assign(n_images(), Grid{});
  bin_start_.clear();
  gx_.resize(N);
  gy_.resize(N);
  gtype_.resize(N);

  std::vector<long long> bin_of;
  for (int img = 0; img < n_images(); ++img) {
    int start = image_offsets_[img], n = image_offsets_[img + 1] - start;
    Grid& g = grids_[img];
    g.bin_offset = bin_start_.size();
    if (n == 0) {
      bin_start_.push_back(start);
      continue;
    }
    auto [xmin, xmax] = std::minmax_element(x_.begin() + start, x_.begin() + start + n);
    auto [ymin, ymax] = std::minmax_element(y_.begin() + start, y_.begin() + start + n);
    double ex = *xmax - *xmin, ey = *ymax - *ymin;
    // side >= r keeps every neighbour in the 3x3 block; the other terms cap
    // the number of bins at about 9n
    g.side = std::max({r, std::sqrt(ex * ey / n), std::max(ex, ey) / (4.0 * n)});
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
    for (int c : counts) bin_start_.push_back(start + c);
    for (int i = type_start_[first]; i < type_start_[first + n_types_]; ++i) {
      int row = type_rows_[i];
      int pos = start + counts[bin_of[row - start]]++;
      gx_[pos] = x_[row];
      gy_[pos] = y_[row];
      gtype_[pos] = type_[row];
    }
  }
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

}  // namespace spicyglm

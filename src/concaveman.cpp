// A C++ port of concaveman (Vladimir Agafonkin, https://github.com/mapbox/concaveman, ISC licence) as run by
// the R package concaveman 1.2.0: the JavaScript bundle it evaluates in V8, with the R-tree (rbush), priority
// queue (tinyqueue), monotone convex hull and exact orientation test it uses, and the rounding its R interface
// applies on the way in (coordinates to 4 decimals by jsonlite::toJSON(), the length threshold to 15
// significant digits by sprintf("%s")). Every step follows the JavaScript, including its order of floating
// point operations, so the polygons are identical, point for point, to concaveman::concaveman().
// Authors and licences of the ported code: the file COPYRIGHTS (inst/COPYRIGHTS in the R packages).

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#include "concaveman.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <limits>
#include <stdexcept>
#include <vector>

namespace concaveman {
namespace {

// ---- V8's Math.exp and Math.log (fdlibm, from V8 11.3 src/base/ieee754.cc) -----------------------------

inline uint64_t bits_of(double d) { uint64_t b; std::memcpy(&b, &d, 8); return b; }
inline double double_of(uint64_t b) { double d; std::memcpy(&d, &b, 8); return d; }
inline int32_t high_word(double d) { return static_cast<int32_t>(bits_of(d) >> 32); }
inline uint32_t low_word(double d) { return static_cast<uint32_t>(bits_of(d) & 0xFFFFFFFFu); }
inline double with_high_word(double d, uint32_t v) {
  return double_of((bits_of(d) & 0xFFFFFFFFu) | (static_cast<uint64_t>(v) << 32));
}
inline double from_words(uint32_t hi, uint32_t lo) { return double_of((static_cast<uint64_t>(hi) << 32) | lo); }

double v8_exp(double x) {
  static const double one = 1.0, halF[2] = {0.5, -0.5}, o_threshold = 7.09782712893383973096e+02,
                      u_threshold = -7.45133219101941108420e+02,
                      ln2HI[2] = {6.93147180369123816490e-01, -6.93147180369123816490e-01},
                      ln2LO[2] = {1.90821492927058770002e-10, -1.90821492927058770002e-10},
                      invln2 = 1.44269504088896338700e+00, P1 = 1.66666666666666019037e-01,
                      P2 = -2.77777777770155933842e-03, P3 = 6.61375632143793436117e-05,
                      P4 = -1.65339022054652515390e-06, P5 = 4.13813679705723846039e-08, E = 2.718281828459045;
  static volatile double huge = 1.0e+300, twom1000 = 9.33263618503218878990e-302, two1023 = 8.988465674311579539e307;
  double y, hi = 0.0, lo = 0.0, c, t, twopk;
  int32_t k = 0, xsb;
  uint32_t hx = static_cast<uint32_t>(high_word(x));
  xsb = (hx >> 31) & 1;
  hx &= 0x7FFFFFFF;
  if (hx >= 0x40862E42) {
    if (hx >= 0x7FF00000) {
      if (((hx & 0xFFFFF) | low_word(x)) != 0) return x + x;
      return (xsb == 0) ? x : 0.0;
    }
    if (x > o_threshold) return huge * huge;
    if (x < u_threshold) return twom1000 * twom1000;
  }
  if (hx > 0x3FD62E42) {
    if (hx < 0x3FF0A2B2) {
      if (x == 1.0) return E;
      hi = x - ln2HI[xsb];
      lo = ln2LO[xsb];
      k = 1 - xsb - xsb;
    } else {
      k = static_cast<int>(invln2 * x + halF[xsb]);
      t = k;
      hi = x - t * ln2HI[0];
      lo = t * ln2LO[0];
    }
    x = hi - lo;
  } else if (hx < 0x3E300000) {
    if (huge + x > one) return one + x;
  } else {
    k = 0;
  }
  t = x * x;
  if (k >= -1021) {
    twopk = from_words(0x3FF00000 + static_cast<int32_t>(static_cast<uint32_t>(k) << 20), 0);
  } else {
    twopk = from_words(0x3FF00000 + (static_cast<uint32_t>(k + 1000) << 20), 0);
  }
  c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
  if (k == 0) return one - ((x * c) / (c - 2.0) - x);
  y = one - ((lo - (x * c) / (2.0 - c)) - hi);
  if (k >= -1021) {
    if (k == 1024) return y * 2.0 * two1023;
    return y * twopk;
  }
  return y * twopk * twom1000;
}

double v8_log(double x) {
  static const double ln2_hi = 6.93147180369123816490e-01, ln2_lo = 1.90821492927058770002e-10,
                      two54 = 1.80143985094819840000e+16, Lg1 = 6.666666666666735130e-01,
                      Lg2 = 3.999999999940941908e-01, Lg3 = 2.857142874366239149e-01,
                      Lg4 = 2.222219843214978396e-01, Lg5 = 1.818357216161805012e-01,
                      Lg6 = 1.531383769920937332e-01, Lg7 = 1.479819860511658591e-01;
  static const double zero = 0.0;
  double hfsq, f, s, z, R, w, t1, t2, dk;
  int32_t k, hx, i, j;
  uint32_t lx;
  hx = high_word(x);
  lx = low_word(x);
  k = 0;
  if (hx < 0x00100000) {
    if (((hx & 0x7FFFFFFF) | lx) == 0) return -std::numeric_limits<double>::infinity();
    if (hx < 0) return std::numeric_limits<double>::quiet_NaN();
    k -= 54;
    x *= two54;
    hx = high_word(x);
  }
  if (hx >= 0x7FF00000) return x + x;
  k += (hx >> 20) - 1023;
  hx &= 0x000FFFFF;
  i = (hx + 0x95F64) & 0x100000;
  x = with_high_word(x, static_cast<uint32_t>(hx | (i ^ 0x3FF00000)));
  k += (i >> 20);
  f = x - 1.0;
  if ((0x000FFFFF & (2 + hx)) < 3) {
    if (f == zero) {
      if (k == 0) return zero;
      dk = static_cast<double>(k);
      return dk * ln2_hi + dk * ln2_lo;
    }
    R = f * f * (0.5 - 0.33333333333333333 * f);
    if (k == 0) return f - R;
    dk = static_cast<double>(k);
    return dk * ln2_hi - ((R - dk * ln2_lo) - f);
  }
  s = f / (2.0 + f);
  dk = static_cast<double>(k);
  z = s * s;
  i = hx - 0x6147A;
  w = z * z;
  j = 0x6B851 - hx;
  t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
  t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
  i |= j;
  R = t2 + t1;
  if (i > 0) {
    hfsq = 0.5 * f * f;
    if (k == 0) return f - (hfsq - s * (hfsq + R));
    return dk * ln2_hi - ((hfsq - (s * (hfsq + R) + dk * ln2_lo)) - f);
  }
  if (k == 0) return f - s * (f - R);
  return dk * ln2_hi - ((s * (f - R) - dk * ln2_lo) - f);
}

// Math.min / Math.max, with -0 < +0
inline double js_min(double a, double b) {
  if (a < b) return a;
  if (b < a) return b;
  return std::signbit(a) ? a : b;
}
inline double js_max(double a, double b) {
  if (a > b) return a;
  if (b > a) return b;
  return std::signbit(a) ? b : a;
}

// ---- robust-orientation: the sign of orient3(a, b, c) = (ay - cy)(bx - cx) - (ax - cx)(by - cy) ---------

inline void two_product(double a, double b, double& x, double& y) {
  static const double splitter = 134217729.0;  // 2^27 + 1
  x = a * b;
  double c = splitter * a, abig = c - a, ahi = c - abig, alo = a - ahi;
  double d = splitter * b, bbig = d - b, bhi = d - bbig, blo = b - bhi;
  double err1 = x - (ahi * bhi), err2 = err1 - (alo * bhi), err3 = err2 - (ahi * blo);
  y = alo * blo - err3;
}
inline void two_sum(double a, double b, double& x, double& y) {
  x = a + b;
  double bv = x - a, av = x - bv, br = b - bv, ar = a - av;
  y = ar + br;
}

// The exact sign of a sum of products, by growing a nonoverlapping expansion (Shewchuk).
int exact_sign(const double (*terms)[2], int n) {
  std::vector<double> e;
  e.reserve(2 * n);
  for (int t = 0; t < n; ++t) {
    double hi, lo;
    two_product(terms[t][0], terms[t][1], hi, lo);
    for (double b : {lo, hi}) {
      double q = b;
      for (double& ei : e) {
        double x, y;
        two_sum(q, ei, x, y);
        ei = y;
        q = x;
      }
      e.push_back(q);
    }
  }
  for (std::size_t i = e.size(); i-- > 0;) {
    if (e[i] > 0) return 1;
    if (e[i] < 0) return -1;
  }
  return 0;
}

int orient_sign(double ax, double ay, double bx, double by, double cx, double cy) {
  static const double EPSILON = 1.1102230246251565e-16;
  static const double ERRBOUND3 = (3.0 + 16.0 * EPSILON) * EPSILON;
  double l = (ay - cy) * (bx - cx);
  double r = (ax - cx) * (by - cy);
  double det = l - r;
  double s;
  auto sgn = [](double v) { return (v > 0) - (v < 0); };
  if (l > 0) {
    if (r <= 0) return sgn(det);
    s = l + r;
  } else if (l < 0) {
    if (r >= 0) return sgn(det);
    s = -(l + r);
  } else {
    return sgn(det);
  }
  double tol = ERRBOUND3 * s;
  if (det >= tol || det <= -tol) return sgn(det);
  // ay*bx - ay*cx - cy*bx + ax*cy + cx*by - ax*by (the cy*cx terms cancel)
  const double terms[6][2] = {{ay, bx}, {-ay, cx}, {-cy, bx}, {ax, cy}, {cx, by}, {-ax, by}};
  return exact_sign(terms, 6);
}

// ---- rbush (the 2.x version in the bundle), over integer item ids ------------------------------------

using BBox = std::array<double, 4>;
inline BBox empty_bbox() {
  const double inf = std::numeric_limits<double>::infinity();
  return {inf, inf, -inf, -inf};
}
inline void extend(BBox& a, const BBox& b) {
  a[0] = js_min(a[0], b[0]);
  a[1] = js_min(a[1], b[1]);
  a[2] = js_max(a[2], b[2]);
  a[3] = js_max(a[3], b[3]);
}
inline double bbox_area(const BBox& a) { return (a[2] - a[0]) * (a[3] - a[1]); }
inline double bbox_margin(const BBox& a) { return (a[2] - a[0]) + (a[3] - a[1]); }
inline double enlarged_area(const BBox& a, const BBox& b) {
  return (js_max(b[2], a[2]) - js_min(b[0], a[0])) * (js_max(b[3], a[3]) - js_min(b[1], a[1]));
}
inline double intersection_area(const BBox& a, const BBox& b) {
  double minX = js_max(a[0], b[0]), minY = js_max(a[1], b[1]), maxX = js_min(a[2], b[2]), maxY = js_min(a[3], b[3]);
  return js_max(0, maxX - minX) * js_max(0, maxY - minY);
}
inline bool contains(const BBox& a, const BBox& b) {
  return a[0] <= b[0] && a[1] <= b[1] && b[2] <= a[2] && b[3] <= a[3];
}
inline bool intersects(const BBox& a, const BBox& b) {
  return b[0] <= a[2] && b[1] <= a[3] && b[2] >= a[0] && b[3] >= a[1];
}

struct RNode {
  std::vector<RNode*> nodes;  // children of an internal node
  std::vector<int> items;     // children of a leaf
  int height = 1;
  bool leaf = true;
  BBox bbox = empty_bbox();
  std::size_t count() const { return leaf ? items.size() : nodes.size(); }
};

// Traits: BBox to_bbox(int), double compare_min_x(int, int), double compare_min_y(int, int)
template <class Traits>
class RBush {
 public:
  explicit RBush(const Traits& traits, int max_entries = 9) : t_(traits) {
    max_entries_ = std::max(4, max_entries);
    min_entries_ = std::max(2, static_cast<int>(std::ceil(max_entries_ * 0.4)));
    clear();
  }

  RNode* root() const { return data_; }

  std::vector<int> search(const BBox& bbox) const {
    RNode* node = data_;
    std::vector<int> result;
    if (!intersects(bbox, node->bbox)) return result;
    std::vector<RNode*> to_search;
    while (node) {
      if (node->leaf) {
        for (int item : node->items)
          if (intersects(bbox, t_.to_bbox(item))) result.push_back(item);
      } else {
        for (RNode* child : node->nodes) {
          if (intersects(bbox, child->bbox)) {
            if (contains(bbox, child->bbox)) all(child, result);
            else to_search.push_back(child);
          }
        }
      }
      node = pop(to_search);
    }
    return result;
  }

  void load(std::vector<int> data) {
    if (data.empty()) return;
    if (static_cast<int>(data.size()) < min_entries_) {
      for (int item : data) insert(item);
      return;
    }
    RNode* node = build(data, 0, static_cast<int>(data.size()) - 1, 0);
    if (data_->count() != 0) throw std::logic_error("concaveman: load() into a non-empty tree");
    data_ = node;
  }

  void insert(int item) { insert_item(item, data_->height - 1); }

  void clear() { data_ = new_node(1, true); }

  void remove(int item) {
    RNode* node = data_;
    BBox bbox = t_.to_bbox(item);
    std::vector<RNode*> path;
    std::vector<long> indexes;
    long i = -1;
    RNode* parent = nullptr;
    bool going_up = false;
    while (node || !path.empty()) {
      if (!node) {  // go up
        node = path.back();
        path.pop_back();
        parent = path.empty() ? nullptr : path.back();
        i = indexes.back();
        indexes.pop_back();
        going_up = true;
      }
      if (node->leaf) {
        auto it = std::find(node->items.begin(), node->items.end(), item);
        if (it != node->items.end()) {
          node->items.erase(it);
          path.push_back(node);
          condense(path);
          return;
        }
      }
      if (!going_up && !node->leaf && contains(node->bbox, bbox)) {  // go down
        path.push_back(node);
        indexes.push_back(i);
        i = 0;
        parent = node;
        node = node->nodes.empty() ? nullptr : node->nodes[0];
      } else if (parent) {  // go right
        i++;
        node = i < static_cast<long>(parent->nodes.size()) ? parent->nodes[i] : nullptr;
        going_up = false;
      } else {
        node = nullptr;
      }
    }
  }

 private:
  Traits t_;
  int max_entries_, min_entries_;
  RNode* data_ = nullptr;
  std::deque<RNode> pool_;

  RNode* new_node(int height, bool leaf) {
    pool_.emplace_back();
    RNode* n = &pool_.back();
    n->height = height;
    n->leaf = leaf;
    return n;
  }
  static RNode* pop(std::vector<RNode*>& v) {
    if (v.empty()) return nullptr;
    RNode* n = v.back();
    v.pop_back();
    return n;
  }

  void all(RNode* node, std::vector<int>& result) const {
    std::vector<RNode*> to_search;
    while (node) {
      if (node->leaf) result.insert(result.end(), node->items.begin(), node->items.end());
      else to_search.insert(to_search.end(), node->nodes.begin(), node->nodes.end());
      node = pop(to_search);
    }
  }

  BBox child_bbox(const RNode* node, std::size_t i) const {
    return node->leaf ? t_.to_bbox(node->items[i]) : node->nodes[i]->bbox;
  }
  BBox dist_bbox(const RNode* node, std::size_t k, std::size_t p) const {
    BBox bbox = empty_bbox();
    for (std::size_t i = k; i < p; ++i) extend(bbox, child_bbox(node, i));
    return bbox;
  }
  void calc_bbox(RNode* node) const { node->bbox = dist_bbox(node, 0, node->count()); }

  RNode* build(std::vector<int>& items, int left, int right, int height) {
    double N = right - left + 1, M = max_entries_;
    if (N <= M) {
      RNode* node = new_node(1, true);
      node->items.assign(items.begin() + left, items.begin() + right + 1);
      calc_bbox(node);
      return node;
    }
    if (!height) {
      height = static_cast<int>(std::ceil(v8_log(N) / v8_log(M)));
      M = std::ceil(N / std::pow(M, height - 1));
    }
    RNode* node = new_node(height, false);
    double N2 = std::ceil(N / M), N1 = N2 * std::ceil(std::sqrt(M));
    auto cx = [this](int a, int b) { return t_.compare_min_x(a, b); };
    auto cy = [this](int a, int b) { return t_.compare_min_y(a, b); };
    multi_select(items, left, right, N1, cx);
    for (double i = left; i <= right; i += N1) {
      double right2 = std::min(i + N1 - 1, static_cast<double>(right));
      multi_select(items, static_cast<int>(i), static_cast<int>(right2), N2, cy);
      for (double j = i; j <= right2; j += N2) {
        double right3 = std::min(j + N2 - 1, right2);
        node->nodes.push_back(build(items, static_cast<int>(j), static_cast<int>(right3), height - 1));
      }
    }
    calc_bbox(node);
    return node;
  }

  RNode* choose_subtree(const BBox& bbox, RNode* node, int level, std::vector<RNode*>& path) const {
    while (true) {
      path.push_back(node);
      if (node->leaf || static_cast<int>(path.size()) - 1 == level) break;
      double min_area = std::numeric_limits<double>::infinity(), min_enlargement = min_area;
      RNode* target = nullptr;
      for (RNode* child : node->nodes) {
        double area = bbox_area(child->bbox);
        double enlargement = enlarged_area(bbox, child->bbox) - area;
        if (enlargement < min_enlargement) {
          min_enlargement = enlargement;
          min_area = area < min_area ? area : min_area;
          target = child;
        } else if (enlargement == min_enlargement) {
          if (area < min_area) {
            min_area = area;
            target = child;
          }
        }
      }
      node = target;
    }
    return node;
  }

  void insert_item(int item, int level) {
    BBox bbox = t_.to_bbox(item);
    std::vector<RNode*> insert_path;
    RNode* node = choose_subtree(bbox, data_, level, insert_path);
    node->items.push_back(item);
    extend(node->bbox, bbox);
    while (level >= 0) {
      if (static_cast<int>(insert_path[level]->count()) > max_entries_) {
        split(insert_path, level);
        level--;
      } else {
        break;
      }
    }
    for (int i = level; i >= 0; i--) extend(insert_path[i]->bbox, bbox);
  }

  void split(std::vector<RNode*>& insert_path, int level) {
    RNode* node = insert_path[level];
    int M = static_cast<int>(node->count()), m = min_entries_;
    choose_split_axis(node, m, M);
    int split_index = choose_split_index(node, m, M);
    RNode* nn = new_node(node->height, node->leaf);
    if (node->leaf) {
      nn->items.assign(node->items.begin() + split_index, node->items.end());
      node->items.resize(split_index);
    } else {
      nn->nodes.assign(node->nodes.begin() + split_index, node->nodes.end());
      node->nodes.resize(split_index);
    }
    calc_bbox(node);
    calc_bbox(nn);
    if (level) insert_path[level - 1]->nodes.push_back(nn);
    else split_root(node, nn);
  }

  void split_root(RNode* node, RNode* nn) {
    data_ = new_node(node->height + 1, false);
    data_->nodes = {node, nn};
    calc_bbox(data_);
  }

  int choose_split_index(const RNode* node, int m, int M) const {
    double min_overlap = std::numeric_limits<double>::infinity(), min_area = min_overlap;
    int index = 0;
    for (int i = m; i <= M - m; i++) {
      BBox b1 = dist_bbox(node, 0, i), b2 = dist_bbox(node, i, M);
      double overlap = intersection_area(b1, b2), area = bbox_area(b1) + bbox_area(b2);
      if (overlap < min_overlap) {
        min_overlap = overlap;
        index = i;
        min_area = area < min_area ? area : min_area;
      } else if (overlap == min_overlap) {
        if (area < min_area) {
          min_area = area;
          index = i;
        }
      }
    }
    return index;
  }

  // Array.prototype.sort in V8 is stable, so std::stable_sort gives the same order
  void sort_children(RNode* node, bool by_x) const {
    if (node->leaf) {
      if (by_x) std::stable_sort(node->items.begin(), node->items.end(), [this](int a, int b) { return t_.compare_min_x(a, b) < 0; });
      else std::stable_sort(node->items.begin(), node->items.end(), [this](int a, int b) { return t_.compare_min_y(a, b) < 0; });
    } else {
      int k = by_x ? 0 : 1;
      std::stable_sort(node->nodes.begin(), node->nodes.end(),
                       [k](const RNode* a, const RNode* b) { return a->bbox[k] - b->bbox[k] < 0; });
    }
  }

  void choose_split_axis(RNode* node, int m, int M) const {
    double x_margin = all_dist_margin(node, m, M, true);
    double y_margin = all_dist_margin(node, m, M, false);
    if (x_margin < y_margin) sort_children(node, true);
  }

  double all_dist_margin(RNode* node, int m, int M, bool by_x) const {
    sort_children(node, by_x);
    BBox left = dist_bbox(node, 0, m), right = dist_bbox(node, M - m, M);
    double margin = bbox_margin(left) + bbox_margin(right);
    for (int i = m; i < M - m; i++) {
      extend(left, child_bbox(node, i));
      margin += bbox_margin(left);
    }
    for (int i = M - m - 1; i >= m; i--) {
      extend(right, child_bbox(node, i));
      margin += bbox_margin(right);
    }
    return margin;
  }

  void condense(std::vector<RNode*>& path) {
    for (int i = static_cast<int>(path.size()) - 1; i >= 0; i--) {
      if (path[i]->count() == 0) {
        if (i > 0) {
          auto& siblings = path[i - 1]->nodes;
          siblings.erase(std::find(siblings.begin(), siblings.end(), path[i]));
        } else {
          clear();
        }
      } else {
        calc_bbox(path[i]);
      }
    }
  }

  // sort so that items come in groups of n unsorted items, with groups sorted between each other
  template <class Cmp>
  static void multi_select(std::vector<int>& arr, int left, int right, double n, const Cmp& compare) {
    std::vector<double> stack = {static_cast<double>(left), static_cast<double>(right)};
    while (!stack.empty()) {
      double r = stack.back();
      stack.pop_back();
      double l = stack.back();
      stack.pop_back();
      if (r - l <= n) continue;
      double mid = l + std::ceil((r - l) / n / 2) * n;
      select(arr, static_cast<long>(l), static_cast<long>(r), static_cast<long>(mid), compare);
      stack.insert(stack.end(), {l, mid, mid, r});
    }
  }

  // Floyd-Rivest selection: the smallest k - left + 1 items between left and right come first
  template <class Cmp>
  static void select(std::vector<int>& arr, long left, long right, long k, const Cmp& compare) {
    while (right > left) {
      if (right - left > 600) {
        double n = right - left + 1, i = k - left + 1;
        double z = v8_log(n);
        double s = 0.5 * v8_exp(2 * z / 3);
        double sd = 0.5 * std::sqrt(z * s * (n - s) / n) * (i - n / 2 < 0 ? -1 : 1);
        long new_left = static_cast<long>(js_max(left, std::floor(k - i * s / n + sd)));
        long new_right = static_cast<long>(js_min(right, std::floor(k + (n - i) * s / n + sd)));
        select(arr, new_left, new_right, k, compare);
      }
      int t = arr[k];
      long i = left, j = right;
      std::swap(arr[left], arr[k]);
      if (compare(arr[right], t) > 0) std::swap(arr[left], arr[right]);
      while (i < j) {
        std::swap(arr[i], arr[j]);
        i++;
        j--;
        while (compare(arr[i], t) < 0) i++;
        while (compare(arr[j], t) > 0) j--;
      }
      if (compare(arr[left], t) == 0) {
        std::swap(arr[left], arr[j]);
      } else {
        j++;
        std::swap(arr[j], arr[right]);
      }
      if (j <= k) left = j + 1;
      if (k <= j) right = j - 1;
    }
  }
};

// ---- the concave hull ----------------------------------------------------------------------------------

struct ListNode {
  int p, prev, next;
  double minX = 0, minY = 0, maxX = 0, maxY = 0;
};

struct Hull {
  const std::vector<double>& x;
  const std::vector<double>& y;
  std::vector<ListNode> list;

  struct PointTraits {
    const Hull* h;
    BBox to_bbox(int i) const { return {h->x[i], h->y[i], h->x[i], h->y[i]}; }
    double compare_min_x(int a, int b) const { return h->x[a] - h->x[b]; }
    double compare_min_y(int a, int b) const { return h->y[a] - h->y[b]; }
  };
  struct SegTraits {
    const Hull* h;
    BBox to_bbox(int i) const {
      const ListNode& n = h->list[i];
      return {n.minX, n.minY, n.maxX, n.maxY};
    }
    double compare_min_x(int a, int b) const { return h->list[a].minX - h->list[b].minX; }
    double compare_min_y(int a, int b) const { return h->list[a].minY - h->list[b].minY; }
  };

  Hull(const std::vector<double>& x_, const std::vector<double>& y_) : x(x_), y(y_) {}

  double sq_dist(int a, int b) const {
    double dx = x[a] - x[b], dy = y[a] - y[b];
    return dx * dx + dy * dy;
  }

  // square distance from point p to the segment (p1, p2)
  double sq_seg_dist(int p, int p1, int p2) const {
    double px = x[p], py = y[p];
    double sx = x[p1], sy = y[p1], dx = x[p2] - sx, dy = y[p2] - sy;
    if (dx != 0 || dy != 0) {
      double t = ((px - sx) * dx + (py - sy) * dy) / (dx * dx + dy * dy);
      if (t > 1) {
        sx = x[p2];
        sy = y[p2];
      } else if (t > 0) {
        sx += dx * t;
        sy += dy * t;
      }
    }
    dx = px - sx;
    dy = py - sy;
    return dx * dx + dy * dy;
  }

  // segment to segment distance (Dan Sunday)
  static double sq_seg_seg_dist(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3) {
    double ux = x1 - x0, uy = y1 - y0, vx = x3 - x2, vy = y3 - y2, wx = x0 - x2, wy = y0 - y2;
    double a = ux * ux + uy * uy, b = ux * vx + uy * vy, c = vx * vx + vy * vy, d = ux * wx + uy * wy,
           e = vx * wx + vy * wy;
    double D = a * c - b * b;
    double sc, sN, tc, tN, sD = D, tD = D;
    if (D == 0) {
      sN = 0;
      sD = 1;
      tN = e;
      tD = c;
    } else {
      sN = b * e - c * d;
      tN = a * e - b * d;
      if (sN < 0) {
        sN = 0;
        tN = e;
        tD = c;
      } else if (sN > sD) {
        sN = sD;
        tN = e + b;
        tD = c;
      }
    }
    if (tN < 0.0) {
      tN = 0.0;
      if (-d < 0.0) sN = 0.0;
      else if (-d > a) sN = sD;
      else {
        sN = -d;
        sD = a;
      }
    } else if (tN > tD) {
      tN = tD;
      if ((-d + b) < 0.0) sN = 0;
      else if (-d + b > a) sN = sD;
      else {
        sN = -d + b;
        sD = a;
      }
    }
    sc = sN == 0 ? 0 : sN / sD;
    tc = tN == 0 ? 0 : tN / tD;
    double cx = (1 - sc) * x0 + sc * x1, cy = (1 - sc) * y0 + sc * y1;
    double cx2 = (1 - tc) * x2 + tc * x3, cy2 = (1 - tc) * y2 + tc * y3;
    double dx = cx2 - cx, dy = cy2 - cy;
    return dx * dx + dy * dy;
  }

  static bool inside(double ax, double ay, const BBox& b) {
    return ax >= b[0] && ax <= b[2] && ay >= b[1] && ay <= b[3];
  }

  // square distance from the segment (a, b) to a bounding box
  double sq_seg_box_dist(int a, int b, const BBox& bb) const {
    double ax = x[a], ay = y[a], bx = x[b], by = y[b];
    if (inside(ax, ay, bb) || inside(bx, by, bb)) return 0;
    double d1 = sq_seg_seg_dist(ax, ay, bx, by, bb[0], bb[1], bb[2], bb[1]);
    if (d1 == 0) return 0;
    double d2 = sq_seg_seg_dist(ax, ay, bx, by, bb[0], bb[1], bb[0], bb[3]);
    if (d2 == 0) return 0;
    double d3 = sq_seg_seg_dist(ax, ay, bx, by, bb[2], bb[1], bb[2], bb[3]);
    if (d3 == 0) return 0;
    double d4 = sq_seg_seg_dist(ax, ay, bx, by, bb[0], bb[3], bb[2], bb[3]);
    if (d4 == 0) return 0;
    return js_min(js_min(d1, d2), js_min(d3, d4));
  }

  int orient(int a, int b, int c) const { return orient_sign(x[a], y[a], x[b], y[b], x[c], y[c]); }

  // do the edges (p1, q1) and (p2, q2) intersect
  bool edges_intersect(int p1, int q1, int p2, int q2) const {
    return p1 != q2 && q1 != p2 && (orient(p1, q1, p2) > 0) != (orient(p1, q1, q2) > 0) &&
           (orient(p2, q2, p1) > 0) != (orient(p2, q2, q1) > 0);
  }

  // the edge (a, b) crosses no edge of the hull
  bool no_intersections(int a, int b, const RBush<SegTraits>& seg_tree) const {
    BBox bb = {js_min(x[a], x[b]), js_min(y[a], y[b]), js_max(x[a], x[b]), js_max(y[a], y[b])};
    for (int e : seg_tree.search(bb))
      if (edges_intersect(list[e].p, list[list[e].next].p, a, b)) return false;
    return true;
  }

  int update_bbox(int n) {
    ListNode& node = list[n];
    int p1 = node.p, p2 = list[node.next].p;
    node.minX = js_min(x[p1], x[p2]);
    node.minY = js_min(y[p1], y[p2]);
    node.maxX = js_max(x[p1], x[p2]);
    node.maxY = js_max(y[p1], y[p2]);
    return n;
  }

  int insert_node(int p, int prev) {
    int n = static_cast<int>(list.size());
    list.push_back(ListNode{p, n, n});
    if (prev >= 0) {
      list[n].next = list[prev].next;
      list[n].prev = prev;
      list[list[prev].next].prev = n;
      list[prev].next = n;
    }
    return n;
  }

  // ray casting
  bool point_in_polygon(int p, const std::vector<int>& vs) const {
    double px = x[p], py = y[p];
    bool in = false;
    for (std::size_t i = 0, j = vs.size() - 1; i < vs.size(); j = i++) {
      double xi = x[vs[i]], yi = y[vs[i]], xj = x[vs[j]], yj = y[vs[j]];
      bool intersect = ((yi > py) != (yj > py)) && (px < (xj - xi) * (py - yi) / (yj - yi) + xi);
      if (intersect) in = !in;
    }
    return in;
  }

  // monotone-convex-hull-2d, returning positions in `points`
  std::vector<int> monotone_hull(const std::vector<int>& points) const {
    int n = static_cast<int>(points.size());
    if (n < 3) {
      std::vector<int> result(n);
      for (int i = 0; i < n; ++i) result[i] = i;
      if (n == 2 && x[points[0]] == x[points[1]] && y[points[0]] == y[points[1]]) return {0};
      return result;
    }
    std::vector<int> sorted(n);
    for (int i = 0; i < n; ++i) sorted[i] = i;
    std::stable_sort(sorted.begin(), sorted.end(), [&](int a, int b) {
      double d = x[points[a]] - x[points[b]];
      if (d != 0) return d < 0;
      return y[points[a]] - y[points[b]] < 0;
    });
    std::vector<int> lower = {sorted[0], sorted[1]}, upper = {sorted[0], sorted[1]};
    for (int i = 2; i < n; ++i) {
      int idx = sorted[i], p = points[idx];
      while (lower.size() > 1 && orient(points[lower[lower.size() - 2]], points[lower.back()], p) <= 0) lower.pop_back();
      lower.push_back(idx);
      while (upper.size() > 1 && orient(points[upper[upper.size() - 2]], points[upper.back()], p) >= 0) upper.pop_back();
      upper.push_back(idx);
    }
    std::vector<int> result(lower);
    for (int j = static_cast<int>(upper.size()) - 2; j > 0; --j) result.push_back(upper[j]);
    return result;
  }

  // convex hull, first dropping the points inside the quadrilateral of the 4 extreme points
  std::vector<int> fast_convex_hull() const {
    int left = 0, top = 0, right = 0, bottom = 0;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
      if (x[i] < x[left]) left = i;
      if (x[i] > x[right]) right = i;
      if (y[i] < y[top]) top = i;
      if (y[i] > y[bottom]) bottom = i;
    }
    std::vector<int> cull = {left, top, right, bottom};
    std::vector<int> filtered(cull);
    for (int i = 0; i < static_cast<int>(x.size()); ++i)
      if (!point_in_polygon(i, cull)) filtered.push_back(i);
    std::vector<int> hull;
    for (int i : monotone_hull(filtered)) hull.push_back(filtered[i]);
    return hull;
  }

  // TinyQueue of points and R-tree nodes by distance
  struct QItem {
    RNode* node;
    int point;  // >= 0 for a point
    double dist;
  };
  struct TinyQueue {
    std::vector<QItem> data;
    static bool less(const QItem& a, const QItem& b) { return a.dist - b.dist < 0; }
    void push(const QItem& item) {
      data.push_back(item);
      std::size_t pos = data.size() - 1;
      while (pos > 0) {
        std::size_t parent = (pos - 1) / 2;
        if (less(data[pos], data[parent])) {
          std::swap(data[parent], data[pos]);
          pos = parent;
        } else {
          break;
        }
      }
    }
    QItem pop() {
      QItem top = data[0];
      data[0] = data.back();
      data.pop_back();
      std::size_t pos = 0, len = data.size();
      while (true) {
        std::size_t left = 2 * pos + 1, right = left + 1, min = pos;
        if (left < len && less(data[left], data[min])) min = left;
        if (right < len && less(data[right], data[min])) min = right;
        if (min == pos) break;
        std::swap(data[min], data[pos]);
        pos = min;
      }
      return top;
    }
  };

  // the point to flex the edge (b, c) inward to, searching the point tree in order of distance to the edge
  int find_candidate(const RBush<PointTraits>& tree, int a, int b, int c, int d, double max_dist,
                     const RBush<SegTraits>& seg_tree) const {
    TinyQueue queue;
    RNode* node = tree.root();
    while (node) {
      if (node->leaf) {
        for (int p : node->items) {
          double dist = sq_seg_dist(p, b, c);
          if (dist > max_dist) continue;
          queue.push({nullptr, p, dist});
        }
      } else {
        for (RNode* child : node->nodes) {
          double dist = sq_seg_box_dist(b, c, child->bbox);
          if (dist > max_dist) continue;
          queue.push({child, -1, dist});
        }
      }
      while (!queue.data.empty() && queue.data[0].point >= 0) {
        QItem item = queue.pop();
        int p = item.point;
        double d0 = sq_seg_dist(p, a, b), d1 = sq_seg_dist(p, c, d);
        if (item.dist < d0 && item.dist < d1 && no_intersections(b, p, seg_tree) && no_intersections(c, p, seg_tree))
          return p;
      }
      node = queue.data.empty() ? nullptr : queue.pop().node;
    }
    return -1;
  }

  std::vector<int> run(double concavity, double length_threshold) {
    std::vector<int> hull = fast_convex_hull();
    RBush<PointTraits> tree(PointTraits{this}, 16);
    std::vector<int> ids(x.size());
    for (std::size_t i = 0; i < ids.size(); ++i) ids[i] = static_cast<int>(i);
    tree.load(ids);

    std::deque<int> queue;
    int last = -1;
    for (int p : hull) {
      tree.remove(p);
      last = insert_node(p, last);
      queue.push_back(last);
    }
    RBush<SegTraits> seg_tree(SegTraits{this}, 16);
    for (int n : queue) seg_tree.insert(update_bbox(n));

    double sq_concavity = concavity * concavity;
    double sq_len_threshold = length_threshold * length_threshold;
    while (!queue.empty()) {
      int node = queue.front();
      queue.pop_front();
      int a = list[node].p, b = list[list[node].next].p;
      double sq_len = sq_dist(a, b);
      if (sq_len < sq_len_threshold) continue;
      double max_sq_len = sq_len / sq_concavity;
      int p = find_candidate(tree, list[list[node].prev].p, a, b, list[list[list[node].next].next].p, max_sq_len, seg_tree);
      if (p >= 0 && js_min(sq_dist(p, a), sq_dist(p, b)) <= max_sq_len) {
        queue.push_back(node);
        queue.push_back(insert_node(p, node));
        tree.remove(p);
        seg_tree.remove(node);
        seg_tree.insert(update_bbox(node));
        seg_tree.insert(update_bbox(list[node].next));
      }
    }
    std::vector<int> ring;
    int node = last;
    do {
      ring.push_back(list[node].p);
      node = list[node].next;
    } while (node != last);
    ring.push_back(list[node].p);
    return ring;
  }
};

// a number as jsonlite::toJSON(digits = 4) writes it and JavaScript reads it back
double json_number(double val) {
  char buf[64];
  if (std::fabs(val) < 2147483647 && std::fabs(val) > 1e-5) {
    // modp_dtoa2(val, buf, 4)
    static const double pow10_4 = 10000;
    bool neg = val < 0;
    double value = neg ? -val : val;
    int whole = static_cast<int>(value);
    double tmp = (value - whole) * pow10_4;
    uint32_t frac = static_cast<uint32_t>(tmp);
    double diff = tmp - frac;
    if (diff > 0.5) {
      ++frac;
      if (frac >= pow10_4) {
        frac = 0;
        ++whole;
      }
    } else if (diff == 0.5 && (frac & 1)) {
      ++frac;
      if (frac >= pow10_4) {
        frac = 0;
        ++whole;
      }
    }
    int count = 4;
    while (count > 0 && (frac % 10) == 0) {
      count--;
      frac /= 10;
    }
    if (count > 0) std::snprintf(buf, sizeof buf, "%s%d.%0*u", neg ? "-" : "", whole, count, frac);
    else std::snprintf(buf, sizeof buf, "%s%d", neg ? "-" : "", whole);
  } else {
    int decimals = static_cast<int>(std::ceil(std::fmin(17, std::fmax(1, std::log10(std::fabs(val))) + 4)));
    std::snprintf(buf, sizeof buf, "%.*g", decimals, val);
  }
  return std::strtod(buf, nullptr);
}

// a number as R's sprintf("%s") writes it (15 significant digits) and JavaScript reads it back
double r_string_number(double val) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.15g", val);
  return std::strtod(buf, nullptr);
}

}  // namespace

std::vector<std::array<double, 2>> concave_hull(const std::vector<double>& x, const std::vector<double>& y,
                                                double concavity, double length_threshold) {
  if (x.size() != y.size()) throw std::invalid_argument("concave_hull: x and y differ in length");
  if (x.empty()) throw std::invalid_argument("concave_hull: no points");
  for (std::size_t i = 0; i < x.size(); ++i)
    if (!std::isfinite(x[i]) || !std::isfinite(y[i])) throw std::invalid_argument("concave_hull: non-finite coordinates");
  // concavity = Math.max(0, concavity); lengthThreshold = lengthThreshold || 0
  concavity = js_max(0, concavity);
  if (std::isnan(length_threshold) || length_threshold == 0) length_threshold = 0;
  Hull hull(x, y);
  std::vector<int> ring = hull.run(concavity, length_threshold);
  std::vector<std::array<double, 2>> out(ring.size());
  for (std::size_t i = 0; i < ring.size(); ++i) out[i] = {x[ring[i]], y[ring[i]]};
  return out;
}

std::vector<std::array<double, 2>> concaveman_r(const std::vector<double>& x, const std::vector<double>& y,
                                                double concavity, double length_threshold) {
  std::vector<double> xr(x.size()), yr(y.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    xr[i] = json_number(x[i]);
    yr[i] = json_number(y[i]);
  }
  return concave_hull(xr, yr, r_string_number(concavity), r_string_number(length_threshold));
}

namespace detail {
double js_exp(double x) { return v8_exp(x); }
double js_log(double x) { return v8_log(x); }
int orient3_sign(double ax, double ay, double bx, double by, double cx, double cy) {
  return orient_sign(ax, ay, bx, by, cx, cy);
}
double json4(double v) { return json_number(v); }
}  // namespace detail

}  // namespace concaveman

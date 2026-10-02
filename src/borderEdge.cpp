#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

#include "spicyCore.h"

namespace {

using spicy::Pt;

struct Seg {
  double x0, y0, x1, y1;
};

// One Sutherland-Hodgman step: keep the part of `in` to the left of the
// directed line a -> b. The subject polygon may be non-convex; the resulting
// polygon can contain zero-width slivers, which do not change its area.
void clipHalfPlane(const std::vector<Pt>& in, std::vector<Pt>& out,
                   const Pt& a, const Pt& b) {
  out.clear();
  const std::size_t n = in.size();
  if (n == 0) return;
  const double ex = b.x - a.x, ey = b.y - a.y;
  Pt prev = in[n - 1];
  double sPrev = ex * (prev.y - a.y) - ey * (prev.x - a.x);
  for (std::size_t i = 0; i < n; ++i) {
    const Pt& cur = in[i];
    const double sCur = ex * (cur.y - a.y) - ey * (cur.x - a.x);
    if ((sCur >= 0) != (sPrev >= 0)) {
      const double t = sPrev / (sPrev - sCur);
      out.push_back({prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y)});
    }
    if (sCur >= 0) out.push_back(cur);
    prev = cur;
    sPrev = sCur;
  }
}

bool isConvex(const std::vector<Pt>& p) {
  const std::size_t n = p.size();
  if (n < 3) return false;
  for (std::size_t i = 0; i < n; ++i) {
    const Pt &a = p[i], &b = p[(i + 1) % n], &c = p[(i + 2) % n];
    if ((b.x - a.x) * (c.y - b.y) - (b.y - a.y) * (c.x - b.x) < 0) return false;
  }
  return true;
}

double signedArea(const std::vector<Pt>& p) {
  const std::size_t n = p.size();
  double a = 0;
  for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
    a += (p[j].x + p[i].x) * (p[i].y - p[j].y);
  }
  return a / 2;
}

}  // namespace

// Area of the intersection of the window with an npoly-gon approximating the
// disc of radius r around each centre. Reproduces spatstat.geom's
// area(intersect.owin(disc(r, npoly = npoly), W)). `rings` are the window's
// boundary polygons in spatstat orientation (outer anticlockwise, holes
// clockwise), so their clipped signed areas sum to the intersection area.
void spicy::discAreas(const double* x, const double* y, std::size_t n, double r, int npoly,
                       const std::vector<spicy::Ring>& rings, double* out) {
  // Disc vertices exactly as spatstat.geom::disc() places them.
  std::vector<double> dx(npoly), dy(npoly);
  const double by = 2 * M_PI / npoly;
  double minX = 0, maxX = 0, minY = 0, maxY = 0;
  for (int k = 0; k < npoly; ++k) {
    const double theta = k * by;
    dx[k] = r * std::cos(theta);
    dy[k] = r * std::sin(theta);
    minX = std::min(minX, dx[k]); maxX = std::max(maxX, dx[k]);
    minY = std::min(minY, dy[k]); maxY = std::max(maxY, dy[k]);
  }

  const std::vector<spicy::Ring>& win = rings;

  std::vector<Pt> a, b, disc(npoly);

  // Convex window (the convex hull and square windows): the window is the
  // intersection of its edges' half-planes, so clip the disc by only those
  // edges whose line comes within r of the centre. Every vertex of the
  // inscribed npoly-gon lies within r of the centre, so the others cannot cut it.
  if (win.size() == 1 && isConvex(win[0])) {
    const std::vector<Pt>& w = win[0];
    const std::size_t m = w.size();
    std::vector<double> len(m);
    for (std::size_t j = 0; j < m; ++j) {
      const Pt &p = w[j], &q = w[(j + 1) % m];
      len[j] = std::hypot(q.x - p.x, q.y - p.y);
    }
    // Area of the whole npoly-gon, the answer when no window edge cuts it.
    double fullArea = 0;
    for (int k = 0, l = npoly - 1; k < npoly; l = k++) fullArea += (dx[l] + dx[k]) * (dy[k] - dy[l]);
    fullArea /= 2;
    std::vector<double> dist(m), sides;
    std::vector<Seg> cuts;
    for (std::size_t c = 0; c < n; ++c) {
      const double cx = x[c], cy = y[c];
      int nCut = 0, cut = -1;
      for (std::size_t j = 0; j < m; ++j) {
        const Pt &p = w[j], &q = w[(j + 1) % m];
        dist[j] = len[j] == 0 ? r : ((q.x - p.x) * (cy - p.y) - (q.y - p.y) * (cx - p.x)) / len[j];
        if (dist[j] < r) { ++nCut; cut = static_cast<int>(j); }
      }
      if (nCut == 0) {
        out[c] = fullArea;
        continue;
      }
      if (nCut == 1 && dist[cut] >= 0) {
        // One edge cuts the disc: subtract the cap beyond it. The vertices
        // outside the edge form one run around the outward normal, so only
        // that run is visited.
        const Pt &p = w[cut], &q = w[(cut + 1) % m];
        const double ex = q.x - p.x, ey = q.y - p.y;
        auto side = [&](int k) { return ex * (dy[k] + cy - p.y) - ey * (dx[k] + cx - p.x); };
        auto idx = [&](int k) { return ((k % npoly) + npoly) % npoly; };
        const double normal = std::atan2(ey, ex) - M_PI / 2;
        const double half = std::acos(std::min(dist[cut] / r, 1.0));
        const int kc = static_cast<int>(std::lround(normal / by));
        const int span = static_cast<int>(std::ceil(half / by)) + 2;
        int ta = span + 1, tb = -span - 1;
        for (int t = -span; t <= span; ++t) {
          if (side(idx(kc + t)) < 0) { ta = std::min(ta, t); tb = std::max(tb, t); }
        }
        if (ta > tb) {
          out[c] = fullArea;
          continue;
        }
        // Cap polygon, in coordinates relative to the centre: the crossing
        // into the run, the outside vertices, and the crossing back out.
        auto crossing = [&](int k0, int k1) {
          const double s0 = side(k0), s1 = side(k1), t = s0 / (s0 - s1);
          return Pt{dx[k0] + t * (dx[k1] - dx[k0]), dy[k0] + t * (dy[k1] - dy[k0])};
        };
        a.clear();
        a.push_back(crossing(idx(kc + ta - 1), idx(kc + ta)));
        for (int t = ta; t <= tb; ++t) a.push_back({dx[idx(kc + t)], dy[idx(kc + t)]});
        a.push_back(crossing(idx(kc + tb), idx(kc + tb + 1)));
        out[c] = std::max(fullArea - signedArea(a), 0.0);
        continue;
      }
      // Several edges: integrate over the boundary of the intersection
      // (Green's theorem), in coordinates relative to the centre. It consists
      // of the parts of the polygon's edges inside the window and the parts of
      // the cutting window edges inside the polygon.
      cuts.clear();
      for (std::size_t j = 0; j < m; ++j)
        if (dist[j] < r) cuts.push_back({w[j].x - cx, w[j].y - cy, w[(j + 1) % m].x - cx, w[(j + 1) % m].y - cy});
      const std::size_t nc = cuts.size();
      sides.resize(static_cast<std::size_t>(npoly) * nc);
      for (int k = 0; k < npoly; ++k)
        for (std::size_t j = 0; j < nc; ++j) {
          const Seg& e = cuts[j];
          sides[k * nc + j] = (e.x1 - e.x0) * (dy[k] - e.y0) - (e.y1 - e.y0) * (dx[k] - e.x0);
        }
      double area = 0;
      for (int k = 0; k < npoly; ++k) {
        const int l = (k + 1) % npoly;
        double lo = 0, hi = 1;
        for (std::size_t j = 0; j < nc && lo < hi; ++j) {
          const double s0 = sides[k * nc + j], s1 = sides[l * nc + j];
          if (s0 >= 0 && s1 >= 0) continue;
          if (s0 < 0 && s1 < 0) { hi = lo; break; }
          const double t = s0 / (s0 - s1);
          if (s0 < 0) lo = std::max(lo, t); else hi = std::min(hi, t);
        }
        if (lo >= hi) continue;
        const double ax = dx[k] + lo * (dx[l] - dx[k]), ay = dy[k] + lo * (dy[l] - dy[k]);
        const double bx = dx[k] + hi * (dx[l] - dx[k]), by_ = dy[k] + hi * (dy[l] - dy[k]);
        area += ax * by_ - bx * ay;
      }
      for (const Seg& e : cuts) {
        double lo = 0, hi = 1;
        for (int k = 0; k < npoly && lo < hi; ++k) {
          const int l = (k + 1) % npoly;
          const double ex = dx[l] - dx[k], ey = dy[l] - dy[k];
          const double s0 = ex * (e.y0 - dy[k]) - ey * (e.x0 - dx[k]);
          const double s1 = ex * (e.y1 - dy[k]) - ey * (e.x1 - dx[k]);
          if (s0 >= 0 && s1 >= 0) continue;
          if (s0 < 0 && s1 < 0) { hi = lo; break; }
          const double t = s0 / (s0 - s1);
          if (s0 < 0) lo = std::max(lo, t); else hi = std::min(hi, t);
        }
        if (lo >= hi) continue;
        const double ax = e.x0 + lo * (e.x1 - e.x0), ay = e.y0 + lo * (e.y1 - e.y0);
        const double bx = e.x0 + hi * (e.x1 - e.x0), by_ = e.y0 + hi * (e.y1 - e.y0);
        area += ax * by_ - bx * ay;
      }
      out[c] = std::max(area / 2, 0.0);
    }
    return;
  }

  for (std::size_t c = 0; c < n; ++c) {
    const double cx = x[c], cy = y[c];
    for (int k = 0; k < npoly; ++k) disc[k] = {dx[k] + cx, dy[k] + cy};
    // The disc's bounding box contains the disc, so clipping to it first is
    // exact and leaves only the few window vertices near this centre.
    const Pt box[4] = {{cx + minX, cy + minY}, {cx + maxX, cy + minY},
                       {cx + maxX, cy + maxY}, {cx + minX, cy + maxY}};
    double area = 0;
    for (const std::vector<Pt>& ring : win) {
      a = ring;
      for (int k = 0; k < 4 && !a.empty(); ++k) {
        clipHalfPlane(a, b, box[k], box[(k + 1) % 4]);
        a.swap(b);
      }
      for (int k = 0; k < npoly && !a.empty(); ++k) {
        clipHalfPlane(a, b, disc[k], disc[(k + 1) % npoly]);
        a.swap(b);
      }
      if (a.size() >= 3) area += signedArea(a);
    }
    out[c] = std::max(area, 0.0);
  }
}

// [[Rcpp::export]]
Rcpp::NumericVector discWindowArea(Rcpp::NumericVector x, Rcpp::NumericVector y,
                                   double r, int npoly, Rcpp::List rings) {
  std::vector<spicy::Ring> win(rings.size());
  for (R_xlen_t i = 0; i < rings.size(); ++i) {
    Rcpp::List ring = rings[i];
    Rcpp::NumericVector rx = ring["x"], ry = ring["y"];
    for (R_xlen_t j = 0; j < rx.size(); ++j) win[i].push_back({rx[j], ry[j]});
  }
  Rcpp::NumericVector out(x.size());
  spicy::discAreas(x.begin(), y.begin(), x.size(), r, npoly, win, out.begin());
  return out;
}

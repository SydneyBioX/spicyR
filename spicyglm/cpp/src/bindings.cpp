#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <limits>

#include "spicyglm/core.hpp"

namespace py = pybind11;
using namespace spicyglm;

namespace {

using IntArray = py::array_t<int, py::array::c_style | py::array::forcecast>;
using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

// NumPy -> std::vector with one memcpy (pybind11's list caster converts element by element)
template <typename T>
std::vector<T> vec(const py::array_t<T, py::array::c_style | py::array::forcecast>& a) {
  if (a.ndim() != 1) throw std::invalid_argument("expected a one-dimensional array");
  return std::vector<T>(a.data(), a.data() + a.size());
}

template <typename T>
py::array_t<T> to_array(const std::vector<T>& v) {
  py::array_t<T> a(v.size());
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

Window parse_window(const std::string& w) {
  if (w == "convex") return Window::Convex;
  if (w == "rectangle") return Window::Rectangle;
  throw std::invalid_argument("window must be 'convex' or 'rectangle'");
}

py::dict pair_fit_dict(const PairFit& p) {
  py::dict out;
  out["beta"] = py::make_tuple(p.fit.beta[0], p.fit.beta[1]);
  out["v_hat"] = p.naive ? p.v_naive : p.cr2.v_hat;
  out["df"] = p.naive ? std::numeric_limits<double>::infinity() : p.cr2.df;  // naive: normal reference
  if (!p.has_diagnostics) return out;
  const PairDiagnostics& D = p.diagnostics;
  py::dict patient, image;
  patient["cluster_id"] = to_array(p.cr2.cluster_id);
  patient["group"] = to_array(p.cr2.group);
  patient["n_i"] = to_array(D.n_i);
  patient["T_i"] = to_array(D.T);
  patient["S_g"] = py::make_tuple(D.S_leverage[0], D.S_leverage[1]);
  patient["l_i"] = to_array(D.l);
  patient["raw_residual_sum"] = to_array(D.raw_sum);
  patient["adjusted_residual_sum"] = to_array(D.adjusted_sum);
  patient["e_i"] = to_array(p.cr2.e);
  patient["influence_i"] = to_array(D.influence);
  patient["y_i"] = to_array(D.y);
  patient["d_i"] = to_array(D.d);
  patient["delta_i"] = to_array(D.delta);
  image["cluster"] = to_array(D.image_cluster);
  image["image_id"] = to_array(D.image_id);
  image["n_ij"] = to_array(D.n_ij);
  image["density_ij"] = to_array(D.density_ij);
  image["l_ij"] = to_array(D.l_ij);
  image["l_ij_group_share"] = to_array(D.l_ij_group_share);
  image["raw_residual_sum_ij"] = to_array(D.raw_sum_ij);
  image["adjusted_residual_sum_ij"] = to_array(D.adjusted_sum_ij);
  image["e_ij"] = to_array(D.e_ij);
  image["e_ij_share_within_patient"] = to_array(D.e_share_ij);
  image["influence_ij"] = to_array(D.influence_ij);
  out["patient"] = patient;
  out["image"] = image;
  return out;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "spicyGLM numeric core";

  m.def("window_area", [](const DoubleArray& x, const DoubleArray& y, const std::string& window) {
    if (x.size() != y.size()) throw std::invalid_argument("x and y differ in length");
    return window_area(x.data(), y.data(), x.size(), parse_window(window));
  }, py::arg("x"), py::arg("y"), py::arg("window"));

  m.def("knn_indices", [](const DoubleArray& x, const DoubleArray& y, const IntArray& image_offsets, int k,
                          int n_threads) {
    std::vector<double> xv = vec(x), yv = vec(y);
    std::vector<int> offsets = vec(image_offsets), idx;
    {
      py::gil_scoped_release release;
      idx = knn_indices(xv, yv, offsets, k, n_threads);
    }
    py::array_t<int> a({static_cast<py::ssize_t>(xv.size()), static_cast<py::ssize_t>(k)});
    std::copy(idx.begin(), idx.end(), a.mutable_data());
    return a;
  }, py::arg("x"), py::arg("y"), py::arg("image_offsets"), py::arg("k"), py::arg("n_threads") = 1);

  py::class_<Dataset>(m, "Dataset")
      .def(py::init([](const DoubleArray& x, const DoubleArray& y, const IntArray& cell_type,
                       const IntArray& image_offsets, int n_types) {
             return Dataset(vec(x), vec(y), vec(cell_type), vec(image_offsets), n_types);
           }),
           py::arg("x"), py::arg("y"), py::arg("cell_type"), py::arg("image_offsets"), py::arg("n_types"))
      .def("image_areas", [](const Dataset& d, const std::string& window) {
        return to_array(d.image_areas(parse_window(window)));
      }, py::arg("window"))
      .def("build_radius_index", &Dataset::build_radius_index, py::arg("r"),
           py::call_guard<py::gil_scoped_release>())
      .def("build_knn", &Dataset::build_knn, py::arg("k"), py::arg("n_threads") = 1,
           py::call_guard<py::gil_scoped_release>())
      .def("poisson_model_data", [](const Dataset& d, const DoubleArray& image_area, int from, int to) {
        std::vector<double> areas = vec(image_area);
        ModelData md;
        {
          py::gil_scoped_release release;
          md = d.poisson_model_data(areas, from, to);
        }
        py::dict out;
        out["row"] = to_array(md.row);
        out["image"] = to_array(md.image);
        out["n"] = to_array(md.n);
        out["density"] = to_array(md.density);
        return out;
      }, py::arg("image_area"), py::arg("from_"), py::arg("to"))
      .def("binomial_model_data", [](const Dataset& d, int from, int to) {
        BinomialModelData md;
        {
          py::gil_scoped_release release;
          md = d.binomial_model_data(from, to);
        }
        py::dict out;
        out["row"] = to_array(md.row);
        out["image"] = to_array(md.image);
        out["n"] = to_array(md.n);
        out["p0"] = to_array(md.p0);
        return out;
      }, py::arg("from_"), py::arg("to"));

  m.def("fit_poisson", [](const IntArray& n, const DoubleArray& density, const IntArray& group,
                          const std::string& estimator) {
    GlmFit f = fit_poisson(vec(n), vec(density), vec(group), estimator);
    py::dict out;
    out["beta"] = py::make_tuple(f.beta[0], f.beta[1]);
    out["mu"] = to_array(f.mu);
    return out;
  }, py::arg("n"), py::arg("density"), py::arg("group"), py::arg("estimator"));

  m.def("fit_binomial", [](const IntArray& n, int k, const DoubleArray& p0, const IntArray& group,
                           const std::string& estimator) {
    GlmFit f = fit_binomial(vec(n), k, vec(p0), vec(group), estimator);
    py::dict out;
    out["beta"] = py::make_tuple(f.beta[0], f.beta[1]);
    out["mu"] = to_array(f.mu);
    return out;
  }, py::arg("n"), py::arg("k"), py::arg("p0"), py::arg("group"), py::arg("estimator"));

  m.def("cr2_wald", [](const IntArray& cluster, const IntArray& image, const IntArray& group,
                       const DoubleArray& var, const DoubleArray& resid) {
    CR2Result r = cr2_wald(vec(cluster), vec(image), vec(group), vec(var), vec(resid));
    py::dict out;
    out["S"] = py::make_tuple(r.S[0], r.S[1]);
    out["cluster_id"] = to_array(r.cluster_id);
    out["group"] = to_array(r.group);
    out["e"] = to_array(r.e);
    out["image_id"] = r.image_id;
    out["adjusted_image_sum"] = r.adjusted_image_sum;
    out["raw_image_sum"] = r.raw_image_sum;
    out["v_hat"] = r.v_hat;
    out["df"] = r.df;
    return out;
  }, py::arg("cluster"), py::arg("image"), py::arg("group"), py::arg("var"), py::arg("resid"));

  m.def("fit_pair_poisson", [](const IntArray& cluster, const IntArray& image, const IntArray& group,
                               const IntArray& n, const DoubleArray& density, const std::string& estimator,
                               const std::string& variance, bool diagnostics) {
    std::vector<int> cv = vec(cluster), iv = vec(image), gv = vec(group), nv = vec(n);
    std::vector<double> dv = vec(density);
    PairFit p;
    {
      py::gil_scoped_release release;
      p = fit_pair_poisson(cv, iv, gv, nv, dv, estimator, variance, diagnostics);
    }
    return pair_fit_dict(p);
  }, py::arg("cluster"), py::arg("image"), py::arg("group"), py::arg("n"), py::arg("density"),
     py::arg("estimator"), py::arg("variance") = "fast", py::arg("diagnostics") = false);

  m.def("fit_pair_binomial", [](const IntArray& cluster, const IntArray& image, const IntArray& group,
                                const IntArray& n, int k, const DoubleArray& p0, const std::string& estimator,
                                const std::string& variance) {
    std::vector<int> cv = vec(cluster), iv = vec(image), gv = vec(group), nv = vec(n);
    std::vector<double> pv = vec(p0);
    PairFit p;
    {
      py::gil_scoped_release release;
      p = fit_pair_binomial(cv, iv, gv, nv, k, pv, estimator, variance);
    }
    return pair_fit_dict(p);
  }, py::arg("cluster"), py::arg("image"), py::arg("group"), py::arg("n"), py::arg("k"),
     py::arg("p0"), py::arg("estimator"), py::arg("variance") = "fast");
}

#ifndef OVMWRAP_IMPL
#define OVMWRAP_IMPL

#include "ovmwrap.h"

namespace {
const int hex2tet[32] = {0, 1, 3, 4, 1, 2, 0, 5, 2, 3, 1, 6, 3, 0, 2, 7,
                         4, 7, 5, 0, 7, 6, 4, 3, 6, 5, 7, 2, 5, 4, 6, 1};

bool check_vector_increase(Eigen::VectorXi &v) {
  for (int i = 1; i < v.size(); ++i) {
    if (v(i) < v(i - 1)) {
      return false;
    }
  }
  return true;
}
}; // namespace

inline std::vector<OpenVolumeMesh::VertexHandle> opposite_vertex_in_cell(
    const OpenVolumeMesh::HexahedralMeshTopologyKernel &mesh,
    const OpenVolumeMesh::CellHandle &cell,
    const OpenVolumeMesh::HalfFaceHandle &hf,
    const std::vector<OpenVolumeMesh::VertexHandle> &v) {
  auto opposite_hf = mesh.opposite_halfface_handle_in_cell(hf, cell);
  auto range_t = mesh.halfface_vertices(opposite_hf);
  std::set<OpenVolumeMesh::VertexHandle> topset(range_t.first, range_t.second);
  std::vector<OpenVolumeMesh::VertexHandle> res;
  for (auto x : v) {
    for (auto vv_it = mesh.vv_iter(x); vv_it.valid(); ++vv_it) {
      if (topset.count(*vv_it)) {
        res.push_back(*vv_it);
        break;
      }
    }
  }
  return res;
}

inline void transform_hex_to_matrix(
    Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXi &b,
    Eigen::MatrixXd &bc, const OpenVolumeMesh::GeometricHexahedralMeshV3d &mesh,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        &fixed) {
  using namespace Eigen;
  // using namespace OpenVolumeMesh;
  int _vnum = mesh.n_vertices();
  int _hnum = mesh.n_cells();
  int _tnum = _hnum * 8;
  V.resize(_vnum, 3);
  T.resize(_tnum, 4);

  b.resize(fixed.size());
  bc.resize(fixed.size(), 3);

  int bnumber = 0;
  for (const auto &_vh : mesh.vertices()) {
    const int _vn = _vh.idx();
    const auto &_p = mesh.vertex(_vh);
    V.row(_vn) = Vector3d{_p[0], _p[1], _p[2]};
    if (fixed.find(_vh) != fixed.end()) {
      b(bnumber) = _vh;
      bc.row(bnumber) = Vector3d{fixed[_vh][0], fixed[_vh][1], fixed[_vh][2]};
      bnumber++;
    }
  }

#ifndef NDEBUG
  assert(check_vector_increase);
#endif

  for (const auto &_ch : mesh.cells()) {
    const int _hn = _ch.idx();
    const auto _xb = mesh.xback_halfface(_ch);
    const auto _xb_range = mesh.halfface_vertices(_xb);
    std::vector<OpenVolumeMesh::VertexHandle> _vertices(_xb_range.first,
                                                        _xb_range.second);
    const auto _xb_range_opposite =
        opposite_vertex_in_cell(mesh, _ch, _xb, _vertices);
    _vertices.insert(_vertices.end(), _xb_range_opposite.begin(),
                     _xb_range_opposite.end());
    for (size_t i = 0; i < 8; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        int _tn = _hn * 8 + i;
        T(_tn, j) = _vertices.at(hex2tet[i * 4 + j]).idx();
      }
    }
  }
}

inline void transform_hex_to_matrix(
    Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXi &b,
    Eigen::MatrixXd &bc, Eigen::MatrixXi &surface,
    const OpenVolumeMesh::GeometricHexahedralMeshV3d &mesh,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        &fixed) {
  using namespace Eigen;
  // using namespace OpenVolumeMesh;
  int _vnum = mesh.n_vertices();
  int _hnum = mesh.n_cells();
  int _tnum = _hnum * 8;
  V.resize(_vnum, 3);
  T.resize(_tnum, 4);

  b.resize(fixed.size());
  bc.resize(fixed.size(), 3);

  int bnumber = 0;
  for (const auto &_vh : mesh.vertices()) {
    const int _vn = _vh.idx();
    const auto &_p = mesh.vertex(_vh);
    V.row(_vn) = Vector3d{_p[0], _p[1], _p[2]};
    if (fixed.find(_vh) != fixed.end()) {
      b(bnumber) = _vh;
      bc.row(bnumber) = Vector3d{fixed[_vh][0], fixed[_vh][1], fixed[_vh][2]};
      bnumber++;
    }
  }

#ifndef NDEBUG
  assert(check_vector_increase);
#endif

  for (const auto &_ch : mesh.cells()) {
    const int _hn = _ch.idx();
    const auto _xb = mesh.xback_halfface(_ch);
    const auto _xb_range = mesh.halfface_vertices(_xb);
    std::vector<OpenVolumeMesh::VertexHandle> _vertices(_xb_range.first,
                                                        _xb_range.second);
    const auto _xb_range_opposite =
        opposite_vertex_in_cell(mesh, _ch, _xb, _vertices);
    _vertices.insert(_vertices.end(), _xb_range_opposite.begin(),
                     _xb_range_opposite.end());
    for (size_t i = 0; i < 8; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        int _tn = _hn * 8 + i;
        T(_tn, j) = _vertices.at(hex2tet[i * 4 + j]).idx();
      }
    }
  }

  int _sfn = 0;
  for (const auto &_hfh : mesh.halffaces()) {
    if (mesh.is_boundary(_hfh)) {
      _sfn += 1;
    }
  }
  surface.resize(_sfn, 4);
  _sfn = 0;
  for (const auto &_hfh : mesh.halffaces()) {
    if (mesh.is_boundary(_hfh)) {
      auto _v_range = mesh.halfface_vertices(_hfh);
      std::vector<OpenVolumeMesh::VertexHandle> _vertices(_v_range.first,
                                                          _v_range.second);
      surface.row(_sfn) = Vector4i{_vertices[0].idx(), _vertices[1].idx(),
                                   _vertices[2].idx(), _vertices[3].idx()};
      _sfn += 1;
    }
  }
}

inline void
transform_matrix_to_hex(const Eigen::MatrixXd &V,
                        OpenVolumeMesh::GeometricHexahedralMeshV3d &mesh) {
  int _vnum = V.rows();
  for (int i = 0; i < _vnum; ++i) {
    mesh.set_vertex(OpenVolumeMesh::VertexHandle(i),
                    OpenVolumeMesh::Vec3d(V(i, 0), V(i, 1), V(i, 2)));
  }
}

inline void write_matrix_to_vtk_tet(const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &T,
                                    std::ostream &stream) {
  stream << "# vtk DataFile Version 3.0\ngenerate by dongmo\nASCII\nDATASET "
            "UNSTRUCTURED_GRID"
         << std::endl;
  stream << "POINTS " << V.rows() << " double" << std::endl;
  for (int i = 0; i < V.rows(); ++i) {
    stream << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << std::endl;
  }
  stream << "CELLS " << T.rows() << " " << 5 * T.rows() << std::endl;
  for (int i = 0; i < T.rows(); ++i) {
    stream << "4 " << T(i, 0) << " " << T(i, 1) << " " << T(i, 2) << " "
           << T(i, 3) << std::endl;
  }
  stream << "CELL_TYPES " << T.rows() << std::endl;
  for (int i = 0; i < T.rows(); ++i) {
    stream << "10 " << std::endl;
  }
  return;
}
#endif
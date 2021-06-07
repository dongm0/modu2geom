#include "msqoptimizer.h"
//#include "ovmwrap.h"

using namespace Mesquite;
/*
ArrayMesh
MsqOperator::OVMmesh2MSQmesh(OpenVolumeMesh::GeometricHexahedralMeshV3d _ovm) {
    std::vector<double> coords;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        auto c = _ovm.vertex(*v_it);
        coords.push_back(c[0]);
        coords.push_back(c[1]);
        coords.push_back(c[2]);
    }
    std::vector<int> fixedflag;
    std::vector<int> connection;
    //可能销毁，想想办法
    ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(),
_ovm.n_cells(), HEXAHEDRON, connection.data());
}
*/
namespace {
std::vector<OpenVolumeMesh::VertexHandle> opposite_vertex_in_cell1(
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
} // namespace

void MsqOperator::Ovm2MsqOut(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                             std::string outname) {
  std::vector<double> coords;
  std::map<OpenVolumeMesh::VertexHandle, int> v_idx;
  int i = 0;
  for (auto v_it = _ovm.vertices_begin(); v_it != _ovm.vertices_end(); ++v_it) {
    auto c = _ovm.vertex(*v_it);
    v_idx[*v_it] = i++;
    coords.push_back(c[0]);
    coords.push_back(c[1]);
    coords.push_back(c[2]);
  }
  // std::vector<bool> fixedflag;
  std::vector<int> connection;
  std::vector<OpenVolumeMesh::VertexHandle> _xb_vertices(4);
  std::vector<OpenVolumeMesh::VertexHandle> _xf_vertices(4);
  for (auto _ch : _ovm.cells()) {
    auto xf = _ovm.xfront_halfface(_ch);
    auto xb = _ovm.xback_halfface(_ch);
    auto _xb_range = _ovm.halfface_vertices(xb);
    auto _xf_range = _ovm.halfface_vertices(xf);
    //_xb_vertices.clear();
    int _num = 0;
    for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
      _xb_vertices.at(_num++) = *_vh;
    }
    _xf_vertices = opposite_vertex_in_cell1(_ovm, _ch, xb, _xb_vertices);
    connection.push_back(v_idx[_xb_vertices[0]]);
    connection.push_back(v_idx[_xb_vertices[1]]);
    connection.push_back(v_idx[_xb_vertices[2]]);
    connection.push_back(v_idx[_xb_vertices[3]]);
    connection.push_back(v_idx[_xf_vertices[0]]);
    connection.push_back(v_idx[_xf_vertices[1]]);
    connection.push_back(v_idx[_xf_vertices[2]]);
    connection.push_back(v_idx[_xf_vertices[3]]);
  }
  // std::shared_ptr<bool> fixed = std::make_shared(new
  // bool(_ovm.n_vertices()));
  bool *fixed = new bool[_ovm.n_vertices()];
  for (int i = 0; i < _ovm.n_vertices(); ++i) {
    fixed[i] = 0;
  }

  // ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(),
  // _ovm.n_cells(), HEXAHEDRON, connection.data());
  MeshImpl msqmesh((int)_ovm.n_vertices(), (int)_ovm.n_cells(), HEXAHEDRON,
                   fixed, coords.data(), connection.data());
  MsqError err;
  msqmesh.write_vtk(outname.c_str(), err);
  // msqmesh.clear();
  delete[] fixed;
}

void MsqOperator::Ovm2MsqOut(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                             std::vector<OpenVolumeMesh::VertexHandle> &_tagged,
                             std::string outname) {
  std::vector<double> coords;
  std::map<OpenVolumeMesh::VertexHandle, int> v_idx;
  int i = 0;
  for (auto v_it = _ovm.vertices_begin(); v_it != _ovm.vertices_end(); ++v_it) {
    auto c = _ovm.vertex(*v_it);
    v_idx[*v_it] = i++;
    coords.push_back(c[0]);
    coords.push_back(c[1]);
    coords.push_back(c[2]);
  }
  // std::vector<bool> fixedflag;
  std::vector<int> connection;
  std::vector<OpenVolumeMesh::VertexHandle> _xb_vertices(4);
  std::vector<OpenVolumeMesh::VertexHandle> _xf_vertices(4);
  for (auto _ch : _ovm.cells()) {
    auto xf = _ovm.xfront_halfface(_ch);
    auto xb = _ovm.xback_halfface(_ch);
    auto _xb_range = _ovm.halfface_vertices(xb);
    auto _xf_range = _ovm.halfface_vertices(xf);
    //_xb_vertices.clear();
    int _num = 0;
    for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
      _xb_vertices.at(_num++) = *_vh;
    }
    _xf_vertices = opposite_vertex_in_cell1(_ovm, _ch, xb, _xb_vertices);
    connection.push_back(v_idx[_xb_vertices[0]]);
    connection.push_back(v_idx[_xb_vertices[1]]);
    connection.push_back(v_idx[_xb_vertices[2]]);
    connection.push_back(v_idx[_xb_vertices[3]]);
    connection.push_back(v_idx[_xf_vertices[0]]);
    connection.push_back(v_idx[_xf_vertices[1]]);
    connection.push_back(v_idx[_xf_vertices[2]]);
    connection.push_back(v_idx[_xf_vertices[3]]);
  }
  // std::shared_ptr<bool> fixed = std::make_shared(new
  // bool(_ovm.n_vertices()));
  bool *fixed = new bool[_ovm.n_vertices()];
  for (int i = 0; i < _ovm.n_vertices(); ++i) {
    fixed[i] = 0;
  }
  for (auto _vh : _tagged) {
    fixed[v_idx[_vh]] = 1;
  }

  // ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(),
  // _ovm.n_cells(), HEXAHEDRON, connection.data());
  MeshImpl msqmesh((int)_ovm.n_vertices(), (int)_ovm.n_cells(), HEXAHEDRON,
                   fixed, coords.data(), connection.data());
  MsqError err;
  msqmesh.write_vtk(outname.c_str(), err);
  // msqmesh.clear();
  delete[] fixed;
}

void MsqOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm) {

  // transform_ovm_2_matrix_scaf(_ovm);
  MeshImpl msqmesh;
  MsqError err;
  msqmesh.read_vtk("ttmp.vtk", err);
  // auto msqmesh = OVMmesh2MSQmesh(_ovm);

  QualityAssessor qa(false, false);
  InstructionQueue inv_dect;
  inv_dect.add_quality_assessor(&qa, err);
  inv_dect.run_instructions(&msqmesh, err);
  int inv_ele = 0, inv_sam = 0;
  qa.get_inverted_element_count(inv_ele, inv_sam, err);

  {
    IdealShapeTarget target;
    // TShapeSizeB1 m1;
    TShapeSizeB1 m2;
    // TSum mymetric(&m1, &m2);

    TQualityMetric metric_0(&target, &m2);
    ElementPMeanP metric(1.0, &metric_0);
    PMeanPTemplate obj_func_opt(1.0, &metric);
    QuasiNewton improver(&obj_func_opt);

    improver.use_global_patch();
    // improver.set_inner_termination_criterion(&e);
    InstructionQueue queue;
    queue.set_master_quality_improver(&improver, err);
    queue.run_instructions(&msqmesh, err);
  }

  std::vector<Mesquite::Mesh::ElementHandle> vertices;
  std::vector<MsqVertex> coordinates;
  msqmesh.get_all_vertices(vertices, err);
  coordinates.resize(vertices.size());
  msqmesh.vertices_get_coordinates(vertices.data(), coordinates.data(),
                                   coordinates.size(), err);
  int i = 0;
  for (auto v_it = _ovm.vertices_begin(); v_it != _ovm.vertices_end(); ++v_it) {
    _ovm.set_vertex(*v_it, OpenVolumeMesh::Geometry::Vec3d(coordinates[i].x(),
                                                           coordinates[i].y(),
                                                           coordinates[i].z()));
    ++i;
  }
}

void MsqOperator::Optimize_old(
    OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm) {

  std::vector<double> coords;
  std::map<OpenVolumeMesh::VertexHandle, int> v_idx;
  int i = 0;
  for (auto v_it = _ovm.vertices_begin(); v_it != _ovm.vertices_end(); ++v_it) {
    auto c = _ovm.vertex(*v_it);
    v_idx[*v_it] = i++;
    coords.push_back(c[0]);
    coords.push_back(c[1]);
    coords.push_back(c[2]);
  }
  // std::vector<bool> fixedflag;
  std::vector<int> connection;
  std::vector<OpenVolumeMesh::VertexHandle> _xb_vertices(4);
  std::vector<OpenVolumeMesh::VertexHandle> _xf_vertices(4);
  for (auto _ch : _ovm.cells()) {
    auto xf = _ovm.xfront_halfface(_ch);
    auto xb = _ovm.xback_halfface(_ch);
    auto _xb_range = _ovm.halfface_vertices(xb);
    auto _xf_range = _ovm.halfface_vertices(xf);
    //_xb_vertices.clear();
    int _num = 0;
    for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
      _xb_vertices.at(_num++) = *_vh;
    }
    _xf_vertices = opposite_vertex_in_cell1(_ovm, _ch, xb, _xb_vertices);
    connection.push_back(v_idx[_xb_vertices[0]]);
    connection.push_back(v_idx[_xb_vertices[1]]);
    connection.push_back(v_idx[_xb_vertices[2]]);
    connection.push_back(v_idx[_xb_vertices[3]]);
    connection.push_back(v_idx[_xf_vertices[0]]);
    connection.push_back(v_idx[_xf_vertices[1]]);
    connection.push_back(v_idx[_xf_vertices[2]]);
    connection.push_back(v_idx[_xf_vertices[3]]);
  }
  bool *fixed = new bool[_ovm.n_vertices()];
  for (int i = 0; i < _ovm.n_vertices(); ++i) {
    if (i <= 1)
      fixed[i] = 1;
    else
      fixed[i] = 0;
  }
  int _test_cnum = _ovm.n_cells();
  // ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(),
  // _ovm.n_cells(), HEXAHEDRON, connection.data());
  MeshImpl msqmesh((int)_ovm.n_vertices(), (int)_ovm.n_cells(), HEXAHEDRON,
                   fixed, coords.data(), connection.data());
  // auto msqmesh = OVMmesh2MSQmesh(_ovm);

  MsqError err;

  QualityAssessor qa(false, false);
  InstructionQueue inv_dect;
  inv_dect.add_quality_assessor(&qa, err);
  inv_dect.run_instructions(&msqmesh, err);
  int inv_ele = 0, inv_sam = 0;
  qa.get_inverted_element_count(inv_ele, inv_sam, err);

  if (false) {
    IdealShapeTarget target;
    // TShapeSizeB1 m1;
    TShapeSizeNB3 m2;
    // TSum mymetric(&m1, &m2);

    TQualityMetric metric_0(&target, &m2);
    ElementPMeanP metric(1.0, &metric_0);
    PMeanPTemplate obj_func_opt(1.0, &metric);
    QuasiNewton improver(&obj_func_opt);
    ;
    improver.use_global_patch();
    // improver.set_inner_termination_criterion(&e);
    InstructionQueue queue;
    queue.set_master_quality_improver(&improver, err);
    queue.run_instructions(&msqmesh, err);
  } else {
    IdealShapeTarget target;
    // TShapeSizeB1 m1;
    TShapeSizeB1 m2;
    // TSum mymetric(&m1, &m2);

    TQualityMetric metric_0(&target, &m2);
    ElementPMeanP metric(1.0, &metric_0);
    PMeanPTemplate obj_func_opt(1.0, &metric);
    QuasiNewton improver(&obj_func_opt);

    improver.use_global_patch();
    // improver.set_inner_termination_criterion(&e);
    InstructionQueue queue;
    queue.set_master_quality_improver(&improver, err);
    queue.run_instructions(&msqmesh, err);
  }

  std::vector<Mesh::ElementHandle> vertices;
  std::vector<MsqVertex> coordinates;
  msqmesh.get_all_vertices(vertices, err);
  coordinates.resize(vertices.size());
  msqmesh.vertices_get_coordinates(vertices.data(), coordinates.data(),
                                   coordinates.size(), err);
  i = 0;
  for (auto v_it = _ovm.vertices_begin(); v_it != _ovm.vertices_end(); ++v_it) {
    _ovm.set_vertex(*v_it, OpenVolumeMesh::Geometry::Vec3d(coordinates[i].x(),
                                                           coordinates[i].y(),
                                                           coordinates[i].z()));
    ++i;
  }
  delete[] fixed;
}

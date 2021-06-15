
#ifndef MYMESH_H
#define MYMESH_H

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
//#include <OpenVolumeMesh/FileManager/FileManagerT_impl.hh>
#include "arapoperator.h"
#include "utils.h"

struct untangleData {
  OpenVolumeMesh::Geometry::Vec3d norm_d;
  OpenVolumeMesh::Geometry::Vec3d face_d;
  OpenVolumeMesh::Geometry::Vec3d bline_d;
  OpenVolumeMesh::Geometry::Vec3d assistV_d;
  std::vector<OpenVolumeMesh::Geometry::Vec3d> vertices;
};

struct controlPara {
  // debug用
  std::vector<int> vint = {};
  std::vector<double> vdouble = {};
  std::vector<char> vchar = {};
};

class MyMesh {
public:
  //可以用>>操作符
  bool ReadTopoFromFile(const std::string &filename);
  bool WriteGeomToVTKFile(const std::string &filename);
  bool WriteGeomToVTKFile(const std::string &filename,
                          std::vector<OpenVolumeMesh::VertexHandle> &_tagged);
  bool GenerateOrder();
  bool GenerateOneCell(const OpenVolumeMesh::CellHandle &_ch);

  bool Optimize();

  int GetTopoVnum() { return m_topomesh.n_vertices(); }
  int GetTopoCnum() { return m_topomesh.n_cells(); }
  int GetGeomVnum() { return m_mesh.n_vertices(); }
  int GetGeomCnum() { return m_mesh.n_cells(); }

  OpenVolumeMesh::CellHandle GetCurrentCellHandle() {
    return m_generate_order[m_cur_idx++];
  }
  bool checkTopo();

private:
  //加入cell的一系列函数
  bool
  AddOneCellCase0(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  bool
  AddOneCellCase1(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  bool
  AddOneCellCase2(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  bool
  AddOneCellCase3(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  bool
  AddOneCellCase4(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  bool
  AddOneCellCase5(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  bool
  AddOneCellCase6(const OpenVolumeMesh::CellHandle &_ch,
                  const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
                  const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

  //一系列操作函数，形式上简洁一些
  OpenVolumeMesh::CellHandle getGeomC(const OpenVolumeMesh::CellHandle &_ch) {
    return m_tm2m_mapping[_ch];
  }
  OpenVolumeMesh::CellHandle getTopoC(const OpenVolumeMesh::CellHandle &_ch) {
    return m_m2tm_mapping[_ch];
  }
  OpenVolumeMesh::VertexHandle
  getTopoV(const OpenVolumeMesh::VertexHandle &_vh) {
    return m_m2tm_v_mapping[_vh];
  }
  OpenVolumeMesh::VertexHandle
  getGeomV(const OpenVolumeMesh::VertexHandle &_vh) {
    return m_tm2m_v_mapping[_vh];
  }
  OpenVolumeMesh::Geometry::Vec3d
  getCoord_topo(const OpenVolumeMesh::VertexHandle &_vh) {
    auto _geomv = getGeomV(_vh);
    return m_mesh.vertex(_geomv);
  }
  void setCoord_topo(const OpenVolumeMesh::VertexHandle &_vh,
                     const OpenVolumeMesh::Geometry::Vec3d _vec) {
    auto _geomv = getGeomV(_vh);
    m_mesh.set_vertex(_geomv, _vec);
  }

  // std::vector<std::vector<OpenVolumeMesh::Vec3d>> untangleBottomFace(const
  // std::vector<untangleData> &uData);

private:
  //私有变量
  OpenVolumeMesh::GeometricHexahedralMeshV3d m_mesh;
  OpenVolumeMesh::TopologicHexahedralMesh m_topomesh;
  std::vector<OpenVolumeMesh::VertexHandle> m_vertices;
  std::vector<OpenVolumeMesh::VertexHandle> m_topo_vertices;
  std::vector<std::vector<uint32_t>> m_cells;
  std::map<OpenVolumeMesh::CellHandle, OpenVolumeMesh::CellHandle>
      m_m2tm_mapping;
  std::map<OpenVolumeMesh::CellHandle, OpenVolumeMesh::CellHandle>
      m_tm2m_mapping;
  std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::VertexHandle>
      m_m2tm_v_mapping;
  std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::VertexHandle>
      m_tm2m_v_mapping;

  std::vector<OpenVolumeMesh::CellHandle> m_generate_order;
  int m_cur_idx = 0;
};
#endif
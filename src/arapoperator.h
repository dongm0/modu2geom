#pragma once

#include "my_scaf.h"
#include "utils.h"

#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

class ArapOperator {
public:
  static ArapOperator &Instance() {
    static ArapOperator m_instance;
    return m_instance;
  }
  void Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                std::map<OpenVolumeMesh::VertexHandle,
                         OpenVolumeMesh::Geometry::Vec3d> &fixed);
  void Deformation(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                   std::map<OpenVolumeMesh::VertexHandle,
                            OpenVolumeMesh::Geometry::Vec3d> &fixed);

private:
  ArapOperator() {}
  ArapOperator(const ArapOperator &_rhs) = delete;
  ArapOperator operator=(const ArapOperator *_rhs) = delete;
};
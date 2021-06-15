#include "my_scaf.h"
#include "utils.h"
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

class ScafOperator {
public:
  static ScafOperator &Instance() {
    static ScafOperator m_instance;
    return m_instance;
  }
  void Deformation(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                   std::vector<OpenVolumeMesh::VertexHandle> fixed);
  void Optimize(
      OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
      std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
          fixed);

private:
  ScafOperator() {}
  ScafOperator(const ScafOperator &_rhs) {}
  ScafOperator operator=(const ScafOperator *_rhs) {}
};
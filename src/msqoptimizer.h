#include "utils.h"
#include <Mesquite_all_headers.hpp>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

using namespace Mesquite;

class MsqOperator {
public:
  static MsqOperator &Instance() {
    static MsqOperator m_instance;
    return m_instance;
  }
  void Ovm2MsqOut(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                  std::string outname);
  void Ovm2MsqOut(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                  std::vector<OpenVolumeMesh::VertexHandle> &_tagged,
                  std::string outname);
  void Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm);
  void Optimize_old(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm);
  void Msq2Ovm(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
               MeshImpl &mesh);

private:
  MsqOperator() {}
  MsqOperator(const MsqOperator &_rhs) {}
  MsqOperator operator=(const MsqOperator *_rhs) {}

  // ArrayMesh OVMmesh2MSQmesh(OpenVolumeMesh::GeometricHexahedralMeshV3d _ovm);
};

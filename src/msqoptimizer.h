#include "utils.h"
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <Mesquite/Mesquite_all_headers.hpp>

using namespace Mesquite;

class MsqOperator {
public:
    static MsqOperator& Instance() {
        static MsqOperator m_instance;
        return m_instance;
    }
    Mesquite::MeshImpl Ovm2Msq(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm);
    void Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm);
    void Msq2Ovm(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, MeshImpl &mesh);
private:
    MsqOperator(){}
    MsqOperator(const MsqOperator &_rhs){}
    MsqOperator operator=(const MsqOperator *_rhs){}

    //ArrayMesh OVMmesh2MSQmesh(OpenVolumeMesh::GeometricHexahedralMeshV3d _ovm);
};
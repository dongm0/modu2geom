#include "utils.h"
#include <Mesquite/Mesquite_all_headers.hpp>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

using namespace Mesquite;

class MsqOperator {
public:
    static MsqOperator& Instance() {
        static MsqOperator m_instance;
        return m_instance;
    }
    void Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm);
private:
    MsqOperator(){}
    MsqOperator(const MsqOperator &_rhs){}
    MsqOperator operator=(const MsqOperator *_rhs){}

    //ArrayMesh OVMmesh2MSQmesh(OpenVolumeMesh::GeometricHexahedralMeshV3d _ovm);
};
#include "utils.h"
#include <igl/arap.h>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

class ArapOperator {
public:
    static ArapOperator& Instance() {
        static ArapOperator m_instance;
        return m_instance;
    }
    void Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, std::vector<OpenVolumeMesh::VertexHandle> fixed);
private:
    ArapOperator(){}
    ArapOperator(const ArapOperator &_rhs){}
    ArapOperator operator=(const ArapOperator *_rhs){}

};}
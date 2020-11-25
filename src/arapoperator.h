#include "utils.h"
#include <igl/arap.h>

class ArapOperator {
public:
    static ArapOperator& Instance() {
        static ArapOperator m_instance;
        return m_instance;
    }
    void Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, );
private:
    ArapOperator(){}
    ArapOperator(const ArapOperator &_rhs){}
    ArapOperator operator=(const ArapOperator *_rhs){}

};}
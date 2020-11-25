#include "arapoperator.h"
#include <Eigen/Dense>

void ArapOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, std::vector<OpenVolumeMesh::VertexHandle> fixed) {
    //找表面并建立映射
    std::vector<OpenVolumeMesh::HalfFaceHandle> surface_hf;
    std::vector<OpenVolumeMesh::VertexHandle> vertices;
    Eigen::Matrix<double, Eigen::Dynamic, 3> V;
    Eigen::Matrix<int, Eigen::Dynamic, 3> F;
    int surface_vnum = 0;
    {
        for (auto hf_it=_ovm.halffaces_begin(); hf_it!=_ovm.halffaces_end(); ++hf_it) {
            if (_ovm.is_boundary(*hf_it)) {
                surface_hf.push_back(*hf_it);
            }
        }
        std::set<OpenVolumeMesh::VertexHandle> vset;
        for (const auto &x : surface_hf) {
            auto range = _ovm.halfface_vertices(x);
            for (auto v_it=range.first; v_it!=range.second; ++v_it) {
                vset.insert(*v_it);
            }
        }
        surface_vnum = vset.size();
        vertices.assign(vset.begin(), vset.end());

        V.resize(surface_vnum, 3);
        F.resize(surface_hf.size()*2, 3);

        for (int i=0; i<surface_vnum; ++i) {
            
        }
    }

}
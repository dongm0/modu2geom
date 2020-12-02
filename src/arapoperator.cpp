#include "arapoperator.h"
#include <Eigen/Dense>

void ArapOperator::Deformation(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, std::vector<OpenVolumeMesh::VertexHandle> fixed) {
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
        std::map<OpenVolumeMesh::VertexHandle, int> inv_mapping;
        for (int i=0; i<vertices.size(); ++i) {
            inv_mapping[vertices[i]] = i;
        }

        V.resize(surface_vnum, 3);
        F.resize(surface_hf.size()*2, 3);

        for (int i=0; i<surface_vnum; ++i) {
            auto coord = vertices.at(i);
            V(i, 0) = _ovm.vertex(coord)[0];
            V(i, 1) = _ovm.vertex(coord)[1];
            V(i, 2) = _ovm.vertex(coord)[2];
        }

        for (int i=0; i<surface_hf.size(); ++i) {
            auto range = _ovm.halfface_vertices(surface_hf[i]);
            int v1 = inv_mapping[*range.first];
            int v2 = inv_mapping[*(range.first+1)];
            int v3 = inv_mapping[*(range.first+2)];
            int v4 = inv_mapping[*(range.first+3)];
            F(i*2, 0) = v1;
            F(i*2, 1) = v2;
            F(i*2, 2) = v3;
            F(i*2+1, 0) = v3;
            F(i*2+1, 1) = v4;
            F(i*2+1, 2) = v1;
        }
    }

}
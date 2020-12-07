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
        std::unordered_map<OpenVolumeMesh::VertexHandle, int> inv_mapping;
        for (int i=0; i<vertices.size(); ++i) {
            inv_mapping[vertices[i]] = i;
        }

        V.resize(surface_vnum, 3);
        F.resize(surface_hf.size()*2, 3);

        for (int i=0; i<surface_vnum; ++i) {
            auto coord = vertices.at(i);
            V(i, 0) = coord[0];
            V(i, 1) = coord[1];
            V(i, 2) = coord[2];
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


void ArapOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d fixed) {
    using namespace Eigen;
    using namespace OpenVolumeMesh;
    Matrix<double, Dynamic, Dynamic> V;
    Matrix<int, Dynamic, Dynamic> T;
    Matrix<int, 1, 1> b;
    Matrix<double, 1, 1> bc;
    std::map<VertexHandle, int> mapping;
    int i=0;
    int j=0;
    std::vector<int> boundary;
    for (auto _vh = _ovm.vertices_begin(); _vh != _ovm.vertices.end(); ++_vh) {
        auto point = _ovm.vertex(_vh);
        V(i, 0) = point[0];
        V(i, 1) = point[1];
        V(i, 2) = point[2];
        if (fixed.find()) {
            b(j) = i;
            bc(j, 0) = fixed[_vh][0];
            bc(j, 1) = fixed[_vh][1];
            bc(j, 2) = fixed[_vh][2];
        }
        mapping[_vh] = i;
        ++i;
    }
    auto V0 = V;
    for (int k=0; k<b.size(); ++k) {
        V0(b(k), 0) = bc(k, 0);
        V0(b(k), 1) = bc(k, 1);
        V0(b(k), 2) = bc(k, 2);
    }
    i=0;
    for (auto _ch = _ovm.cells_begin(); _ch != _ovm.cells_end(); ++_ch) {
        for (auto _vh = _ovm.cv_iter(); _vh)
    }
}
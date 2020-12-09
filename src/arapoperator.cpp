#include "arapoperator.h"
#include <Eigen/Dense>
/*
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
*/

const int hex2tet[32] = {0, 1, 3, 4, 
                         1, 2, 0, 5, 
                         2, 3, 6, 1, 
                         3, 0, 2, 7, 
                         4, 7, 5, 0, 
                         6, 5, 7, 2, 
                         5, 4, 6, 1};

void ArapOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d> fixed) {
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
    for (auto _vh = _ovm.vertices_begin(); _vh != _ovm.vertices_end(); ++_vh) {
        auto point = _ovm.vertex(*_vh);
        V(i, 0) = point[0];
        V(i, 1) = point[1];
        V(i, 2) = point[2];
        if (fixed.find(*_vh) != fixed.end()) {
            b(j) = i;
            bc(j, 0) = fixed[*_vh][0];
            bc(j, 1) = fixed[*_vh][1];
            bc(j, 2) = fixed[*_vh][2];
        }
        mapping[*_vh] = i;
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
        auto _xb = _ovm.xback_halfface(*_ch);
        auto _xf = _ovm.xfront_halfface(*_ch);
        auto _xb_range = _ovm.halfface_vertices(_xb);
        auto _xf_range = _ovm.halfface_vertices(_xf);
        std::vector<int> _cell_number;
        std::vector<int> _xf_number;
        for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
            _cell_number.push_back(mapping[*_vh]);
        }
        int _cnt = 0;
        for (auto _vh = _xf_range.first; _vh != _xf_range.second; ++_vh) {
            _xf_number.push_back(mapping[*_vh]);
        }
        for (auto _vh = _ovm.vv_iter(*_xb_range.first); _vh.valid(); ++_vh) {
            for (int ii=0; ii<4; ++ii) {
                if (mapping[*_vh] == _xf_number[ii]) {
                    for (int jj=0; jj<4; ++jj) {
                        _cell_number.push_back(_xf_number[(ii-jj+3)%4]);
                    }
                    break;
                }
            }
            if (_cell_number.size() > 4) break;
        }
        for (int j=0; j<8; ++j) {
            for (int k=0; k<4; ++k) {
                T(i*8+j, k) = _cell_number[hex2tet[j*4+k]];
            }
        }
    }

    //slim完了
    std::vector<VertexHandle> inmapping(mapping.size());
    for (auto x : mapping) {
        inmapping[x.second] = x.first;
    }
    for (int i=0; i<mapping.size(); ++i) {
        _ovm.set_vertex(inmapping[i], Vec3d(V(i, 0), V(i, 1), V(i, 2)));
    }

    
}
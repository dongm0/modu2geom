#include "arapoperator.h"
#include <Eigen/Dense>
#include <igl/slim.h>

const int hex2tet[32] = {0, 1, 3, 4, 
                         1, 2, 0, 5, 
                         2, 3, 1, 6, 
                         3, 0, 2, 7, 
                         4, 7, 5, 0, 
                         7, 6, 4, 3,
                         6, 5, 7, 2, 
                         5, 4, 6, 1};

void ArapOperator::Deformation(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm, std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d> fixed) {
    using namespace Eigen;
    using namespace OpenVolumeMesh;
    Matrix<double, Dynamic, Dynamic> V;
    V.resize(_ovm.n_vertices(), 3);
    Matrix<int, Dynamic, Dynamic> T;
    T.resize(_ovm.n_cells()*8, 4);
    //Matrix<int, 1, 1> b;
    VectorXi b;
    b.resize(fixed.size());
    Matrix<double, Dynamic, Dynamic> bc;
    bc.resize(fixed.size(), 3);
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
            j++;
        }
        mapping[*_vh] = i;
        ++i;
    }
    auto V0 = V;

    std::vector<std::vector<double>> V_tmp;
    for (int vn=0; vn<V.rows(); ++vn) {
        V_tmp.push_back({V(vn, 0), V(vn, 1), V(vn, 2)});
    }
    std::vector<VertexHandle> inmapping(mapping.size());
    for (auto x : mapping) {
        inmapping[x.second] = x.first;
    }
    int cnumber=0;
    for (auto _ch = _ovm.cells_begin(); _ch != _ovm.cells_end(); ++_ch) {
        auto _xb = _ovm.xback_halfface(*_ch);
        auto _xf = _ovm.xfront_halfface(*_ch);
        auto _xb_range = _ovm.halfface_vertices(_xb);
        auto _xf_range = _ovm.halfface_vertices(_xf);
        std::vector<int> _cell_number;
        std::vector<int> _xf_number;
        std::set<int> _xf_set;
        for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
            _cell_number.push_back(mapping[*_vh]);
        }
        int _cnt = 0;
        for (auto _vh = _xf_range.first; _vh != _xf_range.second; ++_vh) {
            _xf_number.push_back(mapping[*_vh]);
            _xf_set.insert(mapping[*_vh]);
        }
        for (j=0; j<4; ++j) {
            auto _vh1 = inmapping[_cell_number[j]];
            for (auto vv : _ovm.vertex_vertices(_vh1)) {
                if (_xf_set.count(mapping[vv]) != 0) {
                    _cell_number.push_back(mapping[vv]);
                    break;
                }
            }
        }

        for (int j=0; j<8; ++j) {
            for (int k=0; k<4; ++k) {
                T(cnumber*8+j, k) = _cell_number[hex2tet[j*4+k]];
            }
        }
        cnumber++;
    }
    std::vector<std::vector<int>> T_tmp;
    for (int cn=0; cn<T.rows(); ++cn) {
        T_tmp.push_back({T(cn, 0), T(cn, 1), T(cn, 2), T(cn, 3)});
    }
    //slim
    igl::SLIMData sData;
    igl::slim_precompute(V, T, V0, sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 1e6);
    igl::slim_solve(sData, 5);
    std::vector<std::vector<double>> V1_tmp;
    for (int vn=0; vn<V.rows(); ++vn) {
        V1_tmp.push_back({sData.V_o(vn, 0), sData.V_o(vn, 1), sData.V_o(vn, 2)});
    }

    //slim完了
    
    for (int i=0; i<mapping.size(); ++i) {
        _ovm.set_vertex(inmapping[i], Vec3d(sData.V_o(i, 0), sData.V_o(i, 1), sData.V_o(i, 2)));
    }

    
}
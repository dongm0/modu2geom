#include "msqoptimizer.h"

using namespace Mesquite;
/*
ArrayMesh MsqOperator::OVMmesh2MSQmesh(OpenVolumeMesh::GeometricHexahedralMeshV3d _ovm) {
    std::vector<double> coords;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        auto c = _ovm.vertex(*v_it);
        coords.push_back(c[0]);
        coords.push_back(c[1]);
        coords.push_back(c[2]);
    }
    std::vector<int> fixedflag;
    std::vector<int> connection;
    //可能销毁，想想办法
    ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(), _ovm.n_cells(), HEXAHEDRON, connection.data());
}
*/

MeshImpl MsqOperator::Ovm2Msq(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm) {
        std::vector<double> coords;
    std::map<OpenVolumeMesh::VertexHandle, int> v_idx;
    int i=0;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        auto c = _ovm.vertex(*v_it);
        v_idx[*v_it] = i++;
        coords.push_back(c[0]);
        coords.push_back(c[1]);
        coords.push_back(c[2]);
    }
    //std::vector<bool> fixedflag;
    std::vector<int> connection;
    for (auto _ch : _ovm.cells()) {
        auto xf = _ovm.xfront_halfface(_ch);
        auto xb = _ovm.xback_halfface(_ch);
        auto _xb_range = _ovm.halfface_vertices(xb);
        auto _xf_range = _ovm.halfface_vertices(xf);
        std::vector<int> _cell_number;
        std::vector<int> _xf_number;
        for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
            _cell_number.push_back(v_idx[*_vh]);
        }
        int _cnt = 0;
        for (auto _vh = _xf_range.first; _vh != _xf_range.second; ++_vh) {
            _xf_number.push_back(v_idx[*_vh]);
        }
        for (auto _vh = _ovm.vv_iter(*_xb_range.first); _vh.valid(); ++_vh) {
            for (int ii=0; ii<4; ++ii) {
                if (v_idx[*_vh] == _xf_number[ii]) {
                    for (int jj=0; jj<4; ++jj) {
                        _cell_number.push_back(_xf_number[(ii-jj+3)%4]);
                    }
                    break;
                }
            }
            if (_cell_number.size() > 4) break;
        };
        for (int ii=0; ii<8; ++ii) {
            connection.push_back(_cell_number[ii]);
        }
    }
    bool *fixed = new bool(_ovm.n_vertices());
    for (int i=0; i<_ovm.n_vertices(); ++i) {
        if (i<=1) fixed[i] = 1;
        else fixed[i] = 0;
    }
    
    //ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(), _ovm.n_cells(), HEXAHEDRON, connection.data());
    MeshImpl msqmesh((int)_ovm.n_vertices(), (int)_ovm.n_cells(), HEXAHEDRON, fixed, coords.data(), connection.data());
    delete[] fixed;
    return msqmesh;
}

void MsqOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm) {
    
    std::vector<double> coords;
    std::map<OpenVolumeMesh::VertexHandle, int> v_idx;
    int i=0;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        auto c = _ovm.vertex(*v_it);
        v_idx[*v_it] = i++;
        coords.push_back(c[0]);
        coords.push_back(c[1]);
        coords.push_back(c[2]);
    }
    //std::vector<bool> fixedflag;
    std::vector<int> connection;
    for (auto _ch : _ovm.cells()) {
        auto xf = _ovm.xfront_halfface(_ch);
        auto xb = _ovm.xback_halfface(_ch);
        auto _xb_range = _ovm.halfface_vertices(xb);
        auto _xf_range = _ovm.halfface_vertices(xf);
        std::vector<int> _cell_number;
        std::vector<int> _xf_number;
        for (auto _vh = _xb_range.first; _vh != _xb_range.second; ++_vh) {
            _cell_number.push_back(v_idx[*_vh]);
        }
        int _cnt = 0;
        for (auto _vh = _xf_range.first; _vh != _xf_range.second; ++_vh) {
            _xf_number.push_back(v_idx[*_vh]);
        }
        for (auto _vh = _ovm.vv_iter(*_xb_range.first); _vh.valid(); ++_vh) {
            for (int ii=0; ii<4; ++ii) {
                if (v_idx[*_vh] == _xf_number[ii]) {
                    for (int jj=0; jj<4; ++jj) {
                        _cell_number.push_back(_xf_number[(ii-jj+3)%4]);
                    }
                    break;
                }
            }
            if (_cell_number.size() > 4) break;
        };
        for (int ii=0; ii<8; ++ii) {
            connection.push_back(_cell_number[ii]);
        }
    }
    bool *fixed = new bool(_ovm.n_vertices());
    for (int i=0; i<_ovm.n_vertices(); ++i) {
        if (i<=1) fixed[i] = 1;
        else fixed[i] = 0;
    }
    
    //ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(), _ovm.n_cells(), HEXAHEDRON, connection.data());
    MeshImpl msqmesh((int)_ovm.n_vertices(), (int)_ovm.n_cells(), HEXAHEDRON, fixed, coords.data(), connection.data());
    //auto msqmesh = OVMmesh2MSQmesh(_ovm);


    MsqError err;
    
    IdealShapeTarget target;
    //TShapeSizeB1 m1;
    TShapeSizeB3 m2;
    //TSum mymetric(&m1, &m2);
    TQualityMetric metric_0(&target, &m2);
    ElementPMeanP metric(1.0, &metric_0);
    PMeanPTemplate obj_func_opt(1.0, &metric);
    ConjugateGradient improver(&obj_func_opt);;
    improver.use_global_patch();
    //improver.set_inner_termination_criterion(&e);
    InstructionQueue queue;
    queue.set_master_quality_improver(&improver, err);
    queue.run_instructions(&msqmesh, err);

    i = 0;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        _ovm.set_vertex(*v_it, OpenVolumeMesh::Geometry::Vec3d(coords[i*3], coords[i*3+1], coords[i*3+2]));
        ++i;
    }
    delete[] fixed;
}

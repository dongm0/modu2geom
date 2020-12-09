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

void MsqOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm) {
    std::vector<double> coords;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        auto c = _ovm.vertex(*v_it);
        coords.push_back(c[0]);
        coords.push_back(c[1]);
        coords.push_back(c[2]);
    }
    std::vector<int> fixedflag;
    std::vector<unsigned long> connection;
    
    ArrayMesh msqmesh(3, _ovm.n_vertices(), coords.data(), fixedflag.data(), _ovm.n_cells(), HEXAHEDRON, connection.data());
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

    int i = 0;
    for (auto v_it = _ovm.vertices_begin(); v_it!=_ovm.vertices_end(); ++v_it) {
        _ovm.set_vertex(*v_it, OpenVolumeMesh::Geometry::Vec3d(coords[i*3], coords[i*3+1], coords[i*3+2]));
    }
}
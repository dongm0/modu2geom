#include <Mesquite/Mesquite_all_headers.hpp>

using namespace Mesquite;

int main(int argc, char** argv) {
    MeshImpl mesh;
    MsqError err;
    mesh.read_vtk(argv[1], err);


    //ShapeImprover smoother;
    //Squared extra_metric;
    //JSquared extra;
    //JSquared extra_metric;
    //IdealWeightMeanRatio extra_metric;
    //smoother.quality_assessor().add_quality_assessment(&extra_metric);
    //MeshDomainAssoc mesh_and_domain = MeshDomainAssoc(&mesh, &plane);
    //smoother.run_instructions( &mesh, err );


    QualityAssessor qa(false, false);
    InstructionQueue q;
    q.add_quality_assessor(&qa, err);
    int inv_ele = 0, inv_sam = 0;
    q.run_instructions(&mesh, err);
    qa.get_inverted_element_count(inv_ele, inv_sam, err);
    std::cout << "inverted elements: " << inv_ele << std::endl;
    std::cout << "inverted samples: " << inv_sam << std::endl;

    //ShapeImprover imp;
    //imp.run_instructions(mesh, err);


    IdealShapeTarget target;
    //TShapeSizeB1 m1;
    TShapeSizeB3 m2;
    //TShapeSizeOrientNB1 m2;
    //TShapeSizeOrientB1 m2;
    //TSum mymetric(&m1, &m2);
    TQualityMetric metric_0(&target, &m2);
    ElementPMeanP metric(1.0, &metric_0);
    PMeanPTemplate obj_func_opt(1.0, &metric);
    
    qa.add_quality_assessment(&metric);
    QuasiNewton improver(&obj_func_opt);;
    improver.use_global_patch();
    //improver.set_inner_termination_criterion(&e);
    InstructionQueue queue;
    //queue.set_master_quality_improver(&improver, err);
    queue.add_quality_assessor(&qa, err);

    queue.run_instructions(&mesh, err);
    auto res = qa.get_all_results();
    for (auto x : res) {
        std::cout << "Max:" << x->get_maximum() << std::endl;
        std::cout << "Ave:" << x->get_average() << std::endl;
        std::cout << "Min:" << x->get_minimum() << std::endl;
    }


    mesh.write_vtk("out.vtk", err);


    return 0;
}
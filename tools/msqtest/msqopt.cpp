#include <Mesquite/Mesquite_all_headers.hpp>

using namespace Mesquite;

int main(int argc, char **argv) {
  MeshImpl mesh;
  MsqError err;
  mesh.read_vtk(argv[1], err);

  IdealShapeTarget target;
  TShapeSizeB3 m2;
  TQualityMetric metric_0(&target, &m2);
  ElementPMeanP metric(1.0, &metric_0);
  PMeanPTemplate obj_func_opt(1.0, &metric);

  // qa.add_quality_assessment(&metric);
  QuasiNewton improver(&obj_func_opt);
  ;
  improver.use_global_patch();
  InstructionQueue queue;
  queue.set_master_quality_improver(&improver, err);
  queue.run_instructions(&mesh, err);

  // ShapeImprovementWrapper improver;

  mesh.write_vtk("out.vtk", err);

  return 0;
}
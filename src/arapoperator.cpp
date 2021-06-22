#include "arapoperator.h"
#include "ovmwrap.h"
#include <Eigen/Dense>
#include <igl/slim.h>

void ArapOperator::Deformation(
    OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        &fixed) {
  using namespace Eigen;
  using namespace OpenVolumeMesh;
  MatrixXd V, bc;
  MatrixXi T, surface;
  VectorXi b;
  transform_hex_to_matrix(V, T, b, bc, surface, _ovm, fixed);
  // printEigenMatrix(bc, std::cout);
  MatrixXd V0 = V;
  igl::SLIMData sData;
  igl::slim_precompute(V, T, V0, sData, igl::SYMMETRIC_DIRICHLET, b, bc, 1e6);
  igl::slim_solve(sData, 5);
  transform_matrix_to_hex(sData.V_o, _ovm);
  // slim
  /*
  igl::my_scaf::SCAFData sData;
  igl::my_scaf::scaf_precompute(V, T, V0, surface, sData,
                                igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b,
                                bc, 1e9);
  V0 = igl::my_scaf::scaf_solve(sData, 5);
  */
  // slim完了
}

void ArapOperator::Optimize(OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
                            std::map<OpenVolumeMesh::VertexHandle,
                                     OpenVolumeMesh::Geometry::Vec3d> &fixed) {
  using namespace Eigen;
  using namespace OpenVolumeMesh;
  MatrixXd V, bc;
  MatrixXi T, surface;
  VectorXi b;
  transform_hex_to_matrix(V, T, b, bc, surface, _ovm, fixed);
  igl::my_scaf::SLIMData data;
  MatrixXd refV;
  igl::my_scaf::slim_precompute(refV, T, V, surface, data,
                                igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b,
                                bc, 1e6);
  igl::my_scaf::slim_solve(data, 3);
  transform_matrix_to_hex(data.V_o, _ovm);
}
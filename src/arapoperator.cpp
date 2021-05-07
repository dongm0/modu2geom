#include "arapoperator.h"
#include "ovmwrap.h"
#include <Eigen/Dense>
#include <igl/slim.h>

void ArapOperator::Deformation(
    OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        fixed) {
  using namespace Eigen;
  using namespace OpenVolumeMesh;
  MatrixXd V, bc;
  MatrixXi T;
  VectorXi b;
  transform_hex_to_matrix(V, T, b, bc, _ovm, fixed);
  MatrixXd V0 = V;
  // slim
  igl::SLIMData sData;
  igl::slim_precompute(V, T, V0, sData,
                       igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 1e6);
  igl::slim_solve(sData, 5);
  // slim完了
  transform_matrix_to_hex(sData.V_o, _ovm);
}

void ArapOperator::Optimize(
    OpenVolumeMesh::GeometricHexahedralMeshV3d &_ovm,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        fixed) {
  using namespace Eigen;
  using namespace OpenVolumeMesh;
  MatrixXd V, bc;
  MatrixXi T, surface;
  VectorXi b;
  transform_hex_to_matrix(V, T, b, bc, surface, _ovm, fixed);
  MatrixXd V_o = V;
  igl::my_scaf::SCAFData data;
  igl::my_scaf::scaf_precompute(V, T, V_o, surface, data,
                                igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b,
                                bc, 1e6);
  V_o = igl::my_scaf::scaf_solve(data, 3);
  transform_matrix_to_hex(V_o, _ovm);
}
#include "arapoperator.h"
#include "myslim.h"
#include "ovmwrapper.h"
#include <Eigen/Dense>
//#include <igl/slim.h>

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
  Eigen::MatrixXi fix_mat(bc.rows(), 3);
  fix_mat.setZero();
  MatrixXd V0 = V;
  SLIMData sData;
  myslim_precompute(sData, std::move(V), std::move(T), false, std::move(V0),
                    Eigen::MatrixXd(), igl::MappingEnergyType::EXP_CONFORMAL, b,
                    bc, fix_mat, 1e6);
  // myslim_precompute(V, T, V0, sData, igl::SYMMETRIC_DIRICHLET, b, bc, 1e6);
  myslim_solve(sData, 5);
  transform_matrix_to_hex(sData.V, _ovm);
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
  // printEigenMatrix(bc, std::cout);
  Eigen::MatrixXd std_ele;
  std_ele.resize(4, 3);
  std_ele << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1;
  Eigen::MatrixXi fix_mat(bc.rows(), 3);
  fix_mat.setZero();

  // MatrixXd V0 = V;
  SLIMData sData;
  myslim_precompute(sData, std::move(V), std::move(T), true, Eigen::MatrixXd(),
                    std_ele, igl::MappingEnergyType::EXP_CONFORMAL, b, bc,
                    fix_mat, 1e6);
  // myslim_precompute(V, T, V0, sData, igl::SYMMETRIC_DIRICHLET, b, bc, 1e6);
  myslim_solve(sData, 5);
  transform_matrix_to_hex(sData.V, _ovm);
}
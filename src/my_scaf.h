#pragma once

//#include "OVMVtkHexIO.h"
//#include "ovmwrap.h"
#include "utils.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/MappingEnergyType.h>
#include <igl/igl_inline.h>
//#include <igl/slim.h>

namespace igl {
namespace my_scaf {
struct SCAFData {
  // parameters
  double scaffold_factor = 10;
  igl::MappingEnergyType scaf_energy =
      igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
  igl::MappingEnergyType slim_energy =
      igl::MappingEnergyType::SYMMETRIC_DIRICHLET;

  int dim = 3;         // 3D only
  double total_energy; // scaffold + isometric
  double energy;       // objective value

  long mv_num = 0, mf_num = 0;
  long sv_num = 0, sf_num = 0;

  double mesh_measure; // area or volume
  double proximal_p = 0;

  std::map<int, Eigen::RowVectorXd> soft_cons;
  double soft_const_p = 1e4;

  // inner mesh
  Eigen::MatrixXd m_V;
  Eigen::MatrixXi m_T;
  Eigen::MatrixXd m_Vref;
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  Eigen::VectorXd m_M;
  Eigen::MatrixXi m_surface;

  uint32_t m_surface_fn, m_surface_vn;
  std::vector<uint32_t> mapping_tetgen2scaf;
  std::unordered_map<uint32_t, uint32_t> mapping_scaf2tetgen;
  std::vector<uint32_t> mapping_t2s;
  std::unordered_map<uint32_t, uint32_t> mapping_s2t;

  // scaffold
  Eigen::MatrixXi s_T;
  Eigen::VectorXd s_M;

  Eigen::MatrixXd w_V;

  // pre_calc data
  Eigen::MatrixXd m_GradRef;

  //
  Eigen::MatrixXd s_Grad;

  bool has_pre_calc = false;
};
// Compute necessary information to start using SCAF
// Inputs:
//		V           #V by 3 list of mesh vertex positions
//		F           #F by 3/3 list of mesh faces (triangles/tets)
//    data          igl::SCAFData
//    slim_energy Energy type to minimize
//    b           list of boundary indices into V (soft constraint)
//    bc          #b by dim list of boundary conditions (soft constraint)
//    soft_p      Soft penalty factor (can be zero)
IGL_INLINE void
scaf_precompute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::MatrixXd &V_init, const Eigen::MatrixXi &surface,
                igl::my_scaf::SCAFData &data,
                igl::MappingEnergyType slim_energy, Eigen::VectorXi &b,
                Eigen::MatrixXd &bc, double soft_p);

// Run iter_num iterations of SCAF, with precomputed data
// Outputs:
//    V_o (in SLIMData): #V by dim list of mesh vertex positions
IGL_INLINE Eigen::MatrixXd scaf_solve(SCAFData &data, int iter_num);

// Set up the SCAF system L * uv = rhs, without solving it.
// Inputs:
//    s:   igl::SCAFData. Will be modified by energy and Jacobian computation.
// Outputs:
//    L:   m by m matrix
//    rhs: m by 1 vector
//         with m = dim * (#V_mesh + #V_scaf - #V_frame)
IGL_INLINE void scaf_system(SCAFData &s, Eigen::SparseMatrix<double> &L,
                            Eigen::VectorXd &rhs);

// Compute SCAF energy
// Inputs:
//    s:     igl::SCAFData
//    w_uv:  (#V_mesh + #V_scaf) by dim matrix
//    whole: Include scaffold if true
IGL_INLINE double compute_energy(SCAFData &s, const Eigen::MatrixXd &w_V,
                                 bool whole);

} // namespace my_scaf
} // namespace igl

#include "my_scaf.cpp"
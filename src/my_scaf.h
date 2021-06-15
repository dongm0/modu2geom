// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MY_SCAF_H
#define MY_SCAF_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/MappingEnergyType.h>
#include <igl/igl_inline.h>

// This option makes the iterations faster (all except the first) by caching the
// sparsity pattern of the matrix involved in the assembly. It should be on if
// you plan to do many iterations, off if you have to change the matrix
// structure at every iteration.
#define SLIM_CACHED

#ifdef SLIM_CACHED
#include <igl/AtA_cached.h>
#endif

namespace igl {
namespace my_scaf { // Compute a SLIM map as derived in "Scalable Locally
                    // Injective Maps" [Rabinovich et al. 2016].
struct SLIMData {
  // Input
  Eigen::MatrixXd V; // #V by 3 list of mesh vertex positions
  Eigen::MatrixXi F; // #F by 3/3 list of mesh faces (triangles/tets)
  MappingEnergyType slim_energy;

  // Optional Input
  // soft constraints
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  double soft_const_p;

  double exp_factor;        // used for exponential energies, ignored otherwise
  bool mesh_improvement_3d; // only supported for 3d

  // Output
  Eigen::MatrixXd V_o; // #V by dim list of mesh vertex positions (dim = 2 for
                       // parametrization, 3 otherwise)
  double energy;       // objective value

  // INTERNAL
  Eigen::VectorXd M;
  double mesh_area;
  double avg_edge_length;
  int v_num;
  int f_num;
  double proximal_p;

  Eigen::VectorXd WGL_M;
  Eigen::VectorXd rhs;
  Eigen::MatrixXd Ri, Ji;
  Eigen::MatrixXd W;
  Eigen::SparseMatrix<double> Dx, Dy, Dz;
  int f_n, v_n;
  bool first_solve;
  bool has_pre_calc = false;
  int dim;

#ifdef SLIM_CACHED
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXi A_data;
  Eigen::SparseMatrix<double> AtA;
  igl::AtA_cached_data AtA_data;
#endif

  // my_scaf
  Eigen::MatrixXd w_V;
  Eigen::MatrixXi w_T;
  Eigen::MatrixXi m_surface;

  uint32_t m_surface_fn, m_surface_vn;
  std::vector<uint32_t> mapping_tetgen2scaf;
  std::unordered_map<uint32_t, uint32_t> mapping_scaf2tetgen;
};

// Compute necessary information to start using SLIM
// Inputs:
//		V           #V by 3 list of mesh vertex positions
//		F           #F by 3/3 list of mesh faces (triangles/tets)
//    b           list of boundary indices into V
//    bc          #b by dim list of boundary conditions
//    soft_p      Soft penalty factor (can be zero)
//    slim_energy Energy to minimize
IGL_INLINE void slim_precompute(const Eigen::MatrixXd &V,
                                const Eigen::MatrixXi &F,
                                const Eigen::MatrixXd &V_init,
                                const Eigen::MatrixXi &surface, SLIMData &data,
                                MappingEnergyType slim_energy,
                                const Eigen::VectorXi &b,
                                const Eigen::MatrixXd &bc, double soft_p);

// Run iter_num iterations of SLIM
// Outputs:
//    V_o (in SLIMData): #V by dim list of mesh vertex positions
IGL_INLINE Eigen::MatrixXd slim_solve(SLIMData &data, int iter_num);

// Internal Routine. Exposed for Integration with SCAF
IGL_INLINE void slim_update_weights_and_closest_rotations_with_jacobians(
    const Eigen::MatrixXd &Ji, igl::MappingEnergyType slim_energy,
    double exp_factor, Eigen::MatrixXd &W, Eigen::MatrixXd &Ri);

IGL_INLINE void slim_buildA(const Eigen::SparseMatrix<double> &Dx,
                            const Eigen::SparseMatrix<double> &Dy,
                            const Eigen::SparseMatrix<double> &Dz,
                            const Eigen::MatrixXd &W,
                            std::vector<Eigen::Triplet<double>> &IJV);
} // namespace my_scaf
} // namespace igl

//#ifndef IGL_STATIC_LIBRARY
#include "my_scaf.cpp"
//#endif

#endif // SLIM_H

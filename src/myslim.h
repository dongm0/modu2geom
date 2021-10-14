#ifndef MYSLIM_H
#define MYSLIM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/MappingEnergyType.h>

struct SLIMData {
  // major
  Eigen::MatrixXd V; // #V by 3 list of mesh vertex positions
  Eigen::MatrixXi T; // #T by 3/3 list of mesh faces (triangles/tets)
  // Eigen::VectorXi weight; // weight of every element;
  igl::MappingEnergyType slim_energy;

  // Optional Input
  bool use_ideal;
  Eigen::MatrixXd V_ref;         // standard mesh reference
  Eigen::MatrixXd ideal_element; // standard element reference

  // soft constraint
  double soft_constraint_factor;
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  Eigen::MatrixXi fixed;

  // internal variables
  Eigen::MatrixXd line_coeff;
  Eigen::MatrixXd plain_coeff;
  // Eigen::MatrixXd sphere_coeff;
  Eigen::MatrixXd V_init;
  double exp_factor; // used for exponential energies, ignored otherwise
  double energy;     // objective value

  Eigen::VectorXd M;
  double mesh_area;

  int v_num;
  int t_num;
  double proximal_p = 0.0001;

  std::vector<Eigen::Triplet<double>> IJV;
  Eigen::VectorXd rhs;
  Eigen::MatrixXd P_ref;
  Eigen::MatrixXd J;

  bool has_pre_calc = false;
  int dim;
};

void myslim_precompute(SLIMData &data, Eigen::MatrixXd &&V, Eigen::MatrixXi &&T,
                       bool use_ideal, Eigen::MatrixXd &&V_ref,
                       Eigen::MatrixXd ideal_element,
                       igl::MappingEnergyType energy_type,
                       const Eigen::VectorXi &b, const Eigen::MatrixXd &bc,
                       const Eigen::MatrixXi &fixed, double soft_p);

void myslim_solve(SLIMData &data, int iter_num);

void myslim_newton_solve(SLIMData &data, int iter_num);
/*
// Internal Routine. Exposed for Integration with SCAF
IGL_INLINE void slim_update_weights_and_closest_rotations_with_jacobians(
    const Eigen::MatrixXd &Ji, igl::MappingEnergyType slim_energy,
    double exp_factor, Eigen::MatrixXd &W, Eigen::MatrixXd &Ri);

IGL_INLINE void slim_buildA(const Eigen::SparseMatrix<double> &Dx,
                            const Eigen::SparseMatrix<double> &Dy,
                            const Eigen::SparseMatrix<double> &Dz,
                            const Eigen::MatrixXd &W,
                            std::vector<Eigen::Triplet<double>> &IJV);
*/

//#include "myslim.cpp"

#endif // MYSLIM_H
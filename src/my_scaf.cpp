
// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Michael Rabinovich
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef MY_SCAF_IMPL
#define MY_SCAF_IMPL

#include "my_scaf.h"
#include "tetgen.h"
#include "utils.h"

#include <igl/arap.h>
#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/edge_lengths.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/grad.h>
#include <igl/local_basis.h>
#include <igl/mapping_energy_with_jacobians.h>
#include <igl/per_face_normals.h>
#include <igl/polar_svd.h>
#include <igl/repdiag.h>
#include <igl/slice_into.h>
#include <igl/vector_area_matrix.h>
#include <igl/volume.h>

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

#include <igl/AtA_cached.h>
#include <igl/Timer.h>
#include <igl/sparse_cached.h>

#ifdef CHOLMOD
#include <Eigen/CholmodSupport>
#endif

namespace igl {
namespace my_scaf {
// static
IGL_INLINE void build_linear_system(SLIMData &s,
                                    Eigen::SparseMatrix<double> &L);
IGL_INLINE void buildRhs(SLIMData &s, const Eigen::SparseMatrix<double> &A);
IGL_INLINE void add_soft_constraints(SLIMData &s,
                                     Eigen::SparseMatrix<double> &L);
IGL_INLINE double compute_soft_const_energy(SLIMData &s,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F,
                                            const Eigen::MatrixXd &V_o);

IGL_INLINE void compute_jacobians(igl::my_scaf::SLIMData &s,
                                  const Eigen::MatrixXd &uv) {
  if (s.F.cols() == 3) {
    // Ji=[D1*u,D2*u,D1*v,D2*v];
    s.Ji.col(0) = s.Dx * uv.col(0);
    s.Ji.col(1) = s.Dy * uv.col(0);
    s.Ji.col(2) = s.Dx * uv.col(1);
    s.Ji.col(3) = s.Dy * uv.col(1);
  } else /*tet mesh*/ {
    // Ji=[D1*u,D2*u,D3*u, D1*v,D2*v, D3*v, D1*w,D2*w,D3*w];
    s.Ji.col(0) = s.Dx * uv.col(0);
    s.Ji.col(1) = s.Dy * uv.col(0);
    s.Ji.col(2) = s.Dz * uv.col(0);
    s.Ji.col(3) = s.Dx * uv.col(1);
    s.Ji.col(4) = s.Dy * uv.col(1);
    s.Ji.col(5) = s.Dz * uv.col(1);
    s.Ji.col(6) = s.Dx * uv.col(2);
    s.Ji.col(7) = s.Dy * uv.col(2);
    s.Ji.col(8) = s.Dz * uv.col(2);
  }
}

IGL_INLINE void update_weights_and_closest_rotations(SLIMData &s,
                                                     Eigen::MatrixXd &uv) {
  compute_jacobians(s, uv);
  igl::my_scaf::slim_update_weights_and_closest_rotations_with_jacobians(
      s.Ji, s.slim_energy, s.exp_factor, s.W, s.Ri);
}

IGL_INLINE void solve_weighted_arap(SLIMData &s, const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F,
                                    Eigen::MatrixXd &uv,
                                    Eigen::VectorXi &soft_b_p,
                                    Eigen::MatrixXd &soft_bc_p) {
  using namespace Eigen;

  Eigen::SparseMatrix<double> L;
  build_linear_system(s, L);

  igl::Timer t;

  // t.start();
  // solve
  Eigen::VectorXd Uc;
#ifndef CHOLMOD
  if (s.dim == 2) {
    SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Uc = solver.compute(L).solve(s.rhs);
  } else { // seems like CG performs much worse for 2D and way better for 3D
    Eigen::VectorXd guess(uv.rows() * s.dim);
    for (int i = 0; i < s.v_num; i++)
      for (int j = 0; j < s.dim; j++)
        guess(uv.rows() * j + i) = uv(i, j); // flatten vector
    ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper> cg;
    cg.setTolerance(1e-8);
    cg.compute(L);
    Uc = cg.solveWithGuess(s.rhs, guess);
  }
#else
  CholmodSimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  Uc = solver.compute(L).solve(s.rhs);
#endif
  for (int i = 0; i < s.dim; i++)
    uv.col(i) = Uc.block(i * s.v_n, 0, s.v_n, 1);

  // t.stop();
  // std::cerr << "solve: " << t.getElapsedTime() << std::endl;
}

IGL_INLINE void pre_calc(SLIMData &s) {
  if (!s.has_pre_calc) {
    s.v_n = s.v_num;
    s.f_n = s.f_num;

    if (s.F.cols() == 3) {
      s.dim = 2;
      Eigen::MatrixXd F1, F2, F3;
      igl::local_basis(s.V, s.F, F1, F2, F3);
      Eigen::SparseMatrix<double> G;
      igl::grad(s.V, s.F, G);
      Eigen::SparseMatrix<double> Face_Proj;

      auto face_proj = [](Eigen::MatrixXd &F) {
        std::vector<Eigen::Triplet<double>> IJV;
        int f_num = F.rows();
        for (int i = 0; i < F.rows(); i++) {
          IJV.push_back(Eigen::Triplet<double>(i, i, F(i, 0)));
          IJV.push_back(Eigen::Triplet<double>(i, i + f_num, F(i, 1)));
          IJV.push_back(Eigen::Triplet<double>(i, i + 2 * f_num, F(i, 2)));
        }
        Eigen::SparseMatrix<double> P(f_num, 3 * f_num);
        P.setFromTriplets(IJV.begin(), IJV.end());
        return P;
      };

      s.Dx = face_proj(F1) * G;
      s.Dy = face_proj(F2) * G;
    } else {
      s.dim = 3;
      Eigen::SparseMatrix<double> G;
      igl::grad(s.V, s.F, G,
                    s.mesh_improvement_3d /*use normal gradient, or one from a "regular" tet*/);
      s.Dx = G.block(0, 0, s.F.rows(), s.V.rows());
      s.Dy = G.block(s.F.rows(), 0, s.F.rows(), s.V.rows());
      s.Dz = G.block(2 * s.F.rows(), 0, s.F.rows(), s.V.rows());
    }

    s.W.resize(s.f_n, s.dim * s.dim);
    s.Dx.makeCompressed();
    s.Dy.makeCompressed();
    s.Dz.makeCompressed();
    s.Ri.resize(s.f_n, s.dim * s.dim);
    s.Ji.resize(s.f_n, s.dim * s.dim);
    s.rhs.resize(s.dim * s.v_num);

    // flattened weight matrix
    s.WGL_M.resize(s.dim * s.dim * s.f_n);
    for (int i = 0; i < s.dim * s.dim; i++)
      for (int j = 0; j < s.f_n; j++)
        s.WGL_M(i * s.f_n + j) = s.M(j);

    s.first_solve = true;
    s.has_pre_calc = true;
  }
}

IGL_INLINE void build_linear_system(SLIMData &s,
                                    Eigen::SparseMatrix<double> &L) {
  // formula (35) in paper
  std::vector<Eigen::Triplet<double>> IJV;

#ifdef SLIM_CACHED
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  if (s.A.rows() == 0) {
    s.A = Eigen::SparseMatrix<double>(s.dim * s.dim * s.f_n, s.dim * s.v_n);
    igl::sparse_cached_precompute(IJV, s.A_data, s.A);
  } else
    igl::sparse_cached(IJV, s.A_data, s.A);
#else
  Eigen::SparseMatrix<double> A(s.dim * s.dim * s.f_n, s.dim * s.v_n);
  slim_buildA(s.Dx, s.Dy, s.Dz, s.W, IJV);
  A.setFromTriplets(IJV.begin(), IJV.end());
  A.makeCompressed();
#endif

#ifdef SLIM_CACHED
#else
  Eigen::SparseMatrix<double> At = A.transpose();
  At.makeCompressed();
#endif

#ifdef SLIM_CACHED
  Eigen::SparseMatrix<double> id_m(s.A.cols(), s.A.cols());
#else
  Eigen::SparseMatrix<double> id_m(A.cols(), A.cols());
#endif

  id_m.setIdentity();

// add proximal penalty
#ifdef SLIM_CACHED
  s.AtA_data.W = s.WGL_M;
  if (s.AtA.rows() == 0)
    igl::AtA_cached_precompute(s.A, s.AtA_data, s.AtA);
  else
    igl::AtA_cached(s.A, s.AtA_data, s.AtA);

  L = s.AtA + s.proximal_p * id_m; // add also a proximal
  L.makeCompressed();

#else
  L = At * s.WGL_M.asDiagonal() * A +
      s.proximal_p * id_m; // add also a proximal term
  L.makeCompressed();
#endif

#ifdef SLIM_CACHED
  buildRhs(s, s.A);
#else
  buildRhs(s, A);
#endif

  Eigen::SparseMatrix<double> OldL = L;
  add_soft_constraints(s, L);
  L.makeCompressed();
}

IGL_INLINE void add_soft_constraints(SLIMData &s,
                                     Eigen::SparseMatrix<double> &L) {
  int v_n = s.v_num;
  for (int d = 0; d < s.dim; d++) {
    for (int i = 0; i < s.b.rows(); i++) {
      int v_idx = s.b(i);
      s.rhs(d * v_n + v_idx) += s.soft_const_p * s.bc(i, d); // rhs
      L.coeffRef(d * v_n + v_idx, d * v_n + v_idx) +=
          s.soft_const_p; // diagonal of matrix
    }
  }
}

IGL_INLINE double compute_slim_energy(SLIMData &s,
                                      const Eigen::MatrixXd &V_new) {
  compute_jacobians(s, V_new);
  return mapping_energy_with_jacobians(s.Ji, s.M, s.slim_energy, s.exp_factor) +
         compute_soft_const_energy(s, s.V, s.F, V_new);
}

IGL_INLINE double compute_soft_const_energy(SLIMData &s,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F,
                                            const Eigen::MatrixXd &V_o) {
  double e = 0;
  for (int i = 0; i < s.b.rows(); i++) {
    e += s.soft_const_p * (s.bc.row(i) - V_o.row(s.b(i))).squaredNorm();
  }
  return e;
}

IGL_INLINE void calculate_surface(SLIMData &s, const Eigen::MatrixXi &surface) {
  using namespace std;
  using namespace Eigen;
  s.m_surface.conservativeResize(surface.rows(), 4);
  s.m_surface.bottomRows(surface.rows()) = surface;
  s.m_surface_fn = s.m_surface.rows();
  std::unordered_set<uint32_t> uniqueset;
  for (uint32_t i = 0; i < s.m_surface_fn; ++i) {
    for (uint32_t j = 0; j < 4; ++j) {
      uniqueset.insert(surface(i, j));
    }
  }
  s.m_surface_vn = uniqueset.size();
  s.mapping_tetgen2scaf.assign(uniqueset.begin(), uniqueset.end());
  for (uint32_t i = 0; i < s.mapping_tetgen2scaf.size(); ++i) {
    s.mapping_scaf2tetgen.insert({s.mapping_tetgen2scaf[i], i});
  }
  return;
}

IGL_INLINE void calculate_scaf_tetgen(SLIMData &s) {
  using namespace Eigen;
  MatrixXd m_uv = s.V_o;
  auto mapping_t2s = s.mapping_tetgen2scaf;
  auto mapping_s2t = s.mapping_scaf2tetgen;

  tetgenio in, out;

  {
    VectorXd uv_max = m_uv.colwise().maxCoeff();
    VectorXd uv_min = m_uv.colwise().minCoeff();
    uv_max += Vector3d(1, 1, 1), uv_min -= Vector3d(1, 1, 1);

    in.numberofpoints = s.m_surface_vn + 8;
    in.pointlist = new double[in.numberofpoints * 3];
    for (uint32_t i = 0; i < s.m_surface_vn; ++i) {
      int scaf_num = mapping_t2s[i];
      in.pointlist[i * 3 + 0] = s.V_o(scaf_num, 0);
      in.pointlist[i * 3 + 1] = s.V_o(scaf_num, 1);
      in.pointlist[i * 3 + 2] = s.V_o(scaf_num, 2);
    }
    {
      int _tn = s.m_surface_vn;
      in.pointlist[(_tn + 0) * 3 + 0] = uv_min(0);
      in.pointlist[(_tn + 0) * 3 + 1] = uv_min(1);
      in.pointlist[(_tn + 0) * 3 + 2] = uv_min(2);

      in.pointlist[(_tn + 1) * 3 + 0] = uv_max(0);
      in.pointlist[(_tn + 1) * 3 + 1] = uv_min(1);
      in.pointlist[(_tn + 1) * 3 + 2] = uv_min(2);

      in.pointlist[(_tn + 2) * 3 + 0] = uv_max(0);
      in.pointlist[(_tn + 2) * 3 + 1] = uv_max(1);
      in.pointlist[(_tn + 2) * 3 + 2] = uv_min(2);

      in.pointlist[(_tn + 3) * 3 + 0] = uv_min(0);
      in.pointlist[(_tn + 3) * 3 + 1] = uv_max(1);
      in.pointlist[(_tn + 3) * 3 + 2] = uv_min(2);

      in.pointlist[(_tn + 4) * 3 + 0] = uv_min(0);
      in.pointlist[(_tn + 4) * 3 + 1] = uv_min(1);
      in.pointlist[(_tn + 4) * 3 + 2] = uv_max(2);

      in.pointlist[(_tn + 5) * 3 + 0] = uv_max(0);
      in.pointlist[(_tn + 5) * 3 + 1] = uv_min(1);
      in.pointlist[(_tn + 5) * 3 + 2] = uv_max(2);

      in.pointlist[(_tn + 6) * 3 + 0] = uv_max(0);
      in.pointlist[(_tn + 6) * 3 + 1] = uv_max(1);
      in.pointlist[(_tn + 6) * 3 + 2] = uv_max(2);

      in.pointlist[(_tn + 7) * 3 + 0] = uv_min(0);
      in.pointlist[(_tn + 7) * 3 + 1] = uv_max(1);
      in.pointlist[(_tn + 7) * 3 + 2] = uv_max(2);
    }
    in.numberoffacets = s.m_surface_fn + 6;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    for (uint32_t i = 0; i < in.numberoffacets; ++i) {
      in.facetlist[i].numberofpolygons = 1;
      in.facetlist[i].numberofholes = 0;
      in.facetlist[i].holelist = nullptr;
      in.facetlist[i].polygonlist = new tetgenio::polygon;
      auto p = in.facetlist[i].polygonlist;
      p->numberofvertices = 4;
      p->vertexlist = new int[4];
      if (i < s.m_surface_fn) {
        for (uint32_t j = 0; j < 4; ++j) {
          p->vertexlist[j] = mapping_s2t[s.m_surface(i, j)];
        }
      }
    }
    {
      int _fn = s.m_surface_fn;
      int _tn = s.m_surface_vn;

      in.facetlist[_fn + 0].polygonlist->vertexlist[0] = _tn + 0;
      in.facetlist[_fn + 0].polygonlist->vertexlist[1] = _tn + 1;
      in.facetlist[_fn + 0].polygonlist->vertexlist[2] = _tn + 2;
      in.facetlist[_fn + 0].polygonlist->vertexlist[3] = _tn + 3;

      in.facetlist[_fn + 1].polygonlist->vertexlist[0] = _tn + 0;
      in.facetlist[_fn + 1].polygonlist->vertexlist[1] = _tn + 4;
      in.facetlist[_fn + 1].polygonlist->vertexlist[2] = _tn + 5;
      in.facetlist[_fn + 1].polygonlist->vertexlist[3] = _tn + 1;

      in.facetlist[_fn + 2].polygonlist->vertexlist[0] = _tn + 1;
      in.facetlist[_fn + 2].polygonlist->vertexlist[1] = _tn + 5;
      in.facetlist[_fn + 2].polygonlist->vertexlist[2] = _tn + 6;
      in.facetlist[_fn + 2].polygonlist->vertexlist[3] = _tn + 2;

      in.facetlist[_fn + 3].polygonlist->vertexlist[0] = _tn + 2;
      in.facetlist[_fn + 3].polygonlist->vertexlist[1] = _tn + 6;
      in.facetlist[_fn + 3].polygonlist->vertexlist[2] = _tn + 7;
      in.facetlist[_fn + 3].polygonlist->vertexlist[3] = _tn + 3;

      in.facetlist[_fn + 4].polygonlist->vertexlist[0] = _tn + 3;
      in.facetlist[_fn + 4].polygonlist->vertexlist[1] = _tn + 7;
      in.facetlist[_fn + 4].polygonlist->vertexlist[2] = _tn + 4;
      in.facetlist[_fn + 4].polygonlist->vertexlist[3] = _tn + 0;

      in.facetlist[_fn + 5].polygonlist->vertexlist[0] = _tn + 4;
      in.facetlist[_fn + 5].polygonlist->vertexlist[1] = _tn + 7;
      in.facetlist[_fn + 5].polygonlist->vertexlist[2] = _tn + 6;
      in.facetlist[_fn + 5].polygonlist->vertexlist[3] = _tn + 5;
    }
    in.numberofholes = 1;
    in.holelist = new double[3];
    {
      auto _t = (s.V_o.row(s.F(0, 0)) + s.V_o.row(s.F(0, 1)) +
                 s.V_o.row(s.F(0, 2)) + s.V_o.row(s.F(0, 3))) /
                4;
      in.holelist[0] = _t(0);
      in.holelist[1] = _t(1);
      in.holelist[2] = _t(2);
    }
  }

  /*
  std::ofstream ftmp;
  ftmp.open("ftmp_mesh.vtk");
  write_matrix_to_vtk_tet(s.m_V, s.m_T, ftmp);
  ftmp.close();
  */
  tetgenbehavior beh;
  beh.plc = 1;
  beh.nobisect = 1;
  tetrahedralize(&beh, &in, &out);

  // generate by tetgen
  int addinVnum1 = out.numberofpoints - in.numberofpoints;
  // generate by me 确定现在加还是之前加
  int addinVnum2 = 8;

  for (uint32_t i = 0; i < addinVnum2 + addinVnum1; ++i) {
    mapping_t2s.push_back(s.V_o.rows() + i);
    mapping_s2t[s.V_o.rows() + i] = s.m_surface_vn + i;
  }

  s.w_V.conservativeResize(s.V_o.rows() + addinVnum1 + addinVnum2, 3);
  s.w_V.topRows(s.V_o.rows()) = s.V_o;
  for (uint32_t i = 0; i < addinVnum2 + addinVnum1; ++i) {
    int _t = s.m_surface_vn + i;
    s.w_V.row(s.V_o.rows() + i) =
        Vector3d(out.pointlist[_t * 3 + 0], out.pointlist[_t * 3 + 1],
                 out.pointlist[_t * 3 + 2]);
  }

  s.w_T.conservativeResize(out.numberoftetrahedra, 4);
  for (uint32_t i = 0; i < out.numberoftetrahedra; ++i) {
    s.w_T.row(i) = Vector4i(mapping_t2s[out.tetrahedronlist[i * 4 + 0]],
                            mapping_t2s[out.tetrahedronlist[i * 4 + 1]],
                            mapping_t2s[out.tetrahedronlist[i * 4 + 2]],
                            mapping_t2s[out.tetrahedronlist[i * 4 + 3]]);
  }

  // std::ofstream ftmp;

  /*
  ftmp.open("ftmp_scaf.vtk");
  write_matrix_to_vtk_tet(s.w_V, s.s_T, ftmp);
  ftmp.close();

  update_scaffold(s);

  igl::my_scaf::compute_my_mesh_grad_inv_matrix(s.w_V, s.s_T, s.s_Grad);
  */
}

IGL_INLINE void buildRhs(SLIMData &s, const Eigen::SparseMatrix<double> &A) {
  Eigen::VectorXd f_rhs(s.dim * s.dim * s.f_n);
  f_rhs.setZero();
  if (s.dim == 2) {
    /*b = [W11*R11 + W12*R21; (formula (36))
         W11*R12 + W12*R22;
         W21*R11 + W22*R21;
         W21*R12 + W22*R22];*/
    for (int i = 0; i < s.f_n; i++) {
      f_rhs(i + 0 * s.f_n) = s.W(i, 0) * s.Ri(i, 0) + s.W(i, 1) * s.Ri(i, 1);
      f_rhs(i + 1 * s.f_n) = s.W(i, 0) * s.Ri(i, 2) + s.W(i, 1) * s.Ri(i, 3);
      f_rhs(i + 2 * s.f_n) = s.W(i, 2) * s.Ri(i, 0) + s.W(i, 3) * s.Ri(i, 1);
      f_rhs(i + 3 * s.f_n) = s.W(i, 2) * s.Ri(i, 2) + s.W(i, 3) * s.Ri(i, 3);
    }
  } else {
    /*b = [W11*R11 + W12*R21 + W13*R31;
         W11*R12 + W12*R22 + W13*R32;
         W11*R13 + W12*R23 + W13*R33;
         W21*R11 + W22*R21 + W23*R31;
         W21*R12 + W22*R22 + W23*R32;
         W21*R13 + W22*R23 + W23*R33;
         W31*R11 + W32*R21 + W33*R31;
         W31*R12 + W32*R22 + W33*R32;
         W31*R13 + W32*R23 + W33*R33;];*/
    for (int i = 0; i < s.f_n; i++) {
      f_rhs(i + 0 * s.f_n) = s.W(i, 0) * s.Ri(i, 0) + s.W(i, 1) * s.Ri(i, 1) +
                             s.W(i, 2) * s.Ri(i, 2);
      f_rhs(i + 1 * s.f_n) = s.W(i, 0) * s.Ri(i, 3) + s.W(i, 1) * s.Ri(i, 4) +
                             s.W(i, 2) * s.Ri(i, 5);
      f_rhs(i + 2 * s.f_n) = s.W(i, 0) * s.Ri(i, 6) + s.W(i, 1) * s.Ri(i, 7) +
                             s.W(i, 2) * s.Ri(i, 8);
      f_rhs(i + 3 * s.f_n) = s.W(i, 3) * s.Ri(i, 0) + s.W(i, 4) * s.Ri(i, 1) +
                             s.W(i, 5) * s.Ri(i, 2);
      f_rhs(i + 4 * s.f_n) = s.W(i, 3) * s.Ri(i, 3) + s.W(i, 4) * s.Ri(i, 4) +
                             s.W(i, 5) * s.Ri(i, 5);
      f_rhs(i + 5 * s.f_n) = s.W(i, 3) * s.Ri(i, 6) + s.W(i, 4) * s.Ri(i, 7) +
                             s.W(i, 5) * s.Ri(i, 8);
      f_rhs(i + 6 * s.f_n) = s.W(i, 6) * s.Ri(i, 0) + s.W(i, 7) * s.Ri(i, 1) +
                             s.W(i, 8) * s.Ri(i, 2);
      f_rhs(i + 7 * s.f_n) = s.W(i, 6) * s.Ri(i, 3) + s.W(i, 7) * s.Ri(i, 4) +
                             s.W(i, 8) * s.Ri(i, 5);
      f_rhs(i + 8 * s.f_n) = s.W(i, 6) * s.Ri(i, 6) + s.W(i, 7) * s.Ri(i, 7) +
                             s.W(i, 8) * s.Ri(i, 8);
    }
  }
  Eigen::VectorXd uv_flat(s.dim * s.v_n);
  for (int i = 0; i < s.dim; i++)
    for (int j = 0; j < s.v_n; j++)
      uv_flat(s.v_n * i + j) = s.V_o(j, i);

  s.rhs = (f_rhs.transpose() * s.WGL_M.asDiagonal() * A).transpose() +
          s.proximal_p * uv_flat;
}

IGL_INLINE void slim_update_weights_and_closest_rotations_with_jacobians(
    const Eigen::MatrixXd &Ji, MappingEnergyType slim_energy, double exp_factor,
    Eigen::MatrixXd &W, Eigen::MatrixXd &Ri) {
  const double eps = 1e-8;
  double exp_f = exp_factor;
  const int dim = (Ji.cols() == 4 ? 2 : 3);

  if (dim == 2) {
    for (int i = 0; i < Ji.rows(); ++i) {
      typedef Eigen::Matrix2d Mat2;
      typedef Eigen::Matrix<double, 2, 2, Eigen::RowMajor> RMat2;
      typedef Eigen::Vector2d Vec2;
      Mat2 ji, ri, ti, ui, vi;
      Vec2 sing;
      Vec2 closest_sing_vec;
      RMat2 mat_W;
      Vec2 m_sing_new;
      double s1, s2;

      ji(0, 0) = Ji(i, 0);
      ji(0, 1) = Ji(i, 1);
      ji(1, 0) = Ji(i, 2);
      ji(1, 1) = Ji(i, 3);

      polar_svd(ji, ri, ti, ui, sing, vi);

      s1 = sing(0);
      s2 = sing(1);

      // Update Weights according to energy
      switch (slim_energy) {
      case MappingEnergyType::ARAP: {
        m_sing_new << 1, 1;
        break;
      }
      case MappingEnergyType::SYMMETRIC_DIRICHLET: {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      case MappingEnergyType::LOG_ARAP: {
        double s1_g = 2 * (log(s1) / s1);
        double s2_g = 2 * (log(s2) / s2);
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      case MappingEnergyType::CONFORMAL: {
        double s1_g = 1 / (2 * s2) - s2 / (2 * pow(s1, 2));
        double s2_g = 1 / (2 * s1) - s1 / (2 * pow(s2, 2));

        double geo_avg = sqrt(s1 * s2);
        double s1_min = geo_avg;
        double s2_min = geo_avg;

        m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))),
            sqrt(s2_g / (2 * (s2 - s2_min)));

        // change local step
        closest_sing_vec << s1_min, s2_min;
        ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
        break;
      }
      case MappingEnergyType::EXP_CONFORMAL: {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));

        double geo_avg = sqrt(s1 * s2);
        double s1_min = geo_avg;
        double s2_min = geo_avg;

        double in_exp = exp_f * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
        double exp_thing = exp(in_exp);

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      case MappingEnergyType::EXP_SYMMETRIC_DIRICHLET: {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));

        double in_exp =
            exp_f * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2));
        double exp_thing = exp(in_exp);

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
        break;
      }
      default:
        assert(false);
      }

      if (std::abs(s1 - 1) < eps)
        m_sing_new(0) = 1;
      if (std::abs(s2 - 1) < eps)
        m_sing_new(1) = 1;
      mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

      W.row(i) = Eigen::Map<Eigen::Matrix<double, 1, 4, Eigen::RowMajor>>(
          mat_W.data());
      // 2) Update local step (doesn't have to be a rotation, for instance in
      // case of conformal energy)
      Ri.row(i) =
          Eigen::Map<Eigen::Matrix<double, 1, 4, Eigen::RowMajor>>(ri.data());
    }
  } else {
    typedef Eigen::Matrix<double, 3, 1> Vec3;
    typedef Eigen::Matrix<double, 3, 3, Eigen::ColMajor> Mat3;
    typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> RMat3;
    Mat3 ji;
    Vec3 m_sing_new;
    Vec3 closest_sing_vec;
    const double sqrt_2 = sqrt(2);
    for (int i = 0; i < Ji.rows(); ++i) {
      ji << Ji(i, 0), Ji(i, 1), Ji(i, 2), Ji(i, 3), Ji(i, 4), Ji(i, 5),
          Ji(i, 6), Ji(i, 7), Ji(i, 8);

      Mat3 ri, ti, ui, vi;
      Vec3 sing;
      polar_svd(ji, ri, ti, ui, sing, vi);

      double s1 = sing(0);
      double s2 = sing(1);
      double s3 = sing(2);

      // 1) Update Weights
      switch (slim_energy) {
      case MappingEnergyType::ARAP: {
        m_sing_new << 1, 1, 1;
        break;
      }
      case MappingEnergyType::LOG_ARAP: {
        double s1_g = 2 * (log(s1) / s1);
        double s2_g = 2 * (log(s2) / s2);
        double s3_g = 2 * (log(s3) / s3);
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
            sqrt(s3_g / (2 * (s3 - 1)));
        break;
      }
      case MappingEnergyType::SYMMETRIC_DIRICHLET: {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));
        double s3_g = 2 * (s3 - pow(s3, -3));
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
            sqrt(s3_g / (2 * (s3 - 1)));
        break;
      }
      case MappingEnergyType::EXP_SYMMETRIC_DIRICHLET: {
        double s1_g = 2 * (s1 - pow(s1, -3));
        double s2_g = 2 * (s2 - pow(s2, -3));
        double s3_g = 2 * (s3 - pow(s3, -3));
        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
            sqrt(s3_g / (2 * (s3 - 1)));

        double in_exp = exp_f * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) +
                                 pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
        double exp_thing = exp(in_exp);

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;
        s3_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
            sqrt(s3_g / (2 * (s3 - 1)));

        break;
      }
      case MappingEnergyType::CONFORMAL: {
        double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

        double s1_g =
            (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) /
            common_div;
        double s2_g =
            (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) /
            common_div;
        double s3_g =
            (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) /
            common_div;

        double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
        double s1_min = closest_s;
        double s2_min = closest_s;
        double s3_min = closest_s;

        m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))),
            sqrt(s2_g / (2 * (s2 - s2_min))), sqrt(s3_g / (2 * (s3 - s3_min)));

        // change local step
        closest_sing_vec << s1_min, s2_min, s3_min;
        ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
        break;
      }
      case MappingEnergyType::EXP_CONFORMAL: {
        // E_conf = (s1^2 + s2^2 + s3^2)/(3*(s1*s2*s3)^(2/3) )
        // dE_conf/ds1 = (-2*(s2*s3)*(s2^2+s3^2 -2*s1^2) ) /
        // (9*(s1*s2*s3)^(5/3)) Argmin E_conf(s1): s1 = sqrt(s1^2+s2^2)/sqrt(2)
        double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

        double s1_g =
            (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) /
            common_div;
        double s2_g =
            (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) /
            common_div;
        double s3_g =
            (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) /
            common_div;

        double in_exp = exp_f * ((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) /
                                 (3 * pow((s1 * s2 * s3), 2. / 3)));
        ;
        double exp_thing = exp(in_exp);

        double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
        double s1_min = closest_s;
        double s2_min = closest_s;
        double s3_min = closest_s;

        s1_g *= exp_thing * exp_f;
        s2_g *= exp_thing * exp_f;
        s3_g *= exp_thing * exp_f;

        m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))),
            sqrt(s2_g / (2 * (s2 - s2_min))), sqrt(s3_g / (2 * (s3 - s3_min)));

        // change local step
        closest_sing_vec << s1_min, s2_min, s3_min;
        ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
        break;
      }
      default:
        assert(false);
      }
      if (std::abs(s1 - 1) < eps)
        m_sing_new(0) = 1;
      if (std::abs(s2 - 1) < eps)
        m_sing_new(1) = 1;
      if (std::abs(s3 - 1) < eps)
        m_sing_new(2) = 1;
      RMat3 mat_W;
      mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

      W.row(i) = Eigen::Map<Eigen::Matrix<double, 1, 9, Eigen::RowMajor>>(
          mat_W.data());
      // 2) Update closest rotations (not rotations in case of conformal energy)
      Ri.row(i) =
          Eigen::Map<Eigen::Matrix<double, 1, 9, Eigen::RowMajor>>(ri.data());
    } // for loop end

  } // if dim end
}

IGL_INLINE void slim_buildA(const Eigen::SparseMatrix<double> &Dx,
                            const Eigen::SparseMatrix<double> &Dy,
                            const Eigen::SparseMatrix<double> &Dz,
                            const Eigen::MatrixXd &W,
                            std::vector<Eigen::Triplet<double>> &IJV) {
  const int dim = (W.cols() == 4) ? 2 : 3;
  const int f_n = W.rows();
  const int v_n = Dx.cols();

  // formula (35) in paper
  if (dim == 2) {
    IJV.reserve(4 * (Dx.outerSize() + Dy.outerSize()));

    /*A = [W11*Dx, W12*Dx;
          W11*Dy, W12*Dy;
          W21*Dx, W22*Dx;
          W21*Dy, W22*Dy];*/
    for (int k = 0; k < Dx.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it) {
        int dx_r = it.row();
        int dx_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * W(dx_r, 0)));
        IJV.push_back(
            Eigen::Triplet<double>(dx_r, v_n + dx_c, val * W(dx_r, 1)));

        IJV.push_back(
            Eigen::Triplet<double>(2 * f_n + dx_r, dx_c, val * W(dx_r, 2)));
        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dx_r, v_n + dx_c,
                                             val * W(dx_r, 3)));
      }
    }

    for (int k = 0; k < Dy.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it) {
        int dy_r = it.row();
        int dy_c = it.col();
        double val = it.value();

        IJV.push_back(
            Eigen::Triplet<double>(f_n + dy_r, dy_c, val * W(dy_r, 0)));
        IJV.push_back(
            Eigen::Triplet<double>(f_n + dy_r, v_n + dy_c, val * W(dy_r, 1)));

        IJV.push_back(
            Eigen::Triplet<double>(3 * f_n + dy_r, dy_c, val * W(dy_r, 2)));
        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dy_r, v_n + dy_c,
                                             val * W(dy_r, 3)));
      }
    }
  } else {

    /*A = [W11*Dx, W12*Dx, W13*Dx;
            W11*Dy, W12*Dy, W13*Dy;
            W11*Dz, W12*Dz, W13*Dz;
            W21*Dx, W22*Dx, W23*Dx;
            W21*Dy, W22*Dy, W23*Dy;
            W21*Dz, W22*Dz, W23*Dz;
            W31*Dx, W32*Dx, W33*Dx;
            W31*Dy, W32*Dy, W33*Dy;
            W31*Dz, W32*Dz, W33*Dz;];*/
    IJV.reserve(9 * (Dx.outerSize() + Dy.outerSize() + Dz.outerSize()));
    for (int k = 0; k < Dx.outerSize(); k++) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it) {
        int dx_r = it.row();
        int dx_c = it.col();
        double val = it.value();

        IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, val * W(dx_r, 0)));
        IJV.push_back(
            Eigen::Triplet<double>(dx_r, v_n + dx_c, val * W(dx_r, 1)));
        IJV.push_back(
            Eigen::Triplet<double>(dx_r, 2 * v_n + dx_c, val * W(dx_r, 2)));

        IJV.push_back(
            Eigen::Triplet<double>(3 * f_n + dx_r, dx_c, val * W(dx_r, 3)));
        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r, v_n + dx_c,
                                             val * W(dx_r, 4)));
        IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r, 2 * v_n + dx_c,
                                             val * W(dx_r, 5)));

        IJV.push_back(
            Eigen::Triplet<double>(6 * f_n + dx_r, dx_c, val * W(dx_r, 6)));
        IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r, v_n + dx_c,
                                             val * W(dx_r, 7)));
        IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r, 2 * v_n + dx_c,
                                             val * W(dx_r, 8)));
      }
    }

    for (int k = 0; k < Dy.outerSize(); k++) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it) {
        int dy_r = it.row();
        int dy_c = it.col();
        double val = it.value();

        IJV.push_back(
            Eigen::Triplet<double>(f_n + dy_r, dy_c, val * W(dy_r, 0)));
        IJV.push_back(
            Eigen::Triplet<double>(f_n + dy_r, v_n + dy_c, val * W(dy_r, 1)));
        IJV.push_back(Eigen::Triplet<double>(f_n + dy_r, 2 * v_n + dy_c,
                                             val * W(dy_r, 2)));

        IJV.push_back(
            Eigen::Triplet<double>(4 * f_n + dy_r, dy_c, val * W(dy_r, 3)));
        IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r, v_n + dy_c,
                                             val * W(dy_r, 4)));
        IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r, 2 * v_n + dy_c,
                                             val * W(dy_r, 5)));

        IJV.push_back(
            Eigen::Triplet<double>(7 * f_n + dy_r, dy_c, val * W(dy_r, 6)));
        IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r, v_n + dy_c,
                                             val * W(dy_r, 7)));
        IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r, 2 * v_n + dy_c,
                                             val * W(dy_r, 8)));
      }
    }

    for (int k = 0; k < Dz.outerSize(); k++) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, k); it; ++it) {
        int dz_r = it.row();
        int dz_c = it.col();
        double val = it.value();

        IJV.push_back(
            Eigen::Triplet<double>(2 * f_n + dz_r, dz_c, val * W(dz_r, 0)));
        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r, v_n + dz_c,
                                             val * W(dz_r, 1)));
        IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r, 2 * v_n + dz_c,
                                             val * W(dz_r, 2)));

        IJV.push_back(
            Eigen::Triplet<double>(5 * f_n + dz_r, dz_c, val * W(dz_r, 3)));
        IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r, v_n + dz_c,
                                             val * W(dz_r, 4)));
        IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r, 2 * v_n + dz_c,
                                             val * W(dz_r, 5)));

        IJV.push_back(
            Eigen::Triplet<double>(8 * f_n + dz_r, dz_c, val * W(dz_r, 6)));
        IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r, v_n + dz_c,
                                             val * W(dz_r, 7)));
        IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r, 2 * v_n + dz_c,
                                             val * W(dz_r, 8)));
      }
    }
  }
}

IGL_INLINE void slim_precompute(const Eigen::MatrixXd &V,
                                const Eigen::MatrixXi &F,
                                const Eigen::MatrixXd &V_init,
                                const Eigen::MatrixXi &surface, SLIMData &data,
                                MappingEnergyType slim_energy,
                                const Eigen::VectorXi &b,
                                const Eigen::MatrixXd &bc, double soft_p) {

  data.V = V;
  data.F = F;
  data.V_o = V_init;

  data.v_num = V.rows();
  data.f_num = F.rows();

  data.slim_energy = slim_energy;

  data.b = b;
  data.bc = bc;
  data.soft_const_p = soft_p;

  data.proximal_p = 0.0001;

  igl::doublearea(V, F, data.M);
  data.M /= 2.;
  data.mesh_area = data.M.sum();
  data.mesh_improvement_3d =
      false; // whether to use a jacobian derived from a real mesh or an
             // abstract regular mesh (used for mesh improvement)
  data.exp_factor = 1.0; // param used only for exponential energies (e.g
                         // exponential symmetric dirichlet)

  assert(F.cols() == 3 || F.cols() == 4);
  calculate_surface(data, surface);
  calculate_scaf_tetgen(data);

  pre_calc(data);
  data.energy = compute_slim_energy(data, data.V_o) / data.mesh_area;
}

IGL_INLINE double flip_avoiding_line_search_scaf(
    const Eigen::MatrixXd &w_V, const Eigen::MatrixXi &w_T,
    const Eigen::MatrixXi F, Eigen::MatrixXd &cur_v,
    const Eigen::MatrixXd &dst_v,
    std::function<double(Eigen::MatrixXd &)> energy, double cur_energy) {
  using namespace std;
  Eigen::MatrixXd d = dst_v - cur_v;
  Eigen::MatrixXd w_d(d.rows(), d.cols());
  w_d.topRows(d.rows()) = d;

  double min_step_to_singularity =
      igl::flip_avoiding::compute_max_step_from_singularities(w_V, w_T, w_d);
  double max_step_size = std::min(1., min_step_to_singularity * 0.8);

  return igl::line_search(cur_v, d, max_step_size, energy, cur_energy);
}

IGL_INLINE Eigen::MatrixXd slim_solve(SLIMData &data, int iter_num) {
  for (int i = 0; i < iter_num; i++) {
    Eigen::MatrixXd dest_res;
    dest_res = data.V_o;

    calculate_scaf_tetgen(data);

    // Solve Weighted Proxy
    update_weights_and_closest_rotations(data, dest_res);
    solve_weighted_arap(data, data.V, data.F, dest_res, data.b, data.bc);

    double old_energy = data.energy;

    std::function<double(Eigen::MatrixXd &)> compute_energy =
        [&](Eigen::MatrixXd &aaa) { return compute_slim_energy(data, aaa); };

    data.energy = flip_avoiding_line_search_scaf(
                      data.w_V, data.w_T, data.F, data.V_o, dest_res,
                      compute_energy, data.energy * data.mesh_area) /
                  data.mesh_area;
  }
  return data.V_o;
}
} // namespace my_scaf
} // namespace igl

#endif

// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "my_scaf.h"
#include "ovmwrap.h"
#include "printEigen.h"
//#include "triangulate.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/SparseQR>
#include <igl/PI.h>
#include <igl/Timer.h>
#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/doublearea.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/flipped_triangles.h>
#include <igl/grad.h>
#include <igl/harmonic.h>
#include <igl/local_basis.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/mapping_energy_with_jacobians.h>
#include <igl/polar_svd.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/slim.h>
#include <igl/volume.h>

#include <algorithm>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>

#include <tetgen.h>

namespace igl {
namespace my_scaf {
static IGL_INLINE void compute_my_mesh_grad_inv_matrix(const Eigen::MatrixXd &V,
                                                       const Eigen::MatrixXi &T,
                                                       Eigen::MatrixXd &G) {
  using namespace Eigen;
  uint32_t cnum = T.rows(), vnum = V.rows();
  G.conservativeResize(cnum * 3, 3);

  Matrix3d P;
  Vector3d dp1, dp2, dp3;
  for (uint32_t i = 0; i < cnum; ++i) {
    dp1 = V.row(T(i, 1)) - V.row(T(i, 0));
    dp2 = V.row(T(i, 2)) - V.row(T(i, 0));
    dp3 = V.row(T(i, 3)) - V.row(T(i, 0));
    P.row(0) = dp1, P.row(1) = dp2, P.row(2) = dp3;
    P.transposeInPlace();
#ifndef NDEBUG
    assert(P.determinant() >= 0);
#endif
    G.block(i * 3, 0, 3, 3) = P.inverse(); // P.block(0, 0, 3, 3);
  }
}

static IGL_INLINE void compute_my_mesh_grad_inv_matrix_standard(
    const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, Eigen::MatrixXd &G) {
  using namespace Eigen;
  uint32_t cnum = T.rows(), vnum = V.rows();
  G.conservativeResize(cnum * 3, 3);

  // Matrix3d P;
  Vector3d dp1, dp2, dp3;
  for (uint32_t i = 0; i < cnum; ++i) {
    dp1 = V.row(T(i, 1)) - V.row(T(i, 0));
    dp2 = V.row(T(i, 2)) - V.row(T(i, 0));
    dp3 = V.row(T(i, 3)) - V.row(T(i, 0));
    auto _len = (dp1.norm() + dp2.norm() + dp3.norm()) /
                3; // pow(fabs(dp1.cross(dp2).dot(dp3)), 1.0 / 3);
    Matrix3d P = (1.0 / _len) * MatrixXd::Identity(3, 3);

    G.block(i * 3, 0, 3, 3) = P; // P.block(0, 0, 3, 3);
  }
}

static IGL_INLINE void update_scaffold(igl::my_scaf::SCAFData &s) {
  s.mv_num = s.m_V.rows();
  s.mf_num = s.m_T.rows();

  // s.v_num = s.w_uv.rows();
  s.sv_num = s.w_V.rows() - s.m_V.rows();
  s.sf_num = s.s_T.rows();

  // s.sv_num = s.s_V.rows();
  // s.f_num = s.sf_num + s.mf_num;

  s.s_M = Eigen::VectorXd::Constant(s.sf_num, s.scaffold_factor);
}

static IGL_INLINE void mesh_improve(igl::my_scaf::SCAFData &s) {
  using namespace Eigen;
  MatrixXd m_uv = s.m_V;
  s.mapping_t2s = s.mapping_tetgen2scaf;
  s.mapping_s2t = s.mapping_scaf2tetgen;

  tetgenio in, out;

  {
    VectorXd uv_max = m_uv.colwise().maxCoeff();
    VectorXd uv_min = m_uv.colwise().minCoeff();
    uv_max += Vector3d(1, 1, 1), uv_min -= Vector3d(1, 1, 1);

    in.numberofpoints = s.m_surface_vn + 8;
    in.pointlist = new double[in.numberofpoints * 3];
    for (uint32_t i = 0; i < s.m_surface_vn; ++i) {
      int scaf_num = s.mapping_t2s[i];
      in.pointlist[i * 3 + 0] = s.m_V(scaf_num, 0);
      in.pointlist[i * 3 + 1] = s.m_V(scaf_num, 1);
      in.pointlist[i * 3 + 2] = s.m_V(scaf_num, 2);
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
          p->vertexlist[j] = s.mapping_s2t[s.m_surface(i, j)];
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
      auto _t = (s.m_V.row(s.m_T(0, 0)) + s.m_V.row(s.m_T(0, 1)) +
                 s.m_V.row(s.m_T(0, 2)) + s.m_V.row(s.m_T(0, 3))) /
                4;
      in.holelist[0] = _t(0);
      in.holelist[1] = _t(1);
      in.holelist[2] = _t(2);
    }
  }

  std::ofstream ftmp;
  ftmp.open("ftmp_mesh.vtk");
  write_matrix_to_vtk_tet(s.m_V, s.m_T, ftmp);
  ftmp.close();

  tetgenbehavior beh;
  beh.plc = 1;
  beh.nobisect = 1;
  tetrahedralize(&beh, &in, &out);

  // generate by tetgen
  int addinVnum1 = out.numberofpoints - in.numberofpoints;
  // generate by me 确定现在加还是之前加
  int addinVnum2 = 8;

  for (uint32_t i = 0; i < addinVnum2 + addinVnum1; ++i) {
    s.mapping_t2s.push_back(s.mv_num + i);
    s.mapping_s2t[s.mv_num + i] = s.m_surface_vn + i;
  }

  s.w_V.conservativeResize(s.mv_num + addinVnum1 + addinVnum2, 3);
  s.w_V.topRows(s.mv_num) = s.m_V;
  for (uint i = 0; i < addinVnum2 + addinVnum1; ++i) {
    int _t = s.m_surface_vn + i;
    s.w_V.row(s.mv_num + i) =
        Vector3d(out.pointlist[_t * 3 + 0], out.pointlist[_t * 3 + 1],
                 out.pointlist[_t * 3 + 2]);
  }

  s.s_T.conservativeResize(out.numberoftetrahedra, 4);
  for (uint32_t i = 0; i < out.numberoftetrahedra; ++i) {
    s.s_T.row(i) = Vector4i(s.mapping_t2s[out.tetrahedronlist[i * 4 + 0]],
                            s.mapping_t2s[out.tetrahedronlist[i * 4 + 1]],
                            s.mapping_t2s[out.tetrahedronlist[i * 4 + 2]],
                            s.mapping_t2s[out.tetrahedronlist[i * 4 + 3]]);
  }

  // std::ofstream ftmp;

  ftmp.open("ftmp_scaf.vtk");
  write_matrix_to_vtk_tet(s.w_V, s.s_T, ftmp);
  ftmp.close();

  update_scaffold(s);

  igl::my_scaf::compute_my_mesh_grad_inv_matrix(s.w_V, s.s_T, s.s_Grad);
}

static IGL_INLINE void add_new_patch(igl::my_scaf::SCAFData &s,
                                     const Eigen::MatrixXd &V_ref,
                                     const Eigen::MatrixXi &F_ref,
                                     const Eigen::MatrixXi &surface,
                                     const Eigen::MatrixXd &uv_init) {
  using namespace std;
  using namespace Eigen;

  assert(uv_init.rows() != 0);
  Eigen::VectorXd M;
  igl::volume(V_ref, F_ref, M);
  s.mesh_measure += M.sum();

  s.m_M.conservativeResize(s.mf_num + M.size());
  s.m_M.bottomRows(M.size()) = M;

  s.m_T.conservativeResize(s.mf_num + F_ref.rows(), 4);
  s.m_T.bottomRows(F_ref.rows()) = F_ref;
  s.mf_num += F_ref.rows();

  s.m_Vref.conservativeResize(s.mv_num + V_ref.rows(), 3);
  s.m_Vref.bottomRows(V_ref.rows()) = V_ref;
  s.mv_num += V_ref.rows();

  s.m_V.conservativeResize(s.m_Vref.rows(), 3);
  s.m_V.bottomRows(s.m_Vref.rows()) = uv_init;

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

  mesh_improve(s);
}

static IGL_INLINE double
compute_soft_constraint_energy(const SCAFData &s, const Eigen::MatrixXd &newV) {
  double e = 0;
  for (auto const &x : s.soft_cons)
    e += s.soft_const_p * (x.second - newV.row(x.first)).squaredNorm();

  return e;
}

IGL_INLINE double compute_energy(SCAFData &s, const Eigen::MatrixXd &w_V,
                                 bool whole) {
  using namespace Eigen;
  double energy = 0;
  if (w_V.rows() != s.sv_num + s.mv_num)
    assert(!whole);
  Matrix3d P, Q, J;
  Vector3d dp1, dp2, dp3;
  for (uint32_t i = 0; i < s.mf_num; ++i) {
    dp1 = w_V.row(s.m_T(i, 1)) - w_V.row(s.m_T(i, 0));
    dp2 = w_V.row(s.m_T(i, 2)) - w_V.row(s.m_T(i, 0));
    dp3 = w_V.row(s.m_T(i, 3)) - w_V.row(s.m_T(i, 0));
    Q.row(0) = dp1, Q.row(1) = dp2, Q.row(2) = dp3;
    Q.transposeInPlace();
    P = s.m_GradRef.block(i * 3, 0, 3, 3);
    J = Q * P;
    auto det = J.determinant();
    Matrix3d ri, vi, ui, ti;
    Vector3d sing;
    polar_svd(J, ri, ti, ui, sing, vi);
    double s1 = sing(0);
    double s2 = sing(1);
    double s3 = sing(2);

    energy += s.m_M(i) * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) +
                          pow(s3, 2) + pow(s3, -2));
  }

  if (whole) {
    for (uint32_t i = 0; i < s.sf_num; ++i) {
      dp1 = w_V.row(s.s_T(i, 1)) - w_V.row(s.s_T(i, 0));
      dp2 = w_V.row(s.s_T(i, 2)) - w_V.row(s.s_T(i, 0));
      dp3 = w_V.row(s.s_T(i, 3)) - w_V.row(s.s_T(i, 0));
      Q.row(0) = dp1, Q.row(1) = dp2, Q.row(2) = dp3;
      Q.transposeInPlace();
      P = s.s_Grad.block(i * 3, 0, 3, 3);
      J = Q * P;
      auto det = J.determinant();
      Matrix3d ri, vi, ui, ti;
      Vector3d sing;
      polar_svd(J, ri, ti, ui, sing, vi);
      double s1 = sing(0);
      double s2 = sing(1);
      double s3 = sing(2);

      energy += s.s_M(i) * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) +
                            pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
    }
  }
  energy += compute_soft_constraint_energy(s, w_V);
  return energy;
}

IGL_INLINE void scaf_build_element_linear_system(
    Eigen::MatrixXd &tmpA, Eigen::VectorXd &tmpb, const Eigen::Matrix3d &P1,
    const Eigen::Vector3d &v0, const Eigen::Vector3d &v1,
    const Eigen::Vector3d &v2, const Eigen::Vector3d &v3) {
  using namespace Eigen;
  tmpA.resize(12, 12);
  tmpb.resize(12);
  Matrix3d Q;
  Q.row(0) = v1 - v0;
  Q.row(1) = v2 - v0;
  Q.row(2) = v3 - v0;
  Q.transposeInPlace();
  auto J = Q * P1;
  MatrixXd P(4, 3);
  P.block(1, 0, 3, 3) = P1;
  {
    auto s0 = -P(1, 0) - P(2, 0) - P(3, 0), s1 = -P(1, 1) - P(2, 1) - P(3, 1),
         s2 = -P(1, 2) - P(2, 2) - P(3, 2);
    P.row(0) = Vector3d(s0, s1, s2);
  }
  Matrix3d ui, vi, ri, ti;
  Vector3d sing;
  polar_svd(J, ri, ti, ui, sing, vi);
  Vector3d new_sing;
  double new_s0, new_s1, new_s2;
  new_s0 = sqrt(1 + 1 / sing(0) + 1 / pow(sing(0), 2) + 1 / pow(sing(0), 3));
  new_s1 = sqrt(1 + 1 / sing(1) + 1 / pow(sing(1), 2) + 1 / pow(sing(1), 3));
  new_s2 = sqrt(1 + 1 / sing(2) + 1 / pow(sing(2), 2) + 1 / pow(sing(2), 3));
  new_sing << new_s0, new_s1, new_s2;
  Matrix3d Wi;
  Wi = ui * new_sing.asDiagonal() * ui.transpose();
  MatrixXd WW, PP;
  WW = Wi.transpose() * Wi;
  PP = P * P.transpose();

  for (uint32_t i0 = 0; i0 < 12; ++i0) {
    for (uint32_t i1 = 0; i1 < 12; ++i1) {
      tmpA(i0, i1) = WW(i0 / 4, i1 / 4) * PP(i0 % 4, i1 % 4);
    }
  }
  // matlab
  MatrixXd B = (P * (Wi * (J - ri)).transpose()) * Wi;
  tmpb << B(0, 0), B(1, 0), B(2, 0), B(3, 0), B(0, 1), B(1, 1), B(2, 1),
      B(3, 1), B(0, 2), B(1, 2), B(2, 2), B(3, 2);
}

static IGL_INLINE double perform_iteration(SCAFData &s) {
  using namespace Eigen;
  Eigen::MatrixXd V_out = s.w_V;
  MatrixXd A;
  // SparseMatrix<double> A;
  VectorXd b;
  A.resize(s.w_V.rows() * 3, s.w_V.rows() * 3);
  A.setZero();
  b.resize(s.w_V.rows() * 3);
  b.setZero();
  uint32_t _disp = s.w_V.rows();

  Matrix3d Q;
  Vector3d dp1, dp2, dp3;
  for (uint32_t i = 0; i < s.mf_num; ++i) {
    MatrixXd tmpA;
    VectorXd tmpb;
    scaf_build_element_linear_system(
        tmpA, tmpb, s.m_GradRef.block(i * 3, 0, 3, 3), V_out.row(s.m_T(i, 0)),
        V_out.row(s.m_T(i, 1)), V_out.row(s.m_T(i, 2)), V_out.row(s.m_T(i, 3)));
    for (uint32_t i0 = 0; i0 < 12; ++i0) {
      for (uint32_t i1 = 0; i1 < 12; ++i1) {
        A(s.m_T(i, i0 % 4) + (i0 / 4) * _disp,
          s.m_T(i, i1 % 4) + (i1 / 4) * _disp) += tmpA(i0, i1) * s.m_M(i);
      }
      b(s.m_T(i, i0 % 4) + (i0 / 4) * _disp) += tmpb(i0) * s.m_M(i);
    }
  }

  for (uint32_t i = 0; i < s.sf_num; ++i) {
    MatrixXd tmpA;
    VectorXd tmpb;
    scaf_build_element_linear_system(
        tmpA, tmpb, s.s_Grad.block(i * 3, 0, 3, 3), V_out.row(s.s_T(i, 0)),
        V_out.row(s.s_T(i, 1)), V_out.row(s.s_T(i, 2)), V_out.row(s.s_T(i, 3)));
    for (uint32_t i0 = 0; i0 < 12; ++i0) {
      for (uint32_t i1 = 0; i1 < 12; ++i1) {
        A(s.s_T(i, i0 % 4) + (i0 / 4) * _disp,
          s.s_T(i, i1 % 4) + (i1 / 4) * _disp) += tmpA(i0, i1) * s.s_M(i);
        if (A(s.s_T(i, i0 % 4) + (i0 / 4) * _disp,
              s.s_T(i, i1 % 4) + (i1 / 4) * _disp) > 1e20) {

          std::cout << A(s.s_T(i, i0 % 4) + (i0 / 4) * _disp,
                         s.s_T(i, i1 % 4) + (i1 / 4) * _disp)
                    << std::endl;
          int asdfhahskjfdl = 0;
        }
        if (tmpA(i0, i1) > 1e20) {
          int alksdjh = 0;
        }
        if (s.s_M(i) > 1e20) {
          int jghasdfadsfjhkg = 0;
        }
      }
      b(s.s_T(i, i0 % 4) + (i0 / 4) * _disp) += tmpb(i0) * s.s_M(i);
    }
  }

  for (auto x : s.soft_cons) {
    int _vn = x.first;
    int _disp = s.w_V.rows();
    Vector3d _vcoord = x.second;
    A(_vn, _vn) += s.soft_const_p;
    b(_vn) += s.soft_const_p * _vcoord(0);
    A(_vn + _disp, _vn + _disp) += s.soft_const_p;
    b(_vn + _disp) += s.soft_const_p * _vcoord(1);
    A(_vn + _disp * 2, _vn + _disp * 2) += s.soft_const_p;
    b(_vn + _disp * 2) += s.soft_const_p * _vcoord(2);
  }

  SimplicialLDLT<SparseMatrix<double>> solver;

  VectorXd d;
  d.resize(b.size());
  d = solver.compute(A.sparseView()).solve(b);
  std::ofstream fout;
  fout.open("A.txt", std::ios::app);
  printEigenMatrix(A, fout);
  printEigenVector(b, fout);
  fout.close();

  MatrixXd disp = s.w_V;
  disp.col(0) = d.block(0, 0, V_out.rows(), 1);
  disp.col(1) = d.block(V_out.rows(), 0, V_out.rows(), 1);
  disp.col(2) = d.block(V_out.rows() * 2, 0, V_out.rows(), 1);
  V_out = s.w_V + disp;
  auto whole_E = [&s](Eigen::MatrixXd &uv) {
    return compute_energy(s, uv, true);
  };

  Eigen::MatrixXi w_T;
  if (s.m_T.cols() == s.s_T.cols())
    igl::cat(1, s.m_T, s.s_T, w_T);
  else
    w_T = s.s_T;
  return igl::flip_avoiding_line_search(w_T, s.w_V, V_out, whole_E, -1) /
         s.mesh_measure;
}

IGL_INLINE void
scaf_precompute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                const Eigen::MatrixXd &V_init, const Eigen::MatrixXi &surface,
                igl::my_scaf::SCAFData &data,
                igl::MappingEnergyType slim_energy, Eigen::VectorXi &b,
                Eigen::MatrixXd &bc, double soft_p) {
  igl::my_scaf::add_new_patch(data, V, F, surface, V_init);
  data.soft_const_p = soft_p;
  for (int i = 0; i < b.rows(); i++)
    data.soft_cons[b(i)] = bc.row(i);
  data.slim_energy = slim_energy;

  // auto &s = data;

  if (!data.has_pre_calc) {
    //
    if (data.m_use_standard) {
      igl::my_scaf::compute_my_mesh_grad_inv_matrix_standard(
          data.m_Vref, data.m_T, data.m_GradRef);
    } else {
      igl::my_scaf::compute_my_mesh_grad_inv_matrix(data.m_Vref, data.m_T,
                                                    data.m_GradRef);
    }
    data.has_pre_calc = true;
  }
}

IGL_INLINE Eigen::MatrixXd scaf_solve(igl::my_scaf::SCAFData &s, int iter_num) {
  using namespace std;
  using namespace Eigen;
  s.energy = igl::my_scaf::compute_energy(s, s.w_V, false) / s.mesh_measure;

  for (int it = 0; it < iter_num; it++) {
    s.total_energy =
        igl::my_scaf::compute_energy(s, s.w_V, true) / s.mesh_measure;
    // s.rect_frame_V = Eigen::MatrixXd();
    igl::my_scaf::mesh_improve(s);

    double new_weight = s.mesh_measure * s.energy / (s.sf_num * 100);
    s.scaffold_factor = new_weight * 1e-4;
    igl::my_scaf::update_scaffold(s);

    s.total_energy = igl::my_scaf::perform_iteration(s);

    s.energy = igl::my_scaf::compute_energy(s, s.w_V, false) / s.mesh_measure;
  }

  return s.w_V.topRows(s.mv_num);
}
} // namespace my_scaf
} // namespace igl

#ifdef IGL_STATIC_LIBRARY
#endif

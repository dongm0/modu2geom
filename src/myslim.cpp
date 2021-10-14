

#ifndef MYSLIM_IMPL
#define MYSLIM_IMPL
#include "myslim.h"

#include <igl/doublearea.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/mapping_energy_with_jacobians.h>
#include <igl/polar_svd.h>
#include <igl/slice_into.h>
#include <igl/volume.h>

#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

//#include <igl/sparse_cached.h>

namespace {
// declare
void buildSLIMSystem(SLIMData &s);
void buildNewtonSystem(SLIMData &s);
// void preCalcSoftConstraints(SLIMData &s);
void preCalcRefJacobian(SLIMData &s);
double calEnergy(SLIMData &s);
void solveLinearSystem(SLIMData &s, Eigen::MatrixXd &dest);
void calJacobian(SLIMData &s);
void addSoftConstraint(SLIMData &s, Eigen::SparseMatrix<double> &L);
/*
Eigen::MatrixXd calc_line_constraint_coeff(double n1, double n2, double n3,
                                           double x0, double y0, double z0);
Eigen::MatrixXd calc_plain_constraint_coeff(double n1, double n2, double n3,
                                            double x0, double y0, double z0);
                                            */
std::pair<Eigen::MatrixXd, Eigen::VectorXd>
calElementSystem(const Eigen::Matrix<double, 1, Eigen::Dynamic> &J,
                 const Eigen::MatrixXd &P, igl::MappingEnergyType energy_type,
                 double exp_factor);
/*
void preCalcSoftConstraints(SLIMData &data) {
  data.line_coeff.resize(data.b_line.size() * 3, 4);
  for (int i = 0; i < data.b_line.size(); ++i) {
    Eigen::Vector3d X = data.bc_line.row(i).leftCols(3);
    Eigen::Vector3d n = data.bc_line.row(i).rightCols(3);
    data.line_coeff.block(i * 3, 0, 3, 4) =
        calc_line_constraint_coeff(n(0), n(1), n(2), X(0), X(1), X(2));
  }
  data.plain_coeff.resize(data.b_plain.size() * 3, 4);
  for (int i = 0; i < data.b_plain.size(); ++i) {
    Eigen::Vector3d X = data.bc_plain.row(i).leftCols(3);
    Eigen::Vector3d n = data.bc_plain.row(i).rightCols(3);
    data.plain_coeff.block(i * 3, 0, 3, 4) =
        calc_plain_constraint_coeff(n(0), n(1), n(2), X(0), X(1), X(2));
  }
}
Eigen::MatrixXd calc_line_constraint_coeff(double n1, double n2, double n3,
                                           double x0, double y0, double z0) {
  auto coeffx = [=](double n1, double n2, double n3) {
    return pow(n1, 4) + pow(n1, 2) * pow(n2, 2) + pow(n1, 2) * pow(n3, 2) -
           2 * pow(n1, 2) - 1;
  };
  auto coeffy = [=](double n1, double n2, double n3) {
    return n1 * n2 * (pow(n1, 2) + pow(n2, 2) + pow(n3, 2) - 2);
  };
  auto coeffz = [=](double n1, double n2, double n3) {
    return n1 * n3 * (pow(n1, 2) + pow(n2, 2) + pow(n3, 2) - 2);
  };
  auto coeffrhs = [=](double n1, double n2, double n3, double x0, double y0,
                      double z0) {
    double d = n1 * x0 + n2 * y0 + n3 * z0;
    return (pow(n1, 2) - 1) * (x0 - n1 * d) + n1 * n2 * (y0 - n2 * d) +
           n1 * n3 * (z0 - n3 * d);
  };
  Eigen::MatrixXd res;
  res.resize(3, 4);
  res.row(0) << coeffx(n1, n2, n3), coeffy(n1, n2, n3), coeffz(n1, n2, n3),
      coeffrhs(n1, n2, n3, x0, y0, z0);
  res.row(1) << coeffx(n2, n3, n1), coeffy(n2, n3, n1), coeffz(n2, n3, n1),
      coeffrhs(n2, n3, n1, y0, z0, x0);
  res.row(2) << coeffx(n3, n1, n2), coeffy(n3, n1, n2), coeffz(n3, n1, n2),
      coeffrhs(n3, n1, n2, z0, x0, y0);
  return res;
}
Eigen::MatrixXd calc_plain_constraint_coeff(double n1, double n2, double n3,
                                            double x0, double y0, double z0) {
  auto coeffx = [=](double n1, double n2, double n3) {
    return pow(n1, 2) * (pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
  };
  auto coeffy = [=](double n1, double n2, double n3) {
    return n1 * n2 * (pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
  };
  auto coeffz = [=](double n1, double n2, double n3) {
    return n1 * n3 * (pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
  };
  auto coeffrhs = [=](double n1, double n2, double n3, double x0, double y0,
                      double z0) {
    double d = n1 * x0 + n2 * y0 + n3 * z0;
    return -d * n1 * (pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
  };
  Eigen::MatrixXd res;
  res.resize(3, 4);
  res.row(0) << coeffx(n1, n2, n3), coeffy(n1, n2, n3), coeffz(n1, n2, n3),
      coeffrhs(n1, n2, n3, x0, y0, z0);
  res.row(1) << coeffx(n2, n3, n1), coeffy(n2, n3, n1), coeffz(n2, n3, n1),
      coeffrhs(n2, n3, n1, y0, z0, x0);
  res.row(2) << coeffx(n3, n1, n2), coeffy(n3, n1, n2), coeffz(n3, n1, n2),
      coeffrhs(n3, n1, n2, z0, x0, y0);
  return res;
}
*/

void preCalcRefJacobian(SLIMData &data) {
  if (data.use_ideal) {
    if (data.dim == 2) {
      data.P_ref.resize(2, 2);
      Eigen::Matrix2d localRefJ;
      auto v1 =
          (data.ideal_element.row(1) - data.ideal_element.row(0)).leftCols(2);
      auto v2 =
          (data.ideal_element.row(2) - data.ideal_element.row(0)).leftCols(2);
      localRefJ.col(0) = v1, localRefJ.col(1) = v2;
      data.P_ref = localRefJ.inverse();
    } else if (data.dim == 3) {
      data.P_ref.resize(3, 3);
      Eigen::Matrix3d localRefJ;
      auto v1 =
          (data.ideal_element.row(1) - data.ideal_element.row(0)).leftCols(3);
      auto v2 =
          (data.ideal_element.row(2) - data.ideal_element.row(0)).leftCols(3);
      auto v3 =
          (data.ideal_element.row(3) - data.ideal_element.row(0)).leftCols(3);
      localRefJ.col(0) = v1, localRefJ.col(1) = v2, localRefJ.col(2) = v3;
      data.P_ref = localRefJ.inverse();
    } else {
      assert(false);
    }
  } else {
    if (data.dim == 2) {
      data.P_ref.resize(data.t_num * 2, 2);
      for (int i = 0; i < data.t_num; ++i) {
        Eigen::Vector3d _x =
            data.V_ref.row(data.T(i, 1)) - data.V_ref.row(data.T(i, 0));
        Eigen::Vector3d _tmp =
            data.V_ref.row(data.T(i, 2)) - data.V_ref.row(data.T(i, 0));
        double x1 = _x.norm();
        _x /= x1;
        Eigen::Vector3d _n = _x.cross(_tmp);
        _n.normalize();
        Eigen::Vector3d _y = _n.cross(_x);
        double x2 = _tmp.dot(_x);
        double y2 = _tmp.dot(_y);
        Eigen::Matrix2d localRefJ;
        localRefJ << 1.0 / x1, -x2 / (x1 * y2), 0, 1.0 / y2;
        data.P_ref.block(i * 2, 0, 2, 2) = localRefJ;
      }
    } else if (data.dim == 3) {
      data.P_ref.resize(data.t_num * 3, 3);
      for (int i = 0; i < data.t_num; ++i) {
        Eigen::Vector3d _v1 =
            data.V_ref.row(data.T(i, 1)) - data.V_ref.row(data.T(i, 0));
        Eigen::Vector3d _v2 =
            data.V_ref.row(data.T(i, 2)) - data.V_ref.row(data.T(i, 0));
        Eigen::Vector3d _v3 =
            data.V_ref.row(data.T(i, 3)) - data.V_ref.row(data.T(i, 0));
        Eigen::Matrix3d localRefJ;
        localRefJ.col(0) = _v1;
        localRefJ.col(1) = _v2;
        localRefJ.col(3) = _v3;
        data.P_ref.block(i * 3, 0, 3, 3) = localRefJ.inverse();
      }
    } else {
      assert(false);
    }
  }
}
void calJacobian(SLIMData &data) {
  if (data.J.rows() == 0) {
    data.J.resize(data.t_num, data.dim * data.dim);
  }
  if (data.dim == 2) {
    for (int i = 0; i < data.t_num; ++i) {
      Eigen::Vector2d _v1 = data.V.row(data.T(i, 1)) - data.V.row(data.T(i, 0));
      Eigen::Vector2d _v2 = data.V.row(data.T(i, 2)) - data.V.row(data.T(i, 0));
      Eigen::Matrix2d localJ;
      localJ.row(0) = _v1;
      localJ.row(1) = _v2;
      localJ.transposeInPlace();
      if (data.use_ideal) {
        auto _tmp = localJ * data.P_ref;
        data.J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(1, 0), _tmp(1, 1);
      } else {
        auto _tmp =
            localJ * data.P_ref.block(i * data.dim, 0, data.dim, data.dim);
        data.J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(1, 0), _tmp(1, 1);
      }
    }
  } else if (data.dim == 3) {
    for (int i = 0; i < data.t_num; ++i) {
      Eigen::Vector3d _v1 = data.V.row(data.T(i, 1)) - data.V.row(data.T(i, 0));
      Eigen::Vector3d _v2 = data.V.row(data.T(i, 2)) - data.V.row(data.T(i, 0));
      Eigen::Vector3d _v3 = data.V.row(data.T(i, 3)) - data.V.row(data.T(i, 0));
      Eigen::Matrix3d localJ;
      localJ.row(0) = _v1;
      localJ.row(1) = _v2;
      localJ.row(2) = _v3;
      localJ.transposeInPlace();
      if (data.use_ideal) {
        auto _tmp = localJ * data.P_ref;
        data.J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(0, 2), _tmp(1, 0),
            _tmp(1, 1), _tmp(1, 2), _tmp(2, 0), _tmp(2, 1), _tmp(2, 2);
      } else {
        auto _tmp =
            localJ * data.P_ref.block(i * data.dim, 0, data.dim, data.dim);
        data.J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(0, 2), _tmp(1, 0),
            _tmp(1, 1), _tmp(1, 2), _tmp(2, 0), _tmp(2, 1), _tmp(2, 2);
      }
    }
  } else {
    assert(false);
  }
}

double calEnergy(SLIMData &data) {
  double e = 0;
  e += igl::mapping_energy_with_jacobians(data.J, data.M, data.slim_energy,
                                          data.exp_factor);
  for (int i = 0; i < data.b.size(); i++) {
    Eigen::Vector3d _d = data.bc.row(i) - data.V.row(data.b(i));
    for (int j = 0; j < 3; ++j) {
      if (data.fixed(i, j) == 1) {
        _d(j) = 0;
      }
    }
    e += data.soft_constraint_factor * _d.squaredNorm();
  }
  /*
  for (int i = 0; i < data.b_line.size(); ++i) {
    auto x = data.V.row(data.b_line(i));
    auto x0 = data.bc_line.row(i).leftCols(3);
    auto n = data.bc_line.row(i).rightCols(3);
    e += data.soft_constraint_factor *
         (x0 - x - n * ((x0 - x).dot(n))).squaredNorm();
  }
  for (int i = 0; i < data.b_plain.size(); ++i) {
    auto x = data.V.row(data.b_plain(i));
    auto x0 = data.bc_plain.row(i).leftCols(3);
    auto n = data.bc_plain.row(i).rightCols(3);
    e += data.soft_constraint_factor * (n * ((x0 - x).dot(n))).squaredNorm();
  }
  */
  return e;
}

void buildSLIMSystem(SLIMData &data) {
  data.IJV.clear();
  if (data.rhs.size() == 0) {
    data.rhs.resize(data.v_num * data.dim);
  }
  data.rhs.setZero();
  assert(data.dim == 2 || data.dim == 3);
  for (int i = 0; i < data.t_num; ++i) {
    Eigen::MatrixXd pref;
    if (data.use_ideal) {
      pref = data.P_ref;
    } else {
      pref = data.P_ref.block(i * data.dim, 0, data.dim, data.dim);
    }
    auto [localA, localb] = calElementSystem(data.J.row(i), pref,
                                             data.slim_energy, data.exp_factor);
    for (int j = 0; j < (data.dim + 1) * data.dim; ++j) {
      int _row =
          data.T(i, j % (data.dim + 1)) + (j / (data.dim + 1)) * data.v_num;
      for (int k = 0; k < (data.dim + 1) * data.dim; ++k) {
        int _col =
            data.T(i, k % (data.dim + 1)) + (k / (data.dim + 1)) * data.v_num;
        data.IJV.push_back(
            Eigen::Triplet<double>(_row, _col, localA(j, k) * data.M(i)));
      }
      data.rhs(_row) -= localb(j) * data.M(i);
    }
  }
}

std::pair<Eigen::MatrixXd, Eigen::VectorXd>
calElementSystem(const Eigen::Matrix<double, 1, Eigen::Dynamic> &J,
                 const Eigen::MatrixXd &P, igl::MappingEnergyType energy_type,
                 double exp_factor) {
  const int dim = P.rows();
  const double eps = 1e-8;
  const double sqrt_2 = sqrt(2);
  if (dim == 2) {
    Eigen::Matrix2d ji, ri, ti, ui, vi;
    Eigen::Vector2d sing, closest_sing_vec, m_sing_new;
    ji(0, 0) = J(0);
    ji(0, 1) = J(1);
    ji(1, 0) = J(2);
    ji(1, 1) = J(3);

    igl::polar_svd(ji, ri, ti, ui, sing, vi);
    double s1 = sing(0);
    double s2 = sing(1);
    switch (energy_type) {
    case igl::MappingEnergyType::ARAP: {
      m_sing_new << 1, 1;
      break;
    }
    case igl::MappingEnergyType::SYMMETRIC_DIRICHLET: {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
      break;
    }
    case igl::MappingEnergyType::LOG_ARAP: {
      double s1_g = 2 * (log(s1) / s1);
      double s2_g = 2 * (log(s2) / s2);
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
      break;
    }
    case igl::MappingEnergyType::CONFORMAL: {
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
    case igl::MappingEnergyType::EXP_CONFORMAL: {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));

      double geo_avg = sqrt(s1 * s2);
      double s1_min = geo_avg;
      double s2_min = geo_avg;

      double in_exp = exp_factor * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
      double exp_thing = exp(in_exp);

      s1_g *= exp_thing * exp_factor;
      s2_g *= exp_thing * exp_factor;

      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
      break;
    }
    case igl::MappingEnergyType::EXP_SYMMETRIC_DIRICHLET: {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));

      double in_exp =
          exp_factor * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2));
      double exp_thing = exp(in_exp);

      s1_g *= exp_thing * exp_factor;
      s2_g *= exp_thing * exp_factor;

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
    Eigen::Matrix2d Wi = ui * m_sing_new.asDiagonal() * ui.transpose();
    Eigen::Matrix<double, 3, 2> P1;
    P1.bottomRows(2) = P;
    P1.row(0) = Eigen::Vector2d(-P.col(0).sum(), -P.col(1).sum());
    Eigen::MatrixXd WW, PP;
    WW = Wi.transpose() * Wi;
    PP = P1 * P1.transpose();
    Eigen::MatrixXd resA;
    resA.resize((dim + 1) * dim, (dim + 1) * dim);
    Eigen::VectorXd resb;
    resb.resize((dim + 1) * dim);

    for (uint32_t i0 = 0; i0 < 6; ++i0) {
      for (uint32_t i1 = 0; i1 < 6; ++i1) {
        resA(i0, i1) = WW(i0 / 3, i1 / 3) * PP(i0 % 3, i1 % 3);
      }
    }
    // matlab
    Eigen::MatrixXd B = (P1 * (Wi * (ji - ri)).transpose()) * Wi;
    resb << B(0, 0), B(1, 0), B(2, 0), B(0, 1), B(1, 1), B(2, 1);
    return {resA, resb};
  } else if (dim == 3) {
    Eigen::Matrix3d ji, ri, ti, ui, vi;
    Eigen::Vector3d sing, closest_sing_vec, m_sing_new;
    // double s1, s2, s3;
    ji(0, 0) = J(0), ji(0, 1) = J(1), ji(0, 2) = J(2);
    ji(1, 0) = J(3), ji(1, 1) = J(4), ji(1, 2) = J(5);
    ji(2, 0) = J(6), ji(2, 1) = J(7), ji(2, 2) = J(8);
    igl::polar_svd(ji, ri, ti, ui, sing, vi);

    double s1 = sing(0);
    double s2 = sing(1);
    double s3 = sing(2);

    // 1) Update Weights
    switch (energy_type) {
    case igl::MappingEnergyType::ARAP: {
      m_sing_new << 1, 1, 1;
      break;
    }
    case igl::MappingEnergyType::LOG_ARAP: {
      double s1_g = 2 * (log(s1) / s1);
      double s2_g = 2 * (log(s2) / s2);
      double s3_g = 2 * (log(s3) / s3);
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
          sqrt(s3_g / (2 * (s3 - 1)));
      break;
    }
    case igl::MappingEnergyType::SYMMETRIC_DIRICHLET: {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));
      double s3_g = 2 * (s3 - pow(s3, -3));
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
          sqrt(s3_g / (2 * (s3 - 1)));
      break;
    }
    case igl::MappingEnergyType::EXP_SYMMETRIC_DIRICHLET: {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));
      double s3_g = 2 * (s3 - pow(s3, -3));
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
          sqrt(s3_g / (2 * (s3 - 1)));

      double in_exp = exp_factor * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) +
                                    pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
      double exp_thing = exp(in_exp);

      s1_g *= exp_thing * exp_factor;
      s2_g *= exp_thing * exp_factor;
      s3_g *= exp_thing * exp_factor;

      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1))),
          sqrt(s3_g / (2 * (s3 - 1)));

      break;
    }
    case igl::MappingEnergyType::CONFORMAL: {
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
    case igl::MappingEnergyType::EXP_CONFORMAL: {
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

      double in_exp = exp_factor * ((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) /
                                    (3 * pow((s1 * s2 * s3), 2. / 3)));
      ;
      double exp_thing = exp(in_exp);

      double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
      double s1_min = closest_s;
      double s2_min = closest_s;
      double s3_min = closest_s;

      s1_g *= exp_thing * exp_factor;
      s2_g *= exp_thing * exp_factor;
      s3_g *= exp_thing * exp_factor;

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
    // RMat3 mat_W;
    Eigen::Matrix3d Wi = ui * m_sing_new.asDiagonal() * ui.transpose();
    Eigen::Matrix<double, 4, 3> P1;
    P1.bottomRows(3) = P;
    P1.row(0) =
        Eigen::Vector3d(-P.col(0).sum(), -P.col(1).sum(), -P.col(2).sum());
    Eigen::MatrixXd WW, PP;
    WW = Wi.transpose() * Wi;
    PP = P1 * P1.transpose();
    Eigen::MatrixXd resA;
    resA.resize((dim + 1) * dim, (dim + 1) * dim);
    Eigen::VectorXd resb;
    resb.resize((dim + 1) * dim);

    for (uint32_t i0 = 0; i0 < 12; ++i0) {
      for (uint32_t i1 = 0; i1 < 12; ++i1) {
        resA(i0, i1) = WW(i0 / 4, i1 / 4) * PP(i0 % 4, i1 % 4);
      }
    }
    // matlab
    Eigen::MatrixXd B = (P1 * (Wi * (ji - ri)).transpose()) * Wi;
    // Eigen::MatrixXd B = (P1 * (Wi * (ji - ri)).transpose()) * Wi;
    resb << B(0, 0), B(1, 0), B(2, 0), B(3, 0), B(0, 1), B(1, 1), B(2, 1),
        B(3, 1), B(0, 2), B(1, 2), B(2, 2), B(3, 2);
    return {resA, resb};
  } else {
    assert(false);
  }
  return {Eigen::MatrixXd(), Eigen::VectorXd()};
}
void solveLinearSystem(SLIMData &s, Eigen::MatrixXd &dest) {
  Eigen::SparseMatrix<double> L(s.v_num * s.dim, s.v_num * s.dim);
  L.setFromTriplets(s.IJV.begin(), s.IJV.end());
  Eigen::VectorXd uv_flat(s.dim * s.v_num);
  for (int i = 0; i < s.dim; ++i) {
    uv_flat.block(i * s.v_num, 0, s.v_num, 1) = s.V.block(0, i, s.v_num, 1);
  }

  Eigen::SparseMatrix<double> id_m(s.v_num * s.dim, s.v_num * s.dim);
  id_m.setIdentity();
  s.rhs += L * uv_flat;
  addSoftConstraint(s, L);
  s.rhs += s.proximal_p * uv_flat;
  L = L + id_m * s.proximal_p;

  Eigen::VectorXd Uc;
  if (s.dim == 2) {
    using namespace Eigen;
    SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Uc = solver.compute(L).solve(s.rhs);
  } else { // seems like CG performs much worse for 2D and way better for 3D
    using namespace Eigen;
    Eigen::VectorXd guess;
    guess = std::move(uv_flat);
    BiCGSTAB<SparseMatrix<double>> solver;
    solver.setTolerance(1e-8);
    solver.compute(L);
    Uc = solver.solveWithGuess(s.rhs, guess);
    // ConjugateGradient<Eigen::SparseMatrix<double>, Lower | Upper> cg;
    // cg.setTolerance(1e-8);
    // cg.compute(L);
    // Uc = cg.solveWithGuess(s.rhs, guess);
  }

  for (int i = 0; i < s.dim; i++)
    dest.col(i) = Uc.block(i * s.v_num, 0, s.v_num, 1);
  // dest += s.V;
}

void addSoftConstraint(SLIMData &s, Eigen::SparseMatrix<double> &L) {
  int v_n = s.v_num;
  for (int i = 0; i < s.b.rows(); ++i) {
    for (int d = 0; d < s.dim; ++d) {
      int v_idx = s.b(i);
      if (s.fixed(i, d) == 1) {
        continue;
      }
      s.rhs(d * v_n + v_idx) += s.soft_constraint_factor * s.bc(i, d); // rhs
      L.coeffRef(d * v_n + v_idx, d * v_n + v_idx) +=
          s.soft_constraint_factor; // diagonal of matrix
    }
  }

  // custom constraints
  /*
  for (int i = 0; i < s.b_line.size(); i++) {
    for (int d = 0; d < s.dim; d++) {
      int v_idx = s.b_line(i);
      L.coeffRef(d * v_n + v_idx, v_idx) +=
          s.soft_constraint_factor * s.line_coeff(i * 3 + d, 0);
      L.coeffRef(d * v_n + v_idx, v_idx + v_n) +=
          s.soft_constraint_factor * s.line_coeff(i * 3 + d, 1);
      L.coeffRef(d * v_n + v_idx, v_idx + 2 * v_n) +=
          s.soft_constraint_factor * s.line_coeff(i * 3 + d, 2);
      s.rhs(d * v_n + v_idx) -=
          s.soft_constraint_factor * s.line_coeff(i * 3 + d, 3); // rhs
    }
  }
  for (int i = 0; i < s.b_plain.size(); i++) {
    for (int d = 0; d < s.dim; d++) {
      int v_idx = s.b_plain(i);
      L.coeffRef(d * v_n + v_idx, v_idx) +=
          s.soft_constraint_factor * s.plain_coeff(i * 3 + d, 0);
      L.coeffRef(d * v_n + v_idx, v_idx + v_n) +=
          s.soft_constraint_factor * s.plain_coeff(i * 3 + d, 1);
      L.coeffRef(d * v_n + v_idx, v_idx + 2 * v_n) +=
          s.soft_constraint_factor * s.plain_coeff(i * 3 + d, 2);
      s.rhs(d * v_n + v_idx) -=
          s.soft_constraint_factor * s.plain_coeff(i * 3 + d, 3); // rhs
    }
  }
  */
}

} // namespace

void myslim_precompute(SLIMData &data, Eigen::MatrixXd &&V, Eigen::MatrixXi &&T,
                       bool use_ideal, Eigen::MatrixXd &&V_ref,
                       Eigen::MatrixXd ideal_element,
                       igl::MappingEnergyType energy_type,
                       const Eigen::VectorXi &b, const Eigen::MatrixXd &bc,
                       const Eigen::MatrixXi &fixed, double soft_p) {
  data.V = V;
  data.T = T;
  data.use_ideal = use_ideal;
  if (data.use_ideal) {
    data.ideal_element = ideal_element;
  } else {
    data.V_ref = V_ref;
  }

  data.v_num = V.rows(), data.t_num = T.rows();
  data.slim_energy = energy_type;
  data.dim = T.cols() - 1;
  if (data.dim != 2 && data.dim != 3) {
    throw std::runtime_error("only support 2d or 3d data!");
  }

  data.soft_constraint_factor = soft_p;
  data.b = b, data.bc = bc, data.fixed = fixed;
  // data.b_line = b_line, data.bc_line = bc_line;
  // data.b_plain = b_plain, data.bc_plain = bc_plain;

  // preCalcSoftConstraints(data);

  igl::doublearea(V, T, data.M);
  data.M /= 2.;
  data.mesh_area = data.M.sum();
  data.exp_factor = 1.8;

  preCalcRefJacobian(data);
  calJacobian(data);
  data.energy = calEnergy(data) / data.mesh_area;

  data.has_pre_calc = true;
}

void myslim_solve(SLIMData &data, int iter_num) {
  for (int i = 0; i < iter_num; ++i) {
    Eigen::MatrixXd dest_res = data.V;
    calJacobian(data);
    buildSLIMSystem(data);
    solveLinearSystem(data, dest_res);

    double old_energy = data.energy;

    std::function<double(Eigen::MatrixXd &)> compute_energy =
        [&](Eigen::MatrixXd &cur) {
          Eigen::MatrixXd J;
          J.resize(data.t_num, data.dim * data.dim);
          if (data.dim == 2) {
            for (int i = 0; i < data.t_num; ++i) {
              Eigen::Vector2d _v1 =
                  cur.row(data.T(i, 1)) - cur.row(data.T(i, 0));
              Eigen::Vector2d _v2 =
                  cur.row(data.T(i, 2)) - cur.row(data.T(i, 0));
              Eigen::Matrix2d localJ;
              localJ.row(0) = _v1;
              localJ.row(1) = _v2;
              localJ.transposeInPlace();
              if (data.use_ideal) {
                auto _tmp = localJ * data.P_ref;
                J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(1, 0), _tmp(1, 1);
              } else {
                auto _tmp = localJ * data.P_ref.block(i * data.dim, 0, data.dim,
                                                      data.dim);
                J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(1, 0), _tmp(1, 1);
              }
            }
          } else if (data.dim == 3) {
            for (int i = 0; i < data.t_num; ++i) {
              Eigen::Vector3d _v1 =
                  cur.row(data.T(i, 1)) - cur.row(data.T(i, 0));
              Eigen::Vector3d _v2 =
                  cur.row(data.T(i, 2)) - cur.row(data.T(i, 0));
              Eigen::Vector3d _v3 =
                  cur.row(data.T(i, 3)) - cur.row(data.T(i, 0));
              Eigen::Matrix3d localJ;
              localJ.row(0) = _v1;
              localJ.row(1) = _v2;
              localJ.row(2) = _v3;
              localJ.transposeInPlace();
              if (data.use_ideal) {
                auto _tmp = localJ * data.P_ref;
                J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(0, 2), _tmp(1, 0),
                    _tmp(1, 1), _tmp(1, 2), _tmp(2, 0), _tmp(2, 1), _tmp(2, 2);
              } else {
                auto _tmp = localJ * data.P_ref.block(i * data.dim, 0, data.dim,
                                                      data.dim);
                J.row(i) << _tmp(0, 0), _tmp(0, 1), _tmp(0, 2), _tmp(1, 0),
                    _tmp(1, 1), _tmp(1, 2), _tmp(2, 0), _tmp(2, 1), _tmp(2, 2);
              }
            }
          } else {
            assert(false);
          }

          double e = 0;
          e += igl::mapping_energy_with_jacobians(J, data.M, data.slim_energy,
                                                  data.exp_factor);
          for (int i = 0; i < data.b.size(); i++) {
            Eigen::Vector3d _d = data.bc.row(i) - data.V.row(data.b(i));
            for (int j = 0; j < 3; ++j) {
              if (data.fixed(i, j) == 1) {
                _d(j) = 0;
              }
            }
            e += data.soft_constraint_factor * _d.squaredNorm();
          }

          return e;
        };

    data.energy =
        igl::flip_avoiding_line_search(data.T, data.V, dest_res, compute_energy,
                                       data.energy * data.mesh_area) /
        data.mesh_area;
    std::cout << "iter: " << i << " energy: " << data.energy << std::endl;

    if (fabs(old_energy - data.energy) / fabs(data.energy) < 1e-9) {
      std::cout << "already converge. stop iteration." << std::endl;
      break;
    }
  }
}

#endif
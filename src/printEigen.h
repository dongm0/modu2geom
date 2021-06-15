
#ifndef PRINTEIGEN_H
#define PRINTEIGEN_H
#include <Eigen/Eigen>
#include <iostream>

template <typename Matrix>
void printEigenMatrix(const Matrix &matrix, std::ostream &stream) {
  for (int i = 0; i < matrix.rows(); ++i) {
    for (int j = 0; j < matrix.cols(); ++j) {
      stream << matrix(i, j) << " ";
    }
    stream << std::endl;
  }
  stream << std::endl;
}

template <typename Vec>
void printEigenVector(const Vec &v, std::ostream &stream) {
  for (int i = 0; i < v.size(); ++i) {
    stream << v(i) << " ";
  }
  stream << std::endl;
  stream << std::endl;
}
#endif
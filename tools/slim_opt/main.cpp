#include <Eigen/Eigen>
#include <fstream>
#include <igl/slim.h>
#include <string>

using namespace Eigen;

const int hex2tet[32] = {0, 1, 3, 4, 1, 2, 0, 5, 2, 3, 1, 6, 3, 0, 2, 7,
                         4, 7, 5, 0, 7, 6, 4, 3, 6, 5, 7, 2, 5, 4, 6, 1};

std::pair<MatrixXi, MatrixXd> readVtk(std::string filename) {
  std::ifstream fin;
  fin.open(filename);
  char buf[512];
  std::string tmp;
  int i;
  double l;
  fin.getline(buf, 512);
  fin.getline(buf, 512);
  fin.getline(buf, 512);
  fin.getline(buf, 512);
  fin >> tmp >> i >> tmp;
  MatrixXd V(i, 3);
  for (int kk = 0; kk < i; ++kk) {
    double i1, i2, i3;
    fin >> i1 >> i2 >> i3;
    V(kk, 0) = i1;
    V(kk, 1) = i2;
    V(kk, 2) = i3;
  }
  fin >> tmp >> i >> tmp;
  MatrixXd T(i * 8, 4);
  for (int kk = 0; kk < i; ++kk) {
    int j[8];
    fin >> tmp;
    for (int kkk = 0; kkk < 8; ++kkk)
      fin >> j[kkk];
    for (int kkk = 0; kkk < 8; ++kkk) {
      T(kk * 8 + kkk, 0) = j[hex2tet[kkk * 4]];
      T(kk * 8 + kkk, 1) = j[hex2tet[kkk * 4 + 1]];
      T(kk * 8 + kkk, 2) = j[hex2tet[kkk * 4 + 2]];
      T(kk * 8 + kkk, 3) = j[hex2tet[kkk * 4 + 3]];
    }
  }
  fin.close();
  return {T, V};
}

int main(int argc, char **argv) {
  MatrixXd V;
  MatrixXi T;
  auto p = readVtk(std::string(argv[1]));
}
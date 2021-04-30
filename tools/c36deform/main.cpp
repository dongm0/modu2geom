#include <eigen3/Eigen/Eigen>
#include <fstream>
#include <igl/slim.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace Eigen;

const int hex2tet[32] = {0, 1, 3, 4, 1, 2, 0, 5, 2, 3, 1, 6, 3, 0, 2, 7,
                         4, 7, 5, 0, 7, 6, 4, 3, 6, 5, 7, 2, 5, 4, 6, 1};

std::vector<int> bt1 = {33, 47, 37, 49, 43, 50, 38, 44, 42,
                        34, 36, 46, 41, 48, 39, 45, 40, 35};
std::vector<std::vector<double>> bmap = {
    {1.9501, -1.4861, 0.4262},  {-0.2029, 0.9559, -0.0952},
    {3.0429, 0.0544, 3.1287},   {4.1954, 0.1059, -1.3891},
    {0.8574, -3.0265, -2.2763}, {-0.2951, -3.0781, 2.2415},
    {3.6191, 0.0802, 0.8698},   {2.5264, -1.4603, -1.8327},
    {0.2812, -3.0523, -0.0174}, {1.3739, -1.5118, 2.6851},
    {1.4200, 0.5051, 1.5168},   {1.9962, 0.5309, -0.7421},
    {0.3273, -1.0353, -1.1858}, {-0.2490, -1.0611, 1.0731},
    {2.3451, 0.3721, 0.5481},   {1.6166, -0.6549, -1.2535},
    {0.1198, -1.7162, -0.0434}, {0.8483, -0.6893, 1.75831}};

std::pair<MatrixXi, MatrixXd> readvtk(const std::string &filename) {
  std::ifstream fin;
  fin.open(filename);
  char tmpchar[256];
  for (int i = 0; i < 4; ++i) {
    fin.getline(tmpchar, 256);
  }
  int cnum, vnum;
  std::string tmpstr;
  fin >> tmpstr >> vnum;
  fin.getline(tmpchar, 256);
  MatrixXd V(vnum, 3);
  for (int i = 0; i < vnum; ++i) {
    fin >> V(i, 0) >> V(i, 1) >> V(i, 2);
  }
  fin >> tmpstr >> cnum;
  fin.getline(tmpchar, 256);
  MatrixXi T(cnum * 8, 4);
  for (int i = 0; i < cnum; ++i) {
    std::vector<int> ttt(8);
    int tt;
    fin >> tt;
    for (int j = 0; j < 8; ++j) {
      fin >> ttt[j];
    }
    for (int j = 0; j < 8; ++j) {
      for (int k = 0; k < 4; ++k) {
        T(i * 8 + j, k) = ttt[hex2tet[j * 4 + k]];
      }
    }
  }
  fin.close();
  return {T, V};
}

int main(int argc, char **argv) {
  std::map<int, std::vector<double>> bm;
  for (int i = 0; i < 18; ++i) {
    bm.insert({bt1[i], bmap[i]});
  }
  VectorXi b(18);
  MatrixXd bc(18, 3);
  int i = 0;
  for (auto &x : bm) {
    b(i) = x.first;
    bc(i, 0) = x.second[0], bc(i, 1) = x.second[1], bc(i, 2) = x.second[2];
    ++i;
  }
  auto [T, V] = readvtk("../msqtest/out_my.vtk");
  igl::SLIMData sData;
  auto V0 = V;
  /*
  for (int i = 0; i < 18; ++i) {
    V0(b(i), 0) = bc(i, 0);
    V0(b(i), 1) = bc(i, 1);
    V0(b(i), 2) = bc(i, 2);
  }
  */
  igl::slim_precompute(V, T, V0, sData,
                       igl::MappingEnergyType::SYMMETRIC_DIRICHLET, b, bc, 1e6);
  igl::slim_solve(sData, 20);
  for (int i = 0; i < V0.rows(); ++i) {
    // std::cout << i << std::endl;
    std::cout << sData.V_o(i, 0) << " " << sData.V_o(i, 1) << " "
              << sData.V_o(i, 2) << std::endl;
  }

  return 0;
}

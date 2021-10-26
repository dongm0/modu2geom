#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "ovmwrapper.h"
#include "utils.h"
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <algorithm>
#include <vector>

using OVM = OpenVolumeMesh::HexahedralMeshTopologyKernel;

struct Topology {
  size_t c_size = 0;
  size_t v_size = 0;
  std::vector<std::vector<unsigned char>> data;

  OVM ToOVM() const {
    OVM mesh;
    if (c_size == 0) {
      throw std::runtime_error("topology needs init!");
    }
    for (int i = 0; i < v_size; ++i) {
      mesh.add_vertex();
    }
    for (int i = 0; i < c_size; ++i) {
      std::vector<OpenVolumeMesh::VertexHandle> _tmp_v_vec;
      std::for_each(data[i].begin(), data[i].end(), [&](unsigned char idx) {
        _tmp_v_vec.push_back(OpenVolumeMesh::VertexHandle(idx));
      });
      std::swap(_tmp_v_vec[1], _tmp_v_vec[3]);
      mesh.add_cell(_tmp_v_vec);
    }
    return mesh;
  }
  void FromOVM(const OVM &ovm) {
    data.clear();
    data.resize(c_size);
    int c = 0;
    for (auto cit : ovm.cells()) {
      for (auto vit : ovm.cell_vertices(cit)) {
        data[c].push_back(vit.idx());
      }
      std::swap(data[c][1], data[c][3]);
      ++c;
    }
  }

  Topology(const std::string &filename) {
    std::ifstream fin;
    fin.open(filename);
    fin >> *this;
  }
  Topology(const OVM &ovm,
           const std::vector<OpenVolumeMesh::CellHandle> &toposeq)
      : c_size(ovm.n_cells()), v_size(ovm.n_vertices()) {
    data.resize(c_size);
    int c = 0;
    for (auto cit : toposeq) {
      for (auto vit : ovm.cell_vertices(cit)) {
        data[c].push_back(vit.idx());
      }
      std::swap(data[c][1], data[c][3]);
      ++c;
    }
  }

  friend std::istream &operator>>(std::istream &in, Topology &A);
  friend std::ostream &operator<<(std::ostream &out, const Topology &A);
};

std::istream &operator>>(std::istream &in, Topology &A) {
  in >> A.c_size;
  A.data.resize(A.c_size);
  for (int i = 0; i < A.c_size; ++i) {
    A.data[i].resize(8);
    for (int j = 0; j < 8; ++j) {
      size_t _tmp;
      in >> _tmp;
      A.data[i][j] = _tmp;
      A.v_size = std::max(A.v_size, _tmp);
    }
  }
  A.v_size += 1;
  return in;
}

std::ostream &operator<<(std::ostream &out, const Topology &A) {
  out << A.c_size << std::endl;
  for (int i = 0; i < A.c_size; ++i) {
    for (int j = 0; j < 8; ++j) {
      out << (int)A.data[i][j] << " ";
    }
    out << std::endl;
  }
  return out;
}

#endif
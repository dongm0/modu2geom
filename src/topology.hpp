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

  Topology() {}
  Topology(const std::string &filename) {
    std::ifstream fin;
    fin.open(filename);
    fin >> *this;
  }
  // Topology(const OVM &ovm,
  //         const std::vector<OpenVolumeMesh::CellHandle> &toposeq)
  //    : c_size(ovm.n_cells()), v_size(ovm.n_vertices()) {
  //  data.resize(c_size);
  //  int c = 0;
  //  for (auto cit : ovm.cells()) {
  //    auto xf_face = ovm.get_oriented_halfface(0, cit);
  //    auto [xf_begin, xf_end] = ovm.halfface_vertices(xf_face);
  //    std::vector<OpenVolumeMesh::VertexHandle> xf(xf_begin, xf_end);
  //    auto xb = opposite_vertex_in_cell(ovm, cit, xf_face, xf);
  //    std::for_each(xf.begin(), xf.end(), [&](OpenVolumeMesh::VertexHandle vh)
  //    {
  //      data[c].push_back(vh.idx());
  //    });
  //    std::for_each(xb.begin(), xb.end(), [&](OpenVolumeMesh::VertexHandle vh)
  //    {
  //      data[c].push_back(vh.idx());
  //    });
  //    ++c;
  //  }
  //}

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
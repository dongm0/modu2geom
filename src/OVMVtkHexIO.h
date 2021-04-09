#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenVolumeMesh;

template <typename HexMesh>
void OVMReadHexVtkStream(HexMesh &mesh, std::istream &stream,
                         bool _topologyCheck = true,
                         bool _computeBottomUpIncidences = true) {
  std::stringstream ss;
  std::string line;
  std::string stmp;

  typedef typename HexMesh::PointT Point;
  Point v(0.f, 0.f, 0.f);

  uint32_t vnum = 0, cnum = 0;

  mesh.clear(false);
  mesh.enable_bottom_up_incidences(false);

  std::getline(stream, line);
  std::getline(stream, line);
  std::getline(stream, line);
  if (line != "ASCII") {
    throw std::runtime_error("ascii only!");
  }
  std::getline(stream, line);
  std::getline(stream, line);
  ss.clear();
  ss.str(line);
  ss >> stmp >> vnum >> stmp;

  for (uint32_t i = 0; i < vnum; ++i) {
    std::getline(stream, line);
    ss.clear();
    ss.str(line);
    ss >> v[0] >> v[1] >> v[2];
    mesh.add_vertex(v);
  }

  std::getline(stream, line);
  ss.clear();
  ss.str(line);
  ss >> stmp >> cnum >> stmp;

  uint32_t _r[8];

  for (uint32_t i = 0; i < cnum; ++i) {
    std::getline(stream, line);
    ss.clear();
    ss.str(line);
    ss >> stmp >> _r[0] >> _r[1] >> _r[2] >> _r[3] >> _r[4] >> _r[5] >> _r[6] >>
        _r[7];
    std::vector<VertexHandle> _vhs = {VertexHandle(_r[0]), VertexHandle(_r[3]),
                                      VertexHandle(_r[2]), VertexHandle(_r[1]),
                                      VertexHandle(_r[4]), VertexHandle(_r[5]),
                                      VertexHandle(_r[6]), VertexHandle(_r[7])};
    mesh.add_cell(_vhs, _topologyCheck);
  };
}

template <typename HexMesh>
void OVMWriteHexMesh(HexMesh &mesh, std::ofstream &stream) {
  stream << "# vtk DataFile Version 3.0\ngenerate by dongmo\nASCII\nDATASET "
            "UNSTRUCTURED_GRID"
         << std::endl;
  stream << "POINTS " << mesh.n_vertices() << " double" << std::endl;
  for (auto &v : mesh.vertices()) {
    auto p = mesh.vertex(v);
    stream << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }
  stream << "CELLS " << mesh.n_cells() << " " << 9 * mesh.n_cells()
         << std::endl;
  std::vector<VertexHandle> _vtmp;
  for (auto &c : mesh.cells()) {
    stream << "8 ";
    _vtmp.clear();
    for (auto v : mesh.hex_vertices(c)) {
      _vtmp.push_back(v);
    }
    stream << _vtmp[0].idx() << " " << _vtmp[1].idx() << " " << _vtmp[2].idx()
           << " " << _vtmp[3].idx() << " " << _vtmp[4].idx() << " "
           << _vtmp[7].idx() << " " << _vtmp[6].idx() << " " << _vtmp[5].idx();
    stream << std::endl;
  }
  stream << "CELL_TYPES " << mesh.n_cells() << std::endl;
  for (uint32_t i = 0; i < mesh.n_cells(); ++i) {
    stream << "12" << std::endl;
  }
}

template <typename HexMesh>
void OVMWriteHexMesh(HexMesh &mesh, std::ofstream &stream,
                     std::vector<VertexHandle> &_tagged) {
  OVMWriteHexMesh(mesh, stream);
  stream << "POINT_DATA " << mesh.n_vertices()
         << "SCALARS fixed int\nLOOKUP_TABLE default" << std::endl;
  std::vector<uint8_t> _fixed(mesh.n_vertices(), 0);
  for (const auto x : _tagged) {
    _fixed.at(x.idx()) = 1;
  }
  for (const auto x : _fixed) {
    stream << x << std::endl;
  }
}
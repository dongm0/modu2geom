#include <OVMVtkHexIO.h>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <fstream>
#include <iostream>
#include <string>

using Mymesh = OpenVolumeMesh::GeometricHexahedralMeshV3d;

template <typename HexMesh>
static void OVMWriteHexMesh_p(HexMesh &mesh, std::ofstream &stream) {
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
    stream << _vtmp[0].idx() << " " << _vtmp[3].idx() << " " << _vtmp[2].idx()
           << " " << _vtmp[1].idx() << " " << _vtmp[4].idx() << " "
           << _vtmp[5].idx() << " " << _vtmp[6].idx() << " " << _vtmp[7].idx();
    stream << std::endl;
  }
  stream << "CELL_TYPES " << mesh.n_cells() << std::endl;
  for (uint32_t i = 0; i < mesh.n_cells(); ++i) {
    stream << "12" << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    std::cout << "invalid argument number." << std::endl;
    return 1;
  }
  Mymesh mesh;
  std::ifstream fin;
  fin.open(std::string(argv[1]));
  OVMReadHexVtkStream(mesh, fin);
  fin.close();
  std::ofstream fout;
  if (argc == 3) {
    fout.open(std::string(argv[2]));
  } else {
    fout.open(std::string("out_") + std::string(argv[1]));
  }
  OVMWriteHexMesh_p(mesh, fout);
  return 0;
}
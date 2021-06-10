#include "mymesh.h"
#include "utils.h"

std::string filepath("/home/dm/github/modu2geom/modudata/modu/");

std::string ifilename("con36-1.txt");
std::string ofilename("res.vtk");

int main(int argc, char **argv) {
#ifdef NDEBUG
  if (argc != 2) {
    std::cerr << "Usage: modu2geom <file>" << std::endl;
    return 1;
  }
#endif
  MyMesh mesh;
#ifdef NDEBUG
  mesh.ReadTopoFromFile(std::string(argv[1]));
#else
  mesh.ReadTopoFromFile(std::string(argv[1]));
  // mesh.ReadTopoFromFile(filepath + ifilename);
#endif
  mesh.GenerateOrder();
  mesh.checkTopo();
  for (int i = 0; i < mesh.GetTopoCnum(); ++i) {
    if (i >= 30) {
      int _tmpppp = 0;
    }
    mesh.GenerateOneCell(mesh.GetCurrentCellHandle());
    mesh.Optimize();
    // mesh.WriteGeomToVTKFile(std::to_string(i) + ".vtk");
    std::cout << i << std::endl;
  }
  mesh.WriteGeomToVTKFile(ofilename);
  return 0;
}
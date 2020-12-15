#include "utils.h"
#include "mymesh.h"

std::string ifilename("con171.txt");
std::string ofilename("");

int main() {
    MyMesh mesh;
    mesh.ReadTopoFromFile(ifilename);
    mesh.GenerateOrder();
    for (int i=0; i<mesh.GetTopoCnum(); ++i) {
        mesh.GenerateOneCell(mesh.GetCurrentCellHandle());
        mesh.WriteGeomToVTKFile(ofilename+"./src"+std::to_string(i)+".vtk");
    }
    mesh.WriteGeomToVTKFile(ofilename);
    return 0;
    
}
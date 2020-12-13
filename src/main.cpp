#include "utils.h"
#include "mymesh.h"

std::string ifilename("");
std::string ofilename("");

int main() {
    MyMesh mesh;
    mesh.ReadTopoFromFile(ifilename);
    mesh.GenerateOrder();
    for (int i=0; i<mesh.GetTopoCnum(); ++i) {
        mesh.GenerateOneCell(mesh.GetCurrentCellHandle());
    }
    mesh.WriteGeomToVTKFile(ofilename);
    return 0;
    
}
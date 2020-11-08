#include "utils.h"
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

class GeomMesh {
private:
    OpenVolumeMesh::GeometricHexahedralMeshV3d m_mesh;
    std::vector<OpenVolumeMesh::VertexHandle> m_vertices;
    
};
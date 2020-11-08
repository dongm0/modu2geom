/*topomesh的头文件，记录了拓扑结构的网格，用ovm volume mesh 存储，并计算生成顺序*/

#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include "utils.h"

class MyMesh {
public:

    //可以用>>操作符
    ErrorCode ReadTopoFromFile(const std::string &filename);
    ErrorCode GenerateOrder(std::vector<uint8_t> &_order);
    ErrorCode GenerateOneCell(const OpenVolumeMesh::CellHandle &_ch);

private:
    OpenVolumeMesh::GeometricHexahedralMeshV3d m_mesh;
    OpenVolumeMesh::TopologicHexahedralMesh m_topomesh;
    std::vector<OpenVolumeMesh::VertexHandle> m_vertices;
    std::vector<OpenVolumeMesh::VertexHandle> m_topo_vertices;
    std::vector<std::vector<uint8_t>> m_cells;
    std::map<OpenVolumeMesh::CellHandle, OpenVolumeMesh::CellHandle> m_m2tm_mapping;
    std::map<OpenVolumeMesh::CellHandle, OpenVolumeMesh::CellHandle> m_tm2m_mapping;
};
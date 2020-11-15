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

    ErrorCode Optimize();

private:

    //加入cell的一系列函数
    ErrorCode AddOneCellCase0(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    ErrorCode AddOneCellCase1(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    ErrorCode AddOneCellCase2(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    ErrorCode AddOneCellCase3(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    ErrorCode AddOneCellCase4(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    ErrorCode AddOneCellCase5(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    ErrorCode AddOneCellCase6(const OpenVolumeMesh::CellHandle &_ch, 
        const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
        const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec);

    //一系列操作函数，形式上简洁一些
    OpenVolumeMesh::CellHandle getGeomC(const OpenVolumeMesh::CellHandle &_ch) {return m_tm2m_mapping[_ch];}
    OpenVolumeMesh::CellHandle getTopoC(const OpenVolumeMesh::CellHandle &_ch) {return m_m2tm_mapping[_ch];}
    OpenVolumeMesh::VertexHandle getTopoV(const OpenVolumeMesh::VertexHandle &_vh) {return m_m2tm_v_mapping[_vh];}
    OpenVolumeMesh::VertexHandle getGeomV(const OpenVolumeMesh::VertexHandle &_vh) {return m_tm2m_v_mapping[_vh];}

    std::vector<OpenVolumeMesh::VertexHandle> opposite_vertex_in_cell(const OpenVolumeMesh::HexahedralMeshTopologyKernel &mesh, 
    const OpenVolumeMesh::CellHandle &cell, const OpenVolumeMesh::HalfFaceHandle &hf, 
    const std::vector<OpenVolumeMesh::VertexHandle> &v) {
        auto opposite_hf = mesh.opposite_halfface_handle_in_cell(hf, cell);
        auto range_t = mesh.halfface_vertices(opposite_hf);
        std::set<OpenVolumeMesh::VertexHandle> topset(range_t.first, range_t.second);
        std::vector<OpenVolumeMesh::VertexHandle> res;
        for (auto x : v) {
            for (auto vv_it=mesh.vv_iter(x); vv_it.valid(); ++vv_it) {
                if (topset.count(*vv_it)) {
                    res.push_back(*vv_it);
                    break;
                }
            }
        }
        return res;
    }

    OpenVolumeMesh::GeometricHexahedralMeshV3d m_mesh;
    OpenVolumeMesh::TopologicHexahedralMesh m_topomesh;
    std::vector<OpenVolumeMesh::VertexHandle> m_vertices;
    std::vector<OpenVolumeMesh::VertexHandle> m_topo_vertices;
    std::vector<std::vector<uint8_t>> m_cells;
    std::map<OpenVolumeMesh::CellHandle, OpenVolumeMesh::CellHandle> m_m2tm_mapping;
    std::map<OpenVolumeMesh::CellHandle, OpenVolumeMesh::CellHandle> m_tm2m_mapping;
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::VertexHandle> m_m2tm_v_mapping;
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::VertexHandle> m_tm2m_v_mapping;
};
#include "mymesh.h"

ErrorCode MyMesh::ReadTopoFromFile(const std::string &filename) {
    using namespace OpenVolumeMesh;
    
    std::ifstream fin(filename.c_str());

    uint8_t _cnum = 0;
    uint8_t _vnum = 0;
    fin >> _cnum;
    m_cells.assign(_cnum, std::vector<uint8_t>(8));
    for (uint8_t i=0; i<_cnum; ++i) {
        for (uint8_t j=0; j<8; ++j) {
            fin >> m_cells.at(i).at(j);
            _vnum = std::max(_vnum, m_cells.at(i).at(j));
        }
    }
    fin.close();

    //add vertices
    for (uint8_t i=0; i<_vnum; ++i) {
        m_vertices.push_back(m_mesh.add_vertex(Vec3d(0, 0, 0)));
        m_topo_vertices.push_back(m_topomesh.add_vertex());
    }

    //add topo cells
    for (uint8_t i=0; i<_cnum; ++i) {
        std::vector<VertexHandle> _cell;
        for (uint8_t j=0; j<8; ++j) {
            _cell.push_back(m_topo_vertices.at(m_cells.at(i).at(j)));
        }
        #ifdef OVM_TOPOLOGY_CHECK
            m_topomesh.add_cell(_cell, true);
        #else
            m_topomesh.add_cell(_cell, false);
        #endif
    }
    return ErrorCode::succeed;
}

ErrorCode MyMesh::GenerateOneCell(const OpenVolumeMesh::CellHandle &_ch) {
    using namespace OpenVolumeMesh;

    if (m_tm2m_mapping.find(_ch) != m_tm2m_mapping.end()) {
        auto _tmpch = m_tm2m_mapping[_ch];
        if (m_m2tm_mapping.find(_tmpch) != m_m2tm_mapping.end() && m_m2tm_mapping[_tmpch]==_ch)
            return ErrorCode::succeed;
        else
            return ErrorCode::failed;
    }
    // 6 cases
    int num_nbh = 0, casenum = -1;
    std::vector<CellHandle> _nbh;
    {
        m_topomesh.cell(_ch).halffaces();
        for (auto cc_iter=m_topomesh.cc_iter(_ch); cc_iter->is_valid(); ++cc_iter) {
            if (m_tm2m_mapping.find(_ch) == m_tm2m_mapping.end()) {
                num_nbh++;
                _nbh.push_back(*cc_iter);
            }
        }
        if (num_nbh == 3) {
            
        }
    }
}
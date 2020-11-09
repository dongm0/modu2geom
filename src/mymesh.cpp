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
        auto _topo_v = m_mesh.add_vertex(Vec3d(0, 0, 0));
        auto _v = m_topomesh.add_vertex();
        m_vertices.push_back(_topo_v);
        m_topo_vertices.push_back(_v);
        m_m2tm_v_mapping[_v] = _topo_v;
        m_tm2m_v_mapping[_topo_v] = _v;
    }

    //add topo cells
    for (uint8_t i=0; i<_cnum; ++i) {
        std::vector<VertexHandle> _cell;
        for (uint8_t j=0; j<8; ++j) {
            _cell.push_back(m_topo_vertices.at(m_cells.at(i).at(j)));
        }
        auto tmp = _cell[5];
        _cell[5] = _cell[7];
        _cell[7] = tmp;
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
    std::vector<HalfFaceHandle> _nbh_hf;
    {
        for (auto hf_handle : m_topomesh.cell(_ch).halffaces()) {
            auto opposite_hf = m_topomesh.opposite_halfface_handle(hf_handle);
            if (opposite_hf.is_valid() && opposite_hf.idx()!=hf_handle.idx()) {
                auto _nbch = m_topomesh.incident_cell(opposite_hf);
                if (m_tm2m_mapping.find(_ch) == m_tm2m_mapping.end()) {
                    num_nbh++;
                    _nbh.push_back(_nbch);
                    _nbh_hf.push_back(hf_handle);
                }
            }
        }
        for (auto cc_iter=m_topomesh.cc_iter(_ch); cc_iter->is_valid(); ++cc_iter) {
            if (m_tm2m_mapping.find(_ch) == m_tm2m_mapping.end()) {
                num_nbh++;
                _nbh.push_back(*cc_iter);
            }
        }
        if (num_nbh == 3) {
            if (m_topomesh.opposite_halfface_handle_in_cell(_nbh_hf[0], _ch) == _nbh_hf[1]
            or m_topomesh.opposite_halfface_handle_in_cell(_nbh_hf[0], _ch) == _nbh_hf[2]
            or m_topomesh.opposite_halfface_handle_in_cell(_nbh_hf[1], _ch) == _nbh_hf[2])
                casenum = 3;
            else
                casenum = 4;
        }
        else {
            if (num_nbh==0 or num_nbh==1 or num_nbh==2)
                casenum = num_nbh;
            else if (num_nbh==4 or num_nbh==5)
                casenum = num_nbh+1;
            else
                return ErrorCode::failed;
        }
    }

}

ErrorCode MyMesh::AddOneCellCase0(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
    auto bottomface = m_topomesh.xfront_halfface(_ch);
    auto topface = m_topomesh.xback_halfface(_ch);
    auto b_range = m_topomesh.halfface_vertices(bottomface);
    auto t_range = m_topomesh.halfface_vertices(topface);
    std::set<OpenVolumeMesh::VertexHandle> topset(t_range.first, t_range.second);
    std::vector<OpenVolumeMesh::VertexHandle> bottomvec(b_range.first, b_range.second);
    std::vector<OpenVolumeMesh::VertexHandle> topvec();

    for (int i=0; i<4; ++i) {
        auto bv = bottomvec[i];
        for (auto vv_it=m_topomesh.vv_iter(bv); vv_it.valid(); ++vv_it) {
            if (topset.find(*vv_it) != topset.end()) {
                topvec.push_back(*vv_it);
                break;
            }
        }
    }

    
    
}
ErrorCode MyMesh::AddOneCellCase1(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {

}
ErrorCode MyMesh::AddOneCellCase2(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {

}
ErrorCode MyMesh::AddOneCellCase3(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {

}
ErrorCode MyMesh::AddOneCellCase4(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {

}
ErrorCode MyMesh::AddOneCellCase5(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {

}
ErrorCode MyMesh::AddOneCellCase6(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {

}


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
                casenum = 4;
            else
                casenum = 3;
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
    #ifndef NDEBUG
        assert(_nbc_vec.size()==0 && _nbhf_vec.size()==0);
    #endif
    using namespace OpenVolumeMesh;
    auto bottomface = m_topomesh.zback_halfface(_ch);
    auto b_range = m_topomesh.halfface_vertices(bottomface);
    std::vector<VertexHandle> bottomvec(b_range.first, b_range.second);
    std::vector<VertexHandle> topvec = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], bottomvec);

    #ifndef NDEBUG
        assert(topvec.size()==4);
    #endif
    m_mesh.set_vertex(getGeomV(bottomvec[0]), Vec3d(0, 0, 0));
    m_mesh.set_vertex(getGeomV(bottomvec[1]), Vec3d(1, 0, 0));
    m_mesh.set_vertex(getGeomV(bottomvec[2]), Vec3d(1, 1, 0));
    m_mesh.set_vertex(getGeomV(bottomvec[3]), Vec3d(0, 1, 0));
    m_mesh.set_vertex(getGeomV(topvec[0]), Vec3d(0, 0, 1));
    m_mesh.set_vertex(getGeomV(topvec[1]), Vec3d(1, 0, 1));
    m_mesh.set_vertex(getGeomV(topvec[2]), Vec3d(1, 1, 1));
    m_mesh.set_vertex(getGeomV(topvec[3]), Vec3d(0, 1, 1));
    std::vector<VertexHandle> cell_vertices{getGeomV(bottomvec[0]), getGeomV(bottomvec[1]), getGeomV(bottomvec[2]), 
    getGeomV(bottomvec[3]), getGeomV(topvec[0]), getGeomV(topvec[3]), getGeomV(topvec[2]), getGeomV(topvec[1])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
}

ErrorCode MyMesh::AddOneCellCase1(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
    using namespace OpenVolumeMesh;
    #ifndef NDEBUG
        assert(_nbhf_vec.size()==1 && _nbc_vec.size()==1);
    #endif
    
    std::vector<VertexHandle> bottomvec(m_topomesh.halfface_vertices(_nbhf_vec[0]).first, m_topomesh.halfface_vertices(_nbhf_vec[0]).second);
    std::vector<VertexHandle> topvec = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], bottomvec);
    #ifndef NDEBUG
        assert(topvec.size()==4);
    #endif

    for (int i=0; i<4; ++i) {
        VertexHandle p0 = getGeomV(bottomvec[(i+3)%4]);
        VertexHandle p1 = getGeomV(bottomvec[i]);
        VertexHandle p2 = getGeomV(bottomvec[(i+1)%4]);
        auto p3_mid = cross(m_mesh.vertex(p1)-m_mesh.vertex(p0), m_mesh.vertex(p2)-m_mesh.vertex(p1));
        p3_mid.normalize_cond();
        p3_mid -= m_mesh.vertex(p1);
        m_mesh.set_vertex(getGeomV(topvec[i]), p3_mid);
    }
    std::vector<VertexHandle> cell_vertices{getGeomV(bottomvec[0]), getGeomV(bottomvec[1]), getGeomV(bottomvec[2]), 
    getGeomV(bottomvec[3]), getGeomV(topvec[0]), getGeomV(topvec[3]), getGeomV(topvec[2]), getGeomV(topvec[1])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
}
ErrorCode MyMesh::AddOneCellCase2(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
    using namespace OpenVolumeMesh;
    #ifndef NDEBUG
        assert(_nbhf_vec.size()==2 && _nbc_vec.size()==2);
    #endif
    std::vector<VertexHandle> wing1, wing2;
    VertexHandle s1, s2;
    //find common edge
    {
        std::vector<VertexHandle> bf1(m_topomesh.halfface_vertices(_nbhf_vec[0]).first, m_topomesh.halfface_vertices(_nbhf_vec[0]).second);
        std::vector<VertexHandle> bf2(m_topomesh.halfface_vertices(_nbhf_vec[1]).first, m_topomesh.halfface_vertices(_nbhf_vec[1]).second);
        std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
        std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
        int st1, st2;
        for (int i=0; i<4; ++i) {
            if (bf2_set.count(bf1[i]) && bf2_set.count(bf1[(i+1)%4])) {
                st1 = i;
                break;
            }
        }
        for (int i=0; i<4; ++i) {
            if (bf1_set.count(bf2[i]) && bf1_set.count(bf2[(i+1)%4])) {
                st2 = i;
                break;
            }
        }
        wing1 = {bf1[st1], bf1[(st1+1)%4], bf1[(st1+2)%4], bf1[(st1+3)%4]};
        wing2 = {bf2[st2], bf2[(st2+1)%4], bf2[(st2+2)%4], bf2[(st2+3)%4]};
        auto _top = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], {wing1[3], wing1[2]});
        s1 = _top[0], s2 = _top[1];
    }
    //deform

    //end
    VertexHandle p0 = getGeomV(wing1[0]), p1 = getGeomV(wing1[3]), p2 = getGeomV(wing2[2]);
    VertexHandle p3 = getGeomV(wing2[0]), p4 = getGeomV(wing2[3]), p5 = getGeomV(wing1[2]);
    auto s1_mid = m_mesh.vertex(p1)+m_mesh.vertex(p2)-2*m_mesh.vertex(p0);
    s1_mid.normalize_cond();
    s1_mid -= m_mesh.vertex(p0);
    m_mesh.set_vertex(getGeomV(s1), s1_mid);
    auto s2_mid = m_mesh.vertex(p4)+m_mesh.vertex(p5)-2*m_mesh.vertex(p3);
    s2_mid.normalize_cond();
    s2_mid -= m_mesh.vertex(p3);
    m_mesh.set_vertex(getGeomV(s2), s2_mid);

    std::vector<VertexHandle> cell_vertices{getGeomV(wing1[0]), getGeomV(wing1[1]), getGeomV(wing1[2]), 
    getGeomV(wing1[3]), getGeomV(wing2[2]), getGeomV(s1), getGeomV(s2), getGeomV(wing2[3])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
}
ErrorCode MyMesh::AddOneCellCase3(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
    using namespace OpenVolumeMesh;
    #ifndef NDEBUG
        assert(_nbhf_vec.size()==3 && _nbc_vec.size()==3);
    #endif
    std::vector<VertexHandle> fan1, fan2, fan3;
    VertexHandle jade;
    //find common edge
    {
        std::vector<VertexHandle> bf1(m_topomesh.halfface_vertices(_nbhf_vec[0]).first, m_topomesh.halfface_vertices(_nbhf_vec[0]).second);
        std::vector<VertexHandle> bf2(m_topomesh.halfface_vertices(_nbhf_vec[1]).first, m_topomesh.halfface_vertices(_nbhf_vec[1]).second);
        std::vector<VertexHandle> bf3(m_topomesh.halfface_vertices(_nbhf_vec[2]).first, m_topomesh.halfface_vertices(_nbhf_vec[2]).second);
        std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
        std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
        std::set<VertexHandle> bf3_set(bf3.begin(), bf3.end());
        int st1, st2, st3;
        for (int i=0; i<4; ++i) {
            if (bf2_set.count(bf1[i]) && bf3_set.count(bf1[i])) {
                st1 = i;
                break;
            }
        }
        for (int i=0; i<4; ++i) {
            if (bf1_set.count(bf2[i]) && bf3_set.count(bf2[i])) {
                st2 = i;
                break;
            }
        }
        for (int i=0; i<4; ++i) {
            if (bf1_set.count(bf3[i]) && bf2_set.count(bf3[i])) {
                st3 = i;
                break;
            }
        }
        fan1 = {bf1[st1], bf1[(st1+1)%4], bf1[(st1+2)%4], bf1[(st1+3)%4]};
        fan2 = {bf2[st2], bf2[(st2+1)%4], bf2[(st2+2)%4], bf2[(st2+3)%4]};
        fan3 = {bf3[st2], bf3[(st2+1)%4], bf3[(st2+2)%4], bf3[(st2+3)%4]};
        auto _top = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], {fan1[2]});
        jade = _top[0];
    }
    //deform

    //end
    VertexHandle p10 = getGeomV(fan1[2]), p11 = getGeomV(fan1[1]), p12 = getGeomV(fan1[3]);
    VertexHandle p20 = getGeomV(fan2[2]), p21 = getGeomV(fan2[1]), p22 = getGeomV(fan2[3]);
    VertexHandle p30 = getGeomV(fan3[2]), p31 = getGeomV(fan3[1]), p32 = getGeomV(fan3[3]);
    auto s1_mid = (m_mesh.vertex(p11)+m_mesh.vertex(p12)-2*m_mesh.vertex(p10)).normalize_cond()-m_mesh.vertex(p10);
    auto s2_mid = (m_mesh.vertex(p21)+m_mesh.vertex(p22)-2*m_mesh.vertex(p20)).normalize_cond()-m_mesh.vertex(p20);
    auto s3_mid = (m_mesh.vertex(p31)+m_mesh.vertex(p32)-2*m_mesh.vertex(p30)).normalize_cond()-m_mesh.vertex(p30);
    m_mesh.set_vertex(getGeomV(jade), (s1_mid+s2_mid+s3_mid)/3);

    auto _fan1_top = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], fan1);

    std::vector<VertexHandle> cell_vertices{getGeomV(fan1[0]), getGeomV(fan1[1]), getGeomV(fan1[2]), 
    getGeomV(fan1[3]), getGeomV(_fan1_top[0]), getGeomV(_fan1_top[3]), getGeomV(_fan1_top[2]), getGeomV(_fan1_top[1])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
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


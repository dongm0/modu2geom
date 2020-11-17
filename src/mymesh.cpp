#include "mymesh.h"
#include <Mesquite/Mesquite_ArrayMesh.hpp>
#include <Mesquite_IdealShapeTarget.hpp>
#include <Mesquite_TMetric.hpp>
#include <Mesquite_TQualityMetric.hpp>
#include <Mesquite_ElementPMeanP.hpp>
#include <Mesquite_InstructionQueue.hpp>
#include <Mesquite_ConjugateGradient.hpp>

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
    {
        //求外向法向方向
        {
            
        }
    }
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
    //find bottom structure and center, jade is the undefined vertex
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
        fan3 = {bf3[st3], bf3[(st3+1)%4], bf3[(st3+2)%4], bf3[(st3+3)%4]};
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
    using namespace OpenVolumeMesh;
    #ifndef NDEBUG
        assert(_nbhf_vec.size()==3 && _nbc_vec.size()==3);
    #endif
    std::vector<VertexHandle> bottom_vec, wing1, wing2;
    
    //find bottom structure
    {
        int num_b, num_w1, num_w2;
        if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) == _nbhf_vec[1]) {
            num_b = 2, num_w1 = 0, num_w2 = 1;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) == _nbhf_vec[2]) {
            num_b = 1, num_w1 = 0, num_w2 = 2;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[1], _ch) == _nbhf_vec[2]) {
            num_b = 0, num_w1 = 1, num_w2 = 2;
        }
        else {
            return ErrorCode::failed;
        }
        std::vector<VertexHandle> bf1(m_topomesh.halfface_vertices(_nbhf_vec[num_b]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_b]).second);
        std::vector<VertexHandle> bf2(m_topomesh.halfface_vertices(_nbhf_vec[num_w1]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_w1]).second);
        std::vector<VertexHandle> bf3(m_topomesh.halfface_vertices(_nbhf_vec[num_w2]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_w2]).second);
        std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
        std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
        std::set<VertexHandle> bf3_set(bf3.begin(), bf3.end());
        int st1, st2, st3;
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
        for (int i=0; i<4; ++i) {
            if (bf1_set.count(bf3[i]) && bf1_set.count(bf3[(i+1)%4])) {
                st3 = i;
                break;
            }
        }
        bottom_vec = {bf1[st1], bf1[(st1+1)%4], bf1[(st1+2)%4], bf1[(st1+3)%4]};
        wing1 = {bf2[st2], bf2[(st2+1)%4], bf2[(st2+2)%4], bf2[(st2+3)%4]};
        wing2 = {bf3[st3], bf3[(st3+1)%4], bf3[(st3+2)%4], bf3[(st3+3)%4]};
    }
    //deform

    //end
    std::vector<VertexHandle> cell_vertices{getGeomV(bottom_vec[0]), getGeomV(bottom_vec[1]), getGeomV(bottom_vec[2]), 
    getGeomV(bottom_vec[3]), getGeomV(wing1[2]), getGeomV(wing2[3]), getGeomV(wing2[2]), getGeomV(wing1[1])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
}
ErrorCode MyMesh::AddOneCellCase5(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
    using namespace OpenVolumeMesh;
    #ifndef NDEBUG
        assert(_nbhf_vec.size()==4 && _nbc_vec.size()==4);
    #endif
    std::vector<VertexHandle> wing1, wing2, side1, side2;
    //find bottom structure
    {
        int num_w1, num_w2, num_s1, num_s2;
        if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) == _nbhf_vec[1]) {
            num_w1 = 0, num_w2 = 1, num_s1 = 2, num_s2 = 3;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) == _nbhf_vec[2]) {
            num_w1 = 0, num_w2 = 2, num_s1 = 1, num_s2 = 3;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) == _nbhf_vec[3]) {
            num_w1 = 0, num_w2 = 3, num_s1 = 1, num_s2 = 2;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[1], _ch) == _nbhf_vec[2]) {
            num_w1 = 1, num_w2 = 2, num_s1 = 0, num_s2 = 3;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[1], _ch) == _nbhf_vec[3]) {
            num_w1 = 1, num_w2 = 3, num_s1 = 0, num_s2 = 2;
        }
        else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[2], _ch) == _nbhf_vec[3]) {
            num_w1 = 2, num_w2 = 3, num_s1 = 0, num_s2 = 1;
        }
        std::vector<VertexHandle> bf1(m_topomesh.halfface_vertices(_nbhf_vec[num_w1]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_w1]).second);
        std::vector<VertexHandle> bf2(m_topomesh.halfface_vertices(_nbhf_vec[num_w2]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_w2]).second);
        std::vector<VertexHandle> sf1(m_topomesh.halfface_vertices(_nbhf_vec[num_s1]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_s1]).second);
        std::vector<VertexHandle> sf2(m_topomesh.halfface_vertices(_nbhf_vec[num_s2]).first, m_topomesh.halfface_vertices(_nbhf_vec[num_s2]).second);
        std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
        std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
        int st1 = -1, st2 = -1, st3 = -1, st4 = -1;
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
        for (int i=0; i<4; ++i) {
            if (bf2_set.count(sf1[i]) && bf2_set.count(sf1[(i+1)%4]) && bf1_set.count(sf1[(i+2)%4])) {
                sf1.swap(sf2);
                break;
            }
        }
        for (int i=0; i<4; ++i) {
            if (sf1[i] == bf1[(st1+2)%4]) {
                st3 = i;
                break;
            }
        }
        for (int i=0; i<4; ++i) {
            if (sf2[i] == bf2[(st2+2)%4]) {
                st4 = i;
                break;
            }
        }
        wing1 = {bf1[st1], bf1[(st1+1)%4], bf1[(st1+2)%4], bf1[(st1+3)%4]};
        wing2 = {bf2[st2], bf2[(st2+1)%4], bf2[(st2+2)%4], bf2[(st2+3)%4]};
        side1 = {sf1[st3], sf1[(st3+1)%4], sf1[(st3+2)%4], sf1[(st3+3)%4]};
        side2 = {sf2[st4], sf2[(st4+1)%4], sf2[(st4+2)%4], sf2[(st4+3)%4]};
    }
    //deform

    //end
    std::vector<VertexHandle> cell_vertices{getGeomV(wing1[0]), getGeomV(wing1[1]), getGeomV(wing1[2]), 
    getGeomV(wing1[3]), getGeomV(wing1[3]), getGeomV(wing2[2]), getGeomV(side2[3]), getGeomV(side1[3])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
}
ErrorCode MyMesh::AddOneCellCase6(const OpenVolumeMesh::CellHandle &_ch, 
const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec, 
const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
    using namespace OpenVolumeMesh;
    #ifndef NDEBUG
        assert(_nbhf_vec.size()==5 && _nbc_vec.size()==5);
    #endif
    std::vector<VertexHandle> bottom_vec;
    std::vector<std::vector<VertexHandle>> wings;
    //find bottom structure
    {
        HalfFaceHandle bottom_hf;
        std::set<HalfFaceHandle> _nbhf_set(_nbhf_vec.begin(), _nbhf_vec.end());
        for (const auto &x : _nbhf_vec) {
            if (_nbhf_set.count(m_topomesh.opposite_halfface_handle_in_cell(x, _ch))==0) {
                bottom_vec.assign(m_topomesh.halfface_vertices(x).first, m_topomesh.halfface_vertices(x).second);
                bottom_hf = x;
            }
        }
        auto hes = m_topomesh.halfface(bottom_hf).halfedges();
        auto b_range = m_topomesh.halfface_vertices(bottom_hf);
        std::vector<VertexHandle> vs(b_range.first, b_range.second); 
        for (int i=0; i<4; ++i) {
            auto _wing_hf = m_topomesh.adjacent_halfface_in_cell(bottom_hf, hes[i]);
            int st = -1;
            auto s_range = m_topomesh.halfface_vertices(_wing_hf);
            std::vector<VertexHandle> _side_vs(s_range.first, s_range.second);
            for (int j=0; j<4; ++j) {
                if (_side_vs[j] == vs[(i+1)%4]) {
                    st = j;
                    break;
                }
            }
            wings.push_back({_side_vs[st], _side_vs[(st+1)%4], _side_vs[(st+2)%4], _side_vs[(st+3)%4]});
        }
    }
    //deform

    //end
    std::vector<VertexHandle> cell_vertices{getGeomV(bottom_vec[0]), getGeomV(bottom_vec[1]), getGeomV(bottom_vec[2]), 
    getGeomV(bottom_vec[3]), getGeomV(wings[0][2]), getGeomV(wings[3][2]), getGeomV(wings[2][2]), getGeomV(wings[1][2])};
    auto _geomch = m_mesh.add_cell(cell_vertices);
    m_tm2m_mapping[_ch] = _geomch;
    m_m2tm_mapping[_geomch] = _ch;

    return ErrorCode::succeed;
}

ErrorCode MyMesh::Optimize() {
    using namespace Mesquite;
    std::vector<double> coords;
    for (auto v_it = m_mesh.vertices_begin(); v_it!=m_mesh.vertices_end(); ++v_it) {
        auto c = m_mesh.vertex(*v_it);
        coords.push_back(c[0]);
        coords.push_back(c[1]);
        coords.push_back(c[2]);
    }
    std::vector<int> fixedflag;
    std::vector<int> connection;

    ArrayMesh msqmesh(3, m_mesh.n_vertices(), coords.data(), fixedflag.data(), m_mesh.n_cells(), HEXAHEDRON, connection.data());
    MsqError err;
    
    IdealShapeTarget target;
    //TShapeSizeB1 m1;
    TShapeSizeB3 m2;
    //TSum mymetric(&m1, &m2);
    TQualityMetric metric_0(&target, &m2);
    ElementPMeanP metric(1.0, &metric_0);
    PMeanPTemplate obj_func_opt(1.0, &metric);
    ConjugateGradient improver(&obj_func_opt);;
    improver.use_global_patch();
    //improver.set_inner_termination_criterion(&e);
    InstructionQueue queue;
    queue.set_master_quality_improver(&improver, err);
    queue.run_instructions(&msqmesh, err);

    int i = 0;
    for (auto v_it = m_mesh.vertices_begin(); v_it!=m_mesh.vertices_end(); ++v_it) {
        m_mesh.set_vertex(*v_it, OpenVolumeMesh::Geometry::Vec3d(coords[i*3], coords[i*3+1], coords[i*3+2]));
    }

    return ErrorCode::succeed;
}
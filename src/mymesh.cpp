
#include "mymesh.h"
#include "OVMVtkHexIO.h"
#include "msqoptimizer.h"
#include "ovmwrap.h"
#include <algorithm>

namespace {

static OpenVolumeMesh::Vec3d
projectVectoDirection(const OpenVolumeMesh::Vec3d &v,
                      const OpenVolumeMesh::Vec3d &ref) {
  auto ref_norm = ref.normalized();
  return ref_norm * (v | ref_norm);
}

static OpenVolumeMesh::Vec3d
projectVectoPlane(const OpenVolumeMesh::Vec3d &v,
                  const OpenVolumeMesh::Vec3d &ref) {
  return v - projectVectoDirection(v, ref);
}

static inline std::vector<std::vector<OpenVolumeMesh::Vec3d>>
untangleBottomFace(const std::vector<untangleData> &uData,
                   double minangle = 9) {
  using namespace OpenVolumeMesh;
  const double my_pi = atan(1) * 4;

  //可以考虑都单位化
  int fsize = uData.size();
  std::vector<std::vector<Vec3d>> _res;
  for (int i = 0; i < fsize; ++i) {
    std::vector<Vec3d> _res1c;
    for (int j = 0; j < uData[i].vertices.size(); ++j) {
      Vec3d corner = uData.at(i).vertices.at(j) - uData.at(i).assistV_d;
      corner = corner -
               ((corner | uData.at(i).bline_d) / uData.at(i).bline_d.length()) *
                   (uData.at(i).bline_d / uData.at(i).bline_d.length());
      Vec3d base = uData.at(i).vertices.at(j) - corner;

      double alpha = acos((corner | uData.at(i).face_d) /
                          (corner.length() * uData.at(i).face_d.length()));
      if ((corner | uData.at(i).norm_d) < 0)
        alpha = 2 * my_pi - alpha;

      if (alpha < (my_pi / (2 * minangle)) or alpha > (5 * my_pi / 4)) {
        double theta = my_pi / (2 * minangle) - alpha;
        if (theta < 0)
          theta += 2 * my_pi;
        corner = cos(theta) * corner +
                 (1 - cos(theta)) * (corner | uData.at(i).bline_d) *
                     uData.at(i).bline_d -
                 sin(theta) * (uData.at(i).bline_d % corner);
      } else if (alpha > (my_pi / 3) and alpha <= (5 * my_pi / 4)) {
        double theta = alpha - my_pi / 3;
        corner = cos(theta) * corner +
                 (1 - cos(theta)) * (corner | uData.at(i).bline_d) *
                     uData.at(i).bline_d +
                 sin(theta) * (uData.at(i).bline_d % corner);
      }

      _res1c.push_back(base + corner);
    }
    _res.push_back(_res1c);
  }
  return _res;
}
} // namespace

bool MyMesh::ReadTopoFrom(const std::string &filename) {
  int dotpos = -1;
  for (int i = 0; i < filename.size(); ++i) {
    if (filename.at(filename.size() - 1 - i) == '.') {
      dotpos = filename.size() - 1 - i;
      break;
    }
  }
  if (dotpos == -1 || dotpos == 0)
    return false;
  auto ext = filename.substr(dotpos);
  if (ext == ".vtk") {
    return ReadTopoFromVTKFile(filename);
  } else if (ext == ".txt") {
    return ReadTopoFromFile(filename);
  } else {
    return false;
  }
}

bool MyMesh::ReadTopoFromFile(const std::string &filename) {
  using namespace OpenVolumeMesh;

  std::ifstream fin;
  fin.open(filename);
  if (fin.fail()) {
    return false;
  }

  uint32_t _cnum = 0;
  uint32_t _vnum = 0;
  int _tmp;
  fin >> _tmp;
  _cnum = _tmp;
  m_cells.resize(_cnum);
  for (uint32_t i = 0; i < _cnum; ++i) {
    for (uint32_t j = 0; j < 8; ++j) {
      fin >> _tmp;
      m_cells[i][j] = _tmp;
      if (m_vids.count(_tmp) == 0) {
        auto vid = m_topomesh.add_vertex();
        m_vids[_tmp] = vid;
      }
    }
  }
  fin.close();
  _vnum = m_vids.size();

  // add topo cells
  std::vector<VertexHandle> _cell(8);
  for (uint32_t i = 0; i < _cnum; ++i) {
    for (uint32_t j = 0; j < 8; ++j) {
      _cell[j] = m_vids.at(m_cells[i][j]);
    }
    std::swap(_cell[3], _cell[1]);
#ifdef OVM_TOPOLOGY_CHECK
    m_topomesh.add_cell(_cell, true);
#else
    m_topomesh.add_cell(_cell, false); // false
#endif
  }
  m_inner_p0[0] = 0, m_inner_p0[1] = 0, m_inner_p0[2] = 0, m_inner_p1[0] = 1,
  m_inner_p1[1] = 0, m_inner_p1[2] = 0, m_inner_p2[0] = 1, m_inner_p2[1] = 1,
  m_inner_p2[2] = 0;
  m_inner_length = 1;
  return true;
}

bool MyMesh::ReadTopoFromVTKFile(const std::string &filename) {
  using namespace OpenVolumeMesh;

  std::ifstream fin;
  fin.open(filename);
  if (fin.fail()) {
    return false;
  }
  std::stringstream ss;
  std::string line;
  std::string stmp;

  uint32_t _cnum = 0;
  uint32_t _vnum = 0;

  std::getline(fin, line);
  std::getline(fin, line);
  std::getline(fin, line);
  if (line != "ASCII") {
    throw std::runtime_error("ascii only!");
  }
  std::getline(fin, line);
  std::getline(fin, line);
  ss.clear();
  ss.str(line);
  ss >> stmp >> _vnum >> stmp;

  for (uint32_t i = 0; i < _vnum; ++i) {
    std::getline(fin, line);
    m_topomesh.add_vertex();
    if (i == 0) {
      ss.clear();
      ss.str(line);
      double _x, _y, _z;
      ss >> _x >> _y >> _z;
      m_inner_p0[0] = _x, m_inner_p0[1] = _y, m_inner_p0[2] = _z;
    } else if (i == 1) {
      ss.clear();
      ss.str(line);
      double _x, _y, _z;
      ss >> _x >> _y >> _z;
      m_inner_p1[0] = _x, m_inner_p1[1] = _y, m_inner_p1[2] = _z;
      m_inner_length = (m_inner_p0 - m_inner_p1).length();
    } else if (i == 2) {
      ss.clear();
      ss.str(line);
      double _x, _y, _z;
      ss >> _x >> _y >> _z;
      OpenVolumeMesh::Geometry::Vec3d _p2;
      _p2[0] = _x, _p2[1] = _y, _p2[2] = _z;
      m_inner_p2 = m_inner_p1 + (projectVectoPlane(_p2 - m_inner_p1,
                                                   m_inner_p1 - m_inner_p0))
                                        .normalized() *
                                    m_inner_length;
    }
  }
  std::getline(fin, line);
  ss.clear();
  ss.str(line);
  ss >> stmp >> _cnum >> stmp;

  int _tmp;
  m_cells.resize(_cnum);
  for (uint32_t i = 0; i < _cnum; ++i) {
    fin >> _tmp;
    for (uint32_t j = 0; j < 8; ++j) {
      fin >> _tmp;
      m_cells[i][j] = _tmp;
      if (m_vids.count(_tmp) == 0) {
        // auto vid = m_topomesh.add_vertex();
        m_vids[_tmp] = OpenVolumeMesh::VertexHandle(_tmp);
      }
    }
  }
  fin.close();
  _vnum = m_vids.size();

  // add topo cells
  std::vector<VertexHandle> _cell(8);
  for (uint32_t i = 0; i < _cnum; ++i) {
    for (uint32_t j = 0; j < 8; ++j) {
      _cell[j] = m_vids.at(m_cells[i][j]);
    }
    std::swap(_cell[3], _cell[1]);
#ifdef OVM_TOPOLOGY_CHECK
    m_topomesh.add_cell(_cell, true);
#else
    m_topomesh.add_cell(_cell, false); // false
#endif
  }
  return true;
}

bool MyMesh::checkTopo() {
  std::vector<int> tmp(10, 0);
  for (auto e : m_topomesh.edges()) {
    if (!m_topomesh.is_boundary(e)) {
      tmp[m_topomesh.valence(e)] += 1;
    }
  }
  for (auto x : tmp)
    std::cout << x << " ";
  std::cout << std::endl;
  return true;
}

bool MyMesh::GenerateOrder() {
  for (auto _ch : m_topomesh.cells()) {
    m_generate_order.push_back(_ch);
  }
  return true;
}

bool MyMesh::WriteGeomToVTKFile(const std::string &filename) {
  // MsqOperator::Instance().Ovm2MsqOut(m_mesh, filename);
  std::ofstream fout;
  fout.open(filename);
  OVMWriteHexMesh(m_mesh, fout);
  fout.close();
  return true;
}
bool MyMesh::WriteGeomToVTKFile(
    const std::string &filename,
    std::vector<OpenVolumeMesh::VertexHandle> &_tagged) {
  std::ofstream fout;
  fout.open(filename);
  OVMWriteHexMesh(m_mesh, fout, _tagged);
  fout.close();
  return true;
}
bool MyMesh::WriteGeomToVTKFileUseTopoMesh(const std::string &filename) {
  std::ofstream stream;
  stream.open(filename);
  stream << "# vtk DataFile Version 3.0\ngenerate by dongmo\nASCII\nDATASET "
            "UNSTRUCTURED_GRID"
         << std::endl;
  stream << "POINTS " << m_topomesh.n_vertices() << " double" << std::endl;
  for (auto &v : m_topomesh.vertices()) {
    auto p = m_mesh.vertex(m_tm2m_v_mapping[v]);
    stream << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }
  stream << "CELLS " << m_topomesh.n_cells() << " " << 9 * m_topomesh.n_cells()
         << std::endl;
  for (int i = 0; i < m_topomesh.n_cells(); ++i) {
    stream << "8 ";
    for (int j = 0; j < 8; ++j) {
      stream << m_vids.at(m_cells[i][j]).idx() << " ";
    }
    stream << std::endl;
  }
  stream << "CELL_TYPES " << m_topomesh.n_cells() << std::endl;
  for (uint32_t i = 0; i < m_topomesh.n_cells(); ++i) {
    stream << "12" << std::endl;
  }
  for (auto x : m_vids) {
    std::cout << x.first << " " << x.second.idx() << std::endl;
  }
}

bool MyMesh::GenerateOneCell(const OpenVolumeMesh::CellHandle &_ch) {
  using namespace OpenVolumeMesh;

  if (m_tm2m_mapping.find(_ch) != m_tm2m_mapping.end()) {
    auto _tmpch = m_tm2m_mapping[_ch];
    if (m_m2tm_mapping.find(_tmpch) != m_m2tm_mapping.end() &&
        m_m2tm_mapping[_tmpch] == _ch)
      return true;
    else
      return false;
  }
  // 6 cases
  int num_nbh = 0, casenum = -1;
  std::vector<CellHandle> _nbh;
  std::vector<HalfFaceHandle> _nbh_hf;
  {
    for (auto hf_handle : m_topomesh.cell(_ch).halffaces()) {
      auto opposite_hf = m_topomesh.opposite_halfface_handle(hf_handle);
      if (opposite_hf.is_valid() && opposite_hf.idx() != hf_handle.idx()) {
        auto _nbch = m_topomesh.incident_cell(opposite_hf);
        if (m_tm2m_mapping.find(_nbch) != m_tm2m_mapping.end()) {
          num_nbh++;
          _nbh.push_back(_nbch);
          _nbh_hf.push_back(hf_handle);
        }
      }
    }

    if (num_nbh == 3) {
      if (m_topomesh.opposite_halfface_handle_in_cell(_nbh_hf[0], _ch) ==
              _nbh_hf[1] or
          m_topomesh.opposite_halfface_handle_in_cell(_nbh_hf[0], _ch) ==
              _nbh_hf[2] or
          m_topomesh.opposite_halfface_handle_in_cell(_nbh_hf[1], _ch) ==
              _nbh_hf[2])
        casenum = 4;
      else
        casenum = 3;
    } else {
      if (num_nbh == 0 or num_nbh == 1 or num_nbh == 2)
        casenum = num_nbh;
      else if (num_nbh == 4 or num_nbh == 5)
        casenum = num_nbh + 1;
      else
        return false;
    }

    if (casenum == 0)
      AddOneCellCase0(_ch, _nbh, _nbh_hf);
    else if (casenum == 1)
      AddOneCellCase1(_ch, _nbh, _nbh_hf);
    else if (casenum == 2)
      AddOneCellCase2(_ch, _nbh, _nbh_hf);
    else if (casenum == 3)
      AddOneCellCase3(_ch, _nbh, _nbh_hf);
    else if (casenum == 4)
      AddOneCellCase4(_ch, _nbh, _nbh_hf);
    else if (casenum == 5)
      AddOneCellCase5(_ch, _nbh, _nbh_hf);
    else if (casenum == 6)
      AddOneCellCase6(_ch, _nbh, _nbh_hf);
    // WriteGeomToVTKFile("tmp.vtk");
  }
  return true;
}

bool MyMesh::AddOneCellCase0(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
#ifndef NDEBUG
  assert(_nbc_vec.size() == 0 && _nbhf_vec.size() == 0);
#endif
  using namespace OpenVolumeMesh;
  auto bottomface = m_topomesh.xfront_halfface(_ch);
  auto b_range = m_topomesh.halfface_vertices(bottomface);
  std::vector<VertexHandle> bottomvec(b_range.first, b_range.second);
  std::vector<VertexHandle> topvec =
      opposite_vertex_in_cell(m_topomesh, _ch, bottomface, bottomvec);

#ifndef NDEBUG
  assert(topvec.size() == 4);
#endif
  std::vector<VertexHandle> geomvertices(8);
  geomvertices[0] = m_mesh.add_vertex(m_inner_p0);
  geomvertices[1] = m_mesh.add_vertex(m_inner_p1);
  geomvertices[2] = m_mesh.add_vertex(m_inner_p2);
  geomvertices[3] = m_mesh.add_vertex(m_inner_p2 + m_inner_p0 - m_inner_p1);

  auto vertical =
      (m_inner_p2 - m_inner_p1).cross(m_inner_p0 - m_inner_p1).normalized();
  geomvertices[4] = m_mesh.add_vertex(m_mesh.vertex(geomvertices[0]) +
                                      m_inner_length * vertical);
  geomvertices[5] = m_mesh.add_vertex(m_mesh.vertex(geomvertices[1]) +
                                      m_inner_length * vertical);
  geomvertices[6] = m_mesh.add_vertex(m_mesh.vertex(geomvertices[2]) +
                                      m_inner_length * vertical);
  geomvertices[7] = m_mesh.add_vertex(m_mesh.vertex(geomvertices[3]) +
                                      m_inner_length * vertical);
  /*
  geomvertices[4] = m_mesh.add_vertex(
      m_inner_p0 +
      m_inner_length *
          (m_inner_p2 + m_inner_p0 - m_inner_p1 - m_inner_p1).normalized());
  geomvertices[5] = m_mesh.add_vertex(
      m_inner_p1 +
      m_inner_length *
          (m_inner_p2 + m_inner_p0 - m_inner_p1 - m_inner_p1).normalized());
  geomvertices[6] = m_mesh.add_vertex(
      m_inner_p2 +
      m_inner_length *
          (m_inner_p2 + m_inner_p0 - m_inner_p1 - m_inner_p1).normalized());
  geomvertices[7] = m_mesh.add_vertex(
      m_inner_p2 + m_inner_p0 - m_inner_p1 +
      m_inner_length *
          (m_inner_p2 + m_inner_p0 - m_inner_p1 - m_inner_p1).normalized());
          */
  m_tm2m_v_mapping[bottomvec[0]] = geomvertices[0];
  m_tm2m_v_mapping[bottomvec[1]] = geomvertices[1];
  m_tm2m_v_mapping[bottomvec[2]] = geomvertices[2];
  m_tm2m_v_mapping[bottomvec[3]] = geomvertices[3];
  m_tm2m_v_mapping[topvec[0]] = geomvertices[4];
  m_tm2m_v_mapping[topvec[1]] = geomvertices[5];
  m_tm2m_v_mapping[topvec[2]] = geomvertices[6];
  m_tm2m_v_mapping[topvec[3]] = geomvertices[7];
  m_m2tm_v_mapping[geomvertices[0]] = bottomvec[0];
  m_m2tm_v_mapping[geomvertices[1]] = bottomvec[1];
  m_m2tm_v_mapping[geomvertices[2]] = bottomvec[2];
  m_m2tm_v_mapping[geomvertices[3]] = bottomvec[3];
  m_m2tm_v_mapping[geomvertices[4]] = topvec[0];
  m_m2tm_v_mapping[geomvertices[5]] = topvec[1];
  m_m2tm_v_mapping[geomvertices[6]] = topvec[2];
  m_m2tm_v_mapping[geomvertices[7]] = topvec[3];

  std::swap(geomvertices[1], geomvertices[3]);
  auto _geomch = m_mesh.add_cell(geomvertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;

  return true;
}

bool MyMesh::AddOneCellCase1(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
  using namespace OpenVolumeMesh;
#ifndef NDEBUG
  assert(_nbhf_vec.size() == 1 && _nbc_vec.size() == 1);
#endif

  std::vector<VertexHandle> bottomvec(
      m_topomesh.halfface_vertices(_nbhf_vec[0]).first,
      m_topomesh.halfface_vertices(_nbhf_vec[0]).second);
  std::vector<VertexHandle> topvec =
      opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], bottomvec);

  for (int i = 0; i < 4; ++i) {
    VertexHandle p0 = getGeomV(bottomvec[(i + 3) % 4]);
    VertexHandle p1 = getGeomV(bottomvec[i]);
    VertexHandle p2 = getGeomV(bottomvec[(i + 1) % 4]);
    auto test_p0_coord = m_mesh.vertex(p0);
    auto test_p1_coord = m_mesh.vertex(p1);
    auto test_p2_coord = m_mesh.vertex(p2);
    auto p3_mid = cross(m_mesh.vertex(p1) - m_mesh.vertex(p0),
                        m_mesh.vertex(p2) - m_mesh.vertex(p1));
    p3_mid.normalize_cond();
    p3_mid = p3_mid * m_inner_length + m_mesh.vertex(p1);
    auto geomvhandle = m_mesh.add_vertex(p3_mid);
    m_m2tm_v_mapping[geomvhandle] = topvec[i];
    m_tm2m_v_mapping[topvec[i]] = geomvhandle;
    // m_mesh.set_vertex(getGeomV(topvec[i]), p3_mid);
  }
  std::vector<VertexHandle> cell_vertices{
      getGeomV(bottomvec[0]), getGeomV(bottomvec[3]), getGeomV(bottomvec[2]),
      getGeomV(bottomvec[1]), getGeomV(topvec[0]),    getGeomV(topvec[1]),
      getGeomV(topvec[2]),    getGeomV(topvec[3])};
  auto _geomch = m_mesh.add_cell(cell_vertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;

  return true;
}
bool MyMesh::AddOneCellCase2(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
  using namespace OpenVolumeMesh;
#ifndef NDEBUG
  assert(_nbhf_vec.size() == 2 && _nbc_vec.size() == 2);
#endif
  std::vector<VertexHandle> wing1, wing2;
  VertexHandle s1, s2;
  //找底部结构
  {
    std::vector<VertexHandle> bf1(
        m_topomesh.halfface_vertices(_nbhf_vec[0]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[0]).second);
    std::vector<VertexHandle> bf2(
        m_topomesh.halfface_vertices(_nbhf_vec[1]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[1]).second);
    std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
    std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
    int st1, st2;
    for (int i = 0; i < 4; ++i) {
      if (bf2_set.count(bf1[i]) && bf2_set.count(bf1[(i + 1) % 4])) {
        st1 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf1_set.count(bf2[i]) && bf1_set.count(bf2[(i + 1) % 4])) {
        st2 = i;
        break;
      }
    }
    wing1 = {bf1[st1], bf1[(st1 + 1) % 4], bf1[(st1 + 2) % 4],
             bf1[(st1 + 3) % 4]};
    wing2 = {bf2[st2], bf2[(st2 + 1) % 4], bf2[(st2 + 2) % 4],
             bf2[(st2 + 3) % 4]};
    auto _top = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0],
                                        {wing1[3], wing1[2]});
    s1 = _top[0], s2 = _top[1];
  }
  // deform
  {
    //求外向法向方向
    Vec3d normdir, wing1dir, wing2dir,
        bottoml; //估计的法向量，以及两边正交的方向，bottoml是两个底面交线的方向
    {
      VertexHandle wbv1 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[0],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[0]), {wing1[0]})[0];
      VertexHandle wbv2 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[0],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[0]), {wing1[1]})[0];
      VertexHandle wbv3 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[1],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[1]), {wing2[0]})[0];
      VertexHandle wbv4 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[1],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[1]), {wing2[1]})[0];
      auto normdir1 = (getCoord_topo(wing1[0]) - getCoord_topo(wbv1) +
                       getCoord_topo(wing1[1]) - getCoord_topo(wbv2))
                          .normalize_cond() +
                      (getCoord_topo(wing2[0]) - getCoord_topo(wbv3) +
                       getCoord_topo(wing2[1]) - getCoord_topo(wbv4))
                          .normalize_cond();
      auto normdir2 = (getCoord_topo(wing1[3]) - getCoord_topo(wing1[0]) +
                       getCoord_topo(wing1[2]) - getCoord_topo(wing1[1]))
                          .normalize_cond() +
                      (getCoord_topo(wing2[3]) - getCoord_topo(wing2[0]) +
                       getCoord_topo(wing2[2]) - getCoord_topo(wing2[1]))
                          .normalize_cond();
      if (normdir1.length() > normdir2.length()) {
        normdir = normdir1;
      } else {
        normdir = normdir2;
      }
      bottoml = getCoord_topo(wing1[0]) - getCoord_topo(wing1[1]);
      // project
      normdir = normdir - ((bottoml | normdir) / bottoml.length()) *
                              (bottoml / bottoml.length());
      wing1dir = (normdir | bottoml) * bottoml + bottoml % normdir;
      wing2dir = -1 * wing1dir;
      normdir.normalize_cond();
      wing1dir.normalize_cond();
      wing2dir.normalize_cond();
      bottoml.normalize_cond();
    }
    // wing的位置调节
    std::map<VertexHandle, Vec3d> fixed;
    {
      std::vector<untangleData> uData(2);
      uData[0] = {normdir,
                  wing1dir,
                  bottoml,
                  getCoord_topo(wing1[0]),
                  {getCoord_topo(wing1[2]), getCoord_topo(wing1[3])}};
      uData[1] = {normdir,
                  wing2dir,
                  -1 * bottoml,
                  getCoord_topo(wing2[0]),
                  {getCoord_topo(wing2[2]), getCoord_topo(wing2[3])}};
      auto new_pos = untangleBottomFace(uData);
      fixed[getGeomV(wing1[0])] = getCoord_topo(wing1[0]);
      fixed[getGeomV(wing1[1])] = getCoord_topo(wing1[1]);
      fixed[getGeomV(wing1[2])] = new_pos[0][0];
      fixed[getGeomV(wing1[3])] = new_pos[0][1];
      fixed[getGeomV(wing2[2])] = new_pos[1][0];
      fixed[getGeomV(wing2[3])] = new_pos[1][1];
      // std::cout << fixed.size() << std::endl;
    }
    //变形
    { ArapOperator::Instance().Deformation(m_mesh, fixed); }
    //WriteGeomToVTKFile("tmp.vtk");
  }
  // end
  VertexHandle p0 = getGeomV(wing1[0]), p1 = getGeomV(wing1[3]),
               p2 = getGeomV(wing2[2]);
  VertexHandle p3 = getGeomV(wing2[0]), p4 = getGeomV(wing2[3]),
               p5 = getGeomV(wing1[2]);
  auto s1_mid = m_mesh.vertex(p1) + m_mesh.vertex(p2) - 2 * m_mesh.vertex(p0);
  s1_mid = projectVectoPlane(s1_mid, m_mesh.vertex(p0) - m_mesh.vertex(p3));
  s1_mid = s1_mid.normalized() * m_inner_length * sqrt(2);
  s1_mid += m_mesh.vertex(p0);
  auto geoms1 = m_mesh.add_vertex(s1_mid);
  m_tm2m_v_mapping[s1] = geoms1;
  m_m2tm_v_mapping[geoms1] = s1;
  // m_mesh.set_vertex(getGeomV(s1), s1_mid);
  auto s2_mid = m_mesh.vertex(p4) + m_mesh.vertex(p5) - 2 * m_mesh.vertex(p3);
  s2_mid = projectVectoPlane(s2_mid, m_mesh.vertex(p0) - m_mesh.vertex(p3));
  s2_mid = s2_mid.normalized() * m_inner_length * sqrt(2);
  s2_mid += m_mesh.vertex(p3);
  auto geoms2 = m_mesh.add_vertex(s2_mid);
  m_tm2m_v_mapping[s2] = geoms2;
  m_m2tm_v_mapping[geoms2] = s2;
  // m_mesh.set_vertex(getGeomV(s2), s2_mid);

  std::vector<VertexHandle> cell_vertices{
      getGeomV(wing1[0]), getGeomV(wing1[3]), getGeomV(wing1[2]),
      getGeomV(wing1[1]), getGeomV(wing2[2]), getGeomV(wing2[3]),
      getGeomV(s2),       getGeomV(s1)};
  auto _geomch = m_mesh.add_cell(cell_vertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;

  return true;
}
bool MyMesh::AddOneCellCase3(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
  using namespace OpenVolumeMesh;
#ifndef NDEBUG
  assert(_nbhf_vec.size() == 3 && _nbc_vec.size() == 3);
#endif
  std::vector<VertexHandle> fan1, fan2, fan3;
  VertexHandle jade;
  // find bottom structure and center, jade is the undefined vertex
  {
    std::vector<VertexHandle> bf1(
        m_topomesh.halfface_vertices(_nbhf_vec[0]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[0]).second);
    std::vector<VertexHandle> bf2(
        m_topomesh.halfface_vertices(_nbhf_vec[1]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[1]).second);
    std::vector<VertexHandle> bf3(
        m_topomesh.halfface_vertices(_nbhf_vec[2]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[2]).second);
    std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
    std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
    std::set<VertexHandle> bf3_set(bf3.begin(), bf3.end());
    int st1, st2, st3;
    for (int i = 0; i < 4; ++i) {
      if (bf2_set.count(bf1[i]) && bf3_set.count(bf1[i])) {
        st1 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf1_set.count(bf2[i]) && bf3_set.count(bf2[i])) {
        st2 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf1_set.count(bf3[i]) && bf2_set.count(bf3[i])) {
        st3 = i;
        break;
      }
    }
    fan1 = {bf1[st1], bf1[(st1 + 1) % 4], bf1[(st1 + 2) % 4],
            bf1[(st1 + 3) % 4]};
    fan2 = {bf2[st2], bf2[(st2 + 1) % 4], bf2[(st2 + 2) % 4],
            bf2[(st2 + 3) % 4]};
    fan3 = {bf3[st3], bf3[(st3 + 1) % 4], bf3[(st3 + 2) % 4],
            bf3[(st3 + 3) % 4]};
    if (fan1[3] != fan2[1]) {
      fan2.swap(fan3);
    }
    jade = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], {fan1[2]})[0];
  }
  // deform
  {
    //结构解扭前的准备
    Vec3d normdir, fan11dir, fan21dir, fan31dir, bottoml1, bottoml2,
        bottoml3; //估计的法向量，以及两边正交的方向，bottoml是两个底面交线的方向
    {
      auto v13 = getCoord_topo(fan1[1]) - getCoord_topo(fan1[0]);
      auto v12 = getCoord_topo(fan2[1]) - getCoord_topo(fan2[0]);
      auto v23 = getCoord_topo(fan3[1]) - getCoord_topo(fan3[0]);
      auto fan1n = v13 % v12;
      auto fan2n = v12 % v23;
      auto fan3n = v23 % v13;
      normdir = fan1n.normalize_cond() + fan2n.normalize_cond() +
                fan3n.normalize_cond();
      normdir.normalize_cond();
      fan11dir = projectVectoPlane(v13, normdir);
      fan11dir.normalize_cond();
      fan21dir = projectVectoPlane(v12, normdir);
      fan21dir.normalize_cond();
      fan31dir = projectVectoPlane(v23, normdir);
      fan31dir.normalize_cond();
      bottoml1 = normdir % fan11dir;
      bottoml2 = normdir % fan21dir;
      bottoml3 = normdir % fan31dir;
    }

    std::map<VertexHandle, Vec3d> fixed;
    {

      std::vector<untangleData> uData(3);
      uData[0] = {normdir,
                  fan11dir,
                  bottoml1,
                  getCoord_topo(fan1[0]),
                  {getCoord_topo(fan1[1])}};
      uData[1] = {normdir,
                  fan21dir,
                  bottoml2,
                  getCoord_topo(fan2[0]),
                  {getCoord_topo(fan2[1])}};
      uData[2] = {normdir,
                  fan31dir,
                  bottoml3,
                  getCoord_topo(fan3[0]),
                  {getCoord_topo(fan3[1])}};
      auto new_pos = untangleBottomFace(uData);

      fixed[getGeomV(fan1[0])] = getCoord_topo(fan1[0]);
      fixed[getGeomV(fan1[1])] = new_pos[0][0];
      fixed[getGeomV(fan2[1])] = new_pos[1][0];
      fixed[getGeomV(fan3[1])] = new_pos[2][0];

      fixed[getGeomV(fan1[2])] =
          new_pos[0][0] + new_pos[1][0] - getCoord_topo(fan1[0]);
      fixed[getGeomV(fan2[2])] =
          new_pos[1][0] + new_pos[2][0] - getCoord_topo(fan1[0]);
      fixed[getGeomV(fan3[2])] =
          new_pos[2][0] + new_pos[0][0] - getCoord_topo(fan1[0]);
    }
    { ArapOperator::Instance().Deformation(m_mesh, fixed); }
  }

  // end

  auto s1_mid =
      m_mesh.vertex(getGeomV(fan1[1])) - m_mesh.vertex(getGeomV(fan1[0]));
  auto s2_mid =
      m_mesh.vertex(getGeomV(fan2[1])) - m_mesh.vertex(getGeomV(fan2[0]));
  auto s3_mid =
      m_mesh.vertex(getGeomV(fan3[1])) - m_mesh.vertex(getGeomV(fan3[0]));
  auto jade_coord =
      (s1_mid + s2_mid + s3_mid).normalized() * sqrt(3) * m_inner_length +
      m_mesh.vertex(getGeomV(fan1[0]));
  auto geomjade = m_mesh.add_vertex(jade_coord);
  m_m2tm_v_mapping[geomjade] = jade;
  m_tm2m_v_mapping[jade] = geomjade;
  // m_mesh.set_vertex(getGeomV(jade), (s1_mid+s2_mid+s3_mid)/3);

  auto _fan1_top = opposite_vertex_in_cell(m_topomesh, _ch, _nbhf_vec[0], fan1);

  std::vector<VertexHandle> cell_vertices{
      getGeomV(fan1[0]), getGeomV(fan1[3]),      getGeomV(fan1[2]),
      getGeomV(fan1[1]), getGeomV(_fan1_top[0]), getGeomV(_fan1_top[1]),
      geomjade,          getGeomV(_fan1_top[3])};
  auto _geomch = m_mesh.add_cell(cell_vertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;
  // WriteGeomToVTKFile("tmp.vtk");

  return true;
}

bool MyMesh::AddOneCellCase4(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
  using namespace OpenVolumeMesh;
#ifndef NDEBUG
  assert(_nbhf_vec.size() == 3 && _nbc_vec.size() == 3);
#endif
  std::vector<VertexHandle> bottom_vec, wing1, wing2;

  // find bottom structure
  {
    int num_b, num_w1, num_w2;
    if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) ==
        _nbhf_vec[1]) {
      num_b = 2, num_w1 = 0, num_w2 = 1;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) ==
               _nbhf_vec[2]) {
      num_b = 1, num_w1 = 0, num_w2 = 2;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[1], _ch) ==
               _nbhf_vec[2]) {
      num_b = 0, num_w1 = 1, num_w2 = 2;
    } else {
      return false;
    }
    std::vector<VertexHandle> bf1(
        m_topomesh.halfface_vertices(_nbhf_vec[num_b]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_b]).second);
    std::vector<VertexHandle> bf2(
        m_topomesh.halfface_vertices(_nbhf_vec[num_w1]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_w1]).second);
    std::vector<VertexHandle> bf3(
        m_topomesh.halfface_vertices(_nbhf_vec[num_w2]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_w2]).second);
    std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
    std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
    std::set<VertexHandle> bf3_set(bf3.begin(), bf3.end());
    int st1, st2, st3;
    for (int i = 0; i < 4; ++i) {
      if (bf2_set.count(bf1[i]) && bf2_set.count(bf1[(i + 1) % 4])) {
        st1 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf1_set.count(bf2[i]) && bf1_set.count(bf2[(i + 1) % 4])) {
        st2 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf1_set.count(bf3[i]) && bf1_set.count(bf3[(i + 1) % 4])) {
        st3 = i;
        break;
      }
    }
    bottom_vec = {bf1[st1], bf1[(st1 + 1) % 4], bf1[(st1 + 2) % 4],
                  bf1[(st1 + 3) % 4]};
    wing1 = {bf2[st2], bf2[(st2 + 1) % 4], bf2[(st2 + 2) % 4],
             bf2[(st2 + 3) % 4]};
    wing2 = {bf3[st3], bf3[(st3 + 1) % 4], bf3[(st3 + 2) % 4],
             bf3[(st3 + 3) % 4]};
  }
  // deform

  {
    //解扭之前的准备
    Vec3d normdir, wing1dir, wing2dir, bottoml1,
        bottoml2; //估计的法向量，以及两边正交的方向
    {
      Vec3d v1 = getCoord_topo(bottom_vec[1]) - getCoord_topo(bottom_vec[0]);
      Vec3d v2 = getCoord_topo(bottom_vec[2]) - getCoord_topo(bottom_vec[1]);
      Vec3d v3 = getCoord_topo(bottom_vec[3]) - getCoord_topo(bottom_vec[2]);
      Vec3d v4 = getCoord_topo(bottom_vec[0]) - getCoord_topo(bottom_vec[3]);
      normdir = v1 % v2 + v2 % v3 + v3 % v4 + v4 % v1;

      bottoml1 = v1, bottoml2 = v3;
      wing1dir = bottoml1 % normdir, wing2dir = bottoml2 % normdir;
    }
    //底面解扭
    std::map<VertexHandle, Vec3d> fixed;
    {
      std::vector<untangleData> uData(2);
      uData[0] = {normdir,
                  wing1dir,
                  bottoml1,
                  getCoord_topo(wing1[0]),
                  {getCoord_topo(wing1[2]), getCoord_topo(wing1[3])}};
      uData[1] = {normdir,
                  wing2dir,
                  bottoml2,
                  getCoord_topo(wing2[0]),
                  {getCoord_topo(wing2[2]), getCoord_topo(wing2[3])}};
      auto new_pos = untangleBottomFace(uData);
      fixed[getGeomV(wing1[0])] = getCoord_topo(wing1[0]);
      fixed[getGeomV(wing1[1])] = getCoord_topo(wing1[1]);
      fixed[getGeomV(wing1[2])] = new_pos[0][0];
      fixed[getGeomV(wing1[3])] = new_pos[0][1];
      fixed[getGeomV(wing2[0])] = getCoord_topo(wing2[0]);
      fixed[getGeomV(wing2[1])] = getCoord_topo(wing2[1]);
      fixed[getGeomV(wing2[2])] = new_pos[1][0];
      fixed[getGeomV(wing2[3])] = new_pos[1][1];
    }
    { ArapOperator::Instance().Deformation(m_mesh, fixed); }
    std::vector<VertexHandle> tagged;
    for (auto x : fixed) {
      tagged.push_back(x.first);
    }

    // WriteGeomToVTKFile("tmp.vtk", tagged);
  }

  // end
  std::vector<VertexHandle> cell_vertices{
      getGeomV(bottom_vec[0]), getGeomV(bottom_vec[3]), getGeomV(bottom_vec[2]),
      getGeomV(bottom_vec[1]), getGeomV(wing1[2]),      getGeomV(wing1[3]),
      getGeomV(wing2[2]),      getGeomV(wing2[3])};
  auto _geomch = m_mesh.add_cell(cell_vertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;

  return true;
}
bool MyMesh::AddOneCellCase5(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
  using namespace OpenVolumeMesh;
#ifndef NDEBUG
  assert(_nbhf_vec.size() == 4 && _nbc_vec.size() == 4);
#endif
  std::vector<VertexHandle> side1, side2, bottom1,
      bottom2; //两个side是opposite的
  int num_s1, num_s2, num_b1, num_b2;
  // find bottom structure
  {
    //这里曾经校正过一遍，可能会出问题，测试的时候多检查
    if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) ==
        _nbhf_vec[1]) {
      num_s1 = 0, num_s2 = 1, num_b1 = 2, num_b2 = 3;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) ==
               _nbhf_vec[2]) {
      num_s1 = 0, num_s2 = 2, num_b1 = 1, num_b2 = 3;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[0], _ch) ==
               _nbhf_vec[3]) {
      num_s1 = 0, num_s2 = 3, num_b1 = 1, num_b2 = 2;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[1], _ch) ==
               _nbhf_vec[2]) {
      num_s1 = 1, num_s2 = 2, num_b1 = 0, num_b2 = 3;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[1], _ch) ==
               _nbhf_vec[3]) {
      num_s1 = 1, num_s2 = 3, num_b1 = 0, num_b2 = 2;
    } else if (m_topomesh.opposite_halfface_handle_in_cell(_nbhf_vec[2], _ch) ==
               _nbhf_vec[3]) {
      num_s1 = 2, num_s2 = 3, num_b1 = 0, num_b2 = 1;
    }
    std::vector<VertexHandle> bf1(
        m_topomesh.halfface_vertices(_nbhf_vec[num_b1]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_b1]).second);
    std::vector<VertexHandle> bf2(
        m_topomesh.halfface_vertices(_nbhf_vec[num_b2]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_b2]).second);
    std::vector<VertexHandle> sf1(
        m_topomesh.halfface_vertices(_nbhf_vec[num_s1]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_s1]).second);
    std::vector<VertexHandle> sf2(
        m_topomesh.halfface_vertices(_nbhf_vec[num_s2]).first,
        m_topomesh.halfface_vertices(_nbhf_vec[num_s2]).second);
    std::set<VertexHandle> bf1_set(bf1.begin(), bf1.end());
    std::set<VertexHandle> bf2_set(bf2.begin(), bf2.end());
    int st1 = -1, st2 = -1, st3 = -1, st4 = -1;
    for (int i = 0; i < 4; ++i) {
      if (bf2_set.count(bf1[i]) && bf2_set.count(bf1[(i + 1) % 4])) {
        st1 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf1_set.count(bf2[i]) && bf1_set.count(bf2[(i + 1) % 4])) {
        st2 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (bf2_set.count(sf1[i]) && bf2_set.count(sf1[(i + 1) % 4]) &&
          bf1_set.count(sf1[(i + 2) % 4])) {
        sf1.swap(sf2);
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (sf1[i] == bf1[(st1 + 2) % 4]) {
        st3 = i;
        break;
      }
    }
    for (int i = 0; i < 4; ++i) {
      if (sf2[i] == bf2[(st2 + 2) % 4]) {
        st4 = i;
        break;
      }
    }
    bottom1 = {bf1[st1], bf1[(st1 + 1) % 4], bf1[(st1 + 2) % 4],
               bf1[(st1 + 3) % 4]};
    bottom2 = {bf2[st2], bf2[(st2 + 1) % 4], bf2[(st2 + 2) % 4],
               bf2[(st2 + 3) % 4]};
    side1 = {sf1[st3], sf1[(st3 + 1) % 4], sf1[(st3 + 2) % 4],
             sf1[(st3 + 3) % 4]};
    side2 = {sf2[st4], sf2[(st4 + 1) % 4], sf2[(st4 + 2) % 4],
             sf2[(st4 + 3) % 4]};
  }
  // deform
  {
    Vec3d normdir, bottom1dir, bottom2dir, side1dir, side2dir, bottoml,
        side1bdir,
        side2bdir; //估计的法向量，以及两边正交的方向，bottoml是两个底面交线的方向
    {
      VertexHandle wbv1 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[num_b1],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[num_b1]),
          {bottom1[0]})[0];
      VertexHandle wbv2 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[num_b1],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[num_b1]),
          {bottom1[1]})[0];
      VertexHandle wbv3 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[num_b2],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[num_b2]),
          {bottom2[0]})[0];
      VertexHandle wbv4 = opposite_vertex_in_cell(
          m_topomesh, _nbc_vec[num_b2],
          m_topomesh.opposite_halfface_handle(_nbhf_vec[num_b2]),
          {bottom2[1]})[0];
      normdir = (getCoord_topo(bottom1[0]) - getCoord_topo(wbv1) +
                 getCoord_topo(bottom1[1]) - getCoord_topo(wbv2))
                    .normalize_cond() +
                (getCoord_topo(bottom2[0]) - getCoord_topo(wbv3) +
                 getCoord_topo(bottom2[1]) - getCoord_topo(wbv4))
                    .normalize_cond();
      bottoml = getCoord_topo(bottom1[0]) - getCoord_topo(bottom1[1]);
      // project
      normdir = normdir - ((bottoml | normdir) / bottoml.length()) *
                              (bottoml / bottoml.length());
      bottom1dir = bottoml % normdir;
      bottom2dir = -1 * bottom1dir;
      side1dir = -bottoml, side2dir = bottoml;
      side1bdir = normdir % side1dir, side2bdir = normdir % side2bdir;
    }
    // wing的位置调节
    std::map<VertexHandle, Vec3d> fixed;
    {
      std::vector<untangleData> uData(2);
      uData[0] = {normdir,
                  bottom1dir,
                  bottoml,
                  getCoord_topo(bottom1[0]),
                  {getCoord_topo(bottom1[2]), getCoord_topo(bottom1[3])}};
      uData[1] = {normdir,
                  bottom2dir,
                  -1 * bottoml,
                  getCoord_topo(bottom2[0]),
                  {getCoord_topo(bottom2[2]), getCoord_topo(bottom2[3])}};
      auto new_pos = untangleBottomFace(uData);
      fixed[getGeomV(bottom1[0])] = getCoord_topo(bottom1[0]);
      fixed[getGeomV(bottom1[1])] = getCoord_topo(bottom1[1]);
      fixed[getGeomV(bottom1[2])] = new_pos[0][0];
      fixed[getGeomV(bottom1[3])] = new_pos[0][1];
      fixed[getGeomV(bottom2[2])] = new_pos[1][0];
      fixed[getGeomV(bottom2[3])] = new_pos[1][1];
      fixed[getGeomV(side1[3])] =
          new_pos[0][0] + new_pos[1][1] - getCoord_topo(bottom1[1]);
      fixed[getGeomV(side2[3])] =
          new_pos[0][1] + new_pos[1][0] - getCoord_topo(bottom1[0]);
    }
    { ArapOperator::Instance().Deformation(m_mesh, fixed); }
    // WriteGeomToVTKFile("tmp.vtk");
  }
  // end
  std::vector<VertexHandle> cell_vertices{
      getGeomV(bottom1[0]), getGeomV(bottom1[3]), getGeomV(bottom1[2]),
      getGeomV(bottom1[1]), getGeomV(bottom2[2]), getGeomV(bottom2[3]),
      getGeomV(side1[3]),   getGeomV(side2[3])};
  auto _geomch = m_mesh.add_cell(cell_vertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;

  return true;
}
bool MyMesh::AddOneCellCase6(
    const OpenVolumeMesh::CellHandle &_ch,
    const std::vector<OpenVolumeMesh::CellHandle> &_nbc_vec,
    const std::vector<OpenVolumeMesh::HalfFaceHandle> &_nbhf_vec) {
  using namespace OpenVolumeMesh;
#ifndef NDEBUG
  assert(_nbhf_vec.size() == 5 && _nbc_vec.size() == 5);
#endif
  std::vector<VertexHandle> bottom_vec;
  std::vector<std::vector<VertexHandle>> wings;
  // find bottom structure

  {
    HalfFaceHandle bottom_hf;
    std::set<HalfFaceHandle> _nbhf_set(_nbhf_vec.begin(), _nbhf_vec.end());
    for (const auto &x : _nbhf_vec) {
      if (_nbhf_set.count(
              m_topomesh.opposite_halfface_handle_in_cell(x, _ch)) == 0) {
        bottom_vec.assign(m_topomesh.halfface_vertices(x).first,
                          m_topomesh.halfface_vertices(x).second);
        bottom_hf = x;
      }
    }
    auto hes = m_topomesh.halfface(bottom_hf).halfedges();
    auto b_range = m_topomesh.halfface_vertices(bottom_hf);
    std::vector<VertexHandle> vs(b_range.first, b_range.second);
    for (int i = 0; i < 4; ++i) {
      auto _wing_hf = m_topomesh.adjacent_halfface_in_cell(bottom_hf, hes[i]);
      int st = -1;
      auto s_range = m_topomesh.halfface_vertices(_wing_hf);
      std::vector<VertexHandle> _side_vs(s_range.first, s_range.second);
      for (int j = 0; j < 4; ++j) {
        if (_side_vs[j] == vs[(i + 1) % 4]) {
          st = j;
          break;
        }
      }
      wings.push_back({_side_vs[st], _side_vs[(st + 1) % 4],
                       _side_vs[(st + 2) % 4], _side_vs[(st + 3) % 4]});
    }
  }
  // deform
  {
    Vec3d normdir, wingdir[4],
        bottoml
            [4]; //估计的法向量，以及两边正交的方向，bottoml是两个底面交线的方向
    {
      Vec3d v1 = getCoord_topo(bottom_vec[1]) - getCoord_topo(bottom_vec[0]);
      Vec3d v2 = getCoord_topo(bottom_vec[2]) - getCoord_topo(bottom_vec[1]);
      Vec3d v3 = getCoord_topo(bottom_vec[3]) - getCoord_topo(bottom_vec[2]);
      Vec3d v4 = getCoord_topo(bottom_vec[0]) - getCoord_topo(bottom_vec[3]);
      normdir = v1 % v2 + v2 % v3 + v3 % v4 + v4 % v1;

      bottoml[0] = v1, bottoml[1] = v2, bottoml[2] = v3, bottoml[3] = v4;
      wingdir[0] = bottoml[0] % normdir, wingdir[1] = bottoml[1] % normdir,
      wingdir[2] = bottoml[2] % normdir, wingdir[3] = bottoml[3] % normdir;
    }
    // wing的位置调节
    std::map<VertexHandle, Vec3d> fixed;
    {
      std::vector<untangleData> uData(2);
      uData[0] = {normdir,
                  wingdir[0],
                  bottoml[0],
                  getCoord_topo(wings[0][0]),
                  {getCoord_topo(wings[0][2]), getCoord_topo(wings[0][3])}};
      uData[1] = {normdir,
                  wingdir[2],
                  bottoml[2],
                  getCoord_topo(wings[2][0]),
                  {getCoord_topo(wings[2][2]), getCoord_topo(wings[2][3])}};
      std::vector<std::vector<OpenVolumeMesh::Vec3d>> new_pos;
      if (_ch.idx() > 28) {
        new_pos = untangleBottomFace(uData, 2);
      } else {
        new_pos = untangleBottomFace(uData);
      }
      fixed[getGeomV(wings[0][0])] = getCoord_topo(wings[0][0]);
      fixed[getGeomV(wings[0][1])] = getCoord_topo(wings[0][1]);
      fixed[getGeomV(wings[0][2])] = new_pos[0][0];
      fixed[getGeomV(wings[0][3])] = new_pos[0][1];
      fixed[getGeomV(wings[2][0])] = getCoord_topo(wings[2][0]);
      fixed[getGeomV(wings[2][1])] = getCoord_topo(wings[2][1]);
      fixed[getGeomV(wings[2][2])] = new_pos[1][0];
      fixed[getGeomV(wings[2][3])] = new_pos[1][1];
    }
    { ArapOperator::Instance().Deformation(m_mesh, fixed); }
    std::vector<OpenVolumeMesh::VertexHandle> _tmp;
    for (auto &x : fixed) {
      _tmp.push_back(x.first);
    }
    // WriteGeomToVTKFile("tmp.vtk", _tmp);
  }

  // end
  std::vector<VertexHandle> cell_vertices{
      getGeomV(bottom_vec[0]), getGeomV(bottom_vec[3]), getGeomV(bottom_vec[2]),
      getGeomV(bottom_vec[1]), getGeomV(wings[0][2]),   getGeomV(wings[1][2]),
      getGeomV(wings[2][2]),   getGeomV(wings[3][2])};
  auto _geomch = m_mesh.add_cell(cell_vertices);
  m_tm2m_mapping[_ch] = _geomch;
  m_m2tm_mapping[_geomch] = _ch;

  return true;
}

bool MyMesh::Optimize() {
  /*
  auto _h0 = OpenVolumeMesh::VertexHandle(0);
  auto _h1 = OpenVolumeMesh::VertexHandle(1);
  auto _p0 = m_mesh.vertex(_h0);
  auto _p1 = m_mesh.vertex(_h1);
  std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d> fixed;
  // fixed.insert({_h0, _p0});
  // fixed.insert({_h1, _p1});
  ArapOperator::Instance().Optimize(m_mesh, fixed);
  */
  MsqOperator::Instance().Optimize_old(m_mesh);
  return true;
}

#ifndef FINDTOPOSEQ_HPP
#define FINDTOPOSEQ_HPP

#include "ovmwrapper.h"
#include "topology.hpp"
#include <Eigen/Dense>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <fstream>
#include <map>
#include <set>
#include <vector>
using Mesh = OpenVolumeMesh::HexahedralMeshTopologyKernel;

namespace FindTopoSeqPrivate {
std::multimap<int, OpenVolumeMesh::CellHandle>
exposedCells(const Mesh &mesh,
             const std::set<OpenVolumeMesh::CellHandle> &chosen) {
  std::multimap<int, OpenVolumeMesh::CellHandle> res;
  for (auto _ch : mesh.cells()) {
    int innernum = 0;
    if (chosen.count(_ch)) {
      continue;
    }
    for (auto _cch : mesh.cell_cells(_ch)) {
      if (chosen.count(_cch) == 0) {
        innernum += 1;
      }
    }
    res.insert({6 - innernum, _ch});
  }
  return res;
}

std::vector<std::vector<OpenVolumeMesh::CellHandle>>
getTopoSeq_Layers(const Mesh &mesh) {
  // 这里假定没有分离的单元
  std::vector<std::vector<OpenVolumeMesh::CellHandle>> layers;
  std::set<OpenVolumeMesh::CellHandle> chosen;
  while (chosen.size() != mesh.n_cells()) {
    std::set<OpenVolumeMesh::CellHandle> tmp;
    auto ex = exposedCells(mesh, chosen);
    int maxoutnum = ex.rbegin()->first;
    auto range = ex.equal_range(maxoutnum);
    for (auto ch = range.first; ch != range.second; ++ch) {
      chosen.insert(ch->second);
      tmp.insert(ch->second);
    }
    layers.push_back(
        std::vector<OpenVolumeMesh::CellHandle>(tmp.begin(), tmp.end()));
  }

  std::reverse(layers.begin(), layers.end());
  return layers;
}

void dfs(const Mesh &mesh, const OpenVolumeMesh::CellHandle &ch,
         const OpenVolumeMesh::HalfFaceHandle &ori,
         std::set<OpenVolumeMesh::CellHandle> &sheet,
         std::map<OpenVolumeMesh::EdgeHandle, int> &singularties) {
  for (auto _heh : mesh.halfface_halfedges(ori)) {
    auto _nhf = mesh.adjacent_halfface_on_sheet(ori, _heh);
    if (!_nhf.is_valid()) {
      continue;
    }
    auto _ch = mesh.incident_cell(_nhf);
    if (!_ch.is_valid()) {
      continue;
    }
    if (sheet.count(_ch) == 0) {
      sheet.insert(_ch);
      dfs(mesh, _ch, _nhf, sheet, singularties);
    }
  }
  auto _range = mesh.halfface_vertices(ori);
  std::vector<OpenVolumeMesh::VertexHandle> f1(_range.first, _range.second);
  std::vector<OpenVolumeMesh::VertexHandle> f2 =
      opposite_vertex_in_cell(mesh, ch, ori, f1);
  for (int i = 0; i < 4; ++i) {
    auto _v1 = f1[i], _v2 = f2[i];
    for (auto _heit = mesh.voh_iter(_v1); _heit->is_valid(); ++_heit) {
      if (mesh.halfedge(*_heit).to_vertex() == _v2) {
        auto _eh = mesh.edge_handle(*_heit);
        if (!mesh.is_boundary(_eh) && (mesh.valence(_eh) != 3)) {
          singularties[_eh] = mesh.valence(_eh);
        }
        break;
      }
    }
  }
}

static std::pair<std::map<OpenVolumeMesh::EdgeHandle, int>,
                 std::set<OpenVolumeMesh::CellHandle>>
findSheet(const Mesh &mesh, const OpenVolumeMesh::CellHandle &ch,
          const OpenVolumeMesh::HalfFaceHandle &ori) {
  std::set<OpenVolumeMesh::CellHandle> res;
  std::map<OpenVolumeMesh::EdgeHandle, int> singularties;
  dfs(mesh, ch, ori, res, singularties);
  return {singularties, res};
}

std::vector<std::vector<OpenVolumeMesh::CellHandle>>
getTopoSeq_Sheets(const Mesh &mesh) {
  auto cmp = [&](const std::map<OpenVolumeMesh::EdgeHandle, int> &sheet1,
                 const std::map<OpenVolumeMesh::EdgeHandle, int> &sheet2) {
    if (sheet1.size() > sheet2.size()) {
      return true;
    } else if (sheet1.size() < sheet2.size()) {
      return false;
    }
    int s1 = 0, s2 = 0;
    for (auto x : sheet1) {
      s1 += (x.second == 3) ? -1 : 1;
    }
    for (auto x : sheet2) {
      s2 += (x.second == 3) ? -1 : 1;
    }
    return abs(s1) > abs(s2);
  };
  std::map<std::map<OpenVolumeMesh::EdgeHandle, int>,
           std::set<OpenVolumeMesh::CellHandle>, decltype(cmp)>
      sheets(cmp);
  for (auto _ch : mesh.cells()) {
    auto p = findSheet(mesh, _ch, mesh.xback_halfface(_ch));
    sheets.insert(p);
    p = findSheet(mesh, _ch, mesh.yback_halfface(_ch));
    sheets.insert(p);
    p = findSheet(mesh, _ch, mesh.zback_halfface(_ch));
    sheets.insert(p);
  }
  std::vector<std::vector<OpenVolumeMesh::CellHandle>> res;
  for (auto x : sheets) {
    res.push_back(std::vector<OpenVolumeMesh::CellHandle>(x.second.begin(),
                                                          x.second.end()));
  }
  return res;
}

std::set<OpenVolumeMesh::CellHandle>
findOuterRangeCell(const Mesh &mesh,
                   const std::set<OpenVolumeMesh::CellHandle> &part) {
  std::set<OpenVolumeMesh::CellHandle> res;
  for (auto _ch : part) {
    res.insert(_ch);
    for (auto _cch : mesh.cell_cells(_ch)) {
      res.insert(_cch);
    }
  }
  return res;
}

std::vector<OpenVolumeMesh::CellHandle> findTopoSeq(const Mesh &mesh) {
  auto layerseq = getTopoSeq_Layers(mesh);
  auto sheetseq = getTopoSeq_Sheets(mesh);
  std::map<OpenVolumeMesh::CellHandle, int> badsheet;
  for (int i = 0; i < layerseq.size(); ++i) {
    for (auto _ch : layerseq[i]) {
      badsheet[_ch] = i;
    }
  }
  auto cmp = [&](const OpenVolumeMesh::CellHandle &c1,
                 const OpenVolumeMesh::CellHandle &c2) {
    return badsheet[c1] < badsheet[c2];
  };
  for (auto &layer : layerseq) {
    std::sort(layer.begin(), layer.end(), cmp);
  }

  std::vector<OpenVolumeMesh::CellHandle> seq_res;
  std::set<OpenVolumeMesh::CellHandle> chosen;
  seq_res.push_back(layerseq[0][0]);
  chosen.insert(layerseq[0][0]);
  auto outer = findOuterRangeCell(mesh, chosen);
  for (int i = 0; i < layerseq.size(); ++i) {
    while (layerseq[i].size() > 0) {
      auto it = layerseq[i].begin();
      if (chosen.count(*it) != 0) {
        layerseq[i].erase(it);
        continue;
      } else if (outer.count(*it) == 0) {
        while (it != layerseq[i].end() && outer.count(*it) == 0) {
          ++it;
        }
      }
      if (it != layerseq[i].end()) {
        if (chosen.count(*it) == 0) {
          seq_res.push_back(*it);
          chosen.insert(*it);
        }
        layerseq[i].erase(it);
        outer = findOuterRangeCell(mesh, chosen);
        continue;
      }
    }
  }
  return seq_res;
}

// original 17moducore.py
const std::array<std::array<int, 4>, 6> faceposition = {
    std::array<int, 4>{0, 1, 2, 3}, std::array<int, 4>{0, 4, 5, 1},
    std::array<int, 4>{1, 5, 6, 2}, std::array<int, 4>{2, 6, 7, 3},
    std::array<int, 4>{3, 7, 4, 0}, std::array<int, 4>{4, 7, 6, 5}};

const std::array<std::array<int, 8>, 6> cellposition = {
    std::array<int, 8>{0, 1, 2, 3, 4, 5, 6, 7},
    std::array<int, 8>{0, 4, 5, 1, 3, 7, 6, 2},
    std::array<int, 8>{1, 5, 6, 2, 0, 4, 7, 3},
    std::array<int, 8>{2, 6, 7, 3, 1, 5, 4, 0},
    std::array<int, 8>{3, 7, 4, 0, 2, 6, 5, 1},
    std::array<int, 8>{4, 7, 6, 5, 0, 3, 2, 1}};

Topology justifySeq(const Topology &topo) {
  std::vector<int> mapping(topo.v_size, -1);
  std::vector<int> inmapping(topo.v_size, -1);
  int cnt = 0;
  for (auto cell : topo.data) {
    for (auto node : cell) {
      if (mapping[node] == -1) {
        mapping[node] = cnt;
        inmapping[cnt++] = node;
      }
    }
  }
  auto topo_data_cpy = topo.data;
  for (auto &c : topo_data_cpy) {
    std::for_each(c.begin(), c.end(),
                  [&](unsigned char &c) { c = mapping[c]; });
  }
  std::set<std::set<int>> faces;

  for (auto cnum = 0; cnum < topo.c_size; ++cnum) {
    std::vector<std::set<int>> cellfaces = [&]() {
      std::vector<std::set<int>> res;
      res.resize(6);
      for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 4; ++j) {
          res[i].insert(topo_data_cpy[cnum][faceposition[i][j]]);
        }
      }
      return res;
    }();

    for (int i = 0; i < 6; ++i) {
      if (faces.find(cellfaces[i]) != faces.end()) {
        auto tmp = [&]() {
          std::vector<unsigned char> res(8);
          for (int j = 0; j < 8; ++j)
            res[j] = topo_data_cpy[cnum][cellposition[i][j]];
          return res;
        }();
        topo_data_cpy[cnum].swap(tmp);
        break;
      }
    }

    for (int i = 0; i < 6; ++i) {
      faces.insert(cellfaces[i]);
    }
  }
  Topology res(topo);
  res.c_size = topo.c_size, res.v_size = topo.v_size,
  res.data.swap(topo_data_cpy);

  return res;
}

} // namespace FindTopoSeqPrivate

Topology modifyTopoSeq(const Topology &topo) {
  auto mesh = topo.ToOVM();
  auto toposeq = FindTopoSeqPrivate::findTopoSeq(mesh);
  Topology topo_mid;
  {
    topo_mid.c_size = topo.c_size;
    topo_mid.v_size = topo.v_size;
    for (int i = 0; i < toposeq.size(); ++i) {
      topo_mid.data.push_back(topo.data[toposeq[i].idx()]);
    }
  }
  return FindTopoSeqPrivate::justifySeq(topo_mid);
  // auto &ret = FindTopoSeqPrivate::justifySeq(topo_mid);
  // return *ret;
}

/*
int main() {
  std::ifstream fin;
  fin.open("out.vtk");
  Mesh mesh;
  OVMReadHexVtkStream(mesh, fin);
  fin.close();
  auto seq = findTopoSeq(mesh);
  std::ofstream fout;
  fout.open("out1.vtk");
  OVMWriteHexMeshBySeq(mesh, seq, fout);
  fout.close();
  return 0;
}
*/

#endif
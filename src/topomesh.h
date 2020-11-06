/*topomesh的头文件，记录了拓扑结构的网格，用ovm volume mesh 存储，并计算生成顺序*/

#include <OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh>
#include "utils.h"

class TopoMesh {
public:
    OpenVolumeMesh::HexahedralMeshTopologyKernel topomesh;

    //可以用>>操作符
    ErrorCode ReadTopoFromFile(const std::string &filename);
    
};
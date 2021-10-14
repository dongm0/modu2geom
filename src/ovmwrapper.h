#ifndef OVMWRAPPER_H
#define OVMWRAPPER_H

#include "utils.h"
#include <Eigen/Dense>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

using namespace OpenVolumeMesh;

std::vector<VertexHandle>
opposite_vertex_in_cell(const HexahedralMeshTopologyKernel &mesh,
                        const CellHandle &cell, const HalfFaceHandle &hf,
                        const std::vector<OpenVolumeMesh::VertexHandle> &v);

void transform_hex_to_matrix(Eigen::MatrixXd &V, Eigen::MatrixXi &T,
                             Eigen::VectorXi &b, Eigen::MatrixXd &bc,
                             const GeometricHexahedralMeshV3d &mesh,
                             std::map<VertexHandle, Geometry::Vec3d> &fixed);

void transform_hex_to_matrix(Eigen::MatrixXd &V, Eigen::MatrixXi &T,
                             Eigen::VectorXi &b, Eigen::MatrixXd &bc,
                             Eigen::MatrixXi &surface,
                             const GeometricHexahedralMeshV3d &mesh,
                             std::map<VertexHandle, Geometry::Vec3d> &fixed);

void transform_matrix_to_hex(const Eigen::MatrixXd &V,
                             GeometricHexahedralMeshV3d &mesh);

void write_matrix_to_vtk_tet(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T,
                             std::ostream &stream);

template <typename HexMesh>
void OVMReadHexVtkStream(HexMesh &mesh, std::istream &stream,
                         bool _topologyCheck = true,
                         bool _computeBottomUpIncidences = true);

template <typename HexMesh>
OVMWriteHexMesh(HexMesh &mesh, std::ofstream &stream);

#include "owmwrapper.impl"

#endif
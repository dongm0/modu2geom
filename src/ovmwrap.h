
#ifndef OVMWRAP_H
#define OVMWRAP_H

#include "utils.h"
#include <Eigen/Dense>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

std::vector<OpenVolumeMesh::VertexHandle> opposite_vertex_in_cell(
    const OpenVolumeMesh::HexahedralMeshTopologyKernel &mesh,
    const OpenVolumeMesh::CellHandle &cell,
    const OpenVolumeMesh::HalfFaceHandle &hf,
    const std::vector<OpenVolumeMesh::VertexHandle> &v);

void transform_hex_to_matrix(
    Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXi &b,
    Eigen::MatrixXd &bc, const OpenVolumeMesh::GeometricHexahedralMeshV3d &mesh,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        &fixed);

void transform_hex_to_matrix(
    Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::VectorXi &b,
    Eigen::MatrixXd &bc, Eigen::MatrixXi &surface,
    const OpenVolumeMesh::GeometricHexahedralMeshV3d &mesh,
    std::map<OpenVolumeMesh::VertexHandle, OpenVolumeMesh::Geometry::Vec3d>
        &fixed);

void transform_matrix_to_hex(const Eigen::MatrixXd &V,
                             OpenVolumeMesh::GeometricHexahedralMeshV3d &mesh);

void write_matrix_to_vtk_tet(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T,
                             std::ostream &stream);

#include "ovmwrap.cpp"

#endif

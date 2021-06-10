cmake_minimum_required(VERSION 3.10)

project(modu2geom)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

set(IGL_INCLUDE_DIR "/home/dm/github/libigl/include/" CACHE PATH igl)
set(EIGEN_INCLUDE_DIR "/usr/include/eigen3/" CACHE PATH eigen)
set(MESQUITE_INCLUDE_DIR "/home/dm/workspace/trilinos-12.12.1-Source/build/install/include/" CACHE PATH eigen)
set(MESQUITE_LIB_DIR "/home/dm/workspace/trilinos-12.12.1-Source/build/install/lib/" CACHE PATH eigen)
set(OVM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party/include/" CACHE PATH ovm)
set(OVM_LIB_DIR "${CMAKE_SOURCE_DIR}/third-party/lib/" CACHE PATH ovm)
set(TETGEN_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party/include/" CACHE PATH tet)
set(TETGEN_LIB_DIR "${CMAKE_SOURCE_DIR}/third-party/lib/" CACHE PATH tet)
MESSAGE(STATUS ${INCLUDE_DIRECTORIES})

add_subdirectory(./src)



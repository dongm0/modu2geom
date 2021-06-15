cmake_minimum_required(VERSION 3.10)

project(modu2geom)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build.")
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")


if (MSVC_VERSION GREATER_EQUAL "1900")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("/std:c++17" _compiler_supports_cxx17)
    if (_compiler_supports_cxx17)
        add_compile_options("/std:c++17")
    endif()
endif()

set(IGL_INCLUDE_DIR "/home/dm/github/libigl/include/" CACHE PATH igl)
set(EIGEN_INCLUDE_DIR "/usr/include/eigen3/" CACHE PATH eigen)
set(MESQUITE_INCLUDE_DIR "/home/dm/workspace/trilinos-12.12.1-Source/build/install/include/" CACHE PATH eigen)
set(MESQUITE_LIB_DIR "/home/dm/workspace/trilinos-12.12.1-Source/build/install/lib/" CACHE PATH eigen)
set(OVM_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party/include/" CACHE PATH ovm)
set(OVM_LIB_DIR "${CMAKE_SOURCE_DIR}/third-party/lib/" CACHE PATH ovm)
set(TETGEN_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/third-party/include/" CACHE PATH tet)
set(TETGEN_LIB_DIR "${CMAKE_SOURCE_DIR}/third-party/lib/" CACHE PATH tet)
#set(IGL_INCLUDE_DIR "C:\\Users\\Mort_\\Documents\\github\\libigl_x86\\include" CACHE PATH igl)
#set(EIGEN_INCLUDE_DIR "C:\\Users\\Mort_\\Documents\\eigen-3.3.8" CACHE PATH eigen)
#set(MESQUITE_INCLUDE_DIR "C:\\Users\\Mort_\\Documents\\trilinos-12.12.1-Source\\build\\install\\include" CACHE PATH eigen)
#set(MESQUITE_LIB_DIR "C:\\Users\\Mort_\\Documents\\trilinos-12.12.1-Source\\build\\install\\lib" CACHE PATH eigen)
#set(OVM_INCLUDE_DIR "C:\\Users\\Mort_\\Documents\\OpenVolumeMesh-master\\build\\install\\include" CACHE PATH ovm)
#set(OVM_LIB_DIR "C:\\Users\\Mort_\\Documents\\OpenVolumeMesh-master\\build\\install\\lib" CACHE PATH ovm)
#set(TETGEN_INCLUDE_DIR "C:\\Users\\Mort_\\Documents\\github\\tetgen" CACHE PATH tet)
#set(TETGEN_LIB_DIR "C:\\Users\\Mort_\\Documents\\github\\tetgen\\build\\RelWithDebInfo" CACHE PATH tet)
#MESSAGE(STATUS ${INCLUDE_DIRECTORIES})

add_subdirectory(./src)



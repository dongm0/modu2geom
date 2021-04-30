cmake_minimum_required(VERSION 3.10)

project(modu2geom)

add_subdirectory(src)

set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

option(c36 "make tools" OFF)

if (c36)
  add_subdirectory(tools)
endif()


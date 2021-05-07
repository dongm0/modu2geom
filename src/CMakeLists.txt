cmake_minimum_required(VERSION 3.10)

include_directories(${CMAKE_SOURCE_DIR}/third-party/include/ /usr/include/eigen3 /home/dm/github/libigl/include)
link_directories(${CMAKE_SOURCE_DIR}/third-party/lib)

add_executable(modu2geom main.cpp)
add_library(mymesh mymesh.cpp)
add_library(arapoperator arapoperator.cpp)
add_library(ovmwrap ovmwrap.cpp)

target_link_libraries(ovmwrap PUBLIC "OpenVolumeMeshd")
target_link_libraries(arapoperator PUBLIC "OpenVolumeMeshd" PUBLIC "pthread" PRIVATE "tet" PUBLIC ovmwrap)
target_link_libraries(mymesh PUBLIC arapoperator PUBLIC ovmwrap)
target_link_libraries(modu2geom PUBLIC mymesh)
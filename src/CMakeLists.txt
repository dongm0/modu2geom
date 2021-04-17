cmake_minimum_required(VERSION 3.10)

include_directories(${CMAKE_SOURCE_DIR}/third-party/include/ /usr/include/eigen3)
link_directories(${CMAKE_SOURCE_DIR}/third-party/lib)

add_executable(modu2geom main.cpp)
add_library(mymesh mymesh.cpp)
add_library(arapoperator arapoperator.cpp)

target_link_libraries(arapoperator PRIVATE "igl" PRIVATE "pthread")
target_link_libraries(mymesh PUBLIC "OpenVolumeMeshd" PRIVATE arapoperator)
target_link_libraries(modu2geom PRIVATE mymesh)
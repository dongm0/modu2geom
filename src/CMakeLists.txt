cmake_minimum_required(VERSION 3.10)

include_directories(${CMAKE_SOURCE_DIR}/third-party/include/ /usr/include/eigen3)
link_directories(${CMAKE_SOURCE_DIR}/third-party/lib)

add_executable(modu2geom main.cpp)
add_library(mymesh mymesh.cpp)
add_library(arapoperator arapoperator.cpp)
add_library(my_scaf my_scaf.cpp)

target_link_libraries(my_scaf PUBLIC "igl")
target_link_libraries(arapoperator PRIVATE "pthread" PRIVATE my_scaf)
target_link_libraries(mymesh PUBLIC "OpenVolumeMeshd" PRIVATE arapoperator)
target_link_libraries(modu2geom PRIVATE mymesh)
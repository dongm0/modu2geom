INC_DIR = -I/home/dm/github/modu2geom/third-party/include/
INC_DIR += -I/home/dm/github/libigl/include/
INC_DIR += -I/home/dm/github/libigl/external/eigen/
LIB_DIR = -L /home/dm/github/modu2geom/third-party/lib/
STATIC_LIBS = /home/dm/github/libigl/build/libigl.a
LIB_DIR += -lmesquite -lpthread -lOpenVolumeMeshd 

CXX_FLAGS = -g -std=c++11

modu2geom: 
	g++ $(CXX_FLAGS) *.cpp $(STATIC_LIBS) -o modu2geom $(INC_DIR) $(LIB_DIR)
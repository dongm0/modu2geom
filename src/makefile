INC_DIR = -I/home/dm/github/modu2geom/third-party/include/
INC_DIR += -I/home/dm/github/libigl/include/
INC_DIR += -I/home/dm/github/libigl/external/eigen/
LIB_DIR = -L /home/dm/github/modu2geom/third-party/lib/
STATIC_LIBS = /home/dm/github/libigl/build/libigl.a
LIB_DIR += -lmesquite -lpthread -lOpenVolumeMeshd 

CXX_FLAGS = -g -std=c++11

modu2geom: mymesh.o 
	g++ $(CXX_FLAGS) main.cpp mymesh.o $(STATIC_LIBS) -o modu2geom $(INC_DIR) $(LIB_DIR) 

mymesh.o: arapoperator.o msqoptimizer.o
	g++ $(CXX_FLAGS) mymesh.cpp arapoperator.o msqoptimizer.o $(STATIC_LIBS) -c -o mymesh.o $(INC_DIR) $(LIB_DIR)

arapoperator.o:
	g++ $(CXX_FLAGS) arapoperator.cpp $(STATIC_LIBS) -c -o arapoperator.o $(INC_DIR) $(LIB_DIR)

msqoptimizer.o:
	g++ $(CXX_FLAGS) msqoptimizer.cpp -c -o msqoptimizer.o $(INC_DIR) $(LIB_DIR)

clean:
	rm -f modu2geom *.o
WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS =  -Wall
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_DEBUG =  $(INC) -I/home/prusch/gtest-1.6.0/include
CFLAGS_DEBUG =  $(CFLAGS) -Wall -g -fopenmp
RESINC_DEBUG =  $(RESINC)
RCFLAGS_DEBUG =  $(RCFLAGS)
LIBDIR_DEBUG =  $(LIBDIR)
LIB_DEBUG = $(LIB) /home/prusch/gtest-1.6.0/make/gtest_main.a -lgomp
LDFLAGS_DEBUG =  $(LDFLAGS) -lpthread
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/GradientBased

INC_RELEASE =  $(INC) -I/home/prusch/gtest-1.6.0/include
CFLAGS_RELEASE =  $(CFLAGS) -O3 -fopenmp
RESINC_RELEASE =  $(RESINC)
RCFLAGS_RELEASE =  $(RCFLAGS)
LIBDIR_RELEASE =  $(LIBDIR)
LIB_RELEASE = $(LIB) /home/prusch/gtest-1.6.0/make/gtest_main.a -lgomp
LDFLAGS_RELEASE =  $(LDFLAGS) -s -lpthread
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/GradientBased

OBJ_DEBUG = $(OBJDIR_DEBUG)/binaryIO.o $(OBJDIR_DEBUG)/cell.o $(OBJDIR_DEBUG)/lattice.o $(OBJDIR_DEBUG)/main.o $(OBJDIR_DEBUG)/matrix.o $(OBJDIR_DEBUG)/paramset.o $(OBJDIR_DEBUG)/vector.o

OBJ_RELEASE = $(OBJDIR_RELEASE)/binaryIO.o $(OBJDIR_RELEASE)/cell.o $(OBJDIR_RELEASE)/lattice.o $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/matrix.o $(OBJDIR_RELEASE)/paramset.o $(OBJDIR_RELEASE)/vector.o

all: debug release

clean: cleandebug cleanrelease

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/binaryIO.o: binaryIO.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c binaryIO.cpp -o $(OBJDIR_DEBUG)/binaryIO.o

$(OBJDIR_DEBUG)/cell.o: cell.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c cell.cpp -o $(OBJDIR_DEBUG)/cell.o

$(OBJDIR_DEBUG)/lattice.o: lattice.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c lattice.cpp -o $(OBJDIR_DEBUG)/lattice.o

$(OBJDIR_DEBUG)/main.o: main.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c main.cpp -o $(OBJDIR_DEBUG)/main.o

$(OBJDIR_DEBUG)/matrix.o: matrix.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c matrix.cpp -o $(OBJDIR_DEBUG)/matrix.o

$(OBJDIR_DEBUG)/paramset.o: paramset.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c paramset.cpp -o $(OBJDIR_DEBUG)/paramset.o

$(OBJDIR_DEBUG)/vector.o: vector.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c vector.cpp -o $(OBJDIR_DEBUG)/vector.o

cleandebug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/binaryIO.o: binaryIO.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c binaryIO.cpp -o $(OBJDIR_RELEASE)/binaryIO.o

$(OBJDIR_RELEASE)/cell.o: cell.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c cell.cpp -o $(OBJDIR_RELEASE)/cell.o

$(OBJDIR_RELEASE)/lattice.o: lattice.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c lattice.cpp -o $(OBJDIR_RELEASE)/lattice.o

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/matrix.o: matrix.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c matrix.cpp -o $(OBJDIR_RELEASE)/matrix.o

$(OBJDIR_RELEASE)/paramset.o: paramset.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c paramset.cpp -o $(OBJDIR_RELEASE)/paramset.o

$(OBJDIR_RELEASE)/vector.o: vector.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c vector.cpp -o $(OBJDIR_RELEASE)/vector.o

cleanrelease: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug cleandebug before_release after_release cleanrelease


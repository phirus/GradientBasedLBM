WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = -I../../gtest-1.6.0/include
CFLAGS =  -Wall
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_DEBUG =  $(INC) 
CFLAGS_DEBUG =  $(CFLAGS) -Wall -g 
RESINC_DEBUG =  $(RESINC)
RCFLAGS_DEBUG =  $(RCFLAGS)
LIBDIR_DEBUG =  $(LIBDIR)
LIB_DEBUG = $(LIB) ../../gtest-1.6.0/make/gtest_main.a 
LDFLAGS_DEBUG =  $(LDFLAGS) -lpthread
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/GradientBased

OBJ_DEBUG = $(OBJDIR_DEBUG)/binaryIO.o $(OBJDIR_DEBUG)/cell.o $(OBJDIR_DEBUG)/lattice.o $(OBJDIR_DEBUG)/main.o $(OBJDIR_DEBUG)/matrix.o $(OBJDIR_DEBUG)/paramset.o $(OBJDIR_DEBUG)/vector.o


all: debug 
clean: cleandebug

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

.PHONY: before_debug after_debug cleandebug
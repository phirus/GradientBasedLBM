WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = -I../../../gtest-1.6.0/include
CFLAGS =  -Wall -g -fopenmp
LIB = ../../../gtest-1.6.0/make/gtest_main.a
LDFLAGS = -lpthread -lgomp

SRC = ../src
OBJ_DIR = obj
BIN_DIR = bin
OUTPUT = $(BIN_DIR)/GradientBased

# creating OBJs
OBJ_DEBUG = $(OBJ_DIR)/Analyze.o $(OBJ_DIR)/BinaryIO.o $(OBJ_DIR)/Constants2D.o $(OBJ_DIR)/Constants3D.o $(OBJ_DIR)/Cell2D.o $(OBJ_DIR)/Cell3D.o $(OBJ_DIR)/Definitions2D.o $(OBJ_DIR)/Definitions3D.o $(OBJ_DIR)/Lattice2D.o $(OBJ_DIR)/main.o $(OBJ_DIR)/Matrix2D.o $(OBJ_DIR)/Matrix3D.o $(OBJ_DIR)/ParamSet.o $(OBJ_DIR)/Preprocess.o $(OBJ_DIR)/Timetrack.o $(OBJ_DIR)/Vector2D.o $(OBJ_DIR)/Vector3D.o

all: build run

before_debug: 
	test -d $(BIN_DIR) || mkdir -p $(BIN_DIR)
	test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)

after_debug: 

build: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG)
	$(LD) -o $(OUTPUT) $(OBJ_DEBUG)  $(LDFLAGS) $(LIB)

$(OBJ_DIR)/Analyze.o: $(SRC)/Analyze.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Analyze.cc -o $(OBJ_DIR)/Analyze.o

$(OBJ_DIR)/BinaryIO.o: $(SRC)/BinaryIO.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/BinaryIO.cc -o $(OBJ_DIR)/BinaryIO.o

$(OBJ_DIR)/Constants2D.o: $(SRC)/2D/Constants2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Constants2D.cc -o $(OBJ_DIR)/Constants2D.o

$(OBJ_DIR)/Constants3D.o: $(SRC)/3D/Constants3D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/3D/Constants3D.cc -o $(OBJ_DIR)/Constants3D.o

$(OBJ_DIR)/Cell2D.o: $(SRC)/2D/Cell2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Cell2D.cc -o $(OBJ_DIR)/Cell2D.o

$(OBJ_DIR)/Cell3D.o: $(SRC)/3D/Cell3D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/3D/Cell3D.cc -o $(OBJ_DIR)/Cell3D.o

$(OBJ_DIR)/Definitions2D.o: $(SRC)/2D/Definitions2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Definitions2D.cc -o $(OBJ_DIR)/Definitions2D.o	

$(OBJ_DIR)/Definitions3D.o: $(SRC)/3D/Definitions3D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/3D/Definitions3D.cc -o $(OBJ_DIR)/Definitions3D.o	

$(OBJ_DIR)/Lattice2D.o: $(SRC)/2D/Lattice2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Lattice2D.cc -o $(OBJ_DIR)/Lattice2D.o

$(OBJ_DIR)/main.o: main.cc
	$(CXX) $(CFLAGS) $(INC) -c main.cc -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/Matrix2D.o: $(SRC)/2D/Matrix2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Matrix2D.cc -o $(OBJ_DIR)/Matrix2D.o

$(OBJ_DIR)/Matrix3D.o: $(SRC)/3D/Matrix3D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/3D/Matrix3D.cc -o $(OBJ_DIR)/Matrix3D.o

$(OBJ_DIR)/ParamSet.o: $(SRC)/ParamSet.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/ParamSet.cc -o $(OBJ_DIR)/ParamSet.o

$(OBJ_DIR)/Preprocess.o: $(SRC)/Preprocess.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Preprocess.cc -o $(OBJ_DIR)/Preprocess.o

$(OBJ_DIR)/Timetrack.o: $(SRC)/Timetrack.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Timetrack.cc -o $(OBJ_DIR)/Timetrack.o

$(OBJ_DIR)/Vector2D.o: $(SRC)/2D/Vector2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Vector2D.cc -o $(OBJ_DIR)/Vector2D.o

$(OBJ_DIR)/Vector3D.o: $(SRC)/3D/Vector3D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/3D/Vector3D.cc -o $(OBJ_DIR)/Vector3D.o

clean: 
	rm -f $(OBJ_DEBUG) $(OUTPUT)
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
run: 
	./$(OUTPUT)

.PHONY: before_debug after_debug clean
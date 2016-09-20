WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

#INC = -I../../../gtest-1.6.0/include
CFLAGS =  -Wall -O3 -fopenmp -std=c++11
#LIB = gomp
LDFLAGS = -lboost_program_options -lgomp -lpthread

SRC = ../src
OBJ_DIR = obj
BIN_DIR = bin
OUTPUT = BubbleTest

# creating OBJs
OBJ_DEBUG = $(OBJ_DIR)/Analyze_basic.o  $(OBJ_DIR)/Analyze2D.o $(OBJ_DIR)/BasicIO.o $(OBJ_DIR)/BinaryIO2D.o $(OBJ_DIR)/Constants2D.o $(OBJ_DIR)/Cell2D.o $(OBJ_DIR)/Definitions2D.o $(OBJ_DIR)/Lattice2D.o $(OBJ_DIR)/main.o $(OBJ_DIR)/Matrix2D.o $(OBJ_DIR)/ParamSet.o $(OBJ_DIR)/Preprocess.o $(OBJ_DIR)/Timetrack.o $(OBJ_DIR)/Vector2D.o

all: clean build

Release: build

before_debug: 
	test -d $(BIN_DIR) || mkdir -p $(BIN_DIR)
	test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)

after_debug: 

build: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG)
	$(LD) -o $(OUTPUT) $(OBJ_DEBUG)  $(LDFLAGS) $(LIB)

$(OBJ_DIR)/Analyze_basic.o: $(SRC)/Analyze_basic.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Analyze_basic.cc -o $(OBJ_DIR)/Analyze_basic.o

$(OBJ_DIR)/Analyze2D.o: $(SRC)/2D/Analyze2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Analyze2D.cc -o $(OBJ_DIR)/Analyze2D.o

$(OBJ_DIR)/BasicIO.o: $(SRC)/BasicIO.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/BasicIO.cc -o $(OBJ_DIR)/BasicIO.o

$(OBJ_DIR)/BinaryIO2D.o: $(SRC)/2D/BinaryIO2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/BinaryIO2D.cc -o $(OBJ_DIR)/BinaryIO2D.o

$(OBJ_DIR)/Constants2D.o: $(SRC)/2D/Constants2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Constants2D.cc -o $(OBJ_DIR)/Constants2D.o

$(OBJ_DIR)/Cell2D.o: $(SRC)/2D/Cell2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Cell2D.cc -o $(OBJ_DIR)/Cell2D.o

$(OBJ_DIR)/Definitions2D.o: $(SRC)/2D/Definitions2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Definitions2D.cc -o $(OBJ_DIR)/Definitions2D.o	

$(OBJ_DIR)/Lattice2D.o: $(SRC)/2D/Lattice2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Lattice2D.cc -o $(OBJ_DIR)/Lattice2D.o

$(OBJ_DIR)/main.o: main.cc
	$(CXX) $(CFLAGS) $(INC) -c main.cc -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/Matrix2D.o: $(SRC)/2D/Matrix2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Matrix2D.cc -o $(OBJ_DIR)/Matrix2D.o

$(OBJ_DIR)/ParamSet.o: $(SRC)/ParamSet.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/ParamSet.cc -o $(OBJ_DIR)/ParamSet.o

$(OBJ_DIR)/Preprocess.o: $(SRC)/Preprocess.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Preprocess.cc -o $(OBJ_DIR)/Preprocess.o

$(OBJ_DIR)/Timetrack.o: $(SRC)/Timetrack.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Timetrack.cc -o $(OBJ_DIR)/Timetrack.o

$(OBJ_DIR)/Vector2D.o: $(SRC)/2D/Vector2D.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/2D/Vector2D.cc -o $(OBJ_DIR)/Vector2D.o

clean: 
	rm -f $(OBJ_DEBUG) $(OUTPUT)
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
	rm -rf *.vtk
	rm -rf *.csv
	rm -rf *.dat
	rm -rf log
	rm -rf BubbleTest
run: 
	./$(OUTPUT)

.PHONY: before_debug after_debug clean

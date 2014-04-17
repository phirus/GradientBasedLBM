WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

#INC = -I../../../gtest-1.6.0/include
CFLAGS =  -Wall -O3 -fopenmp
#LIB = gomp
LDFLAGS = -lboost_program_options -lgomp
#-lpthread

SRC = ../../src
OBJ_DIR = obj
OUTPUT = translate_params

# creating OBJs
OBJ_DEBUG = $(OBJ_DIR)/BinaryIO.o $(OBJ_DIR)/Constants.o $(OBJ_DIR)/Cell.o $(OBJ_DIR)/Definitions.o $(OBJ_DIR)/Lattice.o $(OBJ_DIR)/main.o $(OBJ_DIR)/Matrix.o $(OBJ_DIR)/ParamSet.o $(OBJ_DIR)/Preprocess.o $(OBJ_DIR)/Timetrack.o $(OBJ_DIR)/Vector.o

all: clean build

Release: build

before_debug: 
	test -d $(BIN_DIR) || mkdir -p $(BIN_DIR)
	test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)

after_debug: 

build: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG)
	$(LD) -o $(OUTPUT) $(OBJ_DEBUG)  $(LDFLAGS) $(LIB)

$(OBJ_DIR)/BinaryIO.o: $(SRC)/BinaryIO.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/BinaryIO.cc -o $(OBJ_DIR)/BinaryIO.o

$(OBJ_DIR)/Constants.o: $(SRC)/Constants.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Constants.cc -o $(OBJ_DIR)/Constants.o

$(OBJ_DIR)/Cell.o: $(SRC)/Cell.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Cell.cc -o $(OBJ_DIR)/Cell.o

$(OBJ_DIR)/Definitions.o: $(SRC)/Definitions.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Definitions.cc -o $(OBJ_DIR)/Definitions.o

$(OBJ_DIR)/Lattice.o: $(SRC)/Lattice.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Lattice.cc -o $(OBJ_DIR)/Lattice.o

$(OBJ_DIR)/main.o: main.cc
	$(CXX) $(CFLAGS) $(INC) -c main.cc -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/Matrix.o: $(SRC)/Matrix.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Matrix.cc -o $(OBJ_DIR)/Matrix.o

$(OBJ_DIR)/ParamSet.o: $(SRC)/ParamSet.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/ParamSet.cc -o $(OBJ_DIR)/ParamSet.o

$(OBJ_DIR)/Preprocess.o: $(SRC)/Preprocess.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Preprocess.cc -o $(OBJ_DIR)/Preprocess.o

$(OBJ_DIR)/Timetrack.o: $(SRC)/Timetrack.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Timetrack.cc -o $(OBJ_DIR)/Timetrack.o

$(OBJ_DIR)/Vector.o: $(SRC)/Vector.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Vector.cc -o $(OBJ_DIR)/Vector.o

clean: 
	rm -f $(OBJ_DEBUG) $(OUTPUT)
	rm -rf $(OBJ_DIR)
run: 
	./$(OUTPUT)

.PHONY: before_debug after_debug clean
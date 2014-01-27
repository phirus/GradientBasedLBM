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
OBJ_DEBUG = $(OBJ_DIR)/binaryIO.o $(OBJ_DIR)/constants.o $(OBJ_DIR)/Cell.o $(OBJ_DIR)/Lattice.o $(OBJ_DIR)/main.o $(OBJ_DIR)/matrix.o $(OBJ_DIR)/paramset.o $(OBJ_DIR)/preprocess.o $(OBJ_DIR)/timetrack.o $(OBJ_DIR)/Vector.o

all: build run

before_debug: 
	test -d $(BIN_DIR) || mkdir -p $(BIN_DIR)
	test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)

after_debug: 

build: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG)
	$(LD) -o $(OUTPUT) $(OBJ_DEBUG)  $(LDFLAGS) $(LIB)

$(OBJ_DIR)/binaryIO.o: $(SRC)/binaryIO.cpp
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/binaryIO.cpp -o $(OBJ_DIR)/binaryIO.o

$(OBJ_DIR)/constants.o: $(SRC)/constants.cpp
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/constants.cpp -o $(OBJ_DIR)/constants.o

$(OBJ_DIR)/Cell.o: $(SRC)/Cell.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Cell.cc -o $(OBJ_DIR)/Cell.o

$(OBJ_DIR)/Lattice.o: $(SRC)/Lattice.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Lattice.cc -o $(OBJ_DIR)/Lattice.o

$(OBJ_DIR)/main.o: main.cpp
	$(CXX) $(CFLAGS) $(INC) -c main.cpp -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/matrix.o: $(SRC)/matrix.cpp
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/matrix.cpp -o $(OBJ_DIR)/matrix.o

$(OBJ_DIR)/paramset.o: $(SRC)/paramset.cpp
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/paramset.cpp -o $(OBJ_DIR)/paramset.o

$(OBJ_DIR)/preprocess.o: $(SRC)/preprocess.cpp
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/preprocess.cpp -o $(OBJ_DIR)/preprocess.o

$(OBJ_DIR)/timetrack.o: $(SRC)/timetrack.cpp
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/timetrack.cpp -o $(OBJ_DIR)/timetrack.o

$(OBJ_DIR)/Vector.o: $(SRC)/Vector.cc
	$(CXX) $(CFLAGS) $(INC) -c $(SRC)/Vector.cc -o $(OBJ_DIR)/Vector.o

clean: 
	rm -f $(OBJ_DEBUG) $(OUTPUT)
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
run: 
	./$(OUTPUT)

.PHONY: before_debug after_debug clean
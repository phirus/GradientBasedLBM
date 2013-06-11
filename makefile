WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = -I../../gtest-1.6.0/include
CFLAGS =  -Wall -g
LIB = ../../gtest-1.6.0/make/gtest_main.a
LDFLAGS = -lpthread

OBJ_DIR = obj/
BIN_DIR = bin/
OUTPUT = $(BIN_DIR)/GradientBased

# creating OBJs

OBJ_DEBUG = $(OBJ_DIR)/binaryIO.o $(OBJ_DIR)/cell.o $(OBJ_DIR)/lattice.o $(OBJ_DIR)/main.o $(OBJ_DIR)/matrix.o $(OBJ_DIR)/paramset.o $(OBJ_DIR)/vector.o


all: debug run

before_debug: 
	test -d $(BIN_DIR) || mkdir -p $(BIN_DIR)
	test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG)
	$(LD) -o $(OUTPUT) $(OBJ_DEBUG)  $(LDFLAGS) $(LIB)

$(OBJ_DIR)/binaryIO.o: binaryIO.cpp
	$(CXX) $(CFLAGS) $(INC) -c binaryIO.cpp -o $(OBJ_DIR)/binaryIO.o

$(OBJ_DIR)/cell.o: cell.cpp
	$(CXX) $(CFLAGS) $(INC) -c cell.cpp -o $(OBJ_DIR)/cell.o

$(OBJ_DIR)/lattice.o: lattice.cpp
	$(CXX) $(CFLAGS) $(INC) -c lattice.cpp -o $(OBJ_DIR)/lattice.o

$(OBJ_DIR)/main.o: main.cpp
	$(CXX) $(CFLAGS) $(INC) -c main.cpp -o $(OBJ_DIR)/main.o

$(OBJ_DIR)/matrix.o: matrix.cpp
	$(CXX) $(CFLAGS) $(INC) -c matrix.cpp -o $(OBJ_DIR)/matrix.o

$(OBJ_DIR)/paramset.o: paramset.cpp
	$(CXX) $(CFLAGS) $(INC) -c paramset.cpp -o $(OBJ_DIR)/paramset.o

$(OBJ_DIR)/vector.o: vector.cpp
	$(CXX) $(CFLAGS) $(INC) -c vector.cpp -o $(OBJ_DIR)/vector.o

clean: 
	rm -f $(OBJ_DEBUG) $(OUTPUT)
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
run: 
	./$(OUTPUT)

.PHONY: before_debug after_debug clean
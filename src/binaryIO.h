#ifndef BINARYIO_H
#define BINARYIO_H

#include<fstream>
#include<iostream>
#include"lattice.h"

void binary_output(const Lattice& l, const string& filename = "data.bin");
const bool binary_input(Lattice& outL, const string& filename = "data.bin");

void techplotOutput(const Lattice& l, int iterNum, bool vebose = false);
void vtkOutput(const Lattice& l, int iterNum);

#endif


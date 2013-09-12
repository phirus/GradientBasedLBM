#ifndef BINARYIO_H
#define BINARYIO_H

#include<fstream>
#include<iostream>
#include<map>

#include"lattice.h"
#include"preprocess.h"

void binary_output(const Lattice& l, const string& filename = "data.bin");
const bool binary_input(Lattice& outL, const string& filename = "data.bin");

void restart_file(const Lattice& l, Preprocess p, const string& filename = "restart.bin");
const bool restart_read(Lattice& l, Preprocess p, const string& filename = "restart.bin");

void techplotOutput(const Lattice& l, int iterNum, bool vebose = false);
void vtkOutput(const Lattice& l, int iterNum);

void paramLogOut(const Lattice& l);

const bool inputQuery(const string& filename, const string& query, double& value);

const ParamSet getFileParams(const string& filename);
const Preprocess getFilePreprocess(const string& filename);

#endif


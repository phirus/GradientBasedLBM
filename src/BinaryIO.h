#ifndef BINARYIO_H
#define BINARYIO_H

#include<fstream>
#include<iostream>
#include<map>

#include"Lattice.h"
#include"Preprocess.h"

void binary_output(const Lattice& l, const string& filename = "data.bin");
const bool binary_input(Lattice& outL, const string& filename = "data.bin");

void restart_file(const Lattice& l, const Preprocess& p, const Timetrack time, const string& filename = "restart.bin");
const bool restart_read(Lattice& outL, Preprocess& p, Timetrack& time, const string& filename = "restart.bin");

void techplotOutput(const Lattice& l, int iterNum, bool vebose = false);
void vtkOutput(const Lattice& l, int iterNum);

void paramLogOut(const Lattice& l);

const bool inputQuery(const string& filename, const string& query, double& value);

const Preprocess getFilePreprocess(const string& filename);
const Timetrack getFileTimetrack(const Preprocess& prepro, const string& filename);

#endif


#ifndef BINARYIO2D_H
#define BINARYIO2D_H

#include<fstream>
#include<iostream>
#include<map>

#include"Lattice2D.h"

#include"../BasicIO.h"

/// binary dump
void write_binary(const Lattice2D& l, const string& filename = "data.bin");
const bool read_binary(Lattice2D& outL, const string& filename = "data.bin");

/// restart files
void write_restart_file(const Lattice2D& l, const Preprocess& p, const Timetrack time, const string& filename = "restart.bin");
const bool read_restart_file(Lattice2D& outL, Preprocess& p, Timetrack& time, const string& filename = "restart.bin");

/// write output
void write_techplot_output(const Lattice2D& l, int iterNum);
void write_techplot_output_alternative(const Lattice2D& l, const string& filename = "alternative.dat");

void write_vtk_output(const Lattice2D& l, const string& filename = "test.vtk");
inline void write_vtk_output(const Lattice2D& l, int iterNum){write_vtk_output(l, createFilename("output_", iterNum, ".vtk"));};


#endif
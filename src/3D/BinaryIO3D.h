#ifndef BINARYIO3D_H
#define BINARYIO3D_H

#include<fstream>
#include<iostream>
#include<map>

#include"Lattice3D.h"

#include"../BasicIO.h"

/// binary dump
void write_binary3D(const Lattice3D& l, const string& filename = "data.bin");
const bool read_binary3D(Lattice3D& outL, const string& filename = "data.bin");

/// restart files
void write_restart_file3D(const Lattice3D& l, const Preprocess& p, const Timetrack time, const string& filename = "restart.bin");
const bool read_restart_file3D(Lattice3D& outL, Preprocess& p, Timetrack& time, const string& filename = "restart.bin");

/// write output
void write_techplot_output3D(const Lattice3D& l, int iterNum);
void write_techplot_output_alternative3D(const Lattice3D& l, const string& filename = "alternative.dat");

void write_vtk_output3D(const Lattice3D& l, const string& filename = "test.vtk");
inline void write_vtk_output3D(const Lattice3D& l, int iterNum){write_vtk_output3D(l, createFilename("output_", iterNum, ".vtk"));};

// void writeBubbleFitData(const Lattice3D& l, const string& filename = "bubbleFit.csv");
#endif
#ifndef BINARYIO2D_H
#define BINARYIO2D_H

#include<fstream>
#include<iostream>
#include<map>

#include"Lattice2D.h"
#include"../Preprocess.h"
#include"../Timetrack.h"

/// binary dump
void write_binary(const Lattice2D& l, const string& filename = "data.bin");
const bool read_binary(Lattice2D& outL, const string& filename = "data.bin");

/// restart files
void write_restart_file(const Lattice2D& l, const Preprocess& p, const Timetrack time, const string& filename = "restart.bin");
const bool read_restart_file(Lattice2D& outL, Preprocess& p, Timetrack& time, const string& filename = "restart.bin");

///auxiliary
const string createFilename(const string& name, int iteration, const string& type = ".bin");

/// write output
void write_techplot_output(const Lattice2D& l, int iterNum);
void write_techplot_output_alternative(const Lattice2D& l, const string& filename = "alternative.dat");

void write_vtk_output(const Lattice2D& l, const string& filename = "test.vtk");
inline void write_vtk_output(const Lattice2D& l, int iterNum){write_vtk_output(l, createFilename("output_", iterNum, ".vtk"));};

void write_data_plot(const std::vector<double> x, const std::vector<double> y1, const std::vector<double> y2, const string& filename = "massplot.dat");
void write_data_plot(const std::vector<double> y, double del_x, const string& filename = "ReynoldsPlot.dat");

void write_param_log(const ParamSet& p);

/// read input
const Preprocess read_preprocess_file(const string& filename);
const Timetrack read_timetrack_file(const string& filename);
const ParamSet read_paramset_file(const string& filename = "paramLog");

const bool input_query(const string& filename, const string& query, double& value);

#endif
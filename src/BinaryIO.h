#ifndef BINARYIO_H
#define BINARYIO_H

#include<fstream>
#include<iostream>
#include<map>

#include"Lattice.h"
#include"Preprocess.h"

void write_binary(const Lattice& l, const string& filename = "data.bin");
const bool read_binary(Lattice& outL, const string& filename = "data.bin");

void write_restart_file(const Lattice& l, const Preprocess& p, const Timetrack time, const string& filename = "restart.bin");
const bool read_restart_file(Lattice& outL, Preprocess& p, Timetrack& time, const string& filename = "restart.bin");

void write_techplot_output(const Lattice& l, int iterNum, bool vebose = false);
void write_vtk_output(const Lattice& l, int iterNum);

void write_param_log(const ParamSet& p);

const bool input_query(const string& filename, const string& query, double& value);

const Preprocess read_preprocess_file(const string& filename);
const Timetrack read_timetrack_file(const Preprocess& prepro, const string& filename);

void write_data_plot(const std::vector<double> x, const std::vector<double> y1, const std::vector<double> y2, const string& filename = "massplot.dat");

#endif
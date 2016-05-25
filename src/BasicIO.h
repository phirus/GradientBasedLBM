#ifndef BASICIO_H
#define BASICIO_H

#include<fstream>
#include<iostream>
#include<sstream>
#include<map>
#include<string.h>
#include<vector>

#include"Preprocess.h"
#include"Timetrack.h"

using namespace std;
typedef std::vector<std::vector<double>> nested_vector; 
/// write output
void write_data_plot(const std::vector<double> x, const std::vector<double> y1, const std::vector<double> y2, const string& filename = "massplot.dat");
void write_data_plot(const std::vector<double> y, double del_x, const string& filename = "ReynoldsPlot.dat");
void write_data_plot(const std::vector<double> y1, const std::vector<double> y2, double del_x, const string& filename = "BubbleVeloPlot.dat");
void write_csv(const nested_vector& data, const string& filename = "BubbleVeloPlot.csv", const string& header= "time , Posx , PosY , v_y , v_x , Re");
void write_param_log(const ParamSet& p);

/// read input
const bool input_query(const string& filename, const string& query, double& value);
const Preprocess read_preprocess_file(const string& filename);
const Timetrack read_timetrack_file(const string& filename);
const ParamSet read_paramset_file(const string& filename = "paramLog");

///auxiliary
const string createFilename(const string& name, int iteration, const string& type = ".bin");

#endif
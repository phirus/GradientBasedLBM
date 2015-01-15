#ifndef ANALYZE_BASIC_H
#define ANALYZE_BASIC_H

#include<cmath>

#include"ParamSet.h"

using namespace std;

const double getEotvos(const ParamSet& params, double resolution = 40);
const double getMorton(const ParamSet& params);
const double getReynolds(const ParamSet& params, double velocity, double resolution = 40);

#endif
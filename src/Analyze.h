#ifndef ANALYZE_H
#define ANALYZE_H

#include<cmath>

#include"Lattice.h"

using namespace std;

const double getEotvos(const ParamSet& params, double resolution = 40);
const double getMorton(const ParamSet& params);
const double getReynolds(const ParamSet& params, double velocity, double resolution = 40);
const double getBubbleVelocity(Lattice l);

#endif
#ifndef ANALYZE_H
#define ANALYZE_H

#include<cmath>

#include"2D/Lattice2D.h"
#include"3D/Lattice3D.h"

using namespace std;

const double getBubbleVelocity(const Lattice2D& l);
const double getBubbleVelocity(const Lattice3D& l);

const double getEotvos(const ParamSet& params, double resolution = 40);
const double getMorton(const ParamSet& params);
const double getReynolds(const ParamSet& params, double velocity, double resolution = 40);
const double getReynolds(const Lattice2D& l, double resolution = 40);
const double getReynolds(const Lattice3D& l, double resolution = 40);

#endif
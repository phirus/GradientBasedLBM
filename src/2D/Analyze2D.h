#ifndef ANALYZE2D_H
#define ANALYZE2D_H

#include<cmath>
#include<vector>

#include"Lattice2D.h"
#include"../Analyze_basic.h"

using namespace std;

const double getBubbleVelocity(const Lattice2D& l);
const double getReynolds(const Lattice2D& l, double resolution = 40);

const vector<double> getShearResidual(const Lattice2D& l);
#endif
#ifndef ANALYZE3D_H
#define ANALYZE3D_H

#include<cmath>

#include"Lattice3D.h"
#include"../Analyze_basic.h"

using namespace std;

const double getBubbleVelocity(const Lattice3D& l);
const double getReynolds(const Lattice3D& l, double resolution = 40);

#endif
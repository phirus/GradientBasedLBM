#ifndef ANALYZE3D_H
#define ANALYZE3D_H

#include<cmath>

#include"Lattice3D.h"
#include"../Analyze_basic.h"

using namespace std;

const Vector3D getBubbleVelocity(const Lattice3D& l);
const double getReynolds(const Lattice3D& l, double resolution = 40);
const double getLineShearSum(const Lattice3D& l);

#endif
#ifndef ANALYZE2D_H
#define ANALYZE2D_H

#include<cmath>

#include"Lattice2D.h"
#include"../Analyze_basic.h"

using namespace std;

const Vector2D getBubbleVelocity(const Lattice2D& l);
const Vector2D getBubblePosition(const Lattice2D& l);
const double getReynolds(const Lattice2D& l, double resolution = 40);
const double getLineShearSum(const Lattice2D& l);

#endif
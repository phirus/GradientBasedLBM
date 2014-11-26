#ifndef CONSTANTSBASIC_H
#define CONSTANTSBASIC_H

#include<boost/array.hpp>

//=========================== CONSTANTS ===========================

/// Pi
const double PI = 3.14159265358979323846264338327950288419716939937510;

/// arbitrary  definition
const double MACH_MAX = 0.1; // maximal erlaubte Mach-Zahl

//=========================== TYPES ===========================

typedef boost::array<double,2> ColSet;         /// < simple 2d vector

//=========================== STRUCTS ===========================

/// structure for interpolation paramters
struct Interpol
{
    double chi, eta, kappa, lambda, ny;
};

//=========================== FUNCTIONS ===========================

/// sums up all elements
inline const double sum(const ColSet& vector){return vector[0] + vector[1] ;};

#endif
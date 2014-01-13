#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include<boost/array.hpp>
#include<boost/multi_array.hpp>
#include"vector.h"

/// contains custom typedefs
typedef boost::array<double,9> array;       /// < used to describe single distributions
typedef boost::array<array,2> FSet;         /// < merge two distributions into a single variable
typedef boost::array<Vector,13> direction ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<double,2> ColSet;         /// < simple 2d vector, y = vec[0], x = vec[1]

/// sums up all elements
inline const double sum(const ColSet& vector){return vector[0] + vector[1] ;};

struct RelaxationPar
{
    double s_2,s_3,s_5;
};

#endif

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include<boost/array.hpp>
#include<boost/multi_array.hpp>
#include"Vector.h"

/// contains custom typedefs
typedef boost::array<double,9> array;       /// < used to describe single distributions
typedef boost::array<array,2> DistributionSetType;         /// < merge two distributions into a single variable
typedef boost::array<Vector,13> direction ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<double,2> ColSet;         /// < simple 2d vector, y = vec[0], x = vec[1]

/// sums up all elements
inline const double sum(const ColSet& vector){return vector[0] + vector[1] ;};

struct RelaxationPar
{
    double s_2,s_3,s_5;
    RelaxationPar(double s2=1, double s3=1, double s5=1):s_2(s2),s_3(s3),s_5(s5){};
};

/// structure for interpolation paramters
struct Interpol
{
    double chi, eta, kappa, lambda, ny;
};

//=========================== FUNCTIONS ===========================

/// functions handling basic operations on arrays
const array array_diff(const array &one, const array &two);
const array array_add(const array &one, const array &two);
const array array_times(const array &foo, double factor);

/// translating array functions to distribution sets 
inline const DistributionSetType distro_diff(const DistributionSetType &one, const DistributionSetType &two)
{
    const DistributionSetType diff = {{array_diff(one[0],two[0]), array_diff(one[1],two[1])}};
    return diff;
};

inline const DistributionSetType distro_add(const DistributionSetType &one, const DistributionSetType &two)
{
    const DistributionSetType foo = {{array_add(one[0],two[0]),array_add(one[1],two[1])}};
    return foo;
};

inline const DistributionSetType distro_add_array(const DistributionSetType &one, const array &two)
{
    const DistributionSetType foo = {{array_add(one[0],two), array_add(one[1],two)}};
    return foo;
};

inline const DistributionSetType distro_times(const DistributionSetType &one, double factor)
{
    const DistributionSetType foo = {{array_times(one[0],factor), array_times(one[1],factor)}};    
    return foo;
};

#endif

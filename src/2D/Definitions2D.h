#ifndef DEFINITIONS2D_H
#define DEFINITIONS2D_H

#include<boost/multi_array.hpp>
#include"Vector2D.h"
#include"../Constants_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<double,9> array2D;       /// < used to describe single distributions
typedef boost::array<array2D,2> DistributionSetType2D;         /// < merge two distributions into a single variable
typedef boost::array<Vector2D,13> direction2D ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<Vector2D,2> VeloSet2D;

struct RelaxationPar2D
{
    double s_2,s_3,s_5;
    RelaxationPar2D(double s2=1, double s3=1, double s5=1):s_2(s2),s_3(s3),s_5(s5){};
};

/// structure for interpolation paramters
struct Interpol
{
    double chi, eta, kappa, lambda, ny;
};

//=========================== FUNCTIONS ===========================

/// functions handling basic operations on arrays
const array2D array_diff_2D(const array2D &one, const array2D &two);
const array2D array_add_2D(const array2D &one, const array2D &two);
const array2D array_times_2D(const array2D &foo, double factor);

/// translating array functions to distribution sets 
inline const DistributionSetType2D distro_diff_2D(const DistributionSetType2D &one, const DistributionSetType2D &two)
{
    const DistributionSetType2D diff = {{array_diff_2D(one[0],two[0]), array_diff_2D(one[1],two[1])}};
    return diff;
};

inline const DistributionSetType2D distro_add_2D(const DistributionSetType2D &one, const DistributionSetType2D &two)
{
    const DistributionSetType2D foo = {{array_add_2D(one[0],two[0]),array_add_2D(one[1],two[1])}};
    return foo;
};

inline const DistributionSetType2D distro_add_array_2D(const DistributionSetType2D &one, const array2D &two)
{
    const DistributionSetType2D foo = {{array_add_2D(one[0],two), array_add_2D(one[1],two)}};
    return foo;
};

inline const DistributionSetType2D distro_times_2D(const DistributionSetType2D &one, double factor)
{
    const DistributionSetType2D foo = {{array_times_2D(one[0],factor), array_times_2D(one[1],factor)}};    
    return foo;
};

#endif
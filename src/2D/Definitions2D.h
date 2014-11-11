#ifndef DEFINITIONS2D_H
#define DEFINITIONS2D_H

#include<boost/array.hpp>
#include<boost/multi_array.hpp>
#include"Vector2D.h"
#include"../Constants_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<double,9> array2D;       /// < used to describe single distributions
typedef boost::array<array2D,2> DistributionSetType2D;         /// < merge two distributions into a single variable
typedef boost::array<Vector2D,13> direction2D ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<Vector2D,2> VeloSet2D;

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
const array2D array2D_diff(const array2D &one, const array2D &two);
const array2D array2D_add(const array2D &one, const array2D &two);
const array2D array2D_times(const array2D &foo, double factor);

/// translating array functions to distribution sets 
inline const DistributionSetType2D distro_diff(const DistributionSetType2D &one, const DistributionSetType2D &two)
{
    const DistributionSetType2D diff = {{array2D_diff(one[0],two[0]), array2D_diff(one[1],two[1])}};
    return diff;
};

inline const DistributionSetType2D distro_add(const DistributionSetType2D &one, const DistributionSetType2D &two)
{
    const DistributionSetType2D foo = {{array2D_add(one[0],two[0]),array2D_add(one[1],two[1])}};
    return foo;
};

inline const DistributionSetType2D distro_add_array2D(const DistributionSetType2D &one, const array2D &two)
{
    const DistributionSetType2D foo = {{array2D_add(one[0],two), array2D_add(one[1],two)}};
    return foo;
};

inline const DistributionSetType2D distro_times(const DistributionSetType2D &one, double factor)
{
    const DistributionSetType2D foo = {{array2D_times(one[0],factor), array2D_times(one[1],factor)}};    
    return foo;
};

#endif
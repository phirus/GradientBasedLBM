#ifndef DEFINITIONS3D_H
#define DEFINITIONS3D_H

#include<boost/multi_array.hpp>
#include"Vector3D.h"
#include"../Constants_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<double,19> array3D;       /// < used to describe single distributions
typedef boost::array<array3D,2> DistributionSetType3D;         /// < merge two distributions into a single variable
typedef boost::array<Vector3D,19> direction3D ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<Vector3D,2> VeloSet3D;

// struct RelaxationPar
// {
//     double s_2,s_3,s_5;
//     RelaxationPar(double s2=1, double s3=1, double s5=1):s_2(s2),s_3(s3),s_5(s5){};
// };

// /// structure for interpolation paramters
// struct Interpol
// {
//     double chi, eta, kappa, lambda, ny;
// };

// //=========================== FUNCTIONS ===========================

/// functions handling basic operations on arrays
const array3D array_diff_3D(const array3D &one, const array3D &two);
const array3D array_add_3D(const array3D &one, const array3D &two);
const array3D array_times_3D(const array3D &foo, double factor);

/// translating array functions to distribution sets 
inline const DistributionSetType3D distro_diff_3D(const DistributionSetType3D &one, const DistributionSetType3D &two)
{
    const DistributionSetType3D diff = {{array_diff_3D(one[0],two[0]), array_diff_3D(one[1],two[1])}};
    return diff;
};

inline const DistributionSetType3D distro_add_3D(const DistributionSetType3D &one, const DistributionSetType3D &two)
{
    const DistributionSetType3D foo = {{array_add_3D(one[0],two[0]),array_add_3D(one[1],two[1])}};
    return foo;
};

inline const DistributionSetType3D distro_add_array_3D(const DistributionSetType3D &one, const array3D &two)
{
    const DistributionSetType3D foo = {{array_add_3D(one[0],two), array_add_3D(one[1],two)}};
    return foo;
};

inline const DistributionSetType3D distro_times_3D(const DistributionSetType3D &one, double factor)
{
    const DistributionSetType3D foo = {{array_times_3D(one[0],factor), array_times_3D(one[1],factor)}};    
    return foo;
};

#endif
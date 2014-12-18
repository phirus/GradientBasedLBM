#ifndef DEFINITIONS3D_H
#define DEFINITIONS3D_H

#include<boost/multi_array.hpp>
#include"Vector3D.h"
#include"../Constants_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<double,19> array3D;       /// < used to describe single distributions
typedef boost::array<double,3> DimSet3D;         /// < simple 3d vector
typedef boost::array<array3D,2> DistributionSetType3D;         /// < merge two distributions into a single variable
typedef boost::array<Vector3D,33> direction3D ; /// < collection of 33 direction vectors (D2Q33)
typedef boost::array<Vector3D,2> VeloSet3D;

struct RelaxationPar3D
{
    double s_2,s_3,s_5,s_11,s_17;
    RelaxationPar3D(double s2=1, double s3=1, double s5=1, double s11=1, double s17=1):s_2(s2),s_3(s3),s_5(s5),s_11(s11),s_17(s17){};
};

//=========================== FUNCTIONS ===========================

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
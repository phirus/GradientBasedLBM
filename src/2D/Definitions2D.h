#ifndef DEFINITIONS2D_H
#define DEFINITIONS2D_H

#include"Vector2D.h"
#include"../Definitions_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<Vector2D,13> direction2D ; /// < collection of 13 direction vectors (D2Q13)
typedef boost::array<Vector2D,2> VeloSet2D;
typedef boost::array<int,2> SizeSet2D;

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

struct pressureTensor2D 
{
    double xx,xy,yx,yy;
    pressureTensor2D() : xx(0),xy(0),yx(0),yy(0) {};
    const double getTrace()const{return xx + yy;};
};

#endif
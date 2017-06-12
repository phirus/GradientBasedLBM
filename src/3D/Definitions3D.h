#ifndef DEFINITIONS3D_H
#define DEFINITIONS3D_H

#include"Vector3D.h"
#include"../Definitions_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<Vector3D,33> direction3D ; /// < collection of 33 direction vectors (D2Q33)
typedef boost::array<Vector3D,2> VeloSet3D;
typedef boost::array<int,3> SizeSet3D;

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

struct pressureTensor3D 
{
    double xx,xy,xz,yx,yy,yz,zx,zy,zz;
    pressureTensor3D() : xx(0), xy(0), xz(0), yx(0), yy(0), yz(0), zx(0), zy(0), zz(0) {};
    const double getTrace()const{return xx + yy + zz;};
};
#endif
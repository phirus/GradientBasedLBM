#ifndef DEFINITIONS_BASIC_H
#define DEFINITIONS_BASIC_H

#include<boost/multi_array.hpp>
#include"Constants_basic.h"

/// contains custom typedefs
//=========================== TYPES ===========================

typedef boost::array<double,9> array2D;       /// < used to describe single distributions
typedef boost::array<array2D,2> DistributionSetType2D;         /// < merge two distributions into a single variable

typedef boost::array<double,19> array3D;       /// < used to describe single distributions
typedef boost::array<double,3> DimSet3D;         /// < simple 3d vector
typedef boost::array<array3D,2> DistributionSetType3D;         /// < merge two distributions into a single variable

struct RelaxationPar2D
{
    double s_2,s_3,s_5;
    RelaxationPar2D(double s2=1, double s3=1, double s5=1):s_2(s2),s_3(s3),s_5(s5){};
};

struct RelaxationPar3D
{
    double s_2,s_3,s_5,s_11,s_17;
    RelaxationPar3D(double s2=1, double s3=1, double s5=1, double s11=1, double s17=1):s_2(s2),s_3(s3),s_5(s5),s_11(s11),s_17(s17){};
};

#endif
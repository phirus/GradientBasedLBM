/// organize the available information on boundary conditions

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include"3D/Definitions3D.h"

using namespace std;

enum boundary_type {periodic, bounceback, pressure};

class BoundaryInformation
{
public: 
    /// Lifecycle
    BoundaryInformation(boundary_type bt = periodic, ColSet density = {{0,0}}, VeloSet3D velocity = {{Vector3D(),Vector3D()}}):type(bt),rho(density),u(velocity){};
    BoundaryInformation(int bt, ColSet density = {{0,0}}, VeloSet3D velocity = {{Vector3D(),Vector3D()}}):type(static_cast<boundary_type>(bt)),rho(density),u(velocity){};
    
    /// accessors
    inline const boundary_type getType()const{return type;};
    inline const ColSet getRho()const{return rho;};
    inline const VeloSet3D getVelocity()const{return u;};

private:
    boundary_type type;
    ColSet rho;                    
    VeloSet3D u;
};

struct Boundaries
{
    BoundaryInformation north;
    BoundaryInformation east;
    BoundaryInformation south;
    BoundaryInformation west;
    BoundaryInformation front;
    BoundaryInformation back;
};

#endif 
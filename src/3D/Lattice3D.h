#ifndef LATTICE3D_H
#define LATTICE3D_H
#include<iostream>
#include<omp.h>

#include"Cell3D.h"
#include"../ParamSet.h"
// #include"Timetrack.h"

/// custom typedef for the whole field of cells
typedef boost::multi_array<Cell3D,3> field3D;

/// contains the domain of the simulation with LB operators
class Lattice3D
{
public:
    /// Lifecycle
    Lattice3D(int x_size=10, int y_size=10, int z_size=10, double fzero_red=1, double fzero_blue=1);
    Lattice3D(const Lattice3D& other);
    ~Lattice3D();

    /// operations
    /// calculations
    void equilibriumIni(); /// < replace all distribution functions with the equilibrium distribution
    void balance(double& mass, double& momentum)const; /// < monitor overall mass and momentum
    void mass_balance(double& liquid_mass, double& gas_mass)const;
    void overallRho();
    direction3D directions(int x, int y, int z)const; /// < calculates positions of neighboring sites (periodical)
    const Vector3D getGradient(int x, int y, int z)const; /// < calculates the color gradient on the position (x,y,z)
    void linearIndex(int index, int& x, int& y, int& z)const;
        
    /// LB steps
    void streamAll(int threads = 1); /// < streaming step
    bool collideAll(int threads = 1, bool gravity = false, bool isLimitActive = true); /// < collision step

    /// walls
    void closedBox(); /// < initialize the Lattice3D (set up walls and calculate rho)
    void bottomWall(); /// < turns the bottom Cells into walls (set up walls and calculate rho)
    void genericWall(std::vector<double> x, std::vector<double> y, std::vector<double> z,  const Vector3D& u_w);
    void lidDrivenCavity(const Vector3D& u_w); /// < initialize the Lattice3D with moving top wall
    void shearWall(const Vector3D& u_w);    /// < initialize the Lattice3D with moving left wall

    /// accessors
    const DimSet3D getSize()const; /// < get the extend of the Lattice3D
    const field3D getData()const{return *data;}; /// < get the data field3D
    const Cell3D getCell(int x, int y, int z)const{return (*data)[x][y][z];};  /// < get a Cell
    const ParamSet getParams()const{return param;}; /// < get the paramter set
    const DistributionSetType3D getF(int x, int y, int z)const{return (*data)[x][y][z].getF();};          /// < get F

    void setData(const field3D& ndata, int x, int y, int z); /// < set the data field3D (and size)
    void setCell(int x, int y, int z, const Cell3D& ncell);    /// < set a Cell
    void setF(int x, int y, int z, int color, const array3D& nf);
    void setF(int x, int y, int z, int color, int index, double value);
    void setParams(const ParamSet& newParam){param = newParam;}; /// < set a new parameter set
    
    /// operators
    Lattice3D& operator=(const Lattice3D& other);
    const bool operator==(const Lattice3D& other)const;

private:
    int xsize, ysize, zsize;   /// < extent of the Lattice3D
    field3D * data;    
    ParamSet param;     /// < set of parameters used during the simulation

    void streamAndBouncePull(Cell3D& tCell, const direction3D& dir)const; /// < internal streaming mechanism with bounce back
};

/// calculates the equilibrium distribution based of a cell
const DistributionSetType3D eqDistro(const ColSet& rho_k, const VeloSet3D& u, const DistributionSetType3D& phi);

const DistributionSetType3D calculate_forcing_term(Vector3D G, VeloSet3D u);

#endif // LATTICE3D_H

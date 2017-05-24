#ifndef LATTICE2D_NO_INTER_H
#define LATTICE2D_NO_INTER_H

#include<iostream>
#include<vector>
#include<omp.h>

#include"Lattice2D.h"

/// contains the domain of the simulation with LB operators
class Lattice2D_no_inter
{
public:
    /// Lifecycle
    Lattice2D_no_inter(int x_size=10, int y_size=10,double fzero_red=1, double fzero_blue=1);
    Lattice2D_no_inter(const Lattice2D_no_inter& other);
    ~Lattice2D_no_inter();

    /// operations
    /// calculations
    void equilibriumIni(); /// < replace all distribution functions with the equilibrium distribution
    void balance(double& mass, double& momentum)const; /// < monitor overall mass and momentum
    void mass_balance(double& liquid_mass, double& gas_mass)const;
    void overallRho();
    direction2D directions(int x, int y)const; /// < calculates positions of neighboring sites (periodical)
             
    /// LB steps
    void streamAll(int threads = 1); /// < streaming step
    bool collideAll(int threads = 1, bool gravity = false, bool isLimitActive = true); /// < collision step
    void evaluateBoundaries(int threads = 1);

    /// walls
    void closedBox(); /// < initialize the Lattice2D_no_inter (set up walls and calculate rho)
    void bottomWall(); /// < initialize the Lattice2D_no_inter (set up walls and calculate rho)
    void genericWall(std::vector<double> x, std::vector<double> y,  const Vector2D& u_w);
    void lidDrivenCavity(const Vector2D& u_w); /// < initialize the Lattice2D_no_inter with moving top wall
    void shearWall(const Vector2D& u_w);    /// < initialize the Lattice2D_no_inter with moving left wall
    void setShearProfile(double gradient, double offset); /// < initialize linear shear profile according to V_y = m x + n

    /// accessors
    const DimSet2D getSize()const; /// < get the extend of the Lattice2D_no_inter
    inline const field2D getData()const{return *data;}; /// < get the data field2D
    const field2D getData(int cutoff)const; 
    inline const Cell2D getCell(int x, int y)const{return (*data)[x][y];};  /// < get a Cell
    inline const ParamSet getParams()const{return param;}; /// < get the paramter set
    inline const Boundaries getBoundaries()const{return bound;};
    inline const int getOffset()const{return offset;};
    inline const DistributionSetType2D getF(int x, int y)const{return (*data)[x][y].getF();};          /// < get F

    void setData(const field2D& ndata, int x, int y); /// < set the data field2D (and size)
    void setCell(int y, int x, const Cell2D& ncell);    /// < set a Cell
    void setF(int x, int y, int color, const array2D& nf);
    void setF(int x, int y, int color, int index, double value);
    void setBoundaries(const Boundaries& newBound);
    void setOffset(int o){offset = o;};
    inline void setParams(const ParamSet& newParam){param = newParam;}; /// < set a new parameter set

    void linearIndex(int index, int& x, int& y)const;
      
    /// boundary treatment
    const bool isBoundary(int x, int y)const;
    void buildWalls(); /// 
    
    const Cell2D boundaryNorthPres(const Cell2D& tmp, ColSet rho)const;
    const Cell2D boundaryNorthVelo(const Cell2D& tmp, double uy)const;

    const Cell2D boundarySouthPres(const Cell2D& tmp, ColSet rho)const;
    const Cell2D boundarySouthVelo(const Cell2D& tmp, double uy)const;

    const Cell2D boundaryWestPres(const Cell2D& tmp, ColSet rho)const;
    const Cell2D boundaryWestVelo(const Cell2D& tmp, double ux)const;

    const Cell2D boundaryEastPres(const Cell2D& tmp, ColSet rho)const;
    const Cell2D boundaryEastVelo(const Cell2D& tmp, double ux)const;
    
    // corners
    const Cell2D cornerNorthWest(const Cell2D& tmp, ColSet rho)const;
    const Cell2D cornerNorthEast(const Cell2D& tmp, ColSet rho)const;
    const Cell2D cornerSouthWest(const Cell2D& tmp, ColSet rho)const;
    const Cell2D cornerSouthEast(const Cell2D& tmp, ColSet rho)const;

    /// operators
    Lattice2D_no_inter& operator=(const Lattice2D_no_inter& other);
    const bool operator==(const Lattice2D_no_inter& other)const;

private:
    int xsize, ysize;   /// < extent of the Lattice2D_no_inter
    field2D * data;    
    ParamSet param;     /// < set of parameters used during the simulation
    Boundaries bound;
    int offset;

    inline const double get_shearvelocity_x(int x, double u_0, double u_1)const{return (u_1 - u_0) / xsize * x + u_0;};
    inline const double get_shearvelocity_y(int y, double u_0, double u_1)const{return (u_1 - u_0) / ysize * y + u_0;};
    void streamAndBouncePull(Cell2D& tCell, const direction2D& dir)const; /// < internal streaming mechanism with bounce back
    const array2D eqDistro_no_inter(double rho, const Vector2D& u);
};


#endif // LATTICE2D_NO_INTER_H

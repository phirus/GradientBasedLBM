#ifndef LATTICE_H
#define LATTICE_H
#include<iostream>
#include<omp.h>

#include"2D/Cell2D.h"
#include"ParamSet.h"
// #include"Timetrack.h"

/// custom typedef for the whole field of cells
typedef boost::multi_array<Cell2D,2> field;

/// contains the domain of the simulation with LB operators
class Lattice
{
public:
    /// Lifecycle
    Lattice(int x_size=10, int y_size=10,double fzero_red=1, double fzero_blue=1);
    Lattice(const Lattice& other);
    ~Lattice();

    /// operations
    /// calculations
    void equilibriumIni(); /// < replace all distribution functions with the equilibrium distribution
    void balance(double& mass, double& momentum)const; /// < monitor overall mass and momentum
    void mass_balance(double& liquid_mass, double& gas_mass)const;
    void overallRho();
    direction2D directions(int x, int y)const; /// < calculates positions of neighboring sites (periodical)
    const Vector2D getGradient(int x, int y)const; /// < calculates the color gradient on the position (y,x)
        
    /// LB steps
    void streamAll(int threads = 1); /// < streaming step
    bool collideAll(int threads = 1, bool gravity = false, bool isLimitActive = true); /// < collision step

    /// walls
    void closedBox(); /// < initialize the Lattice (set up walls and calculate rho)
    void bottomWall(); /// < initialize the Lattice (set up walls and calculate rho)

    /// accessors
    const ColSet getSize()const; /// < get the extend of the Lattice
    const field getData()const{return *data;}; /// < get the data field
    const Cell2D getCell(int x, int y)const{return (*data)[x][y];};  /// < get a Cell
    const ParamSet getParams()const{return param;}; /// < get the paramter set
    const DistributionSetType2D getF(int x, int y)const{return (*data)[x][y].getF();};          /// < get F

    void setData(const field& ndata, int x, int y); /// < set the data field (and size)
    void setCell(int y, int x, const Cell2D& ncell);    /// < set a Cell
    void setF(int x, int y, int color, const array2D& nf);
    void setF(int x, int y, int color, int index, double value);
    void setParams(const ParamSet& newParam){param = newParam;}; /// < set a new parameter set
    
    /// operators
    Lattice& operator=(const Lattice& other);
    const bool operator==(const Lattice& other)const;

private:
    int xsize, ysize;   /// < extent of the Lattice
    field * data;    
    ParamSet param;     /// < set of parameters used during the simulation

    inline void linearIndex(int index, int& x, int& y)const;
    void streamAndBouncePull(Cell2D& tCell, const direction2D& dir)const; /// < internal streaming mechanism with bounce back
};

/// calculates the equilibrium distribution based of a cell
const DistributionSetType2D eqDistro(const ColSet& rho_k, const VeloSet2D& u, const DistributionSetType2D& phi);

const DistributionSetType2D calculate_forcing_term(Vector2D G, VeloSet2D u);

#endif // LATTICE_H

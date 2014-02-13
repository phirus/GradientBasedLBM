#ifndef LATTICE_H
#define LATTICE_H

#include<omp.h>

#include"Cell.h"
#include"ParamSet.h"
#include"Timetrack.h"

/// custom typedef for the whole field of cells
typedef boost::multi_array<Cell,2> field;

/// contains the domain of the simulation with LB operators
class Lattice
{
public:
    Lattice(int x_size=10, int y_size=10,double fzero_red=1, double fzero_blue=1);
    Lattice(const Lattice& other);
    ~Lattice();
    Lattice& operator=(const Lattice& other);
    const bool operator==(const Lattice& other)const;

    /// set-methods
    void setData(const field& ndata, int x, int y); /// < set the data field (and size)
    void setCell(int y, int x, const Cell& ncell);    /// < set a Cell
    void setF(int x, int y, int color, const array& nf);
    void setF(int x, int y, int color, int index, double value);
    void setParams(const ParamSet& newParam){param = newParam;}; /// < set a new parameter set

    /// get-methods
    const ColSet getSize()const; /// < get the extend of the Lattice
    const field getData()const{return *data;}; /// < get the data field
    const Cell getCell(int x, int y)const{return (*data)[x][y];};  /// < get a Cell
    const ParamSet getParams()const{return param;}; /// < get the paramter set
    const DistributionSetType getF(int x, int y)const{return (*data)[x][y].getF();};          /// < get F

    /// calculations
    void equilibriumIni(); /// < replace all distribution functions with the equilibrium distribution
    void balance(double& mass, double& momentum)const; /// < monitor overall mass and momentum
    void mass_balance(double& liquid_mass, double& gas_mass)const;
    void overallRho();
    direction directions(int x, int y)const; /// < calculates positions of neighboring sites (periodical)
    const Vector getGradient(int x, int y)const; /// < calculates the color gradient on the position (y,x)

    /// walls
    void closedBox(); /// < initialize the Lattice (set up walls and calculate rho)
    void bottomWall(); /// < initialize the Lattice (set up walls and calculate rho)

    /// LB steps
    void streamAll(int threads = 1); /// < streaming step
    void collideAll(int threads = 1, bool gravity = false, bool isLimitActive = true); /// < collision step

private:
    int xsize, ysize;   /// < extent of the Lattice
    field * data;    
    ParamSet param;     /// < set of parameters used during the simulation

    inline void linearIndex(int index, int& x, int& y)const;
    void streamAndBouncePull(Cell& tCell, const direction& dir)const; /// < internal streaming mechanism with bounce back
};

/// calculates the equilibrium distribution based of a cell
const DistributionSetType eqDistro(const ColSet& rho_k, const VeloSet& u, const DistributionSetType& phi);

const DistributionSetType calculate_forcing_term(Vector G, VeloSet u);

#endif // LATTICE_H

#ifndef LATTICE_H
#define LATTICE_H

#include<omp.h>
#include<fstream>
#include<string>

#include"cell.h"
#include"paramset.h"

/// custom typedef for the whole field of cells
typedef boost::multi_array<Cell,2> field;

/// calculates the equilibrium distribution based of a cell
const FSet eqDistro(const Cell& tmp, const FSet& phi);

/// computes a difference array (needed for MRT)
const array arrayDiff(const array &one, const array &two);

/// contains the domain of the simulation with LB operators
class Lattice
{
public:
    Lattice(int x_size=10, int y_size=10,double fzero_red=1, double fzero_blue=1);

    void initialize(); /// < initialize the Lattice (set up walls and calculate rho)
    void balance(double& mass, double& momentum); /// < monitor overall mass and momentum
    void overallRho();

    void streamAll(int threads = 0); /// < streaming step
    void collideAll(int threads = 0, bool gravity = false); /// < collision step

    direction directions(int x, int y)const; /// < calculates positions of neighboring sites (periodical)
    const Vector getGradient(int x, int y)const; /// < calculates the color gradient on the position (y,x)

    void equilibriumIni(); /// < replace all distribution functions with the equilibrium distribution

    // access internal elements
    const field getData()const; /// < get the data field
    void setData(const field& ndata, int x, int y); /// < set the data field (and size)

    const ColSet getSize()const; /// < get the extend of the Lattice

    void setCell(int y, int x, const Cell& ncell);    /// < set a Cell
    const Cell getCell(int x, int y)const;                 /// < get a Cell

    const ParamSet getParams()const; /// < get the paramter set
    void setParams(const ParamSet& newParam); /// < set a new parameter set

    const FSet getF(int x, int y)const;          /// < get F
    void setF(int x, int y, int color, const array& nf);
    void setF(int x, int y, int color, int index, double value);

    void techplotOutput(int iterNum, bool vebose = false);
    void vtkOutput(int iterNum);

    const bool operator==(const Lattice& oher)const;

//    void collideGravity();

private:

    field data;         /// < conatins nearly all the data
    int xsize, ysize;   /// < extent of the Lattice
    ParamSet param;     /// < set of parameters used during the simulation

    inline void linearIndex(int index, int& x, int& y)const;
    void streamAndBouncePull(Cell& tCell, const direction& dir)const; /// < internal streaming mechanism with bounce back
};

#endif // LATTICE_H

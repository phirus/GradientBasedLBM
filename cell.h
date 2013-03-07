#ifndef CELL_H
#define CELL_H

#include"constants.h"

using namespace std;

/// all distributions are provided as a Cell
/** all distributions are provided as a Cell with an additional element representing
    the density to reduce redundance when calculating rho or values like the color field
    or deltaRho    */

class Cell
{
public:
    Cell(double fzero_red=1, double fzero_blue=1, bool solid = false); // constr
    Cell(const array& finiRed,const array& finiBlue); // constr

    void calcRho();             /// < calculates both densities
    const double calcPsi()const;     /// < calculates the color field based on the densities

    const Vector calcU()const;      /// < calculates the velocities based on the densities

    const double getDeltaRho()const; /// < returns rho_red - rho_blue

    ///  access the internal elements
    const FSet getF()const;
    void setF(const FSet& newF);
    const ColSet getRho()const;
    const bool getIsSolid()const;
    void setIsSolid(bool tmp);

    const bool operator==(const Cell& other)const;

private:
    FSet f;                     /// < set of two distributions
    ColSet rho;                    /// < rho_r = rho[0], rho_b = rho[1]
    bool isSolid;               /// < used to mark solid cells
};

#endif // CELL_H

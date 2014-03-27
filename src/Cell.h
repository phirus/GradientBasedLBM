/** all distributions are provided as a Cell
* with an additional element representing the density to reduce redundance 
* when calculating rho or values like the color field or deltaRho
*/

#ifndef CELL_H
#define CELL_H

#include"Constants.h"
using namespace std;

class Cell
{
public:
    /// Lifecylce
    Cell(double fzero_dense=1, double fzero_dilute=1, bool solid = false); /// < construcor
    Cell(const array& finiDense,const array& finiDilute); // constr
    Cell(const DistributionSetType& newF); // constr

    /// operators
    const bool operator==(const Cell& other)const;

    /// operations
    void calcRho();                                         /// < calculates both densities, the velocity and delta rho
    const double calcPsi()const;                            /// < calculates the color field based on the densities
 
    /// acsessors
    inline void setF(const DistributionSetType& newF){f = newF;};
    inline void setIsSolid(bool tmp){isSolid = tmp;};

    inline const DistributionSetType getF()const{return f;};
    inline const ColSet getRho()const{return rho;};
    inline const bool getIsSolid()const{return isSolid;};
    inline const VeloSet getU()const{return u;};
    inline const double getDeltaRho()const{return delta;};  /// < returns rho_red - rho_blue
   

private:
    DistributionSetType f;         /// < set of two distributions
    ColSet rho;                    /// < rho_r = rho[0], rho_b = rho[1]
    VeloSet u;                      /// < velocity vector resultińg from the distribution
    bool isSolid;                  /// < used to mark solid cells
    double delta;
};

#endif // CELL_H

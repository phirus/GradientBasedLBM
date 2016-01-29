/** all distributions are provided as a Cell
* with an additional element representing the density to reduce redundance 
* when calculating rho or values like the color field or deltaRho
*/

#ifndef CELL3D_H
#define CELL3D_H

#include"Constants3D.h"
using namespace std;

class Cell3D
{
public:
    /// Lifecylce
    Cell3D(double fzero_dense=1, double fzero_dilute=1, bool solid = false); /// < construcor
    Cell3D(const array3D& finiDense,const array3D& finiDilute); // constr
    Cell3D(const DistributionSetType3D& newF); // constr

    /// operators
    const bool operator==(const Cell3D& other)const;

    /// operations
    void calcRho();                                         /// < calculates both densities, the velocity and delta rho
    const double calcPsi()const;                            /// < calculates the color field based on the densities
 
    /// acsessors
    inline void setF(const DistributionSetType3D& newF){f = newF;};
    inline void setIsSolid(bool tmp){isSolid = tmp;};
    inline void setSolidVelocity(const VeloSet3D& newU){if(isSolid == true) u = newU;};

    inline const DistributionSetType3D getF()const{return f;};
    inline const ColSet getRho()const{return rho;};
    inline const bool getIsSolid()const{return isSolid;};
    inline const VeloSet3D getU()const{return u;};
    inline const double getDeltaRho()const{return delta;};  /// < returns rho_red - rho_blue

private:
    DistributionSetType3D f;         /// < set of two distributions
    ColSet rho;                    /// < rho_r = rho[0], rho_b = rho[1]
    VeloSet3D u;                      /// < velocity vector resultińg from the distribution
    bool isSolid;                  /// < used to mark solid cells
    double delta;
};

#endif // CELL3D_H

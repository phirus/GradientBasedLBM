#include "Cell2D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Cell2D::Cell2D(double fzero_dense, double fzero_dilute, bool solid):isSolid(solid),delta(0)
{
    f[0][0] = fzero_dense;
    f[1][0] = fzero_dilute;

    for (int i=1; i<=8; i++)
    {
        f[0][i]=0;
        f[1][i]=0;
    }
    rho[0] = 0;
    rho[1] = 0;

    u[0].x = 0;
    u[0].y = 0;
    u[1].x = 0;
    u[1].y = 0;

}

// like a copy constructor for the bulk phase
Cell2D::Cell2D(const array2D& finiDense, const array2D& finiDilute):isSolid(false),delta(0)
{
    f.at(0) = finiDense;
    f.at(1) = finiDilute;

    rho[0] = 0;
    rho[1] = 0;

    u[0].x = 0;
    u[0].y = 0;
    u[1].x = 0;
    u[1].y = 0;
}

Cell2D::Cell2D(const DistributionSetType2D& newF):isSolid(false),delta(0)
{
    f.at(0) = newF.at(0);
    f.at(1) = newF.at(1);

    rho[0] = 0;
    rho[1] = 0;

    u[0].x = 0;
    u[0].y = 0;
    u[1].x = 0;
    u[1].y = 0;

}

//=========================== OPERATORS ===========================

const bool Cell2D::operator==(const Cell2D& other)const {
    bool exit = true;
    if (isSolid != other.getIsSolid()) exit = false;

    ColSet rhoOther = other.getRho();
    if(rho[0] != rhoOther[0] || rho[1] != rhoOther[1]) exit = false;

    DistributionSetType2D fOther = other.getF();
    for(int color=0; color <2; color++){
        for(int q=0;q<9;q++){
            if(f[color][q] != fOther[color][q]) exit = false;
            }
        }
    return exit;
}

//=========================== OPERATIONS ===========================

void Cell2D::calcRho()
{
    // initialize density
    rho[0] = 0;
    rho[1] = 0;

    if (isSolid == false)
    {
        // initialize velocities
        u[0].x = 0;
        u[0].y = 0;
        u[1].x = 0;
        u[1].y = 0;

        // iterate
        for (int i=0; i<9; i++)
        {
            rho[0] += f[0][i];
            rho[1] += f[1][i];
            u[0].x += ( f[0][i]) * DIRECTION_2D[i].x;
            u[0].y += ( f[0][i]) * DIRECTION_2D[i].y;
            u[1].x += ( f[1][i] ) * DIRECTION_2D[i].x;
            u[1].y += ( f[1][i] ) * DIRECTION_2D[i].y;
        }

        delta = rho[0]-rho[1];
        
        if(rho[0] >0) {
            u[0].x /= rho[0];
            u[0].y /= rho[0];
        }
        if(rho[1] >0) {
            u[1].x /= rho[1];
            u[1].y /= rho[1];
        }
    }
}

const double Cell2D::calcPsi()const
{
    double rhoSum = sum(rho);
    if(rhoSum > 0) return (rho[0] - rho[1])/( rhoSum ); // < prevent division by 0
    else return 0;
}
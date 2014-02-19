#include "Cell.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Cell::Cell(double fzero_red, double fzero_blue, bool solid):isSolid(solid),delta(0)
{
    f[0][0] = fzero_red;
    f[1][0] = fzero_blue;

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
Cell::Cell(const array& finiRed, const array& finiBlue):isSolid(false),delta(0)
{
    f.at(0) = finiRed;
    f.at(1) = finiBlue;

    rho[0] = 0;
    rho[1] = 0;

    u[0].x = 0;
    u[0].y = 0;
    u[1].x = 0;
    u[1].y = 0;
}

Cell::Cell(const DistributionSetType& newF):isSolid(false),delta(0)
{
    f.at(0) = newF.at(0);
    f.at(1) = newF.at(0);

    rho[0] = 0;
    rho[1] = 0;

    u[0].x = 0;
    u[0].y = 0;
    u[1].x = 0;
    u[1].y = 0;

}

//=========================== OPERATORS ===========================

const bool Cell::operator==(const Cell& other)const {
    bool exit = true;
    if (isSolid != other.getIsSolid()) exit = false;

    ColSet rhoOther = other.getRho();
    if(rho[0] != rhoOther[0] || rho[1] != rhoOther[1]) exit = false;

    DistributionSetType fOther = other.getF();
    for(int color=0; color <2; color++){
        for(int q=0;q<9;q++){
            if(f[color][q] != fOther[color][q]) exit = false;
            }
        }
    return exit;
}

//=========================== OPERATIONS ===========================

void Cell::calcRho()
{
    // initialize
    rho[0] = 0;
    rho[1] = 0;

    u[0].x = 0;
    u[0].y = 0;
    u[1].x = 0;
    u[1].y = 0;

    // iterate
    if (isSolid == false)
    {
        for (int i=0; i<9; i++)
        {
            rho[0] += f[0][i];
            rho[1] += f[1][i];
            u[0].x += ( f[0][i]) * DIRECTION[i].x;
            u[0].y += ( f[0][i]) * DIRECTION[i].y;
            u[1].x += ( f[1][i] ) * DIRECTION[i].x;
            u[1].y += ( f[1][i] ) * DIRECTION[i].y;
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

const double Cell::calcPsi()const
{
    double rhoSum = sum(rho);
    if(rhoSum > 0) return (rho[0] - rho[1])/( rhoSum ); // < prevent division by 0
    else return 0;
}
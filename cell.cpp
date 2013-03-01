#include "cell.h"

Cell::Cell(double fzero_red, double fzero_blue, bool solid):isSolid(solid)
{
//    isSolid = solid;
    f[0][0] = fzero_red;
    f[1][0] = fzero_blue;

    for (int i=1; i<=8; i++)
    {
        f[0][i]=0;
        f[1][i]=0;
    }
    rho[0] = 0;
    rho[1] = 0;
}
// like a copy constructor for the bulk phase
Cell::Cell(const array& finiRed, const array& finiBlue):isSolid(false)
{
    f.at(0) = finiRed;
    f.at(1) = finiBlue;

    rho[0] = 0;
    rho[1] = 0;
}

void Cell::calcRho()
{
    // initialize
    rho[0] = 0;
    rho[1] = 0;

    // iterate
    if (isSolid == false)
    {
        for (int i=0; i<9; i++)
        {
            rho[0] += f[0][i];
            rho[1] += f[1][i];
        }
    }
}

const double Cell::calcPsi()const
{
    double rhoSum = sum(rho);
    if(rhoSum > 0) return (rho[0] - rho[1])/( rhoSum ); // < prevent division by 0
    else return 0;
}

const Vector Cell::calcU()const
{
    // initialize
    double rhoSum = sum(rho);
    Vector u(0,0);

    // iterate
    if (isSolid == false && rhoSum > 0) // < prevent division by 0
    {
        for (int i=0; i<9; i++)
        {
            u.x += ( f[0][i] + f[1][i] ) * e[i].x;
            u.y += ( f[0][i] + f[1][i] ) * e[i].y;
        }
        u.x /= rhoSum;
        u.y /= rhoSum;
    }
    return u;
}


const FSet Cell::getF()const
{
    return f;
}

void Cell::setF(const FSet& newF){
    f = newF;
}

const ColSet Cell::getRho()const
{
    return rho;
}

const double Cell::getDeltaRho()const
{
    return (rho[0]-rho[1]);
}

void Cell::setIsSolid(bool tmp)
{
    isSolid = tmp;
}
const bool Cell::getIsSolid()const
{
    return isSolid;
}

const bool Cell::operator==(const Cell& other)const {
    bool exit = true;
    if (isSolid != other.getIsSolid()) exit = false;

    ColSet rhoOther = other.getRho();
    if(rho[0] != rhoOther[0] || rho[1] != rhoOther[1]) exit = false;

    FSet fOther = other.getF();
    for(int color=0; color <2; color++){
        for(int q=0;q<9;q++){
            if(f[color][q] != fOther[color][q]) exit = false;
            }
        }
    return exit;
}

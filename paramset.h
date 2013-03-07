#ifndef PARAMSET_H
#define PARAMSET_H

#include"definitions.h"
#include<string.h>

using namespace std;

/// structure for interpolation paramters
struct Interpol
{
    double chi, eta, kappa, lambda, ny;
};

struct RelaxationPar
{
    double s_2,s_3,s_5;
};

/// collection of all parameters used during the simulation
class ParamSet
{
public:
    ParamSet(double omR=1, double omB =1,double rhoR = 1, double gamma = 1000, double alB=0.2, double d=0.1, double bet=0.99, double sig = 1e-4, double grav = 0.001, double s2 = 1, double s3 = 1, double s5 = 1.2, double sound = 1484, double length = 0.001); /// < consructor

    const FSet getPhi();     /// calculates phi, based on alpha_b and rho (density ratio)
    const double getOmega(double psi);    /// return omega, based on inter and the color field

    /// access internal elements
    void setOmega(double omR, double omB, double d);
    void setAlpha(double alB);
    void setRatio(double rhoR, double ratio);
    void setBeta(double bet);
    void setSigma(double sig);
    const double getBeta()const;
    const ColSet getAk(double omega)const;
    const double getG()const;
    void setG(double grav);
    const double getRhoR()const{return rhoRed;};
    const Interpol getInter()const; // function is only used for unit testing
    void setRelaxation(double s_2, double s_3, double s_5);
    const RelaxationPar getRelaxation()const;

    const bool operator==(const ParamSet& other)const;

private:
    double omegaRed, omegaBlue; /// < relaxation parameters for both colors
    double rhoRed, gamma;       /// < liquid density and density ratio
    double alphaRed, alphaBlue; /// < alpha_b, relevant for equilibrium distribution
    double delta;               /// < thickness of boundary layer, relevant for equilibrium distribution
    double beta;                /// < beta, relevant for recoloring operator
    double sigma;               /// < surface tension
    double g;                   /// < gravity
    Interpol inter;             /// < interpolation parameters for finding omega, relevant for equilibrium distribution
    RelaxationPar relax;

    double c_s;   // speed of sound
    double timestep;
    double spacestep;


    void calcInter();           /// < calculate the interpolation paramters based on omega and delta
    void calcAlR();
};

#endif

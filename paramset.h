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
    ParamSet(double omR=1, double omB =1,double rhoR = 1, double gammaIni = 1000, double alB=0.2, double deltaIni=0.1, double betaIni=0.99, double sigmaIni = 1e-4, double c_sIni = 1484, double length = 0.001); /// < consructor

    const FSet getPhi();                    /// < calculates phi, based on alpha_b and rho (density ratio)
    const double getOmega(double psi);      /// < return omega, based on inter and the color field

    /// access internal elements
    void setOmega(double omR, double omB, double d);
    void setRatio(double rhoR, double ratio);
    void setAlpha(double alB);
    void setBeta(double bet);
    void setSigma(double sig);
    void setDeltaX(double dx);
    void setSoundSpeed(double sos);

    const double getBeta()const;
    const ColSet getAk(double omega)const;
    const double getG()const{return gravity;};
    void setG(double grav);
    const double getRhoR()const{return rhoRed;};

    void setRelaxation(double s_2, double s_3, double s_5);

    const boost::array<double,13> getEverything()const;
    const RelaxationPar getRelaxation()const;
    const Interpol getInter()const;
    const bool operator==(const ParamSet& other)const;

private:
    double omegaRed, omegaBlue; /// < relaxation parameters for both colors
    double rhoRed, gamma;       /// < liquid density and density ratio
    double alphaRed, alphaBlue; /// < alpha_b, relevant for equilibrium distribution
    double delta;               /// < thickness of boundary layer, relevant for equilibrium distribution
    double beta;                /// < beta, relevant for recoloring operator
    double sigma;               /// < surface tension
    Interpol inter;             /// < interpolation parameters for finding omega, relevant for equilibrium distribution
    RelaxationPar relax;

    double c_s ;                /// < speed of sound / m * s^-1
    double timestep;            /// < timestep /s
    double spacestep;           /// < spacestep /m
    double gravity;             /// < gravity /-

    void calcInter();           /// < calculate the interpolation paramters based on omega and delta
    void calcAlR();
    inline void calcTimestep(){timestep = spacestep / (c_s * sqrt(3) );} ;
    inline void calcGravity(double g = 9.81){gravity = g * timestep * timestep / spacestep;};
};

#endif

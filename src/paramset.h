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
    ParamSet(double omR=1, double omB =1,double rhoR = 1, double gammaIni = 1000, double alB=0.2, double deltaIni=0.1, double betaIni=0.99, double sigmaIni = 1e-4, double c_sIni = 1484, double length = 0.001, double g = 9.81); /// < consructor

    /// get-methods, including calculations if necessary
    const FSet getPhi()const;                   /// < calculates phi, based on alpha_b and rho (density ratio)
    const double getOmega(double psi)const;          /// < return omega, based on inter and the color field
    const ColSet getAk(double omega)const;
    const boost::array<double,14> getEverything()const;

    const double getBeta()const{return beta;};
    const double getG()const{return gravity;};
    const double getRhoR()const{return rhoRed;};
    const Interpol getInter()const{return inter;};
    const RelaxationPar getRelaxation()const{return relax;};
    const double getDeltaT()const{return timestep;};

    // needed only for output
    const double getDeltaX()const{return spacestep;};
    const double getSoundspeed()const{return c_s;};
    const double getOmegaRed()const{return omegaRed;};
    const double getOmegaBlue()const{return omegaBlue;};
    const double getInterfaceThickness()const{return delta;};
    const double getGamma()const{return gamma;};
    const double getAlpha()const{return alphaBlue;};
    const double getSigma()const{return sigma;};
    const double getSpeedLimit()const{return speedlimit;};

    /// set-methods, including calculations if necessary
    void setOmega(double omR, double omB, double d);
    void setAlpha(double alB);
    void setRatio(double rhoR, double ratio);

    void setDeltaX(double dx);
    void setSoundSpeed(double sos);
    void setOriginalG(double g);

    void setRelaxation(double s_2, double s_3, double s_5);

    void setBeta(double bet){beta=bet;};
    void setSigma(double sig){sigma=sig;};

    /// overloaded == Operator
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
    double original_g;           /// < gravity / m * s^-2
    double gravity;             /// < gravity /-

    double speedlimit;          /// < maximum allowed velocity

    void calcInter();           /// < calculate the interpolation paramters based on omega and delta
    void calcAlR();
    void calcTimestep(); 

    inline void calcGravity(){gravity = original_g * timestep * timestep / spacestep;};
};

#endif

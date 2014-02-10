#ifndef PARAMSET_H
#define PARAMSET_H

#include"Definitions.h"
#include<string.h>

using namespace std;


/// collection of all parameters used during the simulation
class ParamSet
{
public:
    /// Lifecycle
    ParamSet(double omR = 1, double omB = 1,double rhoR = 1, double gammaIni = 1000,double sigmaIni = 1e-4, double g = 9.81, double c_limit = 1, double t_step = 1e-3, RelaxationPar rel = RelaxationPar(1,1,1), double alB = 0.2, double deltaIni = 0.1, double betaIni = 0.99); /// < consructor

    /// get-methods, including calculations if necessary
    const DistributionSetType getPhi()const;                   /// < calculates phi, based on alpha_b and rho (density ratio)
    const double getOmega(double psi)const;          /// < return omega, based on inter and the color field
    const ColSet getAk(double omega)const;
    const boost::array<double,12> getEverything()const;

    const double getBeta()const{return beta;};
    const double getG()const{return gravity;};
    const double getRhoR()const{return rhoRed;};
    const Interpol getInter()const{return inter;};
    const RelaxationPar getRelaxation()const{return relax;};
    const double getDeltaT()const{return timestep;};

    // needed only for output
    const double getOmegaRed()const{return omegaRed;};
    const double getOmegaBlue()const{return omegaBlue;};
    const double getInterfaceThickness()const{return delta;};
    const double getGamma()const{return gamma;};
    const double getAlpha()const{return alphaBlue;};
    const double getSigma()const{return sigma;};
    const double getSpeedlimit()const{return speedlimit;};

    /// set-methods, including calculations if necessary
    void setOmega(double omR, double omB, double d);
    void setAlpha(double alB);
    void setRatio(double rhoR, double ratio);
    void setRelaxation(double s_2, double s_3, double s_5);
    void setBeta(double bet){beta=bet;};

    /// overloaded == Operator
    const bool operator==(const ParamSet& other)const;

private:
    // given
    double omegaRed, omegaBlue; /// < relaxation parameters for both colors
    double rhoRed, gamma;       /// < liquid density and density ratio
    double alphaRed, alphaBlue; /// < alpha_b, relevant for equilibrium distribution
    double delta;               /// < thickness of boundary layer, relevant for equilibrium distribution
    double beta;                /// < beta, relevant for recoloring operator    
    double sigma;               /// < dimensionless surface tension
    double gravity;             /// < dimensionless gravity 
    // drive through
    double speedlimit;          /// < maximum allowed velocity
    double timestep;            /// < LB timestep
    // deduced
    RelaxationPar relax;
    Interpol inter;             /// < interpolation parameters for finding omega, relevant for equilibrium distribution

    /// operations
    void calcInter();           /// < calculate the interpolation paramters based on omega and delta
    void calcAlR();
};

#endif

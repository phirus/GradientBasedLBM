#include"paramset.h"
ParamSet::ParamSet(double omR, double omB, double rhoR,double g, double alB, double d, double bet, double sig, double grav, double s2, double s3, double s5, double sound, double length):omegaRed(omR),omegaBlue(omB),rhoRed(rhoR),gamma(g),alphaBlue(alB),delta(d),beta(bet), sigma(sig), g(grav), c_s(sound), spacestep(length)
{
    relax.s_2 = s2;
    relax.s_3 = s3;
    relax.s_5 = s5;
    timestep = spacestep / c_s;
    calcInter();
    calcAlR();
}

void ParamSet::setAlpha(double alB)
{
    alphaBlue = alB;
    calcAlR();
}

const FSet ParamSet::getPhi()
{
    FSet phi;
    phi.at(0).at(0) = alphaRed;
    phi.at(1).at(0) = alphaBlue;
    for (int i = 1;i<9; i+=2)
    {
        phi.at(0).at(i) = (1-alphaRed)/5;
        phi.at(1).at(i) = (1-alphaBlue)/5;
    }
    for (int i = 2;i<9; i+=2)
    {
        phi.at(0).at(i) = (1-alphaRed)/20;
        phi.at(1).at(i) = (1-alphaBlue)/20;
    }

    return phi;
}

void ParamSet::setOmega(double omR, double omB, double d)
{
    omegaRed = omR;
    omegaBlue = omB;
    delta = d;
    calcInter();
}

void ParamSet::setRatio(double rhoR, double ratio)
{
    rhoRed = rhoR;
    gamma = ratio;
    calcAlR();
}

void ParamSet::setBeta(double bet){
    beta = bet;
}

void ParamSet::setSigma(double sig){
    sigma = sig;
}

const Interpol ParamSet::getInter()const
{
    return inter;
}

const double ParamSet::getOmega(double psi)
{
    if (psi > delta)
    {
        return omegaRed;
    }
    else if (psi > 0 && psi<= delta)
    {
        return inter.chi + inter.eta * psi + inter.kappa * psi * psi;
    }
    else if (psi<=0 && psi >= -delta)
    {
        return inter.chi + inter.lambda * psi + inter.ny *psi *psi;
    }
    else
    {
        return omegaBlue;
    }
}

const double ParamSet::getBeta()const{return beta;}

const ColSet ParamSet::getAk(double omega)const
{
    double A = (9 * omega * sigma) / (2* rhoRed * (1+1/gamma) );
    ColSet Ak = {{A,A}};
    return Ak;
}

const double ParamSet::getG()const{return g;}
void ParamSet::setG(double grav){g = grav;}

void ParamSet::calcInter()
{
    inter.chi    = (2*omegaRed*omegaBlue)/( omegaRed + omegaBlue );
    inter.eta    = 2*(omegaRed - inter.chi)/delta;
    inter.kappa  = -inter.eta / ( 2*delta );
    inter.lambda = 2*(inter.chi - omegaBlue)/delta;
    inter.ny     = inter.lambda / ( 2*delta );
}

void ParamSet::calcAlR(){
    alphaRed = 1- (1- alphaBlue)/gamma;
}

void ParamSet::setRelaxation(double s_2, double s_3, double s_5)
{
    relax.s_2 = s_2;
    relax.s_3 = s_3;
    relax.s_5 = s_5;
}
const RelaxationPar ParamSet::getRelaxation()const{return relax;}

const bool ParamSet::operator==(const ParamSet& other)const{
    if (memcmp(this, &other,sizeof(ParamSet)) == 0) return true;
    return false;
}

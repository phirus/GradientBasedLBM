#include"paramset.h"
ParamSet::ParamSet(double omR, double omB, double rhoR,double g, double alB, double d, double bet, double sig, double s2, double s3, double s5, double sound, double length,  double grav):omegaRed(omR),omegaBlue(omB),rhoRed(rhoR),gamma(g),alphaBlue(alB),delta(d),beta(bet), sigma(sig), c_s(sound), spacestep(length),g(grav)
{
    relax.s_2 = s2;
    relax.s_3 = s3;
    relax.s_5 = s5;
    g = 1;
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

const boost::array<double,13> ParamSet::getEverything()const{
    boost::array<double,13> pinkie;
    pinkie[0] = omegaRed;
    pinkie[1] = omegaBlue;
    pinkie[2] = rhoRed;
    pinkie[3] = gamma;
    pinkie[4] = alphaRed;
    pinkie[5] = alphaBlue;
    pinkie[6] = delta;
    pinkie[7] = beta;
    pinkie[8] = sigma;
    pinkie[9] = c_s;
    pinkie[10] = timestep;
    pinkie[11] = spacestep;
    pinkie[12] = g;

    return pinkie;
}


const bool ParamSet::operator==(const ParamSet& other)const{
    bool control = true;
    {
        boost::array<double,13> foo, bar;
        foo = getEverything();
        bar = other.getEverything();

        if(foo != bar) control = false;
    }
    {
        Interpol finter, binter;
        finter = getInter();
        binter = other.getInter();

        if(finter.chi != binter.chi || finter.eta != binter.eta || finter.kappa != binter.kappa || finter.lambda != binter.lambda || finter.ny != binter.ny) control = false;
    }
    {
        RelaxationPar frelax,brelax;
        frelax = getRelaxation();
        brelax = other.getRelaxation();

        if(frelax.s_2 != brelax.s_2 || frelax.s_3 != brelax.s_3 || frelax.s_5 != brelax.s_5) control = false;
    }
    return control;
}

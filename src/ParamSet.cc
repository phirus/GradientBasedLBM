#include"ParamSet.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

ParamSet::ParamSet(double omR, double omB, double rhoR ,double gammaIni, double sigmaIni, double g, double t_step, double s_step,RelaxationPar3D rel, double b_visc, double alB, double deltaIni, double betaIni):
omegaRed(omR),
omegaBlue(omB),
rhoRed(rhoR),
gamma(gammaIni),
alphaBlue(alB),
delta(deltaIni),
beta(betaIni), 
sigma(sigmaIni), 
gravity(g),
bulk_visco(b_visc), 
timestep(t_step),
spacestep(s_step),
relax(rel)
{
    calcInter();
    calcAlR();
}

//=========================== OPERATIONS ===========================

const DistributionSetType2D ParamSet::getPhi2D()const
{
    DistributionSetType2D phi;
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

const DistributionSetType3D ParamSet::getPhi3D()const
{
    DistributionSetType3D phi;
    phi.at(0).at(0) = alphaRed;
    phi.at(1).at(0) = alphaBlue;
    for (int i = 1;i<19; i+=1)
    {
        if(i == 1 || i == 3 || i == 5 || i == 7 || i == 9 || i == 14)
        {
            phi.at(0).at(i) = (1-alphaRed)/12;
            phi.at(1).at(i) = (1-alphaBlue)/12;
        }
        else 
        {
            phi.at(0).at(i) = (1-alphaRed)/24;
            phi.at(1).at(i) = (1-alphaBlue)/24;
        }        
    }
    return phi;
}

const double ParamSet::getOmega(double psi)const
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

// const ColSet ParamSet::getAk(double omega)const
// {
//     double A = (9 * omega * sigma) / (2* rhoRed * (1+1/gamma) );
//     ColSet Ak = {{A,A}};
//     return Ak;
// }

const ColSet ParamSet::getAk(double omega)const
{
    // double tmp = pow(spacestep,4) / pow(timestep,3);
    double A = (9 * omega * sigma) / (4);
    ColSet Ak = {{A,A}};
    return Ak;
}

const RelaxationPar2D ParamSet::getRelaxation2D(double psi)const
{
    const double omega = getOmega(psi);
    const RelaxationPar2D rel = RelaxationPar2D(omega, getS2(omega), relax.s_3, relax.s_5);
    return rel;
}

const RelaxationPar3D ParamSet::getRelaxation3D(double psi)const
{
    const double omega = getOmega(psi);
    const RelaxationPar3D rel = RelaxationPar3D(omega, getS2(omega), relax.s_3, relax.s_5, relax.s_11, relax.s_17);
    return rel;
}


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
    pinkie[9] = gravity;
    pinkie[10] = bulk_visco;
    pinkie[11] = timestep;
    pinkie[12] = spacestep;
    

    return pinkie;
}

void ParamSet::setOmega(double omR, double omB, double d)
{
    omegaRed = omR;
    omegaBlue = omB;
    delta = d;
    calcInter();
}

void ParamSet::setAlpha(double alB)
{
    alphaBlue = alB;
    calcAlR();
}

void ParamSet::setRatio(double rhoR, double ratio)
{
    rhoRed = rhoR;
    gamma = ratio;
    calcAlR();
}


void ParamSet::setRelaxation(double s_2, double s_3, double s_5)
{
    relax.s_2 = s_2;
    relax.s_3 = s_3;
    relax.s_5 = s_5;
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

        if(finter.chi != binter.chi ) control = false;
        if(finter.eta != binter.eta ) control = false;
        if(finter.kappa != binter.kappa ) control = false;
        if(finter.lambda != binter.lambda ) control = false;
        if(finter.ny != binter.ny) control = false;
    }
    {
        RelaxationPar3D frelax,brelax;
        frelax = getRelaxation3D();
        brelax = other.getRelaxation3D();

        if(frelax.s_2 != brelax.s_2 || frelax.s_3 != brelax.s_3 || frelax.s_5 != brelax.s_5) control = false;
    }
    return control;
}


///////////////////////////// PRIVATE /////////////////////////////

//=========================== OPERATIONS ===========================

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
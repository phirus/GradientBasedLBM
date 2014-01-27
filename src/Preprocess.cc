#include"Preprocess.h"

Preprocess::Preprocess(double Re, double Mo, double Eo, double res, double rl, double gam, double dia,double mu_rate, double soundspeed, double sig, double grav, double s_three, double s_five):
ReynoldsMax(Re), Morton(Mo), Eotvos(Eo),
resolution(res), rho_l(rl), gamma(gam), 
diameter(dia), mu_ratio(mu_rate), 
c_s(soundspeed), sigma(sig), g(grav),
s_3(s_three), s_5(s_five)
{
    deduceAll();
}

void Preprocess::deduceAll(){
	calcTau();
	calcSpeedlimit();
	calcSpacestep();
	calcTimestep();
	calcNu();
    calcS2();
	calcDelRho();
}

const ParamSet Preprocess::getParamSet()const{
	const double omega = 1/tau;
	const double rho_r = 1;  // normalized
    const RelaxationPar relax(s_2,s_3,s_3);
	ParamSet param(omega, omega, rho_r, gamma, convertSigma(), convertG(), speedlimit, timestep, relax);
	return param;
}

void Preprocess::refine(double factor){
	c_s *= factor;			// raise the speedof sound
	ReynoldsMax *=  factor;		// with constant diameter and nu, Reynolds raises with the Speed of sound

	deduceAll();		// deduces all new parameters
}

const bool Preprocess::operator==(const Preprocess& other)const
{
    bool exit = true;
    if(ReynoldsMax != other.getReynoldsMax()) exit = false;
    if(Morton != other.getMorton()) exit = false;
    if(Eotvos != other.getEotvos()) exit = false;
    if(resolution != other.getResolution()) exit = false;
    if(rho_l != other.getRhoL()) exit = false;
    if(gamma != other.getGamma()) exit = false;
    if(diameter != other.getDiameter()) exit = false;
    if(mu_ratio != other.getMuRatio()) exit = false;
    if(c_s != other.getSoundspeed()) exit = false;
    if(sigma != other.getSigma()) exit = false;
    if(g != other.getGPhys()) exit = false;
    if(s_3 != other.getS_3()) exit = false;
    if(s_5 != other.getS_5()) exit = false;


    if(tau != other.getTau()) exit = false;
    if(s_2 != other.getS2()) exit = false;
    if(speedlimit != other.getSpeedlimit()) exit = false;
    if(spacestep != other.getSpacestep()) exit = false;
    if(timestep != other.getTimestep()) exit = false;   
    if(delRho != other.getDelRho()) exit = false;
    if(nu != other.getNu()) exit = false;
    
    return exit;
}
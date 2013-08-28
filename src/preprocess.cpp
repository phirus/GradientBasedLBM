#include"preprocess.h"
Preprocess::Preprocess(double Re, double Mo, double Eo)
{
    ReynoldsMax = Re;
    Morton = Mo;
    Eotvos = Eo;
}

void Preprocess::setResolution(double val){
	resolution = val;
	calcSpacestep();
}

void Preprocess::setRhoL(double val){
	rho_l = val;	
}

void Preprocess::setRho0(double val){
	rho_0 = val;	
}

void Preprocess::setGamma(double val){
	gamma = val;
}

void Preprocess::setDiameter(double val){
	diameter = val;
}

// rough approximation of the terminalrise velocity (Eq. 7-3 from 'Bubbles, Drops, and Particles' by Clift et al.)
const double Preprocess::estimateVelocity()const{
	double u_t = sqrt(2.14 * sigma / (rho_l * diameter) + 0.505 * g * diameter);
	return u_t;
}

void Preprocess::deduceAll(){
	


}
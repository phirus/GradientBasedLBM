#include"preprocess.h"
Preprocess::Preprocess(double Re, double Mo, double Eo)
{
    Reynolds = Re;
    Morton = Mo;
    Eotvos = Eo;
}

void Preprocess::setDiameter(double val){
	diameter = val;
	calcSpacestep();
}

void Preprocess::setResolution(double val){
	resolution = val;
	calcSpacestep();
}

// rough approximation of the terminalrise velocity (Eq. 7-3 from 'Bubbles, Drops, and Particles' by Clift et al.)
const double Preprocess::estimateVelocity()const{
	double u_t = sqrt(2.14 * sigma / (rho_l * diameter) + 0.505 * g * diameter);
	return u_t;
}
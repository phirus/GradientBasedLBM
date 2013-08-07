#include"preprocess.h"
Preprocess::Preprocess(double Re, double Mo, double Eo)
{
    Reynolds = Re;
    Morton = Mo;
    Eotvos = Eo;
}

const double Preprocess::getTau()const{
	double tau = (resolution / Mach_max   * sqrt(3) / Reynolds ) + 0.5;

	return tau;
}
void Preprocess::setDiameter(double val){
	diameter = val;
	calcSpacestep();
}

void Preprocess::setResolution(double val){
	resolution = val;
	calcSpacestep();
}

void Preprocess::calcTimestep(){
    speedlimit = c_s * sqrt(3);             // calculates lattice sound speed
    timestep = spacestep / speedlimit;      // timestep
    speedlimit /= Mach_max;                 // maximum velocity based on the maximum Mach-number
}
#include"preprocess.h"
Preprocess::Preprocess(double Re, double Mo, double Eo, double res, double rl, double gam, double dia, double soundspeed, double sig, double grav):ReynoldsMax(Re), Morton(Mo),Eotvos(Eo),resolution(res), rho_l(rl), gamma(gam), diameter(dia), c_s(soundspeed), sigma(sig), g(grav)
{
    deduceAll();
}

void Preprocess::deduceAll(){
	calcTau();
	calcSpeedlimit();
	calcSpacestep();
	calcTimestep();
	calcNu();
	calcDelRho();
}
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
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

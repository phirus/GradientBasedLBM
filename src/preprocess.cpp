#include"preprocess.h"
Preprocess::Preprocess(double Re, double Mo, double Eo)
{
    ReynoldsMax = Re;
    Morton = Mo;
    Eotvos = Eo;
}

void Preprocess::setResolution(double val){
	resolution = val;
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


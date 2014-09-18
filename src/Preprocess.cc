#include"Preprocess.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Preprocess::Preprocess(double Re, double Mo, double Eo, double res, double rhol, double gamma_ini, double mu_rate, double s_three, double s_five, int width_ini, int height_ini, int max_iter_ini, int vtk__ini, int tecplot_ini, int restart_ini):
ReynoldsMax(Re), Morton(Mo), Eotvos(Eo),
resolution(res), rho_l(rhol), gamma(gamma_ini), 
muRatio(mu_rate), s_3(s_three), s_5(s_five),
width(width_ini), height(height_ini), max_iter(max_iter_ini), vtk_interval(vtk__ini), 
tecplot_interval(tecplot_ini), restart_interval(restart_ini)
{
    deduceAll();
}

//=========================== OPERATIONS ===========================

const ParamSet Preprocess::getParamSet()const{
    const double omega = 1/tau;
    const double rho_r = 1;  // normalized
    const RelaxationPar relax(s_2,s_3,s_5);
    ParamSet param(omega, omega, rho_r, gamma, convertSigma(), g, 1, timestep, spacestep, relax);
    return param;
}

//=========================== OPERATOR ===========================

const bool Preprocess::operator==(const Preprocess& other)const
{
    bool exit = true;
    if(ReynoldsMax != other.getReynoldsMax()) exit = false;
    if(Morton != other.getMorton()) exit = false;
    if(Eotvos != other.getEotvos()) exit = false;
    if(resolution != other.getResolution()) exit = false;
    if(rho_l != other.getRhoL()) exit = false;
    if(gamma != other.getGamma()) exit = false;
    if(muRatio != other.getMuRatio()) exit = false;
    if(s_3 != other.getS_3()) exit = false;
    if(s_5 != other.getS_5()) exit = false;
    if(spacestep != other.getSpacestep()) exit = false;
    if(timestep != other.getTimestep()) exit = false;   
    if(tau != other.getTau()) exit = false;
    if(delRho != other.getDelRho()) exit = false;
    if(c_s != other.getSoundspeed()) exit = false;
    if(nu != other.getNu()) exit = false;
    if(s_2 != other.getS2()) exit = false;

    if(sigma != other.getSigma()) exit = false;
    if(g != other.getG()) exit = false;
    if(width != other.getWidth()) exit = false;
    if(height != other.getHeight()) exit = false;
    if(max_iter != other.getIterMax()) exit = false;
    if(vtk_interval != other.getVtkInterval()) exit = false;
    if(tecplot_interval != other.getTecplotInterval()) exit = false;
    if(restart_interval != other.getRestartInterval()) exit = false;
    
    return exit;
}
    
///////////////////////////// PRIVATE /////////////////////////////

//=========================== OPERATIONS ===========================

void Preprocess::deduceAll(){
	calcTau();
    calcDelRho();
    calcSpacestep();
    calcTimestep();
    calcSoundspeed();
	calcNu();
    calcS2();	


}
#include"Preprocess.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Preprocess::Preprocess(double Re, double Mo, double Eo, double res, double rhol, double gamma_ini, double mu_rate, double bulk_vi, double s3, double s5, double s11, double s17, bool shear, double shear_rate, int xCells_ini, int yCells_ini ,int zCells_ini):
ReynoldsMax(Re), Morton(Mo), Eotvos(Eo),
resolution(res), rho_l(rhol), gamma(gamma_ini), mu_ratio(mu_rate), 
s_3(s3), s_5(s5), s_11(s11), s_17(s17),
isShearFlow(shear), shearRate(shear_rate),
xCells(xCells_ini), yCells(yCells_ini), zCells(zCells_ini),bulk_visco(bulk_vi)
{
    deduceAll();
}

//=========================== OPERATIONS ===========================

const ParamSet Preprocess::getParamSet()const{
    const double omega_l = 1/tau_l;
    const double omega_g = 1/tau_g;
    const double rho_r = 1;  // normalized
    const RelaxationPar3D relax(1,1,s_3,s_5,s_11,s_17);
    ParamSet param(omega_l, omega_g, rho_r, gamma, convertSigma(), g, timestep, spacestep, relax, bulk_visco);
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
    if(mu_ratio != other.getMuRatio()) exit = false;
    if(bulk_visco != other.getBulkVisco()) exit = false;
    if(s_3 != other.getS_3()) exit = false;
    if(s_5 != other.getS_5()) exit = false;
    if(s_11 != other.getS_11()) exit = false;
    if(s_17 != other.getS_17()) exit = false;
    if(isShearFlow != other.getIsShearFlow()) exit = false;
    if(shearRate != other.getShearRate()) exit = false;
    if(spacestep != other.getSpacestep()) exit = false;
    if(timestep != other.getTimestep()) exit = false;   
    if(c_s != other.getSoundspeed()) exit = false;
    if(tau_l != other.getTauL()) exit = false;
    if(tau_g != other.getTauG()) exit = false;
    if(delRho != other.getDelRho()) exit = false;
    if(nu != other.getNu()) exit = false;

    if(sigma != other.getSigma()) exit = false;
    if(g != other.getG()) exit = false;

    if(xCells != other.getXCells()) exit = false;
    if(yCells != other.getYCells()) exit = false;
    if(zCells != other.getZCells()) exit = false;

    
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
    calcTauG();
	calcNu();
    calcSigma();
    calcG();
}
#include"Preprocess_Drop.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Preprocess_Drop::Preprocess_Drop(double velo_frac, double Oh, double We, double res, double rhol, double gamma_ini, double mu_rate, double bulk_vi, double s3, double s5, double s11, double s17, int xCells_ini, int yCells_ini, int zCells_ini):
VeloFrac(velo_frac), Ohnesorge(Oh), Weber(We),
resolution(res), rho_l(rhol), gamma(gamma_ini), mu_ratio(mu_rate),
s_3(s3), s_5(s5), s_11(s11), s_17(s17),
xCells(xCells_ini), yCells(yCells_ini), zCells(zCells_ini), bulk_visco(bulk_vi)
{
    deduceAll();
}

//=========================== OPERATIONS ===========================

const ParamSet Preprocess_Drop::getParamSet()const{
    const double omega_l = 1/tau_l;
    const double omega_g = 1/tau_g;
    const double rho_r = 1;  // normalized
    const RelaxationPar3D relax(1,1,s_3,s_5,s_11,s_17);
    ParamSet param(omega_l, omega_g, rho_r, gamma, convertSigma(), 0, timestep, spacestep, relax, bulk_visco);
    return param;
}

//=========================== OPERATOR ===========================

const bool Preprocess_Drop::operator==(const Preprocess_Drop& other)const
{
    bool exit = true;

    if(VeloFrac != other.getVeloFrac()) exit = false;
    if(Ohnesorge != other.getOhnesorge()) exit = false;
    if(Weber != other.getWeber()) exit = false;
    if(resolution != other.getResolution()) exit = false;
    if(rho_l != other.getRhoL()) exit = false;
    if(gamma != other.getGamma()) exit = false;
    if(mu_ratio != other.getMuRatio()) exit = false;
    if(bulk_visco != other.getBulkVisco()) exit = false;
    if(s_3 != other.getS_3()) exit = false;
    if(s_5 != other.getS_5()) exit = false;
    if(s_11 != other.getS_11()) exit = false;
    if(s_17 != other.getS_17()) exit = false;
    if(spacestep != other.getSpacestep()) exit = false;
    if(timestep != other.getTimestep()) exit = false;   
    if(c_s != other.getSoundspeed()) exit = false;
    if(tau_l != other.getTauL()) exit = false;
    if(tau_g != other.getTauG()) exit = false;
    if(delRho != other.getDelRho()) exit = false;
    if(nu != other.getNu()) exit = false;
    if(sigma != other.getSigma()) exit = false;
    if(nu_ratio != other.getNuRatio()) exit = false;
    if(u_d != other.getUD()) exit = false;

    if(xCells != other.getXCells()) exit = false;
    if(yCells != other.getYCells()) exit = false;
    if(zCells != other.getZCells()) exit = false;

    return exit;
}
    
///////////////////////////// PRIVATE /////////////////////////////

//=========================== OPERATIONS ===========================

void Preprocess_Drop::deduceAll()
{
    calcSpacestep(); 
    calcTimestep(); 
    calcSoundspeed();
    calcNuRatio();
    calcDelRho();    
    calcUD();
    calcSigma();
    calcNu();
    calcTau();
    calcTauG();
}
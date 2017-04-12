/// preprocesses the physical properties to LB properties that are then stored in a ParamSet
/// preprocesses the physical properties to LB properties that are then stored in a ParamSet

#ifndef PREPROCESS_DROP_H
#define PREPROCESS_DROP_H

#include"ParamSet.h"
#include"Constants_basic.h"

using namespace std;

class Preprocess_Drop
{
public: 
	/// Lifecycle
	Preprocess_Drop(double velo_frac = 0.1, double Oh = 0.5, double We = 10, double res = 30, double rhol = 1, double gamma_ini = 2, double mu_rate = 1, double bulk_vi = 2, double s3 = 1, double s5 = 1, double s11 = 1, double s17 = 1, int xCells_ini = 50, int yCells_ini = 50 ,int zCells_ini = 50);

	/// operations

	// unit conversions
	inline const double convertSigma()const{return sigma / (rho_l) ;};  /// kg/s^2 -> -
	// get the parameter set
	const ParamSet getParamSet()const;

	/// accessors
	inline const double getVeloFrac()const{return VeloFrac;};
	inline const double getOhnesorge()const{return Ohnesorge;};
	inline const double getWeber()const{return Weber;};
	inline const double getResolution()const{return resolution;};
	inline const double getRhoL()const{return rho_l;};
	inline const double getGamma()const{return gamma;};
	inline const double getMuRatio()const{return mu_ratio;};
	inline const double getBulkVisco()const{return bulk_visco;};
	inline const double getS_3()const{return s_3;};
	inline const double getS_5()const{return s_5;};
	inline const double getS_11()const{return s_11;};
	inline const double getS_17()const{return s_17;};
		
	inline const double getSpacestep()const{return spacestep;};
	inline const double getTimestep()const{return timestep;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getNuRatio()const{return nu_ratio;};

	inline const double getUD()const{return u_d;};
	inline const double getTauL()const {return tau_l;};
	inline const double getTauG()const {return tau_g;};
	inline const double getDelRho()const{return delRho;};
	inline const double getNu()const{return nu;};
	inline const double getSigma()const{return sigma;};
	
	inline const int getXCells()const{return xCells;};
	inline const int getYCells()const{return yCells;};
	inline const int getZCells()const{return zCells;};


	/// operators
	const bool operator==(const Preprocess_Drop& other)const;

private:
	// given
	double VeloFrac; 		/// u_1 and u_2 as fraction of maximum velocity
	double Ohnesorge; 		/// Ohnesorge-number
	double Weber;			/// < Weber-Number
	double resolution; 		/// < width of bubble in cells
	double rho_l ;			/// < liquid density
  	double gamma; 			/// < density ratio
  	double mu_ratio;		/// viscosity ratio l/g
   	double s_3, s_5, s_11, s_17;
	
	// stored
	int xCells, yCells, zCells;
	double bulk_visco;		/// ratio of second to first viscosity mu'/mu
	
    // deduced
    double spacestep;  /// < spacestep /m	
	double timestep;   /// < timestep /s
	double c_s;
    double nu_ratio;
   	double delRho;     /// < density difference / kg * m^-3

	double u_d; 	/// speed of left and right drop
    double sigma;   	   	/// < surface tension
    double tau_l,tau_g;
    double nu;
   	
	/// operations
	inline void calcSpacestep(){spacestep = 1;}	
	inline void calcTimestep(){timestep = 1;} 
	inline void calcSoundspeed(){c_s = 1.0/sqrt(3);};
	inline void calcNuRatio(){nu_ratio = mu_ratio / gamma;};
	inline void calcDelRho(){delRho = rho_l * (1.0 - 1.0/gamma);};
	
	inline void calcUD(){u_d = 0.1 * VeloFrac * c_s;};
	inline void calcSigma(){sigma = (rho_l * resolution * 4 * u_d * u_d) / Weber ;};
	inline void calcNu(){nu = Ohnesorge * sqrt(rho_l * sigma * resolution) / rho_l;};
	inline void calcTau(){tau_l = (nu / (c_s * c_s) ) + 0.5;};
	inline void calcTauG(){tau_g = (tau_l-0.5)/nu_ratio + 0.5;};
	
	void deduceAll();
};

#endif

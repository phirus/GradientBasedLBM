/// preprocesses the physical properties to LB properties that are then stored in a ParamSet
/// preprocesses the physical properties to LB properties that are then stored in a ParamSet

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include"ParamSet.h"
#include"Constants_basic.h"

using namespace std;

class Preprocess
{
public: 
	/// Lifecycle
	Preprocess(double Re = 10, double Mo = 100, double Eo = 10, double res = 30, double rhol = 1, double gamma_ini = 2, double mu_rate = 1, double bulk_vi = 2, double s3 = 1, double s5 = 1, double s11 = 1, double s17 = 1, bool shear = false, double shear_rate = 0, int xCells_ini = 50, int yCells_ini = 50 ,int zCells_ini = 50);

	/// operations

	// unit conversions
	inline const double convertSigma()const{return sigma / (rho_l) ;};  /// kg/s^2 -> -
	// get the parameter set
	const ParamSet getParamSet()const;

	/// accessors
	inline const double getReynoldsMax()const{return ReynoldsMax;};
	inline const double getMorton()const{return Morton;};
	inline const double getEotvos()const{return Eotvos;};
	inline const double getResolution()const{return resolution;};
	inline const double getRhoL()const{return rho_l;};
	inline const double getGamma()const{return gamma;};
	inline const double getMuRatio()const{return mu_ratio;};
	inline const double getBulkVisco()const{return bulk_visco;};
	inline const double getS_3()const{return s_3;};
	inline const double getS_5()const{return s_5;};
	inline const double getS_11()const{return s_11;};
	inline const double getS_17()const{return s_17;};
	inline const bool getIsShearFlow()const{return isShearFlow;};
	inline const double getShearRate()const{return shearRate;};
	
	inline const double getSpacestep()const{return spacestep;};
	inline const double getTimestep()const{return timestep;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getTauL()const {return tau_l;};
	inline const double getTauG()const {return tau_g;};
	inline const double getDelRho()const{return delRho;};
	inline const double getNu()const{return nu;};
	inline const double getSigma()const{return sigma;};
	inline const double getG()const{return g;};

	inline const int getXCells()const{return xCells;};
	inline const int getYCells()const{return yCells;};
	inline const int getZCells()const{return zCells;};


	void setReynoldsMax(double val){ReynoldsMax = val;};

	/// operators
	const bool operator==(const Preprocess& other)const;

private:
	// given
	double ReynoldsMax ; 	/// < maximum Reynolds-Number
	double Morton;			/// < Morton-Number
	double Eotvos;			/// < Eotvos-Number
	double resolution; 		/// < width of bubble in cells
	double rho_l ;			/// < liquid density
  	double gamma; 			/// < density ratio
  	double mu_ratio;		/// viscosity ratio l/g
  	double s_3, s_5, s_11, s_17;
	bool isShearFlow; 		/// < switch that determines active shear flow
	double shearRate;		/// < shear rate in (1/s)

	// stored
	int xCells, yCells, zCells;
	double bulk_visco;		/// ratio of second to first viscosity mu'/mu
	

    // deduced
    double spacestep;  /// < spacestep /m	
	double timestep;   /// < timestep /s
	double c_s;

    double tau_l,tau_g;
   	double delRho;     /// < density difference / kg * m^-3
    double nu;

   	double sigma;   	   	/// < surface tension
	double g;          		/// < gravity / m * s^-2

	/// operations
	inline void calcTau(){tau_l = (resolution * MACH_MAX   * sqrt(3) / ReynoldsMax ) + 0.5;};
	inline void calcDelRho(){delRho = rho_l * (1.0 - 1.0/gamma);};
	inline void calcSpacestep(){spacestep = 1;}	//diameter / resolution;};
	inline void calcTimestep(){timestep = 1;} 	//spacestep / (sqrt(3) * c_s);};
	inline void calcSoundspeed(){c_s = 1.0/sqrt(3);};
	inline void calcTauG(){tau_g = (tau_l-0.5)/mu_ratio + 0.5;};
	inline void calcNu(){nu = c_s * c_s * timestep * (tau_l - 0.5);};
	inline void calcSigma(){sigma = sqrt( ( Eotvos * pow((tau_l - 0.5),4) ) / (81 * resolution * resolution * Morton)) * rho_l;};
	inline void calcG(){g = (Eotvos * sigma) / ( rho_l * (1 - 1/gamma) * resolution * resolution );};

	void deduceAll();
};

#endif

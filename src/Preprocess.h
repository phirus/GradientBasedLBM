/// preprocesses the physical properties to LB properties that are then stored in a ParamSet

#ifndef PREPROCESS_H
#define PREPROCESS_H

#include"ParamSet.h"
#include"Constants.h"

using namespace std;

class Preprocess
{
public: 
	/// Lifecycle
	Preprocess(double Re = 10, double Mo = 100, double Eo = 10, double res = 30, double rhol = 1, double gamma_ini = 2, double mu_rate = 2, double s_three = 1, double s_five = 1, int width_ini = 120, int height_ini = 360);

	/// operations
	// unit conversions
	inline const double convertSigma()const{return sigma / (rho_l) ;};  /// kg/s^2 -> -
	// inline const double convertRhoG()const{return 1/gamma;}; 
	// get the parameter set
	const ParamSet getParamSet()const;


	/// accessors
	inline const double getReynoldsMax()const{return ReynoldsMax;};
	inline const double getMorton()const{return Morton;};
	inline const double getEotvos()const{return Eotvos;};
	inline const double getResolution()const{return resolution;};
	inline const double getRhoL()const{return rho_l;};
	inline const double getGamma()const{return gamma;};
	inline const double getMuRatio()const{return muRatio;};
	inline const double getS_3()const{return s_3;};
	inline const double getS_5()const{return s_5;};
	
	inline const double getSpacestep()const{return spacestep;};
	inline const double getTimestep()const{return timestep;};
	inline const double getTau()const {return tau;};
	inline const double getDelRho()const{return delRho;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getNu()const{return nu;};
	inline const double getS2()const {return s_2;};
	inline const double getSigma()const{return sigma;};
	inline const double getG()const{return g;};

	inline const int getWidth()const{return width;};
	inline const int getHeight()const{return height;};

	void setReynoldsMax(double val){ReynoldsMax = val;};
//////FAKE
	inline double getDiameter()const{return 1;};
	inline double getSpeedlimit()const{return 1;};
	inline const double getGPhys()const{return 1;};

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
  	double muRatio;		/// ratio of second to first viscosity mu'/mu
	double s_3, s_5;

	// stored
	int width;
	int height;

    // deduced
    double spacestep;  /// < spacestep /m	
	double timestep;   /// < timestep /s

    double tau;
   	double delRho;     /// < density difference / kg * m^-3
    double c_s;
	double nu;
    double s_2;

   	double sigma;   	   	/// < surface tension
	double g;          		/// < gravity / m * s^-2

	/// operations
	inline void calcTau(){tau = (resolution * MACH_MAX   * sqrt(3) / ReynoldsMax ) + 0.5;};
	inline void calcDelRho(){delRho = rho_l * (1.0 - 1.0/gamma);};
	inline void calcSpacestep(){spacestep = 1;}	//diameter / resolution;};
	inline void calcTimestep(){timestep = 1;} 	//spacestep / (sqrt(3) * c_s);};
	inline void calcSoundspeed(){c_s = 1.0/sqrt(3);};
	inline void calcNu(){nu = c_s * c_s * timestep * (tau - 0.5);};
	inline void calcS2(){s_2 = 1.0/( (nu * muRatio) / (c_s*c_s * timestep) + 0.5);};
	inline void calcSigma(){sigma = sqrt( ( Eotvos * pow((tau - 0.5),4) ) / (81 * resolution * resolution * Morton)) * rho_l;};
	inline void calcG(){g = (Eotvos * sigma) / ( rho_l * (1 - 1/gamma) * resolution * resolution );};


	void deduceAll();
};

#endif

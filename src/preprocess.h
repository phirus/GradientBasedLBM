#ifndef PREPROCESS_H
#define PREPROCESS_H

#include"paramset.h"
#include"constants.h"

using namespace std;


/// preprocesses the physical properties to LB properties that are then stored in a ParamSet
class Preprocess
{
public: 
	// constr
	Preprocess(double Re = 50, double Mo = 1e-3, double Eo = 20, double res = 40, double rl= 1000, double gam = 5, double dia = 0.1, double soundspeed = 10, double sig= 1e-4, double grav = 10);

	// get methods
	inline const double getReynoldsMax()const{return ReynoldsMax;};
	inline const double getMorton()const{return Morton;};
	inline const double getEotvos()const{return Eotvos;};
	inline const double getResolution()const{return resolution;};
	inline const double getRhoL()const{return rho_l;};
	inline const double getGamma()const{return gamma;};
	inline const double getDiameter()const{return diameter;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getSigma()const{return sigma;};
	inline const double getGPhys()const{return g;};
	
	inline const double getTau()const {return tau;};
	inline const double getSpeedlimit()const{return speedlimit;};
	inline const double getSpacestep()const{return spacestep;};
	inline const double getTimestep()const{return timestep;};
	inline const double getNu()const{return nu;};
	inline const double getDelRho()const{return delRho;};
	
	// set methods
	void setReynoldsMax(double val){ReynoldsMax = val;};

	// unit conversions
	inline const double convertG()const{return g * timestep * timestep / spacestep;};		///  m/s^2 -> -
	inline const double convertSigma()const{return sigma * timestep * timestep / (rho_l * spacestep * spacestep * spacestep) ;};  /// kg/s^2 -> -

	inline const double convertRhoL()const{return rho_l/rho_l;};
	inline const double convertRhoG()const{return rho_g/rho_l;};

	// get the parameter set
	const ParamSet getParamSet()const;

	// refine timestep
	void refine(double factor = 1.1); // standard factor: 10%
	
	/// overloaded operators
	const bool operator==(const Preprocess& other)const;

private:
	// given
	double ReynoldsMax ; 	/// < maximum Reynolds-Number
	double Morton;			/// < Morton-Number
	double Eotvos;			/// < Eotvos-Number
	double resolution; 		/// < width of bubble in cells
	double rho_l ;			/// < liquid density
  	double gamma; 			/// < density ratio
  	double diameter;   		/// < bubble diameter /m
	double c_s; 	       	/// < speed of sound / m * s^-1
	double sigma;   	   	/// < surface tension
	double g;          		/// < gravity / m * s^-2	


    // deduced
    double tau;
    double speedlimit; /// < maximum allowed velocity
	double spacestep;  /// < spacestep /m	
	double timestep;   /// < timestep /s
	double delRho;     /// < density difference / kg * m^-3
	double rho_g;
	double nu;

	// methods
	// calculations
	inline void calcTau(){tau = (resolution * MACH_MAX   * sqrt(3) / ReynoldsMax ) + 0.5;};
	inline void calcSpeedlimit(){speedlimit = MACH_MAX * sqrt(3) * c_s ;};
	inline void calcSpacestep(){spacestep = diameter / resolution;};
	inline void calcTimestep(){timestep = spacestep / (sqrt(3) * c_s);};
	inline void calcNu(){nu = c_s * c_s * timestep * (tau - 0.5);};
	inline void calcDelRho(){delRho = rho_l * (1 - 1/gamma);};

	void deduceAll();

};

#endif

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
	Preprocess(double Re, double Mo, double Eo);



	// get methods
	inline const double getReynolds()const{return Reynolds;};
	inline const double getMorton()const{return Morton;};
	inline const double getEotvos()const{return Eotvos;};
	inline const double getResolution()const{return resolution;};
	inline const double getRhoL()const{return rho_l;};
	inline const double getRho0()const{return rho_0;};
	inline const double getGamma()const{return gamma;};

	inline const double getDiameter()const{return diameter;};
	inline const double getGPhys()const{return g;};
	inline const double getSigma()const{return sigma;};
	
	inline const double getTau()const {return tau;};
	inline const double getSpeedlimit()const{return speedlimit;};
	inline const double getTimestep()const{return timestep;};
	inline const double getSpacestep()const{return spacestep;};
	inline const double getNu()const{return nu;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getDelRho()const{return delRho;};
	
	// set methods
	inline void setReynolds(double val){Reynolds = val;};
	inline void setMorton(double val){Morton = val;};
	inline void setEotvos(double val){Eotvos = val;};
	void setResolution(double val);


	inline void setG(double val){g = val;};

	void setDiameter(double val);
	

	// calculations
	const double estimateVelocity()const;


	// unit conversions
	inline const double convertG()const{return g * timestep * timestep / spacestep;};		///  m/s^2 -> -
	inline const double convertSigma()const{return sigma * timestep * timestep / (rho_0 * spacestep * spacestep * spacestep) ;};  /// kg/s^2 -> -

	inline const double convertRhoL()const{return rho_l/rho_0;};
	inline const double convertRhoG()const{return rho_g/rho_0;};
	


private:
	// given
	double Reynolds ; 
	double Morton;
	double Eotvos;
	double resolution; // width of bubble in cells
	double rho_l ;
  	double rho_0; /// reference density
  	double gamma; 	/// < density ratio
  	double diameter;   /// < bubble diameter /m

    // deduced
    double tau;
    double speedlimit; /// < maximum allowed velocity
	double timestep;   /// < timestep /s
	double spacestep;  /// < spacestep /m	
	double nu;
	double c_s;        /// < speed of sound / m * s^-1
	double rho_g;
	double delRho;     /// < density difference / kg * m^-3

	double sigma;      /// < surface tension
	double g;          /// < gravity / m * s^-2	
	

	// methods
	// calculations
	inline void calcTau(){tau = (resolution / MACH_MAX   * sqrt(3) / Reynolds ) + 0.5;};
	inline void calcSoundspeed(){c_s = MACH_MAX * estimateVelocity();};
	inline void calcSpeedlimit(){speedlimit = MACH_MAX * sqrt(3) * c_s ;};

	inline void calcSpacestep(){spacestep = diameter / resolution;};
	inline void calcTimestep(){timestep = spacestep / (sqrt(3) * c_s);};

	inline void calcNu(){nu = c_s * c_s * timestep * (tau - 0.5);};
	inline void calcDelRho(){delRho = rho_l * (1 - 1/gamma);};
	

};

#endif

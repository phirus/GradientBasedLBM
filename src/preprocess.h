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
	const double getReynolds()const{return Reynolds;};
	const double getMorton()const{return Morton;};
	const double getEotvos()const{return Eotvos;};
	const double getResolution()const{return resolution;};
	const double getDiameter()const{return diameter;};

	// set methods
	void setReynolds(double val){Reynolds = val;};
	void setMorton(double val){Morton = val;};
	void setEotvos(double val){Eotvos = val;};
	void setG(double val){g = val;};
	void setSigma(double val){sigma = val;};
	void setDiameter(double val);
	void setResolution(double val);

	// unit conversions
	inline const double convertG()const{return g * timestep * timestep / spacestep;};		///  m/s^2 -> -
	inline const double convertSigma()const{return sigma * timestep * timestep / (rho_0 * spacestep * spacestep * spacestep) ;};  /// kg/s^2 -> -

	inline const double convertRhoR()const{return rho_r/rho_0;};
	inline const double convertRhoB()const{return rho_b/rho_0;};




private:
	// given
	double Reynolds ; 
	double Morton;
	double Eotvos;
	double g;          /// < gravity / m * s^-2
	double sigma;      /// < surface tension
    double diameter;   /// < bubble diameter /m
  	double resolution; // width of bubble in cells
  	double rho_b, rho_r;
  	double rho_0; /// reference density

    // deduced
    double tau;
    double speedlimit; /// < maximum allowed velocity
	double timestep;   /// < timestep /s
	double spacestep;  /// < spacestep /m
	
	double nu;

	double c_s;        /// < speed of sound / m * s^-1

	// methods
	// calculations
	inline void calcTau(){tau = (resolution / Mach_max   * sqrt(3) / Reynolds ) + 0.5;};
	inline void calcSpacestep(){spacestep = diameter / resolution;};
	inline void calcSpeedlimit(){speedlimit = Mach_max * sqrt(3) * c_s ;};
	inline void calcTimestep(){timestep = spacestep / (sqrt(3) * c_s);};
	inline void calcNu(){nu = c_s * c_s * timestep * (tau - 0.5);};

	void calcSoundspeed();		// < calculates the speed of sound based on a rough approximation






};

#endif

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
	inline const double getReynoldsMax()const{return ReynoldsMax;};
	inline const double getMorton()const{return Morton;};
	inline const double getEotvos()const{return Eotvos;};
	inline const double getResolution()const{return resolution;};
	inline const double getRhoL()const{return rho_l;};
	inline const double getRho0()const{return rho_0;};
	inline const double getGamma()const{return gamma;};
	inline const double getDiameter()const{return diameter;};

	inline const double getTau()const {return tau;};
	inline const double getSpeedlimit()const{return speedlimit;};
	inline const double getTimestep()const{return timestep;};
	inline const double getSpacestep()const{return spacestep;};
	inline const double getNu()const{return nu;};
	inline const double getSoundspeed()const{return c_s;};
	inline const double getDelRho()const{return delRho;};
	
	inline const double getGPhys()const{return g;};
	inline const double getSigma()const{return sigma;};


	// set methods
	void setReynoldsMax(double val){ReynoldsMax = val;};
	void setMorton(double val){Morton = val;};
	void setEotvos(double val){Eotvos = val;};
	void setResolution(double val);
	void setRhoL(double val);
	void setRho0(double val);
	void setGamma(double val);
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
	double ReynoldsMax ; 	/// < maximum Reynolds-Number
	double Morton;			/// < Morton-Number
	double Eotvos;			/// < Eotvos-Number
	double resolution; 		/// < width of bubble in cells
	double rho_l ;			/// < liquid density
  	double rho_0; 			/// < reference density
  	double gamma; 			/// < density ratio
  	double diameter;   		/// < bubble diameter /m

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
	inline void calcTau(){tau = (resolution / MACH_MAX   * sqrt(3) / ReynoldsMax ) + 0.5;};
	inline void calcSoundspeed(){c_s = MACH_MAX * estimateVelocity();};
	inline void calcSpeedlimit(){speedlimit = MACH_MAX * sqrt(3) * c_s ;};
	inline void calcSpacestep(){spacestep = diameter / resolution;};
	inline void calcTimestep(){timestep = spacestep / (sqrt(3) * c_s);};
	inline void calcNu(){nu = c_s * c_s * timestep * (tau - 0.5);};
	inline void calcDelRho(){delRho = rho_l * (1 - 1/gamma);};

	inline void calcSigma(){sigma = g * delRho *diameter*diameter / Eotvos;};
	inline void calcG(){g = Morton * pow(rho_l,2) * pow(sigma,3) / (delRho * pow((nu * rho_l),4));};
	
	void deduceAll();

};

#endif

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






	// calculations
	const double getTau()const;
	void calcSpacestep(){spacestep = diameter / resolution};


private:
	double Reynolds ; 
	double Morton;
	double Eotvos;
	double g;          /// < gravity / m * s^-2
	double sigma;      /// < surface tension
	double timestep;   /// < timestep /s
	double spacestep;  /// < spacestep /m
	double diameter;   /// < bubble diameter /m

//	double speedlimit;          /// < maximum allowed velocity
	double resolution; // width of bubble in cells

};

#endif

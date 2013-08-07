#ifndef PREPROCESS_H
#define PREPROCESS_H

#include"paramset.h"
#include"constants.h"

using namespace std;


/// collection of all parameters used during the simulation
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

	// set methods

	// calculations
	const double getTau()const;

private:
	double Reynolds ; 
	double Morton;
	double Eotvos;
	double resolution; // width of bubble in cells

};

#endif

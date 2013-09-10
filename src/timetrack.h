#ifndef TIMETRACK_H
#define TIMETRACK_H

#include<vector>
#include<cmath>

using namespace std;

/// tracks the physical time throughout the simulation
class Timetrack
{
public: 
	Timetrack(double iniTime = 1e-3, double fac = 1.1);
	inline void timestep(){count++;};
	inline void refine(){refinelist.push_back(count-refinelist.back());};
	const double getTime()const;

	// get methods
	const double getDTini()const{return dtIni;};
	const double getFactor()const{return factor;};
	const int getCount()const{return count;};
	const vector<int> getList()const{return refinelist;};

	// set methods
	void setDTini(double iniTime){dtIni = iniTime;};
	void setFactor(double fac){factor = fac;};
	void setCount(int c){count = c;};
	void setVector(const vector<int>& vec){refinelist = vec;};

private:
	double dtIni;		//   not constant, so the standard copy contr. works
	double factor;		//   not constant, so the standard copy contr. works
	int count;
	vector<int> refinelist;
};
#endif
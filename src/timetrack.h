#ifndef TIMETRACK_H
#define TIMETRACK_H

#include<vector>
#include<cmath>

using namespace std;

/// preprocesses the physical properties to LB properties that are then stored in a ParamSet
class Timetrack
{
public: 
	Timetrack(double iniTime = 1e-3, double fac = 1.1);
	inline void timestep(){count++;};
	inline void refine(){refinelist.push_back(count);};
	double getTime()const;

private:
	const double dtIni;
	const double factor;
	int count;
	vector<int> refinelist;
};

#endif

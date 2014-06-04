/// tracks the physical time and termination conditions throughout the simulation 

#ifndef TIMETRACK_H
#define TIMETRACK_H

#include<vector>
#include<cmath>

using namespace std;

class Timetrack
{
public: 
	/// Lifecycle
	Timetrack(double iniTime = 1e-3, double fac = 1.1, int t_c = 1e5, int tech = 1e3, int restart = 1e4, double resi = 1e-3);

	/// operations	
	inline void timestep(){count++;};
	inline void refine(){refinelist.push_back(count-refinelist.back());};
	const double getTime()const;
	const bool proceed()const;

	/// accessors
	inline const double getDTini()const{return dtIni;};
	inline const double getFactor()const{return factor;};
	inline const int getCount()const{return count;};
	inline const vector<int> getList()const{return refinelist;};
	inline const int getMaxCount()const{return terminalCount;};
	inline const int getTechPlotInt()const{return techplotInterval;};
	inline const int getRestartInt()const{return restartInterval;};
	inline const double getResidual()const{return relativeResidual;};

	inline void setDTini(double iniTime){dtIni = iniTime;};
	inline void setFactor(double fac){factor = fac;};
	inline void setCount(int c){count = c;};
	inline void setVector(const vector<int>& vec){refinelist = vec;};
	inline void setMaxCount(int t_c){terminalCount = t_c;};
	inline void setTechPlotInt(int tmp){techplotInterval = tmp;};
	inline void setRestartInt(int tmp){restartInterval = tmp;};
	inline void setResidual(double tmp){relativeResidual = tmp;};

	/// operators
	const bool operator==(const Timetrack& other)const;

private:
	double dtIni;		//   not constant, so the standard copy contr. works
	double factor;		//   not constant, so the standard copy contr. works
	int count;
	vector<int> refinelist;
	// termination conditions
	int terminalCount;			// maximum number of time steps
	int techplotInterval;
	int restartInterval;
	double relativeResidual;
};
#endif
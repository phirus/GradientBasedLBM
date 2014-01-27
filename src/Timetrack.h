#ifndef TIMETRACK_H
#define TIMETRACK_H

#include<vector>
#include<cmath>

using namespace std;

/// tracks the physical time throughout the simulation
class Timetrack
{
public: 
	Timetrack(double iniTime = 1e-3, double fac = 1.1, int t_c = 1e5, double t_t = 5, int tech = 1e3, int rest = 1e4);
	inline void timestep(){count++;};
	inline void refine(){refinelist.push_back(count-refinelist.back());};
	const double getTime()const;
	const bool proceed()const;

	/// get methods
	inline const double getDTini()const{return dtIni;};
	inline const double getFactor()const{return factor;};
	inline const int getCount()const{return count;};
	inline const vector<int> getList()const{return refinelist;};
	inline const int getMaxCount()const{return terminal_count;};
	inline const double getMaxTime()const{return terminal_time;};
	inline const int getTechPlotInt()const{return techplot_interval;};
	inline const int getRestartInt()const{return restart_interval;};

	/// set methods
	inline void setDTini(double iniTime){dtIni = iniTime;};
	inline void setFactor(double fac){factor = fac;};
	inline void setCount(int c){count = c;};
	inline void setVector(const vector<int>& vec){refinelist = vec;};
	inline void setMaxCount(int t_c){terminal_count = t_c;};
	inline void setMaxTime(double t){terminal_time = t;};
	inline void setTechPlotInt(int tmp){techplot_interval = tmp;};
	inline void setRestartInt(int tmp){restart_interval = tmp;};

	/// overloaded operators
	const bool operator==(const Timetrack& other)const;


private:
	double dtIni;		//   not constant, so the standard copy contr. works
	double factor;		//   not constant, so the standard copy contr. works
	int count;
	vector<int> refinelist;
	// termination conditions
	int terminal_count;			// maximum number of time steps
	double terminal_time;		// maximum simulation time
	int techplot_interval;
	int restart_interval;
};
#endif
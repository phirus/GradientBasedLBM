#include"timetrack.h"

Timetrack::Timetrack(double iniTime, double fac):
dtIni(iniTime)
,factor(fac)
,count(0)
,refinelist(1,0)
{}

const double Timetrack::getTime()const{
	double time = 0;
	int j = 0;	
	for(unsigned int i = 1; i<refinelist.size() ; i++){
		time += refinelist[i] * dtIni * pow(factor,i-1);
		j += refinelist[i];
	}
	time += (count - j) * dtIni * pow(factor,refinelist.size()-1);
	return time;
}
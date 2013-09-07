#include"timetrack.h"

Timetrack::Timetrack(double iniTime, double fac):
dtIni(iniTime)
,factor(fac)
,count(0)
,refinelist(1,0)
{}

double Timetrack::getTime()const{
	double time(0);	
	for(unsigned int i = 0; i< refinelist.size(); i++){
		time += refinelist[i] * pow(factor,i);
	}
	return time;
}
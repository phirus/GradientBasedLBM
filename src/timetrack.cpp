#include"timetrack.h"

Timetrack::Timetrack(double iniTime, double fac, int t_c, double t_t):
dtIni(iniTime)
,factor(fac)
,count(0)
,refinelist(1,0)
,terminal_count(t_c)
,terminal_time(t_t)
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

const bool Timetrack::operator==(const Timetrack& other)const
{
    bool exit = true;
    if(dtIni != other.getDTini()) exit = false;
    if(factor != other.getFactor()) exit = false;
    if(count != other.getCount()) exit = false;
    if(refinelist != other.getList()) exit = false;
    if(terminal_count != other.getMaxCount()) exit = false;
    if(terminal_time != other.getMaxTime()) exit = false;

    return exit;
}

const bool Timetrack::proceed()const
{
	bool getGoing = true;
	if(count > terminal_count) getGoing = false;
	if(getTime() > terminal_time) getGoing = false;

	return getGoing;
}
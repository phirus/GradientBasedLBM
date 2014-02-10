#include"Timetrack.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Timetrack::Timetrack(double iniTime, double fac, int t_c, double t_t, int tech, int rest):
dtIni(iniTime)
,factor(fac)
,count(0)
,refinelist(1,0)
,terminalCount(t_c)
,terminalTime(t_t)
,techplotInterval(tech)
,restartInterval(rest)
{}

//=========================== OPERATIONS ===========================

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

const bool Timetrack::proceed()const
{
	bool getGoing = true;
	if(count > terminalCount) getGoing = false;
	if(getTime() > terminalTime) getGoing = false;

	return getGoing;
}

//=========================== OPERATORS ===========================

const bool Timetrack::operator==(const Timetrack& other)const
{
    bool exit = true;
    if(dtIni != other.getDTini()) exit = false;
    if(factor != other.getFactor()) exit = false;
    if(count != other.getCount()) exit = false;
    if(refinelist != other.getList()) exit = false;
    if(terminalCount != other.getMaxCount()) exit = false;
    if(terminalTime != other.getMaxTime()) exit = false;
    if(techplotInterval != other.getTechPlotInt()) exit = false;
    if(restartInterval != other.getRestartInt()) exit = false;

    return exit;
}
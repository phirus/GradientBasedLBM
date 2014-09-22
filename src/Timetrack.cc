#include"Timetrack.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Timetrack::Timetrack(int t_c, int tech, int restart):
count(0)
,terminalCount(t_c)
,techplotInterval(tech)
,restartInterval(restart)
{}

//=========================== OPERATIONS ===========================

const bool Timetrack::proceed()const
{
	bool getGoing = true;
	if(count > terminalCount) getGoing = false;

	return getGoing;
}

//=========================== OPERATORS ===========================

const bool Timetrack::operator==(const Timetrack& other)const
{
    bool exit = true;
    if(count != other.getCount()) exit = false;
    if(terminalCount != other.getMaxCount()) exit = false;
    if(techplotInterval != other.getTechPlotInt()) exit = false;
    if(restartInterval != other.getRestartInt()) exit = false;
   
    return exit;
}
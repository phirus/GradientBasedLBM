#include"Timetrack.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Timetrack::Timetrack(int t_c, int out, int restart):
count(0)
,terminalCount(t_c)
,outputInterval(out)
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
    if(outputInterval != other.getOutputInt()) exit = false;
    if(restartInterval != other.getRestartInt()) exit = false;
   
    return exit;
}
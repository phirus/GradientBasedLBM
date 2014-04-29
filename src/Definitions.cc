#include"Definitions.h"

///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const array array_diff(const array &one, const array &two)
{
    array a;
    for (int i=0; i<9; i++)
    {
        a[i] = one[i]-two[i];
    }
    return a;
}

const array array_add(const array &one, const array &two)
{
    array a;
    for (int i=0; i<9; i++)
    {
        a[i] = one[i]+two[i];
    }
    return a;
}

const array array_times(const array &foo, double factor)
{
    array bar = foo;
    for(boost::array<double,9>::iterator it = bar.begin(); it != bar.end(); ++it)
    {
        *it *= factor;
    }
    return bar;
}

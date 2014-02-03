#include"Definitions.h"

const array array_diff(const array &one, const array &two)
{
    array a;
    for (int i=0; i<9; i++)
    {
        a[i] = one[i]-two[i];
    }
    return a;
}

const DistributionSetType distro_diff(const DistributionSetType &one, const DistributionSetType &two)
{
    DistributionSetType diff;
    diff[0] = array_diff(one[0],two[0]);
    diff[1] = array_diff(one[1],two[1]);
    return diff;
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
        *it = *it * factor;
    }
    return bar;
}
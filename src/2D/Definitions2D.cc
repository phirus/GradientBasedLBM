#include"Definitions2D.h"

///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const array2D array_diff_2D(const array2D &one, const array2D &two)
{
    array2D a;
    for (int i=0; i<9; i++)
    {
        a[i] = one[i]-two[i];
    }
    return a;
}

const array2D array_add_2D(const array2D &one, const array2D &two)
{
    array2D a;
    for (int i=0; i<9; i++)
    {
        a[i] = one[i]+two[i];
    }
    return a;
}

const array2D array_times_2D(const array2D &foo, double factor)
{
    array2D bar = foo;
    for(boost::array<double,9>::iterator it = bar.begin(); it != bar.end(); ++it)
    {
        *it *= factor;
    }
    return bar;
}

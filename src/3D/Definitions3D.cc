#include"Definitions3D.h"

///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const array3D array_diff_3D(const array3D &one, const array3D &two)
{
    array3D a;
    for (int i=0; i<19; i++)
    {
        a[i] = one[i]-two[i];
    }
    return a;
}

const array3D array_add_3D(const array3D &one, const array3D &two)
{
    array3D a;
    for (int i=0; i<19; i++)
    {
        a[i] = one[i]+two[i];
    }
    return a;
}

const array3D array_times_3D(const array3D &foo, double factor)
{
    array3D bar = foo;
    for(boost::array<double,19>::iterator it = bar.begin(); it != bar.end(); ++it)
    {
        *it *= factor;
    }
    return bar;
}

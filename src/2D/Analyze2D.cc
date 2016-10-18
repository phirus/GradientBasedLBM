#include"Analyze2D.h"

const double getLineShearSum(const Lattice2D& l)
{
    field2D data = l.getData();
    const DimSet2D extent = l.getSize();
    
    const int y_m = extent[1] /2;

    double sum = 0;
    Cell2D tmp_cell;
    VeloSet2D tmp_velo;

    for (int x = 0; x<extent[0];x++)
    {
        tmp_cell = data[x][y_m];
        tmp_cell.calcRho();
        tmp_velo = tmp_cell.getU();
        sum += tmp_velo[0].y;
    }

    return sum;
}
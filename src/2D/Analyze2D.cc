#include"Analyze2D.h"

const Vector2D getBubbleVelocity(const Lattice2D& l)
{
    const field2D data = l.getData();
    const DimSet2D extent = l.getSize();

    Cell2D tmp_cell;
    VeloSet2D tmp_velo;
    ColSet tmp_rho;
    Vector2D momentum(0,0);
    double rho_sum(0);

    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            tmp_cell = data[x][y];
            tmp_cell.calcRho();
            tmp_velo = tmp_cell.getU();
            tmp_rho = tmp_cell.getRho();
            if ( tmp_cell.calcPsi() < -0.99) {
                momentum = momentum + (tmp_velo[1] * tmp_rho[1]);
                rho_sum += tmp_rho[1];
            }
        }
    }

    const Vector2D velocity = momentum * (1.0/ rho_sum);
    return velocity;
}

const Vector2D getBubblePosition(const Lattice2D& l)
{
    const field2D data = l.getData();
    const DimSet2D extent = l.getSize();

    Cell2D tmp_cell;
    ColSet tmp_rho;
    Vector2D position(0,0);
    double rho_sum(0);

    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            tmp_cell = data[x][y];
            tmp_cell.calcRho();
            tmp_rho = tmp_cell.getRho();
            if ( tmp_cell.calcPsi() < -0.99) {
                position.x = position.x + (x * tmp_rho[1]);
                position.y = position.y + (y * tmp_rho[1]);
                rho_sum += tmp_rho[1];
            }
        }
    }

    position = position * (1.0/ rho_sum);
    return position;
}

const double getReynolds(const Lattice2D& l, double resolution)
{
    const ParamSet params = l.getParams();
    const Vector2D velocity = getBubbleVelocity(l);
    const double reynolds = getReynolds(params, velocity.y, resolution);

    return reynolds;
}

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
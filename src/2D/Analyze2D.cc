#include"Analyze2D.h"

const double getBubbleVelocity(const Lattice2D& l)
{
    field2D data = l.getData();
    ColSet extent = l.getSize();

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

    Vector2D velocity = momentum * (1.0/ rho_sum);

    double upward_velo = velocity.y;
    return upward_velo;
}

const double getReynolds(const Lattice2D& l, double resolution)
{
    const ParamSet params = l.getParams();
    const double velocity = getBubbleVelocity(l);
    const double reynolds = getReynolds(params, velocity, resolution);

    return reynolds;
}
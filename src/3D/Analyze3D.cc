#include"Analyze3D.h"

const Vector3D getBubbleVelocity(const Lattice3D& l)
{
    const field3D data = l.getData();
    const DimSet3D extent = l.getSize();

    Cell3D tmp_cell;
    VeloSet3D tmp_velo;
    ColSet tmp_rho;
    Vector3D momentum(0,0,0);
    double rho_sum(0);

    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            for (int z = 0; z<extent[2];z++)
            {
                tmp_cell = data[x][y][z];
                tmp_cell.calcRho();
                tmp_velo = tmp_cell.getU();
                tmp_rho = tmp_cell.getRho();
                if ( tmp_cell.calcPsi() < -0.99) 
                {
                    momentum = momentum + (tmp_velo[1] * tmp_rho[1]);
                    rho_sum += tmp_rho[1];
                }
            }
         }
    }

    Vector3D velocity = momentum * (1.0/ rho_sum);
    return velocity;
}

const Vector3D getBubblePosition(const Lattice3D& l)
{
    const field3D data = l.getData();
    const DimSet3D extent = l.getSize();

    Cell3D tmp_cell;
    ColSet tmp_rho;
    Vector3D position(0,0,0);
    double rho_sum(0);

    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            for (int z = 0; z<extent[2];z++)
            {
                tmp_cell = data[x][y][z];
                tmp_cell.calcRho();
                tmp_rho = tmp_cell.getRho();
                if ( tmp_cell.calcPsi() < -0.99) 
                {
                    position.x = position.x + (x * tmp_rho[1]);
                    position.y = position.y + (y * tmp_rho[1]);
                    position.z = position.z + (z * tmp_rho[1]);
                    rho_sum += tmp_rho[1];
                }
            }
         }
    }

    position = position * (1.0/ rho_sum);
    return position;
}

const double getReynolds(const Lattice3D& l, double resolution)
{
    const ParamSet params = l.getParams();
    const Vector3D velocity = getBubbleVelocity(l);
    const double reynolds = getReynolds(params, velocity.z, resolution);

    return reynolds;
}

const double getLineShearSum(const Lattice3D& l)
{
    field3D data = l.getData();
    const DimSet3D extent = l.getSize();
    
    const int y_m = extent[1] /2;
    const int z_m = extent[2] /2;

    double sum = 0;
    Cell3D tmp_cell;
    VeloSet3D tmp_velo;

    for (int x = 0; x<extent[0];x++)
    {
        tmp_cell = data[x][y_m][z_m];
        tmp_cell.calcRho();
        tmp_velo = tmp_cell.getU();
        sum += tmp_velo[0].z;
    }

    return sum;
}
#include"Analyze.h"

const double getBubbleVelocity(const Lattice& l)
{
    field data = l.getData();
    ColSet extent = l.getSize();

    Cell tmp_cell;
    VeloSet tmp_velo;
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

const double getEotvos(const ParamSet& params, double resolution)
{
    // const double sigma = params.getSigma();
    const double g = params.getG();
    const double gamma = params.getGamma();

    const double tau = 1.0 / params.getOmegaRed();
    const double dt = params.getDeltaT();
    const ColSet A_k = params.getAk(params.getOmegaRed());
    const double sigma_alt = (2.0/9.0) * (A_k[0] + A_k[1]) *dt * tau; 

    const double eotvos = (g * (1.0 - 1.0/gamma) * resolution * resolution) / sigma_alt;
    return eotvos;
}

const double getMorton(const ParamSet& params)
{
    // const double sigma = params.getSigma();
    const double g = params.getG();
    const double gamma = params.getGamma();

    const double tau = 1.0 / params.getOmegaRed();
    const double dt = params.getDeltaT();
    const ColSet A_k = params.getAk(params.getOmegaRed());
    const double sigma_alt = (2.0/9.0) * (A_k[0] + A_k[1]) *dt * tau;
    
    const double morton = (g * pow((tau - 0.5),4) * (1.0 - 1.0/gamma) ) / (sigma_alt * sigma_alt * sigma_alt);
    return morton;
}

const double getReynolds(const ParamSet& params, double velocity, double resolution)
{
    const double tau = 1.0 / params.getOmegaRed();
    const double reynolds = (3 * velocity * resolution) / (tau - 0.5);
    return reynolds;
}

const double getReynolds(const Lattice& l, double resolution)
{
    const ParamSet params = l.getParams();
    const double velocity = getBubbleVelocity(l);
    const double reynolds = getReynolds(params, velocity, resolution);

    return reynolds;
}


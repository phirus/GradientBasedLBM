#include"Analyze.h"

const double getBubbleVelocity(const Lattice& l)
{
    field data = l.getData();
    ColSet extent = l.getSize();

    VeloSet tmp_velo;
    ColSet tmp_rho;
    Vector momentum(0,0);
    double rho_sum(0);

    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            tmp_velo = data[x][y].getU();
            tmp_rho = data[x][y].getRho();
            momentum = momentum + (tmp_velo[1] * tmp_rho[1]);
            rho_sum += tmp_rho[1];
        }
    }

    Vector velocity = momentum * (1.0/ rho_sum);

    double upward_velo = velocity.y;
    return upward_velo;
}

const double getEotvos(const ParamSet& params, double resolution)
{
    const double sigma = params.getSigma();
    const double g = params.getG();
    const double gamma = params.getGamma();
    const double eotvos = (g * (1.0 - 1.0/gamma) * resolution * resolution) / sigma;
    return eotvos;
}

const double getMorton(const ParamSet& params)
{
    const double sigma = params.getSigma();
    const double g = params.getG();
    const double gamma = params.getGamma();
    const double tau = 1.0 / params.getOmegaRed();
    const double morton = (g * pow((tau - 0.5),4) * (1.0 - 1.0/gamma) ) / (sigma * sigma * sigma);
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
#include"Analyze.h"

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

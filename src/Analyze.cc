#include"Analyze.h"

const double getEotvos(const ParamSet& params, double resolution)
{
    const double sigma = params.getSigma();
    const double g = params.getG();
    const double gamma = params.getGamma();
    const double eotvos = (g * (1.0 - 1.0/gamma) * resolution * resolution) / sigma;
    return eotvos;
}



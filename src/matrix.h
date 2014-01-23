#ifndef MATRIX_H
#define MATRIX_H

#include"definitions.h"

// Git Test

/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

class Matrix
{
    public:
    Matrix(boost::multi_array<double,2> m);
    Matrix(RelaxationPar relax, double omega = 1);

    void resetOmega(double omega);

    const array operator*(const array &other) const;
    const double linewise(const array &oher, int line) const;

    private:
    boost::multi_array<double,2> matrix;

};

#endif

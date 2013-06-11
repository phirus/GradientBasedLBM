#ifndef MATRIX_H
#define MATRIX_H

#include"definitions.h"

// Git Test

/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

class Matrix
{
    public:
    Matrix(bool standard = true);
    Matrix(double s_2, double s_3, double s_5, bool momentum = false);
    const array operator*(const array &other) const;
    const double linewise(const array &oher, int line) const;
    void addOmega(double omega);

    private:
    boost::multi_array<double,2> matrix;

};

#endif

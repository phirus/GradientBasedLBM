#ifndef MATRIX_H
#define MATRIX_H

#include"definitions.h"

// Git Test

/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

class Matrix
{
    public:
        Matrix(bool identity = false);
        Matrix(const boost::multi_array<double,2> &m);
        Matrix(const RelaxationPar &relax, double omega = 1);
        Matrix(const Matrix &other);

        void resetOmega(double omega);

        // get
        inline const boost::multi_array<double,2> getData()const{return matrix;};

        // multiplication
        const array operator*(const array &other) const;
        const double linewise(const array &oher, int line) const;

        const Matrix operator*(double other)const;

        // addition
        const Matrix operator+(const Matrix &other)const;
        const Matrix operator-(const Matrix &other)const;

        // is equal?
        const bool operator==(const Matrix &other)const; 

    private:
        boost::multi_array<double,2> matrix;
};

#endif

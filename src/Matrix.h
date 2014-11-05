/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

#ifndef MATRIX_H
#define MATRIX_H

#include"2D/Definitions2D.h"

class Matrix
{
    public:
        /// Lifecylce
        Matrix(bool identity = false);
        Matrix(const boost::multi_array<double,2> &m);
        Matrix(const RelaxationPar &relax, double omega = 1);
        Matrix(const Matrix &other);

         /// operators
        const array2D operator*(const array2D &other) const;
        const DistributionSetType2D operator*(const DistributionSetType2D &other) const;
        const Matrix operator*(double other)const;
        const Matrix operator+(const Matrix &other)const;
        const Matrix operator-(const Matrix &other)const;
        const bool operator==(const Matrix &other)const;       
        
        /// operations
        const double linewise(const array2D &other, int line) const;
        
        /// accessors
        inline const boost::multi_array<double,2> getData()const{return matrix;};
        void resetOmega(double omega);

    private:
        boost::multi_array<double,2> matrix;
};

#endif

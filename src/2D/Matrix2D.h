/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

#ifndef MATRIX2D_H
#define MATRIX2D_H

#include"Definitions2D.h"

class Matrix2D
{
    public:
        /// Lifecylce
        Matrix2D(bool identity = false);
        Matrix2D(const boost::multi_array<double,2> &m);
        Matrix2D(const RelaxationPar &relax, double omega = 1);
        Matrix2D(const Matrix2D &other);

         /// operators
        const array2D operator*(const array2D &other) const;
        const DistributionSetType2D operator*(const DistributionSetType2D &other) const;
        const Matrix2D operator*(double other)const;
        const Matrix2D operator+(const Matrix2D &other)const;
        const Matrix2D operator-(const Matrix2D &other)const;
        const bool operator==(const Matrix2D &other)const;       
        
        /// operations
        const double linewise(const array2D &other, int line) const;
        
        /// accessors
        inline const boost::multi_array<double,2> getData()const{return matrix;};
        void resetOmega(double omega);

    private:
        boost::multi_array<double,2> matrix;
};

#endif
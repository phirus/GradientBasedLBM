/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

#ifndef MATRIX3D_H
#define MATRIX3D_H

#include"Definitions3D.h"

class Matrix3D
{
    public:
        /// Lifecylce
        Matrix3D(bool identity = false);
        Matrix3D(double diag);
        Matrix3D(const boost::multi_array<double,2> &m);
        Matrix3D(const RelaxationPar3D &relax, bool forcingterm);
        Matrix3D(const RelaxationPar3D &relax);
        Matrix3D(const Matrix3D &other);

         /// operators
        const Matrix3D operator*(const Matrix3D &other)const;
        const array3D operator*(const array3D &other) const;
        const DistributionSetType3D operator*(const DistributionSetType3D &other) const;
        const Matrix3D operator*(double other)const;
        const Matrix3D operator+(const Matrix3D &other)const;
        const Matrix3D operator-(const Matrix3D &other)const;
        const bool operator==(const Matrix3D &other)const;       
        
        /// operations
        const array3D diagMult(const array3D &other) const;
        const DistributionSetType3D diagMult(const DistributionSetType3D &other) const;
        const double linewise(const array3D &other, int line) const;
        
        /// accessors

        inline const boost::multi_array<double,2> getData()const{return matrix;};
        void resetOmega(double omega);

    private:
        boost::multi_array<double,2> matrix;
};

#endif

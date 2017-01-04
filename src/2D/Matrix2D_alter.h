/// The Matrix-Class is made for the MRT-Colission-Step
/** The Matrix-Class is made for the MRT-Colission-Step , its only operation is multiplication with an 9-entry-array */

#ifndef MATRIX2D_ALTER_H
#define MATRIX2D_ALTER_H

#include"Definitions2D.h"
#include"Eigen/Dense"


class Matrix2D_alter
{
    public:
        /// Lifecylce
        Matrix2D_alter(const Eigen::Matrix<float,9,9> &other):m(other){};
        Matrix2D_alter(const Matrix2D_alter &other):m(other.getData()){};
        Matrix2D_alter(bool identity = false);
        Matrix2D_alter(const RelaxationPar2D &relax, float omega = 1);
        

         /// operators
        inline const Matrix2D_alter operator*(float other)const{return Matrix2D_alter(m* other);};
        const Matrix2D_alter operator+(const Matrix2D_alter &other)const{return Matrix2D_alter(m + other.getData());};
        const Matrix2D_alter operator-(const Matrix2D_alter &other)const{return Matrix2D_alter(m - other.getData());};
        const bool operator==(const Matrix2D_alter &other)const{return (m == other.getData());};       
       
       const array2D operator*(const array2D &other) const;
        const DistributionSetType2D operator*(const DistributionSetType2D &other) const;

        /// accessors
        inline const Eigen::Matrix<float,9,9> getData()const{return m;};
        void resetOmega(float omega);

    private:
        Eigen::Matrix<float,9,9> m;
};

const Eigen::Matrix<float,1,9> boost2Vector2D(array2D);


#endif

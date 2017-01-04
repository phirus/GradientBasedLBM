#include"Matrix2D_alter.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Matrix2D_alter::Matrix2D_alter(bool identity)
{
    if (identity == true) m.setIdentity(); 
    else m.setZero();
}

Matrix2D_alter::Matrix2D_alter(const RelaxationPar2D &relax, float omega)
{
    m.setZero();

    m(0,0) = 1;
    m(1,1) = relax.s_2;
    m(2,2) = relax.s_3;
    m(3,3) = 1;
    m(4,4) = relax.s_5;
    m(5,5) = 1;
    m(6,6) = relax.s_5;
    m(7,7) = omega;
    m(8,8) = omega;
}

//=========================== OPERATORS ===========================

const array2D Matrix2D_alter::operator*(const array2D &other)const {
    array2D a;

    for (int i= 0; i<9;i++){
        float sum = 0;
        for(int j = 0; j<9;j++){
            sum += other[j] * m(i,j);
        }
        a[i] = sum;
    }
    return a;



    // array2D a;

    // for (int i= 0; i<9;i++){
    //     float sum = 0;
    //     for(int j = 0; j<9;j++){
    //         sum += other[j] * m(i,j);
    //     }
    //     a[i] = sum;
    // }
    // return a;
}

const DistributionSetType2D Matrix2D_alter::operator*(const DistributionSetType2D &other)const {
    DistributionSetType2D a;

    a[0] = *this * other[0];
    a[1] = *this * other[1];
 
    return a;
}

//=========================== ACCESSORS ===========================

void Matrix2D_alter::resetOmega(float omega){
    m(7,7) = omega;
    m(8,8) = omega;
}
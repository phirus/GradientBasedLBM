#include"Matrix2D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Matrix2D::Matrix2D(bool identity):matrix(boost::extents[9][9])
{
    if (identity == true) {
        for(int i = 0;i<9;i++){
            for(int j=0;j<9;j++){
                matrix[i][j] = (i == j) ? 1 : 0;
            }
        }
    }
    else {
        for(int i = 0;i<9;i++){
            for(int j=0;j<9;j++){
                matrix[i][j] = 0;
            }
        }
    }
}

Matrix2D::Matrix2D(double diag):matrix(boost::extents[9][9])
{
    for(int i = 0;i<9;i++)
    {
        for(int j=0;j<9;j++)
        {
            matrix[i][j] = (i == j) ? diag : 0;
        }
    }
}


Matrix2D::Matrix2D(const boost::multi_array<double,2> &m):matrix(boost::extents[9][9])
{
    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            matrix[i][j] = m[i][j];
        }
    }
}


Matrix2D::Matrix2D(const RelaxationPar2D &relax, bool forcingterm, double omega):matrix(boost::extents[9][9])
{
    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            matrix[i][j] = 0;
        }
    }

    if (forcingterm == false)
    {
        matrix[0][0] = 1;
        matrix[1][1] = relax.s_2;
        matrix[2][2] = relax.s_3;
        matrix[3][3] = 1;
        matrix[4][4] = relax.s_5;
        matrix[5][5] = 1;
        matrix[6][6] = relax.s_5;
        matrix[7][7] = omega;
        matrix[8][8] = omega;

    }
    else
    {
        matrix[0][0] = 0.5;
        matrix[1][1] = 1 - 0.5 * relax.s_2;
        matrix[2][2] = 1 - 0.5 * relax.s_3;
        matrix[3][3] = 0.5;
        matrix[4][4] = 1 - 0.5 * relax.s_5;
        matrix[5][5] = 0.5;
        matrix[6][6] = 1 - 0.5 * relax.s_5;
        matrix[7][7] = 1 - 0.5 * omega;
        matrix[8][8] = 1 - 0.5 * omega;

    }
}

Matrix2D::Matrix2D(const Matrix2D &other):matrix(boost::extents[9][9])
{
    matrix = other.getData();
}

//=========================== OPERATORS ===========================
const Matrix2D Matrix2D::operator*(const Matrix2D &other)const
{
    boost::multi_array<double,2> m(boost::extents[9][9]);
    boost::multi_array<double,2> mother = other.getData();
    for (int i=0; i<9; i++)
    {
        for (int j=0; j<9; j++)
        {
            double sum = 0;
            for (int k = 0; k<9; k++)
            {
                sum += matrix[i][k] * mother[k][j];
            }
            m[i][j] = sum;
        }
    }
    return Matrix2D(m);
}

const array2D Matrix2D::operator*(const array2D &other)const 
{
    array2D a;

    for (int i= 0; i<9;i++){
        double sum = 0;
        for(int j = 0; j<9;j++){
            sum += other[j] * matrix[i][j];
        }
        a[i] = sum;
    }
    return a;
}

const DistributionSetType2D Matrix2D::operator*(const DistributionSetType2D &other)const 
{
    DistributionSetType2D a;

    a[0] = *this * other[0];
    a[1] = *this * other[1];
 
    return a;
}

const Matrix2D Matrix2D::operator*(double other)const{
    boost::multi_array<double,2> m(boost::extents[9][9]);

    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            m[i][j] = matrix[i][j] * other;
        }
    }
    return Matrix2D(m);
}

const Matrix2D Matrix2D::operator+(const Matrix2D &other)const{
    boost::multi_array<double,2> m(boost::extents[9][9]); 
    boost::multi_array<double,2> mother = other.getData();

    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            m[i][j] = matrix[i][j] + mother[i][j];
        }
    }
    return Matrix2D(m);
}

const Matrix2D Matrix2D::operator-(const Matrix2D &other)const{
    boost::multi_array<double,2> m(boost::extents[9][9]); 
    boost::multi_array<double,2> mother = other.getData();

    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            m[i][j] = matrix[i][j] - mother[i][j];
        }
    }
    return Matrix2D(m);
}

const bool Matrix2D::operator==(const Matrix2D &other)const{
    boost::multi_array<double,2> mother = other.getData();
    bool equal = false;
    if (matrix == mother) equal = true;

    return equal;
}

//=========================== OPERATIONS ===========================

const array2D Matrix2D::diagMult(const array2D &other) const
{
    array2D a;
    for (int i= 0; i<9;i++){
        a[i] = other[i] * matrix[i][i];
    }
    return a;
}

const DistributionSetType2D Matrix2D::diagMult(const DistributionSetType2D &other) const
{
    DistributionSetType2D a;

    a[0] = diagMult(other[0]);
    a[1] = diagMult(other[1]);
 
    return a;
}


// const double Matrix2D::linewise(const array2D &other, int line)const{
//     double sum = 0;
//     for(int j = 0; j<9;j++){
//         sum += other[j] * matrix[line][j];
//     }
//     return sum;
// }

//=========================== ACCESSORS ===========================

void Matrix2D::resetOmega(double omega){
    matrix[7][7] = omega;
    matrix[8][8] = omega;
}
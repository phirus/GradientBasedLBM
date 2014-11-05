#include"Matrix.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Matrix::Matrix(bool identity):matrix(boost::extents[9][9])
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

Matrix::Matrix(const boost::multi_array<double,2> &m):matrix(boost::extents[9][9])
{
    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            matrix[i][j] = m[i][j];
        }
    }
}


Matrix::Matrix(const RelaxationPar &relax, double omega):matrix(boost::extents[9][9])
{
    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            matrix[i][j] = 0;
        }
    }
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

Matrix::Matrix(const Matrix &other):matrix(boost::extents[9][9])
{
    matrix = other.getData();
}

//=========================== OPERATORS ===========================

const array2D Matrix::operator*(const array2D &other)const {
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

const DistributionSetType Matrix::operator*(const DistributionSetType &other)const {
    DistributionSetType a;

    a[0] = *this * other[0];
    a[1] = *this * other[1];
 
    return a;
}

const Matrix Matrix::operator*(double other)const{
    boost::multi_array<double,2> m(boost::extents[9][9]);

    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            m[i][j] = matrix[i][j] * other;
        }
    }
    return Matrix(m);
}

const Matrix Matrix::operator+(const Matrix &other)const{
    boost::multi_array<double,2> m(boost::extents[9][9]); 
    boost::multi_array<double,2> mother = other.getData();

    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            m[i][j] = matrix[i][j] + mother[i][j];
        }
    }
    return Matrix(m);
}

const Matrix Matrix::operator-(const Matrix &other)const{
    boost::multi_array<double,2> m(boost::extents[9][9]); 
    boost::multi_array<double,2> mother = other.getData();

    for(int i = 0;i<9;i++){
        for(int j=0;j<9;j++){
            m[i][j] = matrix[i][j] - mother[i][j];
        }
    }
    return Matrix(m);
}

const bool Matrix::operator==(const Matrix &other)const{
    boost::multi_array<double,2> mother = other.getData();
    bool equal = false;
    if (matrix == mother) equal = true;

    return equal;
}

//=========================== OPERATIONS ===========================

const double Matrix::linewise(const array2D &other, int line)const{
    double sum = 0;
    for(int j = 0; j<9;j++){
        sum += other[j] * matrix[line][j];
    }
    return sum;
}

//=========================== ACCESSORS ===========================

void Matrix::resetOmega(double omega){
    matrix[7][7] = omega;
    matrix[8][8] = omega;
}
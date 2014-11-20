#include"Matrix3D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Matrix3D::Matrix3D(bool identity):matrix(boost::extents[19][19])
{
    if (identity == true) {
        for(int i = 0;i<19;i++){
            for(int j=0;j<19;j++){
                matrix[i][j] = (i == j) ? 1 : 0;
            }
        }
    }
    else {
        for(int i = 0;i<19;i++){
            for(int j=0;j<19;j++){
                matrix[i][j] = 0;
            }
        }
    }
}

Matrix3D::Matrix3D(const boost::multi_array<double,2> &m):matrix(boost::extents[19][19])
{
    for(int i = 0;i<19;i++){
        for(int j=0;j<19;j++){
            matrix[i][j] = m[i][j];
        }
    }
}


Matrix3D::Matrix3D(const RelaxationPar3D &relax, double omega):matrix(boost::extents[19][19])
{
    for(int i = 0;i<19;i++){
        for(int j=0;j<19;j++){
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
    matrix[7][7] = 1;
    matrix[8][8] = relax.s_5;

    matrix[9][9] = omega;
    matrix[10][10] = relax.s_11;
    matrix[11][11] = omega;
    matrix[12][12] = relax.s_11;
    matrix[13][13] = omega;
    matrix[14][14] = omega;
    matrix[15][15] = omega;
    matrix[16][16] = relax.s_17;
    matrix[17][17] = relax.s_17;
    matrix[18][18] = relax.s_17;
}

Matrix3D::Matrix3D(const Matrix3D &other):matrix(boost::extents[19][19])
{
    matrix = other.getData();
}

//=========================== OPERATORS ===========================

const array3D Matrix3D::operator*(const array3D &other)const {
    array3D a;

    for (int i= 0; i<19;i++){
        double sum = 0;
        for(int j = 0; j<19;j++){
            sum += other[j] * matrix[i][j];
        }
        a[i] = sum;
    }
    return a;
}

const DistributionSetType3D Matrix3D::operator*(const DistributionSetType3D &other)const {
    DistributionSetType3D a;

    a[0] = *this * other[0];
    a[1] = *this * other[1];
 
    return a;
}

const Matrix3D Matrix3D::operator*(double other)const{
    boost::multi_array<double,2> m(boost::extents[19][19]);

    for(int i = 0;i<19;i++){
        for(int j=0;j<19;j++){
            m[i][j] = matrix[i][j] * other;
        }
    }
    return Matrix3D(m);
}

const Matrix3D Matrix3D::operator+(const Matrix3D &other)const{
    boost::multi_array<double,2> m(boost::extents[19][19]); 
    boost::multi_array<double,2> mother = other.getData();

    for(int i = 0;i<19;i++){
        for(int j=0;j<19;j++){
            m[i][j] = matrix[i][j] + mother[i][j];
        }
    }
    return Matrix3D(m);
}

const Matrix3D Matrix3D::operator-(const Matrix3D &other)const{
    boost::multi_array<double,2> m(boost::extents[19][19]); 
    boost::multi_array<double,2> mother = other.getData();

    for(int i = 0;i<19;i++){
        for(int j=0;j<19;j++){
            m[i][j] = matrix[i][j] - mother[i][j];
        }
    }
    return Matrix3D(m);
}

const bool Matrix3D::operator==(const Matrix3D &other)const{
    boost::multi_array<double,2> mother = other.getData();
    bool equal = false;
    if (matrix == mother) equal = true;

    return equal;
}

//=========================== OPERATIONS ===========================

const double Matrix3D::linewise(const array3D &other, int line)const{
    double sum = 0;
    for(int j = 0; j<19;j++){
        sum += other[j] * matrix[line][j];
    }
    return sum;
}

//=========================== ACCESSORS ===========================

void Matrix3D::resetOmega(double omega){
    matrix[9][9]   = omega;
    matrix[11][11] = omega;
    matrix[13][13] = omega;
    matrix[14][14] = omega;
    matrix[15][15] = omega;
}
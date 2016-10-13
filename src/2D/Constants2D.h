#ifndef CONSTANTS2D_H
#define CONSTANTS2D_H

#include"Matrix2D.h"
//#include"Matrix2D_alter.h"
#include"Definitions2D.h"

/// contains constants

/// B from Reis et al.
const array2D B_2D = {{-0.1481481481481481, 0.07407407407407407, 0.04629629629629629, 0.07407407407407407, 0.04629629629629629, 0.07407407407407407, 0.04629629629629629, 0.07407407407407407, 0.04629629629629629}};

/// D2Q9 weights
const array2D WEIGHTS_2D = {{0.4444444444444444, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778}};

/// D2Q13 weights for 4th order isotropic gradient
const boost::array<double,13> GRAD_WEIGHTS_2D = {{0, 0.2666666666666667, 0.1, 0.2666666666666667, 0.1, 0.2666666666666667, 0.1, 0.2666666666666667, 0.1, 0.00833333333333333, 0.00833333333333333, 0.00833333333333333, 0.00833333333333333}};

/// D2Q13 directions
const Vector2D e0(0,0),e1(1,0),e2(1,1),e3(0,1),e4(-1,1),e5(-1,0),e6(-1,-1),e7(0,-1),e8(1,-1),e9(2,0),e10(0,2),e11(-2,0),e12(0,-2);
const boost::array<Vector2D,13> DIRECTION_2D = {{e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12}};
const boost::array<double,13> DIRECTION_ABS_2D = {{e0.Abs(),e1.Abs(),e2.Abs(),e3.Abs(),e4.Abs(),e5.Abs(),e6.Abs(),e7.Abs(),e8.Abs(),e9.Abs(),e10.Abs(),e11.Abs(),e12.Abs()}};

/// D2Q9 streaming indices
const boost::array<int,9> PULL_INDEX_2D = {{0,5,6,7,8,1,2,3,4}};

/// Transformation-Matrix2D
const boost::multi_array<double,2> define_trafo_matrix_2D();
const boost::multi_array<double,2> define_inverse_trafo_matrix_2D();

const Matrix2D TRAFO_MATRIX2D(define_trafo_matrix_2D());
const Matrix2D INV_TRAFO_MATRIX2D(define_inverse_trafo_matrix_2D());

/// Transformation-Matrix2D_alter
// const Eigen::Matrix<double,9,9> define_trafo_matrix_2D_alter();
// const Eigen::Matrix<double,9,9> define_inverse_trafo_matrix_2D_alter();

// const Matrix2D_alter TRAFO_MATRIX2D_alter_alter(define_trafo_matrix_2D_alter());
// const Matrix2D_alter INV_TRAFO_MATRIX2D_alter(define_inverse_trafo_matrix_2D_alter());

#endif
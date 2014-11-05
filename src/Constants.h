#ifndef CONSTANTS_H
#define CONSTANTS_H

#include"Matrix.h"
#include"2D/Definitions2D.h"

/// contains constants

/// Pi
const double PI = 3.14159265358979323846264338327950288419716939937510;

/// B from Reis et al.
const array2D B = {{-0.1481481481481481, 0.07407407407407407, 0.04629629629629629, 0.07407407407407407, 0.04629629629629629, 0.07407407407407407, 0.04629629629629629, 0.07407407407407407, 0.04629629629629629}};

/// D2Q9 weights
const array2D WEIGHTS = {{0.4444444444444444, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778, 0.1111111111111111, 0.02777777777777778}};

/// D2Q13 weights for 4th order isotropic gradient
const boost::array<double,13> GRAD_WEIGHTS = {{0, 0.2666666666666667, 0.1, 0.2666666666666667, 0.1, 0.2666666666666667, 0.1, 0.2666666666666667, 0.1, 0.00833333333333333, 0.00833333333333333, 0.00833333333333333, 0.00833333333333333}};

/// D2Q13 directions
const Vector2D e0(0,0),e1(1,0),e2(1,1),e3(0,1),e4(-1,1),e5(-1,0),e6(-1,-1),e7(0,-1),e8(1,-1),e9(2,0),e10(0,2),e11(-2,0),e12(0,-2);
const boost::array<Vector2D,13> DIRECTION = {{e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12}};

/// Transformation-Matrix2D
const boost::multi_array<double,2> define_trafo_matrix();
const boost::multi_array<double,2> define_inverse_trafo_matrix();

const Matrix2D TRAFO_MATRIX(define_trafo_matrix());
const Matrix2D INV_TRAFO_MATRIX(define_inverse_trafo_matrix());

/// arbitrary  definition
const double MACH_MAX = 0.1; // maximal erlaubte Mach-Zahl

#endif

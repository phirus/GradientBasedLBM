/** TheVector class handles typical 3D vectors
* provides functionfor absolute value and angles between two vectors 
*/

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include<cmath>

class Vector3D
{
    public:
    double x,y,z;

    /// Lifecycle
    Vector(double xn = 0, double yn = 0, double xn = 0):x(xn),y(yn),z(zn){};

    /// operators
    inline const Vector operator+(const Vector& other)const{return Vector(x+other.x, y+other.y, z+other.z);};  /// < addition-operator
    inline const Vector operator-(const Vector& other)const{return Vector(x-other.x, y-other.y, z-other.z);};  /// < subtraction-operator
    inline const double operator*(const Vector&other)const{return x*other.x + y*other.y + z*other.z;};   /// < scalar product
    inline const Vector operator*(double c)const{return Vector(x*c,y*c,z*c);};                   /// < multiplication with a number

    /// operations
    inline const double Abs()const{return sqrt(x*x+y*y+z*z);};   /// < absolute value
    inline const double Sum()const{return x+y+z;};             /// < sum over all components
    const double Angle(const Vector& other)const;
};

#endif

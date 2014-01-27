/** TheVector class handles typical 2D vectors
* provides functionfor absolute value and angles between two vectors 
*/

#ifndef VECTOR_H
#define VECTOR_H

#include<cmath>

class Vector
{
    public:
    double x,y;
    Vector(double xn = 0, double yn = 0):x(xn),y(yn){};

    /// overloaded operators
    inline const Vector operator+(const Vector& other)const{return Vector(x+other.x, y+other.y);};  /// < addition-operator
    inline const Vector operator-(const Vector& other)const{return Vector(x-other.x, y-other.y);};  /// < subtraction-operator

    inline const double operator*(const Vector&other)const{return x*other.x + y*other.y;};   /// < scalar product
    inline const Vector operator*(double c)const{return Vector(x*c,y*c);};                   /// < multiplication with a number

    /// important functions
    inline const double Abs()const{return sqrt(x*x+y*y);};   /// < absolute value
    inline const double Sum()const{return x+y;};             /// < sum over all components
    const double Angle(const Vector& other)const;
};

#endif

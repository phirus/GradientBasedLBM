#include"vector.h"

Vector::Vector(double xn, double yn):x(xn),y(yn){}

const double Vector::angle(const Vector& other)const
{
    double av1 = abs();
    double av2 = other.abs();
    if(av1 == 0 || av2 == 0) return 0; // prevent division by 0
    else return (*this * other) / (av1*av2);
}

#include"Vector.h"

const double Vector::Angle(const Vector& other)const
{
    double av1 = Abs();
    double av2 = other.Abs();
    if(0 == av1 || 0 == av2) return 0; // prevent division by 0
    else return (*this * other) / (av1*av2);
}

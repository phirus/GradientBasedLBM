#ifndef BUBBLEBOX3D_H
#define BUBBLEBOX3D_H

class BubbleBox3D
{
public:
    /// Lifecycle
    BubbleBox3D(double xi=0, double yi=0, double zi=0, double dxi=0, double dyi=0, double dzi=0):x(xi),y(yi),z(zi),dx(dxi),dy(dyi),dz(dzi) {};

    /// operations
    const bool isInBox(double xc, double yc, double zc)const;

    /// accessors
    inline const double getX()const{return x;};
    inline const double getY()const{return y;};
    inline const double getZ()const{return z;};
    inline const double getDX()const{return dx;};
    inline const double getDY()const{return dy;};
    inline const double getDZ()const{return dz;};

    void setX(double xi){x = xi;};
    void setY(double yi){y = yi;};
    void setZ(double zi){z = zi;};

    void setDX(double dxi){dx = dxi;};
    void setDY(double dyi){dy = dyi;};
    void setDZ(double dzi){dz = dzi;};

    void setBubble(double xc, double yc, double zc);

    /// operators
    const bool operator==(const BubbleBox3D& other)const;
    BubbleBox3D& operator=(const BubbleBox3D& other);

private:
    double x,y,z;    // coordinates of left, lower corner
    double dx,dy,dz;
};

#endif // BUBBLEBOX3D_H

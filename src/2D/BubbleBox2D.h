#ifndef BUBBLEBOX2D_H
#define BUBBLEBOX2D_H

class BubbleBox2D
{
public:
    /// Lifecycle
    BubbleBox2D(double xi=0, double yi=0, double wi=0, double hi=0):x(xi),y(yi),w(wi),h(hi){};

    /// operations
    const bool isInBox(double xc, double yc)const;

    /// accessors
    inline const double getX()const{return x;};
    inline const double getY()const{return y;};
    inline const double getW()const{return w;};
    inline const double getH()const{return h;};

    void setX(double xi){x = xi;};
    void setY(double yi){y = yi;};
    void setW(double wi){w = wi;};
    void setH(double hi){h = hi;};
    void setBubble(double xc, double yc);

    /// operators
    const bool operator==(const BubbleBox2D& other)const;
    BubbleBox2D& operator=(const BubbleBox2D& other);

private:
    double x,y;    // coordinates of left, lower corner
    double w,h;
};

#endif // BUBBLEBOX2D_H

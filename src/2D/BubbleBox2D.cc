#include"BubbleBox2D.h"

const bool BubbleBox2D::isInBox(double xc, double yc)const
{
    if (yc < y) return false;
    if (yc > y+h) return false;
    if (xc < x) return false;
    if (xc > x+w) return false;
    return true;
}

void BubbleBox2D::setBubble(double xc, double yc)
{
    x = xc - 0.5 * w;
    y = yc - 0.5 * h;
}

const bool BubbleBox2D::operator==(const BubbleBox2D& other)const
{
    return ((x == other.getX()) && (y == other.getY()) && (w == other.getW()) && (h == other.getH()));
}

BubbleBox2D& BubbleBox2D::operator=(const BubbleBox2D& other)
{
    this->setX(other.getX());
    this->setY(other.getY());
    this->setW(other.getW());
    this->setH(other.getH());

    return *this;
}
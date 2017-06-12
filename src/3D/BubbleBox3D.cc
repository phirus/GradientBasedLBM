#include"BubbleBox3D.h"

const bool BubbleBox3D::isInBox(double xc, double yc, double zc)const
{
    if (zc < z) return false;
    if (zc > z + dz) return false;
    if (yc < y) return false;
    if (yc > y + dy) return false;
    if (xc < x) return false;
    if (xc > x + dx) return false;
    return true;
}

void BubbleBox3D::setBubble(double xc, double yc, double zc)
{
    x = xc - 0.5 * dx;
    y = yc - 0.5 * dy;
    z = zc - 0.5 * dz;
}

const bool BubbleBox3D::operator==(const BubbleBox3D& other)const
{
    return ((x == other.getX()) && (y == other.getY()) && (z == other.getZ()) && (dx == other.getDX()) && (dy == other.getDY()) && (dz == other.getDZ()) );
}

BubbleBox3D& BubbleBox3D::operator=(const BubbleBox3D& other)
{
    this->setX(other.getX());
    this->setY(other.getY());
    this->setZ(other.getZ());
    this->setDX(other.getDX());
    this->setDY(other.getDY());
    this->setDZ(other.getDZ());

    return *this;
}
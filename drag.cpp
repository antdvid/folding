#include "drag.h"
#include <FronTier.h>

// Point Drag
bool PointDrag::isPresetPoint(double p[]) {
    return rad*rad > (p[0]-cent[0])*(p[0]-cent[0]) +
                     (p[1]-cent[1])*(p[1]-cent[1]) +
                     (p[2]-cent[2])*(p[2]-cent[2]);
}

void PointDrag::setVel(double p[], double v[]) {
    if (isPresetPoint(p))
        std::copy(m_v,m_v+3,v);
}

void PointDrag::setAccel(double p[],double a[]) {
    if (isPresetPoint(p))
        std::copy(m_a,m_a+3,a);
}    

PointDrag::PointDrag(const double c[], double r, const double v[], const double a[],  double t)
            : Drag(t), rad(r) {
    if (v) 
	std::copy(v,v+3,m_v);
    else 
	std::fill(m_v,m_v+3,0);

    if (a) 
	std::copy(a,a+3,m_a);
    else 
	std::fill(m_a,m_a+3,0);

    if (c)
        std::copy(c,c+3,cent);
}

//Dravity drag
void GravityDrag::setAccel(double p[], double a[]) {
    if (isPresetPoint(p))
        std::fill(a,a+3,0);
    else
        std::copy(m_a,m_a+3,a);
}

GravityDrag::GravityDrag(const double c[], const double r, const double a[], double t)
             : PointDrag(c,r,nullptr,a,t) {}

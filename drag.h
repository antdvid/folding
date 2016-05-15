#ifndef DRAG_H
#define DRAG_H
#include <FronTier.h>
#include <algorithm>
#include "folding_state.h"

struct Drag {
    double m_t;
    virtual bool isPresetPoint(double[]) = 0; 

    virtual void setVel(double[],double[]) = 0;

    virtual void setAccel(double[],double[]) = 0;

    Drag(double t) : m_t(t) {}
};

struct PointDrag : public Drag {
    double cent[3];
    double rad;
    double m_v[3];
    double m_a[3];
    virtual bool isPresetPoint(double p[]);
    virtual void setVel(double p[], double v[]);

    virtual void setAccel(double p[],double a[]);

    PointDrag(const double c[], double r, const double v[], const double a[],  double t);
};

struct GravityDrag : public PointDrag {
    void setAccel(double p[], double a[]); 

    GravityDrag(const double c[], const double r, const double a[], double t);
};

#endif

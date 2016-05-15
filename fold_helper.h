#ifndef FOLD_HELPER_H
#define FOLD_HELPER_H
#include <FronTier.h>
#include <algorithm>
#include "folding_state.h"

struct Drag {
    double m_t;
    virtual bool isDragPoint(double[]) = 0;
    virtual bool isPresetPoint(double[]) {return false;}
    virtual void setVel(double[],double[]) = 0;
    virtual void setAccel(double[],double[]) = 0;
    Drag(double t) : m_t(t) {}
};

struct PointDrag : public Drag {
    double cent[3];
    double rad;
    double m_v[3];
    double m_a[3];
    virtual bool isDragPoint(double p[]) {
	return false;
    }
    virtual bool isPresetPoint(double p[]) {
	return rad*rad > (p[0]-cent[0])*(p[0]-cent[0]) +
                         (p[1]-cent[1])*(p[1]-cent[1]) +
                         (p[2]-cent[2])*(p[2]-cent[2]);
    }
    virtual void setVel(double p[], double v[]) {
	if (isPresetPoint(p))
            std::copy(m_v,m_v+3,v);
    }

    virtual void setAccel(double p[],double a[]) {
	if (isPresetPoint(p))
	    std::copy(m_a,m_a+3,a);
    }

    PointDrag(const double c[], double r, const double v[], const double a[],  double t)
                : Drag(t), rad(r) {
	if (v) std::copy(v,v+3,m_v);
	else std::fill(m_v,m_v+3,0);

	if (a) std::copy(a,a+3,m_a);
	else std::fill(m_a,m_a+3,0);

	if (c)
        std::copy(c,c+3,cent);
    }
};

struct GravityDrag : public Drag {
    double cent[3];
    double rad;
    double m_a[3];

    virtual bool isDragPoint(double p[]) {
	return true;
    }
    virtual bool isPresetPoint(double p[]) {
	return rad*rad > (p[0]-cent[0])*(p[0]-cent[0]) +
                         (p[1]-cent[1])*(p[1]-cent[1]) +
                         (p[2]-cent[2])*(p[2]-cent[2]);
    }

    virtual void setVel(double p[],double v[]) {
	return;
    }

    virtual void getAccel(double p[], double a[]) {
        if (isDragPoint(p))
	    std::copy(m_a,m_a+3,a);
    }

    GravityDrag(const double c[], const double r, const double a[], double t)
                : Drag(t), rad(r) {

        if (a) std::copy(a,a+3,m_a);
        else std::fill(m_a,m_a+3,0);

	if (c)
        std::copy(c,c+3,this->cent);
    }

};

//adapter for FronTier geometry:
//Tri, Bond, Point
class FoldPoint;

class Cell{
public:
    Cell(){}
    virtual int num_of_points() = 0; 
    virtual FoldPoint* Point(int) = 0; 
    virtual ~Cell(){};
};

class FoldTri: public Cell {
public:
    FoldPoint* points[3];
    int num_of_points(){return 3;}
    FoldPoint* Point(int); 
    FoldTri(TRI*);
    ~FoldTri();
    TRI* getTri() {return tri;}
private:
    FoldTri();
    TRI* tri;
};

class FoldBond: public Cell {
public:
    FoldPoint* points[2];
    int num_of_points(){return 2;}
    FoldPoint* Point(int);
    FoldBond(BOND*);
    ~FoldBond();
private:
    FoldBond();
    BOND* bond;
};

class FoldPoint: public Cell{
public:
    int num_of_points(){return 1;}
    FoldPoint* Point(int);
    FoldPoint(POINT* p) : point(p) {}
    double* getCoords();
    POINT* getPoint(){return point;} 
    double* getVel();
    double* getAccel();
    ~FoldPoint(){};
private:
    POINT* point;
    FoldPoint();
};
#endif

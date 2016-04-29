#ifndef FOLD_HELPER_H
#define FOLD_HELPER_H
#include <FronTier.h>
class Slice{
public:
    enum Side{UPWARDS,DOWNWARDS};
    double getDistance(const double*);
    const double* getCenter() {return center;};
    static void setThickness(double h) {s_thick = h;}
    static double getThickness() {return s_thick;}
    void setSide(Side s) {this->m_side = s;}
    Side getSide() {return m_side;}
    void setDirection(int dir) {this->m_dir = dir;}
    int getDirection() {return m_dir;}
    Slice(double[],Side,int);
    static double s_thick;
private:
    Slice(){}
    double center[3];
    Side m_side;
    int m_dir;
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
private:
    FoldTri();
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
};

class FoldPoint: public Cell{
public:
    int num_of_points(){return 1;}
    FoldPoint* Point(int);
    FoldPoint(POINT* p) : point(p) {}
    double* coords();
    POINT* getFronTierPoint(){return point;} 
    ~FoldPoint(){};
private:
    POINT* point;
    FoldPoint();
};
#endif

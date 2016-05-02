#ifndef FOLD_HELPER_H
#define FOLD_HELPER_H
#include <FronTier.h>
class Slice{
public:
    enum Dir{UPWARDS = 0,DOWNWARDS};
    enum Nor{WEST = 0,EAST,SOUTH,NORTH,LOWER,UPPER};
    double getDistance(const double*);
    const double* getCenter() {return center;};
    static void setThickness(double h) {s_thick = h;}
    static double getThickness() {return s_thick;}
    void setNormal(Nor n) {this->m_nor = n;}
    Nor getNormal() {return m_nor;}
    int getNormalDir() {return m_nor/2;}
    int getNormalSide() {return m_nor%2;}
    void setDirection(Dir dir) {this->m_dir = dir;}
    Dir getDirection() {return m_dir;}
    bool isInGap(const double[]);
    Slice(double[],Dir,Nor);
    static double s_thick;
private:
    Slice(){}
    double center[3];
    Nor m_nor;
    Dir m_dir;
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

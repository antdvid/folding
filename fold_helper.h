#ifndef FOLD_HELPER_H
#define FOLD_HELPER_H
#include <FronTier.h>
class Slice{
    double center[3];
    double normal[3];
    enum Side{POSITIVE_SIDE,NEGATIVE_SIDE};
public:
    double getDistance(const double*);
    Side getSide(const double*);
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
    POINT* point;
    int num_of_points(){return 1;}
    FoldPoint* Point(int);
    FoldPoint(POINT* p) : point(p) {}
    ~FoldPoint(){};
private:
    FoldPoint();
};
#endif

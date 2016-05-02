#ifndef FOLDING_H
#define FOLDING_H
#include<vector>
#include "fold_helper.h"
class Folder {
public:
    virtual void doFlatten(int) = 0;
    virtual void doFolding() = 0;
    virtual ~Folder(){}
    void setDirection(int dir){this->m_direction = dir;} 
    void inputFoldingSlices(Slice*);
    void setThickness(double h){this->m_thickness = h;}
protected:
    std::vector<Slice*> slices;
    std::vector<Cell*> cells;
    std::vector<FoldPoint*> pts;
    void computeRotateCenter();
    Folder(){}
    double getThickness(){return m_thickness;}
    int getDirection(){return m_direction;}
    void getBoundBox(double[2][3]);
    void getFoldBoundBox(double[2][3],Slice*,int);
    void snapSliceToGrid(Slice*);
    virtual void movePointsInGap(Slice*, double) = 0;
    virtual void rotatePoints(Slice*,double h) = 0;
private:
    double m_thickness;
    int    m_direction;
};

class Folder3d:public Folder {
public:
    Folder3d(SURFACE*);
    ~Folder3d();
    void doFlatten(int);
    void doFolding();
    Folder3d(){}
    void movePointsInGap(Slice*,double);
    void rotatePoints(Slice*,double);
private:
    void doFolding(Slice*);
};
#endif

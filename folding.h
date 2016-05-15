#ifndef FOLDING_H
#define FOLDING_H
#include<vector>
#include "fold_helper.h"
#include "ft_spring_solver.h"

class Folder {
public:
    virtual void doFlatten(int) = 0;
    virtual void doFolding() = 0;
    virtual ~Folder(){}
    void addDrag(Drag* drag) { drags.push_back(drag);}
    void setThickness(double h){this->m_thickness = h;}
protected:
    std::vector<Cell*> cells;
    std::vector<Drag*> drags;
    std::vector<FoldPoint*> surf_pts;
    std::vector<FoldPoint*> curv_pts;
    Folder(){}
    double getThickness(){return m_thickness;}
    void getBoundBox(double[2][3]);
private:
    double m_thickness;
    int    m_direction;
};

class Folder3d:public Folder {
public:
    Folder3d(INTERFACE*, SURFACE*);
    ~Folder3d();
    void doFlatten(int);
    void doFolding();
    Folder3d(){}
private:
    void doFolding(Drag*,FT_SpringSolver*);
    INTERFACE* m_intfc;
    void recalcNormal();
};
#endif

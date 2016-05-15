#ifndef FOLDING_H
#define FOLDING_H
#include<vector>
#include "drag.h"
#include "ft_spring_solver.h"

class Folder {
public:
    virtual void doFolding() = 0;
    virtual ~Folder(){}
    void addDrag(Drag* drag) { drags.push_back(drag);}
    void setThickness(double h){this->m_thickness = h;}
protected:
    Folder(){}
    double getThickness(){return m_thickness;}
    std::vector<Drag*> drags;
private:
    double m_thickness;
};

class Folder3d:public Folder {
public:
    Folder3d(INTERFACE*, SURFACE*);
    ~Folder3d();
    void doFolding();
    Folder3d(){}
private:
    void doFolding(Drag*,FT_SpringSolver*);
    INTERFACE* m_intfc;
};
#endif

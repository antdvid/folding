#ifndef FOLDING_H
#define FOLDING_H
#include<vector>
#include "drag.h"
#include "ft_spring_solver.h"
#include "../Collision/collid.h"

class Folder {
public:
    virtual void doFolding() = 0;
    virtual ~Folder(){}
    void addDrag(Drag* drag) { drags.push_back(drag);}
    void addDragsFromFile(std::string);
    void setThickness(double h){this->m_thickness = h;}
    void setFrameStepSize(double dt) {max_dt = dt;}
    virtual void setupMovie(std::string,double) = 0;
protected:
    Folder(){}
    double getThickness(){return m_thickness;}
    double getFrameStepSize() {return max_dt;}
    std::vector<Drag*> drags;
private:
    static double max_dt;
    static double m_thickness;
};

struct Movie;
class Folder3d:public Folder {
public:
    Folder3d(INTERFACE*, SURFACE*);
    ~Folder3d();
    void doFolding();
    Folder3d(){}
    void setupMovie(std::string,double);
private:
    void doFolding(Drag*,FT_SpringSolver*,CollisionSolver*);
    void straightenStrings();
    INTERFACE* m_intfc;
    Movie* movie;
};

struct Movie {
    double mv_dt;
    int mv_count;
    std::string mv_dir;
    INTERFACE* mv_intfc;
    bool doMakeMovie;
    void recordMovieFrame();
    bool isMovieTime(double t);
    Movie() : mv_dt(0),mv_count(1),doMakeMovie(false) {}
};

#endif

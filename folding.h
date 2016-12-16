#ifndef FOLDING_H
#define FOLDING_H
#include<vector>
#include "../Collision/collid.h"
#include "drag.h"
#include "spring_solver.h"
#include "bending.h"

class Folder {
public:
    virtual void doFolding() = 0;
    virtual ~Folder(){}
    void addDrag(Drag* drag) { drags.push_back(drag);}
    void addDragsFromFile(std::string);
    void setThickness(double h){this->m_thickness = h;}
    void setFrameStepSize(double dt) {max_dt = dt;}
    virtual void setupMovie(std::string, std::string, double) = 0;
    void setOdeScheme(SpringSolver::ODE_SCHEME scheme) {ode_scheme = scheme;}
    void setSpringParameters(double, double, double);
    SpringSolver::SpringParameter& getSpringParams(){return spring_params;}
    SpringSolver::ODE_SCHEME getOdeScheme() {return ode_scheme;}
protected:
    Folder(){}
    double getThickness(){return m_thickness;}
    double getFrameStepSize() {return max_dt;}
    std::vector<Drag*> drags;
private:
    static double max_dt;
    static double m_thickness;
    SpringSolver::ODE_SCHEME ode_scheme = SpringSolver::EXPLICIT;
    SpringSolver::SpringParameter spring_params;
};

struct Movie;
class Folder3d:public Folder {
public:
    Folder3d(INTERFACE*, SURFACE*);
    ~Folder3d();
    void doFolding();
    Folder3d(){}
    void setupMovie(std::string, std::string, double);
private:
    void doFolding(Drag*,SpringSolver*,CollisionSolver*);
    void straightenStrings();
    INTERFACE* m_intfc;
    Movie* movie;
};

struct Movie {
    double mv_dt;
    int mv_count;
    std::string mv_gv_dir;
    std::string mv_vtk_dir; 
    INTERFACE* mv_intfc;
    bool doMakeMovie;
    void recordMovieFrame();
    bool isMovieTime(double t);
    Movie() : mv_dt(0),mv_count(1),doMakeMovie(false) {}
};

//helper functions
extern void FT_Intfc2SpringMesh(INTERFACE*, std::vector<SpringVertex*>&);
#endif

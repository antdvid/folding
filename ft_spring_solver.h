#ifndef FT_SpringSolver_H
#define FT_SpringSolver_H
#include "fold_helper.h"
#include "spring_solver.h"
#include "folding_state.h"

class FT_SpringVertex: public SpringVertex
{
public:
   void getExternalAccel(double*); 
   FT_SpringVertex(POINT* p):SpringVertex(p){}
};

class FT_SpringSolver: public SpringSolver
{
private:
    INTERFACE* intfc;
    Drag* drag;
    void assemblePointsFromSurf(SURFACE*);
    void assemblePointsFromCurve(CURVE *);
    void assemblePointsFromNode(NODE*);
    void setConnection(POINT*,POINT*,double);
public: 
    void assemblePoints();
    void presetPoints();
    void setDrag(Drag*); 
    FT_SpringSolver(INTERFACE* input): intfc(input){}
};
#endif


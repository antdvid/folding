#ifndef FT_SpringSolver_H
#define FT_SpringSolver_H
#include "spring_solver.h"
#include "../iFluid/ifluid_state.h"
#include <unordered_map>
#include "drag.h"

class FT_SpringVertex: public SpringVertex
{
public:
   FT_SpringVertex(POINT* p);
};

class FT_SpringSolver: public SpringSolver
{
private:
    INTERFACE* intfc;
    Drag* drag;
    std::unordered_map<POINT*,size_t> ht;
    void assemblePointsFromSurf(SURFACE*);
    void assemblePointsFromCurve(CURVE *);
    void assemblePointsFromNode(NODE*);
    void setConnection(POINT*,POINT*,double);
public: 
    void assemblePoints();
    void presetPoints();
    void setPresetVelocity(SpringVertex* sv);
    void setDrag(Drag*); 
    void updateVelocityToState();
    void updateVelocityFromState();
    void resetVelocity();
    FT_SpringSolver(INTERFACE* input): intfc(input){}
};
#endif


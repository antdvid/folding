#include <FronTier.h>
#include "spring_solver.h"

class BendingForce: public SpringSolver::ExtForceInterface
{
    INTERFACE* intfc;
    void calculateBendingForce3d2003(POINT* p1, TRI* t1, TRI* t2);
    void calculateBendingForce3d2006(POINT* p1, TRI* t1, TRI* t2);
    void clear_surf_point_force(SURFACE*);
public: 
    double* getExternalForce(SpringVertex* sv);
    void computeExternalForce();
    BendingForce(INTERFACE* _intfc) : intfc(_intfc) {}
};


#include <FronTier.h>
#include "spring_solver.h"

class BendingForce: public SpringSolver::ExtForceInterface
{
    INTERFACE* intfc;
    // bend stiffness
    double bends; 
    // bend damping
    double bendd; 
    void calculateBendingForce3d2003(POINT* p1, TRI* t1, TRI* t2);
    void calculateBendingForce3d2006(POINT* p1, TRI* t1, TRI* t2);
    void clear_surf_point_force(SURFACE*);
public: 
    double* getExternalForce(SpringVertex* sv);
    void computeExternalForce();
    void getParaFromFile(const char*); 
    double& getBendStiff() { return bends; }
    double& getBendDamp() { return bendd; }
    BendingForce(INTERFACE* _intfc, double s = 0.01, double d = 0.0) : 
		intfc(_intfc) { bends = s; bendd = d; }
};


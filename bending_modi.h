#include <FronTier.h>
#include "spring_solver.h"
#include <unordered_map>
#include <tuple>

class BendingForce: public SpringSolver::ExtForceInterface
{
    void calculateBendingForce3d2003(POINT* p1, TRI* t1, TRI* t2);
    void calculateBendingForce3d2006(POINT* p1, TRI* t1, TRI* t2);
protected: 
    INTERFACE* intfc;
    // bend stiffness
    double bends; 
    // bend damping
    double bendd; 
    void clear_surf_point_force(SURFACE*);
public: 
    double* getExternalForce(SpringVertex* sv);
    virtual void computeExternalForce();
    void getParaFromFile(const char*); 
    double& getBendStiff() { return bends; }
    double& getBendDamp() { return bendd; }
    BendingForce(INTERFACE* _intfc, double s = 0.0, double d = 0.0) : 
		intfc(_intfc) { bends = s; bendd = d; }
};

class particleBending : public BendingForce {
    std::unordered_map<POINT*, std::tuple<double, double, double>>;  
    void buildMap();    
public : 
    particleBending(INTERFACE* intfc, double s = 0.0, double d = 0.0) : 
		BendingForce(intfc, s, d) {
	buildMap(); 
    } 
    virtual void computeExternalForce(); 
}




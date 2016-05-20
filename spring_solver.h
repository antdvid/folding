#ifndef SPRING_SOLVER_H
#define SPRING_SOLVER_H
#include <vector>
#include <FronTier.h>
#include <algorithm>
#include <string>

struct SpringVertex 
{
    POINT* m_p;
    double* x;
    double v[3];
    double f[3];
    double f_ext[3];
    std::vector<size_t> index_nb;
    std::vector<double> length0; 
    bool is_registered;
    
    SpringVertex(POINT* p):m_p(p),x(Coords(p)),
			   is_registered(false){}
    ~SpringVertex(){}
    double* getCoords(){return Coords(m_p);}
    POINT* getPoint(){return m_p;}
    double* getVel() {return v;}
    double* getExternalAccel() {return f_ext;}
    void addNeighbor(size_t,double);
    bool isRegistered(){return is_registered;}
    void setRegistered(){is_registered = true;}
};

class SpringSolver 
{
    Front* m_front;
    double m_dt;
    void checkVertexNeighbors();
    int findPoints(POINT*);
    void rk4SpringSolver(double,int);
    void computeAccel(SpringVertex*);
    struct SpringParameter {
	double k;
	double lambda;
	double m;
	SpringParameter():k(1000),lambda(0.01),m(0.01){}
    } springParameter;
public:
    void setParameters(double,double,double);
    void printAdjacencyList(std::string);
    void printPointList(std::string);
    double getTimeStepSize();
    virtual void doSolve(double);
    SpringSolver(){}
    SpringSolver(Front* front):m_front(front){}
    ~SpringSolver(){}
    //the following functions need to be 
    //implemented in the derived class
    
    //implement a function to create
    //a vector of SpringVertex and its adjacency
    virtual void assemblePoints() = 0;
    //implement a function to find preset points
    //make its is_registered to be true
    virtual void presetPoints() = 0;
protected:
    std::vector<SpringVertex*> pts;
};

#endif

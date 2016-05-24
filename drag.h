#ifndef DRAG_H
#define DRAG_H
#include <algorithm>
#include "spring_solver.h"

class Drag {
    void operator=(const Drag&);
    Drag(const Drag&);
public:
    double m_t;
    virtual bool isPresetPoint(double*) = 0; 

    virtual void setVel(SpringVertex*) = 0;

    virtual void setAccel(SpringVertex*) = 0;

    void setTimeStepSize(double dt) {m_dt = dt;}
    double getTimeStepSize() {return m_dt;}

    //prototype pattern
    //the drag will be created using 
    //dragFactory() by providing appropriate 
    //data, should not call constructor directly
    friend class DragProtoInit;
    class Info {
	std::string _id;
	std::vector<double> _data;
     public:
	Info(std::string ident, std::vector<double> &dat):
		_id(ident), _data(dat) {}
	Info() {}
	std::string &id() {return _id;}
	const std::string& id() const {return _id;}
	std::vector<double>& data() {return _data;}
	const std::vector<double>& data() const {return _data;}
    };
    virtual std::string id() = 0; 
    static std::vector<Drag*> prototypes;
    static Drag* dragFactory(const Info &);
    virtual Drag* clone(const Info &) = 0;
protected:
    double m_dt;
    Drag(double t) : m_t(t) {}
    Drag(){}
};

class PointDrag : public Drag {
public:
    double cent[3];
    double rad;
    double m_v[3];
    double m_a[3];
    virtual bool isPresetPoint(double*);
    virtual void setVel(SpringVertex*);
    virtual void setAccel(SpringVertex*);

    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "PointDrag";}
    PointDrag(const double c[], double r, const double v[], const double a[],  double t);
    PointDrag(){};
};

class GravityDrag : public PointDrag {
public:
    void setAccel(SpringVertex*); 
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "GravityDrag";}
    GravityDrag(const double c[], const double r, const double a[], double t);
    GravityDrag(){}
};

class MultiplePointDrag : public Drag {
public:
    struct DragPoint {
	double cent[3];
	double rad;
	double m_v[3];
	double m_a[3];
	DragPoint(double[],double,double[],double[]);
    };
    std::vector<DragPoint*> dragPoints;
    bool isPresetPoint(double*);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    int getDragPointNumber(double []);
    void addDragPoint(DragPoint*);
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "MultiplePointDrag";}
    MultiplePointDrag(double t);
    MultiplePointDrag(){}
};

class SpinDrag : public Drag {
public:
    double spinOrig[3];
    double spinDir[3];
    double angVel;
    double foldingBox[2][3];
    
    bool isPresetPoint(double*);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    std::string id() {return "SpinDrag";}
    virtual Drag* clone(const Drag::Info&);
    void spinToAxis(double c[], double dir[], double theta, double p[]);
    SpinDrag() {}
};

class TuckDrag : public SpinDrag {
public:
    double radius;
    double band;
    bool isPresetPoint(double*);
    std::string id() {return "TuckDrag";}
    Drag* clone(const Drag::Info&);
    void setVel(SpringVertex*);
    TuckDrag(){}
};

#endif

#ifndef DRAG_H
#define DRAG_H
#include <algorithm>
#include "spring_solver.h"

class Drag {
    void operator=(const Drag&);
    Drag(const Drag&);
    static double m_thickness;
    static double m_tol;
public:
    double m_t;
    virtual bool isPresetPoint(SpringVertex*) = 0; 

    virtual void setVel(SpringVertex*) = 0;

    virtual void setAccel(SpringVertex*) = 0;

    void setTimeStepSize(double dt) {m_dt = dt;}
    double getTimeStepSize() {return m_dt;}
    
    static void setThickness(double t){m_thickness = t;}
    static double getThickness() {return m_thickness;}
    static void setTolerance(double t){m_tol = t;}
    static double getTolerance() {return m_tol;}

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
    enum {FREE_POINT,ROTATE_POINT,STATIC_POINT};
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
    virtual bool isPresetPoint(SpringVertex*);
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
    bool isPresetPoint(SpringVertex*);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    int getDragPointNumber(double []);
    void addDragPoint(DragPoint*);
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "MultiplePointDrag";}
    MultiplePointDrag(double t);
    MultiplePointDrag(){}
};

class FoldDrag : public Drag {
public:
    double spinOrig[3];
    double spinDir[3];
    double angVel;
    double foldingBox[2][3];
    
    bool isPresetPoint(SpringVertex*);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    std::string id() {return "FoldDrag";}
    virtual Drag* clone(const Drag::Info&);
    FoldDrag() {}
};

class TuckDrag : public Drag {
    double tuck_line[2][3];
    double spinOrig[3];
    double spinDir[3];
    double angVel;
    double constrain_box[2][3];
public:
    bool isPresetPoint(SpringVertex*);
    std::string id() {return "TuckDrag";}
    Drag* clone(const Drag::Info&);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    TuckDrag(){}
};

class CloseUmbrellaDrag : public Drag {
    int num_tuck_line;
    double spinOrig[3];
    double spinDir[3];
    double angVel;
public:
    bool isPresetPoint(SpringVertex*);
    std::string id() {return "CloseUmbrellaDrag";}
    Drag* clone(const Drag::Info&);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*){}
    CloseUmbrellaDrag() {}
};

class RelaxDrag : public Drag {
public: 
    bool isPresetPoint(SpringVertex* sv) {return false;}
    std::string id() {return "RelaxDrag";}
    Drag* clone (const Drag::Info&);
    void setVel(SpringVertex* sv) {}
    void setAccel(SpringVertex* sv) {}
};
#endif

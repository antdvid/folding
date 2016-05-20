#ifndef DRAG_H
#define DRAG_H
#include <algorithm>

class Drag {
    void operator=(const Drag&);
    Drag(const Drag&);
public:
    double m_t;
    virtual bool isPresetPoint(double[]) = 0; 

    virtual void setVel(double[],double[]) = 0;

    virtual void setAccel(double[],double[]) = 0;
    Drag(double t) : m_t(t) {}

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
	std::string id() const {return _id;}
	std::vector<double>& data() {return _data;}
	std::vector<double> data() const {return _data;}
	//const std::vector<double>& data() const {return _data;}
    };
    virtual std::string id() = 0; 
    static std::vector<Drag*> prototypes;
    static Drag* dragFactory(const Info &);
    virtual Drag* clone(const Info &) = 0;
protected:
   Drag(){}
};

class PointDrag : public Drag {
public:
    double cent[3];
    double rad;
    double m_v[3];
    double m_a[3];
    virtual bool isPresetPoint(double p[]);
    virtual void setVel(double p[], double v[]);

    virtual void setAccel(double p[],double a[]);
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "PointDrag";}
    PointDrag(const double c[], double r, const double v[], const double a[],  double t);
    PointDrag(){};
};

class GravityDrag : public PointDrag {
public:
    void setAccel(double p[], double a[]); 
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
    bool isPresetPoint(double []);
    void setVel(double [], double v[]);
    void setAccel(double [], double a[]);
    int getDragPointNumber(double []);
    void addDragPoint(DragPoint*);
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "MultiplePointDrag";}
    MultiplePointDrag(double t);
    MultiplePointDrag(){}
};
#endif

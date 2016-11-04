#include "drag.h"
#include <iostream>
const double EPS = 1e-15;
static double pointToLineDistance(double[],double[], double[]);
static void spinToAxis(double c[], double dir[], double theta, double p[]);
static bool isPointInBox(double*,  double[][3]);
static bool isPointOnLine(double[],double[],double[], double);
static bool isPointInBall(double p[], double c[], double r);
double Drag::m_thickness = 0.001;
double Drag::m_tol = 0.001;

static double distance_between(const double p1[], 
		const double p2[], int dim) {
    double dist = 0;
    for (int i = 0; i < dim; ++i) {
	dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return sqrt(dist);
}
//Drag
Drag* Drag::dragFactory(const Drag::Info& info) {
    for (size_t i = 0; i < prototypes.size(); ++i) {
	if (prototypes[i]->id() == info.id()) 
	{
	    std::cout << "Adding drag: " 
		<< prototypes[i]->id() << std::endl;
	    for (size_t j = 0; j < info.data().size(); ++j)
		std::cout << info.data()[j] << " ";
	    std::cout << std::endl;
	    return prototypes[i]->clone(info);
	}
    }
    std::cout << "Warning: unknown drag type: " 
	      << info.id() << std::endl;
    return NULL;
}

static void spinToAxis(double cen[], double dir[], double theta, double p[])
{
    double a,b,c,u,v,w,x,y,z;
    //normalize dir
    double len = Mag3d(dir);
    if (len < EPS)
    {
	std::cout << "Warning: len < EPS" << std::endl;
        printf("dir = %f %f %f\n",dir[0],dir[1],dir[2]);
	return;
    }
    else
    for (int i = 0; i < 3; ++i)
	dir[i] /= len;

    a = cen[0]; b = cen[1]; c = cen[2];
    u = dir[0]; v = dir[1]; w = dir[2];
    x = p[0]; y = p[1]; z = p[2];
    p[0] = (a*(v*v+w*w)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) + 
            x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta);
    p[1] = (b*(u*u+w*w)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) +
            y*cos(theta) + (c*u-a*w+w*x-u*z)*sin(theta);
    p[2] = (c*(u*u+v*v)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) +
            z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta);
}

// Point Drag
Drag* PointDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 11) {
	std::cout << "Point drag should have "
	          << "11 parameters, " 
		  << "but " << info.data().size() 
		  << " parameters are given" 
		  << std::endl; 
	return NULL;
    }
    else {
	const std::vector<double>& v = info.data();
	const double *it = &(v.front());
	return new PointDrag(it,*(it+3),it+4,it+7,*(it+10));
    }
}

bool PointDrag::isPresetPoint(SpringVertex* sv) {
    double *p = sv->getCoords();
    return rad*rad > (p[0]-cent[0])*(p[0]-cent[0]) +
                     (p[1]-cent[1])*(p[1]-cent[1]) +
                     (p[2]-cent[2])*(p[2]-cent[2]);
}

void PointDrag::setVel(SpringVertex* sv) {
    double *v = sv->getVel();
    if (sv->isRegistered())
        std::copy(m_v,m_v+3,v);
}

void PointDrag::setAccel(SpringVertex* sv) {
    double *a = sv->getExternalAccel();
    if (sv->isRegistered())
        std::copy(m_a,m_a+3,a);
}    

PointDrag::PointDrag(const double c[], double r, const double v[], const double a[],  double t)
            : Drag(t), rad(r) {
    if (v) 
	std::copy(v,v+3,m_v);
    else 
	std::fill(m_v,m_v+3,0);

    if (a) 
	std::copy(a,a+3,m_a);
    else 
	std::fill(m_a,m_a+3,0);

    if (c)
        std::copy(c,c+3,cent);
}

//Dravity drag
void GravityDrag::setAccel(SpringVertex* sv) {
    double *a = sv->getExternalAccel();
    if (sv->isRegistered())
        std::fill(a,a+3,0);
    else
        std::copy(m_a,m_a+3,a);
}

GravityDrag::GravityDrag(const double c[], const double r, const double a[], double t)
             : PointDrag(c,r,nullptr,a,t) {}

Drag* GravityDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 8) {
	std::cout << "Gravity drag should have "
	          << "8 parameters, " 
		  << "but " << info.data().size() 
		  << " parameters are given" 
		  << std::endl; 
	return NULL;
    }
    else {
	const std::vector<double>& v = info.data();
	const double *it = &(v.front());
	return new GravityDrag(it,*(it+3),it+4,*(it+7));
    }
}

//MultiplePointDrag
bool MultiplePointDrag::isPresetPoint(SpringVertex* sv) {
    double *p = sv->getCoords();
    if (getDragPointNumber(p) != -1)
	return true;
    else
        return false;
}

int MultiplePointDrag::getDragPointNumber(double p[]) {
    for (size_t i = 0; i < dragPoints.size(); ++i) 
        if (distance_between(p,dragPoints[i]->cent,3) 
		< dragPoints[i]->rad)
	    return i;
    return -1;
}

void MultiplePointDrag::setVel(SpringVertex* sv) {
    double *p = sv->getCoords();
    double *vel = sv->getVel();

    int index = getDragPointNumber(p);
    if ( index  > -1) {
	size_t i = static_cast<size_t>(index);
	std::copy(dragPoints[i]->m_v,dragPoints[i]->m_v+3,vel);
    }
} 

void MultiplePointDrag::setAccel(SpringVertex* sv) {
    double *p = sv->getCoords();
    double *accel = sv->getExternalAccel();
    int index = getDragPointNumber(p);
    if ( index  > -1) {
	size_t i = static_cast<size_t>(index);
        std::copy(dragPoints[i]->m_a,dragPoints[i]->m_a+3,accel);
    }
}

void MultiplePointDrag::addDragPoint(DragPoint* dragPoint) {
    dragPoints.push_back(dragPoint);
}

MultiplePointDrag::MultiplePointDrag(double t): Drag(t) {}

MultiplePointDrag::DragPoint::DragPoint(double c[], double r,
		double v[], double a[]) {
    std::copy(v,v+3,m_v);
    std::copy(a,a+3,m_a);
    std::copy(c,c+3,cent);
    rad = r;
}

Drag* MultiplePointDrag::clone(const Drag::Info& info) {
    if ((info.data().size() - 1)%10 != 0) {
	std::cout << "Multiple point drag should have "
	          << "10x + 1 parameters, " 
		  << "but " << info.data().size() 
		  << " parameters are given" 
		  << std::endl; 
	return NULL;
    }
    else {
	MultiplePointDrag* mpd = new MultiplePointDrag(info.data().back());
	int num_of_point = (info.data().size()-1)/10;
	std::vector<double> v = info.data();
	for (int i = 0; i < num_of_point; ++i) {
	    double *it = &(v.front()) + i*10;
	    mpd->addDragPoint(new DragPoint(it,*(it+3),it+4,it+7));
	}
	return mpd;
    }
}

//FoldDrag
bool FoldDrag::isPresetPoint(SpringVertex* sv) {
    double *p = sv->getCoords();
    double dist = pointToLineDistance(p,spinOrig,spinDir);
    if (dist < getThickness())
    {
	sv->point_type = FREE_POINT;
	return false;
    }
    else if (isPointInBox(p,foldingBox)) {
	sv->point_type = ROTATE_POINT;
	return true;
    }
    else {
	sv->point_type = STATIC_POINT;
        return true;
    }
}

static double pointToLineDistance(double p[],double orig[], double dir[]) {
    double v[3] = {0};
    double proj_v[3] = {0};
    double len = Mag3d(dir);
    for (size_t i = 0; i < 3; ++i) {
	v[i] = p[i] - orig[i];
	dir[i] /= len;
	proj_v[i] = v[i]*dir[i];
	v[i] = v[i] - proj_v[i];
    }
    return Mag3d(v);
}

static bool isPointInBox(double* p, double box[2][3]) {
    for (size_t i = 0; i < 3; ++i) {
	if (p[i] < box[0][i] || p[i] > box[1][i])
	    return false;
    } 
    return true;
}

void FoldDrag::setVel(SpringVertex* sv) {
    double *p0 = sv->getCoords();
    if (sv->point_type == ROTATE_POINT) {
	double p1[3];
	std::copy(p0,p0+3,p1);
	if (m_dt < EPS)
	    std::cout << "Warning: m_dt = " << m_dt
		      << " is too small. "
		      << "Have you called setTimeStepSize() before using it?" 
		      << std::endl;
	double theta = angVel * m_dt;
	spinToAxis(spinOrig,spinDir,theta,p1);
	double *v = sv->getVel();
	for (size_t i = 0; i < 3; ++i) {
	    v[i] = (p1[i] - p0[i])/m_dt;
	}
    }
    else if (sv->point_type == STATIC_POINT) {
	std::fill(sv->getVel(),sv->getVel()+3,0);
    }
    else 
	return;	
}

void FoldDrag::setAccel(SpringVertex* sv) {}

Drag* FoldDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 14) {
        std::cout << "Fold drag should have "
                  << "14 parameters, "
                  << "but " << info.data().size()
                  << " parameters are given"
                  << std::endl;
        return NULL;
    }
    else {
        FoldDrag* td = new FoldDrag();
        const double* it = &(info.data().front());
        std::copy(it,it+3,td->spinOrig);
        std::copy(it+3,it+6,td->spinDir);
        td->angVel = *(it+6);
	std::copy(it+7,it+10,td->foldingBox[0]);
	std::copy(it+10,it+13,td->foldingBox[1]);
        td->m_t = *(it+13);
        return td;
    }
}

//TuckDrag
static bool isPointOnLine(double p[],  double p1[], 
		     double p2[], double tol) {
    double v[3];
    double v1[3], v2[3];
    for (size_t i = 0; i < 3; ++i) {
	v[i] = p2[i] - p1[i];
	v1[i] = p[i]  - p1[i];
	v2[i] = p[i]  - p2[i];
    }
    if (pointToLineDistance(p,p1,v) < tol && 
	Dot3d(v1,v2) < 0)  {
	return true;		
    }     
    else
	return false;
}

static bool isPointInBall(double p[], double c[], double r) {
    double v[3];
    std::transform(p,p+3,c,v,std::minus<double>());
    return (Mag3d(v) < r) ? true : false;
}

bool TuckDrag::isPresetPoint(SpringVertex* sv) {
    double* p = sv->getCoords();
    if (isPointOnLine(p,tuck_line[0],tuck_line[1],getTolerance()))
    {
	sv->point_type = ROTATE_POINT;
	return true;
    }
    else if (isPointInBox(p,constrain_box))
    {
        sv->point_type = STATIC_POINT;
	return true;
    }
    else
    {
	sv->point_type = FREE_POINT;
	return false;
    }
}

Drag* TuckDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 20) {
	std::cout << "Tuck drag should have "
                  << "20 parameters, "
                  << "but " << info.data().size()
                  << " parameters are given"
                  << std::endl;
        return NULL;
    }
    else {
	TuckDrag* td = new TuckDrag();
	const double* it = &(info.data().front());
	std::copy(it,    it+3, td->tuck_line[0]);
	std::copy(it+3,  it+6, td->tuck_line[1]);
	std::copy(it+6,  it+9, td->spinOrig);
	std::copy(it+9, it+12, td->spinDir);
	td->angVel = *(it+12);
	std::copy(it+13, it+16, td->constrain_box[0]);
	std::copy(it+16, it+19, td->constrain_box[1]);
	td->m_t = *(it+19);
	return td;
    }
}

void TuckDrag::setVel(SpringVertex* sv) {
    double *p0 = sv->getCoords();
    double *v = sv->getVel();
    double p1[3];
    std::copy(p0,p0+3,p1);

    if (sv->point_type == ROTATE_POINT) {
    	double theta = angVel * m_dt;
    	spinToAxis(spinOrig,spinDir,theta,p1);
    	for (size_t i = 0; i < 3; ++i)
            v[i] = (p1[i] - p0[i])/m_dt;
    }
    else 
	return;
}

void TuckDrag::setAccel(SpringVertex* sv) {
    return;
}

//closeUmbrellaDrag
bool CloseUmbrellaDrag::isPresetPoint(SpringVertex* sv) {
    double d_theta = 2*PI/num_tuck_line;
    double theta = 0;
    double *p = sv->getCoords();
    double len = distance_between_positions(p,spinOrig,3);
    if (!isPointInBall(p,spinOrig,EPS)) {
        theta = (p[1] > spinOrig[1]) ? acos((p[0]-spinOrig[0])/len) :
				       acos((p[0]-spinOrig[0])/len) + PI;
	double theta_l = floor(theta/d_theta)*d_theta;
	double theta_r = ceil(theta/d_theta)*d_theta;
	double dist = std::min(len*fabs(sin(theta-theta_l)),len*fabs(sin(theta-theta_r)));
	if (dist < getTolerance()) {
	    sv->point_type = ROTATE_POINT;
	    return true;
	}
	else {
	    sv->point_type = FREE_POINT;
	    return false;
	}
    }
    else {
	sv->point_type = STATIC_POINT;
    	return true;
    }
}

void CloseUmbrellaDrag::setVel(SpringVertex* sv) {
    double *p0 = sv->getCoords();
    double *v = sv->getVel();
    double p1[3];
    std::copy(p0,p0+3,p1);
    double rotateDir[3];
    double p_to_o[3];
    std::transform(p0,p0+3,spinOrig,p_to_o,std::minus<double>());

    if (sv->point_type == ROTATE_POINT) {
        double theta = angVel * m_dt;
	Cross3d(p_to_o,spinDir,rotateDir);
	if (Mag3d(rotateDir) < EPS) return;
        spinToAxis(spinOrig,rotateDir,theta,p1);
        for (size_t i = 0; i < 3; ++i)
            v[i] = (p1[i] - p0[i])/m_dt;
    }
    else
        return;
}

Drag* CloseUmbrellaDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 9) {
        std::cout << "closeUmbrellaDrag should have "
                  << "9 parameters, "
                  << "but " << info.data().size()
                  << " parameters are given"
                  << std::endl;
        return NULL;
    }
    else {
        CloseUmbrellaDrag* cud = new CloseUmbrellaDrag();
        const double* it = &(info.data().front());
	cud->num_tuck_line = *it;
	std::copy(it+1,it+4,cud->spinOrig);	
	std::copy(it+4,it+7,cud->spinDir);	
	cud->angVel = *(it+7);
        cud->m_t = *(it+8);
        return cud;
    }
}

//RelaxDrag
Drag* RelaxDrag::clone(const Drag::Info& info) {
    RelaxDrag* rd = new RelaxDrag();
    const double* it = &(info.data().front());
    rd->m_t = *it;
    return rd;
}

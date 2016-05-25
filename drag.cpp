#include "drag.h"
#include <iostream>
const double EPS = 1e-10;

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

bool PointDrag::isPresetPoint(double* p) {
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
bool MultiplePointDrag::isPresetPoint(double *p) {
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

//spinDrag
void SpinDrag::spinToAxis(double cen[], double dir[], double theta, double p[])
{
    double a,b,c,u,v,w,x,y,z;
    //normalize dir
    double len = 0;
    for (int i = 0; i < 3; ++i)
	len += dir[i]*dir[i];
    if (len < EPS)
	return;
    else
    for (int i = 0; i < 3; ++i)
	dir[i] /= sqrt(len);

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

bool SpinDrag::isPresetPoint(double *p) {
    for (size_t i = 0; i < 3; ++i)
	if (p[i] < foldingBox[0][i] || p[i] > foldingBox[1][i])
	    return false;
    return true;
}

void SpinDrag::setVel(SpringVertex* sv) {
    if (sv->isRegistered()) {
	double *p0 = sv->getCoords();
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
    else 
	return;
}

void SpinDrag::setAccel(SpringVertex* sv) {}

Drag* SpinDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 14) {
        std::cout << "Spin drag should have "
                  << "9 parameters, "
                  << "but " << info.data().size()
                  << " parameters are given"
                  << std::endl;
        return NULL;
    }
    else {
        SpinDrag* td = new SpinDrag();
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

bool TuckDrag::isPresetPoint(double* p) {
    double v1[3], dot3d = 0;
    for (size_t i = 0; i < 3; ++i) {
	v1[i] = p[i] - spinOrig[i];
	dot3d += v1[i] * spinDir[i];
    }    
    if (fabs(dot3d) < EPS && fabs(Mag3d(v1)-radius) < band) {
	return true;
    }
    //including apex and load node
    else if (Mag3d(v1) < band)
	return true;
    else if (sqrt(v1[0]*v1[0]+v1[1]*v1[1]) < band &&
	     p[2] < spinOrig[2])
	return true;
    else
	return false;
}

Drag* TuckDrag::clone(const Drag::Info& info) {
    if (info.data().size() != 10) {
	std::cout << "Tuck drag should have "
                  << "10 parameters, "
                  << "but " << info.data().size()
                  << " parameters are given"
                  << std::endl;
        return NULL;
    }
    else {
	TuckDrag* td = new TuckDrag();
	const double* it = &(info.data().front());
	std::copy(it,it+3,td->spinOrig);
	std::copy(it+3,it+6,td->spinDir);
	td->angVel = *(it+6);
	td->radius = *(it+7);
	td->band = *(it+8);
	td->m_t = *(it+9);
	td->setLoadNodeVel();
	return td;
    }
}

void TuckDrag::setVel(SpringVertex* sv) {
    double *p0 = sv->getCoords();
    double *v = sv->getVel();
    double p1[3];
    std::copy(p0,p0+3,p1);

    double vec[3];
    for (size_t i = 0; i < 3; ++i) {
	vec[i] = p0[i] - spinOrig[i];
    }

    //make the apex fixed
    if (Mag3d(vec) < band) {
	std::fill(v,v+3,0);
	return;
    }
    
    double crx[3];
    Cross3d(vec,spinDir,crx);
    if (m_dt < EPS)
        std::cout << "Warning: m_dt = " << m_dt
                  << " is too small. "
                  << "Have you called setTimeStepSize() before using it?"
                  << std::endl;
    double theta = angVel * m_dt;
    spinToAxis(spinOrig,crx,theta,p1);
    for (size_t i = 0; i < 3; ++i) {
        v[i] = (p1[i] - p0[i])/m_dt;
    }
    //let the load node go down 
    if (sqrt(vec[0]*vec[0]+vec[1]*vec[1]) < band &&
	p0[2] < spinOrig[2]) {
	v[0] = v[1] = 0.0;
	v[2] = load_v;
	return;
    }
}

void TuckDrag::setLoadNodeVel() {
    double t = angVel*m_t;
    load_v = radius*sin(t)/m_t;
}

#include "drag.h"
#include <iostream>

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
	    std::cout << "Add drag " << info.id() << std::endl;
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
	std::vector<double> v = info.data();
	double *it = &(v.front());
	return new PointDrag(it,*(it+3),it+4,it+7,*(it+10));
    }
}

bool PointDrag::isPresetPoint(double p[]) {
    return rad*rad > (p[0]-cent[0])*(p[0]-cent[0]) +
                     (p[1]-cent[1])*(p[1]-cent[1]) +
                     (p[2]-cent[2])*(p[2]-cent[2]);
}

void PointDrag::setVel(double p[], double v[]) {
    if (isPresetPoint(p))
        std::copy(m_v,m_v+3,v);
}

void PointDrag::setAccel(double p[],double a[]) {
    if (isPresetPoint(p))
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
void GravityDrag::setAccel(double p[], double a[]) {
    if (isPresetPoint(p))
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
	std::vector<double> v = info.data();
	double *it = &(v.front());
	return new GravityDrag(it,*(it+3),it+4,*(it+7));
    }
}

//MultiplePointDrag
bool MultiplePointDrag::isPresetPoint(double p[]) {
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

void MultiplePointDrag::setVel(double p[], double vel[]) {
    int index = getDragPointNumber(p);
    if ( index  > -1) {
	size_t i = static_cast<size_t>(index);
	std::copy(dragPoints[i]->m_v,dragPoints[i]->m_v+3,vel);
    }
} 

void MultiplePointDrag::setAccel(double p[], double accel[]) {
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

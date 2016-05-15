#include "fold_helper.h"
#include <FronTier.h>

FoldTri::FoldTri(TRI* t) {
	for (int i = 0; i < 3; ++i)
	    points[i] = new FoldPoint(Point_of_tri(t)[i]);
	this->tri = t;
}

FoldTri::~FoldTri() {
	for (int i = 0; i < 3; ++i)
	    delete points[i];
}

FoldPoint* FoldTri::Point(int i) {return points[i];}

FoldBond::FoldBond(BOND* b){
	points[0] = new FoldPoint(b->start);
	points[1] = new FoldPoint(b->end);
	this->bond = b;
}

FoldBond::~FoldBond() {
	for (int i = 0; i < 2; ++i)
	    delete points[i];
}

FoldPoint* FoldBond::Point(int i) {return points[i];}

FoldPoint* FoldPoint::Point(int i) {return this;}

double* FoldPoint::getCoords() { return Coords(point);}

#include "fold_helper.h"
#include <FronTier.h>
double Slice::getDistance(const double* p) {
	return 0;
}

Slice::Side Slice::getSide(const double* p) {
	return POSITIVE_SIDE;
}

FoldTri::FoldTri(TRI* t) {
	for (int i = 0; i < 3; ++i)
	    points[i] = new FoldPoint(Point_of_tri(t)[i]);
}

FoldTri::~FoldTri() {
	for (int i = 0; i < 3; ++i)
	    delete points[i];
}

FoldPoint* FoldTri::Point(int i) {return points[i];}

FoldBond::FoldBond(BOND* b){
	points[0] = new FoldPoint(b->start);
	points[1] = new FoldPoint(b->end);
}

FoldBond::~FoldBond() {
	for (int i = 0; i < 2; ++i)
	    delete points[i];
}

FoldPoint* FoldBond::Point(int i) {return points[i];}

FoldPoint* FoldPoint::Point(int i) {return this;}

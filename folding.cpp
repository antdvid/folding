#include "folding.h"

static void unsortSurface(SURFACE*);

void Folder::computeRotateCenter(){

}

void Folder::inputFoldingSlices(Slice* s) {
    slices.push_back(s);
}

void Folder::getBoundBox(double bb[2][3]) {
    for (int i = 0; i < 3; i++)
    {
	bb[0][i] = HUGE;
	bb[1][i] = -HUGE;
    }
    for (std::vector<Cell*>::iterator it = cells.begin();
	 it != cells.end(); it ++) 
    for (int i = 0; i < (*it)->num_of_points(); ++i) 
    {
	FoldPoint* p = (*it)->Point(i);
	double* crds = p->coords();
	for (int j = 0; j < 3; ++j) {
	    bb[0][j] = std::min(bb[0][j],crds[j]);
	    bb[1][j] = std::max(bb[1][j],crds[j]);
	}
    } 
}

void Folder3d::doFlatten(int dir) {
    double bb[2][3];
    getBoundBox(bb);
    double h = getThickness();
    double flatDist = 0.5*(bb[1][dir] - bb[0][dir] - h);
    double flatLine[2];
    flatLine[0] = bb[0][dir] + flatDist;
    flatLine[1] = bb[1][dir] - flatDist;
    for (std::vector<Cell*>::iterator it = cells.begin();
         it != cells.end(); it ++)
    for (int i = 0; i < (*it)->num_of_points(); ++i) 
    {
        FoldPoint* p = (*it)->Point(i);
        double* crds = p->coords();
        crds[dir] = std::max(crds[dir],flatLine[0]);
        crds[dir] = std::min(crds[dir],flatLine[1]);
    }
}

void Folder3d::doFolding() {
    doFlatten(getDirection());
    printf("slice number = %lu\n",slices.size());
    for (std::vector<Slice*>::iterator it = slices.begin();
	 it != slices.end(); ++it) 
    {
	printf("fold base on %lu slice\n",it-slices.begin());
	doFolding(*it);
    }
}

void Folder3d::doFolding(Slice* s) {
    double bb[2][3];
    getBoundBox(bb);
    double h = getThickness()*0.5;
    int dir = getDirection();
    double rot_cent[3];
    movePointsInGap(s,h);
    if (s->getSide() == Slice::UPWARDS)
        rotatePoints(s,bb[1][dir]+h);
    else
	rotatePoints(s,bb[0][dir]-h);
}

void Folder3d::movePointsInGap(Slice* s, double h) {
    int slice_dir = s->getDirection();
    Slice::Side side = s->getSide();
    int fold_dir = this->getDirection();
    double gap = s->getThickness();
    const double* center = s->getCenter();

    for (std::vector<FoldPoint*>::iterator it = pts.begin();
	 it != pts.end(); ++it) 
    {
	double* crds = (*it)->coords();
	if (crds[slice_dir] > center[slice_dir]-gap && 
	    crds[slice_dir] < center[slice_dir]+gap) {
	    crds[fold_dir] = (side == Slice::UPWARDS) ? 
			     crds[fold_dir] + h :
			     crds[fold_dir] - h;
        }
    }
}

void Folder3d::rotatePoints(Slice* s,double h) {
    int dir = s->getDirection();
    const double* center = s->getCenter();
    for (std::vector<FoldPoint*>::iterator it = pts.begin();
         it != pts.end(); ++it) 
    {
	
    }
}

Folder3d::Folder3d(SURFACE* s) {
    TRI* t;
    unsortSurface(s);
    surf_tri_loop(s,t) {
	cells.push_back(new FoldTri(t));
	Cell* c = cells.back();
	for (int i = 0; i < 3; ++i) {
	    FoldPoint* p = c->Point(i);
	    POINT* fp = p->getFronTierPoint();
	    if (!sorted(fp)) {
		sorted(fp) = YES;
		pts.push_back(p);
	    }	    
	}
    }
}

static void unsortSurface(SURFACE* s) {
    TRI* t;
    surf_tri_loop(s,t) {
	for (int i = 0; i < 3; ++i) {
	    POINT* p = Point_of_tri(t)[i];
	    sorted(p) = NO;
	}
    }
}

Folder3d::~Folder3d() {
    cells.clear();	
}

#include "folding.h"

static void unsortSurface(SURFACE*);

void Folder::computeRotateCenter(){

}

void Folder::inputFoldingSlices(Slice* s) {
    snapSliceToGrid(s);
    slices.push_back(s);
}

void Folder::snapSliceToGrid(Slice* s) {
   int dir = s->getNormalDir();
   const double* cen = s->getCenter();
   double min_dist = HUGE;
   double best_center[3] = {0.0};
   memcpy(best_center,cen,3*sizeof(double));

   for (std::vector<FoldPoint*>::iterator it = pts.begin();
	it != pts.end(); ++it) 
   {
	double* crds = (*it)->coords();
	double dist = fabs(cen[dir]-crds[dir]);
	if (dist < min_dist) {
	    min_dist = dist;
	    best_center[dir] = crds[dir];
	}
   }
   s->setCenter(best_center);
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

void Folder::getFoldBoundBox(double bb[2][3], Slice* s, int side) {
    //side = -1 count lower part;
    //side = 0  count inner part;
    //side = 1  count higher part;
    double gap = s->getThickness()*0.5;
    int dir = s->getNormalDir();
    const double* center = s->getCenter();
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
	if (side == 0 && 
	   (crds[dir] > center[dir]+gap ||
	    crds[dir] < center[dir]-gap))
	    continue;
	else if (side == 1 &&
	    crds[dir] <= center[dir]+gap)
	    continue;
	else if (side == -1 &&
	    crds[dir] >= center[dir]-gap)
	    continue;
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
    double flatDist = std::max(0.5*(bb[1][dir] - bb[0][dir] - h),0.0);
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
    for (std::vector<Slice*>::iterator it = slices.begin();
	 it != slices.end(); ++it) 
    {
	printf("folding base on %lu slice ...\n",it-slices.begin());
	doFolding(*it);
    }
}

void Folder3d::doFolding(Slice* s) {
    int slice_dir = s->getNormalDir();
    int folder_dir = this->getDirection();
    s->setThickness(this->getThickness());
    if (slice_dir == folder_dir) return;

    double bb[2][3];
    getFoldBoundBox(bb,s,0);
    double fbb[2][3];
    int fside = (s->getNormalSide() == 0) ? 1 : -1;
    getFoldBoundBox(fbb,s,fside);
    double h = 0.0;
    if (s->getDirection() == Slice::UPWARDS)
	h = 0.501*fabs(fbb[1][folder_dir]-bb[1][folder_dir]);
    else
	h = 0.501*fabs(fbb[0][folder_dir]-bb[0][folder_dir]);

    h = std::max(h,this->getThickness());

    movePointsInGap(s,h);

    if (s->getDirection() == Slice::UPWARDS)
        rotatePoints(s,bb[1][folder_dir]+h);
    else
	rotatePoints(s,bb[0][folder_dir]-h);
}

void Folder3d::movePointsInGap(Slice* s, double h) {
    int slice_dir = s->getNormalDir();
    Slice::Dir push_dir = s->getDirection();
    double gap = s->getThickness()/2.0;
    const double *center = s->getCenter();
    int fold_dir = this->getDirection();

    for (std::vector<FoldPoint*>::iterator it = pts.begin();
	 it != pts.end(); ++it) 
    {
	double* crds = (*it)->coords();
	if (crds[slice_dir] > center[slice_dir]-gap && 
	    crds[slice_dir] < center[slice_dir]+gap) {
	    crds[fold_dir] = (push_dir == Slice::UPWARDS) ? 
			     crds[fold_dir] + h :
			     crds[fold_dir] - h;
        }
    }
}

void Folder3d::rotatePoints(Slice* s,double h) {
    
/*
 *before rotation:

        *  rotate center   
   _____/\_____
   _____/\_____  
         
	^      folder direction + push direction
	|
        |
        -----> slice direction + side 

 *after rotation:
    ____________
   /  __________
  /  /
  \  \__________
   \____________

*/
    int slice_dir = s->getNormalDir();
    int side = s->getNormalSide();
    int folder_dir = this->getDirection();
    double gap = s->getThickness()*0.5;
    const double* center = s->getCenter();
    double rot_cen[3];

    for (std::vector<FoldPoint*>::iterator it = pts.begin();
         it != pts.end(); ++it) 
    {
	double* crds = (*it)->coords();
	memcpy(rot_cen,crds,3*sizeof(double));
	rot_cen[folder_dir] = h;
	rot_cen[slice_dir] = center[slice_dir];
	if ((side == 0 && 
	    crds[slice_dir] > rot_cen[slice_dir] + gap) ||
	    (side == 1 &&
	    crds[slice_dir] < rot_cen[slice_dir] - gap))
	{
	    for (int j = 0; j < 3; ++j)
	        crds[j] = 2*rot_cen[j] - crds[j];    
	}
        else if (s->isInGap(crds))
	{
	    crds[folder_dir] = rot_cen[folder_dir];
	    double dist = distance_between_positions(crds,rot_cen,3);
	    if (side == 0)
		crds[slice_dir] = rot_cen[slice_dir] + dist;
	    else
		crds[slice_dir] = rot_cen[slice_dir] - dist;
	}
   	else 
	    continue;
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

#include "folding.h"
#include "spring_solver.h"
#include "ft_spring_solver.h"

static void unsortSurface(SURFACE*);
static void unsortCurve(CURVE*);

//abstract base class Folder
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
	double* crds = p->getCoords();
	for (int j = 0; j < 3; ++j) {
	    bb[0][j] = std::min(bb[0][j],crds[j]);
	    bb[1][j] = std::max(bb[1][j],crds[j]);
	}
    } 
}

//derived class Folder3d
void Folder3d::doFlatten(int dir) {
    double bb[2][3];
    getBoundBox(bb);
    double th = getThickness();
    double centerLine = 0.5*(bb[1][dir] + bb[0][dir]);
    double compRatio = th/(bb[1][dir]-bb[0][dir]);
    for (std::vector<FoldPoint*>::iterator it = surf_pts.begin();
	it != surf_pts.end(); ++it) 
    {
	double *crds = (*it)->getCoords();
	crds[dir] = (crds[dir]-centerLine)*compRatio+centerLine;	
    }
    recalcNormal();
}

void Folder3d::doFolding() {
    FT_SpringSolver* sp_solver = new FT_SpringSolver(m_intfc);
    sp_solver->assemblePoints();
    sp_solver->setParameters(100,0.1,0.01);
    for (std::vector<Drag*>::iterator it = drags.begin();
	 it != drags.end(); ++it) 
    {
	printf("folding base on %lu drag ...\n",it-drags.begin());
	doFolding(*it,sp_solver);
	
    }
    recalcNormal();
}

void Folder3d::recalcNormal() {
    typedef TRI FT_TRI;
    for (std::vector<Cell*>::iterator it = cells.begin();
	it != cells.end(); ++it) {
	FoldTri* foldtri = dynamic_cast<FoldTri*>(*it);
	if (foldtri) {
	    FT_TRI* t = foldtri->getTri();
	    set_normal_of_tri(t);
	}
    }
}

void Folder3d::doFolding(Drag* drag, FT_SpringSolver* sp_solver) {
    sp_solver->setDrag(drag);
    sp_solver->doSolve(drag->m_t);
}

Folder3d::Folder3d(INTERFACE* intfc, SURFACE* s) : m_intfc(intfc)
{
    TRI* t;
    BOND* b;
    CURVE** cc;
    unsortSurface(s);
    intfc_curve_loop(intfc,cc) 
	unsortCurve(*cc);

    surf_tri_loop(s,t) {
	cells.push_back(new FoldTri(t));
	Cell* c = cells.back();
	for (int i = 0; i < c->num_of_points(); ++i) {
	    FoldPoint* p = c->Point(i);
	    POINT* fp = p->getPoint();
	    if (!sorted(fp)) {
		sorted(fp) = YES;
		surf_pts.push_back(p);
	    }	    
	}
    }
    intfc_curve_loop(intfc,cc) {
	if (hsbdry_type(*cc) != STRING_HSBDRY)
	    continue;
	curve_bond_loop(*cc,b) {
	    cells.push_back(new FoldBond(b));
	    Cell* c = cells.back();
	    for (int i = 0; i < c->num_of_points(); ++i) {
	        FoldPoint* p = c->Point(i);
	        POINT* fp = p->getPoint();
	        if (!sorted(fp)) {
		    sorted(fp) = YES;
	    	    curv_pts.push_back(p);
	        }
	    }
	}
    }
    printf("In folding:\n");
    printf("%lu number of points in surf\n%lu number of points in curve\n",
	    surf_pts.size(), curv_pts.size());
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

static void unsortCurve(CURVE* c) {
    BOND* b;
    curve_bond_loop(c,b) {
	sorted(b->start) = NO;
	sorted(b->end) = NO;
    }
}

Folder3d::~Folder3d() {
    cells.clear();	
}

#include "ft_spring_solver.h"

static void unsortIntfc(INTERFACE*);
static void unsortSurface(SURFACE*);
static void unsortCurve(CURVE*);
static bool isElasticCurve(CURVE* c) {
    return (hsbdry_type(c) == STRING_HSBDRY ||
            hsbdry_type(c) == MONO_COMP_HSBDRY ||
            hsbdry_type(c) == GORE_HSBDRY);
}

FT_SpringVertex::FT_SpringVertex(POINT* p) : SpringVertex(p){
    STATE* state = (STATE*)left_state(p);
    std::copy(state->vel,state->vel+3,this->v);
    for (int i = 0; i < 3; ++i)
	f_ext[i] = state->fluid_accel[i] +
		   state->other_accel[i];
}

void FT_SpringSolver::assemblePoints() {
    //assemble tris list from input intfc
    //this function should be called before
    //spring interior dynamics computed
    unsortIntfc(intfc);
    pts.clear();

    SURFACE** s;
    intfc_surface_loop(intfc,s) {
        if (wave_type(*s) == ELASTIC_BOUNDARY)
            assemblePointsFromSurf(*s);
    }

    CURVE** c;
    intfc_curve_loop(intfc,c) {
	if (isElasticCurve(*c))
            assemblePointsFromCurve(*c);
    }

    NODE** n;
    intfc_node_loop(intfc,n) {
	assemblePointsFromNode(*n);
    }
}

void FT_SpringSolver::assemblePointsFromSurf(SURFACE* s) {
	TRI* tris[20];
	TRI* tri;
        surf_tri_loop(s,tri)
        for (int i = 0; i < 3; ++i)
        {
            POINT* p = Point_of_tri(tri)[i];
            if(Boundary_point(p)) continue;
            int nt = 0;
            PointAndFirstRingTris(p,Hyper_surf_element(tri),
                            Hyper_surf(s),&nt,tris);
            for (int j = 0; j < nt; ++j)
            for (int l = 0; l < 3; ++l)
            if (Point_of_tri(tris[j])[l] == p)
            {
                POINT* p_nb = Point_of_tri(tris[j])[(l+1)%3];
		setConnection(p,p_nb,tris[j]->side_length0[(l+1)%3]);
            }
        }
}

void FT_SpringSolver::assemblePointsFromCurve(CURVE *curve) {
 	//assemble from curve
	for (BOND* b = curve->first; b != curve->last; b = b->next)
        for (int i = 0; i < 2; ++i)
        {
            POINT* p = b->end;
            POINT* p_nb = (i == 0) ? b->start : b->next->end;
            double dist = (i == 0) ? b->length0 : b->next->length0;
	    setConnection(p,p_nb,dist);
        }

	//connect with btris 
	TRI** tris;
	for (BOND* b = curve->first; b != curve->last; b = b->next)
	{
	    POINT* p = b->end;
	    for (BOND_TRI** btris = Btris(b); btris && *btris; ++btris)
            {
                int nt = I_FirstRingTrisAroundPoint(p,(*btris)->tri,&tris);
		for (int j = 0; j < nt; ++j)
                {
                    for (int side = 0; side < 3; ++side)
                    {
                        if (p == Point_of_tri(tris[j])[side])
                        {
                            if (is_side_bdry(tris[j],side))
                                continue;
                            POINT* p_nb = Point_of_tri(tris[j])[(side+1)%3];
			    setConnection(p,p_nb,tris[j]->side_length0[(side+1)%3]);
                        }
                    }
                }
	    }
	}
}

void FT_SpringSolver::assemblePointsFromNode(NODE* n) {
    for (CURVE** c = n->in_curves; c && *c; ++c)
    {
	if (!isElasticCurve(*c)) continue;
	BOND* b = (*c)->last;
	POINT* p_nb = b->start;
	setConnection(n->posn,p_nb,b->length0);	
    }

    for (CURVE** c = n->out_curves; c && *c; ++c)
    {
	if (!isElasticCurve(*c)) continue;
	BOND* b = (*c)->first;
	POINT* p_nb = b->end;
	setConnection(n->posn,p_nb,b->length0);	
    }
}

void FT_SpringSolver::setConnection(POINT* p1, POINT* p2, double length0) {
    POINT* points[2] = {p1,p2};
    for (int i = 0; i < 2; ++i) {
	if (sorted(points[i])) continue;
	sorted(points[i]) = YES;
	points[i]->indx = pts.size();
	pts.push_back(new FT_SpringVertex(points[i]));
    }
    if (length0 < 0)
	length0 = separation(p1,p2,3);
    size_t ind = static_cast<size_t>(p1->indx);
    size_t ind_nb = static_cast<size_t>(p2->indx);
    pts[ind]->addNeighbor(ind_nb,length0);
    pts[ind_nb]->addNeighbor(ind,length0);
}

void FT_SpringSolver::presetPoints() {
    printf("Entering presetPoints()\n");
    for (size_t i = 0; i < pts.size(); ++i) 
    {
	drag->setVel(pts[i]->getCoords(),pts[i]->getVel());
	drag->setAccel(pts[i]->getCoords(),pts[i]->getExternalAccel());
	if (drag->isPresetPoint(pts[i]->getCoords())) {
	    printf("Registered point = [%f %f %f]\n",
		pts[i]->getCoords()[0],pts[i]->getCoords()[1],
		pts[i]->getCoords()[2]);
	    pts[i]->setRegistered();
	}
	
    }	
}

void FT_SpringSolver::setDrag(Drag* dg) {
    drag = dg;
    presetPoints();
}

static void unsortIntfc(INTERFACE* intfc) {
    SURFACE** s;
    intfc_surface_loop(intfc,s)
        unsortSurface(*s);
    CURVE** c;
    intfc_curve_loop(intfc,c)
        unsortCurve(*c);
}

static void unsortSurface(SURFACE* s) {
    TRI* t;
    surf_tri_loop(s,t) {
        for (int i = 0; i < 3; ++i) {
            POINT* p = Point_of_tri(t)[i];
            sorted(p) = NO;
            p->indx = -1;
        }
    }
}

static void unsortCurve(CURVE* c) {
    BOND* b;
    curve_bond_loop(c,b) {
        sorted(b->start) = NO;
        sorted(b->end) = NO;
        b->start->indx = b->end->indx = -1;
    }
}


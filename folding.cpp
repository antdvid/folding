#include "folding.h"
#include "spring_solver.h"
#include <fstream>
#include <sstream>

static void setCollisionFreePoints3d(INTERFACE*, Drag*);

double Folder::m_thickness = 0.001;
double Folder::max_dt = 0.01;

void Folder::addDragsFromFile(std::string fname) {
    std::ifstream ifs(fname);
    if (!ifs.is_open()) {
	std::cout << "File \"" << fname << 
	"\" doesn't exist" << std::endl;
    }

    std::string line;
    std::string temp; 
    double tmp;
    Drag::Info info;
    std::set<std::string> foldset;
    
    for (std::vector<Drag*>::iterator it = Drag::prototypes.begin();
                it != Drag::prototypes.end(); it++)
         foldset.insert((*it) -> id());
    getline(ifs, line);
    while (true) {
        if (ifs.eof())
        {
            if (!info.empty())
            {
                addDrag(Drag::dragFactory(info));
                info.clear();
            }
            break;
        }
        if (line[0] == '\0')
        {
            addDrag(Drag::dragFactory(info));
            info.clear();
            getline(ifs, line);
            continue;
        }

        std::istringstream ss(line);
        ss >> temp;
        if (foldset.find(temp) != foldset.end())
        {
                info.id() = temp;
                while (ss >> tmp)
                        info.data().push_back(tmp);
                addDrag(Drag::dragFactory(info));
                info.clear();
        }
        else if (temp == "Foldingplan")
        {
            ss >> info.id();
            while (ss >> tmp)
                info.data().push_back(tmp);
        }
	else
        {
            while (temp.back() != ':')
                ss >> temp;
            while (ss >> tmp)
                info.data().push_back(tmp);
        }
        getline(ifs, line);
    }
    ifs.clear();
    ifs.close();
}

void Folder::setSpringParameters(double k, double lambda, double m) {
    spring_params.k = k;
    spring_params.lambda = lambda;
    spring_params.m = m;
}

void Folder3d::doFolding() {
    SpringSolver* sp_solver = SpringSolver::createSpringSolver(
						getOdeScheme()); 
    CollisionSolver* cd_solver = new CollisionSolver3d();
    
    //configure collision solver
    cd_solver->assembleFromInterface(m_intfc,getFrameStepSize());
    cd_solver->setSpringConstant(getSpringParams().k);
    cd_solver->setFrictionConstant(getSpringParams().lambda);
    cd_solver->setPointMass(getSpringParams().m);
    cd_solver->setFabricThickness(getThickness());
  
    //configure spring solver
    FT_Intfc2SpringMesh(m_intfc, sp_solver->getSpringMesh());

    if (getSpringParams().k == 0 || getSpringParams().m == 0)
	throw std::invalid_argument("tensile stiffness and mass cannot zero");
    sp_solver->setParameters(getSpringParams());
    //sp_solver->ext_forces.push_back(new BendingForce(m_intfc));

    Drag::setTolerance(m_intfc->table->rect_grid.h[0]*0.5);
    Drag::setThickness(0.001);

    for (std::vector<Drag*>::iterator it = drags.begin();
	 it != drags.end(); ++it) 
    {
	if (*it == nullptr) continue;
	printf("folding base on %lu drag: %s ...\n",
		it-drags.begin(),(*it)->id().c_str());
	doFolding(*it,sp_solver,cd_solver);
	sp_solver->resetVelocity();
    }
    straightenStrings();
}

double Folder3d::computePotentialEnergy()
{
    SURFACE** s;
    double E = 0;
    intfc_surface_loop(m_intfc, s)
    {
	if (wave_type(*s) != ELASTIC_BOUNDARY)
	    continue;
	TRI* t;
	surf_tri_loop(*s, t)
	{
	    for (int i = 0; i < 3; ++i)
	    {
		double len = separation(Point_of_tri(t)[i],
				Point_of_tri(t)[(i+1)%3], 3); 
		E += 0.5*spring_params.k*pow(len - t->side_length0[i], 2.0);	
	    }
	}
    } 
    return E*0.5;
}

void Folder3d::doFolding(
     Drag* drag, 
     SpringSolver* sp_solver,
     CollisionSolver* cd_solver) 
{
    double min_dt = sp_solver->getTimeStepSize();
    double max_dt = std::max(min_dt,getFrameStepSize());
    
    sp_solver->setTimeStepSize(max_dt);
    sp_solver->setDrag(drag);

    setCollisionFreePoints3d(m_intfc,drag);

    //call spring solver and collision solver
    //based on adaptive dt
    static double t = 0;
    double dt = max_dt;
    while (t < drag->m_t - MACH_EPS) {
	printf("--------------------------------\n");
	printf("dt = %f, t = %f, total time = %f\n",
		dt,t,drag->m_t);
	printf("--------------------------------\n");
	if (t+dt > drag->m_t) 
	    dt = drag->m_t-t;

	cd_solver->recordOriginPosition();
	cd_solver->setTimeStepSize(dt);

    	sp_solver->solve(dt);
	
	recordData(t,movie->out_dir);

	cd_solver->resolveCollision();

	t = t + dt;
	if (movie->isMovieTime(t))
	    movie->recordMovieFrame();
    }
}

void Folder3d::appendDataToFile(double x, double y, std::string fname)
{
    //if new file, clear it and add to list
    std::ofstream ofs;
    if (dataFileSet.find(fname) == dataFileSet.end())
    {
	dataFileSet.insert(fname);
	ofs.open(fname, std::ofstream::trunc);
	ofs.close();
    }
    ofs.open(fname, std::ofstream::app);
    ofs << x << " " << y << std::endl;
    ofs.close();
}

void Folder3d::recordData(double t, std::string out_dir)
{
    appendDataToFile(t, computePotentialEnergy(), out_dir+"/"+"energy");
}

void Folder3d::setupMovie(std::string dname, std::string oname, 
			double dt) {
    movie = new Movie();
    movie->mv_dt = dt;
    std::string pathgv = oname + '/' + dname + "/gview";
    std::string pathvtk = oname + '/' + dname + "/vtk"; 
    movie->mv_gv_dir = pathgv;
    movie->mv_vtk_dir = pathvtk; 
    movie->out_dir = oname;
    movie->doMakeMovie = true;
    movie->mv_intfc = this->m_intfc;
    std::string sys_cmd = "mkdir -p " + pathgv;
    system(sys_cmd.c_str());
    sys_cmd = "mkdir -p " + pathvtk; 
    system(sys_cmd.c_str());
}

void Movie::recordMovieFrame() {
    std::string fname = mv_gv_dir + "/intfc-" + std::to_string(mv_count);
    gview_plot_interface(fname.c_str(),mv_intfc);

    // keep marking number to have fixed number of digits
    std::string count; 
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << mv_count;
    count = ss.str();
    fname = mv_vtk_dir + "/vtk-" + count;
    std::string sys_cmd = "mkdir -p " + fname; 
    system(sys_cmd.c_str());
    vtk_interface_plot(fname.c_str(),mv_intfc, YES,
                                        mv_dt * mv_count ,mv_count);
    mv_count++; 
}

bool Movie::isMovieTime(double t) {
    if (doMakeMovie)
        return fabs(t-mv_count*mv_dt) < 1e-10;
    else
	return false;
}

static void setCollisionFreePoints3d(INTERFACE* intfc, Drag* drag)
{
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        assert (intfc->dim == 3);

        next_point(intfc,NULL,NULL,NULL);
        while(next_point(intfc,&p,&hse,&hs)){
            STATE* sl = (STATE*)left_state(p);
            sl->is_fixed = false;
            if (wave_type(hs) == NEUMANN_BOUNDARY ||
                wave_type(hs) == MOVABLE_BODY_BOUNDARY)
            {
                sl->is_fixed = true;
            }
        }

        CURVE **c;
        BOND* b;
        intfc_curve_loop(intfc,c){
            if (hsbdry_type(*c) != FIXED_HSBDRY)
                continue;
            for (b = (*c)->first; b != (*c)->last; b = b->next)
            {
                STATE* sl = (STATE*)left_state(b->end);
                sl->is_fixed = true;
            }
        }
}       /* setCollisionFreePoints3d() */

Folder3d::Folder3d(INTERFACE* intfc, SURFACE* s) : m_intfc(intfc)
{}

Folder3d::~Folder3d() {
    if (movie) delete movie;
}

static int numConnectCurves(NODE* n) {
    int ans = 0;
    for (CURVE** c = n->in_curves; c && *c; ++c)
	ans++;

    for (CURVE** c = n->out_curves; c && *c; ++c)
	ans++;
    return ans;
}

void Folder3d::straightenStrings() {
    std::cout << "Straighten the string" << std::endl;
    //straighten the string curves
    CURVE** c;
    double z_avg = 0;
    int c_count = 1;
    intfc_curve_loop(m_intfc,c) {
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	//count how many bonds in a curve
        int count = 0;
        BOND* b;
	curve_bond_loop(*c,b)
	    count++;

	//adjust load node
	POINT* p1, *p2;
	p1 = (*c)->start->posn;
	p2 = (*c)->end->posn;
	if (numConnectCurves((*c)->start) < numConnectCurves((*c)->end))
	    std::swap(p1,p2); //make sure p1 is load node;

	double l = (*c)->first->length0 * count;
	double new_z  =  Coords(p2)[2] - sqrt(sqr(l) 
		       - sqr(Coords(p1)[0]-Coords(p2)[0]) 
		       - sqr(Coords(p1)[1]-Coords(p2)[1]));
	z_avg = z_avg + (new_z - z_avg)/(c_count++);
	Coords(p1)[2] = z_avg;
    }

    intfc_curve_loop(m_intfc,c) {
        if (hsbdry_type(*c) != STRING_HSBDRY) continue;
	//count how many bonds in a curve
        int count = 0;
        BOND* b;
	curve_bond_loop(*c,b)
	    count++;

        double dir[3] = {0};
        for (size_t i = 0; i < 3; ++i)
            dir[i] = Coords((*c)->end->posn)[i] -
                     Coords((*c)->start->posn)[i];
        //normalize dir
        double len = Mag3d(dir);
        for (size_t i = 0; i < 3; ++i)
            dir[i] /= len;

        //move bonds
        double len0 = len/count;
        count = 1;
        for (b = (*c)->first; b != (*c)->last; b = b->next) {
            POINT* p = b->end;
            for (size_t i = 0; i < 3; ++i)
                Coords(p)[i] = Coords((*c)->start->posn)[i]
                             + (count)*len0 * dir[i];
            count++;
        }
   }
}

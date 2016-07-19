#include "spring_solver.h"
#include <assert.h>
#include <iterator>
#include <fstream>
#include <iomanip>

//class SpringVeterx
void SpringVertex::addNeighbor(size_t i_nb, double len0) 
{
    index_nb.push_back(i_nb);
    length0.push_back(len0);
}

//class SpringSolver
void SpringSolver::printAdjacencyList(std::string fname) {
    std::ofstream ofs(fname.c_str());
    for (size_t i = 0; i < pts.size(); ++i) {
	ofs << i << " ";
	std::copy(pts[i]->index_nb.begin(),pts[i]->index_nb.end(),
		std::ostream_iterator<size_t>(ofs," "));
	ofs << "\n";
    }
}

void SpringSolver::printPointList(std::string fname) {
    printf("%lu number of points in list\n",pts.size());
    std::ofstream ofs(fname.c_str());
    ofs << std::setprecision(9);
    for (size_t i = 0; i < pts.size(); ++i) {
	double* crds = pts[i]->getCoords();
	std::copy(crds,crds+3,std::ostream_iterator<double>(ofs," "));
	ofs << pts[i]->getPoint()<< "\n";
    }
}

void SpringSolver::checkVertexNeighbors() {
    for (size_t i = 0; i < pts.size(); ++i) {
	SpringVertex* sv = pts[i];
	for (size_t j = 0; j < sv->index_nb.size(); ++j) {
	    SpringVertex* sv_nb = pts[static_cast<size_t>(sv->index_nb[j])];
	    sv_nb->addNeighbor(i,sv->length0[j]);
	}
    }
}
void SpringSolver::setParameters(double k, double lambda, double m) {
    springParameter.k = k;
    springParameter.lambda = lambda;
    springParameter.m = m;
}

void SpringSolver::computeAccel(SpringVertex* sv) {
	double k = springParameter.k;
	double m = springParameter.m;
	double lambda = springParameter.lambda;
	const int dim = 3;
	double vec[3] = {0};
	double v_rel[3] = {0};

	if (sv->isRegistered())
	{
	    setPresetVelocity(sv);
	    std::copy(sv->f_ext,sv->f_ext+3,sv->f);
	    return;
	}

        for (int i = 0; i < dim; ++i)
            sv->f[i] = 0.0;

        for (size_t i = 0; i < sv->index_nb.size(); ++i)
        {
            double len = 0.0;
	    SpringVertex* sv_nb = pts[sv->index_nb[i]];
            for (int j = 0; j < dim; ++j)
            {
                vec[j] = sv_nb->x[j] - sv->x[j];
                len += sqr(vec[j]);
                v_rel[j] = sv_nb->v[j] - sv->v[j];
            }
            len = sqrt(len);
	    for (int j = 0; j < dim; ++j)
            {
		if (len > MACH_EPS)
                    vec[j] /= len;
		else
		    vec[j] = 0.0;
                sv->f[j] += k*((len - sv->length0[i])*vec[j])/m;
                sv->f[j] += lambda*v_rel[j]/m;
            }
	    if (sv->length0[i] < 0) clean_up(ERROR);
        }
	
        for (int j = 0; j < dim; ++j)
        {
            sv->f[j] += (sv->getExternalAccel())[j];
        }
}



void SpringSolver::rk4SpringSolver(double dt, int n_loop) {
    const size_t size = pts.size();
    static std::vector<std::vector<double> > x_old(size,std::vector<double>(3,0));
    static std::vector<std::vector<double> > x_new(x_old);
    static std::vector<std::vector<double> > v_old(x_old);
    static std::vector<std::vector<double> > v_new(x_old);

    size_t i, j;
    const int dim = 3;

    printf("Starting spring solver:\n");
    for (int n = 0; n < n_loop; ++n)
    {
	printf("    #sub_step  = %d/%d\n",n+1,n_loop);
    	for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
	for (j = 0; j < dim; ++j)
        {
	    x_old[i][j] = pts[i]->x[j];
	    v_old[i][j] = pts[i]->v[j];
        }

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
    	    x_new[i][j] = x_old[i][j] + dt*pts[i]->v[j]/6.0;
            v_new[i][j] = v_old[i][j] + dt*pts[i]->f[j]/6.0;
            pts[i]->x[j] = x_old[i][j] + 0.5*pts[i]->v[j]*dt;
            pts[i]->v[j] = v_old[i][j] + 0.5*pts[i]->f[j]*dt;
        }

        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
    	    x_new[i][j] += dt*pts[i]->v[j]/3.0;
            v_new[i][j] += dt*pts[i]->f[j]/3.0;
            pts[i]->x[j] = x_old[i][j] + 0.5*pts[i]->v[j]*dt;
            pts[i]->v[j] = v_old[i][j] + 0.5*pts[i]->f[j]*dt;
        }
        
        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
    	    x_new[i][j] += dt*pts[i]->v[j]/3.0;
            v_new[i][j] += dt*pts[i]->f[j]/3.0;
            pts[i]->x[j] = x_old[i][j] + pts[i]->v[j]*dt;
            pts[i]->v[j] = v_old[i][j] + pts[i]->f[j]*dt; 
        }

        for (i = 0; i < size; ++i)
            computeAccel(pts[i]);

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
    	    x_new[i][j] += dt*pts[i]->v[j]/6.0;
            v_new[i][j] += dt*pts[i]->f[j]/6.0;
        }

        for (i = 0; i < size; ++i)
        for (j = 0; j < dim; ++j)
        {
    	    pts[i]->x[j] = x_new[i][j];
            pts[i]->v[j] = v_new[i][j];
        }
    }
}

void SpringSolver::doSolve(double t) {
    double dt = getTimeStepSize(); 
    int nt = t/dt+1;
    dt = t/nt;
    rk4SpringSolver(dt,nt);
}

double SpringSolver::getTimeStepSize() {
    return 0.1*sqrt(springParameter.m
		/springParameter.k);
}

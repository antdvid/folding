#include "folding.h"
#include "spring_solver.h"
#include "ft_spring_solver.h"

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
}

void Folder3d::doFolding(Drag* drag, FT_SpringSolver* sp_solver) {
    sp_solver->setDrag(drag);
    sp_solver->doSolve(drag->m_t);
}

Folder3d::Folder3d(INTERFACE* intfc, SURFACE* s) : m_intfc(intfc)
{}

Folder3d::~Folder3d() {}

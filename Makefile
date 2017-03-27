CXX=mpicxx -g -rdynamic -fopenmp -pedantic -Wno-long-long -Wno-unused-result -std=c++11

incs =   -I.. -I../util   -I/usr/local/pkg/hdf/include -DUSE_HDF  -I/usr/local/pkg/mpich2/include -D__MPI__   -DHAS_FENV    -D__GD__ -I/usr/include  -D__HYPRE__   -I../include -D__DAMPING__
libincs =  -L/usr/lib  /usr/local/pkg/petsc/lib/libfblas.a /usr/lib/x86_64-linux-gnu/libgfortran.so.3 -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib  -L../lib/x86_64 -L.
libs = -lgd  /usr/local/pkg/petsc/lib/libfblas.a /usr/lib/x86_64-linux-gnu/libgfortran.so.3  -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread -lmpich -lpthread -lmpfr

CGAL_DIR=
CVODE_DIR=/usr/local/pkg/cvode

CGAL_Include=
CGAL_Lib= -lCGAL_Core -lCGAL_ImageIO -lCGAL

CVODE_Include=-I/usr/local/pkg/cvode/include
CVODE_Lib=-L/usr/local/pkg/cvode/lib -lsundials_cvode -lsundials_nvecserial

NLOPT_Include=-I./lib
NLOPT_Lib=-L./lib -lnlopt

test: test.o folding.o folding_helper.o drag.o dcollid3d.o dcollid.o spring_solver.o drag_proto.o ex_spring_solver.o im_spring_solver.o bending.o cgal.o origami.o
	$(CXX) $^ -lFronTier -lm -o test $(libincs) $(libs) $(incs) $(CGAL_Include) $(CGAL_Lib) $(CVODE_Include) $(CVODE_Lib) $(NLOPT_Lib) -lgmp -lmpfr -frounding-math
dcollid3d.o: ../Collision/dcollid3d.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
dcollid.o: ../Collision/dcollid.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
cgal.o: cgal.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -L/usr/lib64 -lgmp -lmpfr -frounding-math
test.o: test.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
im_spring_solver.o: im_spring_solver.cpp
	$(CXX) $< -c $(incs) $(CVODE_Include) -frounding-math
origami.o: origami.cpp
	$(CXX) $< -c $(incs) $(NLOPT_Include)
%.o: %.cpp
	$(CXX) $< -c $(incs) $(NLOPT_Include) $(CGAL_Include) -frounding-math

-include ../devel-deps.inc

clean:
	rm -rf *.o test 
tagsfile:
	ctags *.h *.cpp ../Collision/*.h ../Collision/*.cpp ../src/*/*.[chf]


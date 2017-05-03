CXX=mpicxx -g -rdynamic -fopenmp -pedantic -Wno-long-long -Wno-unused-result -std=c++11

incs =   -I.. -I../util   -I/usr/local/pkg/hdf/include -DUSE_HDF   -DHAS_FENV    -D__MPI__  -D__GD__ -I/usr/local/pkg/gd/include  -D__HYPRE__   -I../include -D__DAMPING__
libincs =  -L/usr/local/pkg/gd/lib  /lib64/libblas.so /lib64/libgfortran.so.3 -L/usr/lib -L/usr/local/pkg/hdf/lib  -L../lib/x86_64 -L.
libs = -lgd  /lib64/libblas.so /lib64/libgfortran.so.3  -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread -lmpfr

CGAL_DIR=/usr/local/pkg/cgal
CVODE_DIR=/usr/local/pkg/cvode

CGAL_Include=-I$(CGAL_DIR)/include
CGAL_Lib=-L$(CGAL_DIR)/lib -lCGAL_Core -lCGAL_ImageIO -lCGAL

CVODE_Include=-I$(CVODE_DIR)/include
CVODE_Lib=-L$(CVODE_DIR)/lib -lsundials_cvode -lsundials_nvecserial

NLOPT_Include=-I./lib
NLOPT_Lib=-L./lib -lnlopt

test: test.o folding.o folding_helper.o drag.o dcollid3d.o dcollid.o spring_solver.o drag_proto.o ex_spring_solver.o im_spring_solver.o bending.o cgal.o origami.o
	$(CXX) $^ -lFronTier -lm -o test $(libincs) $(libs) $(CGAL_Lib) $(CVODE_Lib) $(NLOPT_Lib) -lgmp -lmpfr -frounding-math
dcollid3d.o: ../Collision/dcollid3d.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
dcollid.o: ../Collision/dcollid.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
cgal.o: cgal.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
test.o: test.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) -frounding-math
im_spring_solver.o: im_spring_solver.cpp
	$(CXX) $< -c $(incs) $(CVODE_Include) -frounding-math
origami.o: origami.cpp
	$(CXX) $< -c $(incs) $(NLOPT_Include)
%.o: %.cpp
	$(CXX) $< -c $(incs) $(CGAL_Include) $(NLOPT_Include) -frounding-math

-include ../devel-deps.inc

clean:
	rm -rf *.o test 
tagsfile:
	ctags *.h *.cpp ../Collision/*.h ../Collision/*.cpp ../src/*/*.[chf]


CXX=mpicxx -g -rdynamic -std=c++11
#uncomment the following if debugging needed
CXX=mpicxx -g -pg -std=c++11 -pedantic -Wall -Wextra -Wno-undef -Wno-comment -Wno-unused-parameter -Wno-long-long -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wno-old-style-cast -Wno-redundant-decls

incs =   -I/usr/local/pkg/mpich2/include -D__MPI__  -I/usr/local/pkg/hdf/include -DUSE_HDF   -D__HYPRE__ -I../include -I/usr/local/pkg/cvode/include
libincs =  -L/usr/lib -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib -L../lib/x86_64 -L/usr/local/pkg/cvode/lib
libs = -lgd -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread -lsundials_cvode -lsundials_nvecserial

test: test.o folding.o folding_helper.o drag.o dcollid3d.o dcollid.o spring_solver.o drag_proto.o ex_spring_solver.o im_spring_solver.o bending.o
	$(CXX) $^ -lFronTier -lm -o test $(libincs) $(libs) -lCGAL_Core -lCGAL_ImageIO -lCGAL -lgmp -frounding-math

dcollid3d.o: ../Collision/dcollid3d.cpp
	$(CXX) $< -c -I../include $(incs) -frounding-math
dcollid.o: ../Collision/dcollid.cpp
	$(CXX) $< -c -I../include $(incs) -frounding-math
%.o: %.cpp
	$(CXX) $< -c -I../include $(incs) -frounding-math

-include ../devel-deps.inc

clean:
	rm -rf *.o test
tagsfile:
	ctags *.cpp ../src/*/*.[chf]

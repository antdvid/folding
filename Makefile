#CXX=mpicxx -g -rdynamic -std=c++11
#uncomment the following if debugging needed
CXX=mpicxx -g -std=c++11 -pedantic -Wall -Wextra -Wno-undef -Wno-comment -Wno-unused-parameter -Wno-long-long -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Werror -Wno-old-style-cast -Wno-redundant-decls

incs =   -I/usr/local/pkg/mpich2/include -D__MPI__  -I/usr/local/pkg/hdf/include -DUSE_HDF   -D__HYPRE__ -I../include
libincs =  -L/usr/lib -L/usr/lib -L/usr/local/pkg/hdf/lib  -L/usr/local/pkg/mpich2/lib -L../lib/x86_64
libs = -lgd -lmfhdf -ldf -ljpeg -lz  -lmpich -lpthread

test: test.o folding.o drag.o dcollid3d.o dcollid.o spring_solver.o ft_spring_solver.o
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

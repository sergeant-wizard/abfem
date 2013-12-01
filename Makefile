EXE			= ABFEM
CC			= g++
#CC			= icpc
OBJS_OTHERS	= matrix_wiz.o sparse.o gid.o global.o matrix_algebraic.o
OBJS_ELEM	= tetra1.o tetra8.o elemface.o hexa27.o
OBJS_MARKER	= marker.o tetra1_marker.o face.o
OBJS		= $(OBJS_OTHERS) $(OBJS_ELEM) $(OBJS_MARKER)
OBJS_DIR	= objs

ALL_CPP_FILES=elemface.cpp gid.cpp global.cpp main.cpp matrix_wiz.cpp sparse.cpp tetra1.cpp tetra1_marker.cpp tetra8.cpp marker.cpp matrix_algebraic.cpp hexa27.cpp

#FLAGS		= -openmp
#MAIN_FLAGS	= -L$(MKLPATH) $(MKLPATH)/libmkl_solver_lp64.a -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread -gxx-name=g++-4.3

#MKLPATH		= /opt/intel/mkl/10.2.5.035/lib/em64t
#MKLPATH		= /opt/intel-local/mkl/lib/intel64/
#MKLPATH		= /opt/intel/mkl/10.2.5.035/lib/64

ABFEM: $(OBJS)
	$(CC) -o ../$@ $(OBJS) $(MAIN_FLAGS) main.cpp;

.cpp.o:
	$(CC) -c $(<F) $(FLAGS)

clean:
	-rm ./objs/*.o ./*.o ../$(EXE)

depend:
	makedepend -- $(FLAGS) -- $(ALL_CPP_FILES)

# DO NOT DELETE

elemface.o: elemface.h /usr/include/math.h /usr/include/sys/reent.h
elemface.o: /usr/include/_ansi.h /usr/include/newlib.h
elemface.o: /usr/include/sys/config.h /usr/include/machine/ieeefp.h
elemface.o: /usr/include/sys/features.h /usr/include/sys/_types.h
elemface.o: /usr/include/machine/_types.h
elemface.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
elemface.o: common.h tetra1_up.hpp tetra1.h element.h sparse.h
elemface.o: /usr/include/stdlib.h /usr/include/machine/stdlib.h
elemface.o: /usr/include/alloca.h matrix.h node.h material.h
elemface.o: matrix_algebraic.h
gid.o: gid.h prepost.h element.h sparse.h /usr/include/stdlib.h
gid.o: /usr/include/machine/ieeefp.h /usr/include/_ansi.h
gid.o: /usr/include/newlib.h /usr/include/sys/config.h
gid.o: /usr/include/sys/features.h /usr/include/sys/reent.h
gid.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
gid.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
gid.o: /usr/include/machine/stdlib.h /usr/include/alloca.h
gid.o: /usr/include/math.h matrix.h node.h tetra8.h tetra1_up.hpp tetra1.h
gid.o: common.h material.h matrix_algebraic.h gid.hpp
global.o: global.h common.h matrix_wiz.h /usr/include/math.h
global.o: /usr/include/sys/reent.h /usr/include/_ansi.h /usr/include/newlib.h
global.o: /usr/include/sys/config.h /usr/include/machine/ieeefp.h
global.o: /usr/include/sys/features.h /usr/include/sys/_types.h
global.o: /usr/include/machine/_types.h /usr/include/machine/_default_types.h
global.o: /usr/include/sys/lock.h matrix.h sparse.h /usr/include/stdlib.h
global.o: /usr/include/machine/stdlib.h /usr/include/alloca.h element.h
global.o: node.h prepost.h
main.o: /usr/include/stdlib.h /usr/include/machine/ieeefp.h
main.o: /usr/include/_ansi.h /usr/include/newlib.h /usr/include/sys/config.h
main.o: /usr/include/sys/features.h /usr/include/sys/reent.h
main.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
main.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
main.o: /usr/include/machine/stdlib.h /usr/include/alloca.h
main.o: /usr/include/math.h incompressible.hpp material.h common.h
main.o: matrix_algebraic.h matrix.h elastic.hpp global_pstab.hpp
main.o: global_up.hpp global.h matrix_wiz.h sparse.h element.h node.h
main.o: prepost.h elemface.h tetra1_up.hpp tetra1.h gid.h tetra8.h
main.o: gid_marker.hpp marker.h face.h elem_gauss.h inc_marker.hpp
main.o: tetra1_marker.h global_marker.hpp gid.hpp hexa27.h elem_gauss.hpp
main.o: elem_gauss_up.hpp elem_gauss_marker.hpp bodyforce.hpp
matrix_wiz.o: matrix_wiz.h /usr/include/math.h /usr/include/sys/reent.h
matrix_wiz.o: /usr/include/_ansi.h /usr/include/newlib.h
matrix_wiz.o: /usr/include/sys/config.h /usr/include/machine/ieeefp.h
matrix_wiz.o: /usr/include/sys/features.h /usr/include/sys/_types.h
matrix_wiz.o: /usr/include/machine/_types.h
matrix_wiz.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
matrix_wiz.o: matrix.h
sparse.o: sparse.h /usr/include/stdlib.h /usr/include/machine/ieeefp.h
sparse.o: /usr/include/_ansi.h /usr/include/newlib.h
sparse.o: /usr/include/sys/config.h /usr/include/sys/features.h
sparse.o: /usr/include/sys/reent.h /usr/include/sys/_types.h
sparse.o: /usr/include/machine/_types.h /usr/include/machine/_default_types.h
sparse.o: /usr/include/sys/lock.h /usr/include/machine/stdlib.h
sparse.o: /usr/include/alloca.h /usr/include/math.h
tetra1.o: tetra1.h common.h element.h sparse.h /usr/include/stdlib.h
tetra1.o: /usr/include/machine/ieeefp.h /usr/include/_ansi.h
tetra1.o: /usr/include/newlib.h /usr/include/sys/config.h
tetra1.o: /usr/include/sys/features.h /usr/include/sys/reent.h
tetra1.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
tetra1.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
tetra1.o: /usr/include/machine/stdlib.h /usr/include/alloca.h
tetra1.o: /usr/include/math.h matrix.h node.h material.h matrix_algebraic.h
tetra1_marker.o: tetra1_marker.h common.h tetra1.h element.h sparse.h
tetra1_marker.o: /usr/include/stdlib.h /usr/include/machine/ieeefp.h
tetra1_marker.o: /usr/include/_ansi.h /usr/include/newlib.h
tetra1_marker.o: /usr/include/sys/config.h /usr/include/sys/features.h
tetra1_marker.o: /usr/include/sys/reent.h /usr/include/sys/_types.h
tetra1_marker.o: /usr/include/machine/_types.h
tetra1_marker.o: /usr/include/machine/_default_types.h
tetra1_marker.o: /usr/include/sys/lock.h /usr/include/machine/stdlib.h
tetra1_marker.o: /usr/include/alloca.h /usr/include/math.h matrix.h node.h
tetra1_marker.o: material.h matrix_algebraic.h tetra1_up.hpp marker.h face.h
tetra1_marker.o: elem_gauss.h
tetra8.o: gid.h prepost.h element.h sparse.h /usr/include/stdlib.h
tetra8.o: /usr/include/machine/ieeefp.h /usr/include/_ansi.h
tetra8.o: /usr/include/newlib.h /usr/include/sys/config.h
tetra8.o: /usr/include/sys/features.h /usr/include/sys/reent.h
tetra8.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
tetra8.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
tetra8.o: /usr/include/machine/stdlib.h /usr/include/alloca.h
tetra8.o: /usr/include/math.h matrix.h node.h tetra8.h tetra1_up.hpp tetra1.h
tetra8.o: common.h material.h matrix_algebraic.h gid.hpp
marker.o: marker.h element.h sparse.h /usr/include/stdlib.h
marker.o: /usr/include/machine/ieeefp.h /usr/include/_ansi.h
marker.o: /usr/include/newlib.h /usr/include/sys/config.h
marker.o: /usr/include/sys/features.h /usr/include/sys/reent.h
marker.o: /usr/include/sys/_types.h /usr/include/machine/_types.h
marker.o: /usr/include/machine/_default_types.h /usr/include/sys/lock.h
marker.o: /usr/include/machine/stdlib.h /usr/include/alloca.h
marker.o: /usr/include/math.h matrix.h node.h matrix_algebraic.h face.h
marker.o: common.h material.h elem_gauss.h
matrix_algebraic.o: matrix_algebraic.h matrix.h
hexa27.o: hexa27.h elem_gauss.hpp elem_gauss.h common.h element.h sparse.h
hexa27.o: /usr/include/stdlib.h /usr/include/machine/ieeefp.h
hexa27.o: /usr/include/_ansi.h /usr/include/newlib.h
hexa27.o: /usr/include/sys/config.h /usr/include/sys/features.h
hexa27.o: /usr/include/sys/reent.h /usr/include/sys/_types.h
hexa27.o: /usr/include/machine/_types.h /usr/include/machine/_default_types.h
hexa27.o: /usr/include/sys/lock.h /usr/include/machine/stdlib.h
hexa27.o: /usr/include/alloca.h /usr/include/math.h matrix.h node.h
hexa27.o: material.h matrix_algebraic.h elem_gauss_up.hpp marker.h face.h
hexa27.o: elem_gauss_marker.hpp

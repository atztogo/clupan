cpp = g++

#cflags_dndebug = -O3 -DBOOST_UBLAS_NO_STD_CERR -DNDEBUG
#cflags = -O3 -DBOOST_UBLAS_NO_STD_CERR
cflags = -O3 -DBOOST_UBLAS_NO_STD_CERR -I./src 
ompflags = -fopenmp

loader = $(cpp)


# spglib library
libspglib = ~/usr/local/spglib-1.9.9/lib/libsymspg.a
# openmp library
libomp = -lgomp
# GSL library
dirgsl = /usr/lib/x86_64-linux-gnu
libgsl = $(dirgsl)/libgsl.a $(dirgsl)/libgslcblas.a

VPATH = src

source_io = input.o output.o error_check.o parse_input.o \
    parse_structure.o structure.o smith.o output_structure.o \
    derivative_lattice.o
source_derivative = sym.o hermite.o permutation.o \
    permutation_atom.o labeling.o main_derivative.o 
source_derivative_to_structure = derivative_to_structure.o 
source_derivative_to_poscar = derivative_to_poscar.o 
source_cluster_search =  sym.o combination.o permutation_atom.o \
    main_cluster_search.o
source_correlation =  sym.o permutation_atom.o \
    hamiltonian.o permutation.o correlation.o main_correlation.o
source_lsf = least_squares.o cv.o main_lsf.o
source_cv_casp = least_squares.o cv.o main_cv_casp.o
source_ga = least_squares.o cv.o combination.o ga.o main_ga.o
source_gss = main_gss.o
source_mc = sym.o permutation_atom.o hamiltonian.o permutation.o \
    correlation.o ewald.o mc_initialize.o mc_common.o \
    mc_binary.o gcmc_binary.o mc.o main_mc.o
source_ewald = ewald.o main_ewald.o
source_ti = main_ti.o
source_sqs = sym.o permutation_atom.o hamiltonian.o permutation.o \
    correlation.o parse_poscar.o mc_initialize.o mc_common.o \
    mc_binary.o mc.o ewald.o main_sqs.o

source_poscar_to_structure = parse_poscar.o poscar_to_structure.o
source_structure_to_poscar = structure_to_poscar.o
source_find_primitive = find_primitive.o
source_refine_cell = refine_cell.o
source_mc_cell_check = mc_cell_check.o
source_correlation_for_casp = correlation_for_casp.o
source_correlation_for_convex = correlation_for_convex.o
source_structure_selection = structure_selection.o
source_modify_structure = modify_structure.o
source_opt_ewald = ewald.o least_squares.o cv.o opt_ewald.o

all: mkdir $(addprefix lib/,$(source_io)) \
    $(addprefix lib/,$(source_derivative)) \
    $(addprefix lib/,$(source_derivative_to_structure)) \
    $(addprefix lib/,$(source_derivative_to_poscar)) \
    $(addprefix lib/,$(source_cluster_search)) \
    $(addprefix lib/,$(source_correlation)) \
    $(addprefix lib/,$(source_lsf)) \
    $(addprefix lib/,$(source_cv_casp)) \
    $(addprefix lib/,$(source_ga)) \
    $(addprefix lib/,$(source_gss)) \
    $(addprefix lib/,$(source_mc)) \
    $(addprefix lib/,$(source_ewald)) \
    $(addprefix lib/,$(source_ti)) \
    $(addprefix lib/,$(source_sqs)) \
    $(addprefix lib/,$(source_poscar_to_structure)) \
    $(addprefix lib/,$(source_structure_to_poscar)) \
    $(addprefix lib/,$(source_find_primitive)) \
    $(addprefix lib/,$(source_refine_cell)) \
    $(addprefix lib/,$(source_mc_cell_check)) \
    $(addprefix lib/,$(source_correlation_for_casp)) \
    $(addprefix lib/,$(source_correlation_for_convex)) \
    $(addprefix lib/,$(source_structure_selection)) \
    $(addprefix lib/,$(source_modify_structure)) \
    $(addprefix lib/,$(source_opt_ewald))
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_derivative)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/derivative
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_derivative_to_structure)) \
        $(libboost) $(libgsl) $(libomp) -o bin/derivative_to_structure
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_derivative_to_poscar)) \
        $(libboost) $(libgsl) $(libomp) -o bin/derivative_to_poscar
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_cluster_search)) \
        $(libboost) $(libgsl) $(libspglib) $(libomp) -o bin/cluster_search
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_correlation)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/correlation
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_lsf)) \
        $(libboost) $(libgsl) $(libomp) -o bin/lsf
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_cv_casp)) \
        $(libboost) $(libgsl) $(libomp) -o bin/cv_casp
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_ga)) \
        $(libboost) $(libgsl) $(libomp) -o bin/ga
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_gss)) \
        $(libboost) $(libgsl) $(libomp) -o bin/gss
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_mc)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/mc
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_ewald)) \
        $(libboost) $(libgsl) $(libomp) -o bin/ewald
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_ti)) \
        $(libboost) $(libgsl) $(libomp) -o bin/ti
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_sqs)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/sqs
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_poscar_to_structure)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/poscar_to_structure
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_structure_to_poscar)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/structure_to_poscar
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_find_primitive)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/find_primitive
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_refine_cell)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/refine_cell
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_mc_cell_check)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/mc_cell_check
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_correlation_for_casp)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/correlation_for_casp
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_correlation_for_convex)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/correlation_for_convex
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_structure_selection)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/structure_selection
	$(loader) $(ldflags) $(addprefix lib/,$(source_io)) \
        $(addprefix lib/,$(source_modify_structure)) \
        $(libboost) $(libspglib) $(libgsl) $(libomp) -o bin/modify_structure

mkdir: 
	@if [ ! -d bin ]; then \
        echo "mkdir bin"; mkdir bin; \
    fi 
	@if [ ! -d lib ]; then \
        echo "mkdir lib"; mkdir lib; \
    fi 

lib/%.o: %.cpp
	$(cpp) -c $(cflags) $(ompflags) -o $@ $<

clean:
	rm -rf bin lib


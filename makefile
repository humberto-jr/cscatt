# How to:
#
# 1) If both dependencies, Intel MKL and GSL, are installed in their
# default location, the user defined potential energy surface (PES)
# is pre-compiled and named pot.o (default), then type
#
# make
#
# to build and link all codes using GNU gcc compiler (default). Notice,
# the environment variable MKL_DIR must be properly defined before make.
#
# 2) If Intel icc compiler is the C compiler of choice instead, type
#
# make CC=icc
#
# 3) If MPI functionalities are to be included in the build, type
#
# make CC=icc with_mpi
#
# or, "make with_mpi" to use the default GNU compiler.
#
# 4) If the user defined PES happens to have a different filename, then
# type
#
# make CC=icc PES_LIB=[filename.o]
#
# or, "make PES_LIB=[filename.o]". And, if the PES also depends on another
# set of object files,
#
# make PES_LIB="[filename.o] [dependence_a.o] [dependence_b.o] [etc.o]"
#
# 5) If the GSL library is not installed in the default directory, type
#
# make GSL_INC=[path]/include GSL_LIB=[path]/bin/libgsl.a
#
# where, [path] points the location where the library is installed. As
# usual the include directory with all header files and the library as
# an .a binary file are needed. Often when working at machines in which
# no root permissions are available one is forced to install GSL
# somewhere else or asking the administrator to install it.
#
# 7) For fresh builds type
#
# make clean
#
# In particular a good practice would be to run "make clean" between the
# building process of routines with and without MPI.
#
# Humberto Jr
# Dec, 2018
#

SHELL = /bin/sh
LINEAR_ALGEBRA = GSL
GSL_DIR = /usr/local

#
# C compilers (GNU gcc by default):
#

CC = gcc
CFLAGS = -W -Wall -std=c99 -pedantic -fopenmp -O3 -I$(GSL_DIR)/include
LDFLAGS = -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm

#
# GNU or Intel MPI wrappers:
#

USE_MPI = no

ifeq ($(CC), mpicc)
	USE_MPI = yes
endif

ifeq ($(CC), mpiicc)
	USE_MPI = yes
endif

#
# IBM xlc compiler:
#

XLF_DIR = /opt/ibmcmp/xlf/bg/14.1

ifeq ($(CC), xlc)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
	override LDFLAGS = -lxlf90_r -lxl -lxlfmath -L$(XLF_DIR)/bglib64
endif

#
# IBM xlc_r compiler:
#

ifeq ($(CC), xlc_r)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
	override LDFLAGS = -lxlf90_r -lxl -lxlfmath -L$(XLF_DIR)/bglib64
endif

#
# IBM mpixlc compiler:
#

ifeq ($(CC), mpixlc)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
	override LDFLAGS = -lxlf90_r -lxl -lxlfmath -L$(XLF_DIR)/bglib64
	USE_MPI = yes
endif

#
# IBM mpixlc_r compiler:
#

ifeq ($(CC), mpixlc_r)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
	override LDFLAGS = -lxlf90_r -lxl -lxlfmath -L$(XLF_DIR)/bglib64
	USE_MPI = yes
endif

#
# PGI pgcc compiler:
#

ifeq ($(CC), pgcc)
	override CFLAGS = -mp -O3 -I$(GSL_DIR)/include
endif

#
# Fortran compilers (using GNU gfortran as default option):
#

FC =

ifeq ($(FC), gfortran)
	LDFLAGS += -lgfortran
endif

ifeq ($(FC), ifort)
	LDFLAGS += -lifcore
endif

ifeq ($(FC), )
	FC = gfortran
endif

#
# User defined potential energy surface (PES):
#

PES_NAME =
PES_OBJECT =

ifeq ($(PES_NAME), )
	PES_MACRO =
else
	ifeq ($(PES_OBJECT), )
		PES_OBJECT = $(PES_NAME).o
	endif

	PES_MACRO = -DEXTERNAL_PES_NAME=$(PES_NAME)
endif

#
# MPI library (using OpenMPI as default option):
#

MPI_DIR = /usr/include/openmpi

ifeq ($(USE_MPI), yes)
	MPI_INC = -DUSE_MPI
endif

#
# Intel MKL (for GNU gcc or Intel icc compilers):
#

MKL_DIR = $(MKLROOT)

ifeq ($(LINEAR_ALGEBRA), MKL)
	LINEAR_ALGEBRA_INC = -DMKL_ILP64 -m64 -I$(MKL_DIR)/include -DUSE_MKL

	ifeq ($(CC), icc)
		LINEAR_ALGEBRA_LIB = -parallel -Wl,--start-group $(MKL_DIR)/lib/intel64/libmkl_intel_ilp64.a $(MKL_DIR)/lib/intel64/libmkl_intel_thread.a $(MKL_DIR)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
	endif

	ifeq ($(CC), gcc)
		LINEAR_ALGEBRA_LIB = -Wl,--start-group $(MKL_DIR)/lib/intel64/libmkl_intel_ilp64.a $(MKL_DIR)/lib/intel64/libmkl_gnu_thread.a $(MKL_DIR)/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
	endif
endif

#
# LAPACKE, LAPACK and BLAS libraries:
#

LAPACKE_DIR = /usr/local
LAPACK_DIR = /usr/local
BLAS_DIR = /usr/local

ifeq ($(LINEAR_ALGEBRA), LAPACKE)
	LINEAR_ALGEBRA_INC = -I$(LAPACKE_DIR)/include -DUSE_LAPACKE
	LINEAR_ALGEBRA_LIB = -L$(LAPACKE_DIR)/lib -llapacke -llapack -lblas -lgfortran -lm
endif

#
# IBM ESSL (for Blue Gene machines):
#

ESSL_DIR = /usr

ifeq ($(LINEAR_ALGEBRA), ESSL)
	LINEAR_ALGEBRA_INC = -I$(ESSL_DIR)/include -DUSE_ESSL
	LINEAR_ALGEBRA_LIB = -L$(ESSL_DIR)/lib64 -lesslbg -lm
endif

#
# MAGMA library:
#

MAGMA_DIR = /usr/local/magma
CUDA_DIR = /usr/local/cuda

ifeq ($(LINEAR_ALGEBRA), MAGMA)
	CFLAGS += -DUSE_MAGMA
	LINEAR_ALGEBRA_INC = -I$(CUDA_DIR)/include -I$(MAGMA_DIR)/include -DADD_
	LINEAR_ALGEBRA_LIB = -L$(MAGMA_DIR)/lib -L$(CUDA_DIR)/lib64 -lmagma -lm
endif

#
# ATLAS library:
#

ATLAS_DIR = /usr/local/atlas

ifeq ($(LINEAR_ALGEBRA), ATLAS)
	LINEAR_ALGEBRA_INC = -I$(ATLAS_DIR)/include -DUSE_ATLAS
	LINEAR_ALGEBRA_LIB = -L$(ATLAS_DIR)/lib -latlas -lptcblas -lm
endif

#
# PETSc library (optionally used by mpi_lib module and requires BLAS/LAPACK):
#

USE_PETSC = no
PETSC_INC =
PETSC_DIR =
PETSC_PREFIX = /usr/local/petsc

ifeq ($(USE_PETSC), yes)
	ifeq ($(PETSC_DIR), )
		PETSC_DIR = $(PETSC_PREFIX)
	endif

	CFLAGS += -std=c11

	PETSC_INC = -I$(PETSC_DIR)/include -DUSE_PETSC
	PETSC_LIB = -Wl,-rpath,$(PETSC_DIR)/lib -L$(PETSC_DIR)/lib -lpetsc -lquadmath -ldl
endif

ifeq ($(USE_MPI), yes)
	PETSC_CONFIG = --prefix=$(PETSC_PREFIX) --with-cc=$(CC) --with-cxx=0 --with-fc=0 --with-shared-libraries=0 --with-x=0
else
	PETSC_CONFIG = --prefix=$(PETSC_PREFIX) --with-cc=$(CC) --with-cxx=0 --with-fc=0 --with-shared-libraries=0 --with-x=0 --with-mpi=0
endif

#
# SLEPc library (optionally used by mpi_lib module and requires PETSc):
#

USE_SLEPC = no
SLEPC_INC =
SLEPC_DIR =
SLEPC_PREFIX = /usr/local/slepc

ifeq ($(USE_SLEPC), yes)
	ifeq ($(SLEPC_DIR), )
		SLEPC_DIR = $(SLEPC_PREFIX)
	endif

	SLEPC_INC = -I$(SLEPC_DIR)/include -DUSE_SLEPC
	SLEPC_LIB = -Wl,-rpath,$(SLEPC_DIR)/lib -L$(SLEPC_DIR)/lib -lslepc
endif

#
# Extra macros, if any, in order to tune the building:
#

USE_MACRO = DUMMY_MACRO

#
# General rules:
#

all: modules drivers
modules: utils matrix nist johnson pes file math mpi_lib
drivers: d_basis pes_print b_print c_print m_print m_basis a+d_multipole a+d_cmatrix pec_print b_resize about

#
# Rules for modules:
# TODO: few modules are not checking if c_lib.h exists before compiling.
#

MODULES_DIR = modules

matrix: $(MODULES_DIR)/matrix.c $(MODULES_DIR)/matrix.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $(LINEAR_ALGEBRA_INC) -D$(USE_MACRO) -c $<
	@echo

nist: $(MODULES_DIR)/nist.c $(MODULES_DIR)/nist.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -c $<
	@echo

johnson: $(MODULES_DIR)/johnson.c $(MODULES_DIR)/johnson.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -c $<
	@echo

pes: $(MODULES_DIR)/pes.c $(MODULES_DIR)/pes.h $(MODULES_DIR)/math.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/cartesian.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $(PES_MACRO) -c $<
	@echo

file: $(MODULES_DIR)/file.c $(MODULES_DIR)/file.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -c $<
	@echo

math: $(MODULES_DIR)/math.c $(MODULES_DIR)/math.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -c $<
	@echo

mpi_lib: $(MODULES_DIR)/mpi_lib.c $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/c_lib.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $(MPI_INC) $(PETSC_INC) $(SLEPC_INC) -c $<
	@echo

blas_lib: $(MODULES_DIR)/blas_lib.c $(MODULES_DIR)/blas_lib.h $(MODULES_DIR)/gsl_lib.h $(MODULES_DIR)/c_lib.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -c $<
	@echo

#network: $(MODULES_DIR)/network.c $(MODULES_DIR)/network.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h
#	@echo "\033[31m$<\033[0m"
#	$(CC) $(CFLAGS) -c $<
#	@ech

#
# Rules for drivers:
#

d_basis: d_basis.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h $(MODULES_DIR)/math.h $(MODULES_DIR)/pes.h $(PES_OBJECT) nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o matrix.o file.o math.o pes.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

about: about.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h $(MODULES_DIR)/math.h $(MODULES_DIR)/nist.h $(MODULES_DIR)/pes.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out matrix.o mpi_lib.o pes.o file.o nist.o math.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

pes_print: pes_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

b_print: b_print.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h utils.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o matrix.o file.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

c_print: c_print.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o matrix.o file.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

m_basis_std: m_basis.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o mpi_lib.o matrix.o file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

m_basis: m_basis.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o mpi_lib.o matrix.o file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(SLEPC_LIB) $(PETSC_LIB) $(LINEAR_ALGEBRA_LIB)
	@echo

network: network.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h matrix.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out matrix.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

test_suit: test_suit.c $(MODULES_DIR)/matrix.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out matrix.o file.o math.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

a+d_multipole: a+d_multipole.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h $(PES_OBJECT)
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out file.o pes.o nist.o math.o mpi_lib.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

m_print: m_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

utils: utils.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -c
	@echo

a+d_cmatrix: a+d_cmatrix.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o mpi_lib.o matrix.o math.o file.o pes.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

pec_print: pec_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

b_resize: b_resize.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o matrix.o file.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

a+t_multipole: a+t_multipole.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h $(PES_OBJECT)
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out utils.o file.o pes.o nist.o math.o mpi_lib.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

numerov: numerov.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h $(PES_OBJECT) math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out utils.o matrix.o file.o pes.o nist.o math.o mpi_lib.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

#
# Rules for tools:
#

TOOLS_DIR = tools

3j: $(TOOLS_DIR)/3j.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

6j: $(TOOLS_DIR)/6j.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

9j: $(TOOLS_DIR)/9j.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

sphe_harmonics: $(TOOLS_DIR)/sphe_harmonics.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

percival_seaton: $(TOOLS_DIR)/percival_seaton.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

sparse_eigen: $(TOOLS_DIR)/sparse_eigen.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out mpi_lib.o matrix.o file.o $(LDFLAGS) $(SLEPC_LIB) $(PETSC_LIB) $(LINEAR_ALGEBRA_LIB)
	@echo

#
# Rules to build/install external libraries:
#

LIB_DIR = lib

gsl: $(LIB_DIR)/gsl-2.5.tar.gz $(GSL_DIR)
	tar -zxvf $<
	cd gsl-2.5/; ./configure CC=$(CC) --prefix=$(GSL_DIR); make; make install
	rm -rf gsl-2.5/

lapacke: $(LIB_DIR)/lapack-3.5.0.tar $(LAPACKE_DIR)/lib $(LAPACKE_DIR)/include
	tar -xvf $<
	cd lapack-3.5.0/; cp make.inc.example make.inc
	cd lapack-3.5.0/BLAS/SRC; make
	cd lapack-3.5.0; make; make lapackelib
	cd lapack-3.5.0; mv librefblas.a libblas.a; cp lib*.a $(LAPACKE_DIR)/lib/; cp lapacke/include/*.h $(LAPACKE_DIR)/include/
	rm -rf lapack-3.5.0/

magma: $(LIB_DIR)/magma-2.5.1-alpha1.tar.gz $(CUDA_DIR) $(MAGMA_DIR)
	tar -zxvf $<
	cp magma-2.5.1-alpha1/make.inc-examples/make.inc.mkl-$(CC) magma-2.5.1-alpha1/make.inc
	cd magma-2.5.1-alpha1/; export CUDADIR=$(CUDA_DIR); export GPU_TARGET="Kepler Maxwell Pascal"; make; make install prefix=$(MAGMA_DIR);
	rm -rf magma-2.5.1-alpha1/

arpack: $(LIB_DIR)/ARPACK.tar.xz $(ARPACKROOT)
	tar -xf $<
	cd ARPACK/; make lib home=$(PWD)/ARPACK PLAT=linux FC=$(FC) FFLAGS=-O3 MAKE=make SHELL=/bin/bash; mv libarpack_*.a $(ARPACKROOT)/
	rm -rf ARPACK

arpack-ng: $(LIB_DIR)/arpack-ng.tar.xz $(ARPACKROOT)
	tar -xf $<
	cd arpack-ng/; ./bootstrap
	cd arpack-ng/; export FFLAGS="-DMKL_ILP64 -I$(MKL_DIR)"; export FCLAGS="-DMKL_ILP64 -I$(MKL_DIR)"; export INTERFACE64=1; ./configure --with-blas=mkl_gf_ilp64 --with-lapack=mkl_gf_ilp64 --enable-icb F77=$(FC) FC=$(FC) CC=$(CC) CXX=g++ --prefix=$(ARPACKROOT)
	cd arpack-ng/; make; make install
	rm -rf arpack-ng

petsc: $(LIB_DIR)/petsc-3.13.4.tar.xz
	tar -xf $<
	mkdir -p $(PETSC_PREFIX)
	cd petsc-3.13.4/; ./configure $(PETSC_CONFIG)
	cd petsc-3.13.4/; make PETSC_DIR=$(PWD)/petsc-3.13.4 PETSC_ARCH=arch-linux2-c-debug all
	cd petsc-3.13.4/; make PETSC_DIR=$(PWD)/petsc-3.13.4 PETSC_ARCH=arch-linux2-c-debug install
	rm -rf petsc-3.13.4

slepc: $(LIB_DIR)/slepc-3.13.4.tar.gz $(PETSC_DIR)
	tar -zxvf $(LIB_DIR)/slepc-3.13.4.tar.gz
	mkdir -p $(SLEPC_PREFIX)
	cd slepc-3.13.4/; export PETSC_DIR=$(PETSC_DIR); export SLEPC_DIR=$(PWD)/slepc-3.13.4; ./configure --prefix=$(SLEPC_PREFIX)
	cd slepc-3.13.4/; make PETSC_DIR=$(PETSC_DIR) SLEPC_DIR=$(PWD)/slepc-3.13.4
	cd slepc-3.13.4/; make PETSC_DIR=$(PETSC_DIR) SLEPC_DIR=$(PWD)/slepc-3.13.4 install
	rm -rf slepc-3.13.4

clean:
	rm -f *.o *.out

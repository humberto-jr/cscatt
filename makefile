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

ifeq ($(CC), xlc)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
endif

#
# IBM xlc_r compiler:
#

ifeq ($(CC), xlc_r)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
endif

#
# IBM mpixlc compiler:
#

ifeq ($(CC), mpixlc)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
	USE_MPI = yes
endif

#
# IBM mpixlc_r compiler:
#

ifeq ($(CC), mpixlc_r)
	override CFLAGS = -std=c99 -q64 -qstrict -qsmp=omp -qthreaded -O5 -I$(GSL_DIR)/include
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
XLF_DIR = /opt/ibmcmp/xlf

ifeq ($(FC), xlf)
	LDFLAGS += -L$(XLF_DIR)/lib -lxlf90 -lxl -lxlfmath
endif

ifeq ($(FC), xlf90)
	LDFLAGS += -L$(XLF_DIR)/lib -lxlf90 -lxl -lxlfmath
endif

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
# User defined potential energy surface (required by the pes module):
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

	ifeq ($(CC), mpiicc)
		LINEAR_ALGEBRA_LIB = -parallel -Wl,--start-group $(MKL_DIR)/lib/intel64/libmkl_intel_ilp64.a $(MKL_DIR)/lib/intel64/libmkl_intel_thread.a $(MKL_DIR)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
	endif

	ifeq ($(CC), gcc)
		LINEAR_ALGEBRA_LIB = -Wl,--start-group $(MKL_DIR)/lib/intel64/libmkl_intel_ilp64.a $(MKL_DIR)/lib/intel64/libmkl_gnu_thread.a $(MKL_DIR)/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
	endif

	ifeq ($(CC), mpicc)
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
	LINEAR_ALGEBRA_LIB = -L$(LAPACKE_DIR)/lib -L$(LAPACKE_DIR)/lib64 -llapacke -llapack -lblas -lm
endif

BLAS_CONFIG = FORTRAN=$(FC)
LAPACK_CONFIG = FORTRAN=$(FC) LOADER=$(FC) CC=$(CC) NOOPT=-O0 TIMER=NONE

ifeq ($(FC), gfortran)
	BLAS_CONFIG += OPTS=-O3
	LAPACK_CONFIG += CFLAGS=-O3 OPTS=-O3
endif

ifeq ($(FC), ifort)
	BLAS_CONFIG += OPTS=-O3
	LAPACK_CONFIG += CFLAGS=-O3 OPTS=-O3
endif

ifeq ($(FC), xlf)
	BLAS_CONFIG += OPTS="-qstrict -O5"
	LAPACK_CONFIG += CFLAGS="-qstrict -O5" OPTS="-qstrict -O5"
endif

ifeq ($(FC), xlf90)
	BLAS_CONFIG += OPTS="-qstrict -O5"
	LAPACK_CONFIG += CFLAGS="-qstrict -O5" OPTS="-qstrict -O5"
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
# CUDA library
#

USE_CUDA = no
CUDA_DIR = /usr/local/cuda

ifeq ($(USE_CUDA), yes)
	LINEAR_ALGEBRA_INC = -I$(CUDA_DIR)/include
	LINEAR_ALGEBRA_LIB = -L$(CUDA_DIR)/lib64
endif

#
# MAGMA library (requires CUDA):
#

MAGMA_DIR = /usr/local/magma

ifeq ($(LINEAR_ALGEBRA), MAGMA)
	CFLAGS += -DUSE_MAGMA
	LINEAR_ALGEBRA_INC += -I$(MAGMA_DIR)/include -DADD_
	LINEAR_ALGEBRA_LIB += -L$(MAGMA_DIR)/lib -lmagma -lm
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

PETSC_CONFIG = --prefix=$(PETSC_PREFIX) --with-cc=$(CC) --with-cxx=0 --with-fc=0 --with-shared-libraries=0 --with-x=0

ifeq ($(USE_MPI), no)
	PETSC_CONFIG += --with-mpi=0
endif

ifeq ($(LINEAR_ALGEBRA), LAPACKE)
	PETSC_CONFIG += --with-blas-lib=$(LAPACKE_DIR)
	PETSC_CONFIG += --with-lapack-lib=$(LAPACKE_DIR)
endif

ifeq ($(USE_CUDA), yes)
	PETSC_CONFIG += --with-cuda --with-cuda-dir=$(CUDA_DIR)
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
modules: utils matrix nist johnson pes file math mpi_lib fgh spline
drivers: d_basis pes_print b_print c_print m_print a+d_basis a+d_multipole a+d_cmatrix pec_print b_resize about

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

fgh: $(MODULES_DIR)/fgh.c $(MODULES_DIR)/fgh.h $(MODULES_DIR)/file.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/globals.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -c $<
	@echo

spline: $(MODULES_DIR)/spline.c $(MODULES_DIR)/spline.h $(MODULES_DIR)/gsl_lib.h $(MODULES_DIR)/globals.h
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

d_basis: d_basis.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h $(MODULES_DIR)/math.h $(MODULES_DIR)/fgh.h $(MODULES_DIR)/pes.h $(PES_OBJECT) nist.o mpi_lib.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out matrix.o file.o math.o fgh.o pes.o nist.o mpi_lib.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

about: about.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h $(MODULES_DIR)/math.h $(MODULES_DIR)/nist.h $(MODULES_DIR)/pes.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out matrix.o mpi_lib.o pes.o file.o nist.o math.o $(PES_OBJECT) $(LDFLAGS) $(SLEPC_LIB) $(PETSC_LIB) $(LINEAR_ALGEBRA_LIB)
	@echo

pes_print: pes_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

b_print: b_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h $(MODULES_DIR)/fgh.h mpi_lib.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out matrix.o file.o fgh.o mpi_lib.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

c_print: c_print.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o matrix.o file.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

m_basis_std: m_basis.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o mpi_lib.o matrix.o file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

a+d_basis: a+d_basis.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
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

m_print: m_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h nist.o math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o pes.o nist.o math.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

utils: utils.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -c
	@echo

# TODO: abcd_pes_print driver is a temporary driver.
abcd_pes_print: abcd_pes_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/fgh.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

a+d_cmatrix: a+d_cmatrix.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out mpi_lib.o matrix.o math.o file.o fgh.o pes.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

pec_print: pec_print.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out file.o pes.o math.o nist.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

b_resize: b_resize.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) -D$(USE_MACRO) $< -o $@.out utils.o matrix.o file.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

a+t_multipole: a+t_multipole.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/math.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h $(PES_OBJECT)
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out file.o pes.o nist.o math.o mpi_lib.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

numerov: numerov.c utils.h $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h $(MODULES_DIR)/file.h $(MODULES_DIR)/pes.h $(PES_OBJECT) math.o nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out utils.o matrix.o file.o pes.o nist.o math.o mpi_lib.o $(PES_OBJECT) $(LDFLAGS) $(LINEAR_ALGEBRA_LIB) $(FORT_LIB)
	@echo

#
# Rules for debug:
#

setup:
	@echo "                CC = $(CC)"
	@echo "                FC = $(FC)"
	@echo "            CFLAGS = $(CFLAGS)"
	@echo "           LDFLAGS = $(LDFLAGS)"
	@echo "           USE_MPI = $(USE_MPI)"
	@echo "          USE_CUDA = $(USE_CUDA)"
	@echo "         USE_PETSC = $(USE_PETSC)"
	@echo "         USE_SLEPC = $(USE_SLEPC)"
	@echo "         PETSC_DIR = $(PETSC_DIR)"
	@echo "         PETSC_INC = $(PETSC_INC)"
	@echo "         PETSC_LIB = $(PETSC_LIB)"
	@echo "         SLEPC_DIR = $(SLEPC_DIR)"
	@echo "         SLEPC_INC = $(SLEPC_INC)"
	@echo "         SLEPC_LIB = $(SLEPC_LIB)"
	@echo "          PES_NAME = $(PES_NAME)"
	@echo "        PES_OBJECT = $(PES_OBJECT)"
	@echo "    LINEAR_ALGEBRA = $(LINEAR_ALGEBRA)"
	@echo "LINEAR_ALGEBRA_INC = $(LINEAR_ALGEBRA_INC)"
	@echo "LINEAR_ALGEBRA_LIB = $(LINEAR_ALGEBRA_LIB)"

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

dgemm_timer: $(TOOLS_DIR)/dgemm_timer.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out matrix.o $(LDFLAGS) $(LINEAR_ALGEBRA_LIB)
	@echo

gaunt: $(TOOLS_DIR)/gaunt.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

mass: $(TOOLS_DIR)/mass.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/nist.h nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out nist.o $(LDFLAGS)
	@echo

simpson_timer: $(TOOLS_DIR)/simpson_timer.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

gauss_legendre_timer: $(TOOLS_DIR)/gauss_legendre_timer.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

legendre_poly: $(TOOLS_DIR)/legendre_poly.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/math.h math.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out math.o $(LDFLAGS)
	@echo

mpi_tester: $(TOOLS_DIR)/mpi_tester.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/mpi_lib.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out mpi_lib.o $(LDFLAGS)
	@echo

multipole_integrand: $(TOOLS_DIR)/multipole_integrand.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/file.h $(MODULES_DIR)/math.h $(MODULES_DIR)/pes.h $(PES_OBJECT) nist.o
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out file.o pes.o nist.o math.o $(PES_OBJECT) $(LDFLAGS) $(FORT_LIB)
	@echo

simpson: $(TOOLS_DIR)/simpson.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/file.h $(MODULES_DIR)/math.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out matrix.o file.o math.o $(LDFLAGS)
	@echo

file_spline: $(TOOLS_DIR)/file_spline.c $(MODULES_DIR)/globals.h $(MODULES_DIR)/matrix.h $(MODULES_DIR)/spline.h $(MODULES_DIR)/file.h
	@echo "\033[31m$<\033[0m"
	$(CC) $(CFLAGS) $< -o $@.out matrix.o spline.o file.o $(LDFLAGS)
	@echo

#
# Rules to build/install external libraries:
#

LIB_DIR = lib

gsl: $(LIB_DIR)/gsl-2.5.tar.gz $(GSL_DIR)
	tar -zxvf $<
	cd gsl-2.5/; ./configure CC=$(CC) --prefix=$(GSL_DIR); make; make install
	rm -rf gsl-2.5/

lapacke: $(LIB_DIR)/lapack-3.5.0.tar
	tar -xvf $<
	mkdir -p $(LAPACKE_DIR)
	mkdir -p $(LAPACKE_DIR)/lib
	mkdir -p $(LAPACKE_DIR)/include
	cd lapack-3.5.0/; cp make.inc.example make.inc
	cd lapack-3.5.0/BLAS/SRC; make $(BLAS_CONFIG)
	cd lapack-3.5.0; make $(LAPACK_CONFIG); make lapackelib
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

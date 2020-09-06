/******************************************************************************

 About
 -----

 This module is a wrapper to MPI functions that encapsulates most of the actual
 calls of MPI routines. The wrapper is thread-safe and the OpenMP master thread
 is responsible for all operations. A macro USE_MPI is used to switch on/off
 the library usage. If off, all functions from this module becomes a set of
 valid compilable dummy calls.

******************************************************************************/

#include "mpi_lib.h"
#include "globals.h"

#if defined(USE_MPI)
	#include "mpi.h"
#endif

#if defined(USE_PETSC)
	#include <petscvec.h>
	#include <petscmat.h>
#endif

#if defined(USE_PETSC) && defined(USE_SLEPC)
	#include "slepceps.h"
#endif

FILE *stream = NULL;

static int rank = 0, size = 1, level = 0, chunk_size = 1, tasks = 1, extra_tasks = 0;
static int *index_min = NULL, *index_max = NULL;

/******************************************************************************

 Type mpi_matrix: defines a matrix object, often very large and sparse, whose
 non-zero elements are allocated among all MPI processes available. It needs
 the PETSc library by defining the USE_PETSC macro during compilation.

******************************************************************************/

struct mpi_matrix
{
	#if defined(USE_PETSC)
		Mat data;
	#endif

	#if defined(USE_SLEPC)
		EPS solver;
	#endif

	int max_row, max_col, first, last;
};

/******************************************************************************

 Type mpi_vector: defines a vector object whose elements are allocated among
 all MPI processes available. It needs the PETSc library by defining the
 USE_PETSC macro during compilation.

******************************************************************************/

struct mpi_vector
{
	#if defined(USE_PETSC)
		Vec data;
	#endif

	int length, first, last;
};

/******************************************************************************

 Macro CHECK_PETSC_ERROR(): checks the error code of PETSc calls and writes an
 error message in the C stderr with the name of the failed function, if the
 code corresponds to an error. When stop = true, the exection is terminated.

******************************************************************************/

#define CHECK_PETSC_ERROR(name, code, stop)                                  \
{                                                                            \
  if (code != 0)                                                             \
  {                                                                          \
    PRINT_ERROR("rank %d, %s failed with error code %d\n", rank, name, code) \
    if (stop) exit(EXIT_FAILURE);                                            \
  }                                                                          \
}

/******************************************************************************

 Macro ASSERT_ROW_INDEX(): check if the p-th element is within the row bounds
 of a matrix pointed by a given pointer.

******************************************************************************/

#define ASSERT_ROW_INDEX(pointer, p) \
{                                    \
  ASSERT((p) > -1)                   \
  ASSERT((p) < pointer->max_row)     \
}

/******************************************************************************

 Macro ASSERT_COL_INDEX(): check if the q-th element is within the column bounds
 of a matrix pointed by a given pointer.

******************************************************************************/

#define ASSERT_COL_INDEX(pointer, q) \
{                                    \
  ASSERT((q) > -1)                   \
  ASSERT((q) < pointer->max_col)     \
}

/******************************************************************************

 Macro ASSERT_RANK(): check if the n-th process is among those in the MPI
 communicator.

******************************************************************************/

#define ASSERT_RANK(n) \
{                      \
  ASSERT((n) > -1)     \
  ASSERT((n) < size)   \
}

/******************************************************************************

 Function mpi_type(): an auxiliary routine that returns the MPI datatype for a
 given C datatype (int, char, float and double).

 NOTE: only few types implemented and should be completed as needed.

******************************************************************************/

#if defined(USE_MPI)
	static MPI_Datatype mpi_type(const c_type c)
	{
		switch (c)
		{
			case type_int:    return MPI_INT;
			case type_char:   return MPI_CHAR;
			case type_float:  return MPI_FLOAT;
			case type_double: return MPI_DOUBLE;

			default:
				PRINT_ERROR("invalid C type %d\n", (int) c)
				exit(EXIT_FAILURE);
		}
	}
#endif

/******************************************************************************

 Function mpi_sizeof(): an auxiliary routine that returns the output of sizeof
 for a given C type.

 NOTE: only few types implemented and should be completed as needed.

******************************************************************************/

static size_t mpi_sizeof(const c_type c)
{
	switch (c)
	{
		case type_int:    return sizeof(int);
		case type_char:   return sizeof(char);
		case type_float:  return sizeof(float);
		case type_double: return sizeof(double);

		default:
			PRINT_ERROR("invalid C type %d\n", (int) c)
			exit(EXIT_FAILURE);
	}
}

/******************************************************************************

 Function mpi_init(): initializes the MPI library and should be invoked before
 any MPI call. It also initializes PETSc and/or SLEPc libraries if USE_PETSC
 and/or USE_SLEPC macros are defined during compilation.

******************************************************************************/

void mpi_init(int argc, char *argv[])
{
	ASSERT(argc > 0)
	ASSERT(argv != NULL)

	#if defined(USE_MPI)
	{
		MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &level);

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
	}
	#endif

	#if defined(USE_PETSC)
	{
		const int info = PetscInitialize(&argc, &argv, NULL, NULL);
		CHECK_PETSC_ERROR("PetscInitialize()", info, true)
	}
	#endif

	#if defined(USE_SLEPC)
	{
		const int info = SlepcInitialize(&argc, &argv, NULL, NULL);
		CHECK_PETSC_ERROR("SlepcInitialize()", info, true)
	}
	#endif

	return;
}

/******************************************************************************

 Function mpi_end(): finalizes the use of MPI and no calls to MPI functions
 should be done after. It also finalizes PETSc and/or SLEPc libraries if
 USE_PETSC and/or USE_SLEPC macros are defined during compilation.

******************************************************************************/

void mpi_end()
{
	#if defined(USE_MPI)
	{
		#pragma omp master
		{
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
		}
	}
	#endif

	#if defined(USE_SLEPC)
	{
		const int info = SlepcFinalize();
		CHECK_PETSC_ERROR("SlepcFinalize()", info, false)
	}
	#endif

	#if defined(USE_PETSC)
	{
		const int info = PetscFinalize();
		CHECK_PETSC_ERROR("PetscFinalize()", info, false)
	}
	#endif

	if (index_min != NULL) free(index_min);
	if (index_max != NULL) free(index_max);

	return;
}

/******************************************************************************

 Function mpi_rank(): returns the ID of the current process.

******************************************************************************/

int mpi_rank()
{
	return rank;
}

/******************************************************************************

 Function mpi_comm_size(): returns the maximum number of processes in
 communication.

******************************************************************************/

int mpi_comm_size()
{
	return size;
}

/******************************************************************************

 Function mpi_thread_level(): returns the thread level used.

******************************************************************************/

int mpi_thread_level()
{
	return level;
}

/******************************************************************************

 Function mpi_barrier(): blocks until all processes and all its threads in the
 communicator have reached this routine.

******************************************************************************/

void mpi_barrier()
{
	#if defined(USE_MPI)
	{
		#pragma omp critical
		MPI_Barrier(MPI_COMM_WORLD);
	}
	#endif
}

/******************************************************************************

 Function mpi_set_tasks(): divides a number of tasks among all processes by
 setting a minimum and maximum task index later returned by mpi_first_task()
 and mpi_last_task().

******************************************************************************/

void mpi_set_tasks(const int max_task)
{
	ASSERT(max_task > 0)

	tasks = max_task;
	chunk_size = max_task/size;

	index_min = allocate(size, sizeof(int), false);
	index_max = allocate(size, sizeof(int), false);

	for (int cpu_rank = 0; cpu_rank < size; ++cpu_rank)
	{
		index_min[cpu_rank] = cpu_rank*chunk_size;
		index_max[cpu_rank] = index_min[cpu_rank] + chunk_size - 1;
	}

	extra_tasks = (max_task - 1) - index_max[size - 1];

	ASSERT(extra_tasks >= 0)
}

/******************************************************************************

 Function mpi_first_task(): returns the index of a first task for this process.

******************************************************************************/

int mpi_first_task()
{
	return index_min[rank];
}

/******************************************************************************

 Function mpi_first_task(): returns the index of a last task for this process.

******************************************************************************/

int mpi_last_task()
{
	return index_max[rank];
}

/******************************************************************************

 Function mpi_first_task(): returns the index of an extra task for this process
 or zero if no extra tasks are scheduled.

******************************************************************************/

int mpi_extra_task()
{
	const int index = (extra_tasks > 0? index_max[size - 1] + rank + 1 : 0);

	return (index < tasks? index : 0);
}

/******************************************************************************

 Function mpi_check(): returns true if there is a message from a given source
 process and tag.

******************************************************************************/

bool mpi_check(const int from)
{
	ASSERT_RANK(from)

	int message_sent = (int) false;

	#if defined(USE_MPI)
	{
		#pragma omp critical
		MPI_Iprobe(from, 666, MPI_COMM_WORLD, &message_sent, MPI_STATUS_IGNORE);
	}
	#endif

	return (bool) message_sent;
}

/******************************************************************************

 Function mpi_send(): sends n elements, up to n_max, from the current process
 to another.

******************************************************************************/

void mpi_send(const int to, const int n_max, const c_type c, void *data)
{
	ASSERT_RANK(to)
	ASSERT(n_max > 0)
	ASSERT(data != NULL)
	ASSERT(c != type_unknown)

	#if defined(USE_MPI)
	{
		#pragma omp critical
		{
			const MPI_Datatype type_name = mpi_type(c);

			MPI_Send(&n_max, 1, MPI_INT, to, 666, MPI_COMM_WORLD);
			MPI_Send(data, n_max, type_name, to, 667, MPI_COMM_WORLD);
		}
	}
	#endif
}

/******************************************************************************

 Function mpi_receive(): receives n elements, up to n_max, from other process.

******************************************************************************/

void mpi_receive(const int from, const int n_max, const c_type c, void *data)
{
	ASSERT_RANK(from)
	ASSERT(n_max > 0)
	ASSERT(data != NULL)
	ASSERT(c != type_unknown)

	#if defined(USE_MPI)
	{
		#pragma omp critical
		{
			const MPI_Datatype type_name = mpi_type(c);

			int m_max = 0;

			MPI_Recv(&m_max, 1,
			         MPI_INT, from, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			m_max = min(n_max, m_max);

			MPI_Recv(data, m_max,
			         type_name, from, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	#endif
}

/******************************************************************************

 Function mpi_print(): prints a formatted message from a given process (master
 thread) in the C stdout of process 0 (master thread).

******************************************************************************/

void mpi_printf(const char line[])
{
	ASSERT(line != NULL)

	#pragma omp critical
	{
		#if defined(USE_MPI)
			const int tag_a = 803 + thread_id();
			const int tag_b = 997 + thread_id();

			if (rank == 0)
			{
				printf("%s", line);

				for (int source = 1; source < size; ++source)
				{
					int message_sent = (int) false;

					MPI_Iprobe(source, tag_a,
					           MPI_COMM_WORLD, &message_sent, MPI_STATUS_IGNORE);

					if (message_sent)
					{
						int length = 0;

						MPI_Recv(&length, 1, MPI_INT,
						         source, tag_a, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						char *new_line = allocate(length, sizeof(char), false);

						MPI_Recv(new_line, length, MPI_CHAR,
						         source, tag_b, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						printf("%s", new_line);

						free(new_line);
					}
				}
			}
			else
			{
				const int length = strlen(line);

				MPI_Request info;

				MPI_Isend(&length, 1, MPI_INT, 0, tag_a, MPI_COMM_WORLD, &info);

				MPI_Isend(line, length, MPI_CHAR, 0, tag_b, MPI_COMM_WORLD, &info);
			}
		#else
			printf("%s", line);
		#endif
	}

	return;
}

/******************************************************************************

 Function mpi_set_stream(): to set a disk stream for readings and writings done
 by process 0 (master thread) on behalf of all other processes.

******************************************************************************/

void mpi_set_stream(FILE *new_stream)
{
	ASSERT(new_stream != NULL);
	stream = new_stream;
}

/******************************************************************************

 Function mpi_fwrite(): process 0 (master thread) writes in the disk stream
 provided to mpi_set_stream() an array of data with a given lenght from all
 other processes, in ascending order of ranks.

******************************************************************************/

void mpi_fwrite(const int length, const c_type c, const void *array)
{
	ASSERT(length > 0)
	ASSERT(array != NULL)
	ASSERT(stream != NULL)

	const int data_size = mpi_sizeof(c);

	#pragma omp critical
	{
		#if defined(USE_MPI)
			const MPI_Datatype type_name = mpi_type(c);

			if (rank == 0)
			{
				if(fwrite(array, length, data_size, stream) != (size_t) length)
				{
					PRINT_ERROR("unable to write %d elements from process 0\n", length)
				}

				for (int source = 1; source < size; ++source)
				{
					int message_sent = (int) false;

					MPI_Iprobe(source, 731,
					           MPI_COMM_WORLD, &message_sent, MPI_STATUS_IGNORE);

					if (message_sent)
					{
						int length = 0;

						MPI_Recv(&length, 1, MPI_INT,
						         source, 731, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						void *new_array = allocate(length, data_size, false);

						MPI_Recv(new_array, length, type_name,
						         source, 732, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						if(fwrite(array, length, data_size, stream) != (size_t) length)
						{
							PRINT_ERROR("unable to write %d elements from process %d\n", length, source)
						}

						free(new_array);
					}
				}
			}
			else
			{
				MPI_Request info;

				MPI_Isend(&length, 1, MPI_INT, 0, 731, MPI_COMM_WORLD, &info);

				MPI_Isend(array, length, type_name, 0, 732, MPI_COMM_WORLD, &info);
			}
		#else
			if(fwrite(array, length, data_size, stream) != (size_t) length)
			{
				PRINT_ERROR("unable to write %d elements from process 0\n", length)
			}
		#endif
	}

	return;
}

/******************************************************************************

 Function mpi_matrix_alloc(): allocate resources for a matrix of shape max_row-
 by-max_col, where chunks of rows are stored in each MPI process available. The
 approximated number of non-zero elements, in the diagonal blocks, given by
 non_zeros[0], and off-diagonal blocks, non_zeros[1], is provided to
 optimize memory allocation.

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

mpi_matrix *mpi_matrix_alloc(const int max_row,
                             const int max_col, const int non_zeros[])
{
	ASSERT(max_row > 0)
	ASSERT(max_col > 0)
	ASSERT(non_zeros[0] >= 0)
	ASSERT(non_zeros[1] >= 0)

	#if defined(USE_PETSC)
	{
		int info;

		mpi_matrix *pointer = allocate(1, sizeof(struct mpi_matrix), true);

		pointer->max_row = max_row;
		pointer->max_col = max_col;

/*		info = MatCreate(MPI_COMM_WORLD, &pointer->data);

		CHECK_PETSC_ERROR("MatCreate()", info, true)

		info = MatSetType(pointer->data, MATAIJ);

		CHECK_PETSC_ERROR("MatSetType()", info, true)

		info = MatSetSizes(pointer->data, PETSC_DECIDE, PETSC_DECIDE, max_row, max_col);

		CHECK_PETSC_ERROR("MatSetSizes()", info, true)

		info = MatSetUp(pointer->data);

		CHECK_PETSC_ERROR("MatSetUp()", info, true)

		info = MatMPIAIJSetPreallocation(pointer->data,
		                                 non_zeros[0], NULL, non_zeros[1], NULL);
*/
		if (mpi_comm_size() > 1)
		{
			info = MatCreateAIJ(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, max_row,
			                    max_col, non_zeros[0], NULL, non_zeros[1], NULL, &pointer->data);

			CHECK_PETSC_ERROR("MatCreateAIJ()", info, true)
		}
		else
		{
			info = MatCreateSeqAIJ(PETSC_COMM_SELF,
			                       max_row, max_col, non_zeros[0], NULL, &pointer->data);

			CHECK_PETSC_ERROR("MatCreateSeqAIJ()", info, true)
		}

		info = MatSetOption(pointer->data, MAT_NEW_NONZERO_ALLOCATION_ERR, false);

		CHECK_PETSC_ERROR("MatSetOption()", info, true)

		info = MatSetOption(pointer->data, MAT_USE_HASH_TABLE, true);

		CHECK_PETSC_ERROR("MatSetOption()", info, true)

		info = MatGetOwnershipRange(pointer->data, &pointer->first, &pointer->last);

		CHECK_PETSC_ERROR("MatGetOwnershipRange()", info, true)

		return pointer;
	}
	#endif

	return NULL;
}

/******************************************************************************

 Function mpi_matrix_free(): release resources allocated by mpi_matrix_alloc().

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_matrix_free(mpi_matrix *m)
{
	ASSERT(m != NULL)

	#if defined(USE_PETSC)
	{
		const int info = MatDestroy(&m->data);
		CHECK_PETSC_ERROR("MatDestroy()", info, true)
		free(m);
	}
	#endif
}

/******************************************************************************

 Function mpi_matrix_set(): stores all non-zero elements that will be used by
 mpi_matrix_build() to construct a sparse matrix among all MPI processes.

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_matrix_set(mpi_matrix *m, const int p, const int q, const double x)
{
	ASSERT(m != NULL)
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)

	#if defined(USE_PETSC)
	{
		if ((p >= m->first) && (p < m->last))
		{
			const int info = MatSetValue(m->data, p, q, x, INSERT_VALUES);
			CHECK_PETSC_ERROR("MatSetValue()", info, true)
		}
	}
	#else
	{
		ASSERT(x == x)
		PRINT_ERROR("%s\n", "the PETSc library is required")
	}
	#endif
}

/******************************************************************************

 Function mpi_matrix_build(): builds the matrix among all MPI processes after
 its elements have been cached by all needed calls of mpi_matrix_set().

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_matrix_build(mpi_matrix *m)
{
	ASSERT(m != NULL)

	#if defined(USE_PETSC)
	{
		int info;

		info = MatAssemblyBegin(m->data, MAT_FINAL_ASSEMBLY);
		info = MatAssemblyEnd(m->data, MAT_FINAL_ASSEMBLY);

		CHECK_PETSC_ERROR("MatAssemblyEnd()", info, true)
	}
	#endif
}

/******************************************************************************

 Function mpi_vector_alloc(): allocate resources for a vector of a given length
 where chunks of elements are stored in each MPI process available.

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

mpi_vector *mpi_vector_alloc(const int length)
{
	ASSERT(length > 0)

	#if defined(USE_PETSC)
	{
		int info;

		mpi_vector *pointer = allocate(1, sizeof(mpi_vector), true);

		info = VecCreate(MPI_COMM_WORLD, &pointer->data);

		CHECK_PETSC_ERROR("VecCreate()", info, true)

		pointer->length = length;

		info = VecSetSizes(pointer->data, PETSC_DECIDE, pointer->length);

		CHECK_PETSC_ERROR("VecSetSizes()", info, true)

		info = VecGetOwnershipRange(pointer->data, &pointer->first, &pointer->last);

		CHECK_PETSC_ERROR("VecGetOwnershipRange()", info, true)

		return pointer;
	}
	#endif

	return NULL;
}

/******************************************************************************

 Function mpi_vector_free(): release resources allocated by mpi_vector_alloc().

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_vector_free(mpi_vector *v)
{
	ASSERT(v != NULL)

	#if defined(USE_PETSC)
	{
		const int info = VecDestroy(&v->data);
		CHECK_PETSC_ERROR("VecDestroy()", info, true)
		free(v);
	}
	#endif
}

/******************************************************************************

 Function mpi_vector_build(): builds the vector among all MPI processes after
 its elements have been cached by all needed calls of mpi_vector_set().

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

 TODO: the actual routine mpi_vector_set() is not yet implemented.

******************************************************************************/

void mpi_vector_build(mpi_vector *v)
{
	ASSERT(v != NULL)

	#if defined(USE_PETSC)
	{
		int info;

		info = VecAssemblyBegin(v->data);
		info = VecAssemblyEnd(v->data);

		CHECK_PETSC_ERROR("VecAssemblyEnd()", info, true)
	}
	#endif
}

/******************************************************************************

 Function mpi_matrix_sparse_eigen(): computes n eigenvalues and eigenvectors of
 a matrix m, where the upper part of the spectrum is resolved if up = true, and
 the lowest part is computed otherwise. On entry, m is expected hermitian. The
 method used is Krylov-Schur, a variation of Arnoldi with a very effective
 restarting technique. In the case of symmetric problems, this is
 equivalent to the thick-restart Lanczos method.

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

 NOTE: it requires the SLEPc library by defining the USE_SLEPC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_matrix_sparse_eigen(mpi_matrix *m, const int n,
                             const int max_step, const double tol, const bool up)
{
	ASSERT(m != NULL)
	ASSERT(n > 0)

	#if defined(USE_PETSC) && defined(USE_SLEPC)
	{
		int info = EPSCreate(MPI_COMM_WORLD, &m->solver);

		CHECK_PETSC_ERROR("EPSCreate()", info, true)

		info = EPSSetOperators(m->solver, m->data, NULL);

		CHECK_PETSC_ERROR("EPSSetOperators()", info, true)

		info = EPSSetDimensions(m->solver, n, 2*n + 10, n + 2*n + 10);

		CHECK_PETSC_ERROR("EPSSetDimensions()", info, true)

		info = EPSSetProblemType(m->solver, EPS_HEP);

		CHECK_PETSC_ERROR("EPSSetProblemType()", info, true)

		const EPSWhich job = (up? EPS_LARGEST_MAGNITUDE : EPS_SMALLEST_MAGNITUDE);

		info = EPSSetWhichEigenpairs(m->solver, job);

		CHECK_PETSC_ERROR("EPSSetWhichEigenpairs()", info, true)

		info = EPSSetType(m->solver, EPSKRYLOVSCHUR);

		CHECK_PETSC_ERROR("EPSSetType()", info, true)

		info = EPSSetTolerances(m->solver, tol, max_step);

		CHECK_PETSC_ERROR("EPSSetTolerances()", info, true)

		info = EPSSolve(m->solver);

		CHECK_PETSC_ERROR("EPSSolve()", info, true)
	}
	#else
	{
		ASSERT(up == up)
		PRINT_ERROR("%s\n", "both PETSc and SLEPc libraries are required")
	}
	#endif
}

/******************************************************************************

 Function mpi_matrix_eigenpair(): returns the n-th eigenvalue and eigenvector
 of a matrix m after a call to mpi_matrix_sparse_eigen().

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

 NOTE: it requires the SLEPc library by defining the USE_SLEPC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

mpi_vector *mpi_matrix_eigenpair(const mpi_matrix *m,
                                 const int n, double *eigenval)
{
	ASSERT(m != NULL)
	ASSERT_COL_INDEX(m, n)

	#if defined(USE_PETSC) && defined(USE_SLEPC)
	{
		int n_max, info;

		info = EPSGetConverged(m->solver, &n_max);

		CHECK_PETSC_ERROR("EPSGetConverged()", info, true)

		if (n >= n_max)
		{
			PRINT_ERROR("n = %d is out of the %d converged solutions\n", n, n_max)

			*eigenval = 0.0;
			return NULL;
		}

		mpi_vector *eigenvec = allocate(1, sizeof(struct mpi_vector), true);

		info = MatCreateVecs(m->data, NULL, &eigenvec->data);

		CHECK_PETSC_ERROR("MatCreateVecs()", info, true)

		info = VecGetSize(eigenvec->data, &eigenvec->length);

		CHECK_PETSC_ERROR("VecGetSize()", info, true)

		info = VecGetOwnershipRange(eigenvec->data, &eigenvec->first, &eigenvec->last);

		CHECK_PETSC_ERROR("VecGetOwnershipRange()", info, true)

		info = EPSGetEigenpair(m->solver, n, eigenval, NULL, eigenvec->data, NULL);

		CHECK_PETSC_ERROR("EPSGetEigenpair()", info, true)

		return eigenvec;
	}
	#else
	{
		PRINT_ERROR("%s\n", "both PETSc and SLEPc libraries are required")
	}
	#endif

	*eigenval = 0.0;
	return NULL;
}

/******************************************************************************

 Function mpi_vector_write(): writes n elements, with n = [n_min, n_max), from
 vector v to an output file using binary format. On entry, the file stream is
 written by the MPI process zero and may be set NULL for other processes.

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_vector_write(const mpi_vector *v,
                      const int n_min, const int n_max, FILE *stream)
{
	ASSERT(v != NULL)
	ASSERT(n_min > -1)
	ASSERT(n_max > n_min)

	if (mpi_rank() == 0)
	{
		ASSERT(stream != NULL)
	}

	#if defined(USE_PETSC)
	{
		mpi_barrier();

		int info;
		double value;

		if (mpi_rank() == 0)
		{
			if (n_min < v->last)
			{
				for (int n = max(n_min, v->first); n < min(n_max, v->last); ++n)
				{
					info = VecGetValues(v->data, 1, &n, &value);

					CHECK_PETSC_ERROR("VecGetValues()", info, true)

					fwrite(&value, sizeof(double), 1, stream);
				}

				if (mpi_comm_size() == 1) return;
			}
			else
			{
/*
 *				NOTE: when this if-block is false a delay is used to make sure other
 *				processes had enough time to send at least the first data received
 *				by the for-block below. Ten milliseconds are used.
 */
				sleep(0.010);
			}

			for (int from = 1; from < mpi_comm_size(); ++from)
			{
				while (mpi_check(from))
				{
					mpi_receive(from, 1, type_double, &value);
					fwrite(&value, sizeof(double), 1, stream);
				}
			}
		}
		else
		{
			if (n_min >= v->last) return;

			for (int n = max(n_min, v->first); n < min(n_max, v->last); ++n)
			{
				info = VecGetValues(v->data, 1, &n, &value);

				CHECK_PETSC_ERROR("VecGetValues()", info, true)

				mpi_send(0, 1, type_double, &value);
			}
		}
	}
	#endif
}

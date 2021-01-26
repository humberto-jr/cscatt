/******************************************************************************

 About
 -----

 This module is a wrapper to MPI functions that encapsulates most of the actual
 calls of MPI routines. The wrapper is thread-safe and only one OpenMP thread
 is executed at the time. Macro USE_MPI is used during compilation to switch
 on/off the MPI usage. If off, all functions from this module become a set of
 valid compilable dummy calls.

******************************************************************************/

#include "mpi_lib.h"
#include "globals.h"

#if defined(USE_MPI) && !defined(USE_PETSC)
	#include "mpi.h"
#endif

#if defined(USE_PETSC)
	#include <petscvec.h>
	#include <petscmat.h>
#endif

#if defined(USE_PETSC) && defined(USE_SLEPC)
	#include "slepceps.h"
#endif

static int this_rank = 0, comm_size = 1, thread_level = 0;

static int chunk_size = 1, tasks = 1, extra_tasks = 0, last_rank_index = 0;

/******************************************************************************

 Type mpi_matrix: defines a matrix object, often very large and sparse, whose
 non-zero elements are stored among all MPI processes available. It needs the
 PETSc and/or SLEPc libraries by defining the USE_PETSC and USE_SLEPC macros
 during compilation.

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

 Type mpi_vector: defines a vector object, often very large and sparse, whose
 non-zero elements are stored among all MPI processes available. It needs the
 PETSc and/or SLEPc libraries by defining the USE_PETSC and USE_SLEPC macros
 during compilation.

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
 error message in the C stderr reporting the problem if the code corresponds to
 an error. When stop = true, the execution is terminated.

******************************************************************************/

#define CHECK_PETSC_ERROR(name, code, stop)                                       \
{                                                                                 \
  if (code != 0)                                                                  \
  {                                                                               \
    PRINT_ERROR("rank %d, %s failed with error code %d\n", this_rank, name, code) \
    if (stop) exit(EXIT_FAILURE);                                                 \
  }                                                                               \
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

 Macro ASSERT_COL_INDEX(): check if the q-th element is within the column
 bounds of a matrix pointed by a given pointer.

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

#define ASSERT_RANK(n)    \
{                         \
  ASSERT((n) > -1)        \
  ASSERT((n) < comm_size) \
}

/******************************************************************************

 Function mpi_init(): initializes the MPI library and shall be invoked before
 any MPI call. It also initializes PETSc and/or SLEPc libraries if USE_PETSC
 and USE_SLEPC macros are defined during compilation.

******************************************************************************/

void mpi_init(int argc, char *argv[])
{
	ASSERT(argc > 0)
	ASSERT(argv != NULL)

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

	#if defined(USE_MPI)
	{
		#if !defined(USE_PETSC)
		{
			#pragma omp master
			MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_level);
		}
		#endif

		#pragma omp master
		MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

		#pragma omp master
		MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	}
	#endif
}

/******************************************************************************

 Function mpi_using_petsc(): return true if the macro USE_PETSC was used during
 compilation, false otherwise.

******************************************************************************/

bool mpi_using_petsc()
{
	#if defined(USE_PETSC)
		return true;
	#else
		return false;
	#endif
}

/******************************************************************************

 Function mpi_using_slepc(): return true if the macro USE_SLEPC was used during
 compilation, false otherwise.

******************************************************************************/

bool mpi_using_slepc()
{
	#if defined(USE_SLEPC)
		return true;
	#else
		return false;
	#endif
}

/******************************************************************************

 Function mpi_end(): finalizes the use of MPI and calls to MPI functions shall
 not be done after. It also finalizes PETSc and/or SLEPc libraries if USE_PETSC
 and USE_SLEPC macros are defined during compilation.

******************************************************************************/

void mpi_end()
{
	#if defined(USE_MPI) && !defined(USE_PETSC)
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
}

/******************************************************************************

 Function mpi_rank(): returns the rank ID of the current MPI process.

******************************************************************************/

int mpi_rank()
{
	return this_rank;
}

/******************************************************************************

 Function mpi_comm_size(): returns the maximum number of MPI processes in
 communication.

******************************************************************************/

int mpi_comm_size()
{
	return comm_size;
}

/******************************************************************************

 Function mpi_thread_level(): returns the level of thread support used.

******************************************************************************/

int mpi_thread_level()
{
	return thread_level;
}

/******************************************************************************

 Function mpi_barrier(): blocks until all MPI processes in the communicator,
 one thread at the time, have reached this routine.

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

 Function mpi_set_tasks(): divides a given number of tasks among all MPI
 processes by setting a minimum and maximum task index that are later returned
 by mpi_first_task() and mpi_last_task(). Each process shall have a different
 set of indices depending on its own rank.

******************************************************************************/

void mpi_set_tasks(const int max_task)
{
	ASSERT(max_task > 0)

	tasks = max_task;
	chunk_size = max_task/comm_size;
	last_rank_index = (comm_size - 1)*chunk_size + (chunk_size - 1);
	extra_tasks = (max_task - 1) - last_rank_index;

	ASSERT(extra_tasks >= 0)
}

/******************************************************************************

 Function mpi_first_task(): returns the index of a first task for the process
 that calls it.

 NOTE: mpi_first_task() shall return 0 when called by the MPI process 0.

******************************************************************************/

int mpi_first_task()
{
	return this_rank*chunk_size;
}

/******************************************************************************

 Function mpi_first_task(): returns the index of a last task for the process
 that calls it.

 NOTE: mpi_first_task() shall return max_task from mpi_set_tasks() when called
 by the last MPI process.

******************************************************************************/

int mpi_last_task()
{
	return this_rank*chunk_size + (chunk_size - 1);
}

/******************************************************************************

 Function mpi_extra_task(): returns the index of an extra task for the process
 that calls it, if max_task from mpi_set_tasks() is not divisible by the
 number of MPI processes. It returns 0 if there is no extra task.

 NOTE: if extra tasks are available they are given to those lower-rank MPI
 processes starting from 0.

******************************************************************************/

int mpi_extra_task()
{
	const int index = (extra_tasks > 0? last_rank_index + this_rank + 1 : 0);

	return (index < tasks? index : 0);
}

/******************************************************************************

 Function mpi_inbox(): returns true if there is a message from a given MPI
 process.

******************************************************************************/

bool mpi_inbox(const int from)
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

 Function mpi_send(): sends n elements from the current MPI process to another.
 Where, type is one of: 'i' for int, 'c' for char, 'f' for float or 'd' for
 double.

******************************************************************************/

void mpi_send(const int to, const int n, const char type, void *data)
{
	ASSERT(n > 0)
	ASSERT_RANK(to)
	ASSERT(data != NULL)

	/* NOTE: to avoid 'unused parameter' warns during compilation of the dummy version. */
	ASSERT(type == type)

	#if defined(USE_MPI)
	{
		#pragma omp critical
		MPI_Send(&n, 1, MPI_INT, to, 666, MPI_COMM_WORLD);

		#pragma omp critical
		switch (type)
		{
			case 'i':
				MPI_Send(data, n, MPI_INT, to, 667, MPI_COMM_WORLD);
				break;

			case 'c':
				MPI_Send(data, n, MPI_CHAR, to, 667, MPI_COMM_WORLD);
				break;

			case 'f':
				MPI_Send(data, n, MPI_FLOAT, to, 667, MPI_COMM_WORLD);
				break;

			case 'd':
				MPI_Send(data, n, MPI_DOUBLE, to, 667, MPI_COMM_WORLD);
				break;

			default:
				PRINT_ERROR("invalid type %c\n", type)
				exit(EXIT_FAILURE);
		}
	}
	#endif
}

/******************************************************************************

 Function mpi_receive(): receives n elements sent by mpi_send() from another
 MPI process. Where, type is one of: 'i' for int, 'c' for char, 'f' for float
 or 'd' for double.

 NOTE: n from the sender can be different from that in the receiver, n', which
 receives only min(n, n').

******************************************************************************/

void mpi_receive(const int from, const int n, const char type, void *data)
{
	ASSERT(n > 0)
	ASSERT_RANK(from)
	ASSERT(data != NULL)

	/* NOTE: to avoid 'unused parameter' warns during compilation of the dummy version. */
	ASSERT(type == type)

	#if defined(USE_MPI)
	{
		int m = 0;

		#pragma omp critical
		MPI_Recv(&m, 1, MPI_INT, from, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		m = min(n, m);

		#pragma omp critical
		switch (type)
		{
			case 'i':
				MPI_Recv(data, m, MPI_INT, from, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				break;

			case 'c':
				MPI_Recv(data, m, MPI_CHAR, from, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				break;

			case 'f':
				MPI_Recv(data, m, MPI_FLOAT, from, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				break;

			case 'd':
				MPI_Recv(data, m, MPI_DOUBLE, from, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				break;

			default:
				PRINT_ERROR("invalid type %c\n", type)
				exit(EXIT_FAILURE);
		}
	}
	#endif
}

/******************************************************************************

 Function mpi_print(): prints a formatted message from a given process (master
 thread) in the C stdout of process 0 (master thread).

 TODO: re-checking stopped here.

******************************************************************************/

/*
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
*/

/******************************************************************************

 Function mpi_fwrite(): process 0 (master thread) writes in the disk stream
 provided to mpi_set_stream() an array of data with a given lenght from all
 other processes, in ascending order of ranks.

******************************************************************************/

/*
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
*/

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

/*
		info = MatCreate(MPI_COMM_WORLD, &pointer->data);

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

/*
		FIXME: PETSc fails when using hash table with MPI.
		info = MatSetOption(pointer->data, MAT_USE_HASH_TABLE, true);

		CHECK_PETSC_ERROR("MatSetOption()", info, true)
*/

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

 Function mpi_matrix_set(): stores each non-zero pq-element that will be used
 by mpi_matrix_build() to construct a sparse matrix among all MPI processes.

 NOTE: it requires the PETSc library by defining the USE_PETSC macro during
 compilation. Otherwise, this is a dummy function.

******************************************************************************/

void mpi_matrix_set(mpi_matrix *m, const int p, const int q, const double x)
{
	ASSERT(m != NULL)
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)

	/* NOTE: to avoid 'unused parameter' warns during compilation of the dummy version. */
	ASSERT(x == x)

	#if defined(USE_PETSC)
	{
		if ((p >= m->first) && (p < m->last))
		{
			const int info = MatSetValue(m->data, p, q, x, INSERT_VALUES);
			CHECK_PETSC_ERROR("MatSetValue()", info, true)
		}
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

int mpi_matrix_sparse_eigen(mpi_matrix *m, const int n,
                            const int max_step, const double tol, const bool up)
{
	ASSERT(n > 0)
	ASSERT(m != NULL)
	ASSERT(max_step > 0)

	/* NOTE: to avoid 'unused parameter' warns during compilation of the dummy version. */
	ASSERT(up == up)
	ASSERT(tol == tol)

	#if defined(USE_PETSC) && defined(USE_SLEPC)
	{
		int info = EPSCreate(MPI_COMM_WORLD, &m->solver);

		CHECK_PETSC_ERROR("EPSCreate()", info, true)

		info = EPSSetOperators(m->solver, m->data, NULL);

		CHECK_PETSC_ERROR("EPSSetOperators()", info, true)

		info = EPSSetDimensions(m->solver, n, 2*n + 10, PETSC_DECIDE);

		CHECK_PETSC_ERROR("EPSSetDimensions()", info, true)

		info = EPSSetProblemType(m->solver, EPS_HEP);

		CHECK_PETSC_ERROR("EPSSetProblemType()", info, true)

		const EPSWhich job = (up? EPS_LARGEST_REAL : EPS_SMALLEST_REAL);

		info = EPSSetWhichEigenpairs(m->solver, job);

		CHECK_PETSC_ERROR("EPSSetWhichEigenpairs()", info, true)

		info = EPSSetType(m->solver, EPSKRYLOVSCHUR);

		CHECK_PETSC_ERROR("EPSSetType()", info, true)

		info = EPSSetTolerances(m->solver, tol, max_step);

		CHECK_PETSC_ERROR("EPSSetTolerances()", info, true)

		info = EPSSolve(m->solver);

		CHECK_PETSC_ERROR("EPSSolve()", info, true)

		int n_max = 0;
		info = EPSGetConverged(m->solver, &n_max);

		CHECK_PETSC_ERROR("EPSGetConverged()", info, true)

		return n_max;
	}
	#endif

	return 0;
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
		mpi_vector *eigenvec = allocate(1, sizeof(struct mpi_vector), true);

		int info = MatCreateVecs(m->data, NULL, &eigenvec->data);

		CHECK_PETSC_ERROR("MatCreateVecs()", info, true)

		info = VecGetSize(eigenvec->data, &eigenvec->length);

		CHECK_PETSC_ERROR("VecGetSize()", info, true)

		info = VecGetOwnershipRange(eigenvec->data, &eigenvec->first, &eigenvec->last);

		CHECK_PETSC_ERROR("VecGetOwnershipRange()", info, true)

		info = EPSGetEigenpair(m->solver, n, eigenval, NULL, eigenvec->data, NULL);

		CHECK_PETSC_ERROR("EPSGetEigenpair()", info, true)

		return eigenvec;
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
				while (mpi_inbox(from))
				{
					mpi_receive(from, 1, 'd', &value);
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

				mpi_send(0, 1, 'd', &value);
			}
		}
	}
	#endif
}

/******************************************************************************

 Function mpi_about(): prints in a given output file the conditions in which
 the module was compiled.

******************************************************************************/

void mpi_about(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "# build date  = %s\n", __DATE__);
	fprintf(output, "# source code = %s\n", __FILE__);

	#if defined(USE_MPI)
		fprintf(output, "# using MPI   = yes\n");
	#else
		fprintf(output, "# using MPI   = no\n");
	#endif

	#if defined(USE_PETSC)
		fprintf(output, "# using PETSc = yes\n");
	#else
		fprintf(output, "# using PETSc = no\n");
	#endif

	#if defined(USE_SLEPC)
		fprintf(output, "# using SLEPc = yes\n");
	#else
		fprintf(output, "# using SLEPc = no\n");
	#endif
}

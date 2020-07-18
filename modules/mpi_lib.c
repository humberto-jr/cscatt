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

static int rank = 0, size = 1, level = 0;

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

 Function mpi_init(): initializes the MPI library and should be invoked before
 any MPI call.

******************************************************************************/

void mpi_init(int argc, char *argv[])
{
	ASSERT(argc > 1)
	ASSERT(argv != NULL)

	#if defined(USE_MPI)
		MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &level);

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
	#endif

	return;
}

/******************************************************************************

 Function mpi_end(): finalizes the use of MPI and no calls to MPI functions
 should be done after.

******************************************************************************/

void mpi_end()
{
	#if defined(USE_MPI)
		#pragma omp master
		{
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
		}
	#endif

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

 Function mpi_send(): sends n elements, up to n_max, of a given array, from one
 process to that same interval in the array of another process.

******************************************************************************/

void mpi_send(const int to,
              const int from,
              const int n_max,
              const c_type c,
              void *array)
{
	ASSERT(to > -1)
	ASSERT(from > -1)
	ASSERT(n_max > 0)
	ASSERT(array != NULL)
	ASSERT(c != type_unknown)

	if (to == from) return;

	#if defined(USE_MPI)
		#pragma omp master
		{
			const MPI_Datatype typename = mpi_type(c);

			if (from == rank)
			{
				MPI_Request info;

				MPI_Isend(&n_max, 1,
				          MPI_INT, to, 666, MPI_COMM_WORLD, &info);

				MPI_Isend(array, n_max,
				          typename, to, 667, MPI_COMM_WORLD, &info);
			}

			if (to == rank)
			{
				int m_max = 0;

				MPI_Recv(&m_max, 1,
				         MPI_INT, from, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				m_max = min(n_max, m_max);

				MPI_Recv(array, m_max,
				         typename, from, 667, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	#endif
}

/******************************************************************************

 Function mpi_print(): writes a short log from a given CPU in the C stdout of
 CPU 0. Notice, if OpenMP is used, only the thread 0 from each CPU sends data
 whereas only the thread 0 from CPU 0 receives.

******************************************************************************/

void mpi_printf(const char line[])
{
	ASSERT(line != NULL)

	#if defined(USE_MPI)
		#pragma omp master
		{
			if (rank == 0)
			{
				printf("%s", line);

				for (int from = 1; from < size; ++from)
				{
					int length = 0;

					MPI_Request info;

					MPI_Irecv(&length, 1, MPI_INT,
					          from, 998, MPI_COMM_WORLD, &info); // bug here

					if (length > 0)
					{
						char *new_line = allocate(length, sizeof(char), false);

						MPI_Recv(new_line, length, MPI_CHAR,
						         from, 999, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						printf("%s", new_line);

						free(new_line);
					}
				}
			}
			else
			{
				const int length = strlen(line);

				MPI_Request info;

				MPI_Isend(&length, 1, MPI_INT, 0, 998, MPI_COMM_WORLD, &info);

				MPI_Isend(line, length, MPI_CHAR, 0, 999, MPI_COMM_WORLD, &info);
			}
		}
	#else
		printf("%s", line);
	#endif

	return;
}

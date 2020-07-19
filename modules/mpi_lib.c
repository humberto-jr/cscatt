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

FILE *stream = NULL;
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
 any MPI call.

******************************************************************************/

void mpi_init(int argc, char *argv[])
{
	ASSERT(argc > 1)
	ASSERT(argv != NULL)

	#if defined(USE_MPI)
		MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &level);

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

void mpi_fwrite(const int length,
                const c_type c, const void *array)
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

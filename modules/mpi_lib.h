#if !defined(MPI_LIB_HEADER)
	#define MPI_LIB_HEADER
	#include "c_lib.h"

	void mpi_init(int argc, char *argv[]);

	void mpi_end();

	int mpi_rank();

	int mpi_comm_size();

	int mpi_thread_level();

	void mpi_send(const int to,
	              const int from,
	              const int n_max,
	              const c_type c,
	              void *array);

	void mpi_printf(const char line[]);

	void mpi_set_stream(FILE *new_stream);

	void mpi_fwrite(const int length, const c_type c, const void *array);

	/******************************************************************************

	 Macro MPI_PRINTF(): prints a formatted message from a given process (master
	 thread) in the C stdout of process 0 (master thread).

	 NOTE: a maximum of 1024 characters only.

	******************************************************************************/

	#define MPI_PRINTF(format, ...)         \
	{                                       \
	  char line[1024];                      \
	  sprintf(line, format, ##__VA_ARGS__); \
	                                        \
	  mpi_printf(line);                     \
	}
#endif

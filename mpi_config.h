#if !defined(MPI_CONFIG_HEADER)
	#define CONFIG_HEADER

	#if defined(USE_MPI)
		#include <mpi.h>
	#endif

	#if !defined(MPI_STDOUT_FORMAT)
		#define MPI_STDOUT_FORMAT "cpu=%d_arrang=%c_J=%d.log"
	#endif

	/******************************************************************************

	 Function open_mpi_stdout(): opens a C stdout for a given MPI CPU.

	 NOTE: Each file shall have a unique filename based on the arrangement and J of
	 the problem in order to avoid overlaps with outputs of other calculations.

	******************************************************************************/

	inline static void open_mpi_stdout(const int cpu_id, const char arrang,
	                                   const int J, const bool to_append)
	{
		char filename[MAX_LINE_LENGTH];
		sprintf(filename, MPI_STDOUT_FORMAT, cpu_id, arrang, J);

		file_init_stdout(filename, to_append);
	}
#endif

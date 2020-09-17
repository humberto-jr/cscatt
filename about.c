#include "modules/pes.h"
#include "modules/matrix.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

int main()
{
	printf("# build date   = %s\n", __DATE__);
	printf("# source code  = %s\n", __FILE__);

	printf("\n");
	matrix_about(stdout);

	printf("\n");
	pes_about(stdout);

	printf("\n");
	mpi_about(stdout);

	return EXIT_SUCCESS;
}

#include "modules/mpi_lib.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

	srand(time(NULL));

	const int last_rank = mpi_comm_size() - 1;

	if (mpi_rank() == 0)
	{
		double tag = (double) rand()/RAND_MAX;

		printf("# Rank 0 created  tag %f\n", tag);

		if (mpi_comm_size() > 1) mpi_send(1, 1, type_double, &tag);
	}
	else
	{
		double tag = 0.0;

		mpi_receive(mpi_rank() - 1, 1, type_double, &tag);

		printf("# Rank %d received tag %f from rank %d\n", mpi_rank(), tag, mpi_rank() - 1);

		if (mpi_rank() < last_rank) mpi_send(mpi_rank() + 1, 1, type_double, &tag);
	}

	mpi_end();
	return EXIT_SUCCESS;
}

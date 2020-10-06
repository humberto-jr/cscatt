#include "modules/file.h"
#include "modules/mpi_lib.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	mpi_init(argc, argv);

	if (argc != 6)
	{
		PRINT_ERROR("%d arguments given. Usage: ./sparse_eigen.out [filename] [n_max] [max_step] [tolerance] [upper/lower (true/false)]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	FILE *input = file_open(argv[1], "r");

	int max_row = 0, max_col = 0, non_zeros[] = {1, 1}, info = 0;

	info = fscanf(input, "%d %d", &max_row, &max_col);

	ASSERT(info == 2)

	mpi_matrix *a = mpi_matrix_alloc(max_row, max_col, non_zeros);

	while (true)
	{
		int p, q;
		double value;
		info = fscanf(input, "%d %d %lf", &p, &q, &value);

		if (info != 3) break;

		mpi_matrix_set(a, p, q, value);
	}

	file_close(&input);

	mpi_matrix_build(a);

	const int n_max = atoi(argv[2]), max_step = atoi(argv[3]);

	const double tolerance = atof(argv[4]);

	const bool up = (bool) atoi(argv[5]);

	mpi_matrix_sparse_eigen(a, n_max, max_step, tolerance, up);

	for (int n = 0; n < n_max; ++n)
	{
		double eigenval = 0.0;
		mpi_vector *eigenvec = mpi_matrix_eigenpair(a, n, &eigenval);

		printf("% -8e\n", eigenval);

		mpi_vector_free(eigenvec);
	}

	mpi_matrix_free(a);

	mpi_end();
	return EXIT_SUCCESS;
}

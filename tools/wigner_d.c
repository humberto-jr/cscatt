#include "modules/math.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		PRINT_ERROR("%d arguments given. Usage: %s [k] [m] [j] [beta, min] [beta, max] [beta, size]\n", argc - 1, argv[0])
		return EXIT_FAILURE;
	}

	const double k = atof(argv[1]);
	const double m = atof(argv[2]);
	const double j_max = atof(argv[3]);
	const double beta_min = atof(argv[4]);
	const double beta_max = atof(argv[5]);
	const int n_max = atoi(argv[6]);

	const double beta_step = (beta_max - beta_min)/as_double(n_max);

	for (int n = 0; n < n_max; ++n)
	{
		const double beta = beta_min + as_double(n)*beta_step;

		double *d = math_wigner_d(k, m, j_max, beta);

		printf("%06f", beta);

		for (int j = 0; j <= as_int(j_max); ++j)
			printf("\t % -8e", d[j]);

		printf("\n");

		free(d);
	}

	return EXIT_SUCCESS;
}

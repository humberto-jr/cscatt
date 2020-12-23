#include "modules/math.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		PRINT_ERROR("%d arguments given. Usage: ./legendre_poly.out [l] [x, min] [x, max] [points]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const int l = atoi(argv[1]);

	const double x_min = atof(argv[2]);
	const double x_max = atof(argv[3]);

	const int n_max = atoi(argv[4]);

	const double x_step = (x_max - x_min)/as_double(n_max);

	for (int n = 0; n < n_max; ++n)
	{
		const double x = x_min + as_double(n)*x_step;

		printf("%6f\t % -8e\n", x, math_legendre_poly(l, x));
	}

	return EXIT_SUCCESS;
}

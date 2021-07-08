#include "modules/math.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		PRINT_ERROR("%d arguments given. Usage: ./%s [l] [m] [theta] [phi]\n", argc - 1, argv[0])
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const size_t l = atoi(argv[1]);
	const int m = atoi(argv[2]);
	const double theta = atof(argv[3]);
	const double phi = atof(argv[4]);

	const double complex y = math_sphe_harmonics(l, m, theta, phi);

	printf("% -8e + i(% -8e)\n", creal(y), cimag(y));

	return EXIT_SUCCESS;
}

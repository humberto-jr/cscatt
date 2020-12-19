#include "modules/math.h"
#include "modules/globals.h"

double f(const int J,
	      const int j1,
	      const int j2,
	      const int l1,
	      const int l2,
	      const int lambda)
{
	double result = pow(-1.0, j1 + j2 - J);

	result *= math_clebsch_gordan(j1, j2, lambda, 0, 0, 0);
	result *= math_clebsch_gordan(l1, l2, lambda, 0, 0, 0);
	result *= math_racah_coeff(j1, l1, j2, l2, J, lambda);
	result *= sqrt(as_double(2*j1 + 1)*as_double(2*j2 + 1)*as_double(2*l1 + 1)*as_double(2*l2 + 1))/as_double(2*lambda + 1);

	return result;
}

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		PRINT_ERROR("%d arguments given. Usage: ./percival_seaton.out [J] [j1] [j2] [l1] [l2] [lambda]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const int J = atoi(argv[1]);
	const int j1 = atoi(argv[2]);
	const int j2 = atoi(argv[3]);
	const int l1 = atoi(argv[4]);
	const int l2 = atoi(argv[5]);
	const int lambda = atoi(argv[6]);

	printf("% -8e\n", math_percival_seaton(J, j1, j2, l1, l2, lambda));
	printf("% -8e\n", f(J, j1, j2, l1, l2, lambda));
	printf("% -8e\n", math_ps(J, j1, j2, l1, l2, lambda));

	return EXIT_SUCCESS;
}

#include "modules/math.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		PRINT_ERROR("%d arguments given. Usage: ./gaunt.out [j1] [j2] [lambda]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const int j1 = atoi(argv[1]);
	const int j2 = atoi(argv[2]);
	const int lambda = atoi(argv[3]);

	for (int k = -j1; k <= j1; ++k)
	{
		printf("% -d\t % -8e\n", k, math_gaunt_coeff(k, j1, j2, lambda));
	}

	return EXIT_SUCCESS;
}

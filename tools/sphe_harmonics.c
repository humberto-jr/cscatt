#include "modules/math.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 5)
	{
		PRINT_ERROR("%d arguments given. Usage: ./sphe_harmonics.out [l] [m] [theta] [phi]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const int a = atoi(argv[1]);
	const int b = atoi(argv[2]);
	const int c = atof(argv[3]);
	const int d = atof(argv[4]);

	printf("% -8e\n", math_sphe_harmonics(a, b, c, d));

	return EXIT_SUCCESS;
}

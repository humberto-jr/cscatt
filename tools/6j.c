#include "modules/math.h"
#include "modules/globals.h"

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		PRINT_ERROR("%d arguments given. Usage: ./6j.out [a] [b] [c] [d] [e] [f]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	const int a = atoi(argv[1]);
	const int b = atoi(argv[2]);
	const int c = atoi(argv[3]);
	const int d = atoi(argv[4]);
	const int e = atoi(argv[5]);
	const int f = atoi(argv[6]);

	printf("% -8e\n", math_wigner_6j(a, b, c, d, e, f));

	return EXIT_SUCCESS;
}

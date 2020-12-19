#include "modules/math.h"
#include "modules/globals.h"

#define RESULT (sqrt(5.0) + 0.5*log(2.0 + sqrt(5.0)))

double test(const double x, const void *params)
{
	ASSERT(params == NULL)

	return sqrt(1.0 + x*x);
}

int main()
{
	printf("# Integral of sqrt(1 + x*x) from 0 to 2 is sqrt(5) + 0.5*log(2 + sqrt(5))\n");
	printf("# points	   result	      error	 time (s)\n");

	for (int n = 2; n <= 64; ++n)
	{
		const double start_time = wall_time();
		const double result = math_gauss_legendre(0.0, 2.0, n, NULL, test);
		const double end_time = wall_time();

		printf("   %5d\t %f\t % 8e\t %f\n", n, result, RESULT - result, end_time - start_time);
	}

	return EXIT_SUCCESS;
}

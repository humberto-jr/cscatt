#include "modules/math.h"
#include "modules/globals.h"

#define RESULT 2.9422553486075

double test(const double x, const void *params)
{
	ASSERT(params == NULL)

	return 1.0/(x*x + 1.0);
}

int main()
{
	printf("# Integral of 1/(x*x + 1) from -10 to 10 is 2.9422553486075\n");
	printf("# points	   result	      error	 time (s)\n");

	for (int n = 50; n < 100000; n += 250)
	{
		double result = 0.0, start_time = 0.0, end_time = 0.0;

		start_time = wall_time();
		result = math_simpson(-10.0, 10.0, n, NULL, false, test);
		end_time = wall_time();

		printf("   %5d\t %f\t % 8e\t %f\n", n, result, RESULT - result, end_time - start_time);

		start_time = wall_time();
		result = math_simpson(-10.0, 10.0, n + 1, NULL, false, test);
		end_time = wall_time();

		printf("   %5d\t %f\t % 8e\t %f\n", n + 1, result, RESULT - result, end_time - start_time);
	}

	return EXIT_SUCCESS;
}

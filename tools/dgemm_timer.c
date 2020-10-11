#include "modules/matrix.h"
#include "modules/globals.h"

int main()
{
	matrix_init_gpu();

	matrix_about(stdout);
	printf("# size	 time (s)\n");

	for (int n = 128; n < 2560; n += 32)
	{
		matrix *a = matrix_alloc(n, n, false);
		matrix *b = matrix_alloc(n, n, false);
		matrix *c = matrix_alloc(n, n, false);

		matrix_set_all(a, 1.0, false);
		matrix_set_all(b, 2.0, false);
		matrix_set_all(c, 3.0, false);

		const double result = as_double(n)*2.0 + 3.0;

		const double start_time = wall_time();

		matrix_multiply(1.0, a, b, 1.0, c);

		const double end_time = wall_time();

		for (int p = 0; p < n; ++p)
		{
			for (int q = 0; q < n; ++q)
			{
				const double error = fabs(matrix_get(c, p, q) - result);

				if (error > 1.0E-8)
				{
					PRINT_ERROR("%s", "matrix_multiply() failed with abs. error = ")
					PRINT_ERROR("%-8e\n at (%d, %d)\n", error, p, q)
				}
			}
		}

		printf(" %5d\t %f\n", n, end_time - start_time);

		matrix_free(a);
		matrix_free(b);
		matrix_free(c);
	}

	matrix_end_gpu();
	return EXIT_SUCCESS;
}

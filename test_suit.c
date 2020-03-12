#include "modules/integral.h"
#include "modules/spherical.h"
#include "modules/cartesian.h"
#include "modules/matrix.h"
#include "modules/globals.h"

int main()
{
	printf("\nmatrix_init_gpu();\n");
	matrix_init_gpu();

	printf("\nmatrix *a = matrix_alloc(5, 3, false);\n");
	matrix *a = matrix_alloc(5, 3, false);

	printf("\nmatrix_set_random(a, false);\n");
	matrix_set_random(a, false);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix *b = matrix_alloc(2, 3, false);\n");
	matrix *b = matrix_alloc(2, 3, false);

	printf("\nmatrix_set_random(b, false);\n");
	matrix_set_random(b, false);

	printf("\nmatrix_write(b, stdout, 2, 3);\n");
	matrix_write(b, stdout, 2, 3);

	printf("\nmatrix_swap(a, b);\n");
	matrix_swap(a, b);

	printf("\nmatrix_write(a, stdout, 2, 3);\n");
	matrix_write(a, stdout, 2, 3);

	printf("\nmatrix_write(b, stdout, 5, 3);\n");
	matrix_write(b, stdout, 5, 3);

	printf("\nmatrix_swap(a, b);\n");
	matrix_swap(a, b);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_write(b, stdout, 2, 3);\n");
	matrix_write(b, stdout, 2, 3);

	printf("\nmatrix_min(a);\n");
	printf(" %f\n", matrix_min(a));

	printf("\nmatrix_max(a);\n");
	printf(" %f\n", matrix_max(a));

	printf("\nmatrix_set_all(a, 0.0, false);\n");
	matrix_set_all(a, 0.0, false);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_set(a, 3, 2, 3.14);\n");
	matrix_set(a, 3, 2, 3.14);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_sum(a, false);\n");
	printf(" %f\n", matrix_sum(a, false));

	printf("\nmatrix_col_sum(a, 2, false);\n");
	printf(" %f\n", matrix_col_sum(a, 2, false));

	printf("\nmatrix_row_sum(a, 3, false);\n");
	printf(" %f\n", matrix_row_sum(a, 3, false));

	printf("\nmatrix_set(a, 3, 0, 1.0);\n");
	matrix_set(a, 3, 0, 1.0);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_row_sum(a, 3, false);\n");
	printf(" %f\n", matrix_row_sum(a, 3, false));

	printf("\nmatrix_incr_all(a, 1.0, false);\n");
	matrix_incr_all(a, 1.0, false);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_decr_all(a, 1.0, false);\n");
	matrix_decr_all(a, 1.0, false);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_decr(a, 3, 0, 1.0);\n");
	matrix_decr(a, 3, 0, 1.0);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix_decr(a, 3, 2, 3.14);\n");
	matrix_decr(a, 3, 2, 3.14);

	printf("\nmatrix_write(a, stdout, 5, 3);\n");
	matrix_write(a, stdout, 5, 3);

	printf("\nmatrix *c = matrix_alloc(9, 9, false);\n");
	matrix *c = matrix_alloc(9, 9, false);

	printf("\nmatrix_set_all(c, 0.0, false);\n");
	matrix_set_all(c, 0.0, false);

	printf("\nmatrix_write(c, stdout, 9, 9);\n");
	matrix_write(c, stdout, 9, 9);

	printf("\nmatrix_null(c);\n");
	printf(" %s\n", (matrix_null(c)? "true" : "false"));

	printf("\nmatrix_negative(c);\n");
	printf(" %s\n", (matrix_negative(c)? "true" : "false"));

	printf("\nmatrix_positive(c);\n");
	printf(" %s\n", (matrix_positive(c)? "true" : "false"));

	printf("\nmatrix_diag_set_all(c, 1.0, false);\n");
	matrix_diag_set_all(c, 1.0, false);

	printf("\nmatrix_write(c, stdout, 9, 9);\n");
	matrix_write(c, stdout, 9, 9);

	printf("\nmatrix_null(c);\n");
	printf(" %s\n", (matrix_null(c)? "true" : "false"));

	printf("\nmatrix_negative(c);\n");
	printf(" %s\n", (matrix_negative(c)? "true" : "false"));

	printf("\nmatrix_positive(c);\n");
	printf(" %s\n", (matrix_positive(c)? "true" : "false"));

	printf("\nmatrix_trace(c, false);\n");
	printf(" %f\n", matrix_trace(c, false));

	printf("\nmatrix_sum(c, false);\n");
	printf(" %f\n", matrix_sum(c, false));

	printf("\nmatrix_symm_set(c, 5, 7, 3.14);\n");
	matrix_symm_set(c, 5, 7, 3.14);

	printf("\nmatrix_write(c, stdout, 9, 9);\n");
	matrix_write(c, stdout, 9, 9);

	printf("\nmatrix_set_random(c, false);\n");
	matrix_set_random(c, false);

	printf("\nmatrix_write(c, stdout, 9, 9);\n");
	matrix_write(c, stdout, 9, 9);

	printf("\nmatrix *d = matrix_get_row(c, 6, false);\n");
	matrix *d = matrix_get_row(c, 6, false);

	printf("\nmatrix_write(d, stdout, 9, 9);\n");
	matrix_write(d, stdout, 9, 9);

	printf("\nmatrix *e = matrix_get_col(c, 1, false);\n");
	matrix *e = matrix_get_col(c, 1, false);

	printf("\nmatrix_write(e, stdout, 9, 9);\n");
	matrix_write(e, stdout, 9, 9);

	printf("\nmatrix *f = matrix_get_diag(c, false);\n");
	matrix *f = matrix_get_diag(c, false);

	printf("\nmatrix_write(f, stdout, 9, 9);\n");
	matrix_write(f, stdout, 9, 9);

	printf("\nmatrix_row_scale(c, 3, 1.0E-6, false);\n");
	matrix_row_scale(c, 3, 1.0E-6, false);

	printf("\nmatrix_write(c, stdout, 9, 9);\n");
	matrix_write(c, stdout, 9, 9);

	printf("\nmatrix_col_scale(c, 0, 1.0E-6, false);\n");
	matrix_col_scale(c, 0, 1.0E-6, false);

	printf("\nmatrix_write(c, stdout, 9, 9);\n");
	matrix_write(c, stdout, 9, 9);

	printf("\nmatrix *g = matrix_get_block(c, 3, 7, 0, 3);\n");
	matrix *g = matrix_get_block(c, 3, 7, 0, 3);

	printf("\nmatrix_write(g, stdout, 9, 9);\n");
	matrix_write(g, stdout, 9, 9);

	printf("\nmatrix *h = matrix_alloc(3, 3, false);\n");
	matrix *h = matrix_alloc(3, 3, false);

	printf("\nmatrix_set(h, 0, 0, -4.0);\n");
	matrix_set(h, 0, 0, -4.0);

	printf("\nmatrix_set(h, 0, 1, 2.0);\n");
	matrix_set(h, 0, 1, 2.0);

	printf("\nmatrix_set(h, 0, 2, -2.0);\n");
	matrix_set(h, 0, 2, -2.0);

	printf("\nmatrix_set(h, 1, 0, 2.0);\n");
	matrix_set(h, 1, 0, 2.0);

	printf("\nmatrix_set(h, 1, 1, -7.0);\n");
	matrix_set(h, 1, 1, -7.0);

	printf("\nmatrix_set(h, 1, 2, 4.0);\n");
	matrix_set(h, 1, 2, 4.0);

	printf("\nmatrix_set(h, 2, 0, -2.0);\n");
	matrix_set(h, 2, 0, -2.0);

	printf("\nmatrix_set(h, 2, 1, 4.0);\n");
	matrix_set(h, 2, 1, 4.0);

	printf("\nmatrix_set(h, 2, 2, -7.0);\n");
	matrix_set(h, 2, 2, -7.0);

	printf("\nmatrix_write(h, stdout, 3, 3);\n");
	matrix_write(h, stdout, 3, 3);

	printf("\ndouble *eigenval = matrix_symm_eigen(h, 'v');\n");
	double *eigenval = matrix_symm_eigen(h, 'v');

	printf("eigenval[0] = %f, error = %f\n", eigenval[0], eigenval[0] + 12.0);
	printf("eigenval[1] = %f, error = %f\n", eigenval[1], eigenval[1] + 3.0);
	printf("eigenval[2] = %f, error = %f\n", eigenval[2], eigenval[2] + 3.0);

	printf("\nmatrix_free(a);\n");
	matrix_free(a);

	printf("\nmatrix_free(b);\n");
	matrix_free(b);

	printf("\nmatrix_free(c);\n");
	matrix_free(c);

	printf("\nmatrix_free(d);\n");
	matrix_free(d);

	printf("\nmatrix_free(e);\n");
	matrix_free(e);

	printf("\nmatrix_free(f);\n");
	matrix_free(f);

	printf("\nmatrix_free(g);\n");
	matrix_free(g);

	printf("\nmatrix_free(h);\n");
	matrix_free(h);

	free(eigenval);

	printf("\ncartesian x = {.x = 3.0, .y = 4.0, .z = 5.0};\n");
	cartesian x = {.x = 3.0, .y = 4.0, .z = 5.0};

	printf("\nspherical y = {.rho = 0.0, .theta = 0.0, .phi = 0.0};\n");
	spherical y = {.rho = 0.0, .theta = 0.0, .phi = 0.0};

	printf("\ncartesian_to_spherical(&x, &y);\n");
	cartesian_to_spherical(&x, &y);

	printf("y.rho = %f, y.theta = %f, y.phi = %f\n", y.rho, y.theta, y.phi);

	printf("\ncartesian z = {.x = 0.0, .y = 0.0, .z = 0.0};\n");
	cartesian z = {.x = 0.0, .y = 0.0, .z = 0.0};

	printf("\ncartesian_from_spherical(&y, &z);\n");
	cartesian_from_spherical(&y, &z);

	printf("z.x = %f, z.y = %f, z.z = %f\n", z.x, z.y, z.z);

	printf("\n");
	printf("Timing for integral_benchmark(simpson_1st, n, &error, &time) for n points; n vs. error, time\n");

	for (int n = 100; n < 1000; n += 50)
	{
		double error = 0.0, time = 0.0;

		integral_benchmark(simpson_1st, n, &error, &time);
		printf("%4d\t % -8e\t % -8e\t\n", n, error, time);
	}

	printf("\n");
	printf("Timing for integral_benchmark(simpson_2nd, n, &error, &time) for n points; n vs. error, time\n");

	for (int n = 100; n < 1000; n += 50)
	{
		double error = 0.0, time = 0.0;

		integral_benchmark(simpson_2nd, n, &error, &time);
		printf("%4d\t % -8e\t % -8e\t\n", n, error, time);
	}

/*
	printf("\n");
	printf("# Test of matrix_inv() #######################################\n");
	{
		matrix *m = matrix_alloc(1, 3, 3, false);

		matrix_set(m, 0, 0, 1.0);
		matrix_set(m, 0, 1, 2.0);
		matrix_set(m, 0, 2, 3.0);
		matrix_set(m, 1, 0, 2.0);
		matrix_set(m, 1, 1, 1.0);
		matrix_set(m, 1, 2, 4.0);
		matrix_set(m, 2, 0, 3.0);
		matrix_set(m, 2, 1, 4.0);
		matrix_set(m, 2, 2, 0.0);

		printf("# Matrix M:\n");
		write_data(stdout, 3, 3, m->data);

		matrix_inv(m);

		printf("\n");
		printf("# Matrix M^-1:\n");
		write_data(stdout, 3, 3, m->data);

		matrix_decr(m, 0, 0, -0.695652);
		matrix_decr(m, 0, 1,  0.521739);
		matrix_decr(m, 0, 2,  0.217391);
		matrix_decr(m, 1, 0,  0.521739);
		matrix_decr(m, 1, 1, -0.391304);
		matrix_decr(m, 1, 2,  0.086956);
		matrix_decr(m, 2, 0,  0.217391);
		matrix_decr(m, 2, 1,  0.086956);
		matrix_decr(m, 2, 2, -0.130435);

		printf("\n");
		printf("# Approx. error (to six digits):\n");
		write_data(stdout, 3, 3, m->data);

		matrix_free(m);
	}

	printf("\n");
	printf("# Test of multich_renorm_numerov() ###########################\n");
	{
		const double x_min = 0.050;
		const double x_max = 6.000;
		const double x_inc = 0.001;

		const int grid_size = (int) (x_max - x_min)/x_inc;

		const double coll_energy = 2.60566;
*/
		/* NOTE: 4He + 20Ne (atomic units) */
/*		const double mass = 6089.0;

		matrix *ratio = matrix_alloc(1, 2, 2, true);

		for (int n = 0; n < grid_size; ++n)
		{
			const double x = x_min + as_double(n)*x_inc;

			matrix *v = olson_smith_model(x);
			johnson_numerov(x_inc, mass, coll_energy, v, ratio);

			matrix_free(v);
		}
	}
*/
	printf("\n");
	printf("Timing for matrix_multiply(1.0, a, b, 1.0, c); n-by-n; n vs. wall time (s)\n");

	for (int n = 128; n < 2560; n += 32)
	{
		a = matrix_alloc(n, n, false);
		b = matrix_alloc(n, n, false);
		c = matrix_alloc(n, n, false);

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

		printf("%4d\t%f\n", n, end_time - start_time);

		matrix_free(a);
		matrix_free(b);
		matrix_free(c);
	}

/*
	printf("\n");
	printf("# Timing of matrix_inv(): n-by-n; n vs. wall time (s)\n");

	for (int n = 128; n < 2560; n += 32)
	{
		matrix *m = matrix_alloc(1, n, n, false);
		matrix_set_random(m, false);

		const double start_time = wall_time();

		matrix_inv(m);

		const double end_time = wall_time();

		printf("%4d\t%f\n", n, end_time - start_time);

		matrix_free(m);
	}
*/

	printf("\nmatrix_end_gpu();\n");
	matrix_end_gpu();

	return EXIT_SUCCESS;
}

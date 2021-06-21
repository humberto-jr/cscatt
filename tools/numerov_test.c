#include "modules/pes.h"
#include "modules/matrix.h"
#include "modules/johnson.h"
#include "modules/globals.h"

int main()
{
	matrix_init_gpu();

	const double x_min = 0.050;
	const double x_max = 6.000;
	const double x_step = 0.001;
	const size_t grid_size = (size_t) (x_max - x_min)/x_step;

	/* NOTE: collision energy of 70.9 eV (2.60566 a.u.). */
	const double coll_energy = 2.60566;

	/* NOTE: 4He + 20Ne (atomic units). */
	const double mass = 6089.0;

	/* NOTE: list of partial waves used at table I of Ref. [1]. */
	const double l[] =
	{
		  0.0,
		 10.0,
		 20.0,
		200.0,
		202.0,
		205.0,
		209.0,
		214.0,
		290.0,
		300.0,
		310.0,
		320.0,
		330.0,
		350.0
	};

	/* NOTE: from table I of Ref. [1]. */
	const double result[] =
	{
		0.279E-3,
		0.164E-2,
		0.110E-1,
		0.215E-1,
		0.552E-1,
		0.900E-1,
		0.553E-1,
		0.487E-3,
		0.167E-4,
		0.146,
		0.202,
		0.941E-1,
		0.226E-1,
		0.329E-3
	};

	const double levels[] =
	{
		pes_olson_smith_model(0, 0, x_max),
		pes_olson_smith_model(1, 1, x_max)
	};

	matrix *r = matrix_alloc(2, 2, true);
	matrix *v = matrix_alloc(2, 2, false);

	printf("#  l      Numerov        Ref. [1]           Error\n");
	printf("# -----------------------------------------------\n");

	for (size_t m = 0; m < 14; ++m)
	{
		const int l_list[] = {(int) l[m], (int) l[m]};

		for (size_t n = 0; n < grid_size; ++n)
		{
			const double x = x_min + as_double(n)*x_step;

			const double l_term = l[m]*(l[m] + 1.0)/(2.0*mass*x*x);

			matrix_set(v, 0, 0, pes_olson_smith_model(0, 0, x) + l_term);
			matrix_set(v, 0, 1, pes_olson_smith_model(0, 1, x));
			matrix_set(v, 1, 0, pes_olson_smith_model(1, 0, x));
			matrix_set(v, 1, 1, pes_olson_smith_model(1, 1, x) + l_term);

			johnson_jcp78_numerov(x_step, mass, coll_energy, v, r, NULL);
		}

		matrix *k
			= johnson_kmatrix(l_list, x_step, coll_energy, mass, levels, r, x_max);

		smatrix *s = johnson_smatrix(k);

		const double s01_real = matrix_get(s->re_part, 0, 1);
		const double s01_imag = matrix_get(s->im_part, 0, 1);
		const double s01 = s01_real*s01_real + s01_imag*s01_imag;

		printf(" %3d\t %f\t %f\t %f\n", (int) l[m], s01, result[m], fabs(s01 - result[m]));

		matrix_free(s->re_part);
		matrix_free(s->im_part);
		matrix_free(k);
		free(s);
	}

	printf("# -----------------------------------------------\n");
	printf("# [1] B. R. Johnson. Journal of Computational Physics, 13, 445-449 (1973)\n");

	matrix_free(v);
	matrix_free(r);

	matrix_end_gpu();
	return EXIT_SUCCESS;
}

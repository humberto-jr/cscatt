/******************************************************************************

 Function pes_olson_smith_model(): return the 2-by-2 model potential of Olson
 and Smith for the system Ne + He^+ at a given internuclear distance x. See
 Eq. (45) and Table I from Ref. [1] for details.

 NOTE: a handy benchmark for the algorithms.

******************************************************************************/

matrix *pes_olson_smith_model(const double x)
{
	matrix *m = matrix_alloc(2, 2, false);

	matrix_set(m, 0, 0, 21.1*exp(-x/0.678)/x);
	matrix_set(m, 0, 1, 0.170*exp(-x/0.667));
	matrix_set(m, 1, 0, 0.170*exp(-x/0.667));
	matrix_set(m, 1, 1, (21.1/x - 12.1)*exp(-x/0.678) + 16.8/27.2113961);

	return m;
}

/******************************************************************************

 Function pes_tully_model(): return one of the three (n = 1, 2, 3) 2-by-2 model
 potentials of J. C. Tully at a given x. See Eq. (21), (23) and (24) of
 Ref. [2].

******************************************************************************/

matrix *pes_tully_model(const int n, const double x)
{
	matrix *m = matrix_alloc(2, 2, false);

	switch (n)
	{
		case 1:
			if (x > 0.0)
			{
				matrix_set(m, 0, 0,  0.01*(1.0 - exp(-1.6*x)));
				matrix_set(m, 1, 1, -0.01*(1.0 - exp(-1.6*x)));
			}
			else
			{
				matrix_set(m, 0, 0, -0.01*(1.0 - exp(1.6*x)));
				matrix_set(m, 1, 1,  0.01*(1.0 - exp(1.6*x)));
			}

			matrix_set(m, 0, 1, 0.005*exp(-1.0*x*x));
			matrix_set(m, 1, 0, 0.005*exp(-1.0*x*x));
		break;

		case 2:
			matrix_set(m, 0, 0,  0.0);
			matrix_set(m, 0, 1,  0.015*exp(-0.06*x*x));
			matrix_set(m, 1, 0,  0.015*exp(-0.06*x*x));
			matrix_set(m, 1, 1, -0.10*exp(-0.28*x*x) + 0.05);
		break;

		case 3:
			matrix_set(m, 0, 0,  0.0006);
			matrix_set(m, 1, 1, -0.0006);

			matrix_set(m, 1, 0, 0.0);

			if (x < 0.0)
				matrix_set(m, 0, 1, 0.10*exp(0.90*x));

			else
				matrix_set(m, 0, 1, 0.10*(2.0 - exp(-0.90*x)));
		break;

		default:
			PRINT_ERROR("invalid index n = %d\n", n)
			exit(EXIT_FAILURE);
	}

	return m;
}

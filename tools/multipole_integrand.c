#include "modules/pes.h"
#include "modules/math.h"
#include "modules/file.h"
#include "modules/globals.h"

#define THETA_GRID_SIZE 1000

#define HEADER "#   theta\t       V(r, R)\t    V'(r, inf)\t        V - V'\t      Legendre\t           sin\t       product\n"

#define FORMAT "% 06f\t % -8e\t % -8e\t % -8e\t % -8e\t % -8e\t % -8e\n"

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		PRINT_ERROR("%d arguments given. Usage: ./multipole_integrand.out [file]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	file_init_stdin(argv[1]);;

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

	const double r = file_keyword(stdin, "r_value", 0.0, INF, 0.0);

	const double R = file_keyword(stdin, "R_value", 0.0, INF, 0.0);

	const int lambda = (int) file_keyword(stdin, "lambda_value", 0.0, INF, 0.0);

	const char arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0);

	const double multipole = pes_legendre_multipole(arrang, lambda, r, R);

	const double theta_step = M_PI/THETA_GRID_SIZE;

	printf("# arrang = %c, r = %f, R = %f, lambda = %d, result = %f\n",
	       arrang, r, R, lambda, multipole);

	printf(HEADER);

	for (int n = 0; n < THETA_GRID_SIZE; ++n)
	{
		const double theta = as_double(n)*theta_step;

		const double sinx = sin(theta);

		const double v = pes_abc(arrang, r, R, theta*180.0/M_PI);

		const double v_inf = pes_abc(arrang, r, 1000.0, theta*180.0/M_PI);

		const double v_diff = v - v_inf;

		const double legendre = math_legendre_poly(lambda, cos(theta));

		const double result = v_diff*legendre*sinx;

		printf(FORMAT, theta, v, v_inf, v_diff, legendre, sinx, result);
	}

	return EXIT_SUCCESS;
}

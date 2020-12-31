#include "modules/pes.h"
#include "modules/math.h"
#include "modules/file.h"
#include "modules/globals.h"

#define HEADER "#   theta\t       V(r, R)\t    V'(r, inf)\t        V - V'\t      Legendre\t           sin\t       product\n"

#define FORMAT "%06f\t % -8e\t % -8e\t % -8e\t % -8e\t % -8e\t % -8e\n"

#define MAX_ORDER 64

double x[MAX_ORDER];
double v[MAX_ORDER];
double sinx[MAX_ORDER];
double v_inf[MAX_ORDER];
double v_diff[MAX_ORDER];
double result[MAX_ORDER];
double legendre[MAX_ORDER];

struct legendre_params
{
	const char arrang;
	const int lambda;
	const double r;
	const double R;
};

double legendre_integrand(const double theta, const void *params)
{
	static int n = 0;

	ASSERT(n < MAX_ORDER)

	struct legendre_params *p = (struct legendre_params *) params;

	x[n] = theta*180.0/M_PI;

	sinx[n] = sin(theta);

	v[n] = pes_abc(p->arrang, p->r, p->R, theta*180.0/M_PI);

	v_inf[n] = pes_abc(p->arrang, p->r, 1000.0, theta*180.0/M_PI);

	v_diff[n] = v[n] - v_inf[n];

	legendre[n] = math_legendre_poly(p->lambda, cos(theta));

	result[n] = v_diff[n]*legendre[n]*sinx[n];

	++n;
	return result[n];
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		PRINT_ERROR("%d arguments given. Usage: ./multipole_integrand.out [input file]\n", argc - 1)
		return EXIT_FAILURE;
	}

	ASSERT(argv != NULL)

	file_init_stdin(argv[1]);;

	pes_init_mass(stdin, 'a');
	pes_init_mass(stdin, 'b');
	pes_init_mass(stdin, 'c');

	pes_init();

	struct legendre_params p =
	{
		.r = file_keyword(stdin, "r_value", 0.0, INF, 0.0),

		.R = file_keyword(stdin, "R_value", 0.0, INF, 0.0),

		.lambda = (int) file_keyword(stdin, "lambda_value", 0.0, INF, 0.0),

		.arrang = 96 + (int) file_keyword(stdin, "arrang", 1.0, 3.0, 1.0)
	};

	const double multipole
		= math_gauss_legendre(0.0, M_PI, MAX_ORDER, &p, legendre_integrand);

	printf("# arrang = %c, r = %f, R = %f, lambda = %d, result = %f\n",
	       p.arrang, p.r, p.R, p.lambda, as_double(2*p.lambda + 1)*multipole/2.0);

	printf(HEADER);

	int *index = math_bubble_sort(MAX_ORDER, x);

	for (int n = 0; n < MAX_ORDER; ++n)
	{
		const int m = index[n];
		printf(FORMAT, x[m], v[m], v_inf[m], v_diff[m], legendre[m], sin(x[m]), result[m]);
	}

	free(index);
	return EXIT_SUCCESS;
}

#include "modules/math.h"
#include "modules/file.h"
#include "modules/globals.h"

struct data
{
	int m, size;
	double *theta, *sigma;
};

FILE *open_input(const int v_in,
                 const int j_in,
                 const int m_in,
                 const int v_out,
                 const int j_out,
                 const int m_out, const int J)
{
	char filename[MAX_LINE_LENGTH];

	sprintf(filename, "dif_csection_iv=%d_ij=%d_im=%d_fv=%d_fj=%d_fm=%d_J=%d.dat", v_in, j_in, m_in, v_out, j_out, m_out, J);
	printf("# Using %s\n", filename);

	return file_open(filename, "r");
}

void read_csection(FILE *input, const double coll_energy, struct data *d)
{
	ASSERT(input != NULL)
	ASSERT(d->sigma == NULL);
	ASSERT(d->theta == NULL);

	d->size = 0;
	char line[MAX_LINE_LENGTH];

	while (fgets(line, sizeof(line), input) != NULL)
	{
		char *token = strtok(line, " \t");

		if (token == NULL)
		{
			PRINT_ERROR("invalid entry at line '%s' (1)\n", line)
			exit(EXIT_FAILURE);
		}

		if (atof(token) != coll_energy) continue;

		token = strtok(NULL, " \t");

		if (token == NULL)
		{
			PRINT_ERROR("invalid entry at line '%s' (2)\n", line)
			exit(EXIT_FAILURE);
		}

		d->theta = realloc(d->theta, sizeof(double)*(d->size + 1));

		d->theta[d->size] = atof(token);

		token = strtok(NULL, " \t");

		if (token == NULL)
		{
			PRINT_ERROR("invalid entry at line '%s' (3)\n", line)
			exit(EXIT_FAILURE);
		}

		d->sigma = realloc(d->sigma, sizeof(double)*(d->size + 1));

		d->sigma[d->size] = atof(token);

		d->size += 1;
	}

	if (d->size == 0)
	{
		PRINT_ERROR("no data read for E = %f\n", coll_energy)
		exit(EXIT_FAILURE);
	}
}

void sum_csection(struct data *a, const struct data *b)
{
	ASSERT(a->sigma != NULL)
	ASSERT(b->sigma != NULL)
	ASSERT(a->theta != NULL)
	ASSERT(b->theta != NULL)
	ASSERT(a->size == b->size)

	for (int n = 0; n < a->size; ++n)
	{
		a->sigma[n] += b->sigma[n];
		ASSERT(a->theta[n] == b->theta[n])
	}
}

struct data *read_all(const int v_in,
                      const int j_in,
                      const int v_out,
                      const int j_out,
                      const int m_out,
                      const int J_max,
											const double coll_energy)
{
	struct data *d = allocate(2*j_in + 1, sizeof(struct data), true);

	int counter = 0;
	for (int m = -j_in; m <= j_in; ++m)
	{
		/* NOTE: each count is for one m-value summed over all J-values. */

		d[counter].m = m;

		/* J = 0 only */
		FILE *input = open_input(v_in, j_in, m, v_out, j_out, m_out, 0);

		struct data temp = {.theta = NULL, .sigma = NULL};

		read_csection(input, coll_energy, &temp);

		fclose(input);

		d[counter].size = temp.size;
		d[counter].sigma = allocate(temp.size, sizeof(double), true);
		d[counter].theta = allocate(temp.size, sizeof(double), true);

		for (int n = 0; n < temp.size; ++n)
			d[counter].theta[n] = temp.theta[n];

		sum_csection(&d[counter], &temp);

		free(temp.theta);
		free(temp.sigma);

		/* J > 0 */
		for (int J = 1; J <= J_max; ++J)
		{
			temp.sigma = NULL;
			temp.theta = NULL;

			input = open_input(v_in, j_in, m, v_out, j_out, m_out, J);

			read_csection(input, coll_energy, &temp);

			fclose(input);

			sum_csection(&d[counter], &temp);

			free(temp.theta);
			free(temp.sigma);
		}

		++counter;
	}

	ASSERT((2*j_in + 1) == counter)

	return d;
}

int main(int argc, char *argv[])
{
	if (argc != 10)
	{
		PRINT_ERROR("%d arguments given. Usage: %s [v, in] [j, in] [m, exp] [v, out] [j, out] [m, out] [J, max] [coll. energy, eV] [beta]\n", argc - 1, argv[0])
		return EXIT_FAILURE;
	}

	const int v_in = atoi(argv[1]);
	const int j_in = atoi(argv[2]);
	const int m_exp = atoi(argv[3]);
	const int v_out = atoi(argv[4]);
	const int j_out = atoi(argv[5]);
	const int m_out = atoi(argv[6]);
	const int J_max = atoi(argv[7]);
	const double coll_energy = atof(argv[8]);
	const double beta = atof(argv[9]);

	struct data *d = read_all(v_in, j_in, v_out, j_out, m_out, J_max, coll_energy);

	const double a = (double) abs(m_exp);

	for (int n = 0; n < d[0].size; ++n)
	{
		printf("% -8e", d[0].theta[n]);

		double sum = 0.0;
		for (int counter = 0; counter < (2*j_in + 1); ++counter)
		{
			double *wigner_d = NULL;
			const double b = (double) abs(d[counter].m);

			if (a > b)
				wigner_d = math_wigner_d(a, b, as_double(j_in), beta);
			else
				wigner_d = math_wigner_d(b, a, as_double(j_in), beta);

			sum += d[counter].sigma[n]*pow(wigner_d[j_in], 2);

			free(wigner_d);
			printf("\t % -8e", d[counter].sigma[n]/pow(1.0E-8, 2));
		}

		printf("\t % -8e\n", sum/pow(1.0E-8, 2));
	}

	free(d);
	return EXIT_SUCCESS;
}

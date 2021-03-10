#include "modules/math.h"
#include "modules/file.h"
#include "modules/globals.h"

struct data
{
	int m, size;
	double *coll_energy, *sigma;
};

FILE *open_input(const int v_in,
                 const int j_in,
                 const int m_in,
                 const int v_out,
                 const int j_out,
                 const int m_out)
{
	char filename[MAX_LINE_LENGTH];

	sprintf(filename, "int_csection_iv=%d_ij=%d_im=%d_fv=%d_fj=%d_fm=%d.dat", v_in, j_in, m_in, v_out, j_out, m_out);
	printf("# Using %s\n", filename);

	return file_open(filename, "r");
}

void read_csection(FILE *input, struct data *d)
{
	ASSERT(input != NULL)
	ASSERT(d->sigma == NULL);
	ASSERT(d->coll_energy == NULL);

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

		d->coll_energy = realloc(d->coll_energy, sizeof(double)*(d->size + 1));

		d->coll_energy[d->size] = atof(token);

		token = strtok(NULL, " \t");

		if (token == NULL)
		{
			PRINT_ERROR("invalid entry at line '%s' (2)\n", line)
			exit(EXIT_FAILURE);
		}

		d->sigma = realloc(d->sigma, sizeof(double)*(d->size + 1));

		d->sigma[d->size] = atof(token);
		d->size += 1;
	}

	ASSERT(d->size > 0)
}

struct data *read_all(const int v_in,
                      const int j_in,
                      const int v_out,
                      const int j_out,
                      const int m_out)
{
	struct data *d = allocate(2*j_in + 1, sizeof(struct data), true);

	int counter = 0;
	for (int m = -j_in; m <= j_in; ++m)
	{
		FILE *input = open_input(v_in, j_in, m, v_out, j_out, m_out);

		read_csection(input, &d[counter]);
		d[counter].m = m;

		fclose(input);
		++counter;
	}

	ASSERT((2*j_in + 1) == counter)

	return d;
}

int main(int argc, char *argv[])
{
	if (argc != 8)
	{
		PRINT_ERROR("%d arguments given. Usage: %s [v, in] [j, in] [m, in] [v, out] [j, out] [m, out] [beta]\n", argc - 1, argv[0])
		return EXIT_FAILURE;
	}

	const int v_in = atoi(argv[1]);
	const int j_in = atoi(argv[2]);
	const int m_in = atoi(argv[3]);
	const int v_out = atoi(argv[4]);
	const int j_out = atoi(argv[5]);
	const int m_out = atoi(argv[6]);
	const double beta = atof(argv[7]);

	struct data *d = read_all(v_in, j_in, v_out, j_out, m_out);

	const double a = (double) abs(m_in);

	for (int n = 0; n < d[0].size; ++n)
	{
		printf("% -8e", d[0].coll_energy[n]*1.160451812E4);

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
			printf("\t % -8e", d[counter].sigma[n]);
		}

		printf("\t % -8e\n", sum);
	}

	free(d);
	return EXIT_SUCCESS;
}

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "modules/globals.h"
#include "modules/matrix.h"
#include "modules/tools.h"

struct network_layer
{
	double *error;
	double *neuron;
	matrix *weight;
};

struct network
{
	int max_layer, max_neuron;
	struct network_layer *layer;
};

static inline double sigmoid(const double x)
{
	return 1.0/(1.0 + exp(-x));
}
/*
static inline double leaky_rect_lu(const double x)
{
	return (x > 0.0? x : 0.01*x);
}*/
/*
static inline void softmax(const struct network *net, const int l)
{
	register double sum = 0.0;

	for (int n = 0; n < net->max_neuron; ++n)
	{
		sum += exp(net->layer[l].neuron[n]);
	}

	for (int n = 0; n < net->max_neuron; ++n)
	{
		net->layer[l].neuron[n] = exp(net->layer[l].neuron[n])/sum;
	}
}*/
/*
static inline void softmax(double neuron[], const int max_neuron)
{
	register double sum = 0.0;

	for (int n = 0; n < max_neuron; ++n)
	{
		sum += exp(neuron[n]);
	}

	ASSERT(sum > 0.0)

	for (int n = 0; n < max_neuron; ++n)
	{
		neuron[n] = exp(neuron[n])/sum;
	}
}*/
/*
void update(struct network *net, const int l, const double error[])
{
	for (int n = 0; n < net->max_neuron; ++n)
	{
		for (int m = 0; m < net->max_neuron; ++m)
		{
			register const double correction
				= error[m]*net->layer[l - 1].neuron[n];

			matrix_incr(net->layer[l].weight, n, m, correction);
		}
	}
}*/

double *pattern_a(const matrix *a, const double input[])
{
	double *result = allocate(matrix_row(a), false);

	for (int n = 0; n < matrix_row(a); ++n)
	{
		register double sum = 0.0;

		for (int m = 0; m < matrix_col(a); ++m)
		{
			sum += matrix_get(a, n, m)*input[m];
		}

		result[n] = sum;
	}

	return result;
}

void feed_forward(struct network *net,
                  const matrix *input, const int p)
{
	for (int n = 0; n < net->max_layer; ++n)
	{
		for (int j = 0; j < net->max_neuron; ++j)
		{
			register double sum = 0.0;

			if (n > 0)
			{
				for (int k = 0; k < net->max_neuron; ++k)
				{
					register const double activ
						= net->layer[n - 1].neuron[k];

					register const double weight
						= matrix_get(net->layer[n].weight, j, k);

					sum += weight*activ;
				}
			}
			else
			{
				for (int k = 0; k < net->max_neuron; ++k)
				{
					register const double activ
						= matrix_get(input, k, p);

					register const double weight
						= matrix_get(net->layer[0].weight, j, k);

					sum += weight*activ;
				}
			}

			net->layer[n].neuron[j] = sigmoid(sum);
		}
	}
}

void back_propagate(struct network *net,
                    const matrix *input, const int p)
{
	register const int l_max = net->max_layer - 1;

/*
 *	Resolve the error in the last layer:
 */

	for (int n = 0; n < net->max_neuron; ++n)
	{
		register const double output
			= net->layer[l_max].neuron[n];

		register const double result
			= matrix_get(input, n, p);

		net->layer[l_max].error[n]
			= output*(1.0 - output)*(result - output);
	}

/*
 *	Resolve the error in the layer n:
 */

	for (int n = (l_max - 1); n > -1; --n)
	{
		for (int j = 0; j < net->max_neuron; ++j)
		{
			register double sum = 0.0;

			for (int k = 0; k < net->max_neuron; ++k)
			{
				register const double error
					= net->layer[n + 1].error[k];

				register const double weight
					= matrix_get(net->layer[n + 1].weight, k, j);

				/* NOTE: weight is transposed */
				sum += error*weight;
			}

			register const double output
				= net->layer[n].neuron[j];

			net->layer[n].error[j]
				= output*(1.0 - output)*sum;
		}
	}
}

void print_network(const struct network *net, const double input[],
                   const double result[], const double error)
{
	static int step = 0;

	printf("\n====================================================\n");

	printf("step = %d\n", ++step);

	printf("\n");
	for (int j = 0; j < net->max_neuron; ++j)
	{
		printf("input[%d] = %f\n", j, input[j]);
	}

	printf("\n");
	for (int j = 0; j < net->max_neuron; ++j)
	{
		printf("result[%d] = %f\n", j, result[j]);
	}

	printf("\n");
	printf("error = %f\n", error);

	for (int n = 0; n < net->max_layer; ++n)
	{
		printf("\n");
		printf("layer = %d\n", n);


		printf("\nactivations:\n");

		for (int j = 0; j < net->max_neuron; ++j)
		{
			printf(" %f\n", net->layer[n].neuron[j]);
		}

		printf("\nerrors:\n");

		for (int j = 0; j < net->max_neuron; ++j)
		{
			printf(" %f\n", net->layer[n].error[j]);
		}

		printf("\nweights:\n");

		matrix_write(net->layer[n].weight, stdout,
		             net->max_neuron, net->max_neuron);
	}

//	printf("\n====================================================\n");
}

void init_net(struct network *net, const int max_layer, const int max_neuron)
{
	net->layer = calloc(max_layer, sizeof(struct network_layer));
	ASSERT(net->layer != NULL)

	for (int l = 0; l < max_layer; ++l)
	{
		net->layer[l].error  = allocate(max_neuron, false);
		net->layer[l].neuron = allocate(max_neuron, false);
		net->layer[l].weight = matrix_alloc(max_neuron, max_neuron, false);

		matrix_set_random(net->layer[l].weight, false);
	}

	net->max_layer = max_layer;
	net->max_neuron = max_neuron;
}

void test(struct network *net)
{
	static int counter = 0;

	matrix *input
		= matrix_alloc(net->max_neuron, 1, false);

	matrix_set_random(input, false);

	feed_forward(net, input, 0);

	printf("\ngeneration %d\n", ++counter);
	printf("n    input       output         error\n");
	printf("=====================================\n");

	for (int n = 0; n < net->max_neuron; ++n)
	{
		register const double a = matrix_get(input, n, 0);
		register const double b = net->layer[net->max_layer - 1].neuron[n];

		printf("%d    ", n);
		printf("%f    ", a);
		printf("%f    ", b);
		printf("%f    ", fabs(a - b));
		printf("\n");
	}

	matrix_free(input);
}

matrix *build_batch(const int max_neuron, const int batch_size)
{
	matrix *batch
		= matrix_alloc(max_neuron, batch_size, false);

	matrix_set_random(batch, false);

	return batch;
}

int main(int argc, char *argv[])
{
	const double learn_rate = 1.0;
	const int batch_size = 1000;
	const int max_layer = 3;

	FILE *input = fopen(argv[1], "r");

	matrix *a
		= matrix_read(input, tools_row_count(input), tools_col_count(input));

	printf(" %s:\n", argv[1]);
	matrix_write(a, stdout, matrix_row(a), matrix_col(a));

	fclose(input);

	struct network net_a, net_b;
	init_net(&net_a, max_layer, matrix_row(a));
	init_net(&net_b, max_layer, matrix_row(a));

	double error_a[max_layer];
	error_a[0] = 0.0;
	error_a[1] = 0.0;
	error_a[2] = 0.0;
	error_a[3] = 0.0;

	while (true)
	{
		matrix *batch = build_batch(matrix_row(a), batch_size);

		for (int n = 0; n < batch_size; ++n)
		{
			feed_forward(&net_a, batch, n);
			back_propagate(&net_a, batch, n);

			for (int l = 0; l < net_a.max_layer; ++l)
			{
				for (int j = 0; j < net_a.max_neuron; ++j)
				{
					for (int k = 0; k < net_a.max_neuron; ++k)
					{
						register double correction;

						if (l > 0)
						{
							correction
								= net_a.layer[l].error[k]*net_a.layer[l - 1].neuron[j];
						}
						else
						{
							correction
								= net_a.layer[0].error[k]*matrix_get(batch, j, n);
						}

						matrix_incr(net_a.layer[l].weight, k, j, learn_rate*correction/as_double(batch_size));
					}
				}

//				error_a[l] += sum;
			}
		}

		matrix_free(batch);

/*		for (int l = 0; l < net_a.max_layer; ++l)
		{
			matrix_incr_all(net_a.layer[l].weight,
			                learn_rate*error_a[l]/as_double(batch_size), false);
		}
*/
		test(&net_a);
	}

	return EXIT_SUCCESS;
}

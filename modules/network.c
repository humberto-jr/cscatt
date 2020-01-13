/******************************************************************************

 About
 -----


******************************************************************************/

#include "matrix.h"
#include "network.h"

struct network_layer
{
	int max_neuron;
	matrix *activation, *weight, *bias;
};

typedef struct network_layer network_layer;

struct network
{
	int max_layer;
	matrix *input;
	network_layer *layer;
};

/******************************************************************************

 Function network_alloc():

******************************************************************************/

network *network_alloc(const int max_layer)
{
	ASSERT(max_layer > 0)

	network *pointer = calloc(1, sizeof(network));
	ASSERT(pointer != NULL)

	pointer->layer = calloc(max_layer, sizeof(network_layer));
	ASSERT(pointer->layer != NULL)

	pointer->max_layer = max_layer;

	return pointer;
}

/******************************************************************************

 Function network_set_layer():

******************************************************************************/

void network_set_layer(network *net, const int l, const int max_neuron)
{
	ASSERT(l > 0)
	ASSERT(l < net->max_layer)

	ASSERT(max_neuron > 0)

	net->layer[l].max_neuron = max_neuron;
}

/******************************************************************************

 Function network_set_input():

******************************************************************************/

void network_set_input(network *net, const int size, const double input[])
{
	ASSERT(size > 0)
	net->input = matrix_alloc(size, 1, false);

	for (int n = 0; n < size; ++n)
	{
		matrix_set(net->input, n, 0, input[n]);
	}
}

/******************************************************************************

 Function network_init_random():

******************************************************************************/

void network_init_random(network *net)
{
	for (int l = 1; l < net->max_layer; ++l)
	{
		if (net->layer[l].max_neuron < 1)
		{
			PRINT_ERROR("no neurons defined for layer %d\n", net->layer[l].max_neuron)
			continue;
		}

		net->layer[l].weight
			= matrix_alloc(net->layer[l].max_neuron, net->layer[l-1].max_neuron, false);

		net->layer[l].activation
			= matrix_alloc(net->layer[l].max_neuron, 1, false);

		net->layer[l].bias
			= matrix_alloc(net->layer[l].max_neuron, 1, false);

		matrix_set_random(net->layer[l].weight, false);
		matrix_set_random(net->layer[l].bias, false);
	}
}

/******************************************************************************

 Function network_init_random():

******************************************************************************/

static void network_sigmoid(network *net, const int l)
{
	ASSERT(l >= 2)

	matrix *last_w = net->layer[l-1].weight;
	matrix *last_a = (l == 2? net->input : net->layer[l-1].activation);

	matrix *this_a = net->layer[l].activation;
	matrix *this_b = net->layer[l].bias;

	matrix_multiply(1.0, last_w, last_a, 0.0, this_a);

	for (int n = 0; n < net->layer[l].max_neuron; ++n)
	{
		register const double z
			= matrix_get(this_a, n, 0) + matrix_get(this_b, n, 0);

		matrix_set(this_a, n, 0, 1.0/(1.0 + exp(-z)));
	}

	last_w = NULL;
	last_a = NULL;
	this_a = NULL;
	this_b = NULL;
}

/******************************************************************************

 Function network_feed_forward():

******************************************************************************/

void network_feed_forward(network *net)
{
	for (int l = 2; l < net->max_layer; ++l)
	{
		network_sigmoid(net, l);
	}
}

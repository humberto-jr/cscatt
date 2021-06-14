#if !defined(VECTOR_HEADER)
	#define VECTOR_HEADER
	#include "globals.h"
	#include "matrix.h"

	typedef matrix vector;

	inline static void vector_init_gpu()
	{
		matrix_init_gpu();
	}

	inline static void vector_end_gpu()
	{
		matrix_end_gpu();
	}

	inline static vector *vector_alloc(const size_t length, const bool set_zero)
	{
		return matrix_alloc(1, length, set_zero);
	}

	inline static void vector_free(vector *v)
	{
		matrix_free(v);
	}

	inline static void vector_set(vector *v, const size_t n, const double x)
	{
		matrix_data_set(v, n, x);
	}

	inline static void vector_set_all(vector *v, const double x)
	{
		matrix_set_all(v, x);
	}

	inline static void vector_set_block(vector *v,
	                                    const size_t n_min,
	                                    const size_t n_max,
	                                    const double x)
	{
		matrix_set_block(v, 0, 0, n_min, n_max, x);
	}

	inline static void vector_set_first(vector *v, const double x)
	{
		matrix_data_set(v, 0, x);
	}

	inline static void vector_set_last(vector *v, const double x)
	{
		matrix_data_set(v, matrix_cols(x) - 1, x);
	}

	inline static void vector_set_random(vector *v)
	{
		matrix_set_random(v);
	}

	inline static void vector_set_zero(vector *v)
	{
		matrix_set_zero(v);
	}

	inline static double vector_get(const vector *v, const size_t n)
	{
		return matrix_data_get(v, n);
	}

	inline static vector *vector_get_block(const vector *v,
	                                       const size_t n_min,
	                                       const size_t n_max)
	{
		return matrix_get_block(v, 0, 0, n_min, n_max);
	}

	inline static double *vector_get_raw(const vector *v)
	{
		return matrix_get_raw_row(v, 0);
	}

	inline static size_t vector_length(const vector *v)
	{
		return matrix_cols(v);
	}

//	inline static vector_copy(vector *a,
//	                          const vector *b, const double alpha, const double beta);
//	{
//		matrix_copy(a, b, alpha, beta);
//	}

	inline static void vector_swap(vector *a, vector *b)
	{
		matrix_swap(a, b);
	}

	inline static double vector_sum(const vector *v)
	{
		return matrix_sum(v);
	}

	inline static double vector_min(const vector *v)
	{
		return matrix_min(v);
	}

	inline static double vector_max(const vector *v)
	{
		return matrix_max(v);
	}

//	void matrix_multiply(const double alpha, const matrix *a,
//	                     const matrix *b, const double beta, matrix *c);

	inline static void vector_add(const double alpha, const vector *a,
	                              const double beta, const vector *b, vector *c)
	{
		matrix_add(alpha, a, beta, b, c);
	}

	inline static void vector_sub(const double alpha, const vector *a,
	                              const double beta, const vector *b, vector *c)
	{
		matrix_sub(alpha, a, beta, b, c);
	}

	inline static bool vector_is_null(const vector *v)
	{
		return matrix_is_null(v);
	}

	inline static bool vector_is_positive(const vector *v)
	{
		return matrix_is_positive(v);
	}

	inline static bool vector_is_negative(const vector *v)
	{
		return matrix_is_negative(v);
	}

	inline static bool vector_has_nan(const vector *v)
	{
		return matrix_has_nan(v);
	}

	inline static bool vector_using_magma()
	{
		return matrix_using_magma();
	}

	inline static bool vector_using_mkl()
	{
		return matrix_using_mkl();
	}

	inline static bool vector_using_lapacke()
	{
		return vector_using_lapacke();
	}

	inline static void vector_save(const vector *v, const char filename[])
	{
		matrix_save(v, filename);
	}

	inline static vector *vector_load(const char filename[])
	{
		return matrix_load(filename);
	}

	inline static size_t vector_sizeof(const vector *v)
	{
		return matrix_sizeof(v);
	}

	inline static void vector_use_omp(vector *v, const bool use)
	{
		matrix_use_omp(v, use);
	}

	inline static void vector_write(const vector *v, FILE *output)
	{
		ASSERT(output != NULL)

		const size_t n_max = matrix_cols(v);

		for (size_t n = 0; n < n_max; ++n)
			fprintf(output, "% -8e\n", matrix_get(v, 0, n));
	}

	inline static void vector_resize(vector *v, const size_t length)
	{
		matrix_reshape(v, 1, length);
	}

	inline static void vector_enlarge(vector *v, const size_t extra_length)
	{
		matrix_reshape(v, 1, matrix_cols(v) + extra_length);
	}

	inline static void vector_about(FILE *output)
	{
		ASSERT(output != NULL)

		fprintf(output, "# build date   = %s\n", __DATE__);
		fprintf(output, "# source code  = %s\n", __FILE__);

		if (matrix_using_magma())
			fprintf(output, "# lin. algebra = MAGMA\n");

		else if (matrix_using_lapacke())
			fprintf(output, "# lin. algebra = LAPACKE + GSL CBLAS\n");

		else if (matrix_using_mkl())
			fprintf(output, "# lin. algebra = Intel Math Kernel Library (MKL)\n");

		else
			fprintf(output, "# lin. algebra = GNU Scientific Library (GSL)\n");
	}
#endif

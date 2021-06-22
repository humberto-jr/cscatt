#if !defined(MATRIX_HEADER)
	#define MATRIX_HEADER
	#include "globals.h"

	typedef struct matrix matrix;

	struct tensor
	{
		matrix *value;
	};

	typedef struct tensor tensor;

	void matrix_init_gpu();

	void matrix_end_gpu();

	matrix *matrix_alloc(const size_t max_row,
	                     const size_t max_col, const bool set_zero);

	matrix *matrix_alloc_as(const matrix *m, const bool set_zero);

	void matrix_free(matrix *m);

	void matrix_set(matrix *m,
	                const size_t p, const size_t q, const double x);

	void matrix_set_all(matrix *m, const double x);

	void matrix_set_diag(matrix *m, const size_t p, const double x);

	void matrix_set_symm(matrix *m,
	                     const size_t p, const size_t q, const double x);

	void matrix_set_row(matrix *m, const size_t p, const double x);

	void matrix_set_col(matrix *m, const size_t q, const double x);

	void matrix_set_block(matrix *m,
	                      const size_t row_min,
	                      const size_t row_max,
	                      const size_t col_min,
	                      const size_t col_max,
	                      const double x);

	void matrix_set_random(matrix *m);

	void matrix_set_zero(matrix *m);

	double matrix_get(const matrix *m, const size_t p, const size_t q);

	matrix *matrix_get_row(const matrix *m, const size_t p);

	matrix *matrix_get_col(const matrix *m, const size_t q);

	matrix *matrix_get_diag(const matrix *m);

	matrix *matrix_get_block(const matrix *m,
	                         const size_t row_min,
	                         const size_t row_max,
	                         const size_t col_min,
	                         const size_t col_max);

	double *matrix_get_raw_row(const matrix *m, const size_t p);

	double *matrix_get_raw_col(const matrix *m, const size_t q);

	size_t matrix_rows(const matrix *m);

	size_t matrix_cols(const matrix *m);

	void matrix_incr(matrix *m,
	                 const size_t p, const size_t q, const double x);

	void matrix_incr_all(matrix *m, const double x);

	void matrix_decr(matrix *m,
	                 const size_t p, const size_t q, const double x);

	void matrix_decr_all(matrix *m, const double x);

	void matrix_scale(matrix *m,
	                  const size_t p, const size_t q, const double x);

	void matrix_scale_all(matrix *m, const double x);

	void matrix_scale_row(matrix *m, const size_t p, const double x);

	void matrix_scale_col(matrix *m, const size_t q, const double x);

	void matrix_copy_element(matrix *a,
	                         const size_t p,
	                         const size_t q,
	                         const matrix *b,
	                         const size_t l,
	                         const size_t k);

	void matrix_copy(matrix *a,
	                 const matrix *b, const double alpha, const double beta);

	void matrix_swap(matrix *a, matrix *b);

	double matrix_trace(const matrix *m);

	double matrix_sum(const matrix *m);

	double matrix_sum_row(const matrix *m, const size_t p);

	double matrix_sum_col(const matrix *m, const size_t q);

	double matrix_min(const matrix *m);

	double matrix_max(const matrix *m);

	void matrix_multiply(const double alpha, const matrix *a,
	                     const matrix *b, const double beta, matrix *c);

	void matrix_add(const double alpha, const matrix *a,
	                const double beta, const matrix *b, matrix *c);

	void matrix_sub(const double alpha, const matrix *a,
	                const double beta, const matrix *b, matrix *c);

	void matrix_inverse(matrix *m);

	double *matrix_symm_eigen(matrix *m, const char job);

	bool matrix_is_null(const matrix *m);

	bool matrix_is_positive(const matrix *m);

	bool matrix_is_negative(const matrix *m);

	bool matrix_is_square(const matrix *m);

	bool matrix_has_nan(const matrix *m);

	bool matrix_using_magma();

	bool matrix_using_mkl();

	bool matrix_using_lapacke();

	void matrix_save(const matrix *m, const char filename[]);

	matrix *matrix_load(const char filename[]);

	matrix *matrix_read(FILE *input,
	                    const size_t max_row, const size_t max_col);

	void matrix_write(const matrix *m,
	                  FILE *output, const size_t max_row, const size_t max_col);

	size_t matrix_sizeof(const matrix *m);

	void matrix_use_omp(matrix *m, const bool use);

	void matrix_reshape(matrix *m, const size_t max_row, const size_t max_col);

	void matrix_data_set(matrix *m, const size_t n, const double x);

	double matrix_data_get(const matrix *m, const size_t n);

	double *matrix_data_raw(const matrix *m);

	size_t matrix_data_length(const matrix *m);

	void matrix_about(FILE *output);
#endif

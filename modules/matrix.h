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

	matrix *matrix_alloc(const int max_row,
	                     const int max_col, const bool set_zero);

	void matrix_free(matrix *m);

	void matrix_set(matrix *m, const int p, const int q, const double x);

	void matrix_set_all(matrix *m, const double x, const bool use_omp);

	void matrix_diag_set(matrix *m, const int p, const double x);

	void matrix_diag_set_all(matrix *m, const double x, const bool use_omp);

	void matrix_symm_set(matrix *m, const int p, const int q, const double x);

	double matrix_get(const matrix *m, const int p, const int q);

	matrix *matrix_get_row(const matrix *m, const int p, const bool use_omp);

	matrix *matrix_get_col(const matrix *m, const int q, const bool use_omp);

	matrix *matrix_get_diag(const matrix *m, const bool use_omp);

	matrix *matrix_get_block(const matrix *m,
	                         const int row_min, const int row_max,
	                         const int col_min, const int col_max);

	double *matrix_raw_row(const matrix *m, const int p, const bool use_omp);

	double *matrix_raw_col(const matrix *m, const int q, const bool use_omp);

	double matrix_get_pow(const matrix *m,
	                      const int p, const int q, const double power);

	int matrix_row(const matrix *m);

	int matrix_col(const matrix *m);

	void matrix_set_random(matrix *m, const bool use_omp);

	void matrix_incr(matrix *m, const int p, const int q, const double x);

	void matrix_incr_all(matrix *m, const double x, const bool use_omp);

	void matrix_decr(matrix *m, const int p, const int q, const double x);

	void matrix_decr_all(matrix *m, const double x, const bool use_omp);

	void matrix_scale(matrix *m, const int p, const int q, const double x);

	void matrix_scale_all(matrix *m, const double x, const bool use_omp);

	void matrix_row_scale(matrix *m,
	                      const int p, const double x, const bool use_omp);

	void matrix_col_scale(matrix *m,
	                      const int q, const double x, const bool use_omp);

	void matrix_copy_element(matrix *a, const int p, const int q,
	                         const matrix *b, const int l, const int k);

	void matrix_copy(matrix *a, const matrix *b, const double alpha,
	                 const double beta, const bool use_omp);

	void matrix_swap(matrix *a, matrix *b);

	double matrix_trace(const matrix *m, const bool use_omp);

	double matrix_sum(const matrix *m, const bool use_omp);

	double matrix_row_sum(const matrix *m, const int p, const bool use_omp);

	double matrix_col_sum(const matrix *m, const int q, const bool use_omp);

	double matrix_min(const matrix *m);

	double matrix_max(const matrix *m);

	void matrix_multiply(const double alpha, const matrix *a,
	                     const matrix *b, const double beta, matrix *c);

	void matrix_add(const double alpha, const matrix *a, const double beta,
	                const matrix *b, matrix *c, const bool use_omp);

	void matrix_sub(const double alpha, const matrix *a, const double beta,
	                const matrix *b, matrix *c, const bool use_omp);

	void matrix_inverse(matrix *m);

	double *matrix_symm_eigen(matrix *m, const char job);

	double matrix_row_quadr(const matrix *m, const int p, const bool use_omp);

	double matrix_col_quadr(const matrix *m, const int q, const bool use_omp);

	bool matrix_null(const matrix *m);

	bool matrix_positive(const matrix *m);

	bool matrix_negative(const matrix *m);

	bool matrix_has_nan(const matrix *m);

	void matrix_save(const matrix *m, const char filename[]);

	matrix *matrix_load(const char filename[]);

	matrix *matrix_read(FILE *input, const int max_row, const int max_col);

	void matrix_write(const matrix *m,
	                  FILE *output, const int max_row, const int max_col);

	size_t matrix_sizeof(const matrix *m);

	void matrix_about(FILE *output);
#endif

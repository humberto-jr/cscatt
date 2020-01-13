/******************************************************************************

 About
 -----

 This module defines the opaque type matrix and many general purpose functions
 needed to manipulate it. Interfaces for few popular linear algebra libraries
 are builtin (including the use of GPUs) and are tuned during compilation.

******************************************************************************/

#if defined(USE_MKL)
	#include <mkl.h>
	#define LINEAR_ALGEBRA_LIB "Intel Math Kernel Library (MKL)"
#endif

#if defined(USE_LAPACKE)
	#include <lapacke.h>
	#include <gsl/gsl_cblas.h>
	#define LINEAR_ALGEBRA_LIB "LAPACKE + GSL CBLAS"

	/* NOTE: typedef for compatibility with MKL and LAPACKE. */
	typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
#endif

#if defined(USE_MAGMA)
	#include <magma_v2.h>
	#include <magma_lapack.h>
	#define LINEAR_ALGEBRA_LIB "MAGMA"

	/* NOTE: matrix_init_gpu() will setup the queue. */
	static magma_queue_t gpu_queue;
#endif

#if !defined(LINEAR_ALGEBRA_LIB)
	#include <gsl/gsl_cblas.h>
	#include <gsl/gsl_eigen.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_matrix.h>
	#define LINEAR_ALGEBRA_LIB "GNU Scientific Library (GSL)"

	/* NOTE: typedef for compatibility with MKL and LAPACKE. */
	typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
#endif

#include "matrix.h"

#if !defined(LAPACK_DISABLE_NAN_CHECK)
	#define LAPACK_DISABLE_NAN_CHECK
#endif

/******************************************************************************

 Macro DATA_OFFSET(): matrices are always ordinary vectors with elements stored
 in a row-major scheme: vector[p*max_col + q], where p = [0, max_row) and q =
 [0, max_col). Thus, DATA_OFFSET expands to the actual data offset expression
 for a given pq-element pointed by a pointer.

******************************************************************************/

struct matrix
{
	double *data;
	int max_row, max_col;
};

#define DATA_OFFSET(name, p, q) name->data[(p)*(name->max_col) + (q)]

/******************************************************************************

 Macro ASSERT_ROW_INDEX(): check if the p-th element is within the row bounds
 of a matrix pointed by a given pointer.

******************************************************************************/

#define ASSERT_ROW_INDEX(pointer, p) \
{                                    \
	ASSERT((p) > -1)                   \
	ASSERT((p) < pointer->max_row)     \
}

/******************************************************************************

 Macro ASSERT_COL_INDEX(): check if the q-th element is within the column bounds
 of a matrix pointed by a given pointer.

******************************************************************************/

#define ASSERT_COL_INDEX(pointer, q) \
{                                    \
	ASSERT((q) > -1)                   \
	ASSERT((q) < pointer->max_col)     \
}

/******************************************************************************

 Function min(): return the min between two integers a and b.

******************************************************************************/

inline static int min(const int a, const int b)
{
	return (a < b? a : b);
}

/******************************************************************************

 Wrapper call_dgemm(): a general call to dgemm interfacing the many different
 possible libraries. For the meaning of each input parameter, please refer to
 dgemm documentation in the netlib repository:

 http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html

 NOTE: Row-major matrices only.

******************************************************************************/

void call_dgemm(const char trans_a,
                const char trans_b,
                const int m,
                const int n,
                const int k,
                const double alpha,
                const double a[],
                const int lda,
                const double b[],
                const int ldb,
                const double beta,
                double c[],
                const int ldc)
{
	#if defined(USE_MAGMA)
	{
/*
 *		NOTE: a normal matrix in the host (row-major) implies transposed matrix
 *		in the device (col-major).
 */
		const magma_trans_t a_form
			= (trans_a == 't'? MagmaNoTrans : MagmaTrans);

		const magma_trans_t b_form
			= (trans_b == 't'? MagmaNoTrans : MagmaTrans);

		double *a_gpu = NULL, *b_gpu = NULL, *c_gpu = NULL;

		magma_dmalloc(&a_gpu, m*k);
		magma_dmalloc(&b_gpu, k*n);
		magma_dmalloc(&c_gpu, m*n);

		magma_setmatrix_async(m, k, sizeof(double), a, lda, a_gpu, lda, gpu_queue);
		magma_setmatrix_async(k, n, sizeof(double), b, ldb, b_gpu, ldb, gpu_queue);

		if (beta != 0.0)
		{
			magma_setmatrix_async(m, n, sizeof(double), c, ldc, c_gpu, ldc, gpu_queue);
		}

		magma_queue_sync(gpu_queue);

		magmablas_dgemm(a_form, b_form, m, n, k,
		                alpha, a_gpu, lda, b_gpu, ldb, beta, c_gpu, ldc, gpu_queue);

		magma_getmatrix(m, n, sizeof(double), c_gpu, ldc, c, ldc, gpu_queue);

		magma_free(a_gpu);
		magma_free(b_gpu);
		magma_free(c_gpu);
	}
	#else
	{
		const CBLAS_TRANSPOSE a_form
			= (trans_a == 't'? CblasTrans : CblasNoTrans);

		const CBLAS_TRANSPOSE b_form
			= (trans_b == 't'? CblasTrans : CblasNoTrans);

		cblas_dgemm(CblasRowMajor, a_form, b_form,
		            m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	}
	#endif
}

/******************************************************************************

 Wrapper call_dsyev(): a general call to dsyev interfacing the many different
 possible libraries. For the meaning of each input parameter, please refer to
 dsyev documentation in the netlib repository:

 http://www.netlib.org/lapack/explore-html/dd/d4c/dsyev_8f.html

 NOTE: when GSL is used, only uplo = 'l' is available.

******************************************************************************/

void call_dsyev(const char jobz,
                const char uplo,
                const int n,
                double a[],
                const int lda,
                double w[])
{
	ASSERT(a != NULL)
	ASSERT(w != NULL)

	#if defined(USE_MAGMA)
	{
		const magma_vec_t job_mode
			= (jobz == 'n'? MagmaNoVec : MagmaVec);

		const magma_uplo_t fill_mode
			= (uplo == 'u'? MagmaUpper : MagmaLower);

		const magma_int_t nb
			= magma_get_dsytrd_nb((magma_int_t) n);

		const magma_int_t lwork
			= (jobz == 'n'? 2*n + n*nb : GSL_MAX(2*n + n*nb, 1 + 6*n + 2*n*n));

		const magma_int_t liwork
			= (jobz == 'n'? n : 3 + 5*n);

		double *work = NULL;
		magma_dmalloc_cpu(&work, lwork);

		magma_int_t *iwork = NULL;
		magma_imalloc_cpu(&iwork, liwork);

		magma_int_t info = 1;

		magma_dsyevd_m(1, job_mode, fill_mode,
		               n, a, lda, w, work, lwork, iwork, liwork, &info);

		if (info != 0)
		{
			PRINT_ERROR("magma_dsyevd_m() failed with error code %d\n", info)
			exit(EXIT_FAILURE);
		}

		magma_free_cpu(work);
		magma_free_cpu(iwork);
	}
	#elif defined(USE_MKL) || defined(USE_LAPACKE)
	{
		const int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR,
		                      jobz, uplo, n, a, lda, w);

		if (info != 0)
		{
			PRINT_ERROR("LAPACKE_dsyev() failed with error code %d\n", info)
			exit(EXIT_FAILURE);
		}
	}
	#else
	{
		ASSERT(n > 0)
		ASSERT(lda <= n)
		ASSERT(uplo == 'l')
		ASSERT(jobz == 'n' || jobz == 'v')

		gsl_matrix_view A = gsl_matrix_view_array(a, n, n);
		gsl_vector_view W = gsl_vector_view_array(w, n);

		gsl_matrix *eigenvec = gsl_matrix_alloc(n, n);

		gsl_eigen_symmv_workspace *work = gsl_eigen_symmv_alloc(n);
		gsl_eigen_symmv(&A.matrix, &W.vector, eigenvec, work);
		gsl_eigen_symmv_free(work);

		gsl_eigen_symmv_sort(&W.vector, eigenvec, GSL_EIGEN_SORT_VAL_ASC);
		if (jobz == 'v') gsl_matrix_memcpy(&A.matrix, eigenvec);

		gsl_matrix_free(eigenvec);
	}
	#endif
}

/******************************************************************************

 Function matrix_init_gpu(): allocate resources needed by MAGMA and shall be
 called before any call to functions of this module when the macro USE_MAGMA
 is defined during the compilation.

 NOTE: this function will setup properties shared by all matrix objects.

******************************************************************************/

void matrix_init_gpu()
{
	#if defined(USE_MAGMA)
		if (magma_init() != MAGMA_SUCCESS)
		{
			PRINT_ERROR("%s\n", "magma_init() failed")
			exit(EXIT_FAILURE);
		}

		magma_queue_create(0, &gpu_queue);
	#else
		PRINT_ERROR("%s\n", "macro USE_MAGMA not defined")
	#endif
}

/******************************************************************************

 Function matrix_end_gpu(): free resources allocated by matrix_init_gpu().

******************************************************************************/

void matrix_end_gpu()
{
	#if defined(USE_MAGMA)
		if (magma_finalize() != MAGMA_SUCCESS)
		{
			PRINT_ERROR("%s\n", "magma_finalize() failed")
		}

		magma_queue_destroy(gpu_queue);
	#else
		PRINT_ERROR("%s\n", "macro USE_MAGMA not defined")
	#endif
}

/******************************************************************************

 Function matrix_alloc(): allocate resources for a matrix of shape max_row-by
 -max_col.

 NOTE: this function shall return only if all input parameters are good and if
 resources are allocated properly. Thus, neither the caller or other functions
 from the module are required to (re)check returned pointers.

******************************************************************************/

matrix *matrix_alloc(const int max_row, const int max_col, const bool set_zero)
{
	ASSERT(max_row > 0)
	ASSERT(max_col > 0)

	matrix *pointer = calloc(1, sizeof(matrix));
	ASSERT(pointer != NULL)

	pointer->max_row = max_row;
	pointer->max_col = max_col;

	#if defined(USE_MAGMA)
		magma_dmalloc_pinned(pointer->data, max_row*max_col);
		if (set_zero) matrix_set_all(pointer, 0.0, false);
	#else
		pointer->data = allocate(max_row*max_col, sizeof(double), set_zero);
	#endif

	return pointer;
}

/******************************************************************************

 Function matrix_free(): release resources allocated by matrix_alloc().

******************************************************************************/

void matrix_free(matrix *m)
{
	#if defined(USE_MAGMA)
		magma_free_pinned(m->data);
	#else
		free(m->data);
	#endif

	free(m);
}

/******************************************************************************

 Function matrix_set(): set x to the pq-element of matrix m.

******************************************************************************/

void matrix_set(matrix *m, const int p, const int q, const double x)
{
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)
	DATA_OFFSET(m, p, q) = x;
}

/******************************************************************************

 Function matrix_set_all(): the same as matrix_set() but for all elements.

******************************************************************************/

void matrix_set_all(matrix *m, const double x, const bool use_omp)
{
	const int n_max = m->max_row*m->max_col;

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		m->data[n] = x;
	}
}

/******************************************************************************

 Function matrix_diag_set(): set x to the p-th diagonal element of matrix m.

******************************************************************************/

void matrix_diag_set(matrix *m, const int p, const double x)
{
	ASSERT_ROW_INDEX(m, p)
	DATA_OFFSET(m, p, p) = x;
}

/******************************************************************************

 Function matrix_diag_set_all(): the same as matrix_diag_set() but for all
 diagonal elements.

******************************************************************************/

void matrix_diag_set_all(matrix *m, const double x, const bool use_omp)
{
	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int n = 0; n < m->max_row; ++n)
	{
		DATA_OFFSET(m, n, n) = x;
	}
}

/******************************************************************************

 Function matrix_symm_set(): set x to both the pq-element and qp-element of a
 symmetric matrix m.

******************************************************************************/

void matrix_symm_set(matrix *m, const int p, const int q, const double x)
{
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)

	DATA_OFFSET(m, p, q) = x;
	DATA_OFFSET(m, q, p) = x;
}

/******************************************************************************

 Function matrix_get(): return the pq-element of matrix m.

******************************************************************************/

double matrix_get(const matrix *m, const int p, const int q)
{
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)
	return DATA_OFFSET(m, p, q);
}

/******************************************************************************

 Function matrix_get_row(): the same as matrix_get() but return the p-th row of
 matrix m.

 NOTE: the return matrix is a row-matrix.

******************************************************************************/

matrix *matrix_get_row(const matrix *m, const int p, const bool use_omp)
{
	ASSERT_ROW_INDEX(m, p)

	matrix *row = matrix_alloc(1, m->max_col, false);

	#pragma omp parallel for default(none) shared(m, row) schedule(static) if(use_omp)
	for (int q = 0; q < m->max_col; ++q)
	{
		DATA_OFFSET(row, 0, q) = DATA_OFFSET(m, p, q);
	}

	return row;
}

/******************************************************************************

 Function matrix_get_col(): the same as matrix_get() but return the q-th column
 of matrix m.

 NOTE: the return matrix is a column-matrix.

******************************************************************************/

matrix *matrix_get_col(const matrix *m, const int q, const bool use_omp)
{
	ASSERT_COL_INDEX(m, q)

	matrix *col = matrix_alloc(m->max_row, 1, false);

	#pragma omp parallel for default(none) shared(m, col) schedule(static) if(use_omp)
	for (int p = 0; p < m->max_row; ++p)
	{
		DATA_OFFSET(col, p, 0) = DATA_OFFSET(m, p, q);
	}

	return col;
}

/******************************************************************************

 Function matrix_get_diag(): the same as matrix_get() but return the diagonal
 of matrix m.

 NOTE: the return matrix is a column-matrix.

******************************************************************************/

matrix *matrix_get_diag(const matrix *m, const bool use_omp)
{
	matrix *diag = matrix_alloc(m->max_row, 1, false);

	#pragma omp parallel for default(none) shared(m, diag) schedule(static) if(use_omp)
	for (int p = 0; p < m->max_row; ++p)
	{
		DATA_OFFSET(diag, p, 0) = DATA_OFFSET(m, p, p);
	}

	return diag;
}

/******************************************************************************

 Function matrix_get_block(): the same as matrix_get() but return a block of
 matrix m.

******************************************************************************/

matrix *matrix_get_block(const matrix *m,
                         const int row_min,
                         const int row_max,
                         const int col_min,
                         const int col_max)
{
	ASSERT(row_max >= row_min)
	ASSERT(col_max >= col_min)

	ASSERT_ROW_INDEX(m, row_min)
	ASSERT_ROW_INDEX(m, row_max)

	ASSERT_COL_INDEX(m, col_min)
	ASSERT_COL_INDEX(m, col_max)

	matrix *block
		= matrix_alloc(row_max - row_min + 1, col_max - col_min + 1, false);

	int row = 0;
	for (int p = row_min; p <= row_max; ++p)
	{
		int col = 0;
		for (int q = col_min; q <= col_max; ++q)
		{
			DATA_OFFSET(block, row, col) = DATA_OFFSET(m, p, q);
			++col;
		}

		++row;
	}

	return block;
}

/******************************************************************************

 Function matrix_row(): return the number of rows in the matrix m.

******************************************************************************/

int matrix_row(const matrix *m)
{
	return m->max_row;
}

/******************************************************************************

 Function matrix_col(): return the number of columns in the matrix m.

******************************************************************************/

int matrix_col(const matrix *m)
{
	return m->max_col;
}

/******************************************************************************

 Function matrix_set_random(): the same as matrix_set_all() but elements are
 defined randomly within [0, 1].

******************************************************************************/

void matrix_set_random(matrix *m, const bool use_omp)
{
	srand(rand());
	const int n_max = m->max_row*m->max_col;

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		m->data[n] = (double) rand()/RAND_MAX;
	}
}

/******************************************************************************

 Function matrix_incr(): increment the pq-element of matrix m by x (+=).

******************************************************************************/

void matrix_incr(matrix *m, const int p, const int q, const double x)
{
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)
	DATA_OFFSET(m, p, q) += x;
}

/******************************************************************************

 Function matrix_incr_all(): the same as matrix_incr() but for all elements.

******************************************************************************/

void matrix_incr_all(matrix *m, const double x, const bool use_omp)
{
	const int n_max = m->max_row*m->max_col;

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		m->data[n] += x;
	}
}

/******************************************************************************

 Function matrix_decr(): decrement the pq-element of matrix m by x (-=).

******************************************************************************/

void matrix_decr(matrix *m, const int p, const int q, const double x)
{
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)
	DATA_OFFSET(m, p, q) -= x;
}

/******************************************************************************

 Function matrix_decr_all(): the same as matrix_decr() but for all elements.

******************************************************************************/

void matrix_decr_all(matrix *m, const double x, const bool use_omp)
{
	const int n_max = m->max_row*m->max_col;

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		m->data[n] -= x;
	}
}

/******************************************************************************

 Function matrix_scale(): scale the pq-element of matrix m by x (*=).

******************************************************************************/

void matrix_scale(matrix *m, const int p, const int q, const double x)
{
	ASSERT_ROW_INDEX(m, p)
	ASSERT_COL_INDEX(m, q)
	DATA_OFFSET(m, p, q) *= x;
}

/******************************************************************************

 Function matrix_scale_all(): the same as matrix_scale() but for all elements.

******************************************************************************/

void matrix_scale_all(matrix *m, const double x, const bool use_omp)
{
	const int n_max = m->max_row*m->max_col;

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		m->data[n] *= x;
	}
}

/******************************************************************************

 Function matrix_row_scale(): the same as matrix_scale() but for the p-th row
 of matrix m.

******************************************************************************/

void matrix_row_scale(matrix *m,
                      const int p, const double x, const bool use_omp)
{
	ASSERT_ROW_INDEX(m, p)

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int q = 0; q < m->max_col; ++q)
	{
		DATA_OFFSET(m, p, q) *= x;
	}
}

/******************************************************************************

 Function matrix_col_scale(): the same as matrix_scale() but for the q-th
 column of matrix m.

******************************************************************************/

void matrix_col_scale(matrix *m,
                      const int q, const double x, const bool use_omp)
{
	ASSERT_COL_INDEX(m, q)

	#pragma omp parallel for default(none) shared(m) schedule(static) if(use_omp)
	for (int p = 0; p < m->max_row; ++p)
	{
		DATA_OFFSET(m, p, q) *= x;
	}
}

/******************************************************************************

 Function matrix_copy_element(): copy the lk-element of matrix b to the
 pq-element of matrix a.

******************************************************************************/

void matrix_copy_element(matrix *a, const int p, const int q,
                         const matrix *b, const int l, const int k)
{
	ASSERT_ROW_INDEX(a, p)
	ASSERT_COL_INDEX(a, q)

	ASSERT_ROW_INDEX(b, l)
	ASSERT_COL_INDEX(b, k)

	DATA_OFFSET(a, p, q) = DATA_OFFSET(b, l, k);
}

/******************************************************************************

 Function matrix_copy(): copy all elements from b to a, a = b*alpha + beta.

******************************************************************************/

void matrix_copy(matrix *a, const matrix *b,
                 const double alpha, const double beta, const bool use_omp)
{
	const int n_max
		= min(a->max_row*a->max_col, b->max_row*b->max_col);

	#pragma omp parallel for default(none) shared(a, b) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		a->data[n] = b->data[n]*alpha + beta;
	}
}

/******************************************************************************

 Function matrix_swap(): swap the shape and elements of matrices a and b.

******************************************************************************/

void matrix_swap(matrix *a, matrix *b)
{
	const int b_row = b->max_row;
	const int b_col = b->max_col;
	double *b_data = b->data;

	b->max_row = a->max_row;
	b->max_col = a->max_col;
	b->data = a->data;

	a->max_row = b_row;
	a->max_col = b_col;
	a->data = b_data;

	b_data = NULL;
}

/******************************************************************************

 Function matrix_trace(): return the trace of (a square) matrix m, tr(m).

******************************************************************************/

double matrix_trace(const matrix *m, const bool use_omp)
{
	double sum = 0.0;

	#pragma omp parallel for default(none) shared(m) reduction(+:sum) if(use_omp)
	for (int n = 0; n < m->max_row; ++n)
	{
		sum += DATA_OFFSET(m, n, n);
	}

	return sum;
}

/******************************************************************************

 Function matrix_sum(): return the sum of all elements of matrix m.

******************************************************************************/

double matrix_sum(const matrix *m, const bool use_omp)
{
	double sum = 0.0;
	const int n_max = m->max_row*m->max_col;

	#pragma omp parallel for default(none) shared(m) reduction(+:sum) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		sum += m->data[n];
	}

	return sum;
}

/******************************************************************************

 Function matrix_row_sum(): return the sum of all elements in the row p of
 matrix m.

******************************************************************************/

double matrix_row_sum(const matrix *m, const int p, const bool use_omp)
{
	ASSERT_ROW_INDEX(m, p)

	double sum = 0.0;

	#pragma omp parallel for default(none) shared(m) reduction(+:sum) if(use_omp)
	for (int q = 0; q < m->max_col; ++q)
	{
		sum += DATA_OFFSET(m, p, q);
	}

	return sum;
}

/******************************************************************************

 Function matrix_col_sum(): return the sum of all elements in the column q of
 matrix m.

******************************************************************************/

double matrix_col_sum(const matrix *m, const int q, const bool use_omp)
{
	ASSERT_COL_INDEX(m, q)

	double sum = 0.0;

	#pragma omp parallel for default(none) shared(m) reduction(+:sum) if(use_omp)
	for (int p = 0; p < m->max_row; ++p)
	{
		sum += DATA_OFFSET(m, p, q);
	}

	return sum;
}

/******************************************************************************

 Function matrix_min(): return the smallest element of matrix m.

******************************************************************************/

double matrix_min(const matrix *m)
{
	double min = INF;
	const int n_max = m->max_row*m->max_col;

	for (int n = 0; n < n_max; ++n)
	{
		if (m->data[n] < min) min = m->data[n];
	}

	return min;
}

/******************************************************************************

 Function matrix_max(): return the biggest element of matrix m.

******************************************************************************/

double matrix_max(const matrix *m)
{
	double max = -INF;
	const int n_max = m->max_row*m->max_col;

	for (int n = 0; n < n_max; ++n)
	{
		if (m->data[n] > max) max = m->data[n];
	}

	return max;
}

/******************************************************************************

 Function matrix_multiply(): perform the operation c = alpha*a*b + beta*c.

******************************************************************************/

void matrix_multiply(const double alpha,
                     const matrix *a, const matrix *b, const double beta, matrix *c)
{
	call_dgemm('n', 'n', a->max_row, b->max_col, a->max_col, alpha, a->data,
	            a->max_row, b->data, a->max_col, beta, c->data, c->max_row);
}

/******************************************************************************

 Function matrix_add(): perform the operation c = alpha*a + beta*b.

******************************************************************************/

void matrix_add(const double alpha, const matrix *a, const double beta,
                const matrix *b, matrix *c, const bool use_omp)
{
	const int n_max
		= min(a->max_row*a->max_col, b->max_row*b->max_col);

	#pragma omp parallel for default(none) shared(a, b, c) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		c->data[n] = a->data[n]*alpha + b->data[n]*beta;
	}
}

/******************************************************************************

 Function matrix_sub(): perform the operation c = alpha*a - beta*b.

******************************************************************************/

void matrix_sub(const double alpha, const matrix *a, const double beta,
                const matrix *b, matrix *c, const bool use_omp)
{
	const int n_max
		= min(a->max_row*a->max_col, b->max_row*b->max_col);

	#pragma omp parallel for default(none) shared(a, b, c) schedule(static) if(use_omp)
	for (int n = 0; n < n_max; ++n)
	{
		c->data[n] = a->data[n]*alpha - b->data[n]*beta;
	}
}

/******************************************************************************

 Function matrix_inverse(): to invert the matrix m.

******************************************************************************/

void matrix_inverse(matrix *m)
{
	#if defined(USE_MAGMA)
	{
		magma_int_t *ipiv = NULL;
		magma_imalloc_cpu(&ipiv, min(m->max_row, m->max_col));

		double *m_gpu = NULL;
		magma_dmalloc(&m_gpu, m->max_row*m->max_col);

		magma_setmatrix(m->max_row, m->max_col, sizeof(double),
		                m->data, m->max_row, m_gpu, m->max_row, gpu_queue);

		magma_int_t info = 1;
		magma_dgetrf_gpu(m->max_row, m->max_col, m_gpu, m->max_row, ipiv, &info);

		if (info != 0)
		{
			PRINT_ERROR("magma_dgetrf_gpu() failed with error code %d\n", info)
			exit(EXIT_FAILURE);
		}

		const magma_int_t lwork
			= m->max_row*magma_get_dgetri_nb(m->max_row);

		double *work = NULL;
		magma_dmalloc(&work, lwork);

		info = 1;
		magma_dgetri_gpu(m->max_row, m_gpu, m->max_row, ipiv, work, lwork, &info);

		if (info != 0)
		{
			PRINT_ERROR("magma_dgetri_gpu() failed with error code %d\n", info)
			exit(EXIT_FAILURE);
		}

		magma_getmatrix(m->max_row, m->max_col, sizeof(double),
		                m_gpu, m->max_row, m->data, m->max_row, gpu_queue);

		magma_free(ipiv);
		magma_free(work);
		magma_free(m_gpu);
	}
	#elif defined(USE_MKL) || defined(USE_LAPACKE)
	{
		#if defined(USE_LAPACKE)
			int *ipiv = malloc(sizeof(int)*m->max_row);
		#else
			long long int *ipiv = malloc(sizeof(long long int)*m->max_row);
		#endif

		ASSERT(ipiv != NULL)

		int info = LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'u', m->max_row,
		                                 m->data, m->max_row, ipiv);

		if (info != 0)
		{
			PRINT_ERROR("LAPACKE_dsytrf() failed with error code %d\n", info)
			exit(EXIT_FAILURE);
		}

		info = LAPACKE_dsytri(LAPACK_ROW_MAJOR, 'u', m->max_row,
		                             m->data, m->max_row, ipiv);

		if (info != 0)
		{
			PRINT_ERROR("LAPACKE_dsytri() failed with error code %d\n", info)
			exit(EXIT_FAILURE);
		}

/*
 *		NOTE: To this point, only the upper triangular part is actually
 *		inverted. Thus, we have to copy lapacke result into the lower one.
 */

		for (int p = 0; p < m->max_row; ++p)
		{
			for (int q = (p + 1); q < m->max_col; ++q)
			{
				DATA_OFFSET(m, q, p) = DATA_OFFSET(m, p, q);
			}
		}

		free(ipiv);
	}
	#else
	{
		double *buffer = allocate(m->max_row*m->max_col, sizeof(double), false);

		for (int n = 0; n < m->max_row*m->max_col; ++n)
		{
			buffer[n] = m->data[n];
		}

		gsl_permutation *p = gsl_permutation_alloc(m->max_row);
		gsl_matrix_view  a = gsl_matrix_view_array(buffer, m->max_row, m->max_col);
		gsl_matrix_view  b = gsl_matrix_view_array(m->data, m->max_row, m->max_col);

		int sign = 0;
		gsl_linalg_LU_decomp(&a.matrix, p, &sign);
		gsl_linalg_LU_invert(&a.matrix, p, &b.matrix);

		gsl_permutation_free(p);
		free(buffer);
	}
	#endif
}

/******************************************************************************

 Function matrix_symm_eigen(): return the eigenvalues of a symmetric matrix. On
 exit, the original matrix is replaced by the respective eigenvectors if job =
 'v'.

******************************************************************************/

double *matrix_symm_eigen(matrix *m, const char job)
{
	double *eigenval = allocate(m->max_row, sizeof(double), false);

	call_dsyev(job, 'l', m->max_row, m->data, m->max_row, eigenval);
	return eigenval;
}

/******************************************************************************

 Function matrix_row_quadr(): perform a 1/3-Simpson quadrature rule in the p-th
 row of matrix m.

******************************************************************************/

double matrix_row_quadr(const matrix *m, const int p, const bool use_omp)
{
	ASSERT_ROW_INDEX(m, p)

	const int q_max = (m->max_col%2 == 0? m->max_col : m->max_col - 1);

	double sum = DATA_OFFSET(m, p, 0) + DATA_OFFSET(m, p, q_max - 1);

	#pragma omp parallel for default(none) shared(m) reduction(+:sum) if(use_omp)
	for (int q = 1; q < (q_max - 2); q += 2)
	{
		sum += 4.0*DATA_OFFSET(m, p, q) + 2.0*DATA_OFFSET(m, p, q + 1);
	}

	return sum/3.0;
}

/******************************************************************************

 Function matrix_col_quadr(): perform a 1/3-Simpson quadrature rule in the q-th
 column of matrix m.

******************************************************************************/

double matrix_col_quadr(const matrix *m, const int q, const bool use_omp)
{
	ASSERT_COL_INDEX(m, q)

	const int p_max = (m->max_row%2 == 0? m->max_row : m->max_row - 1);

	double sum = DATA_OFFSET(m, 0, q) + DATA_OFFSET(m, p_max - 1, q);

	#pragma omp parallel for default(none) shared(m) reduction(+:sum) if(use_omp)
	for (int p = 1; p < (p_max - 2); p += 2)
	{
		sum += 4.0*DATA_OFFSET(m, p, q) + 2.0*DATA_OFFSET(m, p + 1, q);
	}

	return sum/3.0;
}

/******************************************************************************

 Function matrix_null(): return true if all elements are zero. Return false
 otherwise.

******************************************************************************/

bool matrix_null(const matrix *m)
{
	const int n_max = m->max_row*m->max_col;

	for (int n = 0; n < n_max; ++n)
	{
		if (m->data[n] != 0.0) return false;
	}

	return true;
}

/******************************************************************************

 Function matrix_positive(): return true if all elements are positive. Return
 false otherwise.

******************************************************************************/

bool matrix_positive(const matrix *m)
{
	const int n_max = m->max_row*m->max_col;

	for (int n = 0; n < n_max; ++n)
	{
		if (m->data[n] < 0.0) return false;
	}

	return true;
}

/******************************************************************************

 Function matrix_negative(): return true if all elements are negative. Return
 false otherwise.

******************************************************************************/

bool matrix_negative(const matrix *m)
{
	const int n_max = m->max_row*m->max_col;

	for (int n = 0; n < n_max; ++n)
	{
		if (m->data[n] >= 0.0) return false;
	}

	return true;
}

/******************************************************************************

 Function matrix_has_nan(): return true if at least one element has a NaN.
 Return false otherwise.

******************************************************************************/

bool matrix_has_nan(const matrix *m)
{
	const int n_max = m->max_row*m->max_col;

	for (int n = 0; n < n_max; ++n)
	{
		if (isnan(m->data[n]) != 0) return true;
	}

	return false;
}
/*****************************************************************************

 Function matrix_save(): saves a matrix object in the disk. Binary format is
 used.

******************************************************************************/

void matrix_save(const matrix *m, const char filename[])
{
	FILE *output = fopen(filename, "wb");

	if (output == NULL)
	{
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	int status;

	status = fwrite(&m->max_row, sizeof(int), 1, output);
	ASSERT(status == 1)

	status = fwrite(&m->max_col, sizeof(int), 1, output);
	ASSERT(status == 1)

	status = fwrite(m->data, sizeof(double), m->max_row*m->max_col, output);
	ASSERT(status == m->max_row*m->max_col)

	fclose(output);
}

/******************************************************************************

 Function matrix_load(): load a matrix object from the disk which has been
 written by matrix_save().

******************************************************************************/

matrix *matrix_load(const char filename[])
{
	FILE *input = fopen(filename, "rb");

	if (input == NULL)
	{
		PRINT_ERROR("unable to open %s\n", filename)
		exit(EXIT_FAILURE);
	}

	int status, max_row, max_col;

	status = fread(&max_row, sizeof(int), 1, input);
	ASSERT(status == 1)

	status = fread(&max_col, sizeof(int), 1, input);
	ASSERT(status == 1)

	matrix *m = matrix_alloc(max_row, max_col, false);

	status = fread(m->data, sizeof(double), m->max_row*m->max_col, input);
	ASSERT(status == m->max_row*m->max_col)

	fclose(input);

	return m;
}

/******************************************************************************

 Function matrix_read(): reads max_row lines and max_col columns from a data
 file and return them as a matrix.

 NOTE: column data separated either by spaces or tabs are expected.

******************************************************************************/

matrix *matrix_read(FILE *input, const int max_row, const int max_col)
{
	ASSERT(max_row > 0)
	ASSERT(max_col > 0)
	ASSERT(input != NULL)

	matrix *m = matrix_alloc(max_row, max_col, false);

	rewind(input);
	char line[MAX_LINE_LENGTH] = "\n";

	for (int p = 0; p < max_row; ++p)
	{
		if (fgets(line, sizeof(line), input) == NULL) break;

		if ((line[0] != '#') && (line[0] != '\n') && (line[0] != '\0'))
		{
			char *token = strtok(line, " \t");

			for (int q = 0; q < max_col; ++q)
			{
				if (token != NULL) DATA_OFFSET(m, p, q) = atof(token);
				token = strtok(NULL, " \t");
			}
		}
	}

	return m;
}

/******************************************************************************

 Function matrix_write(): writes to a data file max_row lines and max_col
 columns of a matrix m.

******************************************************************************/

void matrix_write(const matrix *m,
                  FILE *output, const int max_row, const int max_col)
{
	ASSERT(max_row > 0)
	ASSERT(max_col > 0)
	ASSERT(output != NULL)

	const int p_max = min(max_row, m->max_row);
	const int q_max = min(max_col, m->max_col);

	for (int p = 0; p < p_max; ++p)
	{
		for (int q = 0; q < q_max; ++q)
		{
/*
 *			NOTE: "% -8e\t" print numbers left-justified with an invisible
 *			plus sign, if any, 8 digits wide in scientific notation + tab.
 */
			fprintf(output, "% -8e\t", DATA_OFFSET(m, p, q));
		}

		fprintf(output, "\n");
	}
}

/******************************************************************************

 Function matrix_sizeof(): return the total size in bits of the content hold by
 matrix m.

******************************************************************************/

size_t matrix_sizeof(const matrix *m)
{
	return 2*sizeof(int) + m->max_row*m->max_col*sizeof(double);
}

/******************************************************************************

 Function matrix_about(): print in a given output file the conditions in which
 the module was compiled.

******************************************************************************/

void matrix_about(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "# build date   = %s\n", __DATE__);
	fprintf(output, "# source code  = %s\n", __FILE__);
	fprintf(output, "# data layout  = %s\n", "row-major scheme");
	fprintf(output, "# lin. algebra = %s\n", LINEAR_ALGEBRA_LIB);
}

/******************************************************************************

 About
 -----

 This module defines wrappers for popular implementations of the BLAS library,
 including those using GPUs, and selectable during compilation.

******************************************************************************/

#if defined(USE_MKL)
	#include <mkl.h>
	#define LINEAR_ALGEBRA_LIB "Intel Math Kernel Library (MKL)"
#endif

#if defined(USE_MAGMA)
	#include <magma_v2.h>
	#include <magma_lapack.h>
	#define LINEAR_ALGEBRA_LIB "MAGMA"

	/* NOTE: blas_init_gpu() will setup the queue. */
	static magma_queue_t gpu_queue;
#endif

#if !defined(LINEAR_ALGEBRA_LIB)
	#include "gsl_lib.h"
	#define LINEAR_ALGEBRA_LIB "GNU Scientific Library (GSL)"

	/* NOTE: typedef for compatibility with MKL. */
	typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
#endif

#include "blas_lib.h"

/******************************************************************************

 Function blas_init_gpu(): allocate resources needed by MAGMA and shall be
 called before any call to functions of this module when the macro USE_MAGMA
 is defined during the compilation.

******************************************************************************/

void blas_init_gpu()
{
	#if defined(USE_MAGMA)
	{
		if (magma_init() != MAGMA_SUCCESS)
		{
			PRINT_ERROR("%s\n", "magma_init() failed")
			exit(EXIT_FAILURE);
		}

		magma_queue_create(0, &gpu_queue);
	}
	#else
	{
		PRINT_ERROR("%s\n", "macro USE_MAGMA not defined")
	}
	#endif
}

/******************************************************************************

 Function blas_end_gpu(): free resources allocated by blas_init_gpu().

******************************************************************************/

void blas_end_gpu()
{
	#if defined(USE_MAGMA)
	{
		if (magma_finalize() != MAGMA_SUCCESS)
		{
			PRINT_ERROR("%s\n", "magma_finalize() failed")
		}

		magma_queue_destroy(gpu_queue);
	}
	#else
	{
		PRINT_ERROR("%s\n", "macro USE_MAGMA not defined")
	}
	#endif
}

/******************************************************************************

 Wrapper blas_dgemm(): a general call to dgemm interfacing the many different
 possible libraries. For the meaning of each input parameter, please refer to
 dgemm documentation in the netlib repository:

 http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html

 NOTE: Row-major matrices only.

******************************************************************************/

void blas_dgemm(const char trans_a,
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
		const magma_trans_t a_form = (trans_a == 't'? MagmaNoTrans : MagmaTrans);
		const magma_trans_t b_form = (trans_b == 't'? MagmaNoTrans : MagmaTrans);

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
		const CBLAS_TRANSPOSE a_form = (trans_a == 't'? CblasTrans : CblasNoTrans);
		const CBLAS_TRANSPOSE b_form = (trans_b == 't'? CblasTrans : CblasNoTrans);

		cblas_dgemm(CblasRowMajor, a_form, b_form,
		            m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	}
	#endif
}

/******************************************************************************

 Function blas_about(): print in a given output file the conditions in which
 the module was compiled.

******************************************************************************/

void blas_about(FILE *output)
{
	ASSERT(output != NULL)

	fprintf(output, "# build date   = %s\n", __DATE__);
	fprintf(output, "# source code  = %s\n", __FILE__);
	fprintf(output, "# lin. algebra = %s\n", LINEAR_ALGEBRA_LIB);
}

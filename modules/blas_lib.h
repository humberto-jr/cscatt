#if !defined(BLAS_LIB_HEADER)
	#define BLAS_LIB_HEADER
	#include "globals.h"

	void blas_init_gpu();

	void blas_end_gpu();

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
	                const int ldc);

	void blas_about(FILE *output);
#endif

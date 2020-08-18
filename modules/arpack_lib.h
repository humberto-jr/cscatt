#if !defined(ARPACK_LIB_HEADER)
	#define ARPACK_LIB_HEADER
	#include "c_lib.h"

	void dsaupd_(int *ido,
	             char bmat[],
	             int *n,
	             char which[],
	             int *nev,
	             double *tol,
	             double resid[],
	             int *ncv,
	             double v[],
	             int *ldv,
	             int iparam[],
	             int ipntr[],
	             double workd[],
	             double workl[],
	             int *lworkl,
	             int *info);

	void dseupd_(bool *rvec,
	             char howmny[],
	             bool select[],
	             double d[],
	             double z[],
	             int *ldz,
	             double *simga,
	             char bmat[],
	             int *n,
	             char which[],
	             int *nev,
	             double *tol,
	             double resid[],
	             int *ncv,
	             double v[],
	             int *ldv,
	             int iparam[],
	             int ipntr[],
	             double workd[],
	             double workl[],
	             int *lworkl,
	             int *info);
/*
	inline static void dsaupd(int ido,
	                          char bmat[],
	                          int n,
	                          char which[],
	                          int nev,
	                          double tol,
	                          double resid[],
	                          int ncv,
	                          double v[],
	                          int ldv,
	                          int iparam[],
	                          int ipntr[],
	                          double workd[],
	                          double workl[],
	                          int lworkl,
	                          int info)
	{
		dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv,
		        v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
	}*/
#endif

#if !defined(PES_HEADER)
	#define PES_HEADER
	#include "coor.h"

	double pes(const jacobi_coor *x);

	double pec(const char arrang, const double r);

	void pes_about(FILE *output);
#endif

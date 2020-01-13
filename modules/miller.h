#if !defined(MILLER_HEADER)
	#define MILLER_HEADER
	#include "globals.h"
	#include "matrix.h"

	double miller_jcp69_rot_integral(const char arrang, const int lambda,
	                                 const double r, const double R);

	double miller_jcp69_vib_integral(const int lambda,
	                                 const char arrang,
	                                 const matrix *wavef_a,
	                                 const matrix *wavef_b,
	                                 const double r_min,
	                                 const double r_max,
	                                 const double R);
#endif

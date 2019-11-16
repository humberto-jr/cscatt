#if !defined(GLOBALS_HEADER)
	#define GLOBALS_HEADER
	#include <stdbool.h>
	#include <stdio.h>
	#include <float.h>
	#include <math.h>

	/******************************************************************************

	 Type real: defines a general real type number to which precision, double or
	 float, can be tuned during compilation.

	 NOTE: not used.

	******************************************************************************/

	#if defined(USE_SINGLE_PRECISION)
		typedef float real;
	#else
		typedef double real;
	#endif

	/******************************************************************************

	 Macro MAX_LINE_LENGTH(): defines a default size for strings to be tuned during
	 compilation.

	******************************************************************************/

	#if !defined(MAX_LINE_LENGTH)
		#define MAX_LINE_LENGTH 1024
	#endif

	/******************************************************************************

	 Macro INF: use an upper bound for real numbers with double precision to define
	 infinity, and likewise, -infinity.

	 NOTE: DBL_MAX from float.h is chosen, however any other reasonably large
	 number is possible.

	******************************************************************************/

	#define INF DBL_MAX

	/******************************************************************************

	 Macro M_PI(): is the value of pi as defined in the GSL library. Its usage is
	 so common that it is worth redefining it here for compatibility reasons when
	 GSL is not available.

	******************************************************************************/

	#if !defined(M_PI)
		#define M_PI 3.14159265358979323846264338328
	#endif

	/******************************************************************************

	 Macro PRINT_ERROR(): writes an error message in the C stderr using the format
	 "# file.c, function(), line n: formatted message here".

	******************************************************************************/

	#define PRINT_ERROR(format, ...)                                            \
	{                                                                           \
		fprintf(stderr, "# %s, %s(), line %d: ", __FILE__, __func__, __LINE__); \
		fprintf(stderr, format, ##__VA_ARGS__);                                 \
	}

	/******************************************************************************

	 Macro ASSERT(): checks the validity (true/false) of a given expression. If
	 false, it shall exit the runtime with a formatted message in the C stderr.

	******************************************************************************/

	#define ASSERT(expr)                                  \
	{                                                     \
		if (!(expr))                                      \
		{                                                 \
			PRINT_ERROR("assertion '%s' failed\n", #expr) \
			exit(EXIT_FAILURE);                           \
		}                                                 \
	}

	/******************************************************************************

	 Macro AS_STRING(): return a given macro as a string between double quotes.

	******************************************************************************/

	#define AS_STRING(macro) #macro

	/******************************************************************************

	 Macro PRINT_MACRO(): return the value of a given macro as a string.

	******************************************************************************/

	#define PRINT_MACRO(macro) AS_STRING(macro)

	/******************************************************************************

	 Function allocate(): allocate resources for n elements with double precision
	 and return a pointer pointing the first one (initialized to zero if set_zero
	 = true).

	******************************************************************************/

	static inline double *allocate(const int n, const bool set_zero)
	{
		ASSERT(n > 0)

		double *pointer
			= (set_zero? calloc(n, sizeof(double)) : malloc(sizeof(double)*n));

		if (pointer == NULL)
		{
			PRINT_ERROR("unable to allocate resources for %d elements\n", n)
			exit(EXIT_FAILURE);
		}

		return pointer;
	}

	/******************************************************************************

	 Function as_double(): is a simple interface to a very common operation of
	 casting an integer, n, as a double precision number.

	******************************************************************************/

	inline static double as_double(const int n)
	{
		return (double) n;
	}
	/******************************************************************************

	 Function side_c(): uses the law of cosines to resolve the side c opposite to
	 the interior angle C (in degree) of a triangle with sides a, b and c.

	******************************************************************************/

	inline static double side_c(const double side_a,
	                            const double side_b, const double angle_c)
	{
		return sqrt(side_a*side_a + side_b*side_b - 2.0*side_a*side_b*cos(angle_c*M_PI/180.0));
	}

	/******************************************************************************

	 Function angle_c(): uses the law of cosines to resolve the angle C (in degree)
	 opposite to the side c of a triangle with sides a, b and c.

	******************************************************************************/

	inline static double angle_c(const double side_a,
	                             const double side_b, const double side_c)
	{
		return acos((side_c*side_c - side_a*side_a - side_b*side_b)/(-2.0*side_a*side_b))*180.0/M_PI;
	}

	/******************************************************************************

	 Function centr_term(): returns the actual centrifugal potential term at x for
	 a given angular momentum l, provided the mass under consideration. All in
	 atomic units.

	******************************************************************************/

	inline static double centr_term(const int l,
	                                const double mass, const double x)
	{
		return as_double(l*(l + 1))/(2.0*mass*x*x);
	}

	/******************************************************************************

	 Function sinc(): returns sin(x)/x for a given x.

	******************************************************************************/

	inline static double sinc(const double x)
	{
		return sin(x)/x;
	}

	/******************************************************************************

	 Function parity(): return the parity p = (-1)^l for a given l; either p = +1
	 or p = -1.

	******************************************************************************/

	inline static int parity(const int l)
	{
		return (int) pow(-1.0, l);
	}
#endif

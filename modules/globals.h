#if !defined(GLOBALS_HEADER)
	#define GLOBALS_HEADER
	#include "clib.h"
	#include <omp.h>

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

	static inline void *allocate(const int n,
	                             const int data_size, const bool set_zero)
	{
		void *pointer = (set_zero? calloc(n, data_size) : malloc(n*data_size));

		if (pointer == NULL)
		{
			PRINT_ERROR("unable to allocate resources for %d elements of %d bits\n", n, data_size)
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

	 Function as_double(): is a simple interface to a very common operation of
	 casting an double, n, as an integer number.

	******************************************************************************/

	inline static int as_int(const double n)
	{
		return (int) n;
	}
	/******************************************************************************

	 Function tools_wall_time(): return the wall time in units of seconds.

	******************************************************************************/

	inline static double wall_time()
	{
		return omp_get_wtime();
	}

	/******************************************************************************

	 Function sinc(): returns sin(x)/x for a given x.

	******************************************************************************/

	inline static double sinc(const double x)
	{
		return sin(x)/x;
	}

	/******************************************************************************

	 Function sigmoid(): return the sigmoid function at x.

	******************************************************************************/

	inline static double sigmoid(const double x)
	{
		return 1.0/(1.0 + exp(-x));
	}

	/******************************************************************************

	 Function left_trim(): trim a string s in the left side.

	******************************************************************************/

	inline static char *left_trim(char s[])
	{
		while (isspace(*s)) ++s;
		return s;
	}

	/******************************************************************************

	 Function right_trim(): trim a string s in the right side.

	******************************************************************************/

	inline static char *right_trim(char s[])
	{
		char *s_backward = s + strlen(s);

		while (isspace(*(--s_backward)));
		*(s_backward + 1) = '\0';

		return s;
	}

	/******************************************************************************

	 Function trim(): trim a string s in both sides.

	******************************************************************************/

	inline static char *trim(char s[])
	{
		return right_trim(left_trim(s));
	}
#endif

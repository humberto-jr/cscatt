#if !defined(GLOBALS_HEADER)
	#define GLOBALS_HEADER
	#include <omp.h>
	#include "c_lib.h"

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

	#define PRINT_ERROR(format, ...)                                          \
	{                                                                         \
	  fprintf(stderr, "# %s, %s(), line %d: ", __FILE__, __func__, __LINE__); \
	  fprintf(stderr, format, ##__VA_ARGS__);                                 \
	}

	/******************************************************************************

	 Macro ASSERT(): checks the validity (true/false) of a given expression. If
	 false, it shall exit the runtime with a formatted message in the C stderr.

	******************************************************************************/

	#define ASSERT(expr)                              \
	{                                                 \
	  if (!(expr))                                    \
	  {                                               \
	    PRINT_ERROR("assertion '%s' failed\n", #expr) \
	    exit(EXIT_FAILURE);                           \
	  }                                               \
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

	 Function allocate(): allocate resources for n elements with data_size bits and
	 return a pointer pointing the first one (initialized to zero if set_zero =
	 true).

	******************************************************************************/

	static inline void *allocate(const int n,
	                             const int data_size, const bool set_zero)
	{
		void *pointer = (set_zero? calloc(n, data_size) : malloc(n*data_size));

		if (pointer == NULL)
		{
			PRINT_ERROR("unable to allocate resources for %d elements of %d bytes\n", n, data_size)
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

	 Function as_int(): is a simple interface to a very common operation of casting
	 an double, n, as an integer number.

	******************************************************************************/

	inline static int as_int(const double n)
	{
		return (int) n;
	}

	/******************************************************************************

	 Function wall_time(): returns the wall time in units of seconds.

	******************************************************************************/

	inline static double wall_time()
	{
		return omp_get_wtime();
	}

	/******************************************************************************

	 Function thread_id(): returns the current OpenMP thread identification.

	******************************************************************************/

	inline static int thread_id()
	{
		return omp_get_thread_num();
	}

	/******************************************************************************

	 Function max_threads(): returns the total OpenMP threads available.

	******************************************************************************/

	inline static int max_threads()
	{
		return omp_get_max_threads();
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

	/******************************************************************************

	 Function time_stamp(): return the current YMDHMS date as a time stamp.

	 Example: 31 May, 2001 09:45:54 AM.

	******************************************************************************/

	inline static const char *time_stamp()
	{
		time_t now = time(NULL);
		const struct tm *info = localtime(&now);

		static char stamp[50];
		strftime(stamp, 50, "%B %d, %Y %I:%M:%S %p", info);

		return stamp;
	}

	/******************************************************************************

	 Function min(): return the min between two integers a and b.

	******************************************************************************/

	inline static int min(const int a, const int b)
	{
		return (a < b? a : b);
	}

	/******************************************************************************

	 Function max(): return the max between two integers a and b.

	******************************************************************************/

	inline static int max(const int a, const int b)
	{
		return (a > b? a : b);
	}

	/******************************************************************************

	 Function find(): scan a given input file searching for the 1st occurrence of a
	 given pattern. It returns the whole line, if found, or an empty string
	 otherwise.

	 NOTE: Lines starting by '#' are ignored.

	******************************************************************************/

	inline static char *find(FILE *input, const char pattern[])
	{
		ASSERT(input != NULL)
		ASSERT(pattern != NULL)

		char *line = allocate(MAX_LINE_LENGTH + 1, sizeof(char), false);

		rewind(input);

		while (fgets(line, MAX_LINE_LENGTH, input) != NULL)
			if (line[0] != '#' && strstr(line, pattern) != NULL) return line;

		line[0] = '\0';

		return line;
	}

	/******************************************************************************

	 Function read_dbl_keyword(): scan a given input file and search for the 1st
	 occurrence of a keyword in the format "[key] = [value]". Where, [value] is a
	 real number such that min < [value] < max. Return [value] if found or a
	 default value otherwise.

	 NOTE: Lines starting by '#' are ignored.

	******************************************************************************/

	inline static double read_dbl_keyword(FILE *input,
	                                      const char key[],
	                                      const double min,
	                                      const double max,
	                                      const double default_value)
	{
		ASSERT(max >= min)

		char *line = find(input, key), *token = NULL;

		if (line[0] != '\0') token = strtok(line, "=");

		if (token != NULL)
		{
			token = trim(token);

			if (strcmp(token, key) == 0)
			{
				token = strtok(NULL, "=");

				const double value = (token != NULL? atof(token) : default_value);
				free(line);

				if (value < min)
					return min;
				else if (value > max)
					return max;
				else
					return value;
			}
		}

		free(line);
		return default_value;
	}

	/******************************************************************************

	 Function read_int_keyword(): scan a given input file and search for the 1st
	 occurrence of a keyword in the format "[key] = [value]". Where, [value] is an
	 integer number such that min < [value] < max. Return [value] if found or a
	 default value otherwise.

	 NOTE: Lines starting by '#' are ignored.

	******************************************************************************/

	inline static int read_int_keyword(FILE *input,
	                                   const char key[],
	                                   const int min,
	                                   const int max,
	                                   const int default_value)
	{
		ASSERT(max >= min)

		char *line = find(input, key), *token = NULL;

		if (line[0] != '\0') token = strtok(line, "=");

		if (token != NULL)
		{
			token = trim(token);

			if (strcmp(token, key) == 0)
			{
				token = strtok(NULL, "=");

				const int value = (token != NULL? atoi(token) : default_value);
				free(line);

				if (value < min)
					return min;
				else if (value > max)
					return max;
				else
					return value;
			}
		}

		free(line);
		return default_value;
	}

	/******************************************************************************

	 Function read_str_keyword(): scan a given input file and search for the 1st
	 occurrence of a keyword in the format "[key] = [value]". Where, [value] is a
	 string with up to MAX_LINE_LENGTH characters. Return [value] if found or the
	 content pointed by replacement.

	 NOTE: Lines starting by '#' are ignored.

	******************************************************************************/

	inline static char *read_str_keyword(FILE *input,
	                                     const char key[], char replacement[])
	{
		ASSERT(replacement != NULL)

		char *line = find(input, key), *token = NULL;

		if (line[0] != '\0') token = strtok(line, "=");

		if (token != NULL)
		{
			token = trim(token);

			if (strcmp(token, key) == 0)
			{
				token = strtok(NULL, "=");

				if (token != NULL)
					token = trim(token);
				else
					token = replacement;
			}
		}
		else
		{
			token = replacement;
		}

		char *value = allocate(strlen(token) + 1, sizeof(char), false);

		strcpy(value, token);
		free(line);

		return value;
	}

	/******************************************************************************

	 Function globals_about(): prints in a given output file the conditions in
	 which the module was compiled.

	******************************************************************************/

	inline static void globals_about(FILE *output)
	{
		ASSERT(output != NULL)

		fprintf(output, "# build date      = %s\n", __DATE__);
		fprintf(output, "# source code     = %s\n", __FILE__);

		#if defined(USE_SINGLE_PRECISION)
			fprintf(output, "# real precision  = single\n");
		#else
			fprintf(output, "# real precision  = double\n");
		#endif

		fprintf(output, "# MAX_LINE_LENGTH = %d\n", MAX_LINE_LENGTH);
		fprintf(output, "# INF             = %8e\n", INF);
	}
#endif

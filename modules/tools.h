#if !defined(TOOLS_HEADER)
	#define TOOLS_HEADER

	char *tools_find_string(FILE *input, const char pattern[]);

	double tools_get_key(FILE *input, const char key[], const double min,
	                     const double max, const double default_value);

	int tools_row_count(FILE *input);

	int tools_col_count(FILE *input);

	const char *tools_time_stamp();

	double tools_wall_time();
#endif

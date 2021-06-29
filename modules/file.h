#if !defined(FILE_HEADER)
	#define FILE_HEADER
	#include "globals.h"

	FILE *file_open(const char filename[], const char mode[]);

	void file_close(FILE **stream);

	bool file_exist(const char filename[]);

	bool file_end(FILE *stream);

	void file_remove(const char filename[]);

	void file_rename(const char old_name[], const char new_name[]);

	void file_init_stdin(const char filename[]);

	void file_init_stdout(const char filename[]);

	char *file_find(FILE *stream, const char pattern[]);

	double file_keyword(FILE *input, const char key[], const double min,
	                    const double max, const double default_value);

	size_t file_row_count(FILE *stream);

	size_t file_col_count(FILE *stream);

	void file_write(const void *buffer,
	                const size_t size, const size_t length, FILE *stream);

	void file_read(void *buffer, const size_t size,
	               const size_t length, FILE *stream, const size_t offset);

	void file_about(FILE *output);
#endif

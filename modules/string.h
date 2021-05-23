#if !defined(STRING_HEADER)
	#define STRING_HEADER

	typedef struct string string;

	void string_increase_storage(string *s);

	string *string_alloc();

	void string_free(string *s);

	void string_set(string *s, const char text[]);

	size_t string_length(const string *s);

	size_t string_max_length(const string *s);

	void string_swap(string *a, string *b);

	string *string_substring(const string *s,
	                         const size_t start, const size_t end);

	void string_crop(string *s, const size_t start, const size_t end);

	string *string_copy(const string *s);

	void string_append(string *s, const char text[]);

	void string_print(const string *s, FILE *stream, const bool end_line);

	int string_count(const string *s, const char pattern[]);

	void string_insert(string *s, const size_t n, const char text[]);

	void string_remove(string *s, const size_t start, const size_t end);

	void string_replace_all(string *s,
	                        const char pattern[], const char insert[]);

	void string_replace_at(string *s, const char pattern[],
	                       const char insert[], const size_t n);

	void string_right_trim(string *s);

	void string_left_trim(string *s);

	void string_trim(string *s);

	void string_reset(string *s);

	char string_at(string *s, const size_t n);

	void string_read_file(string *s, const char filename[]);

	string *string_file_line(const string *s, const size_t n);

	void string_set_lower(string *s);

	void string_set_upper(string *s);

	bool string_compare(const string *a, const string *b);

	int string_find_first(const string *s, const char pattern[]);

	int string_find_from(const string *s,
	                     const char pattern[], const size_t start);

	void string_concatenate(string *a, const string *b);

	char *string_as_array(const string *s);

	void string_tokenize(string *s, const char delim[]);

	int string_token_count(const string *s);

	void string_token_print(const string *s,
	                        const size_t n, FILE *stream, const bool end_line);

	void string_token_print_all(const string *s,
	                            FILE *stream, const bool end_line);

	void string_token_trim(string *s, const size_t n);

	void string_token_trim_all(string *s);

	void string_token_replace(string *s, const size_t n,
	                          const char pattern[], const char insert[]);

	void string_token_replace_all(string *s,
	                              const char pattern[], const char insert[]);

	size_t string_token_length(string *s, const size_t n);
#endif

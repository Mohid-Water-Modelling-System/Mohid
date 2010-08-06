#ifndef stringH
#define stringH

#include <sstream>
#include <string>
#include <cerrno>

#include "definitions.h"

using namespace std;

#ifndef DEFAULT_STRING_LENGTH
  #define DEFAULT_STRING_LENGTH 1024
#endif

typedef enum { _all_, _left_, _right_ } TRIM_TYPE;

class EStringIndexOutOfRange
{
public:
	EStringIndexOutOfRange() {}
};

class EStringConversionError
{
public:
	EStringConversionError() {}
};

class String
{
public:
	String(long start_capacity = 256);
	String(const String &str, long start_capacity = 256);
	String(const char *str, long start_capacity = 256);
	String(const char c, long start_capacity = 256);
	String(int number, const String str_to_repeat);
	~String();

public:
	char &operator [] (long index);
	const char &operator [] (long index) const;
	const String &operator = (const String &str);
	const String &operator = (const String *str);
	const String &operator = (const char *str);	
	const String &operator = (const char c);
	const String &operator = (const int number);
	const String &operator = (const double number);
	bool operator == (const String &str) const;
	bool operator == (const char *str) const;
	bool operator == (const char c) const;
	bool operator != (const String &str) const;
	bool operator != (const char *str) const;
	bool operator != (const char c) const;
	const String &operator += (const String &str);
	const String &operator += (const char *str);
	const String &operator += (const char c);
	const String &operator += (const int l);
	String operator () (long start, long characters = -1);
	String operator () (long start, const String &to_find, bool no_case = false);
	String operator () (long start, const char *to_find, bool no_case = false);

public:
	void SetString(String *str);
	void SetString(const char *str);
	bool Replace(const String &old_token, const String &new_token); 
	bool Replace(const char *old_token, const char *new_token);
	int Compare(const String &str, bool no_case = false, bool partial = false);
	int Compare(const char *str, bool no_case = false, bool partial = false);
	int Compare(String *str, bool no_case = false, bool partial = false);
	void Copy(char *dest, long characters = -1, long start = 0);
	long Find(const String &token, bool no_case = false);
	long Find(const char *token, bool no_case = false);
	const char *CharPtr(void) const;
	long Length(void);
	void ClearRightZeros(void);
	void Cut(const int start);
	void RemoveSpaces(void);
	void ToUpper(void);
	void ToLower(void);
	void Trim(TRIM_TYPE trim_type = _left_); 
	float AsFloat(void);
	double AsDouble(void);
	int AsInt(void);
	long AsLong(void);
	bool IsEmpty(void);
	bool IsAmong(const int index, const String &str);
	bool IsAmong(const int index, const char *str);
	int FindUpper(long start);

private:
	char *str_ptr;
	long capacity;
	long length;

private:
	void _CheckIndex(long index) const;
	void _SetString(const char *str, const long new_capacity = 0);
	void _TrimLeft(void);
	void _TrimRight(void);
	bool _CheckEmpty(void) const;
	void _ChangeCapacity(long new_capacity, const bool keep = true);
	void _Init(void);
};

#endif
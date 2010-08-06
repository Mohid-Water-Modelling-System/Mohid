#include "string.h"

//=============================================================================
//  Constructors & Destructors
//=============================================================================
String::String(long start_capacity)
{
	_Init();
	_ChangeCapacity(start_capacity, false);
}

String::String(const String &str, long start_capacity)
{
	_Init();
	_SetString(str.str_ptr, start_capacity);
}

String::String(const char *str, long start_capacity)
{
	_Init();
	_SetString(str, start_capacity);
}

String::String(const char c, long start_capacity)
{
	char  str[2];

	str[0] = c;
	str[1] = '\0';

	_Init();
	_SetString(str, start_capacity);
}

String::String(int number, const String str_to_repeat)
{
	_Init();

	String str_repetition("");
	for (int i = 0; i < number; i++)
		str_repetition += str_to_repeat;

	_SetString(str_repetition.CharPtr()); 
}

String::~String()
{
	if (str_ptr != NULL)
		delete [] str_ptr;
}

//=============================================================================
//  Operators
//=============================================================================

char &String::operator [] (long index)
{
	_CheckIndex(index);
	return str_ptr[index];
}

const char &String::operator [] (long index) const
{
	_CheckIndex(index);
	return str_ptr[index];
}

const String &String::operator = (const String &str)
{
	_SetString(str.str_ptr);		
	return *this;
}

const String &String::operator = (const String *str)
{
	_SetString(str->str_ptr);		
	return *this;
}

const String &String::operator = (const char *str)
{
	_SetString(str);		
	return *this;
}

const String &String::operator = (const char c)
{
	char str[2];

	str[0] = c;
	str[1] = '\0';

	_SetString(str);		
	return *this;
}

const String &String::operator = (const int number)
{
	stringstream ss;
  ss << number;

	_SetString(ss.str().c_str());
	return *this;
}

const String &String::operator = (const double number)
{
	stringstream ss;
  ss << number;

	_SetString(ss.str().c_str());
	return *this;
}

bool String::operator == (const String &str) const
{
	for (long i = 0; i < length; i++)
		if (str_ptr[i] != str.str_ptr[i])
			return false;

	return true;
}

bool String::operator == (const char *str) const
{
	bool empty_l = _CheckEmpty();
	bool empty_r = (str == NULL || str[0] == '\0');

	if (empty_l && empty_r)
		return true;
	else if (empty_l || empty_r)
		return false;

	long len = strlen(str);
	if (length != len)
		return false;

	for (long i = 0; i < length; i++)
		if (str_ptr[i] != str[i])
			return false;

	return true;
}

bool String::operator == (const char c) const
{
	if (length <= 0 || length > 1)
		return false;

	return (c == str_ptr[0]);
}

bool String::operator != (const String &str) const
{
	return !(*this == str);
}

bool String::operator != (const char *str) const
{
	return !(*this == str);
}

bool String::operator != (const char c) const
{
	return !(*this == c);
}

const String &String::operator += (const String &str)
{
	return (*this += str.CharPtr());
}

const String &String::operator += (const char *str)
{
	if (str == NULL || str[0] == '\0')
		return *this;

	long len = strlen(str);
	if ((len + length) > capacity)
		_ChangeCapacity(length + len);

	memcpy(&str_ptr[length], str, len + 1);
	length += len;

	return *this;
}

const String &String::operator += (const char c)
{
	if ((length + 1) > capacity)
		_ChangeCapacity(capacity + 50);

	str_ptr[length++] = c;
	str_ptr[length] = '\0';

	return *this;
}

const String &String::operator += (const int number)
{
  stringstream ss;
  ss << number;

	long len = ss.str().length();
	if ((len + length) > capacity)
		_ChangeCapacity(length + len);

	memcpy(&str_ptr[length], ss.str().c_str(), len + 1);
	length += len;

	return *this;
}

String String::operator () (long start, long characters)
{
	_CheckIndex(start);

	long len;
	if (characters < 0 || (start + characters) > length)
		len = length - start;
	else
		len = characters;

	String t = "";
	for (long i = 0; i < len; i++, start++)
		t += str_ptr[start];

	return t;
}

String String::operator () (long start, const String &to_find, bool no_case)
{
	return operator () (start, to_find.CharPtr(), no_case);
}

String String::operator () (long start, const char *to_find, bool no_case)
{
	_CheckIndex(start);

	long p = Find(to_find, no_case);
	if (p == -1)
		return operator ()(start);
	else
		return operator ()(start, p);
}

//=============================================================================
//  Public functions
//=============================================================================
void String::SetString(String *str)
{
	_SetString(str->CharPtr());		
}

void String::SetString(const char *str)
{
	_SetString(str);
}

bool String::Replace(const String &old_token, const String &new_token)
{
	return Replace(old_token.CharPtr(), new_token.CharPtr());
}

bool String::Replace(const char *old_token, const char *new_token)
{
	long l, l1, l2;
	char *ptr_start;

	if (_CheckEmpty())
		return false;

	ptr_start = strstr (str_ptr, old_token);

	if (ptr_start == NULL)
		return false;

	l  = (long) strlen (str_ptr);
	l1 = (long) strlen (old_token);
	l2 = (long) strlen (new_token);

	memmove (ptr_start + l2, ptr_start + l1, l - (ptr_start - str_ptr + l1) + 1);
	memcpy (ptr_start, new_token, l2);

	return true;
}

int String::Compare(const String &str, bool no_case, bool partial)
{
	return Compare(str.CharPtr(), no_case, partial);
}

int String::Compare(const char *str, bool no_case, bool partial)
{
	//
	//Compares str in relation to the current string on this instance
	// if str is lesser return -1, if is grater, return 1, 0 otherwise
	//

	bool empty_l = _CheckEmpty();
	bool empty_r = (str == NULL || str[0] == '\0');

	if (empty_l && empty_r)
		return 0;
	else if (empty_l)
		return 1;
	else if (empty_r)
		return -1;

	int c1, c2;
	char *str1 = str_ptr;
	const char *str2 = str;

	if (no_case)
	{
		do
		{
			c1 = tolower (*str1++);
			c2 = tolower (*str2++);
		}
		while (c1 != '\0' && c2 != '\0' && c1 == c2);
	}
	else
	{
		do
		{
			c1 = *str1++;
			c2 = *str2++;
		}
		while (c1 != '\0' && c2 != '\0' && c1 == c2);
	}

	if (c1 == '\0' || c2 == '\0')
	{
		if (c1 == '\0' && c2 == '\0')
			return 0;

		if (c2 == '\0')
		{
			if (partial)
				return 0;
			else
				return -1;
		}
		else
			return 1;
	}
	else
	{
		if (c1 > c2)
			return -1;
		else
			return 1;
	}
}

int String::Compare(String *str, bool no_case, bool partial)
{
	return Compare(str->CharPtr(), no_case, partial);
}

void String::Copy(char *dest, long characters, long start)
{
	_CheckIndex(start);

	if (characters < 0)
		characters = length;

	if ((start + characters) > length)
		characters = length - start;

	long  i;
	for (i = 0; i < characters; i++)
		dest[i] = str_ptr[i + start];
	dest[i] = '\0';
}

long String::Find(const String &token, bool no_case)
{
	return Find(token.CharPtr(), no_case);
}

long String::Find(const char *token, bool no_case)
{
	if (_CheckEmpty())
		return -1;

	const char *t, *tok, *str = str_ptr;
	long i;

	if (no_case)
	{
		i = 0;
		while (*str != '\0')
		{
			tok = token;
			if (tolower(*str++) == tolower(*tok++))
			{
				t = str;
				while (*t != '\0' && *tok != '\0')
					if (tolower(*t++) != tolower(*tok++))
						break;

					if(*tok == '\0')
						return i;
			}

			i++;
		}
	}
	else
	{
		i = 0;
		while (*str != '\0')
		{
			tok = token;
			if (*str++ == *tok++)
			{
				t = str;
				while (*t != '\0' && *tok != '\0')
					if (*t++ != *tok++)
						break;

					if(*tok == '\0')
						return i;
			}

			i++;
		}
	}
	return -1;
}

const char *String::CharPtr(void) const
{
	return (const char *)str_ptr;
}

void String::ClearRightZeros(void)
{
	if (_CheckEmpty())
		return;

	long old_length = length;
	while (length > 1 && str_ptr[length - 1] == '0') 
		--length;

	if (old_length > length && str_ptr[length] == '.')
		length--;

	str_ptr[length] = '\0';
}

void String::Cut(const int start)
{
	if (_CheckEmpty())
		return;

	_CheckIndex(start);

	str_ptr[start] = '\0';
	length = start;
}

void String::RemoveSpaces(void)
{
	if (_CheckEmpty())
		return;

	long i, j;

  for (i = j = 0; str_ptr[i] != '\0'; i++)
  {
    if (str_ptr[i] != ' ')
      str_ptr[j++] = str_ptr[i];
  }
  str_ptr[j] = '\0';
}

void String::ToUpper(void)
{
	for (long i = 0; i < length; i++)
		str_ptr[i] = (char) toupper((int) str_ptr[i]);
}

void String::ToLower(void)
{
	for (long i = 0; i < length; i++)
		str_ptr[i] = (char) tolower((int) str_ptr[i]);
}

void String::Trim(TRIM_TYPE trim_type)
{
	switch(trim_type)
	{
	case _left_:
		_TrimLeft();
		break;
	case _right_:
		_TrimRight();
		break;
	case _all_:
		_TrimLeft();
		_TrimRight();
		break;
	}
}

float String::AsFloat(void)
{
	return 0.0;
}

double String::AsDouble(void)
{
	char *temp = str_ptr;
	double r;

	r = strtod(temp, NULL);
	if (errno == ERANGE)
		throw EStringConversionError();

	return r;
}

int String::AsInt(void)
{
	char *temp = str_ptr;
	int r;

	r = atoi(temp);
	if (errno == ERANGE)
		throw EStringConversionError();

	return r;
}

long String::AsLong(void)
{
	char *temp = str_ptr;
	long r;

	r = atol(temp);
	if (errno == ERANGE)
		throw EStringConversionError();

	return r;
}

bool String::IsEmpty(void)
{
	return _CheckEmpty();
}

bool String::IsAmong(const int index, const String &str)
{
	if (_CheckEmpty())
		return false;

	for (long i = 0; i < str.length; i++)
		if (str_ptr[index] == str[i])
			return true;

	return false;
}

bool String::IsAmong(const int index, const char *str)
{
	if (_CheckEmpty())
		return false;

	for (long i = 0; str[i] != '\0'; i++)
		if (str_ptr[index] == str[i])
			return true;

	return false;
}

int String::FindUpper(long start)
{
	if (start >= length)
		return -1;

	char *ptr = &str_ptr[start];
	while (*ptr != '\0')
	{
		if (isupper(*ptr))
			return start;

		start++;
		ptr++;
	}

	return -1;
}

long String::Length(void)
{
	return length;
}

//=============================================================================
//  Private functions
//=============================================================================

void String::_CheckIndex(long index) const
{
	if (index > length)
		throw EStringIndexOutOfRange();
}

void String::_SetString(const char *str, long new_capacity)
{
	long len;

	if (str != NULL)
		len = strlen(str);
	else
		len = 0;

	new_capacity = max(len, new_capacity);

	if ((new_capacity - capacity) > 0)
		_ChangeCapacity(max(len, new_capacity), false);

	length = len;

	if (length > 0)
	{
		memcpy(str_ptr, str, length + 1);
		str_ptr[length] = '\0';
	}
	else if(capacity > 0)
	{
		str_ptr[0] = '\0';
	}

}

void String::_TrimLeft(void)
{
	if (_CheckEmpty())
		return;

	long i;
	for (i = 0; str_ptr[i] == ' ' && i < length; i++);
	length -= i;
	memcpy(str_ptr, &(str_ptr[i]), length + 1);
}

void String::_TrimRight(void)
{
	if (_CheckEmpty())
		return;

	for (; str_ptr[length - 1] == ' ' && length > 0; length--);
	str_ptr[length] = '\0';
}

bool String::_CheckEmpty(void) const
{
	if(str_ptr == NULL || length == 0)
		return true;

	return false;
}

void String::_ChangeCapacity(const long new_capacity, const bool keep)
{
	char *new_str = NULL;

	if (new_capacity <= 0)
	{
		delete [] str_ptr;

		str_ptr  = NULL;
		length   = 0;
		capacity = 0;

		return;
	}

	if (new_capacity != capacity)
	{
		length = min(length, new_capacity);
		new_str = new char [new_capacity + 1];

		if (keep && str_ptr != NULL)
		{
			memcpy(new_str, str_ptr, length);
			delete [] str_ptr;
		}

		str_ptr = new_str;
	}

	str_ptr[length] = '\0';
	capacity = new_capacity;
}

void String::_Init(void)
{
	this->str_ptr  = NULL;
	this->length   = 0;
	this->capacity = 0;
}
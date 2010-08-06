#ifndef listofpointersH
#define listofpointersH

#include "string.h"
#include "exceptionhandler.h"

#ifndef LIST_CAPACITY_ADD_STEP
	#define LIST_CAPACITY_ADD_STEP 10
#endif

#ifndef LIST_INITIAL_CAPACITY
	#define LIST_INITIAL_CAPACITY 1
#endif

template <class T>
class ListOfPointers
{
public:
	ListOfPointers(int initial_capacity = LIST_INITIAL_CAPACITY);
	~ListOfPointers();

public:
	inline T *operator[] (const int index) { return Element(index); }

public:
	inline int Count() { return _count; }
	inline int Capacity() { return _count; }
	inline bool OutOfBounds() { return _out_of_bounds; }
	inline int LastIndex() { return _count - 1; }
	inline bool IsEmpty() { if (count <= 0) return true; return false; }
	inline T *Element(int index) 
	{
		if (!_CheckIndex(index))
		{
			_out_of_bounds = true;
			return NULL;
		}

		_out_of_bounds = false;
		return _items[index];
	}
	inline T *Items(int index) { return Element(index); }
	inline T *Last() { return Element(0); }
	inline T *First() { return Element(_count - 1); }

public:
	int AddNew(T *to_copy);
	bool Store(const int index, T * item);

	bool Swap(const int i1, const int i2);

	void Clear();
	void SetNewCapacity(int new_capacity);
	
	void CopyTo(ListOfPointers<T> *copy);

public:
	bool Sort(int start = 0, bool no_case = true);

	T *Search(String *name, bool no_case = false, bool partial = false);
	T *Search(char *name, bool no_case = false, bool partial = false);
	T *Search(String *name, int &index, bool no_case = false, bool partial = false);
	T *Search(char *name, int &index, bool no_case = false, bool partial = false);

protected:
	bool _CheckIndex(int index, bool throw_exception = false);

protected:
	T **_items; //matrix of T types elements
	T *t_item; //temporary item used on many operations

	bool _out_of_bounds; //used to indicate if the last operation was done with a valid index

	int _count, //counts the number of "slots" used
	    _capacity; //counts the number of "slots" on this list

	char message[2049];
	String t_string;
};


//---------------------------------------------------------------------------
// Constructors & Destructors
//---------------------------------------------------------------------------
template <class T>
ListOfPointers<T>::ListOfPointers(int initial_capacity)
{
	_count = 0;
	_capacity = initial_capacity;
	_out_of_bounds = false;

	_items = NULL; 
	_items = new T *[_capacity];
}

template <class T>
ListOfPointers<T>::~ListOfPointers()
{
	delete [] _items;
}

//---------------------------------------------------------------------------
// Add and Store functions
//---------------------------------------------------------------------------
template <class T>
inline int ListOfPointers<T>::AddNew(T *to_copy)
{
	if (_count == _capacity)
		SetNewCapacity(_capacity + LIST_CAPACITY_ADD_STEP);

	_items[_count] = to_copy;

	_count++;
	return _count - 1;
}

template <class T>
bool ListOfPointers<T>::Store(const int index, T * item)
{
	if (!_CheckIndex(index)) 
		return false;

	_items[index] = item;
	return true;
}

//---------------------------------------------------------------------------
// Element operations
//---------------------------------------------------------------------------
template <class T>
inline bool ListOfPointers<T>::Swap(const int i1, const int i2)
{
	if (!_CheckIndex(i1) || !_CheckIndex(i2))
	{
		_out_of_bounds = true;
		return false;
	}

	t_item = _items[i1];
	_items[i1] = _items[i2];
	_items[i2] = t_item;

	return true;
}

//---------------------------------------------------------------------------
// Capacity operations
//---------------------------------------------------------------------------
template <class T>
inline void ListOfPointers<T>::Clear()
{
	_count = 0;
}

template <class T>
inline void ListOfPointers<T>::SetNewCapacity(int new_capacity)
{
	if (new_capacity <= 0)
	{
		Clear();
		return;
	}

	if (new_capacity <= _capacity)
	{
		for (int i = new_capacity; i < _count; i++)
			_items[i]->Reset();

		_count = new_capacity;
		return;
	}

	T **new_items_list;
	int old_capacity = _capacity;

	new_items_list = new T *[new_capacity];
	
	for (int i = 0; i < _count; i++)
			new_items_list[i] = _items[i];

	delete [] _items;

	_capacity = new_capacity;
	_items = new_items_list;
}

//---------------------------------------------------------------------------
// Copy functions
//---------------------------------------------------------------------------
template <class T>
inline void ListOfPointers<T>::CopyTo(ListOfPointers<T> *copy)
{
	copy->Clear();
	copy->SetNewCapacity(_capacity);

	for (int i = 0; i < _count; i++)
		copy->AddNew(_items[i]);
}

//---------------------------------------------------------------------------
// Sort and Search functions
//---------------------------------------------------------------------------
template <class T>
inline bool ListOfPointers<T>::Sort(int start, bool no_case)
{
	bool swap;
	T *a, *b;
	int i;

	start++;
	if (start == 0)
		return false;

	if (start > _count)
		return true;

	do
	{
		swap = false;
		for (i = start; i < _count; i++)
		{
			a = _items[i-1];
			b = _items[i];

			if (a->name.Compare(b->name, no_case) == -1) //-1 => a > b 
			{
				Swap(i-1, i);
				swap = true;
			}
		}
	}
	while(swap);

	return true;
}

template <class T>
inline T *ListOfPointers<T>::Search(String *name, bool no_case, bool partial)
{
	int index;
	return Search(name, index, no_case, partial);
}

template <class T>
inline T *ListOfPointers<T>::Search(char *name, bool no_case, bool partial)
{
	int index;
	return Search(name, index, no_case, partial);
}

template <class T>
inline T *ListOfPointers<T>::Search(String *name, int &index, bool no_case, bool partial)
{
	//ToDo: Make a system to know if list is sorted to use a quicksort search instead a one-by-one comparison
	for (index = 0; index < _count; index++)
		if (_items[index]->name.Compare(name, no_case, partial) == 0)
			return _items[index];

	index = -1;
	return NULL;
}

template <class T>
inline T *ListOfPointers<T>::Search(char *name, int &index, bool no_case, bool partial)
{
	//ToDo: Make a system to know if list is sorted to use a quicksort search instead a one-by-one comparison
	for (index = 0; index < _count; index++)
		if (_items[index]->name.Compare(name, no_case, partial) == 0)
			return _items[index];

	index = -1;
	return NULL;
}

//---------------------------------------------------------------------------
// Internal functions
//---------------------------------------------------------------------------
template <class T>
bool ListOfPointers<T>::_CheckIndex(int index, bool throw_exception)
{
	if (index < 0 || index >= _count)
	{
		if (throw_exception)
		{			
			sprintf(message, "Index (%d) out of bounds. Count on this list: %d.", index, _count);
			throw ExceptionHandler(message);
		}
		else
			return false;
	}

	return true;
}

//OLD Code
//---------------------------------------------------------------------------
// Constructors & Destructors
//---------------------------------------------------------------------------
/*
template <class T>
ListOfPointers<T>::ListOfPointers(bool is_owner)
{
	count    = 0;
	capacity = LIST_INITIAL_CAPACITY;

	owner = is_owner; 

	_ResetPosition();

	items = new T *[capacity];

	for (int i = 0; i < capacity; i++)
		items[i] = NULL;
}

template <class T>
ListOfPointers<T>::ListOfPointers(ListOfPointers<T> *list)
{
	list->CopyTo(this);
}

template <class T>
ListOfPointers<T>::~ListOfPointers()
{
	if (owner)
		for (int i = 0; i < count; i++)
			Delete(i);

	delete [] items;
}

//---------------------------------------------------------------------------
// Operators
//---------------------------------------------------------------------------
template <class T>
T * ListOfPointers<T>::operator[] (const int index)
{
	_CheckIndex(index);
	return items[index];
}

template <class T>
void ListOfPointers<T>::operator++ (void)
{
	p++;
	_CheckPosition();
}

template <class T>
void ListOfPointers<T>::operator-- (void)
{
	if (--p < -1)
	{
		out_of_range = true;	
		p = -1;
	}
	else
		out_of_range = false;
}

//---------------------------------------------------------------------------
// Public Functions
//---------------------------------------------------------------------------
template <class T>
int ListOfPointers<T>::Add(T *new_item)
{
	if (count == capacity)
		Capacity(capacity + LIST_CAPACITY_ADD_STEP);

	items[count] = new_item;
	int index = count++;

	return index;
}

template <class T>
int ListOfPointers<T>::Add(const String &name)
{
	return Add(name.CharPtr());
}

template <class T>
int ListOfPointers<T>::Add(const char *name)
{
	T * new_item = new T;
	new_item->name = name;
	return Add(new_item);
}

template <class T>
T *ListOfPointers<T>::Element(int index)
{
	if (index < 0 || index >= count)
		return NULL;

	return items[index];
}

template <class T>
T *ListOfPointers<T>::AddNew()
{ 
	T *item = NULL; 
	item = new T; 
	Add(item); 
	return item; 
}

template <class T>
bool ListOfPointers<T>::AddNew(int n_items)
{
	for (int i = 0; i < n_items; i++)
		AddNew();

	return true;
}

template <class T>
void ListOfPointers<T>::Clear()
{
	if (owner)
		for (int i = 0; i < count; i++)
			Delete(i);

	count = 0;

	_ResetPosition();
}

template <class T>
void ListOfPointers<T>::Delete(int index, bool compress)
{
	_CheckIndex(index);

	if (owner)
		delete items[index];
	items[index] = NULL;

	if (compress)
		Compress();
}

template <class T>
void ListOfPointers<T>::Compress()
{
	int i, j;

	for (i = 0; i < count; i++)
	{
		if (items[i] == NULL)
		{
			j = i + 1;
			while(j < count && items[j] == NULL) 
				j++;

			if (j < count && items[j] != NULL)
			{
				items[i] = items[j];
				items[j] = NULL;
			}
			else
				break;
		}
	}

	count = i;
	_CheckPosition();
}

template <class T>
void ListOfPointers<T>::Capacity(int new_capacity)
{
	if (new_capacity == capacity)
		return;

	if (new_capacity <= 0)
	{
		Clear();
		return;
	}

	T **new_items;
	int old_capacity = capacity;

	new_items = new T * [new_capacity];

	int i;
	if (count > new_capacity)
	{		
		for (i = 0; i < new_capacity; i++)
			new_items[i] = items[i];

		if (owner)
			for (; i < count; i++)
				Delete(i);

		count = new_capacity;
	}
	else
	{
		for (i = 0; i < count; i++)
			new_items[i] = items[i];

		for (; i < new_capacity; i++)
			new_items[i] = NULL;
	}

	if (capacity > 0)
		delete [] items;

	capacity = new_capacity;
	items = new_items;

	_CheckPosition();
}

template <class T>
void ListOfPointers<T>::Swap(const int i1, const int i2)
{
	_CheckIndex(i1);
	_CheckIndex(i2);

	T *temp = items[i1];
	items[i1] = items[i2];
	items[i2] = temp;
}

template <class T>
T *ListOfPointers<T>::Items(const int index)
{
	return operator[](index);
}

template <class T>
T *ListOfPointers<T>::Item()
{
	if (out_of_range)
		return NULL;

	return operator[](p);
}

template <class T>
void ListOfPointers<T>::MoveFirst() 
{
	if (count > 0) 
	{ 
		out_of_range = false;
		p = 0; 
	}
}

template <class T>
void ListOfPointers<T>::MoveLast() 
{
	if (count > 0) 
	{ 
		p = count - 1; 
		out_of_range = false; 
	}
}

template <class T>
void ListOfPointers<T>::MoveTo(int new_position)
{
	_CheckIndex(new_position);
	p = new_position;
}

template <class T>
T *ListOfPointers<T>::First()
{
	if (count > 0)
		return items[0];

	return NULL;
}

template <class T>
T *ListOfPointers<T>::Last()
{
	if (count > 0)
		return items[count - 1];

	return NULL;
}

template <class T>
void ListOfPointers<T>::CopyTo(ListOfPointers<T> *copy)
{
	copy->Clear();
	copy->Capacity(capacity);
	copy->Owner(owner);

	T *item;
	for (int i = 0; i < count; i++)
	{
		item = items[i];
		if (item == NULL)
			copy->Add((T *)NULL);
		else
			copy->Add(item->Copy());
	}

	copy->ResetPosition();
}

template <class T>
ListOfPointers<T> *ListOfPointers<T>::Copy()
{
	ListOfPointers<T> *new_list;
	new_list = new ListOfPointers<T>;

	CopyTo(new_list);

	return new_list;
}

template <class T>
void ListOfPointers<T>::Fill(T *to_fill_in, int new_capacity)
{
	Clear();
	Capacity(new_capacity);

	for (int i = 0; i < capacity; i++)
		if (to_fill_in == NULL) Add(NULL);
		else
		{
			T *new_i = new T;
			new_i = to_fill_in;
			Add(new_i);
		}
}

template <class T>
bool ListOfPointers<T>::Sort(int start, bool no_case)
{
	bool swap;
	T *a, *b;
	int i;

	start++;
	if (start == 0 || start >= count)
		return false;

	do
	{
		swap = false;
		for (i = start; i < count; i++)
		{
			a = (T *)items[i-1];
			b = (T *)items[i];

			if (a->name.Compare(b->name, no_case) == -1) //-1 => a > b 
			{
				Swap(i-1, i);
				swap = true;
			}
		}
	}
	while(swap);

	return true;
}

template <class T>
bool ListOfPointers<T>::SortItself(int start = 0, bool no_case = true)
{
	bool swap;
	T *a, *b;
	int i;

	start++;
	if (start == 0 || start >= count)
		return false;

	do
	{
		swap = false;
		for (i = start; i < count; i++)
		{
			a = (T *)items[i-1];
			b = (T *)items[i];

			if (a->Compare(b, no_case) == -1) //-1 => a > b 
			{
				Swap(i-1, i);
				swap = true;
			}
		}
	}
	while(swap);

	return true;
}

template <class T>
T *ListOfPointers<T>::Search(const String &name, bool no_case, bool partial)
{
	int index;
	return Search(name.CharPtr(), index, no_case, partial);
}

template <class T>
T *ListOfPointers<T>::Search(const char *name, bool no_case, bool partial)
{
	int index;
	return Search(name, index, no_case, partial);
}

template <class T>
T *ListOfPointers<T>::Search(const String &name, int &index, bool no_case, bool partial)
{
	return Search(name.CharPtr(), index, no_case, partial);
}

template <class T>
T *ListOfPointers<T>::Search(const char *name, int &index, bool no_case, bool partial)
{
	T *item;

	for (index = 0; index < count; index++)
	{
		item = (T *)items[index];
		if (item != NULL && item->name.Compare(name, no_case, partial) == 0)
			return item;
	}

	index = -1;
	return NULL;
}

template <class T>
T *ListOfPointers<T>::Search(T *to_find)
{
	for (int index = 0; index < count; index++)
		if (*(items[index]) == *to_find)
			return items[index];

	return NULL;
}

template <class T>
T *ListOfPointers<T>::Search(T *to_find, int &index)
{
	for (index = 0; index < count; index++)
		if (*(items[index]) == *to_find)
			return items[index];

	index = -1;
	return NULL;
}

template <class T>
T *ListOfPointers<T>::Store(const String &name, bool replace)
{
	int index;
	T *item;

	item = Search(name);

	if (item != NULL && replace)
	{
		item->Reset();
		item->name = name;
	}
	else if (item == NULL)
	{
		index = Add(name);
		item = items[index];
	}
		
	return item;
}

template <class T>
T *ListOfPointers<T>::Store(const char *name, bool replace)
{
	String n(name);
	return Store(n, replace);
}

template <class T>
bool ListOfPointers<T>::Store(const int index, T * item)
{
	if (!_CheckIndex(index, false)) return false;

	if (owner && items[index] != NULL)
		Delete(index);

	items[index] = item;
	return true;
}

//---------------------------------------------------------------------------
// Private Functions
//---------------------------------------------------------------------------
template <class T>
bool ListOfPointers<T>::_CheckIndex(int index, bool throw_exception)
{
	if (index >= 0 && index < count)
		return true;
	else if (throw_exception)
		throw EListOfPointersIndexOutOfRange();

	return false;
}

template <class T>
void ListOfPointers<T>::_ResetPosition()
{
	if (count > 0)
	{
		p = 0;
		out_of_range = false;
	}
	else
	{
		p = -1;
		out_of_range = true;
	}
}

template <class T>
void ListOfPointers<T>::_CheckPosition()
{
	if (p > count)
	{
		out_of_range = true;
		p = count;
	}
	else if (p < 0)
	{
		out_of_range = true;
		p = -1;
	}
	else
		out_of_range = false;
}
*/
#endif
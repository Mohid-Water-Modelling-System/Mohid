#ifndef llistH
#define llistH

#include "string.h"
#include "exceptionhandler.h"

template <class T>
class ListItem
{
public:
	ListItem()
	{
		try
		{
			data = NULL;
			data = new T;

			n = NULL;
			p = NULL;
		}
		catch(...)
		{
			printf("\nPhreeqC Object: Erro de alocação de memória para um item de lista...\n");
			throw;
		}
	}
	~ListItem()
	{
		delete data;
	}

public:
	T *data;
	ListItem *n, *p;
	
	int index;
};

template<class T>
class List
{
public:
	List(int initial_capacity = 1);
	~List();

public:
	inline T *operator[] (const int index) { return Element(index); }

public:
	inline int Count() { return _count; }
	inline int Capacity() { return _capacity; }
	inline int LastIndex() { return _count - 1; }
	inline bool IsEmpty() { if (_count <= 0) return true; return false; }
	inline T *Element(int index) 
	{
		if (!CheckIndex(index))
			return NULL;

		return _index[index]->data;
	}
	inline T *Items(int index) { return Element(index); }
	inline T *Last() { return Element(_count - 1); }
	inline T *First() { return Element(0); }

public:
	T *AddNew();
	T *AddNew(String *name);
	T *AddNew(char *name);
	T *AddNew(T *to_copy);
	//T *AddNew(int &index);
	T *AddNew(int &index, String *name);
	T *AddNew(int &index, char *name);
	T *AddNew(int &index, T *to_copy);
	bool AddNew(int n_items);
	T *Store(String *name, bool replace = true);

	/*
		The Delete function must be used ONLY with the control of instances in the "interface.cpp" code
		Uses the ->index in the ListItem to find the correct "item" to delete, so, do not work with "ordered" lists
	*/
	bool Delete(int index);

public:
	bool Sort(int start = 0, bool no_case = true);
	T *Search(String *name, bool no_case = false, bool partial = false);
	T *Search(char *name, bool no_case = false, bool partial = false);
	T *Search(String *name, int &index, bool no_case = false, bool partial = false);
	T *Search(char *name, int &index, bool no_case = false, bool partial = false);	

public:
	void Clear();
	void Swap(int i1, int i2);
	void CopyTo(List<T> *copy);
	void SetNewCapacity(int new_capacity); //Reset the _index matrix
	void CompressByIndex(int new_count);

protected:
	ListItem<T> *_first,
							 *_last,
							 *_next,
							 *_temp,
							 **_index;

	int _count,
			_capacity;

private:
	void AddListItem();
	T *GetListItem();
	T *GetListItem(int &index);
	ListItem<T> *GetItem(int index);
	void DeleteList();
	void ResetIndex();

private:
	bool CheckIndex(int i);

};

//==========================================================================================================
// Constructors & Destructors
//==========================================================================================================
template <class T>
List<T>::List(int initial_capacity)
{
	_first = NULL;
	_last = NULL;
	_next = NULL;
	_temp = NULL;
	_index = NULL;

	_count = 0;
	_capacity = 0;

	if (initial_capacity <= 0)
		initial_capacity = 1;
		
	SetNewCapacity(initial_capacity);
}

template <class T>
List<T>::~List()
{
	DeleteList();
	delete [] _index;
}

//==========================================================================================================
// Add & Store 
//==========================================================================================================
template <class T>
T *List<T>::AddNew()
{
	return GetListItem();
}

template <class T>
T *List<T>::AddNew(T *to_copy)
{
	T *new_item = GetListItem();
	to_copy->CopyTo(new_item);
	return new_item;
}

template <class T>
T *List<T>::AddNew(String *name)
{
	T *new_item = GetListItem();
	new_item->name = name;
	return new_item;
}

template <class T>
T *List<T>::AddNew(char *name)
{
	T *new_item = GetListItem();
	new_item->name = name;
	return new_item;
}

/*
template <class T>
T *List<T>::AddNew(int &index)
{
	return GetListItem(index);
}
*/

template <class T>
T *List<T>::AddNew(int &index, String *name)
{
	T *new_item = GetListItem(index);
	new_item->name = name;
	return new_item;
}

template <class T>
T *List<T>::AddNew(int &index, char *name)
{
	T *new_item = GetListItem(index);
	new_item->name = name;
	return new_item;
}

template <class T>
T *List<T>::AddNew(int &index, T *to_copy)
{
	T *new_item = GetListItem(index);
	to_copy->CopyTo(new_item);
	return new_item;
}

template <class T>
bool List<T>::AddNew(int n_items)
{
	for (int i = 0; i < n_items; i++)
		GetListItem();

	return true;
}

template <class T>
T *List<T>::Store(String *name, bool replace = true)
{
	T *item = Search(name);

	if (item != NULL && replace)
	{
		item->Reset();
		item->name = name;
	}
	else if (item == NULL)
	{
		item = AddNew(name);
	}
		
	return item;
}

template <class T>
bool List<T>::Delete(int index)
{
	if (_count > 1)
	{
		for (int i = 0; i < _count; i++)
		{
			if (_index[i]->index == index)
			{
				//Store item pointer
				_temp = _index[i];

				int j;
				for (j = i; j < _capacity - 1; j++)
					_index[j] = _index[j + 1];
				_index[j] = _temp;

				if (_temp != _last)
				{
					if (_temp == _first)
					{
						_first = _temp->n;
						_first->p = NULL;
					}

					_temp->p = _last;
					_last->n = _temp;
					_last = _temp;
					_temp->n = NULL;
				}

				break;
			}
		}
	}
	else
	{
		_next = _first;
	}

	_count--;
	return true;
}
//==========================================================================================================
// Search & Sort
//==========================================================================================================
template <class T>
bool List<T>::Sort(int start, bool no_case)
{
	bool swap;
	T *a, *b;

	start++;
	if (start <= 0)
		return false;

	if (start >= _count)
		return true;

	int i;
	do
	{
		swap = false;

		for (i = start; i < _count; i++)
		{
			a = _index[i-1]->data;
			b = _index[i]->data;

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
inline T *List<T>::Search(String *name, bool no_case, bool partial)
{
	int index;
	return Search(name, index, no_case, partial);
}

template <class T>
inline T *List<T>::Search(char *name, bool no_case, bool partial)
{
	int index;
	return Search(name, index, no_case, partial);
}

template <class T>
inline T *List<T>::Search(String *name, int &index, bool no_case, bool partial)
{
	//ToDo: Make a system to know if list is sorted to use a quicksort search instead a one-by-one comparison
	for (index = 0; index < _count; index++)
		if (_index[index]->data->name.Compare(name, no_case, partial) == 0)
			return _index[index]->data;

	index = -1;
	return NULL;
}

template <class T>
inline T *List<T>::Search(char *name, int &index, bool no_case, bool partial)
{
	//ToDo: Make a system to know if list is sorted to use a quicksort search instead a one-by-one comparison
	if (name != NULL) 
	{
		for (index = 0; index < _count; index++)
			if (_index[index]->data->name.Compare(name, no_case, partial) == 0)
				return _index[index]->data;
	}

	index = -1;
	return NULL;
}

//==========================================================================================================
// Utilities
//==========================================================================================================
template <class T>
void List<T>::Clear()
{
	//reset "index" matrix if _count is greater than 0
	ResetIndex();

	//reset _count
	_count = 0;

	//reset _next pointer
	_next = _first;
}

template <class T>
void List<T>::Swap(int i1, int i2)
{
	if (!CheckIndex(i1) || !CheckIndex(i2))
	{
		char message[500];
		sprintf(message, "Index out of bounds.");
		throw ExceptionHandler(message);
	}

	_temp = _index[i1];
	_index[i1] = _index[i2];
	_index[i2] = _temp;
}

template <class T>
void List<T>::CopyTo(List<T> *copy)
{
	int i, index;

	copy->Clear();
	copy->SetNewCapacity(_capacity);
	copy->ResetIndex();

	_temp = _first;
	for (i = 0; i < _count; i++)
	{
		copy->AddNew(_temp->data);
		_temp = _temp->n;
	}

	for (i = 0; i < _count; i++)
	{
		index = _index[i]->index;
		_temp = copy->GetItem(index);
		copy->_index[i] = _temp;
	}
}

template <class T>
void List<T>::CompressByIndex(int new_count)
{
	if (new_count < 0 || new_count > _count)
		throw ExceptionHandler("Index out of bounds.");

	if (new_count == _count)
		return;

	_first = _index[0];
	_temp = _first;
	_temp->index = 0;
	for (int i = 1; i < _capacity; i++)
	{
		_temp->n = _index[i];
		_temp = _temp->n;
		_temp->index = i;
	}
	_temp->n = NULL;

	_count = new_count;
	_next = _index[_count];
}

//==========================================================================================================
// List control
//==========================================================================================================
template <class T>
void List<T>::AddListItem()
{
	ListItem<T> *new_item = new ListItem<T>;

	new_item->n = NULL;
	new_item->p = NULL;
	new_item->index = _capacity++;

	if (_first == NULL)
		_first = new_item;		
	else
		_last->n = new_item;

	new_item->p = _last;
	_last = new_item;
}

template <class T>
T *List<T>::GetListItem()
{
	T *item = _next->data;

	_count++;

	if (_count == _capacity)
		SetNewCapacity(_capacity + 1);
	else
		_next = _index[_count];

	return item;
}

template <class T>
T *List<T>::GetListItem(int &index)
{
	T *item = _next->data;
	index = _count;

	_count++;

	if (_count == _capacity)
		SetNewCapacity(_capacity + 10);
	else
		_next = _index[_count];

	return item;
}

template <class T>
ListItem<T> *List<T>::GetItem(int index)
{
	if (!CheckIndex(index))
		return NULL;

	_temp = _first;
	for (int i = 0; i < index; i++)
		_temp = _temp->n;

	return _temp;
}

template <class T>
void List<T>::DeleteList()
{
	if (_count > 0)
		while (_first != NULL)
		{
			_temp = _first->n;
			delete _first;
			_first = _temp;
		}

	_first = NULL;
	_last = NULL;
}

template <class T>
void List<T>::ResetIndex()
{
	int i;

	_temp = _first;

	for (i = 0; i < _capacity; i++)
	{			
		_index[i] = _temp;
		_temp = _temp->n;
	}
}

template <class T>
void List<T>::SetNewCapacity(int new_capacity)
{
	if (new_capacity == _capacity)
		return;

	if (new_capacity <= 0)
	{
		Clear();
		return;
	}

	if (new_capacity < _capacity)
	{
		//Set "_count" in case new_capacity is lesser than the number of "slots" that were in use
		if (_count > new_capacity)
			_count = new_capacity;

		//reset the index matrix case _count > 0
		ResetIndex();

		//Correct "next" pointer
		_next = _index[_count];

		return;
	}

	int i;
	int old_capacity = _capacity;

	//Store the position of the actual LAST element of the list
	_temp = _last;

	//Add (new_capacity - _capacity) elements to the list
	for (i = 0; i < (new_capacity - old_capacity); i++)
		AddListItem(); //AddListItem change _capacity (+1)

	//Creates a new index matrix
	ListItem<T> **new_index = new ListItem<T> * [new_capacity];

	//Pass the actual "order" of the index to the new index matrix
	for (i = 0; i < old_capacity; i++)
		new_index[i] = _index[i];

	//delete the actual index matrix and associates the index pointer to the new allocated index matrix
	delete [] _index;
	_index = new_index;

	//Insert a pointer to the new elements on the list in the end of index matrix 
	//(after the position of the old last element)
	if (_temp != NULL)
	{
		for (i = old_capacity; i < new_capacity; i++)
		{
			_temp = _temp->n;
			_index[i] = _temp;
		}
	}
	else
		ResetIndex();

	//Correct "next" pointer
	_next = _index[_count];
}

//==========================================================================================================
// Internal
//==========================================================================================================
template <class T>
bool List<T>::CheckIndex(int i)
{
	if (i < 0)
		return false;

	if (i >= _count)
		return false;

	return true;
}

#endif
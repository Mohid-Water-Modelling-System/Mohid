#ifndef arrayH
#define arrayH

#include <algorithm>
#include <stdio.h>

template<class T>
class Array
{
public:
	Array();
	~Array();

public:
	void Capacity(int new_capacity);
	int Capacity() { return capacity; }
	int Count() { return count; }
	void SetPointerToArray(T **pointer);
	T *PointerToArray() { return array_; }
	void Add(T *data);
	void Fill(T *value_to_fill, int new_capacity);
	void PrintToFile(String filename);

private:
	int capacity;
	int count;
	T *array_;
	T **pointer_to_array;
};

template <class T>
Array<T>::Array()
{
	capacity = 0;
	count = 0;
	pointer_to_array = new T *;

	array_ = NULL;
	*pointer_to_array = array_;
}

template <class T>
Array<T>::~Array()
{
	if (array_ != NULL) delete [] array_;
}

template <class T>
void Array<T>::Capacity(int new_capacity)
{
	if (new_capacity <= 0)
		return;

	T *t_array = new T [new_capacity];

	int max = min(capacity, new_capacity);
	for (int i = 0; i < max; i++)
		t_array[i] = array_[i];

	delete [] array_;

	array_ = t_array;
	*pointer_to_array = array_;
	capacity = new_capacity;
	count = min(count, capacity);
}

template <class T>
void Array<T>::SetPointerToArray(T **pointer)
{
	pointer_to_array = pointer;
}

template <class T>
void Array<T>::Add(T *data)
{
	if (count >= capacity)
		Capacity(capacity + 20);

	array_[count] = *data;
	count++;
}

template <class T>
void Array<T>::Fill(T *value_to_fill, int new_capacity)
{
	Capacity(new_capacity);

	for (int i = 0; i < capacity; i++)
		array_[i] = *value_to_fill;
}

template <class T>
void Array<T>::PrintToFile(String filename)
{
	FILE *	f = fopen(filename.CharPtr(), "wt");

	for (int i = 0; i < capacity; i++)
		fprintf(f, "%f\n", array_[i]);

	fclose(f);
}

#endif
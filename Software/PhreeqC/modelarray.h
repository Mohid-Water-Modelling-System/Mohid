#ifndef modelarrayH
#define modelarrayH

#include <stdio.h>
#include "constants.h"


class ModelArray
{
public:
	ModelArray();
	~ModelArray();

public:
	void ArrNewMax(int new_max);
	void ResidualNewMax(int new_max);
	void DeltaNewMax(int new_max);
	void NormalNewMax(int new_max);
	void IneqArrayNewMax(int new_max);
	void BackEqNewMax(int new_max);
	void ZeroNewMax(int new_max);
	void ResNewMax(int new_max);
	void Delta1NewMax(int new_max);
	void CuNewMax(int new_max);
	void IuNewMax(int new_max);
	void IsNewMax(int new_max);
	void XArgENewMax(int new_max);
	void ResArgENewMax(int new_max);
	void ScratchENewMax(int new_max);

	void Fill(LDBLE *_array, LDBLE value, int count);
	void Fill(int *_array, int value, int count);

public:
	inline int arr_max() { return _arr_max; }
	inline int residual_max() { return _residual_max; }
	inline int delta_max() { return _delta_max; }
	inline int ineq_array_max() { return _ineq_array_max; }
	inline int normal_max() { return _normal_max; }
	inline int back_eq_max() { return _back_eq_max; }
	inline int zero_max() { return _zero_max; }
	inline int res_max() { return _res_max; }
	inline int delta1_max() { return _delta1_max; }
	inline int cu_max() { return _cu_max; }
	inline int iu_max() { return _iu_max; }
	inline int is_max() { return _is_max; }
	inline int x_arg_e_max() { return _x_arg_e_max; } 
	inline int res_arg_e_max() { return _res_arg_e_max; } 
	inline int scratch_e_max() { return _scratch_e_max; } 


protected:
	LDBLE *arr,
				*residual,
				*delta;

	LDBLE *normal,
				*ineq_array,
				*zero, 
				*res, 
				*delta1,
				*cu,
				*x_arg_e,
				*res_arg_e,
				*scratch_e;

	int *back_eq,
		  *iu,
			*is;

private:
	int _arr_max,
			_residual_max,
			_delta_max;

	int _normal_max,
			_ineq_array_max,
			_zero_max, 
			_res_max, 
			_delta1_max,
			_cu_max,
			_x_arg_e_max,
			_res_arg_e_max,
			_scratch_e_max;

	int _back_eq_max,
			_iu_max,
			_is_max;

	int _arr_capacity,
			_residual_capacity,
			_delta_capacity;

	int _normal_capacity,
			_ineq_array_capacity,
			_zero_capacity, 
			_res_capacity, 
			_delta1_capacity,
			_cu_capacity,
			_x_arg_e_capacity,
			_res_arg_e_capacity,
			_scratch_e_capacity;

	int _back_eq_capacity,
			_iu_capacity,
			_is_capacity;
};


inline 
void ModelArray::ArrNewMax(int new_max)
{
	_arr_max = new_max;

	if (_arr_max > _arr_capacity)
	{
		_arr_capacity = _arr_max;

		delete [] arr;
		arr = NULL;
		arr = new LDBLE [_arr_capacity];
	}
}

inline 
void ModelArray::ResidualNewMax(int new_max)
{
	_residual_max = new_max;

	if (_residual_max > _residual_capacity)
	{
		_residual_capacity = _residual_max;

		delete [] residual;
		residual = NULL;
		residual = new LDBLE [_residual_capacity];
	}
}

inline 
void ModelArray::DeltaNewMax(int new_max)
{
	_delta_max = new_max;

	if (_delta_max > _delta_capacity)
	{
		_delta_capacity = _delta_max;

		delete [] delta;
		delta = NULL;
		delta = new LDBLE [_delta_capacity];
	}
}

inline
void ModelArray::NormalNewMax(int new_max)
{
	_normal_max = new_max;

	if (_normal_max > _normal_capacity)
	{
		_normal_capacity = _normal_max;

		delete [] normal;
		normal = NULL;
		normal = new LDBLE [_normal_capacity];
	}
}

inline
void ModelArray::IneqArrayNewMax(int new_max)
{
	_ineq_array_max = new_max;

	if (_ineq_array_max > _ineq_array_capacity)
	{
		_ineq_array_capacity = _ineq_array_max;

		delete [] ineq_array;
		ineq_array = NULL;
		ineq_array = new LDBLE [_ineq_array_capacity];
	}
}

inline
void ModelArray::BackEqNewMax(int new_max)
{
	_back_eq_max = new_max;

	if (_back_eq_max > _back_eq_capacity)
	{
		_back_eq_capacity = _back_eq_max;

		delete [] back_eq;
		back_eq = NULL;
		back_eq = new int [_back_eq_capacity];
	}
}

inline
void ModelArray::ZeroNewMax(int new_max)
{
	_zero_max = new_max;

	if (_zero_max > _zero_capacity)
	{
		_zero_capacity = _zero_max;

		delete [] zero;
		zero = NULL;
		zero = new LDBLE [_zero_capacity];
	}
}

inline
void ModelArray::ResNewMax(int new_max)
{
	_res_max = new_max;

	if (_res_max > _res_capacity)
	{
		_res_capacity = _res_max;

		delete [] res;
		res = NULL;
		res = new LDBLE [_res_capacity];
	}
}

inline
void ModelArray::Delta1NewMax(int new_max)
{
	_delta1_max = new_max;

	if (_delta1_max > _delta1_capacity)
	{
		_delta1_capacity = _delta1_max;

		delete [] delta1;
		delta1 = NULL;
		delta1 = new LDBLE [_delta1_capacity];
	}
}

inline
void ModelArray::CuNewMax(int new_max)
{
	_cu_max = new_max;

	if (_cu_max > _cu_capacity)
	{
		_cu_capacity = _cu_max;

		delete [] cu;
		cu = NULL;
		cu = new LDBLE [_cu_capacity];
	}
}

inline
void ModelArray::IuNewMax(int new_max)
{
	_iu_max = new_max;

	if (_iu_max > _iu_capacity)
	{
		_iu_capacity = _iu_max;

		delete [] iu;
		iu = NULL;
		iu = new int [_iu_capacity];
	}
}

inline
void ModelArray::IsNewMax(int new_max)
{
	_is_max = new_max;

	if (_is_max > _is_capacity)
	{
		_is_capacity = _is_max;

		delete [] is;
		is = NULL;
		is = new int [_is_capacity];
	}
}

inline
void ModelArray::XArgENewMax(int new_max)
{
	_x_arg_e_max = new_max;

	if (_x_arg_e_max > _x_arg_e_capacity)
	{
		_x_arg_e_capacity = _x_arg_e_max;

		delete [] x_arg_e;
		x_arg_e = NULL;
		x_arg_e = new LDBLE [_x_arg_e_capacity];
	}
}

inline
void ModelArray::ResArgENewMax(int new_max)
{
	_res_arg_e_max = new_max;

	if (_res_arg_e_max > _res_arg_e_capacity)
	{
		_res_arg_e_capacity = _res_arg_e_max;

		delete [] res_arg_e;
		res_arg_e = NULL;
		res_arg_e = new LDBLE [_res_arg_e_capacity];
	}
}

inline
void ModelArray::ScratchENewMax(int new_max)
{
	_scratch_e_max = new_max;

	if (_scratch_e_max > _scratch_e_capacity)
	{
		_scratch_e_capacity = _scratch_e_max;

		delete [] scratch_e;
		scratch_e = NULL;
		scratch_e = new LDBLE [_scratch_e_capacity];
	}
}

inline 
void ModelArray::Fill(LDBLE *_array, LDBLE value, int count)
{
	for (int i = 0; i < count; i++)
		_array[i] = value;
}

inline 
void ModelArray::Fill(int *_array, int value, int count)
{
	for (int i = 0; i < count; i++)
		_array[i] = value;
}

#endif
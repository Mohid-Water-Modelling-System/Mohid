#include "modelarray.h"

ModelArray::ModelArray()
{
	_arr_max = 0;
	_arr_capacity = 0;
	arr = NULL;

	_residual_max = 0;
	_residual_capacity = 0;
	residual = NULL;

	_delta_max = 0;
	_delta_capacity = 0;
	delta = NULL;

	_normal_max = 0;
	_normal_capacity = 0;
	normal = NULL;

	_ineq_array_max = 0;
	_ineq_array_capacity = 0;
	ineq_array = NULL;

	_back_eq_max = 0;
	_back_eq_capacity = 0;
	back_eq = NULL;

	_zero_max = 0;
	_zero_capacity = 0;
	zero = NULL;

	_res_max = 0;
	_res_capacity = 0;
	res = NULL;

	_delta1_max = 0;
	_delta1_capacity = 0;
	delta1 = NULL;

	_cu_max = 0;
	_cu_capacity = 0;
	cu = NULL;

	_iu_max = 0;
	_iu_capacity = 0;
	iu = NULL;

	_is_max = 0;
	_is_capacity = 0;
	is = NULL;

	_x_arg_e_max = 0;
	_x_arg_e_capacity = 0;
	x_arg_e = NULL;

	_res_arg_e_max = 0;
	_res_arg_e_capacity = 0;
	res_arg_e = NULL;

	_scratch_e_max = 0;
	_scratch_e_capacity = 0;
	scratch_e = NULL;
}

ModelArray::~ModelArray()
{
	delete arr;
	delete residual;
	delete delta;
	delete normal;
	delete ineq_array;
	delete back_eq;
	delete zero;
	delete res;
	delete delta1;
	delete cu;
	delete iu;
	delete is;
	delete x_arg_e;
	delete res_arg_e;
	delete scratch_e;
}
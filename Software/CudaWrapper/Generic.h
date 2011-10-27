#pragma once

class Generic
{
public:
	/// Print a value of any type to cout
	/// Inline because it's a template (export on templates) is not supported in VS2010)
	template <class Type> inline static void Print(Type val)
	{
		cout << "CudaWrapper " << val << endl;
	}
	
	/// Catch an exception, clean up CUDA. Exit after user input
	static void CatchException(exception e);
};
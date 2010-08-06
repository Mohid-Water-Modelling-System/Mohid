#ifndef useH
#define useH

#include "datastructures.h"

//---------------------------------------------------------------------------
// class Use
//---------------------------------------------------------------------------
class Use
{
public:
	Use(); 
	~Use();

public:
	void Reset();
	void ResetPointers();
	void ClearLists();
	bool CheckIn();

public:
	Solution *SearchSolution(int number);
	Exchange *SearchExchange(int number);
	Surface *SearchSurface(int number);
	
public:
	//these must be provided BEFORE SetUse() is called
	int sol_number,
			ppa_number,
			exc_number,
			sur_number,
			gas_number,
			ssa_number,
			irr_number;

	bool ppa_in,
			 exc_in,
			 sol_in,
			 sur_in,
			 gas_in,
			 ssa_in,
			 irr_in;

	class Solution *sol_p;
	class PPAssemblage *ppa_p;
	class Exchange *exc_p;
	class Surface *sur_p;
	class GasPhase *gas_p;
	class SSAssemblage *ssa_p;
	//class Irrev *irr_p;

	List<class Solution> *sol_list;
	List<class PPAssemblage> *ppa_list;
	List<class Exchange> *exc_list;
	List<class Surface> *sur_list;
	List<class GasPhase> *gas_list;
	List<class SSAssemblage> *ssa_list;
	//List<class Irrev> *irr_list;
};


#endif
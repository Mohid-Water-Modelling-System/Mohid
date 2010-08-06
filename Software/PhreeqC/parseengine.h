#ifndef parseengineH
#define parseengineH

#include "basicfunctions.h"
#include "datastructures.h"

class ParseEngine:
	public BasicFunctions
{
public:
	ParseEngine(GlobalData *gd);
	~ParseEngine();

public:
	bool GetSpecies(char **str);
	bool ParseEq(String &eq, bool association = true);
	bool GetElementsInSpecies(char **str, LDBLE coef);
	bool GetSecondaryElementsInSpecies(char **str, LDBLE &coef);
	bool ComputeGFW(String &str, LDBLE &gfw);
	void CopyToTempEOSList(List<class ElementOfSpecies> *list, LDBLE coef); //add_elt_list
	void SaveEOSList(List<class ElementOfSpecies> *list);

protected:
	bool GetCoef(char **str, LDBLE &coef);
	bool CombineElements(void); //elt_list_combine
	ElementOfSpecies *AddElementToList(String &element);
	void SetGlobalData(GlobalData *gd) { this->gd = gd; }

protected:
	int parent_count;
	List<ElementOfSpecies> *eos_list;
	Reaction rt;

protected:
	GlobalData *gd;

private:
	List<ElementOfSpecies> *t_eos_list, *eos_list_ptr;

};

#endif
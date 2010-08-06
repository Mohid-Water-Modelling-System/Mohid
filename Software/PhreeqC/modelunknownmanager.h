#pragma once
#ifndef modelunknownmanagerH
#define modelunknownmanagerH

//#include "unknownlist.h"
#include "modeldata.h"
#include "exceptionhandler.h"

class ModelUnknownManager
{
public:
	ModelUnknownManager();
	~ModelUnknownManager();

public:
	void SetData(ModelData *md_ptr);
	void SetupUnknowns();

protected:
	List<Unknown> *unknown_list;
	//Unknown *unknown_list_ptr;

	int count_unknowns;

private:
	ModelData *md;

	char message[2048];
};

#endif

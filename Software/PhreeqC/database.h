#ifndef databaseH
#define databaseH

#include "definitions.h"
#include "constants.h"
#include "databaseengine.h"
#include "file.h"
#include "debug.h"

class Database:
	public DBEngine
{
public:
	Database(GlobalData *gd);
	~Database();

public:
	bool Load(String &file_name);
	bool TidyDatabase();

private:
	bool LoadSolutionMasterSpeciesBlock();
	bool LoadSolutionSpeciesBlock();
	bool LoadSolutionSpeciesOption();
	bool LoadSolutionSpeciesEquation();
	bool LoadPhasesBlock();
	bool LoadPhaseOption();
	bool LoadPhaseEquation();
	bool LoadExchangeMasterSpeciesBlock();
	bool LoadExchangeSpeciesBlock();
	bool LoadExchangeSpeciesOption();
	bool LoadExchangeSpeciesEquation();
	bool LoadSurfaceMasterSpeciesBlock();
	bool LoadSurfaceSpeciesBlock();
	bool LoadSurfaceSpeciesOption();
	bool LoadSurfaceSpeciesEquation();

private:
	File *file;
	List<Keyword> options_list;
	Species *s;
	Phase *p;
	Debug d;

private: //for debug purposes
	void PrintToFileSpeciesInfo(const char *file_name);
	void PrintToFileSpeciesInfo2(const char *file_name);
	void PrintToFilePhasesInfo(const char *file_name);
};

#endif
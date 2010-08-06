#ifndef modelH
#define modelH

#include "modeldata.h"
#include "database.h"
#include "basicfunctions.h"
#include "ioengine.h"
#include "modelengine.h"
#include "tidy.h"
#include "parseengine.h"
#include "errors.h"
#include "messages.h"

struct ModelOptions
{
	int print_always;
};

class Model:
	public ParseEngine,
	public IOEngine,
	public Messages
{
public:
	Model();
	Model(Model *copy);
	~Model();

public: //List related functions
	void Reset() {}
	void CopyTo(Model *copy) {}

	String name;
	int number;

public:
	bool LoadDatabase(String file_name);

public:
	void SetOptions(int option, bool state);
	bool SetupModel();
	bool RunModel(bool print_input, bool print_result);

	void PrintUserInputToFile();
	void PrintUserInputToFile(FILE *file);

	bool PrintResults(String filename);
	void SetErrorFile(String error_file) 
	{
		error = new ErrorFile(error_file);
		me->error = error;
	}

private: 
	bool RunInitialSolutionCalculations();
	bool RunInitialExchangersCalculations();
	bool RunInitialGasPhaseCalculations();
	bool RunInitialSurfacesCalculations();
	bool RunReactions();

	bool ResetModel();
	void SaveDataStructures();
	bool SetUse();
	void SetupSurface();
	void SetupExchange();

public:
	Solution *solution_r;

public:
	Database *db;
	ModelData *md;
	GlobalData *gd;

	ModelEngine *me;

	Tidy *t;

	int default_pe;

public:
	bool exc_exists,
			 sur_exists,
			 gas_exists,
			 ppa_exists,
			 ssa_exists;

	Solution *solution;
	Exchange *exchange;
	Surface *surface;
	GasPhase *gasphase;
	PPAssemblage *ppassemblage;
	SSAssemblage *ssassemblage;

private:
	ErrorFile *error;
	ModelOptions *options;

	FILE *output_results;
};

#endif
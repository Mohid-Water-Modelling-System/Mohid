#ifndef debugH
#define debugH

#include <stdio.h>

#include "definitions.h"
#include "list.h"
#include "datastructures.h"
#include "array.h"
#include "specieslist.h"
#include "modeldata.h"

class Debug
{
public:
	Debug();
	~Debug();

public:
	void OpenSpeciesReactionFile(String filename);
	void CloseSpeciesReactionFile();
	void PrintSpeciesReaction(Species *s);
	void PrintMasterList(String filename, GlobalData *gd);
	void PrintSpeciesList(String filename, GlobalData *gd);
	void PrintUnknownInfoList(String filename, List<UnknownInfo> *ui);
	void PrintUnknownList(String filename, GlobalData *gd);
	void PrintArray(String filename, LDBLE *array_, int count);
	void PrintPhaseList(String filename, List<Phase> *p);
	void OpenUnknownFile(String filename);
	void CloseUnknownFile();
	void PrintUnknownToFile(List<Unknown> *x, String comment);
	void OpenArrayFile(String filename);
	void CloseArrayFile();
	void PrintArrayToFile(LDBLE *array_, int count, String comment);
	void PrintSpeciesInfoListToFile(String filename, SpeciesList *sl);

	void OpenSpeciesInfoFile(String filename);
	void CommentIntoSpeciesInfoFile(String comment);
	void PrintSpeciesInfoToFile(String comment, SpeciesList *si_list, ModelData *md, GlobalData *gd);
	void CloseSpeciesInfoFile();
	LDBLE Under(LDBLE xval, LDBLE LOG_10);

	void OpenGammasFile();
	void CommentIntoGammasFile(String comment);
	void PrintGammasToFile(String comment, List<Species> *list);
	void CloseGammasFile();

	void OpenDeltaFile();
	void CommentIntoDeltaFile(String comment);
	void PrintDeltaToFile(String comment, LDBLE *delta, int count);
	void CloseDeltaFile();

	void OpenLAFile();
	void CommentIntoLAFile(String comment);
	void PrintLAToFile(String comment, int count, List<Unknown> *x);
	void CloseLAFile();

	void OpenDebugFile();
	void CloseDebugFile();

	void OpenSetFile();
	void CloseSetFile();
	void OpenReviseGuessesFile();
	void CloseReviseGuessesFile();
	void OpenBuildModelFile();
	void CloseBuildModelFile();
	void OpenMolalitiesFile();
	void CloseMolalitiesFile();
	void OpenMBSumsFile();
	void CloseMBSumsFile();
	void OpenBuildSSAssemblageFile();
	void CloseBuildSSAssemblageFile();
	void OpenBuildJacobianSumsFile();
	void CloseBuildJacobianSumsFile();
	void OpenIneqFile();
	void CloseIneqFile();
	void OpenResetFile();
	void CloseResetFile();

	void OpenSetupSolutionFile();
	void CloseSetupSolutionFile();


	//======================================================================================================
public:
	void OpenAddSolutionFile();
	void CloseAddSolutionFile();
	void OpenAddExchangeFile();
	void CloseAddExchangeFile();
	void OpenAddPPAssemblageFile();
	void CloseAddPPAssemblageFile();

public:
	FILE *addsolution_f,
		   *addexchange_f,
			 *addppassemblage_f;

	int addsolution_count,
			addexchange_count,
			addppassemblage_count;

public:
	int si_count;
	int gammas_count;
	int delta_count;
	int arr_count, la_count;

	int debug_count;

	int set_count,
			reviseguesses_count,
			buildmodel_count,
			molalities_count,
			mbsums_count,
			buildssassemblage_count,
			buildjacobiansums_count,
			ineq_count,
			reset_count,
			setup_solution_count;

	bool debug_status;

	FILE *debug_file;

	FILE *srxn_f, 
		   *u_f,
			 *arr_f,
			 *si_f,			 
			 *delta_f, 
			 *la_f;

	FILE *set_f,
			 *gammas_f,
			 *reviseguesses_f,
			 *buildmodel_f,
			 *molalities_f,
			 *mbsums_f,
			 *buildssassemblage_f,
			 *buildjacobiansums_f,
			 *ineq_f,
			 *reset_f,
			 *setup_solution_f;


};

#endif
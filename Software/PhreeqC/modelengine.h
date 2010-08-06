#ifndef modelengineH
#define modelengineH

#include "modelarray.h"
#include "modelunknownmanager.h"

#include "datastructures.h"
#include "parseengine.h"
#include "tokenlist.h"
#include "specieslist.h"
#include "units.h"
#include "array.h"
#include "debug.h"
#include "modeldata.h"
#include "errors.h"
#include "ssassemblagemodel.h"

class ModelEngine:
	public ModelArray,
	public ModelUnknownManager,
	public SSAssemblageModel,
	public Units
{
public:
	ModelEngine(GlobalData *gd, ModelData *md);
	~ModelEngine();

public:
	void SavePPResults(PPAssemblage *ppa);

	void ResetData();
	void FreeModelAllocs();

	void SaveSolution(Solution *dest); //xsolutionsave
	void StepSaveSurf(Surface *dest);
	void StepSaveExch(Exchange *dest);

	void SaveExchange(Exchange *dest);
	void SaveGasPhase(GasPhase *dest);
	void SaveSurface(Surface *dest);

	void CopyUse(int i);
	bool SetInitialMoles(int i);

	bool SetReaction(int i);

	void SaveSolutionResults();

	void PrintResults(String filename);
	void PrintTotals(FILE *f);
	void PrintSpecies(FILE *f);

	bool Prepare();

	bool ConvertUnits(Solution *sol_p); 
	bool Clear(); 
	
	bool SetupSolution();
	bool SetupExchange();
	bool SetupSurface();
	bool SetupPurePhases();
	bool SetupGasPhases();
	bool SetupSSAssemblage();
	bool SetupRelatedSurface();

	bool TidyRedox(); 

	Unknown *FindSurfaceChargeUnknown(String *str, int plane);

	LDBLE SSRoot(LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq, LDBLE xbaq);
	LDBLE SSF(LDBLE xb, LDBLE a0, LDBLE a1, LDBLE kc, LDBLE kb, LDBLE xcaq, LDBLE xbaq);
	LDBLE SSHalve(LDBLE a0, LDBLE a1, LDBLE x0, LDBLE x1, LDBLE kc, LDBLE kb, LDBLE xcaq, LDBLE xbaq);

	LDBLE CalcPSIAvg(LDBLE surf_chrg_eq);
	void SSBinary (SS *s_s_ptr);
	void SSIdeal(SS *s_s_ptr);

	bool CalcInitG();
	bool CalcAllG();
	bool CalcAllDonnanMusic();
	bool CalcInitDonnan();
	bool CalcInitDonnanMusic();

	CONVERGE_RESULT ExecuteModel();
	CONVERGE_RESULT ExecuteSurfaceModel();

	LDBLE QRombMidPnt(LDBLE x1, LDBLE x2);
	LDBLE GFunction (LDBLE x_value);

	LDBLE MidPnt(LDBLE x1, LDBLE x2, int n);
	void Polint(LDBLE * xa, LDBLE * ya, int n, LDBLE xv, LDBLE * yv, LDBLE * dy);

	bool GasPhaseCheck(GasPhase *gas_phase_ptr);
	bool PPAssemblageCheck (PPAssemblage *pp_assemblage_ptr);
	bool CheckPPAssemblage(PPAssemblage *ppa);
	bool SSAssemblageCheck(SSAssemblage *s_s_assemblage_ptr);
	CONVERGE_RESULT SolutionCheck();

	bool AddSolution(Solution *sol_p, LDBLE extensive, LDBLE intensive);
	bool AddPPAssemblage(PPAssemblage *pp_assemblage_ptr);
	bool AddSSAssemblage(SSAssemblage *s_s_assemblage_ptr);
	//bool AddReaction (Irrev *irrev_ptr, int step_number, LDBLE step_fraction);
	bool AddExchange (Exchange *exchange_ptr);
	bool AddSurface(Surface *surface_ptr);
	bool AddGasPhase(GasPhase *gas_phase_ptr);

	Master *MasterPrimarySearch(String str);

	//bool ReactionCalc(Irrev *irrev_ptr);

	bool SumDiffuseLayer(SurfaceCharge *surface_charge_ptr1);

	bool BuildModel(); //build_model
	void PrintCentered(FILE *f, String str);
	bool GetMasterList(String &list, Master *m, ListOfPointers<Master> *copy); //get_list_master_ptrs
	bool SetupMasterReaction(ListOfPointers<Master> *m, Reaction **r); //setup_master_rxn
	bool RewriteMasterToSecondary(Master *m1, Master *m2); //rewrite_master_to_secondary
	bool GetReactionCoef(Reaction *r, String name, LDBLE &coef, int start = 1); //rxn_find_coef & trxn_find_coef
	bool CheckIn(); //inout
	bool WriteMassActionEqnX(); //write_mass_action_eqn_x
	bool AddPotentialFactor(); // add_potential_factor
	bool AddCDMusicFactors(int n); //add_cd_music_factors
	bool WriteMBEqnX(); //write_mb_eqn_x
	bool AddSurfaceChargeBalance(); //add_surface_charge_balance
	bool AddCDMusicChargeBalances(int n); //add_cd_music_charge_balance
	bool MBForSpeciesAQ(int n);
	bool StoreMBUnknowns(Unknown *unknown_ptr, LDBLE *LDBLE_ptr, LDBLE coef, LDBLE *gamma_ptr);
	bool MBForSpeciesEX(int n);
	bool MBForSpeciesSURF(int n);
	bool BuildMBSums();
	bool BuildJacobianSums(int k);
	bool WriteMBForSpeciesList(int n);
	bool BuildSpeciesList(int n);
	bool StoreMB(LDBLE * source, LDBLE * target, LDBLE coef);
	bool WritePhaseSysTotal(int n);
	bool BuildSolutionPhaseBoundaries();
	bool BuildPurePhases();
	bool BuildMinExch();
	bool BuildMinSurface();
	bool BuildGasPhase();
	bool BuildSSAssemblage();
	bool SaveModel();
	bool StoreDN(int k, LDBLE * source, int row, LDBLE coef_in, LDBLE * gamma_source);
	bool StoreJacob(LDBLE * source, LDBLE * target, LDBLE coef);
	bool IsSpecial(Species *spec);
	bool StoreJacob0(int row, int column, LDBLE coef);
	bool ChangeHydrogenInEOSList(LDBLE charge);
	bool StoreSumDeltas(LDBLE * source, LDBLE * target, LDBLE coef);
	Unknown *FindSurfaceChargeUnknown(String &name, int plane); //find_surface_charge_unknown
	
	bool MBGases();
	bool MBSS();
	CONVERGE_RESULT Residuals();
	bool KTemp(LDBLE tempc);
	bool Set(bool initial);
	bool CheckResiduals();
	bool SumSpecies();
	bool NumericalJacobian();
	bool JacobianSums();
	int Ineq(int in_kode);
	bool Reset();
	bool Gammas(LDBLE mu); 
	bool Molalities(bool allow_overflow);
	bool ReviseGuesses();
	bool InitialSurfaceWater();
	bool MBSums();
	bool SwitchBases();
	bool Reprep();
	Master *SurfaceGetPSIMaster(String name, int plane);
	LDBLE Under(LDBLE xval);
	bool SSPrep(LDBLE t, SS *s_s_ptr);
	bool InitialGuesses();
	bool CalculateValues();
	bool IneqInit(int max_row_count, int max_column_count);
	bool CL1(int k, int l, int m, int n,
					 int nklmd, int n2d,
					 LDBLE * q,
					 int *kode, LDBLE toler,
					 int *iter, LDBLE * x, LDBLE * res, LDBLE * error,
					 LDBLE * cu, int *iu, int *s, int check);
	bool CL1Space(int check, int n2d, int klm, int nklmd);
	bool CalcAllDonnan();
	bool CalcGasPressures();
	bool CalcSSFractions();
	bool ResetupMaster();
	LDBLE CalcSC();
	
	bool RunReactions(int i, LDBLE step_fraction);

	CONVERGE_RESULT SetAndRunWrapper(int i, LDBLE step_fraction, int nsaver);
	CONVERGE_RESULT SetAndRun(int i, int nsaver, LDBLE step_fraction);
	CONVERGE_RESULT Step(LDBLE step_fraction);
	
	Debug d;

	ErrorFile *error;

private:
	GlobalData *gd;
	ModelData *md;	
	Reaction *r_temp;

	PPAssemblage ppa_save_1, ppa_save_2;
	SSAssemblage ssa_save_1, ssa_save_2;

	TokenList tl; //Used in GetMasterList
	TokenList token; //Used in TidyRedox



	char error_string[3000];
};

#endif
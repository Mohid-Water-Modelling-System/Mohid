#ifndef datastructuresH
#define datastructuresH

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>

#include "exceptionhandler.h"
#include "definitions.h"
#include "units.h"
#include "basicfunctions.h"
#include "constants.h"
#include "string.h"
#include "list.h"
#include "listofpointers.h"
#include "data.h"

using namespace std;
//============================================================================================================
// class ReactionToken
//============================================================================================================
class ReactionToken: //rxn_token_temp & rxn_token
	public Data
{
public:
	ReactionToken();
	ReactionToken(ReactionToken *copy);
	ReactionToken(class Species *s, const LDBLE coef);
	~ReactionToken();

public:
	void Reset();
	//ReactionToken *Copy();
	void CopyTo(ReactionToken *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	String name;

  LDBLE coef;
  LDBLE z; //Charge on species

	class Unknown *u;
  class Species *s;
};

//============================================================================================================
// class Reaction
//============================================================================================================
class Reaction: //reaction
	public BasicFunctions,
	public Data
{
public:
	Reaction();
	Reaction(Reaction *copy);
	~Reaction();

public:
	void Reset(); 
	Reaction *Copy();
	void CopyTo(class Reaction *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	void AddTRXN(class Reaction *r, LDBLE coef, bool combine);	
	void AddPhaseTRXN(class Reaction *r, LDBLE coef, bool combine);
	void AddPhaseTRXN(class Phase *p, class Reaction *r); //CopyPhaseReactionToTempReaction
	void AddSpeciesTRXN(Species *s);
	void TRXNCombine();
	void TRXNReverseK();
	bool TRXNSwap(String token, int position = 0);
	bool TRXNMultiply(LDBLE coef);

public:
	void CopyReactions(class Reaction *reaction);
	void CopyReactionsToSpecies(class Reaction *reaction);
	void CopyReactionsToPhase(class Reaction *reaction);
	void CopyReactionsToExchange(class Reaction *reaction);

public:
	String name;
  LDBLE logk[8];
  List<ReactionToken> *token_list;
};

//============================================================================================================
// class Element
//============================================================================================================
class Element: //element
	public Data
{
public:
	Element();
	Element(Element *copy);
	~Element();

public:
	void Clear();
	void Reset();
	//Element *Copy();
	void CopyTo(Element *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	String name;
  class Master *master;
  class Master *primary;
  LDBLE gfw;
};

//============================================================================================================
// class Master
//============================================================================================================
class Master: //master
	public Data
{
public:
	Master();
	Master(Master *copy);
	~Master();

public:
	void Clear(); //Clear only modeled data
	void Reset(); //Reset ALL data
	//Master *Copy();
	void CopyTo(Master *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	String name;
	String gfw_formula;	// formula from which to calcuate gfw

	bool in;	// true if in model, FALSE if out, REWRITE if other mb eq
	bool rewrite; 
	bool primary; // true if master species is primary

	int number;	// sequence number in list of masters
	int type;	// AQ or EX
	
	LDBLE coef;	// coefficient of element in master species
	LDBLE total; // total concentration for element or valence state
	LDBLE alk; // alkalinity of species 
	LDBLE gfw; // default gfw for species
	LDBLE total_primary;

	class Element *e;	// element structure
	class Species *s;	// pointer to species structure
	class Unknown *u;	// pointer to unknown structure

	class Reaction *rxn_primary;	// reaction writes master species in terms of primary master species
	class Reaction *rxn_secondary;	// reaction writes master species in terms of secondary master species
	class Reaction **pe_rxn;	// e- written in terms of redox couple (or e-), points to location
};

//============================================================================================================
// class ElementOfSpecies
//============================================================================================================
class ElementOfSpecies //eos_list
{
public:
	ElementOfSpecies();
	ElementOfSpecies(ElementOfSpecies *copy);
	~ElementOfSpecies();

public:
	void Reset();
	//ElementOfSpecies *Copy();
	void CopyTo(ElementOfSpecies *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	String name;
  class Element *e; // pointer to element structure
  LDBLE coef; // number of element e's in eqn
};

//============================================================================================================
// class NameCoef
//============================================================================================================
class NameCoef //name_coef
{
public:
	NameCoef();
	NameCoef(NameCoef *copy);
	~NameCoef();

public:
	void Reset();
	//NameCoef *Copy();
	void CopyTo(NameCoef *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
  String name;
  LDBLE coef;
};

//============================================================================================================
// class Species
//============================================================================================================
class Species //species
{
public:
	Species();
	Species(Species *copy);
	~Species();

public:
	void Clear();
	void Reset();
	//Species *Copy();
	void CopyTo(Species *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	String name;
  String mole_balance; // formula for mole balance

  LDBLE gfw; // gram formula wt of species
  LDBLE z; // charge of species 
  LDBLE dw; // tracer diffusion coefficient in water at 25oC, m2/s 
  LDBLE equiv; // equivalents in exchange species 
  LDBLE alk; // alkalinity of species, used for cec in exchange 
  LDBLE carbon; // stoichiometric coefficient of carbon in species 
  LDBLE co2; // stoichiometric coefficient of C(4) in species 
  LDBLE h; // stoichiometric coefficient of H in species 
  LDBLE o; // stoichiometric coefficient of O in species 
  LDBLE dha, dhb; // WATEQ Debye Huckel a and b-dot 
  LDBLE lk; // log10 k at working temperature 
  LDBLE logk[8]; // log kt0, delh, 6 coefficients analalytical expression 
  LDBLE lg; // log10 activity coefficient, gamma 
  LDBLE lg_pitzer; // log10 activity coefficient, from pitzer calculation 
  LDBLE lm; // log10 molality 
  LDBLE la;	// log10 activity 
  LDBLE dg;	// gamma term for jacobian 
  LDBLE dg_total_g;
  LDBLE moles; // moles in solution; moles/mass_water = molality 
  LDBLE tot_g_moles; // (1 + sum(g)) * moles 
  LDBLE tot_dh2o_moles; // sum(moles*g*Ws/Waq) 

	DELTA_H_UNIT original_units;

  int number;
  int type;	// flag indicating presence in model and types of equations 
  int gflag; // flag for preferred activity coef eqn 
  int exch_gflag; // flag for preferred activity coef eqn 

  bool in; // species used in model if true 
  bool check_equation; // switch to check equation for charge and element balance 

  LDBLE cd_music[5];
  LDBLE dz[3];

  class Master *primary;	// points to master species list, NULL if not primary master 
  class Master *secondary;	// points to master species list, NULL if not secondary master 

	class Reaction *rxn;		// pointer to data base reaction
  class Reaction *rxn_s;	// pointer to reaction converted to secondary and primary master species 
	class Reaction *rxn_x;	// reaction to be used in model 

	List<class ElementOfSpecies> *eos_list; // pointer to next element
  List<class ElementOfSpecies> *e_sec_list; //*next_secondary;	
	List<class ElementOfSpecies> *eos_sys_total; //*next_sys_total;
  List<class NameCoef> *add_logk; 
	List<class SpeciesDiffLayer> *diff_layer;
};

//============================================================================================================
// class SpeciesDiffLayer
//============================================================================================================
class SpeciesDiffLayer
{
public:
	SpeciesDiffLayer();
	SpeciesDiffLayer(SpeciesDiffLayer *copy);
	~SpeciesDiffLayer();

public:
	void Reset();
	//SpeciesDiffLayer *Copy();
	void CopyTo(SpeciesDiffLayer *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	LDBLE g_moles,
		    dg_g_moles,
				dx_moles,
				drelated_moles,
				dh2o_moles;

	class SurfaceCharge *charge; //ToDo: Doubt if this is to allocate or only to "point"...

	int count_g;
};

//============================================================================================================
// class Phase
//============================================================================================================
class Phase //phase
{				
public:
	Phase();
	Phase(Phase *copy);
	~Phase();

public:
	void Clear();
	void Reset();
	//Phase *Copy();
	void CopyTo(Phase *copy);
	void PrintToFile(FILE *file, int spaces);

public:
	String name, // name of species
				 formula;		// chemical formula 

	bool check_equation, // switch to check equation for charge and element balance
			 in, // species used in model if true
			 in_system;

	DELTA_H_UNIT original_units;	// enum with original delta H units

	List<class NameCoef> *add_logk; //class NameCoef
  List<class ElementOfSpecies> *eos_list; // pointer to list of elements in phase //class ElementOfSpecies *next_elt;	
	List<class ElementOfSpecies> *eos_sys_total; //next_sys_total;

	class Reaction *rxn; // pointer to data base reaction
  class Reaction *rxn_s; // pointer to reaction converted to secondary and primary master species
	class Reaction *rxn_x; // reaction to be used in model

  int type;	// flag indicating presence in model and types of equations

  LDBLE logk[8], // log kt0, delh, 6 coefficients analalytical expression
	      lk,	// log10 k at working temperature
        moles_x,
				p_soln_x,
				fraction_x,
				log10_lambda, 
				log10_fraction_x,
				dn, 
				dnb, 
				dnc,
				gn, 
				gntot,
				gn_n, 
				gntot_n;
};

//============================================================================================================
// class Unknown
//============================================================================================================
class Unknown: //unknown
	public Data
{
public:
	Unknown();
	Unknown(Unknown *copy);
	~Unknown();

public:
	void Reset();
	//Unknown *Copy();
	void CopyTo(Unknown *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
	//stores the results to be passed to MOHID (result defined in SaveSolutionResults function)		  
	LDBLE molality_result, moles_result; 

	int type;

	LDBLE moles,
				ln_moles,
	      f,
				sum,
				delta,
				la,
				si,
				related_moles,
				mass_water;

	bool dissolve_only,
			 s_s_in;

	int number,
			s_s_comp_number;

	String name;

	ListOfPointers<class Master> *master;
	ListOfPointers<class Unknown> *comp_unknowns;
	
	class Conc *total;
	class Phase *p;
	class SurfaceComp *surface_comp;
	class PurePhase *pure_phase;
	class GasPhase *gas_phase;
	class Species *s;
	class SurfaceCharge *surface_charge;
	class Unknown *phase_unknown;
	class SSComp *s_s_comp;
	class ExchComp *exch_comp;
	class Unknown *potential_unknown, 
								*potential_unknown1, 
								*potential_unknown2;
	class SS *s_s;
};

//============================================================================================================
// class SpeciesInfo
//============================================================================================================
class SpeciesInfo
{	
public:
	SpeciesInfo();
	SpeciesInfo(SpeciesInfo *copy);
	~SpeciesInfo();

public:
	void Reset();
	//SpeciesInfo *Copy();
	void CopyTo(SpeciesInfo *copy);

public:
	class Species *master_s;
  class Species *s;

  LDBLE coef;
};

//============================================================================================================
// class Surface
//============================================================================================================
struct SurCompData
{
	String formula,
				 phase_name,
				 rate_name;

	LDBLE moles,
				phase_proportion;

	LDBLE area, grams;

	bool use_moles, use_area, use_grams;
};

class Surface
{
public:
	Surface();
	Surface(Surface *copy);
	~Surface();

public:
	void Reset();
	void CopyTo(Surface *copy);
	//Surface *Copy();

public:
	void SetEquilibria(bool on_off = true) { solution_equilibria = on_off; }
	void SetDiffuseLayer(LDBLE thickness = 1e-8) { dl_type = BORKOVEK_DL; this->thickness = thickness; }
	void SetNoElectrostatic() { type = NO_EDL; }
	void SetOnlyCounterIons(bool on_off = true) { only_counter_ions = on_off; }
	void SetSurComp(SurCompData data);

public:
  SURFACE_TYPE type;

	bool new_def;

	int number;

	DIFFUSE_LAYER_TYPE dl_type;

	List<class SurfaceComp> *comps;
	List<class SurfaceCharge> *charge;
  
	bool related_phases,
			 related_rate,
			 only_counter_ions,
			 solution_equilibria;

	LDBLE debye_lengths,
		    thickness,
				DDL_limit,
				DDL_viscosity;

	SITES_UNITS sites_units;


};

//============================================================================================================
// class SurfaceComp
//============================================================================================================
class SurfaceComp
{
public:
	SurfaceComp();
	SurfaceComp(SurfaceComp *copy);
	~SurfaceComp();

public:
	void Reset();
	void CopyTo(SurfaceComp *copy);
	//SurfaceComp *Copy();

public:
	String name; //same as master->name (that is same as master->e->name)
  String formula;
  String phase_name;
  String rate_name;

	List<class ElementOfSpecies> *formula_totals;
	List<class ElementOfSpecies> *totals;

  class Master *master;

  int charge;

	LDBLE formula_z;
  LDBLE moles;
  LDBLE la;
  LDBLE cb;
  LDBLE phase_proportion;
  LDBLE Dw;			// diffusion coefficient in water, used in MCD. No transport if 0

	bool use_moles, use_area, use_grams, related_phases;
	LDBLE area, grams;
};

//============================================================================================================
// class SurfaceCharge
//============================================================================================================
class SurfaceCharge
{
public:
	SurfaceCharge();
	SurfaceCharge(SurfaceCharge *copy);
	~SurfaceCharge();

public:
	void Reset();
	//SurfaceCharge *Copy();
	void CopyTo(SurfaceCharge *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
  String name;

  LDBLE specific_area;
  LDBLE grams;
  LDBLE charge_balance;
  LDBLE mass_water;
  
	List<class ElementOfSpecies> *diffuse_layer_totals;
  List<class SurfaceDiffLayer> *g;	// stores g and dg/dXd for each ionic charge
  //int count_g;

	LDBLE la_psi, la_psi1, la_psi2;
  LDBLE psi, psi1, psi2;

  LDBLE capacitance[2];
  LDBLE sigma0, sigma1, sigma2, sigmaddl;
};

//============================================================================================================
// class SurfaceDiffLayer
//============================================================================================================
class SurfaceDiffLayer
{
public:
	SurfaceDiffLayer();
	SurfaceDiffLayer(SurfaceDiffLayer *copy);
	~SurfaceDiffLayer();

public:
	void Reset();
	//SurfaceDiffLayer *Copy();
	void CopyTo(SurfaceDiffLayer *copy);

public:
  LDBLE charge;
  LDBLE g;
  LDBLE dg;
  LDBLE psi_to_z;
};

//============================================================================================================
// class UnknownInfo
//============================================================================================================
class UnknownInfo
{
public:
	UnknownInfo();
	UnknownInfo(UnknownInfo *copy);
	~UnknownInfo();

public:
	void Reset();
	void CopyTo(UnknownInfo *copy);
	void PrintToFile(FILE *file, int spaces = 0);

public:
  Unknown *u;

	LDBLE *source;
  LDBLE *gamma_source;
  LDBLE coef;
};

//============================================================================================================
// class STCoef
//============================================================================================================
class STCoef //list_1 & list_2
{
public:
	STCoef();
	STCoef(STCoef *copy);
	~STCoef();
	
public:
	void Reset();
	void CopyTo(STCoef *copy);

public:
	void PrintToFile(FILE *file, int spaces = 0);

public:
  LDBLE *source;
  LDBLE *target;
	LDBLE coef;
};

//============================================================================================================
// class PPAssemblage
//============================================================================================================
struct PPData
{
  String name;
  String add_formula;

  LDBLE si;	
	LDBLE moles;

  bool force_equality;
  bool dissolve_only;
};

class PPAssemblage
{
public:
	PPAssemblage(class GlobalData *gd = NULL);
	PPAssemblage(PPAssemblage *copy);
	~PPAssemblage();

public:
	void Reset();
	void CopyTo(PPAssemblage *copy);

public:
	int AddPurePhase(PPData &data, int &PurePhaseID);
	void ChangeMoles(String *name, LDBLE new_moles);
	void Print(FILE *file);

public:
	int number;

	bool new_def;

  List<class ElementOfSpecies> *eos_list;
  List<class PurePhase> *pure_phases;

public:
	class GlobalData *gd;
};

//============================================================================================================
// class PurePhase
//============================================================================================================
class PurePhase
{
public:
	PurePhase();
	PurePhase(PurePhase *copy);
	~PurePhase();

public:
	void Reset();
	void Clear();
	//PurePhase *Copy();
	void CopyTo(PurePhase *copy);

public:
  class Phase *phase;

  String name;
  String add_formula;

  LDBLE si;
  LDBLE moles;
  LDBLE delta;
  LDBLE initial_moles;

  bool force_equality;
  bool dissolve_only;
};

//============================================================================================================
// class Exchange
//============================================================================================================
class Exchange
{
public:
	Exchange();
	Exchange(Exchange *copy);
	~Exchange();

public:
	void Reset();
	void CopyTo(Exchange *copy);
	//Exchange *Copy();
	void Print(FILE *file);

public:
	int AddExchanger(class ExchCompData *exchanger);

public:
	String name;
	bool new_def;
	int number;

	List<class ExchComp> *comps;

	bool solution_equilibria;

	bool related_phases,
			 related_rate;
};

//============================================================================================================
// class ExchCompData
//============================================================================================================
class ExchCompData
{
public:
	ExchCompData()
	{
		formula = "";
		name = "";
		amount = 0.0;
		type = 0;
	}

public:
	String formula, //exchange formula
				 name; //name of the phase or kinect reaction
	LDBLE amount; //can be exchange concentration if type 0 or exchange_per_mole if type 1 or 2
	int type; //0 - none, 1 - equilibrium_phase, 2 - kinetic_reactant

};

//============================================================================================================
// class ExchComp
//============================================================================================================
class ExchComp
{
public:
	ExchComp();
	ExchComp(ExchComp *copy);
	~ExchComp();

public:
	void Reset();
	void CopyTo(ExchComp *copy);
	//ExchComp *Copy();

public:
	String name;
	int type;
	LDBLE amount;

  String formula;
  String phase_name;
  String rate_name;
  
	LDBLE formula_z;
  LDBLE moles;
  LDBLE la;
  LDBLE charge_balance;
  LDBLE phase_proportion;

	List<ElementOfSpecies> *formula_totals;
	List<ElementOfSpecies> *totals;
  
	class Master *master;  
};


//============================================================================================================
// class GasPhase
//============================================================================================================
typedef enum { PRESSURE = 1, VOLUME } GasPhaseType;

class GasPhase
{
public:
	GasPhase();
	GasPhase(class GasPhase * copy);
	~GasPhase();

public:
	void Clear();
	void Reset();
	//GasPhase *Copy();
	void CopyTo(GasPhase *copy);
	void PrintToFile(FILE *file, int spaces);

public:
	void SetPressure(LDBLE value) { total_p = value; }
	void SetVolume(LDBLE value) { volume = value; }
	void SetTemperature(LDBLE value) { temperature = value + (LDBLE)273.15; }
	void SetType(GasPhaseType type) { this->type = type; }
	void SetEquilibrium(bool on_off) { solution_equilibria = on_off; }
	void AddGasComp(String name, LDBLE p_read = NAN);

public:
  //String description;

	bool solution_equilibria,
			 new_def;

	int number;

  GasPhaseType type;

  LDBLE total_p; // pressure
  LDBLE total_moles;
  LDBLE volume;
  LDBLE temperature;

  List<class GasComp> *comps;
};

//============================================================================================================
// class SSAssemblage
//============================================================================================================
class SSAssemblage
{
public:
	SSAssemblage();
	SSAssemblage(SSAssemblage *copy);
	~SSAssemblage();

public:
	void Reset();
	void CopyTo(SSAssemblage *copy);
	//SSAssemblage *Copy();

public:
  List<class SS> *ss_list;

	int number;
};

//============================================================================================================
// class SS
//============================================================================================================
struct SSCompData
{
public:
	String name; //phase name
	LDBLE moles; //NAN if not inputed

public:
	SSCompData() { moles = (LDBLE)NAN; }
};

class SS
{
public:
	SS();
	SS(SS *copy);
	~SS();

public:
	void Reset();
	void CopyTo(SS* copy);
	//SS* Copy();

public:
	bool AddComp(struct SSCompData *comp); //Add comps to solid solution
	void NondimensionalGuggenheim(LDBLE p0 = 0.0, LDBLE p1 = 0.0);
	
public:
  String name;

  List<class SSComp> *comps_list;

  LDBLE total_moles,
				dn,
				a0, 
				a1,
				ag0, 
				ag1,
				tk, 
				xb1, 
				xb2,
				p[4];

  bool s_s_in,
			 ideal,
	     miscibility,
			 spinodal;

	int input_case;
};

//============================================================================================================
// class SSComp
//============================================================================================================
class SSComp
{
public:
	SSComp();
	SSComp(SSComp *copy);
	~SSComp();

public:
	void Reset();
	void CopyTo(SSComp *copy);
	//SSComp *Copy();

public:
  String name;

  class Phase *phase;

  LDBLE initial_moles,
				moles,
				init_moles,
				delta,
				fraction_x,
				log10_lambda,
				log10_fraction_x,
				dn, 
				dnc, 
				dnb;
};

//============================================================================================================
// class GasComp
//============================================================================================================
class GasComp
{
public:
	GasComp();
	GasComp(GasComp *copy);
	~GasComp();

public:
	void Clear();
	void Reset();
	//GasComp *Copy();
	void CopyTo(GasComp *copy);
	void PrintToFile(FILE *file, int spaces);

public:
  class Phase *phase;

  String name;

  LDBLE p_read;
  LDBLE moles;
  LDBLE initial_moles;
};

//============================================================================================================
// class Irrev
//============================================================================================================
/*
class Irrev
{
public:
  List<NameCoef> *list;
	List<ElementOfSpecies> *elts;
  LDBLE *steps;
  char *units;
  int count_steps;
  int count_list;
};
*/

//============================================================================================================
// class GlobalData
//============================================================================================================
class GlobalData
{
public:
	GlobalData();
	~GlobalData();

public:
	void ResetToExecuteModel();

public:
	List<Element> element_list; //Element     
	List<Species> species_list; //Species
	List<Master> master_list; //Master
	List<Phase> phase_list; //Phase

	List<class ChargeGroup> charge_group; 

	Species *s_hplus,
					*s_h3oplus,
					*s_eminus,
					*s_h2o,
					*s_o2,
					*s_co3,
					*s_h2;

	Element *e_h_one; //element_h_one

	PPAssemblage *ppa;
};

//============================================================================================================
// class Solution
//============================================================================================================
struct PXData
{
	LDBLE value,
			  sat_index;

	String phase;
	
	bool charge;
};

struct RedoxData
{
	String element_1,
				 element_2;
	
	LDBLE valence_1,
				valence_2;
};

struct ConcData
{
public:
	ConcData() 
	{
		Reset();
	}

	void Reset()
	{
		element = "";
		as = "";
		phase = "";
		redox.element_1 = "";
		redox.element_2 = "";
		redox.valence_1 = 0.0;
		redox.valence_2 = 0.0;
		for (int i = 0; i < 5; i ++) use[i] = false;
		gfw = 0.0;
		value = 0.0;
		unit = _mol_kgw_;
		id = -1;
		sat_index = 0.0;
		charge = false;
	}

public:
	String element,
				 as,
				 phase;

	UNITS unit;

	LDBLE gfw,
		    value,
				sat_index;

	bool charge;

	struct RedoxData redox;

	int id;

	bool use[5]; //0 - as, 1 - unit, 2 - gfw, 3 - redox, 4 - charge or phase
};

class Solution
{
public:
	Solution(); //used for COPY only, because do not initialize PE
	Solution(Solution *copy);
	Solution(GlobalData *gd);
	~Solution();

public:
	void Reset();
	void CopyTo(Solution *dest);
	void InitPE() { InitPE(gd); }
	void InitPE(GlobalData *gd);

public:
	void SetPE(PXData &pe);
	void SetPH(PXData &ph);
	int SetRedox(RedoxData &rd);
	int SetElementConc(ConcData &cd, int &MasterSpeciesID, LDBLE &gfw);
	void ChangeElementConc(char *element, LDBLE new_conc);
	void ChangeElementConc(int index, LDBLE new_conc);

public:
	void Print(FILE *file);

public:
	int number;

  LDBLE tc;
  LDBLE ph;
  LDBLE solution_pe;
  LDBLE mu;
  LDBLE ah2o;
  LDBLE density;
  LDBLE total_h;
  LDBLE total_o;
  LDBLE cb;
  LDBLE mass_water;
  LDBLE total_alk; //total_alkalinity

  UNITS units;

  int default_pe;

	List<class PEData> *pe; //struct pe_data *e;
	List<class Conc> *totals; //struct conc *totals;
	List<class MasterActivity> *ma,
		                         *species_gamma;

	GlobalData *gd;

	Conc *ph_conc_ptr,
		   *pe_conc_ptr;

private:
	Conc * SetConc(CONC_TYPE type, LDBLE value, String phase, LDBLE sat_index, bool charge);

};

//============================================================================================================
// class MasterActivity
//============================================================================================================
class MasterActivity
{
public:
	MasterActivity();
	MasterActivity(MasterActivity *copy);
	~MasterActivity();

public:
	void Reset();
	//MasterActivity *Copy();
	void CopyTo(MasterActivity * copy);

public:
	String description;
	LDBLE la;
};

//============================================================================================================
// class ChargeGroup
//============================================================================================================
class ChargeGroup
{
public:
	ChargeGroup();
	ChargeGroup(ChargeGroup *copy);
	~ChargeGroup();

public:
	void Reset();
	void CopyTo(ChargeGroup *copy);
	//ChargeGroup *Copy();

public:
  LDBLE z;
  LDBLE eq;
};

//============================================================================================================
// class PEData
//============================================================================================================
class PEData
{
public:
	PEData();
	PEData(PEData *copy);
	~PEData();

public:
	bool operator == (const PEData &right);

public:
	void Reset();
	//PEData *Copy();
	void CopyTo(PEData *copy);

public:
	String name;

	String element_1,
				 element_2;

	LDBLE valence_1,
				valence_2;

	class Reaction *rxn;

private:
	void AllocMemory();

};

//============================================================================================================
// class Conc
//============================================================================================================
class Conc
{
public:
	Conc();
	Conc(Conc *copy);
	~Conc();

public:
	void Reset();
	//Conc *Copy();
	void CopyTo(Conc *copy);

public:
	bool SetDescription(String description);

public:
  String name; //description;

	CONC_TYPE type;
	bool charge;

  LDBLE moles;
  LDBLE input_conc;
	LDBLE phase_si;
  LDBLE gfw;

  UNITS units;

	String equation_name;
  String as;
  
	Phase *p; //phase
  
  int n_pe;

	int id;

	ListOfPointers<Master> *m_list;

	//class Unknown *x;
};

//============================================================================================================
// class Isotope
//============================================================================================================
/*
class Isotope
{
public:
	Isotope();
	Isotope(Isotope *copy);
	~Isotope();

public:
	void Reset();
	Isotope *Copy();
	void CopyTo(Isotope *copy);

public:
  LDBLE isotope_number;
  LDBLE total;
  LDBLE ratio;
  LDBLE ratio_uncertainty;
  LDBLE x_ratio_uncertainty;
  LDBLE coef;			// coefficient of element in phase 

	String elt_name;
  String name; //char *isotope_name;

	//Master *master;
  //Master *primary;
};
*/

#endif
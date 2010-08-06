#ifndef constantsH
#define constantsH

#include "definitions.h"

#ifndef NAN
#   define NAN -99999999
#endif

#define MISSING (LDBLE)-9999.999

#define F_C_MOL (LDBLE)96493.5		/* C/mol or joule/volt-eq */
#define F_KJ_V_EQ  (LDBLE)96.4935	/* kJ/volt-eq */
#define F_KCAL_V_EQ (LDBLE)23.0623	/* kcal/volt-eq */
#define R_LITER_ATM (LDBLE)0.0820597	/* L-atm/deg-mol */
#define R_KCAL_DEG_MOL (LDBLE)0.00198726	/* kcal/deg-mol */
#define R_KJ_DEG_MOL (LDBLE)0.00831470	/* kJ/deg-mol */
#define EPSILON (LDBLE)78.5		/* dialectric constant, dimensionless */
#define EPSILON_ZERO (LDBLE)8.854e-12	/* permittivity of free space, C/V-m = C**2/m-J */
#define JOULES_PER_CALORIE (LDBLE)4.1840
#define AVOGADRO (LDBLE)6.02252e23	/* atoms / mole */

#define DEFAULT_STR_LENGTH 1024
#define TOLERANCE (LDBLE)1e-9

#define TRUE 1
#define FALSE 0
#define OK 1
#define ERROR 0
#define STOP 1
#define CONTINUE 0

#define DISP 2
#define STAG 3
#define NOMIX 4

#define CONVERGED 2
#define MASS_BALANCE 3
/*
  #define OSCILLATE 4
  #define H2O_LIMITS 5
*/
#define REWRITE 2
#define INIT -1

/* check_line values, plus EMPTY, EOF, OK */
#define KEYWORD 3

/* copy_token values */
#define EMPTY 2
#define UPPER 4
#define LOWER 5
#define DIGIT 6
#define UNKNOWN 7
#define OPTION 8

/* species types */
#define AQ 0
#define HPLUS 1
#define H2O 2
#define EMINUS 3
#define SOLID 4
#define EX 5
#define SURF 6
#define SURF_PSI 7
#define SURF_PSI1 8
#define SURF_PSI2 9

/* unknown types */
#define MB 10
#define ALK 11
#define CB 12
#define SOLUTION_PHASE_BOUNDARY 13
#define MU 14
#define AH2O 15
#define MH 16
#define MH2O 17
#define PP 18
#define EXCH 19
#define SURFACE 20
#define SURFACE_CB 21
#define SURFACE_CB1 22
#define SURFACE_CB2 23
#define GAS_MOLES 24
#define S_S_MOLES 25
#define PITZER_GAMMA 26
/* state */
#define INITIALIZE               0
#define INITIAL_SOLUTION   1
#define INITIAL_EXCHANGE   2
#define INITIAL_SURFACE 3
#define INITIAL_GAS_PHASE  4
#define REACTION                   5
#define INVERSE                 6
#define ADVECTION                 7
#define TRANSPORT                 8
#define PHAST                     9

/* constaints in mass balance */
#define EITHER 0
#define DISSOLVE 1
#define PRECIPITATE -1

/* gas phase type */
//#define PRESSURE 1
//#define VOLUME 2

#define MAX_PP_ASSEMBLAGE 10	/* default estimate of the number of phase assemblages */
#define MAX_ADD_EQUATIONS 20	/* maximum number of equations added together to reduce eqn to
				   master species */
#define MAX_ELEMENTS 50		/* default estimate of the number of elements */
#define MAX_LENGTH 256		/* maximum number of characters component name */
#define MAX_LINE 80		/* estimate of maximum line length */
#define MAX_LM (LDBLE)3.0		/* maximum log molality allowed in intermediate iterations */
#define MIN_LM (LDBLE)-30.0		/* minimum log molality allowed before molality set to zero */
#define MAX_MASS_BALANCE 10	/* initial guess of number mass balance equations for a solution */
#define MAX_MASTER 50		/* default estimate of the number of master species */
#define MAX_ELTS 15		/* default estimate for maximum number of times elements occur in
				   an equation */
#define MAX_PHASES 500		/* initial guess of number of phases defined */
#define MAX_SOLUTION 10		/* The maximum number of solutions allowed */
#define MAX_S 500		/* default estimate for maximum number of species in aqueous model */
#define MAX_STRINGS 3000
#define MAX_SUM_JACOB0 50	/* list used to calculate jacobian */
#define MAX_SUM_JACOB1 500	/* list used to calculate jacobian */
#define MAX_SUM_JACOB2 500	/* list used to calculate jacobian */
#define MAX_SUM_MB 500		/* list used to calculate mass balance sums */
#define MAX_TRXN 16		/* default estimate for maximum number of components in an eqn */
#define MAX_UNKNOWNS 15		/* default estimate for maximum number of unknowns in model */
#define TOL (LDBLE)1e-9		/* tolerance for comparisons of double numbers */
#define LOG_ZERO_MOLALITY -30	/* molalities <= LOG_ZERO_MOLALITY are considered equal to zero */
#define MIN_TOTAL 1e-25
#define MIN_TOTAL_SS MIN_TOTAL
#define MIN_RELATED_SURFACE MIN_TOTAL*100
#define MIN_RELATED_LOG_ACTIVITY -30

#define MAX_QUAD 20
#define K_POLY 5


typedef enum 
{ 
	_end_, _keyword_, _empty_, _option_, _none_, _error_, _upper_, _lower_, _digit_, _unknown_ 
} RETURN_TYPE;

/*
typedef enum 
{ 
	_kcal_, _cal_, _kjoules_, _joules_ 
} DELTA_H_UNIT;
*/

#define DELTA_H_UNIT int

#define _kcal_    0
#define _cal_     1
#define _kjoules_ 2
#define _joules_  3

/*
typedef enum 
{ 
	_unknown_unit_,
	_mol_l_, _mmol_l_, _umol_l_, 
	_g_l_, _mg_l_, _ug_l_, 
	_eq_l_, _meq_l_, _ueq_l_, 
	_mol_kgs_, _mmol_kgs_, _umol_kgs_, 
	_g_kgs_, _mg_kgs_, _ug_kgs_, 
	_eq_kgs_, _meq_kgs_, _ueq_kgs_,
	_mol_kgw_, _mmol_kgw_, _umol_kgw_, 
	_g_kgw_, _mg_kgw_, _ug_kgw_,
	_eq_kgw_, _meq_kgw_, _ueq_kgw_
} UNITS;
*/

#define UNITS int
#define _unknown_unit_ 0
#define _mol_l_        1
#define _mmol_l_			 2
#define _umol_l_       3
#define _g_l_          4
#define _mg_l_         5
#define _ug_l_         6
#define _eq_l_         7
#define _meq_l_        8
#define _ueq_l_        9
#define _mol_kgs_      10
#define _mmol_kgs_     11
#define _umol_kgs_     12
#define _g_kgs_        13
#define _mg_kgs_       14
#define _ug_kgs_       15 
#define _eq_kgs_       16
#define _meq_kgs_      17
#define _ueq_kgs_      18
#define _mol_kgw_      19
#define _mmol_kgw_     20
#define _umol_kgw_     21
#define _g_kgw_        22
#define _mg_kgw_       23
#define _ug_kgw_       24
#define _eq_kgw_       25
#define _meq_kgw_      26
#define _ueq_kgw_      27

/*
typedef enum
{
	_unknown_type_, _ph_, _pe_, _conc_
} CONC_TYPE;
*/

#define CONC_TYPE int
#define _unknown_type_ 0
#define _ph_           1
#define _pe_           2
#define _conc_         3

/*
typedef enum 
{ 
	UNKNOWN_DL, NO_EDL, DDL, CD_MUSIC 
} SURFACE_TYPE;
*/

#define SURFACE_TYPE int

#define UNKNOWN_DL 0
#define NO_EDL     1
#define DDL        2
#define CD_MUSIC   3

/*
typedef enum 
{ 
	NO_DL, BORKOVEK_DL, DONNAN_DL 
} DIFFUSE_LAYER_TYPE;
*/

#define DIFFUSE_LAYER_TYPE int

#define NO_DL       0
#define BORKOVEK_DL 1
#define DONNAN_DL   2

/*
typedef enum
{
	CONVERGED_CR, MASS_BALANCE_CR, NOT_CONVERGED_CR, OK_CR, ERROR_CR
} CONVERGE_RESULT;
*/

#define CONVERGE_RESULT int

#define CONVERGED_CR     0
#define MASS_BALANCE_CR  1
#define NOT_CONVERGED_CR 2
#define OK_CR            3
#define ERROR_CR         4

//typedef double realtype;

#define ZERO_TOL (LDBLE)1.0e-30

/*
enum SITES_UNITS
{ SITES_ABSOLUTE, SITES_DENSITY };
*/

#define SITES_UNITS int

#define SITES_ABSOLUTE 0
#define SITES_DENSITY  1

#endif
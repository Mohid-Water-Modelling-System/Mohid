#ifndef definitionsH
#define definitionsH

#define _CRT_SECURE_NO_WARNINGS

#ifdef DEBUG_MOHID
	#define DEBUG_ALL_FUNCTIONS
#endif

#ifdef DEBUG_ALL_FUNCTIONS
	#define DEBUG_AddSolution
	#define DEBUG_AddExchange
	#define DEBUG_AddPPAssemblage
#endif

#ifdef USE_LONG_DOUBLE
	#define LDBLE long double
	#define SCANFORMAT "%Lf"
#else
	#ifdef USE_DOUBLE
		#define LDBLE double
		#define SCANFORMAT "%lf"
	#else		
		#define LDBLE float
		#define SCANFORMAT "%f"
	#endif
#endif

#endif
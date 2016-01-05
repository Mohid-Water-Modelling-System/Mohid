!*MODULE PhreeqcRM PHREEQC Reaction Module for Transport Codes
!> @brief Fortran Documentation for the geochemical reaction module PhreeqcRM. 
!> @par "" 
!> "USE PhreeqcRM" is included in Fortran source code to define the PhreeqcRM functions.
!> For Windows, define the module by including the file RM_interface.F90 in your project.
!> For Linux, configure, compile, and install the PhreeqcRM library and module file. 
!> You will need installed include directory (-I) added to the project) to reference the module file.
!> You will need to link to the library to produce the executable for your code.
!>
    MODULE PhreeqcRM
    IMPLICIT NONE
    !!!SAVE
#if defined(NDEBUG)
    LOGICAL :: rmf_debug=.false.
#else
    LOGICAL :: rmf_debug=.true.
#endif     
    INTEGER, PRIVATE  :: rmf_nxyz=-1
    INTEGER, PRIVATE  :: rmf_ncomps=-1
    PRIVATE :: ChK_Concentrations2Utility
    PRIVATE :: Chk_CreateMapping
    PRIVATE :: Chk_GetConcentrations
    PRIVATE :: Chk_GetDensity
    PRIVATE :: Chk_GetEndCell
    PRIVATE :: Chk_GetGfw
    PRIVATE :: Chk_GetSaturation
    PRIVATE :: Chk_GetSelectedOutput
    PRIVATE :: Chk_GetSolutionVolume
    PRIVATE :: Chk_GetSpeciesConcentrations
    PRIVATE :: Chk_GetSpeciesD25
    PRIVATE :: Chk_GetSpeciesZ
    PRIVATE :: Chk_GetStartCell
    PRIVATE :: Chk_InitialPhreeqc2Concentrations
    PRIVATE :: Chk_InitialPhreeqc2Module
    PRIVATE :: Chk_InitialPhreeqcCell2Module
    PRIVATE :: Chk_InitialPhreeqc2SpeciesConcentrations
    PRIVATE :: Chk_SetConcentrations
    PRIVATE :: Chk_SetDensity
    PRIVATE :: Chk_SetPorosity
    PRIVATE :: Chk_SetPressure
    PRIVATE :: Chk_SetPrintChemistryMask
    PRIVATE :: Chk_SetSaturation
    PRIVATE :: Chk_SetTemperature
    PRIVATE :: Chk_SpeciesConcentrations2Module
    PRIVATE :: Chk_Double1D
    PRIVATE :: Chk_Double2D
    PRIVATE :: Chk_Integer1D
    PRIVATE :: Chk_Integer2D
    CONTAINS
    
!> Abort the program. 
!> @a irm_result will be interpreted as
!> an IRM_RESULT value and decoded; @a err_str will be printed; and the reaction module
!> will be destroyed. If using MPI, an MPI_Abort message will be sent before the reaction
!> module is destroyed. If the @a id is an invalid instance, RM_Abort will return a value of
!> IRM_BADINSTANCE, otherwise the program will exit with a return code of 4.
!> @param id            The instance id returned from @ref RM_Create.
!> @param irm_result        Integer treated as an IRM_RESULT return code.
!> @param err_str       String to be printed as an error message.
!> @retval IRM_RESULT   Program will exit before returning unless @a id is an invalid reaction module id.
!> @see                 @ref RM_Destroy, @ref RM_ErrorMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm)
!> string = "SELECTED_OUTPUT 5; -pH;RUN_CELLS; -cells 1"
!> status = RunString(iphreeqc_id, string)
!> if (status .ne. 0) status = RM_Abort(id, status, "IPhreeqc RunString failed")
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root or workers.

INTEGER FUNCTION RM_Abort(id, irm_result, err_str)
	USE ISO_C_BINDING
    IMPLICIT NONE 
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_Abort(id, result, str) &
			BIND(C, NAME='RMF_Abort')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: result
            CHARACTER(KIND=C_CHAR), INTENT(in) :: str(*)
        END FUNCTION RMF_Abort
    END INTERFACE     
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: irm_result
    CHARACTER(len=*), INTENT(in) :: err_str
    RM_Abort = RMF_Abort(id, irm_result, trim(err_str)//C_NULL_CHAR)
    RETURN    
END FUNCTION RM_Abort

!> Close the output and log files.
!> @param id            The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                 @ref RM_OpenFiles, @ref RM_SetFilePrefix
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_CloseFiles(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called only by root.

INTEGER FUNCTION RM_CloseFiles(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_CloseFiles(id) &
			BIND(C, NAME='RMF_CloseFiles')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_CloseFiles
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_CloseFiles = RMF_CloseFiles(id)
    RETURN    
END FUNCTION RM_CloseFiles

!> @a N sets of component concentrations are converted to SOLUTIONs numbered 1-@a n in the Utility IPhreeqc.
!> The solutions can be reacted and manipulated with the methods of IPhreeqc. If solution concentration units
!> (@ref RM_SetUnitsSolution) are per liter, one liter of solution is created in the Utility instance; if solution
!> concentration units are mass fraction, one kilogram of solution is created in the Utility instance.
!> The motivation for this
!> method is the mixing of solutions in wells, where it may be necessary to calculate solution properties
!> (pH for example) or react the mixture to form scale minerals. The code fragments below make a mixture of
!> concentrations and then calculate the pH of the mixture.
!> @param id            The instance @a id returned from @ref RM_Create.
!> @param c             Array of concentrations to be made SOLUTIONs in Utility IPhreeqc, array size is 
!> (@a n, @a ncomps) where @a ncomps is the number of components (@ref RM_GetComponentCount).
!> @param n             The number of sets of concentrations.
!> @param tc            Array of temperatures to apply to the SOLUTIONs, in degree C. Array of size @a n.
!> @param p_atm         Array of pressures to apply to the SOLUTIONs, in atm. Array of size n.
!> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate (c_well(1,ncomps))
!> do i = 1, ncomps
!>   c_well(1,i) = 0.5 * c(1,i) + 0.5 * c(10,i)
!> enddo
!> allocate(tc(1), p_atm(1))
!> tc(1) = 15.0
!> p_atm(1) = 3.0
!> iphreeqc_id = RM_Concentrations2Utility(id, c_well, 1, tc, p_atm)
!> string = "SELECTED_OUTPUT 5; -pH; RUN_CELLS; -cells 1"
!> status = RunString(iphreeqc_id, string)
!> status = SetCurrentSelectedOutputUserNumber(iphreeqc_id, 5)
!> status = GetSelectedOutputValue(iphreeqc_id, 1, 1, vtype, pH, svalue)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called only by root.

INTEGER FUNCTION RM_Concentrations2Utility(id, c, n, tc, p_atm)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_Concentrations2Utility(id, c, n, tc, p_atm) &
			BIND(C, NAME='RMF_Concentrations2Utility')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: c(*)
            INTEGER(KIND=C_INT), INTENT(in) :: n
            REAL(KIND=C_DOUBLE), INTENT(in) :: tc(*), p_atm(*)
        END FUNCTION RMF_Concentrations2Utility  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: c
    INTEGER, INTENT(in) :: n
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: tc, p_atm
    if (rmf_debug) CALL ChK_Concentrations2Utility(id, c, n, tc, p_atm)
    RM_Concentrations2Utility = RMF_Concentrations2Utility(id, c, n, tc, p_atm)
    return
END FUNCTION RM_Concentrations2Utility

SUBROUTINE ChK_Concentrations2Utility(id, c, n, tc, p_atm)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: c
    INTEGER, INTENT(in) :: n
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: tc, p_atm
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double2D(id, c, n, rmf_ncomps, "Concentration", "RM_Concentrations2Utility")
    errors = errors + Chk_Double1D(id, tc, n, "Temperature", "RM_Concentrations2Utility")
    errors = errors + Chk_Double1D(id, p_atm, n, "Pressure", "RM_Concentrations2Utility")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument(s) in RM_Concentrations2Utility")
    endif
END SUBROUTINE Chk_Concentrations2Utility  

!> Creates a reaction module. If the code is compiled with
!> the preprocessor directive USE_OPENMP, the reaction module is multithreaded.
!> If the code is compiled with the preprocessor directive USE_MPI, the reaction
!> module will use MPI and multiple processes. If neither preprocessor directive is used,
!> the reaction module will be serial (unparallelized).
!> @param nxyz                   The number of grid cells in the user's model.
!> @param nthreads (or @a comm, MPI)       When using OPENMP, the argument (@a nthreads) is the number of worker threads to be used.
!> If @a nthreads <= 0, the number of threads is set equal to the number of processors of the computer.
!> When using MPI, the argument (@a comm) is the MPI communicator to use within the reaction module.
!> @retval Id of the PhreeqcRM instance, negative is failure (See @ref RM_DecodeError).
!> @see                 @ref RM_Destroy
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> nxyz = 40
!> #ifdef USE_MPI
!>   id = RM_Create(nxyz, MPI_COMM_WORLD)
!>   call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
!>   if (status .ne. MPI_SUCCESS) then
!>     stop "Failed to get mpi_myself"
!>   endif
!>   if (mpi_myself > 0) then
!>     status = RM_MpiWorker(id)
!>     status = RM_Destroy(id)
!>     return
!>   endif
!> #else
!>   nthreads = 3
!>   id = RM_Create(nxyz, nthreads)
!> #endif
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and workers. 

INTEGER FUNCTION RM_Create(nxyz, nthreads) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_Create(nxyz, nthreads) &
			BIND(C, NAME='RMF_Create') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: nxyz
			INTEGER(KIND=C_INT), INTENT(in) :: nthreads
        END FUNCTION RMF_Create
	END INTERFACE
    INTEGER, INTENT(in) :: nxyz
	INTEGER, INTENT(in) :: nthreads
    RM_Create = RMF_Create(nxyz, nthreads) 
    rmf_nxyz = nxyz
  
    return
END FUNCTION RM_Create

!> Provides a mapping from grid cells in the user's model to reaction cells in PhreeqcRM.
!> The mapping is used to eliminate inactive cells and to use symmetry to decrease the number of cells for which chemistry must be run.
!> The mapping may be many-to-one to account for symmetry.
!> Default is a one-to-one mapping--all user grid cells are reaction cells (equivalent to @a grid2chem values of 0,1,2,3,...,@a nxyz-1).
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param grid2chem        An array of integers: Nonnegative is a reaction cell number (0 based), negative is an inactive cell. Array of size @a nxyz (number of grid cells).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ! For demonstation, two equivalent rows by symmetry
!> allocate(grid2chem(nxyz))
!> do i = 1, nxyz/2
!>   grid2chem(i) = i - 1
!>   grid2chem(i+nxyz/2) = i - 1
!> enddo
!> status = RM_CreateMapping(id, grid2chem)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_CreateMapping(id, grid2chem)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION RMF_CreateMapping(id, grid2chem) &
			BIND(C, NAME='RMF_CreateMapping')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: grid2chem(*)
        END FUNCTION RMF_CreateMapping
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:) :: grid2chem
    if (rmf_debug) call Chk_CreateMapping(id, grid2chem)
    RM_CreateMapping = RMF_CreateMapping(id, grid2chem)
    return
END FUNCTION RM_CreateMapping

SUBROUTINE Chk_CreateMapping(id, grid2chem)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:) :: grid2chem
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Integer1D(id, grid2chem, rmf_nxyz, "Grid2chem mapping", "RM_CreateMapping")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument(s) in RM_CreateMapping")
    endif
END SUBROUTINE Chk_CreateMapping

!> If @a e is negative, this method prints an error message corresponding to IRM_RESULT @a e. If @a e is non-negative, no action is taken.
!> @param id                   The instance @a id returned from @ref RM_Create.
!> @param e                    An IRM_RESULT value returned by one of the reaction-module methods.
!> @retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
!> @par IRM_RESULT definition:
!> @htmlonly
!> <CODE>
!> <PRE>
!> typedef enum {
!>   IRM_OK            =  0,  //Success
!>   IRM_OUTOFMEMORY   = -1,  //Failure, Out of memory
!>   IRM_BADVARTYPE    = -2,  //Failure, Invalid VAR type
!>   IRM_INVALIDARG    = -3,  //Failure, Invalid argument
!>   IRM_INVALIDROW    = -4,  //Failure, Invalid row
!>   IRM_INVALIDCOL    = -5,  //Failure, Invalid column
!>   IRM_BADINSTANCE   = -6,  //Failure, Invalid rm instance id
!>   IRM_FAIL          = -7,  //Failure, Unspecified
!> } IRM_RESULT;
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_CreateMapping(id, grid2chem)
!> if (status < 0) status = RM_DecodeError(id, status)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Can be called by root and (or) workers.

INTEGER FUNCTION RM_DecodeError(id, e)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
		INTEGER(KIND=C_INT) FUNCTION RMF_DecodeError(id, e) &
			BIND(C, NAME='RMF_DecodeError')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: e
        END FUNCTION RMF_DecodeError
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: e
    RM_DecodeError = RMF_DecodeError(id, e)
    return
END FUNCTION RM_DecodeError

!> Destroys a reaction module.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT   0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_Create
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_Destroy(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and workers.

INTEGER FUNCTION RM_Destroy(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_Destroy(id) &
			BIND(C, NAME='RMF_Destroy')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_Destroy
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_Destroy = RMF_Destroy(id)
    return
END FUNCTION RM_Destroy

!> Writes the contents of all workers to file in _RAW formats, including SOLUTIONs and all reactants.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param dump_on          Signal for writing the dump file: 1 true, 0 false.
!> @param append           Signal to append to the contents of the dump file: 1 true, 0 false.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetDumpFileName
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>      dump_on = 1
!> append = 0
!> status = RM_SetDumpFileName(id, "advection_f90.dmp")
!> status = RM_DumpModule(id, dump_on, append)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root; workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_DumpModule(id, dump_on, append) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_DumpModule(id, dump_on, append) &
			BIND(C, NAME='RMF_DumpModule')
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: dump_on
            INTEGER(KIND=C_INT), INTENT(in) :: append
        END FUNCTION RMF_DumpModule
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: dump_on
    INTEGER, INTENT(in) :: append
    RM_DumpModule = RMF_DumpModule(id, dump_on, append)
    return
END FUNCTION RM_DumpModule

!> Send an error message to the screen, the output file, and the log file.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param errstr           String to be printed.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_OpenFiles, @ref RM_LogMessage, @ref RM_OutputMessage, @ref RM_ScreenMessage, @ref RM_WarningMessage.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_ErrorMessage(id, "Goodbye world")
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers; root writes to output and log files.

INTEGER FUNCTION RM_ErrorMessage(id, errstr)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_ErrorMessage(id, errstr) &
			BIND(C, NAME='RMF_ErrorMessage')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: errstr(*)
        END FUNCTION RMF_ErrorMessage
	END INTERFACE 
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: errstr
    RM_ErrorMessage = RMF_ErrorMessage(id, trim(errstr)//C_NULL_CHAR)
    return
END FUNCTION RM_ErrorMessage

!> Returns the number of items in the list of all elements in the InitialPhreeqc instance.
!> Elements are those that have been defined in a solution or any other reactant (EQUILIBRIUM_PHASE, KINETICS, and others).
!> The method can be called multiple times and the list that is created is cummulative.
!> The list is the set of components that needs to be transported. 
!> By default the list
!> includes water, excess H and excess O (the H and O not contained in water);
!> alternatively, the list may be set to contain total H and total O (@ref RM_SetComponentH2O),
!> which requires transport results to be accurate to eight or nine significant digits.
!> If multicomponent diffusion (MCD) is to be modeled, there is a capability to retrieve aqueous species concentrations
!> (@ref RM_GetSpeciesConcentrations) and to set new solution concentrations after MCD by using individual species concentrations
!> (@ref RM_SpeciesConcentrations2Module). To use these methods the save-species property needs to be turned on (@ref RM_SetSpeciesSaveOn).
!> If the save-species property is on, RM_FindComponents will generate
!> a list of aqueous species (@ref RM_GetSpeciesCount, @ref RM_GetSpeciesName), their diffusion coefficients at 25 C (@ref RM_GetSpeciesD25),
!> their charge (@ref RM_GetSpeciesZ).
!> @param id            The instance @a id returned from @ref RM_Create.
!> @retval              Number of components currently in the list, or IRM_RESULT error code (see @ref RM_DecodeError).
!> @see                 @ref RM_GetComponent, @ref RM_SetSpeciesSaveOn, @ref RM_GetSpeciesConcentrations, 
!> @ref RM_SpeciesConcentrations2Module, @ref RM_GetSpeciesCount, @ref RM_GetSpeciesName, 
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesZ, @ref RM_SetComponentH2O.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ! Get list of components
!> ncomps = RM_FindComponents(id)
!> allocate(components(ncomps))
!> do i = 1, ncomps
!>   status = RM_GetComponent(id, i, components(i))
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_FindComponents(id) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_FindComponents(id) &
			BIND(C, NAME='RMF_FindComponents') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_FindComponents  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_FindComponents = RMF_FindComponents(id)
    rmf_ncomps = RM_FindComponents
    return
END FUNCTION RM_FindComponents  

!> Fills an array with the cell numbers in the user's numbering sytstem that map to a cell in the
!> PhreeqcRM numbering system. The mapping is defined by @ref RM_CreateMapping.

!> @param id            The instance @a id returned from @ref RM_Create.
!> @param n             A cell number in the PhreeqcRM numbering system (0 <= n < @ref RM_GetChemistryCellCount).
!> @param list          Array to store the user cell numbers mapped to PhreeqcRM cell @a n.
!> @param size          Input, the allocated size of @a list; it is an error if the array is too small. 
!>                      Output, the number of cells mapped to cell @a n.
!> @retval              IRM_RESULT error code (see @ref RM_DecodeError).
!> 
!> @see                 @ref RM_CreateMapping, @ref RM_GetChemistryCellCount, @ref RM_GetGridCellCount.
!> 
!> @par C Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> if (RM_GetBackwardMapping(rm_id, rm_cell_number, list, size) .eq. 0) then
!>   if (fstr(1:l) .eq. "HYDRAULIC_K") then
!>     my_basic_fortran_callback = K_ptr(list(1)+1) 
!>   endif
!> endif
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetBackwardMapping(id, n, list, size) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetBackwardMapping(id, n, list, size) &
			BIND(C, NAME='RMF_GetBackwardMapping') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in)    :: id, n
            INTEGER(KIND=C_INT), INTENT(in)    :: list(*)
            INTEGER(KIND=C_INT), INTENT(inout) :: size
        END FUNCTION RMF_GetBackwardMapping  
	END INTERFACE
    INTEGER, INTENT(in)    :: id, n
    INTEGER, INTENT(in)    :: list(*)
    INTEGER, INTENT(inout) :: size
    RM_GetBackwardMapping = RMF_GetBackwardMapping(id, n, list, size)
    return
END FUNCTION RM_GetBackwardMapping  

!> Returns the number of chemistry cells in the reaction module. The number of chemistry cells is defined by
!> the set of non-negative integers in the mapping from user grid cells (@ref RM_CreateMapping).
!> The number of chemistry cells is less than or equal to the number of cells in the user's model.
!> @param id            The instance @a id returned from @ref RM_Create.
!> @retval              Number of chemistry cells, or IRM_RESULT error code (see @ref RM_DecodeError).
!> @see                 @ref RM_CreateMapping, @ref RM_GetGridCellCount.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_CreateMapping(id, grid2chem)
!> nchem = RM_GetChemistryCellCount(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetChemistryCellCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetChemistryCellCount(id) &
			BIND(C, NAME='RMF_GetChemistryCellCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetChemistryCellCount 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetChemistryCellCount = RMF_GetChemistryCellCount(id)
    return
END FUNCTION RM_GetChemistryCellCount 

!> Retrieves an item from the reaction-module component list that was generated by calls to @ref RM_FindComponents.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param num              The number of the component to be retrieved. Fortran, 1 based.
!> @param comp_name        The string value associated with component @a num.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetComponentCount
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ! Get list of components
!> ncomps = RM_FindComponents(id)
!> allocate(components(ncomps))
!> do i = 1, ncomps
!>   status = RM_GetComponent(id, i, components(i))
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetComponent(id, num, comp_name)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetComponent(id, num, comp_name, l) &
			BIND(C, NAME='RMF_GetComponent')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, num, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: comp_name(*)
        END FUNCTION RMF_GetComponent 
	END INTERFACE
    INTEGER, INTENT(in) :: id, num
    CHARACTER(len=*), INTENT(inout) :: comp_name
    RM_GetComponent = RMF_GetComponent(id, num, comp_name, len(comp_name))
    return
END FUNCTION RM_GetComponent 

!> Returns the number of components in the reaction-module component list. 
!> The component list is generated by calls to @ref RM_FindComponents.
!> The return value from the last call to @ref RM_FindComponents is equal to the return value from RM_GetComponentCount.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 The number of components in the reaction-module component list, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetComponent.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ncomps1 = RM_GetComponentCount(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_GetComponentCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetComponentCount(id) &
			BIND(C, NAME='RMF_GetComponentCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetComponentCount 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetComponentCount = RMF_GetComponentCount(id)
END FUNCTION RM_GetComponentCount 

!> Transfer solution concentrations from each reaction cell 
!> to the concentration array given in the argument list (@a c).
!> Units of concentration for @a c are defined by @ref RM_SetUnitsSolution. 
!> For concentration units of per liter, 
!> the solution volume is used to calculate the concentrations for @a c. 
!> For mass fraction concentration units, 
!> the solution mass is used to calculate concentrations for @a c.
!> Two options are available for the volume and mass of solution 
!> that are used in converting to transport concentrations: (1) the volume and mass of solution are
!> calculated by PHREEQC, or 
!> (2) the volume of solution is the product of saturation (@ref RM_SetSaturation),
!> porosity (@ref RM_SetPorosity), and representative volume (@ref RM_SetRepresentativeVolume),
!> and the mass of solution is volume times density as defined by @ref RM_SetDensity.
!> @ref RM_UseSolutionDensityVolume determines which option is used.
!> For option 1, the databases that have partial molar volume definitions needed
!> to accurately calculate solution volume are
!> phreeqc.dat, Amm.dat, and pitzer.dat. 
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param c                Array to receive the concentrations. Dimension of the array is (@a nxyz, @a ncomps),
!> where @a nxyz is the number of user grid cells and @a ncomps is the result of @ref RM_FindComponents or @ref RM_GetComponentCount.
!> Values for inactive cells are set to 1e30.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> 
!> @see                    @ref RM_FindComponents, @ref RM_GetComponentCount, @ref RM_GetSaturation,
!> @ref RM_SetConcentrations, @ref RM_SetDensity, @ref RM_SetRepresentativeVolume, @ref RM_SetSaturation,
!> @ref RM_SetUnitsSolution, @ref RM_UseSolutionDensityVolume.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(c(nxyz, ncomps))
!> status = RM_RunCells(id)
!> status = RM_GetConcentrations(id, c)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetConcentrations(id, c) 
	USE ISO_C_BINDING  
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetConcentrations(id, c) &
			BIND(C, NAME='RMF_GetConcentrations')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out)  :: c(*)
        END FUNCTION RMF_GetConcentrations 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:,:) :: c
    if (rmf_debug) call Chk_GetConcentrations(id, c)  
    RM_GetConcentrations = RMF_GetConcentrations(id, c)   
    return
END FUNCTION RM_GetConcentrations         

SUBROUTINE Chk_GetConcentrations(id, c)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: c
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double2D(id, c, rmf_nxyz, rmf_ncomps, "concentration", "RM_GetConcentrations")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument(s) in RM_GetConcentrations")
    endif
END SUBROUTINE Chk_GetConcentrations
#ifdef SKIP
INTEGER FUNCTION RM_GetConcentrations1D(id, c) 
	USE ISO_C_BINDING  
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetConcentrations(id, c) &
			BIND(C, NAME='RMF_GetConcentrations')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out)  :: c(*)
        END FUNCTION RMF_GetConcentrations 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:,:) :: c
    RM_GetConcentrations1D = RMF_GetConcentrations(id, c)   
    return
END FUNCTION RM_GetConcentrations1D    
#endif
!> Transfer solution densities from the reaction cells to the array given in the argument list (@a density). 
!> Densities are those calculated by the reaction module.
!> Only the following databases distributed with PhreeqcRM have molar volume information needed to accurately calculate density:
!> phreeqc.dat, Amm.dat, and pitzer.dat.
!> 
!> @param id                   The instance @a id returned from @ref RM_Create.
!> @param density              Array to receive the densities. Dimension of the array is @a nxyz,
!> where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount). Values for inactive cells are set to 1e30.
!> @retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(density(nxyz))
!> status = RM_RunCells(id)
!> status = RM_GetDensity(id, density)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetDensity(id, density)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetDensity(id, density) &
			BIND(C, NAME='RMF_GetDensity')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: density(*)
        END FUNCTION RMF_GetDensity 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), dimension(:) :: density
    if (rmf_debug) call Chk_GetDensity(id, density)
    RM_GetDensity = RMF_GetDensity(id, density) 
    return
END FUNCTION RM_GetDensity 

SUBROUTINE Chk_GetDensity(id, density)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: density
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, density, rmf_nxyz, "density", "RM_GetDensity")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetDensity")
    endif
END SUBROUTINE Chk_GetDensity

!> Returns an array with the ending cell numbers from the range of cell numbers assigned to each worker.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param ec               Array to receive the ending cell numbers. Dimension of the array is 
!>                         the number of threads (OpenMP) or the number of processes (MPI).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_Create, @ref RM_GetMpiTasks, @ref RM_GetStartCell, @ref RM_GetThreadCount.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> n = RM_GetThreadCount(id) * RM_GetMpiTasks(id)
!> allocate(ec(n))
!> status = RM_GetEndCell(id, ec)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetEndCell(id, ec)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetEndCell(id, ec) &
			BIND(C, NAME='RMF_GetEndCell')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(out):: ec(*)
        END FUNCTION RMF_GetEndCell 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(out), DIMENSION(:) :: ec
    if (rmf_debug) call Chk_GetEndCell(id, ec)
    RM_GetEndCell = RMF_GetEndCell(id, ec)
    RETURN
END FUNCTION RM_GetEndCell

SUBROUTINE Chk_GetEndCell(id, ec)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:) :: ec
    INTEGER :: errors, n
    errors = 0
    n = RM_GetMpiTasks(id) * RM_GetThreadCount(id)
    errors = errors + Chk_Integer1D(id, ec, n, "EndCell", "RM_GetEndCell")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetEndCell")
    endif
END SUBROUTINE Chk_GetEndCell

!> Returns a string containing error messages related to the last call to a PhreeqcRM method to 
!> the character argument (@a errstr).
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param errstr           The error string related to the last call to a PhreeqcRM method.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> character(len=:), allocatable :: errstr
!> integer l
!> if (status .ne. 0) then
!>   l = RM_GetErrorStringLength(id)
!>   allocate (character(len=l) :: errstr)
!>   status = RM_GetErrorString(id, errstr)
!>   write(*,"(A)") errstr
!>   deallocate(errstr)
!>   status = RM_Destroy(id)
!>   stop
!> endif </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetErrorString(id, errstr)  
	USE ISO_C_BINDING 
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetErrorString(id, errstr, l) &
			BIND(C, NAME='RMF_GetErrorString')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
			INTEGER(KIND=C_INT), INTENT(in) :: l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: errstr(*)
        END FUNCTION RMF_GetErrorString 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: errstr
    RM_GetErrorString = RMF_GetErrorString(id, errstr, len(errstr))   
END FUNCTION RM_GetErrorString 

!> Returns the length of the string that contains error messages related to the last call to a PhreeqcRM method. 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval int             Length of the error message string.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> character(len=:), allocatable :: errstr
!> integer l
!> if (status .ne. 0) then
!>   l = RM_GetErrorStringLength(id)
!>   allocate (character(len=l) :: errstr)
!>   status = RM_GetErrorString(id, errstr)
!>   write(*,"(A)") errstr
!>   deallocate(errstr)
!>   status = RM_Destroy(id)
!>   stop
!> endif 
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetErrorStringLength(id)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetErrorStringLength(id) &
			BIND(C, NAME='RMF_GetErrorStringLength')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetErrorStringLength 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetErrorStringLength = RMF_GetErrorStringLength(id) 
END FUNCTION RM_GetErrorStringLength 
 
!> Returns the reaction-module file prefix to the character argument (@a prefix).
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param prefix           Character string where the prefix is written.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetFilePrefix.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> character(100) :: string
!> character(200) :: string1
!> status = RM_GetFilePrefix(id, string)
!> string1 = "File prefix: "//string
!> status = RM_OutputMessage(id, string1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetFilePrefix(id, prefix)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetFilePrefix(id, prefix, l) &
			BIND(C, NAME='RMF_GetFilePrefix')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: prefix(*)
        END FUNCTION RMF_GetFilePrefix
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: prefix
    integer l
    l = len(prefix)
    RM_GetFilePrefix = RMF_GetFilePrefix(id, prefix, l)
END FUNCTION RM_GetFilePrefix

!> Returns the gram formula weights (g/mol) for the components in the reaction-module component list.
!> @param id               The instance id returned from @ref RM_Create.
!> @param gfw              Array to receive the gram formula weights. Dimension of the array is @a ncomps,
!> where @a ncomps is the number of components in the component list.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetComponentCount, @ref RM_GetComponent.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> character(100),   dimension(:), allocatable   :: components
!> double precision, dimension(:), allocatable   :: gfw
!> ncomps = RM_FindComponents(id)
!> allocate(components(ncomps))
!> allocate(gfw(ncomps))
!> status = RM_GetGfw(id, gfw)
!> do i = 1, ncomps
!>   status = RM_GetComponent(id, i, components(i))
!>   write(string,"(A10, F15.4)") components(i), gfw(i)
!>   status = RM_OutputMessage(id, string)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_GetGfw(id, gfw)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetGfw(id, gfw) &
			BIND(C, NAME='RMF_GetGfw')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: gfw(*)
        END FUNCTION RMF_GetGfw 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: gfw
    if (rmf_debug) call Chk_GetGfw(id, gfw) 
    RM_GetGfw = RMF_GetGfw(id, gfw)   
END FUNCTION RM_GetGfw 

SUBROUTINE Chk_GetGfw(id, gfw) 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: gfw
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, gfw, rmf_ncomps, "gfw", "RM_GetGfw")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetGfw")
    endif
END SUBROUTINE Chk_GetGfw

!> Returns the number of grid cells in the user's model, which is defined in the call to @ref RM_Create.
!> The mapping from grid cells to reaction cells is defined by @ref RM_CreateMapping.
!> The number of reaction cells may be less than the number of grid cells if
!> there are inactive regions or symmetry in the model definition.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 Number of grid cells in the user's model, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_Create,  @ref RM_CreateMapping.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> nxyz = RM_GetGridCellCount(id)
!> write(string1, "(A,I)") "Number of grid cells in the user's model: ", nxyz
!> status = RM_OutputMessage(id, trim(string1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetGridCellCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetGridCellCount(id) &
			BIND(C, NAME='RMF_GetGridCellCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetGridCellCount
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetGridCellCount = RMF_GetGridCellCount(id)
END FUNCTION RM_GetGridCellCount

!> Returns an IPhreeqc id for the @a ith IPhreeqc instance in the reaction module.
!> For the threaded version, there are @a nthreads + 2 IPhreeqc instances, where
!> @a nthreads is defined in the constructor (@ref RM_Create).
!> The number of threads can be determined by @ref RM_GetThreadCount.
!> The first @a nthreads (0 based) instances will be the workers, the
!> next (@a nthreads) is the InitialPhreeqc instance, and the next (@a nthreads + 1) is the Utility instance.
!> Getting the IPhreeqc pointer for one of these instances allows the user to use any of the IPhreeqc methods
!> on that instance.
!> For MPI, each process has exactly three IPhreeqc instances, one worker (number 0),
!> one InitialPhreeqc instance (number 1), and one Utility instance (number 2).
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param i                The number of the IPhreeqc instance to be retrieved (0 based).
!> @retval                 IPhreeqc id for the @a ith IPhreeqc instance, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_Create, @ref RM_GetThreadCount. 
!> See IPhreeqc documentation for descriptions of IPhreeqc methods.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> ! Utility pointer is worker number nthreads + 1
!> iphreeqc_id1 = RM_GetIPhreeqcId(id, RM_GetThreadCount(id) + 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetIPhreeqcId(id, i)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetIPhreeqcId(id, i) &
			BIND(C, NAME='RMF_GetIPhreeqcId')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: i
        END FUNCTION RMF_GetIPhreeqcId
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: i
    RM_GetIPhreeqcId = RMF_GetIPhreeqcId(id, i)
END FUNCTION RM_GetIPhreeqcId

!> Returns the MPI task number. For the OPENMP version, the task number is always
!> zero and the result of @ref RM_GetMpiTasks is one. For the MPI version,
!> the root task number is zero, and all workers have a task number greater than zero.
!> The number of tasks can be obtained with @ref RM_GetMpiTasks. The number of
!> tasks and computer hosts are determined at run time by the mpiexec command, and the
!> number of reaction-module processes is defined by the communicator used in
!> constructing the reaction modules (@ref RM_Create).
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 The MPI task number for a process, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetMpiTasks.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string1, "(A,I)") "MPI task number: ", RM_GetMpiMyself(id)
!> status = RM_OutputMessage(id, string1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetMpiMyself(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetMpiMyself(id) &
			BIND(C, NAME='RMF_GetMpiMyself')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetMpiMyself
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetMpiMyself = RMF_GetMpiMyself(id)
END FUNCTION RM_GetMpiMyself

!> Returns the number of MPI processes (tasks) assigned to the reaction module.
!> For the OPENMP version, the number of tasks is always
!> one (although there may be multiple threads, @ref RM_GetThreadCount),
!> and the task number returned by @ref RM_GetMpiMyself is zero. For the MPI version, the number of
!> tasks and computer hosts are determined at run time by the mpiexec command. An MPI communicator
!> is used in constructing reaction modules for MPI. The communicator may define a subset of the
!> total number of MPI processes.
!> The root task number is zero, and all workers have a task number greater than zero.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 The number of MPI  processes assigned to the reaction module,
!> negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetMpiMyself, @ref RM_Create.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> mpi_tasks = RM_GetMpiTasks(id)
!> write(string1, "(A,I)") "Number of MPI processes: ", mpi_tasks
!> status = RM_OutputMessage(id, string1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetMpiTasks(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetMpiTasks(id) &
			BIND(C, NAME='RMF_GetMpiTasks')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetMpiTasks
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetMpiTasks = RMF_GetMpiTasks(id)
END FUNCTION RM_GetMpiTasks

!> Returns the user number for the @a nth selected-output definition.
!> Definitions are sorted by user number. Phreeqc allows multiple selected-output
!> definitions, each of which is assigned a nonnegative integer identifier by the
!> user. The number of definitions can be obtained by @ref RM_GetSelectedOutputCount.
!> To cycle through all of the definitions, RM_GetNthSelectedOutputUserNumber
!> can be used to identify the user number for each selected-output definition
!> in sequence. @ref RM_SetCurrentSelectedOutputUserNumber is then used to select
!> that user number for selected-output processing.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param n                The sequence number of the selected-output definition for which the user number will be returned.
!> Fortran, 1 based.
!> @retval                 The user number of the @a nth selected-output definition, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetSelectedOutput,
!> @ref RM_GetSelectedOutputColumnCount, @ref RM_GetSelectedOutputCount,
!> @ref RM_GetSelectedOutputHeading,
!> @ref RM_GetSelectedOutputRowCount, @ref RM_SetCurrentSelectedOutputUserNumber, @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   write(*,*) "Selected output sequence number: ", isel)
!>   write(*,*) "Selected output user number:     ", n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   allocate(selected_out(nxyz,col))
!>   status = RM_GetSelectedOutput(id, selected_out)
!>   ! Process results here
!>   deallocate(selected_out)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_GetNthSelectedOutputUserNumber(id, n)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetNthSelectedOutputUserNumber(id, n) &
			BIND(C, NAME='RMF_GetNthSelectedOutputUserNumber')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n
        END FUNCTION RMF_GetNthSelectedOutputUserNumber 
	END INTERFACE
    INTEGER, INTENT(in) :: id, n
    RM_GetNthSelectedOutputUserNumber = RMF_GetNthSelectedOutputUserNumber(id, n)
END FUNCTION RM_GetNthSelectedOutputUserNumber 

!> Returns a vector of saturations (@a sat_calc) as calculated by the reaction module.
!> Reactions will change the volume of solution in a cell.
!> The transport code must decide whether to ignore or account for this change in solution volume due to reactions.
!> Following reactions, the cell saturation is calculated as solution volume (@ref RM_GetSolutionVolume)
!> divided by the product of representative volume (@ref RM_SetRepresentativeVolume) and the porosity (@ref RM_SetPorosity).
!> The cell saturation returned by @a RM_GetSaturation may be less than or greater than the saturation set by the transport code
!> (@ref RM_SetSaturation), and may be greater than or less than 1.0, even in fully saturated simulations.
!> Only the following databases distributed with PhreeqcRM have molar volume information needed
!> to accurately calculate solution volume and saturation: phreeqc.dat, Amm.dat, and pitzer.dat.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param sat_calc              Vector to receive the saturations. Dimension of the array is set to @a nxyz,
!> where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount).
!> Values for inactive cells are set to 1e30.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> 
!> @see                    @ref RM_GetSolutionVolume, @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume, @ref RM_SetSaturation.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(sat_calc(nxyz))
!> status = RM_RunCells(id)
!> status = RM_GetSaturation(id, sat_calc)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetSaturation(id, sat_calc)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSaturation(id, sat_calc) &
			BIND(C, NAME='RMF_GetSaturation')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: sat_calc(*)
        END FUNCTION RMF_GetSaturation
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:) :: sat_calc
    if (rmf_debug) call Chk_GetSaturation(id, sat_calc)
    RM_GetSaturation = RMF_GetSaturation(id, sat_calc)
END FUNCTION RM_GetSaturation

SUBROUTINE Chk_GetSaturation(id, sat)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: sat
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, sat, rmf_nxyz, "saturation", "RM_GetSaturation")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetSaturation")
    endif
END SUBROUTINE Chk_GetSaturation

!> Populates an array with values from the current selected-output definition. @ref RM_SetCurrentSelectedOutputUserNumber
!> determines which of the selected-output definitions is used to populate the array.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param so               An array to contain the selected-output values. Size of the array is (@a nxyz, @a col),
!> where @a nxyz is the number of grid cells in the user's model (@ref RM_GetGridCellCount), and @a col is the number of
!> columns in the selected-output definition (@ref RM_GetSelectedOutputColumnCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetNthSelectedOutputUserNumber,
!> @ref RM_GetSelectedOutputColumnCount, @ref RM_GetSelectedOutputCount, @ref RM_GetSelectedOutputHeading,
!> @ref RM_GetSelectedOutputRowCount, @ref RM_SetCurrentSelectedOutputUserNumber, @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   allocate(selected_out(nxyz,col))
!>   status = RM_GetSelectedOutput(id, selected_out)
!>   ! Process results here
!>   deallocate(selected_out)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetSelectedOutput(id, so)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSelectedOutput(id, so) &
			BIND(C, NAME='RMF_GetSelectedOutput')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: so(*)
        END FUNCTION RMF_GetSelectedOutput
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(out) :: so
    if (rmf_debug) call Chk_GetSelectedOutput(id, so)
    RM_GetSelectedOutput = RMF_GetSelectedOutput(id, so)
END FUNCTION RM_GetSelectedOutput
 
SUBROUTINE Chk_GetSelectedOutput(id, so)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: so
    INTEGER :: errors, ncol
    ncol = RM_GetSelectedOutputColumnCount(id)
    errors = 0
    errors = errors + Chk_Double2D(id, so, rmf_nxyz, ncol, "selected output", "RM_GetSelectedOutput")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetSelectedOutput")
    endif
END SUBROUTINE Chk_GetSelectedOutput

!> Returns the number of columns in the current selected-output definition. @ref RM_SetCurrentSelectedOutputUserNumber
!> determines which of the selected-output definitions is used.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 Number of columns in the current selected-output definition, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetNthSelectedOutputUserNumber, @ref RM_GetSelectedOutput,
!> @ref RM_GetSelectedOutputCount, @ref RM_GetSelectedOutputHeading,
!> @ref RM_GetSelectedOutputRowCount, @ref RM_SetCurrentSelectedOutputUserNumber, @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   allocate(selected_out(nxyz,col))
!>   status = RM_GetSelectedOutput(id, selected_out)
!>   ! Process results here
!>   deallocate(selected_out)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_GetSelectedOutputColumnCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSelectedOutputColumnCount(id) &
			BIND(C, NAME='RMF_GetSelectedOutputColumnCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetSelectedOutputColumnCount
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetSelectedOutputColumnCount = RMF_GetSelectedOutputColumnCount(id)
END FUNCTION RM_GetSelectedOutputColumnCount

!> Returns the number of selected-output definitions. @ref RM_SetCurrentSelectedOutputUserNumber
!> determines which of the selected-output definitions is used.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 Number of selected-output definitions, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetNthSelectedOutputUserNumber, @ref RM_GetSelectedOutput,
!> @ref RM_GetSelectedOutputColumnCount, @ref RM_GetSelectedOutputHeading,
!> @ref RM_GetSelectedOutputRowCount, @ref RM_SetCurrentSelectedOutputUserNumber, @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   allocate(selected_out(nxyz,col))
!>   status = RM_GetSelectedOutput(id, selected_out)
!>   ! Process results here
!>   deallocate(selected_out)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_GetSelectedOutputCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSelectedOutputCount(id) &
			BIND(C, NAME='RMF_GetSelectedOutputCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetSelectedOutputCount
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetSelectedOutputCount = RMF_GetSelectedOutputCount(id)
END FUNCTION RM_GetSelectedOutputCount

!> Returns a selected output heading. The number of headings is determined by @ref RM_GetSelectedOutputColumnCount.
!> @ref RM_SetCurrentSelectedOutputUserNumber
!> determines which of the selected-output definitions is used.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param icol             The sequence number of the heading to be retrieved. Fortran, 1 based.
!> @param heading          A string buffer to receive the heading.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetNthSelectedOutputUserNumber, @ref RM_GetSelectedOutput,
!> @ref RM_GetSelectedOutputColumnCount, @ref RM_GetSelectedOutputCount,
!> @ref RM_GetSelectedOutputRowCount, @ref RM_SetCurrentSelectedOutputUserNumber, @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   do j = 1, col
!>     status = RM_GetSelectedOutputHeading(id, j, heading)
!>     write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ", trim(heading)
!>   enddo
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_GetSelectedOutputHeading(id, icol, heading)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSelectedOutputHeading(id, icol, heading, l) &
			BIND(C, NAME='RMF_GetSelectedOutputHeading')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, icol, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: heading(*)
        END FUNCTION RMF_GetSelectedOutputHeading
	END INTERFACE
    INTEGER, INTENT(in) :: id, icol
    CHARACTER(len=*), INTENT(out) :: heading
    RM_GetSelectedOutputHeading = RMF_GetSelectedOutputHeading(id, icol, heading, len(heading))
END FUNCTION RM_GetSelectedOutputHeading

!> Returns the number of rows in the current selected-output definition. However, the method
!> is included only for convenience; the number of rows is always equal to the number of
!> grid cells in the user's model, and is equal to @ref RM_GetGridCellCount.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 Number of rows in the current selected-output definition, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetNthSelectedOutputUserNumber, @ref RM_GetSelectedOutput, @ref RM_GetSelectedOutputColumnCount,
!> @ref RM_GetSelectedOutputCount, @ref RM_GetSelectedOutputHeading,
!> @ref RM_SetCurrentSelectedOutputUserNumber, @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   allocate(selected_out(nxyz,col))
!>   status = RM_GetSelectedOutput(id, selected_out)
!>   ! Print results
!>   do i = 1, RM_GetSelectedOutputRowCount(id)
!>     write(*,*) "Cell number ", i
!>     write(*,*) "     Selected output: "
!>     do j = 1, col
!>       status = RM_GetSelectedOutputHeading(id, j, heading)
!>       write(*,'(10x,i2,A2,A10,A2,f10.4)') j, " ", trim(heading),": ", selected_out(i,j)
!>     enddo
!>   enddo
!>   deallocate(selected_out)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.
        
INTEGER FUNCTION RM_GetSelectedOutputRowCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSelectedOutputRowCount(id) &
			BIND(C, NAME='RMF_GetSelectedOutputRowCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetSelectedOutputRowCount
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetSelectedOutputRowCount = RMF_GetSelectedOutputRowCount(id)
END FUNCTION RM_GetSelectedOutputRowCount

!> Transfer solution volumes from the reaction cells to the array given in the argument list (@a vol).
!> Solution volumes are those calculated by the reaction module.
!> Only the following databases distributed with PhreeqcRM have molar volume information 
!> needed to accurately calculate solution volume:
!> phreeqc.dat, Amm.dat, and pitzer.dat.
!> 
!> @param id                   The instance @a id returned from @ref RM_Create.
!> @param vol                  Array to receive the solution volumes. Dimension of the array is (@a nxyz),
!> where @a nxyz is the number of user grid cells. Values for inactive cells are set to 1e30. 
!> @retval IRM_RESULT         0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetSaturation.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(volume(nxyz))
!> status = RM_RunCells(id)
!> status = RM_GetSolutionVolume(id, volume)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetSolutionVolume(id, vol)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSolutionVolume(id, vol) &
			BIND(C, NAME='RMF_GetSolutionVolume')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: vol(*)
        END FUNCTION RMF_GetSolutionVolume 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:) :: vol
    if (rmf_debug) call Chk_GetDensity(id, vol)
    RM_GetSolutionVolume = RMF_GetSolutionVolume(id, vol)   
END FUNCTION RM_GetSolutionVolume 

SUBROUTINE Chk_GetSolutionVolume(id, vol)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: vol
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, vol, rmf_nxyz, "vol", "RM_GetSolutionVolume")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetSolutionVolume")
    endif
END SUBROUTINE Chk_GetSolutionVolume

!> Transfer concentrations of aqueous species to the array argument (@a species_conc)
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> The list of aqueous
!> species is determined by @ref RM_FindComponents and includes all
!> aqueous species that can be made from the set of components.
!> Solution volumes used to calculate mol/L are calculated by the reaction module.
!> Only the following databases distributed with PhreeqcRM have molar volume information 
!> needed to accurately calculate solution volume:
!> phreeqc.dat, Amm.dat, and pitzer.dat.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param species_conc     Array to receive the aqueous species concentrations. 
!> Dimension of the array is (@a nxyz, @a nspecies),
!> where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount), 
!> and @a nspecies is the number of aqueous species (@ref RM_GetSpeciesCount).
!> Concentrations are moles per liter.
!> Values for inactive cells are set to 1e30.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesCount, @ref RM_GetSpeciesD25, @ref RM_GetSpeciesZ,
!> @ref RM_GetSpeciesName, @ref RM_SpeciesConcentrations2Module, @ref RM_GetSpeciesSaveOn, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetSpeciesSaveOn(id, 1)
!> ncomps = RM_FindComponents(id)
!> nspecies = RM_GetSpeciesCount(id)
!> nxyz = RM_GetGridCellCount(id)
!> allocate(species_c(nxyz, nspecies))
!> status = RM_RunCells(id)
!> status = RM_GetSpeciesConcentrations(id, species_c)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_GetSpeciesConcentrations(id, species_conc) 
	USE ISO_C_BINDING  
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSpeciesConcentrations(id, species_conc) &
			BIND(C, NAME='RMF_GetSpeciesConcentrations')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: species_conc(*)
        END FUNCTION RMF_GetSpeciesConcentrations 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:,:) :: species_conc
	if (rmf_debug) call Chk_GetSpeciesConcentrations(id, species_conc)
    RM_GetSpeciesConcentrations = RMF_GetSpeciesConcentrations(id, species_conc)
END FUNCTION RM_GetSpeciesConcentrations 

SUBROUTINE Chk_GetSpeciesConcentrations(id, species_conc)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: species_conc
    INTEGER :: errors, nspecies
    nspecies = RM_GetSpeciesCount(id)
    errors = 0
    errors = errors + Chk_Double2D(id, species_conc, rmf_nxyz, nspecies, "species concentration", "RM_GetSpeciesConcentrations")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetSpeciesConcentrations")
    endif
END SUBROUTINE Chk_GetSpeciesConcentrations

!> The number of aqueous species used in the reaction module.
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> The list of aqueous
!> species is determined by @ref RM_FindComponents and includes all
!> aqueous species that can be made from the set of components.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT      The number of aqueous species, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, 
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesZ,
!> @ref RM_GetSpeciesName, @ref RM_SpeciesConcentrations2Module, 
!> @ref RM_GetSpeciesSaveOn, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetSpeciesSaveOn(id, 1)
!> ncomps = RM_FindComponents(id)
!> nspecies = RM_GetSpeciesCount(id)
!> nxyz = RM_GetGridCellCount(id)
!> allocate(species_c(nxyz, nspecies))
!> status = RM_RunCells(id)
!> status = RM_GetSpeciesConcentrations(id, species_c)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetSpeciesCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSpeciesCount(id) &
			BIND(C, NAME='RMF_GetSpeciesCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetSpeciesCount
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetSpeciesCount = RMF_GetSpeciesCount(id)
END FUNCTION RM_GetSpeciesCount

!> Transfers diffusion coefficients at 25C to the array argument (@a diffc).
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> Diffusion coefficients are defined in SOLUTION_SPECIES data blocks, normally in the database file.
!> Databases distributed with the reaction module that have diffusion coefficients defined are
!> phreeqc.dat, Amm.dat, and pitzer.dat.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param diffc            Array to receive the diffusion coefficients at 25 C, m^2/s.
!> Dimension of the array is @a nspecies,
!> where @a nspecies is is the number of aqueous species (@ref RM_GetSpeciesCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, @ref RM_GetSpeciesCount, @ref RM_GetSpeciesZ,
!> @ref RM_GetSpeciesName, @ref RM_SpeciesConcentrations2Module, @ref RM_GetSpeciesSaveOn, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetSpeciesSaveOn(id, 1)
!> ncomps = RM_FindComponents(id)
!> nspecies = RM_GetSpeciesCount(id)
!> allocate(diffc(nspecies))
!> status = RM_GetSpeciesD25(id, diffc)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetSpeciesD25(id, diffc)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSpeciesD25(id, diffc) &
			BIND(C, NAME='RMF_GetSpeciesD25')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: diffc(*)
        END FUNCTION RMF_GetSpeciesD25 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:) :: diffc
	if (rmf_debug) call Chk_GetSpeciesD25(id, diffc)
    RM_GetSpeciesD25 = RMF_GetSpeciesD25(id, diffc)
END FUNCTION RM_GetSpeciesD25 

SUBROUTINE Chk_GetSpeciesD25(id, diffc)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: diffc
    INTEGER :: errors, nspecies
    nspecies = RM_GetSpeciesCount(id)
    errors = 0
    errors = errors + Chk_Double1D(id, diffc, nspecies, "diffusion coefficient", "RM_GetSpeciesD25")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetSpeciesD25")
    endif
END SUBROUTINE Chk_GetSpeciesD25

!> Transfers the name of the @a ith aqueous species to the character argument (@a name).
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> The list of aqueous
!> species is determined by @ref RM_FindComponents and includes all
!> aqueous species that can be made from the set of components.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param i                Sequence number of the species in the species list. Fortran, 1 based.
!> @param name             Character array to receive the species name.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, @ref RM_GetSpeciesCount,
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesZ,
!> @ref RM_SpeciesConcentrations2Module, @ref RM_GetSpeciesSaveOn, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> char*100 name
!> ...
!> status = RM_SetSpeciesSaveOn(id, 1)
!> ncomps = RM_FindComponents(id)
!> nspecies = RM_GetSpeciesCount(id)
!> do i = 1, nspecies
!>   status = RM_GetSpeciesName(id, i, name)
!>   write(*,*) name
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetSpeciesName(id, i, name)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSpeciesName(id, i, name, l) &
			BIND(C, NAME='RMF_GetSpeciesName')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, i, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: name(*)
        END FUNCTION RMF_GetSpeciesName
	END INTERFACE
    INTEGER, INTENT(in) :: id, i
    CHARACTER(len=*), INTENT(out) :: name
    RM_GetSpeciesName = RMF_GetSpeciesName(id, i, name, len(name))
END FUNCTION RM_GetSpeciesName

!> Returns the value of the species-save property.
!> By default, concentrations of aqueous species are not saved. Setting the species-save property to true allows
!> aqueous species concentrations to be retrieved
!> with @ref RM_GetSpeciesConcentrations, and solution compositions to be set with
!> @ref RM_SpeciesConcentrations2Module.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT      0, species are not saved; 1, species are saved; negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, @ref RM_GetSpeciesCount,
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesZ,
!> @ref RM_GetSpeciesName, @ref RM_SpeciesConcentrations2Module, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> save_on = RM_GetSpeciesSaveOn(id)
!> if (save_on .ne. 0) then
!>   write(*,*) "Reaction module is saving species concentrations"
!> else
!>   write(*,*) "Reaction module is not saving species concentrations"
!> end
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetSpeciesSaveOn(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSpeciesSaveOn(id) &
			BIND(C, NAME='RMF_GetSpeciesSaveOn')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetSpeciesSaveOn
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetSpeciesSaveOn = RMF_GetSpeciesSaveOn(id)
END FUNCTION RM_GetSpeciesSaveOn

!> Transfers the charge of each aqueous species to the array argument (@a  z).
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param z                Array that receives the charge for each aqueous species.
!> Dimension of the array is @a nspecies,
!> where @a nspecies is is the number of aqueous species (@ref RM_GetSpeciesCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, @ref RM_GetSpeciesCount,
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesName, @ref RM_SpeciesConcentrations2Module, 
!> @ref RM_GetSpeciesSaveOn, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetSpeciesSaveOn(id, 1)
!> ncomps = RM_FindComponents(id)
!> nspecies = RM_GetSpeciesCount(id)
!> allocate(z(nspecies))
!> status = RM_GetSpeciesZ(id, z)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetSpeciesZ(id, z)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetSpeciesZ(id, z) &
			BIND(C, NAME='RMF_GetSpeciesZ')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(out) :: z(*)
        END FUNCTION RMF_GetSpeciesZ 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(out), DIMENSION(:) :: z
	if (rmf_debug) call Chk_GetSpeciesZ(id, z) 
    RM_GetSpeciesZ = RMF_GetSpeciesZ(id, z)
END FUNCTION RM_GetSpeciesZ 
 
SUBROUTINE Chk_GetSpeciesZ(id, z) 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: z
    INTEGER :: errors, nspecies
    nspecies = RM_GetSpeciesCount(id)
    errors = 0
    errors = errors + Chk_Double1D(id, z, nspecies, "species charge", "RM_GetSpeciesZ")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_GetSpeciesZ")
    endif
END SUBROUTINE Chk_GetSpeciesZ

!> Returns an array with the starting cell numbers from the range of cell numbers assigned to each worker.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param sc               Array to receive the starting cell numbers. Dimension of the array is 
!>                         the number of threads (OpenMP) or the number of processes (MPI).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_Create, @ref RM_GetEndCell, @ref RM_GetMpiTasks, @ref RM_GetThreadCount.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> n = RM_GetThreadCount(id) * RM_GetMpiTasks(id)
!> allocate(sc(n))
!> status = RM_GetStartCell(id, sc)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_GetStartCell(id, sc)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetStartCell(id, sc) &
			BIND(C, NAME='RMF_GetStartCell')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(out):: sc(*)
        END FUNCTION RMF_GetStartCell 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(out), DIMENSION(:) :: sc
    if (rmf_debug) call Chk_GetStartCell(id, sc)
    RM_GetStartCell = RMF_GetStartCell(id, sc)
    RETURN
END FUNCTION RM_GetStartCell

SUBROUTINE Chk_GetStartCell(id, sc)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:) :: sc
    INTEGER :: errors, n
    errors = 0
    n = RM_GetMpiTasks(id) * RM_GetThreadCount(id)
    errors = errors + Chk_Integer1D(id, sc, n, "StartCell", "RM_GetStartCell")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RMF_GetStartCell")
    endif
END SUBROUTINE Chk_GetStartCell

!> Returns the number of threads, which is equal to the number of workers used to run in parallel with OPENMP.
!> For the OPENMP version, the number of threads is set implicitly or explicitly with @ref RM_Create. For the
!> MPI version, the number of threads is always one for each process.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 The number of threads, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetMpiTasks.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string1, "(A,I)") "Number of threads: ", RM_GetThreadCount(id)
!> status = RM_OutputMessage(id, string1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers; result is always 1.

INTEGER FUNCTION RM_GetThreadCount(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_GetThreadCount(id) &
			BIND(C, NAME='RMF_GetThreadCount')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetThreadCount
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetThreadCount = RMF_GetThreadCount(id)
END FUNCTION RM_GetThreadCount

!> Returns the current simulation time in seconds. The reaction module does not change the time value, so the
!> returned value is equal to the default (0.0) or the last time set by @ref RM_SetTime.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 The current simulation time in seconds.
!> @see                    @ref RM_GetTimeConversion, @ref RM_GetTimeStep, @ref RM_SetTime,
!> @ref RM_SetTimeConversion, @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
!>       RM_GetTime(id) * RM_GetTimeConversion(id), " days"
!> status = RM_LogMessage(id, string)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.
        
DOUBLE PRECISION FUNCTION RM_GetTime(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_GetTime(id) &
			BIND(C, NAME='RMF_GetTime')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetTime
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetTime = RMF_GetTime(id)
END FUNCTION RM_GetTime
  
!> Returns a multiplier to convert time from seconds to another unit, as specified by the user.
!> The reaction module uses seconds as the time unit. The user can set a conversion
!> factor (@ref RM_SetTimeConversion) and retrieve it with RM_GetTimeConversion. The
!> reaction module only uses the conversion factor when printing the long version
!> of cell chemistry (@ref RM_SetPrintChemistryOn), which is rare. Default conversion factor is 1.0.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 Multiplier to convert seconds to another time unit.
!> @see                    @ref RM_GetTime, @ref RM_GetTimeStep, @ref RM_SetTime, @ref RM_SetTimeConversion, @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
!>       RM_GetTime(id) * RM_GetTimeConversion(id), " days"
!> status = RM_LogMessage(id, string)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

DOUBLE PRECISION FUNCTION RM_GetTimeConversion(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_GetTimeConversion(id) &
			BIND(C, NAME='RMF_GetTimeConversion')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetTimeConversion
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetTimeConversion = RMF_GetTimeConversion(id)
END FUNCTION RM_GetTimeConversion

!> Returns the current simulation time step in seconds.
!> This is the time over which kinetic reactions are integrated in a call to @ref RM_RunCells.
!> The reaction module does not change the time step value, so the
!> returned value is equal to the default (0.0) or the last time step set by @ref RM_SetTimeStep.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval                 The current simulation time step in seconds.
!> @see                    @ref RM_GetTime, @ref RM_GetTimeConversion, @ref RM_SetTime,
!> @ref RM_SetTimeConversion, @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string, "(A32,F15.1,A)") "          Time step             ", &
!>       RM_GetTimeStep(id) * RM_GetTimeConversion(id), " days"
!> status = RM_LogMessage(id, string)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

DOUBLE PRECISION FUNCTION RM_GetTimeStep(id)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION RMF_GetTimeStep(id) &
			BIND(C, NAME='RMF_GetTimeStep')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_GetTimeStep 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_GetTimeStep = RMF_GetTimeStep(id)
END FUNCTION RM_GetTimeStep 

!> Fills an array (@a bc_conc) with concentrations from solutions in the InitialPhreeqc instance.
!> The method is used to obtain concentrations for boundary conditions. If a negative value
!> is used for a cell in @a bc1, then the highest numbered solution in the InitialPhreeqc instance
!> will be used for that cell. Concentrations may be a mixture of two
!> solutions, @a bc1 and @a bc2, with a mixing fraction for @a bc1 1 of
!> @a f1 and mixing fraction for @a bc2 of (1 - @a f1).
!> A negative value for @a bc2 implies no mixing, and the associated value for @a f1 is ignored.
!> If @a bc2 and @a f1 are omitted,
!> no mixing is used; concentrations are derived from @a bc1 only.
!> 
!> @param id                  The instance @a id returned from @ref RM_Create.
!> @param bc_conc                   Array of concentrations extracted from the InitialPhreeqc instance.
!> The dimension of @a bc_conc is (@a n_boundary, @a ncomp),
!> where @a ncomp is the number of components returned from @ref RM_FindComponents or @ref RM_GetComponentCount.
!> @param n_boundary          The number of boundary condition solutions that need to be filled.
!> @param bc1  Array of solution index numbers that refer to solutions in the InitialPhreeqc instance.
!> Size is @a n_boundary.
!> @param bc2  Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance
!> and are defined to mix with @a bc1.
!> Size is @a n_boundary. Optional in Fortran.
!> @param f1           Fraction of @a bc1 that mixes with (1-@a f1) of @a bc2.
!> Size is (n_boundary). Optional in Fortran.
!> @retval IRM_RESULT         0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                       @ref RM_FindComponents, @ref RM_GetComponentCount.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> nbound = 1
!> allocate(bc1(nbound), bc2(nbound), f1(nbound))
!> allocate(bc_conc(nbound, ncomps))
!> bc1 = 0           ! solution 0 from InitialPhreeqc instance
!> bc2 = -1          ! no bc2 solution for mixing
!> f1 = 1.0          ! mixing fraction for bc1
!> status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, f1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_InitialPhreeqc2Concentrations(id, bc_conc, n_boundary, bc1, bc2, f1) 
	USE ISO_C_BINDING  
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqc2Concentrations(id, bc_conc, n_boundary, bc1) &
			BIND(C, NAME='RMF_InitialPhreeqc2Concentrations')
			USE ISO_C_BINDING   
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(OUT) :: bc_conc(*)
            INTEGER(KIND=C_INT), INTENT(IN) :: n_boundary, bc1(*)
        END FUNCTION RMF_InitialPhreeqc2Concentrations    
        INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqc2Concentrations2(id, bc_conc, n_boundary, bc1, bc2, f1) &
			BIND(C, NAME='RMF_InitialPhreeqc2Concentrations2')
			USE ISO_C_BINDING   
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(OUT) :: bc_conc(*)
            INTEGER(KIND=C_INT), INTENT(IN) :: n_boundary, bc1(*)
            INTEGER(KIND=C_INT), INTENT(IN) :: bc2(*)
            REAL(KIND=C_DOUBLE), INTENT(IN) :: f1(*)
        END FUNCTION RMF_InitialPhreeqc2Concentrations2
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(:,:) :: bc_conc
    INTEGER, INTENT(IN) :: n_boundary 
    INTEGER, INTENT(IN), DIMENSION(:) :: bc1
    INTEGER, INTENT(IN), DIMENSION(:) , OPTIONAL :: bc2
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) , OPTIONAL :: f1
	if (rmf_debug) call Chk_InitialPhreeqc2Concentrations(id, bc_conc, n_boundary, bc1, bc2, f1) 
    if (present(bc2) .and. present(f1)) then
        RM_InitialPhreeqc2Concentrations = RMF_InitialPhreeqc2Concentrations2(id, bc_conc, n_boundary, bc1, bc2, f1)
    else
        RM_InitialPhreeqc2Concentrations = RMF_InitialPhreeqc2Concentrations(id, bc_conc, n_boundary, bc1)
    endif
END FUNCTION RM_InitialPhreeqc2Concentrations    

SUBROUTINE Chk_InitialPhreeqc2Concentrations(id, bc_conc, n_boundary, bc1, bc2, f1) 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: bc_conc
    INTEGER, INTENT(IN) :: n_boundary 
    INTEGER, INTENT(IN), DIMENSION(:) :: bc1
    INTEGER, INTENT(IN), DIMENSION(:) , OPTIONAL :: bc2
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) , OPTIONAL :: f1
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double2D(id, bc_conc, n_boundary, rmf_ncomps, "concentration", "RM_InitialPhreeqc2Concentrations")
    errors = errors + Chk_Integer1D(id, bc1, n_boundary, "bc1", "RM_InitialPhreeqc2Concentrations")
    if (present(bc2)) then
        errors = errors + Chk_Integer1D(id, bc2, n_boundary, "bc2", "RM_InitialPhreeqc2Concentrations")
    endif
    if (present(f1)) then
        errors = errors + Chk_Double1D(id, f1, n_boundary, "f1", "RM_InitialPhreeqc2Concentrations")
    endif
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_InitialPhreeqc2Concentrations")
    endif
END SUBROUTINE Chk_InitialPhreeqc2Concentrations

!> Transfer solutions and reactants from the InitialPhreeqc instance to the reaction-module workers, possibly with mixing.
!> In its simplest form, @a ic1 is used to select initial conditions, including solutions and reactants,
!> for each cell of the model, without mixing.
!> @a ic1 is dimensioned (@a nxyz, 7), where @a nxyz is the number of grid cells in the user's model
!> (@ref RM_GetGridCellCount). The dimension of 7 refers to solutions and reactants in the following order:
!> (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE, (4) SURFACE, (5) GAS_PHASE,
!> (6) SOLID_SOLUTIONS, and (7) KINETICS. In Fortran, ic1(100, 4) = 2, indicates that
!> cell 99 (0 based) contains the SURFACE definition with user number 2 that has been defined in the
!> InitialPhreeqc instance (either by @ref RM_RunFile or @ref RM_RunString).
!> @n@n
!> It is also possible to mix solutions and reactants to obtain the initial conditions for cells. For mixing,
!> @a ic2 contains numbers for a second entity that mixes with the entity defined in @a ic1.
!> @a F1 contains the mixing fraction for @a ic1, whereas (1 - @a f1) is the mixing fraction for
!> @a ic2.
!> In Fortran, ic1(100, 4) = 2, initial_conditions2(100, 4) = 3, f1(100, 4) = 0.25 indicates that
!> cell 99 (0 based) contains a mixture of 0.25 SURFACE 2 and 0.75 SURFACE 3, where the surface 
!> compositions have been defined in the InitialPhreeqc instance. 
!> If the user number in @a ic2 is negative, no mixing occurs.
!> If @a ic2 and @a f1 are omitted,
!> no mixing is used, and initial conditions are derived solely from @a ic1.
!> 
!> @param id                  The instance @a id returned from @ref RM_Create.
!> @param ic1 Array of solution and reactant index numbers that refer to definitions in the InitialPhreeqc instance.
!> Size is (@a nxyz,7). The order of definitions is given above.
!> Negative values are ignored, resulting in no definition of that entity for that cell.
!> @param ic2  Array of solution and reactant index numbers that refer to definitions in the InitialPhreeqc instance.
!> Nonnegative values of @a ic2 result in mixing with the entities defined in @a ic1.
!> Negative values result in no mixing.
!> Size is (@a nxyz,7). The order of definitions is given above.
!> Optional in Fortran; omitting results in no mixing.
!> @param f1           Fraction of ic1 that mixes with (1-@a f1) of ic2.
!> Size is (nxyz,7). The order of definitions is given above.
!> Optional in Fortran; omitting results in no mixing.
!> @retval IRM_RESULT          0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                        @ref RM_InitialPhreeqcCell2Module.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(ic1(nxyz,7), ic2(nxyz,7), f1(nxyz,7))
!> ic1 = -1
!> ic2 = -1
!> f1 = 1.0
!> do i = 1, nxyz
!>   ic1(i,1) = 1       ! Solution 1
!>   ic1(i,2) = -1      ! Equilibrium phases none
!>   ic1(i,3) = 1       ! Exchange 1
!>   ic1(i,4) = -1      ! Surface none
!>   ic1(i,5) = -1      ! Gas phase none
!>   ic1(i,6) = -1      ! Solid solutions none
!>   ic1(i,7) = -1      ! Kinetics none
!> enddo
!> status = RM_InitialPhreeqc2Module(id, ic1, ic2, f1)1))
!> ! No mixing is defined, so the following is equivalent
!> status = RM_InitialPhreeqc2Module(id, ic1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_InitialPhreeqc2Module(id, ic1, ic2, f1)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqc2Module(id, ic1) &
			BIND(C, NAME='RMF_InitialPhreeqc2Module')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: ic1(*)
        END FUNCTION RMF_InitialPhreeqc2Module  
	END INTERFACE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqc2Module2(id, ic1, ic2, f1) &
			BIND(C, NAME='RMF_InitialPhreeqc2Module2')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: ic1(*)
            INTEGER(KIND=C_INT), INTENT(in) :: ic2(*)
            REAL(KIND=C_DOUBLE), INTENT(in) :: f1(*)
        END FUNCTION RMF_InitialPhreeqc2Module2  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:,:) :: ic1
    INTEGER, INTENT(in), DIMENSION(:,:), OPTIONAL :: ic2
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:), OPTIONAL :: f1
	if (rmf_debug) call Chk_InitialPhreeqc2Module(id, ic1, ic2, f1)
    if (present(ic2) .and. present(f1)) then
        RM_InitialPhreeqc2Module = RMF_InitialPhreeqc2Module2(id, ic1, ic2, f1)
    else
        RM_InitialPhreeqc2Module = RMF_InitialPhreeqc2Module(id, ic1)  
    endif    
END FUNCTION RM_InitialPhreeqc2Module    

SUBROUTINE Chk_InitialPhreeqc2Module(id, ic1, ic2, f1) 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(IN), DIMENSION(:,:) :: ic1
    INTEGER, INTENT(IN), DIMENSION(:,:) , OPTIONAL :: ic2
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) , OPTIONAL :: f1
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Integer2D(id, ic1, rmf_nxyz, 7, "ic1", "RM_InitialPhreeqc2Module")
    if (present(ic2)) then
        errors = errors + Chk_Integer2D(id, ic2, rmf_nxyz, 7, "ic2", "RM_InitialPhreeqc2Module")
    endif
    if (present(f1)) then
        errors = errors + Chk_Double2D(id, f1, rmf_nxyz, 7, "f1", "RM_InitialPhreeqc2Module")
    endif
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_InitialPhreeqc2Module")
    endif
END SUBROUTINE Chk_InitialPhreeqc2Module

!> Fills an array (@a bc_conc) with aqueous species concentrations from solutions in the InitialPhreeqc instance.
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> The method is used to obtain aqueous species concentrations for boundary conditions. If a negative value
!> is used for a cell in @a bc1, then the highest numbered solution in the InitialPhreeqc instance
!> will be used for that cell.
!> Concentrations may be a mixture of two
!> solutions, @a bc1 and @a bc2, with a mixing fraction for @a bc1 of
!> @a f1 and mixing fraction for @a bc2 of (1 - @a f1).
!> A negative value for @a bc2 implies no mixing, and the associated value for @a f1 is ignored.
!> If @a bc2 and @a f1 are omitted,
!> no mixing is used; concentrations are derived from @a bc1 only.
!> 
!> @param id                  The instance @a id returned from @ref RM_Create.
!> @param bc_conc           Array of aqueous concentrations extracted from the InitialPhreeqc instance.
!> The dimension of @a species_c is (@a n_boundary, @a nspecies),
!> where @a nspecies is the number of aqueous species returned from @ref RM_GetSpeciesCount.
!> @param n_boundary          The number of boundary condition solutions that need to be filled.
!> @param bc1  Array of solution index numbers that refer to solutions in the InitialPhreeqc instance.
!> Size is @a n_boundary.
!> @param bc2  Array of solution index numbers that that refer to solutions in the InitialPhreeqc instance
!> and are defined to mix with @a bc1.
!> Size is @a n_boundary. Optional in Fortran.
!> @param f1           Fraction of @a bc1 that mixes with (1-@a f1) of @a bc2.
!> Size is @a n_boundary. Optional in Fortran.
!> @retval IRM_RESULT         0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                  @ref RM_FindComponents, @ref RM_GetSpeciesCount, @ref RM_SetSpeciesSaveOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> nbound = 1
!> allocate(bc1(nbound), bc2(nbound), f1(nbound))
!> allocate(bc_conc(nbound, ncomps))
!> bc1 = 0           ! solution 0 from InitialPhreeqc instance
!> bc2 = -1          ! no bc2 solution for mixing
!> f1 = 1.0          ! mixing fraction for bc1
!> status = RM_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, nbound, bc1, bc2, f1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, n_boundary, bc1, bc2, f1) 
	USE ISO_C_BINDING  
    IMPLICIT NONE
    INTERFACE
           INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, n_boundary, bc1) &
			BIND(C, NAME='RMF_InitialPhreeqc2SpeciesConcentrations')   
			USE ISO_C_BINDING
                IMPLICIT NONE
                INTEGER(KIND=C_INT), INTENT(in) :: id
                REAL(KIND=C_DOUBLE), INTENT(OUT) :: bc_conc(*)
                INTEGER(KIND=C_INT), INTENT(IN) :: n_boundary, bc1(*)
        END FUNCTION RMF_InitialPhreeqc2SpeciesConcentrations    
        INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqc2SpeciesConcentrations2(id, bc_conc, n_boundary, bc1, bc2, f1) &
			BIND(C, NAME='RMF_InitialPhreeqc2SpeciesConcentrations2')   
			USE ISO_C_BINDING
                IMPLICIT NONE
                INTEGER(KIND=C_INT), INTENT(in) :: id
                REAL(KIND=C_DOUBLE), INTENT(OUT) :: bc_conc(*)
                INTEGER(KIND=C_INT), INTENT(IN) :: n_boundary, bc1(*)
                INTEGER(KIND=C_INT), INTENT(IN) :: bc2(*)
                REAL(KIND=C_DOUBLE), INTENT(IN) :: f1(*)
        END FUNCTION RMF_InitialPhreeqc2SpeciesConcentrations2  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: bc_conc
    INTEGER, INTENT(IN) :: n_boundary
    INTEGER, INTENT(IN), DIMENSION(:) :: bc1
    INTEGER, INTENT(IN), DIMENSION(:), OPTIONAL :: bc2
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:), OPTIONAL :: f1
	if (rmf_debug) call Chk_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, n_boundary, bc1, bc2, f1) 
    if (present(bc2) .and. present(f1)) then
        RM_InitialPhreeqc2SpeciesConcentrations = &
            RMF_InitialPhreeqc2SpeciesConcentrations2(id, bc_conc, n_boundary, bc1, bc2, f1)
    else
        RM_InitialPhreeqc2SpeciesConcentrations = &
            RMF_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, n_boundary, bc1)
    endif
END FUNCTION RM_InitialPhreeqc2SpeciesConcentrations          

SUBROUTINE Chk_InitialPhreeqc2SpeciesConcentrations(id, bc_conc, n_boundary, bc1, bc2, f1) 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: bc_conc
    INTEGER, INTENT(IN) :: n_boundary 
    INTEGER, INTENT(IN), DIMENSION(:) :: bc1
    INTEGER, INTENT(IN), DIMENSION(:) , OPTIONAL :: bc2
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) , OPTIONAL :: f1
    INTEGER :: errors, nspecies
    nspecies = RM_GetSpeciesCount(id)
    errors = 0
    errors = errors + Chk_Double2D(id, bc_conc, n_boundary, nspecies, "concentration", "RM_InitialPhreeqc2SpeciesConcentrations")
    errors = errors + Chk_Integer1D(id, bc1, n_boundary, "bc1", "RM_InitialPhreeqc2SpeciesConcentrations")
    if (present(bc2)) then
        errors = errors + Chk_Integer1D(id, bc2, n_boundary, "bc2", "RM_InitialPhreeqc2SpeciesConcentrations")
    endif
    if (present(f1)) then
        errors = errors + Chk_Double1D(id, f1, n_boundary, "f1", "RM_InitialPhreeqc2SpeciesConcentrations")
    endif
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_InitialPhreeqc2SpeciesConcentrations")
    endif
END SUBROUTINE Chk_InitialPhreeqc2SpeciesConcentrations

!> A cell numbered @a n_user in the InitialPhreeqc instance is selected to populate a series of cells.
!> All reactants with the number @a n_user are transferred along with the solution.
!> If MIX @a n_user exists, it is used for the definition of the solution.
!> If @a n_user is negative, @a n_user is redefined to be the largest solution or MIX number in the InitialPhreeqc instance.
!> All reactants for each cell in the list @a cell_numbers are removed before the cell
!> definition is copied from the InitialPhreeqc instance to the workers.
!> @param id                 The instance @a id returned from @ref RM_Create.
!> @param n_user                  Cell number refers to a solution or MIX and associated reactants in the InitialPhreeqc instance.
!> A negative number indicates the largest solution or MIX number in the InitialPhreeqc instance will be used.
!> @param cell_numbers     A list of cell numbers in the user's grid-cell numbering system that will be populated with
!> cell @a n_user from the InitialPhreeqc instance.
!> @param n_cell The number of cell numbers in the @a cell_numbers list.
!> @retval IRM_RESULT        0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                      @ref RM_InitialPhreeqc2Module.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate (cell_numbers(2))
!> cell_numbers(1) = 18
!> cell_numbers(2) = 19
!> ! n will be the largest SOLUTION number in InitialPhreeqc instance
!> ! copies solution and reactants to cells 18 and 19
!> status = RM_InitialPhreeqcCell2Module(id, -1, cell_numbers, 2)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_InitialPhreeqcCell2Module(id, n_user, cell_numbers, n_cell)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_InitialPhreeqcCell2Module(id, n_user, cell_numbers, n_cell) &
			BIND(C, NAME='RMF_InitialPhreeqcCell2Module')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: n_user
            INTEGER(KIND=C_INT), INTENT(in) :: cell_numbers(*)
            INTEGER(KIND=C_INT), INTENT(in) :: n_cell
        END FUNCTION RMF_InitialPhreeqcCell2Module  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: n_user
    INTEGER, INTENT(in), DIMENSION(:) :: cell_numbers
    INTEGER, INTENT(in) :: n_cell
	if (rmf_debug) call Chk_InitialPhreeqcCell2Module(id, n_user, cell_numbers, n_cell)
    RM_InitialPhreeqcCell2Module = RMF_InitialPhreeqcCell2Module(id, n_user, cell_numbers, n_cell)
END FUNCTION RM_InitialPhreeqcCell2Module   

SUBROUTINE Chk_InitialPhreeqcCell2Module(id, n_user, cell_numbers, n_cell)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id, n_user, n_cell
    INTEGER, INTENT(in), DIMENSION(:) :: cell_numbers
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Integer1D(id, cell_numbers, n_cell, "cell numbers", "RM_InitialPhreeqcCell2Module")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_InitialPhreeqcCell2Module")
    endif
END SUBROUTINE Chk_InitialPhreeqcCell2Module

!> Load a database for all IPhreeqc instances--workers, InitialPhreeqc, and Utility. All definitions
!> of the reaction module are cleared (SOLUTION_SPECIES, PHASES, SOLUTIONs, etc.), and the database is read.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param db_name          String containing the database name.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_Create.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_LoadDatabase(id, "phreeqc.dat")
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_LoadDatabase(id, db_name) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_LoadDatabase(id, db_name) &
			BIND(C, NAME='RMF_LoadDatabase') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: db_name(*)
        END FUNCTION RMF_LoadDatabase 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: db_name
    RM_LoadDatabase = RMF_LoadDatabase(id, trim(db_name)//C_NULL_CHAR)
END FUNCTION RM_LoadDatabase 
    
!> Print a message to the log file.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param str              String to be printed.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_OpenFiles, @ref RM_ErrorMessage, @ref RM_OutputMessage, @ref RM_ScreenMessage, @ref RM_WarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string, "(A32,F15.1,A)") "Beginning transport calculation ", &
!>       RM_GetTime(id) * RM_GetTimeConversion(id), " days"
!> status = RM_LogMessage(id, string)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_LogMessage(id, str) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_LogMessage(id, str) &
			BIND(C, NAME='RMF_LogMessage') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: str(*)
        END FUNCTION RMF_LogMessage
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: str
    RM_LogMessage = RMF_LogMessage(id, trim(str)//C_NULL_CHAR)
END FUNCTION RM_LogMessage

!> MPI only. Workers (processes with @ref RM_GetMpiMyself > 0) must call RM_MpiWorker to be able to
!> respond to messages from the root to accept data, perform calculations, and
!> (or) return data. RM_MpiWorker contains a loop that reads a message from root, performs a
!> task, and waits for another message from root. @ref RM_SetConcentrations, @ref RM_RunCells, and @ref RM_GetConcentrations
!> are examples of methods that send a message from root to get the workers to perform a task. The workers will
!> respond to all methods that are designated "workers must be in the loop of RM_MpiWorker" in the
!> MPI section of the method documentation.
!> The workers will continue to respond to messages from root until root calls
!> @ref RM_MpiWorkerBreak.
!> @n@n
!> (Advanced) The list of tasks that the workers perform can be extended by using @ref RM_SetMpiWorkerCallback.
!> It is then possible to use the MPI processes to perform other developer-defined tasks, such as transport calculations, without
!> exiting from the RM_MpiWorker loop. Alternatively, root calls @ref RM_MpiWorkerBreak to allow the workers to continue
!> past a call to RM_MpiWorker. The workers perform developer-defined calculations, and then RM_MpiWorker is called again to respond to
!> requests from root to perform reaction-module tasks.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError). RM_MpiWorker returns a value only when
!> @ref RM_MpiWorkerBreak is called by root.
!> @see                    @ref RM_MpiWorkerBreak, @ref RM_SetMpiWorkerCallback.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_MpiWorker(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by all workers.

INTEGER FUNCTION RM_MpiWorker(id) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_MpiWorker(id) &
			BIND(C, NAME='RMF_MpiWorker') 
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_MpiWorker
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_MpiWorker = RMF_MpiWorker(id)
END FUNCTION RM_MpiWorker

!> MPI only. This method is called by root to force workers (processes with @ref RM_GetMpiMyself > 0)
!> to return from a call to @ref RM_MpiWorker.
!> @ref RM_MpiWorker contains a loop that reads a message from root, performs a
!> task, and waits for another message from root. The workers respond to all methods that are designated
!> "workers must be in the loop of RM_MpiWorker" in the
!> MPI section of the method documentation.
!> The workers will continue to respond to messages from root until root calls RM_MpiWorkerBreak.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_MpiWorker, @ref RM_SetMpiWorkerCallback.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_MpiWorkerBreak(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_MpiWorkerBreak(id) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_MpiWorkerBreak(id) &
			BIND(C, NAME='RMF_MpiWorkerBreak') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_MpiWorkerBreak
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_MpiWorkerBreak = RMF_MpiWorkerBreak(id)
END FUNCTION RM_MpiWorkerBreak

!> Opens the output and log files. Files are named prefix.chem.txt and prefix.log.txt 
!> based on the prefix defined by @ref RM_SetFilePrefix.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetFilePrefix, @ref RM_GetFilePrefix, @ref RM_CloseFiles,
!> @ref RM_ErrorMessage, @ref RM_LogMessage, @ref RM_OutputMessage, @ref RM_WarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetFilePrefix(id, "Advect_f90")
!> status = RM_OpenFiles(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_OpenFiles(id) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_OpenFiles(id) &
			BIND(C, NAME='RMF_OpenFiles') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_OpenFiles
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_OpenFiles = RMF_OpenFiles(id)
END FUNCTION RM_OpenFiles

!> Print a message to the output file.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param str              String to be printed.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_ErrorMessage, @ref RM_LogMessage, @ref RM_ScreenMessage, @ref RM_WarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string1, "(A,I10)") "Number of threads:                                ", RM_GetThreadCount(id)
!> status = RM_OutputMessage(id, string1)
!> write(string1, "(A,I10)") "Number of MPI processes:                          ", RM_GetMpiTasks(id)
!> status = RM_OutputMessage(id, string1)
!> write(string1, "(A,I10)") "MPI task number:                                  ", RM_GetMpiMyself(id)
!> status = RM_OutputMessage(id, string1)
!> status = RM_GetFilePrefix(id, string)
!> write(string1, "(A,A)") "File prefix:                                      ", string
!> status = RM_OutputMessage(id, trim(string1))
!> write(string1, "(A,I10)") "Number of grid cells in the user's model:         ", RM_GetGridCellCount(id)
!> status = RM_OutputMessage(id, trim(string1))
!> write(string1, "(A,I10)") "Number of chemistry cells in the reaction module: ", RM_GetChemistryCellCount(id)
!> status = RM_OutputMessage(id, trim(string1))
!> write(string1, "(A,I10)") "Number of components for transport:               ", RM_GetComponentCount(id)
!> status = RM_OutputMessage(id, trim(string1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_OutputMessage(id, str)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_OutputMessage(id, str) &
			BIND(C, NAME='RMF_OutputMessage')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: str(*)
        END FUNCTION RMF_OutputMessage
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: str
    RM_OutputMessage = RMF_OutputMessage(id, trim(str)//C_NULL_CHAR)
END FUNCTION RM_OutputMessage

!> Runs a reaction step for all of the cells in the reaction module.
!> Normally, tranport concentrations are transferred to the reaction cells (@ref RM_SetConcentrations) before
!> reaction calculations are run. The length of time over which kinetic reactions are integrated is set
!> by @ref RM_SetTimeStep. Other properties that may need to be updated as a result of the transport
!> calculations include porosity (@ref RM_SetPorosity), saturation (@ref RM_SetSaturation),
!> temperature (@ref RM_SetTemperature), and pressure (@ref RM_SetPressure).
!> @param id               The instance @a id returned from @ref RM_Create.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetConcentrations,  @ref RM_SetPorosity,
!> @ref RM_SetPressure, @ref RM_SetSaturation, @ref RM_SetTemperature, @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetPorosity(id, por)                ! If pore volume changes
!> status = RM_SetSaturation(id, sat)              ! If saturation changes
!> status = RM_SetTemperature(id, temperature)     ! If temperature changes
!> status = RM_SetPressure(id, pressure)           ! If pressure changes
!> status = RM_SetConcentrations(id, c)          ! Transported concentrations
!> status = RM_SetTimeStep(id, time_step)             ! Time step for kinetic reactions
!> status = RM_RunCells(id)
!> status = RM_GetConcentrations(id, c)          ! Concentrations after reaction
!> status = RM_GetDensity(id, density)             ! Density after reaction
!> status = RM_GetSolutionVolume(id, volume)       ! Solution volume after reaction
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_RunCells(id)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_RunCells(id) &
			BIND(C, NAME='RMF_RunCells')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RMF_RunCells  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    RM_RunCells = RMF_RunCells(id)
END FUNCTION RM_RunCells  

!> Run a PHREEQC input file. The first three arguments determine which IPhreeqc instances will run
!> the file--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
!> files that modify the thermodynamic database should be run by all three sets of instances.
!> Files with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
!> be run by the workers. Files that contain initial conditions or boundary conditions should
!> be run by the InitialPhreeqc instance.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param workers          1, the workers will run the file; 0, the workers will not run the file.
!> @param initial_phreeqc  1, the InitialPhreeqc instance will run the file; 0, the InitialPhreeqc will not run the file.
!> @param utility          1, the Utility instance will run the file; 0, the Utility instance will not run the file.
!> @param chem_name        Name of the file to run.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_RunString.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_RunFile(id, 1, 1, 1, "advect.pqi")
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_RunFile(id, workers, initial_phreeqc, utility, chem_name)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_RunFile(id, workers, initial_phreeqc, utility, chem_name) &
			BIND(C, NAME='RMF_RunFile')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: workers, initial_phreeqc, utility
            CHARACTER(KIND=C_CHAR), INTENT(in) :: chem_name(*)
        END FUNCTION RMF_RunFile   
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: workers, initial_phreeqc, utility
    CHARACTER(len=*), INTENT(in) :: chem_name
    RM_RunFile = RMF_RunFile(id, workers, initial_phreeqc, utility, trim(chem_name)//C_NULL_CHAR)
END FUNCTION RM_RunFile   

!> Run a PHREEQC input string. The first three arguments determine which
!> IPhreeqc instances will run
!> the string--the workers, the InitialPhreeqc instance, and (or) the Utility instance. Input
!> strings that modify the thermodynamic database should be run by all three sets of instances.
!> Strings with SELECTED_OUTPUT definitions that will be used during the time-stepping loop need to
!> be run by the workers. Strings that contain initial conditions or boundary conditions should
!> be run by the InitialPhreeqc instance.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param workers          1, the workers will run the string; 0, the workers will not run the string.
!> @param initial_phreeqc  1, the InitialPhreeqc instance will run the string; 0, the InitialPhreeqc will not run the string.
!> @param utility          1, the Utility instance will run the string; 0, the Utility instance will not run the string.
!> @param input_string     String containing PHREEQC input.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_RunFile.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> string = "DELETE; -all"
!> status = RM_RunString(id, 1, 0, 1, string)  ! workers, initial_phreeqc, utility
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_RunString(id, initial_phreeqc, workers, utility, input_string)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_RunString(id, initial_phreeqc, workers, utility, input_string) &
			BIND(C, NAME='RMF_RunString')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: initial_phreeqc, workers, utility
            CHARACTER(KIND=C_CHAR), INTENT(in) :: input_string(*)
        END FUNCTION RMF_RunString   
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: initial_phreeqc, workers, utility
    CHARACTER(len=*), INTENT(in) :: input_string
    RM_RunString = RMF_RunString(id, initial_phreeqc, workers, utility, trim(input_string)//C_NULL_CHAR)
END FUNCTION RM_RunString   

!> Print message to the screen.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param str              String to be printed.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_ErrorMessage,  @ref RM_LogMessage, @ref RM_OutputMessage, @ref RM_WarningMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> write(string, "(A32,F15.1,A)") "Beginning reaction calculation  ", &
!>       time * RM_GetTimeConversion(id), " days"
!> status = RM_ScreenMessage(id, string)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_ScreenMessage(id, str) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_ScreenMessage(id, str) &
			BIND(C, NAME='RMF_ScreenMessage') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: str(*)
        END FUNCTION RMF_ScreenMessage 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: str
    RM_ScreenMessage = RMF_ScreenMessage(id, trim(str)//C_NULL_CHAR) 
END FUNCTION RM_ScreenMessage   

!> Select whether to include H2O in the component list.
!> The concentrations of H and O must be known
!> accurately (8 to 10 significant digits) for the numerical method of
!> PHREEQC to produce accurate pH and pe values.
!> Because most of the H and O are in the water species,
!> it may be more robust (require less accuracy in transport) to
!> transport the excess H and O (the H and O not in water) and water.
!> The default setting (@a true) is to include water, excess H, and excess O as components.
!> A setting of @a false will include total H and total O as components.
!> @a RM_SetComponentH2O must be called before @ref RM_FindComponents.
!> 
!> @param id               The instance id returned from @ref RM_Create.
!> @param tf               0, total H and O are included in the component list; 1, excess H, excess O, and water
!> are included in the component list.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetComponentH2O(id, 0)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetComponentH2O(id, tf)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetComponentH2O(id, tf) &
			BIND(C, NAME='RMF_SetComponentH2O')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: tf
        END FUNCTION RMF_SetComponentH2O
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: tf
    RM_SetComponentH2O = RMF_SetComponentH2O(id, tf)
END FUNCTION RM_SetComponentH2O

!> Use the vector of concentrations (@a c) to set the moles of components in each reaction cell.
!> The volume of water in a cell is the product of porosity (@ref RM_SetPorosity), saturation (@ref RM_SetSaturation),
!> and reference volume (@ref RM_SetRepresentativeVolume).
!> The moles of each component are determined by the volume of water and per liter concentrations.
!> If concentration units (@ref RM_SetUnitsSolution) are mass fraction, the
!> density (as specified by @ref RM_SetDensity) is used to convert from mass fraction to per mass per liter.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param c                Array of component concentrations. Size of array is (@a nxyz, @a ncomps), where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount), and @a ncomps is the number of components as determined
!> by @ref RM_FindComponents or @ref RM_GetComponentCount.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetDensity, @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume,
!> @ref RM_SetSaturation, @ref RM_SetUnitsSolution.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(c(nxyz, ncomps))
!> ...
!> call advect_f90(c, bc_conc, ncomps, nxyz)
!> status = RM_SetPorosity(id, por)               ! If porosity changes
!> status = RM_SetSaturation(id, sat)             ! If saturation changes
!> status = RM_SetTemperature(id, temperature))   ! If temperature changes
!> status = RM_SetPressure(id, pressure)          ! If pressure changes
!> status = RM_SetConcentrations(id, c)           ! Transported concentrations
!> status = RM_SetTimeStep(id, time_step)         ! Time step for kinetic reactions
!> status = RM_SetTime(id, time)                  ! Current time
!> status = RM_RunCells(id)
!> status = RM_GetConcentrations(id, c)           ! Concentrations after reaction
!> status = RM_GetDensity(id, density)            ! Density after reaction
!> status = RM_GetSolutionVolume(id, volume)      ! Solution volume after reaction
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetConcentrations(id, c)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetConcentrations(id, c) &
			BIND(C, NAME='RMF_SetConcentrations')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: c(*)
        END FUNCTION RMF_SetConcentrations
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: c
	if (rmf_debug) call Chk_SetConcentrations(id, c)
    RM_SetConcentrations = RMF_SetConcentrations(id, c)
END FUNCTION RM_SetConcentrations
	
SUBROUTINE Chk_SetConcentrations(id, c)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: c
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double2D(id, c, rmf_nxyz, rmf_ncomps, "concentration", "RM_SetConcentrations")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetConcentrations")
    endif
END SUBROUTINE Chk_SetConcentrations

#ifdef SKIP
INTEGER FUNCTION RM_SetConcentrations1D(id, c)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetConcentrations(id, c) &
			BIND(C, NAME='RMF_SetConcentrations')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: c(*)
        END FUNCTION RMF_SetConcentrations
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: c
    RM_SetConcentrations1D = RMF_SetConcentrations(id, c)
END FUNCTION RM_SetConcentrations1D
#endif

!> Select the current selected output by user number. The user may define multiple SELECTED_OUTPUT
!> data blocks for the workers. A user number is specified for each data block. The value of
!> the argument @a n_user selects which of the SELECTED_OUTPUT definitions will be used
!> for selected-output operations.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param n_user           User number of the SELECTED_OUTPUT data block that is to be used.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetNthSelectedOutputUserNumber, @ref RM_GetSelectedOutput, @ref RM_GetSelectedOutputColumnCount,
!> @ref RM_GetSelectedOutputCount, @ref RM_GetSelectedOutputRowCount, @ref RM_GetSelectedOutputHeading,
!> @ref RM_SetSelectedOutputOn.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> do isel = 1, RM_GetSelectedOutputCount(id)
!>   n_user = RM_GetNthSelectedOutputUserNumber(id, isel)
!>   status = RM_SetCurrentSelectedOutputUserNumber(id, n_user)
!>   col = RM_GetSelectedOutputColumnCount(id)
!>   allocate(selected_out(nxyz,col))
!>   status = RM_GetSelectedOutput(id, selected_out)
!>   ! Process results here
!>   deallocate(selected_out)
!> enddo
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.
!>  */

INTEGER FUNCTION RM_SetCurrentSelectedOutputUserNumber(id, n_user)  
	USE ISO_C_BINDING 
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetCurrentSelectedOutputUserNumber(id, n_user) &
			BIND(C, NAME='RMF_SetCurrentSelectedOutputUserNumber') 
			USE ISO_C_BINDING  
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: n_user
        END FUNCTION RMF_SetCurrentSelectedOutputUserNumber
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: n_user
    RM_SetCurrentSelectedOutputUserNumber = RMF_SetCurrentSelectedOutputUserNumber(id, n_user)
END FUNCTION RM_SetCurrentSelectedOutputUserNumber

!> Set the density for each reaction cell. These density values are used 
!> when converting from transported mass fraction concentrations (@ref RM_SetUnitsSolution) to
!> produce per liter concentrations during a call to @ref RM_SetConcentrations.
!> They are also used when converting from module concentrations to transport concentrations
!> of mass fraction (@ref RM_GetConcentrations), if @ref RM_UseSolutionDensityVolume is set to @a false.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param density          Array of densities. Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetConcentrations, @ref RM_SetConcentrations, 
!> @ref RM_SetUnitsSolution, @ref RM_UseSolutionDensityVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(density(nxyz))
!> density = 1.0
!> status = RM_SetDensity(id, density(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetDensity(id, density)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetDensity(id, density) &
			BIND(C, NAME='RMF_SetDensity')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: density(*)
        END FUNCTION RMF_SetDensity 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: density
	if (rmf_debug) call Chk_SetDensity(id, density)
    RM_SetDensity = RMF_SetDensity(id, density)
END FUNCTION RM_SetDensity 

SUBROUTINE Chk_SetDensity(id, density)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: density
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, density, rmf_nxyz, "density", "RM_SetDensity")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetDensity")
    endif
END SUBROUTINE Chk_SetDensity

!> Set the name of the dump file. It is the name used by @ref RM_DumpModule.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param dump_name        Name of dump file.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_DumpModule.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetDumpFileName(id, "advection_f90.dmp")
!> dump_on = 1
!> append = 0
!> status = RM_DumpModule(id, dump_on, append)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_SetDumpFileName(id, dump_name) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetDumpFileName(id, dump_name) &
			BIND(C, NAME='RMF_SetDumpFileName') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: dump_name(*)
        END FUNCTION RMF_SetDumpFileName  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: dump_name
    RM_SetDumpFileName = RMF_SetDumpFileName(id, trim(dump_name)//C_NULL_CHAR)
END FUNCTION RM_SetDumpFileName   

!> Set the action to be taken when the reaction module encounters an error.
!> Options are 0, return to calling program with an error return code (default);
!> 1, throw an exception, in C++, the exception can be caught, for C and Fortran, the program will exit; or
!> 2, attempt to exit gracefully.
!> @param id               The instance id returned from @ref RM_Create.
!> @param mode             Error handling mode: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> id = RM_create(nxyz, nthreads)
!> status = RM_SetErrorHandlerMode(id, 2)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetErrorHandlerMode(id, mode)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetErrorHandlerMode(id, mode) &
			BIND(C, NAME='RMF_SetErrorHandlerMode')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: mode
        END FUNCTION RMF_SetErrorHandlerMode    
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: mode
    RM_SetErrorHandlerMode = RMF_SetErrorHandlerMode(id, mode)
END FUNCTION RM_SetErrorHandlerMode        

!> Set the prefix for the output (prefix.chem.txt) and log (prefix.log.txt) files.
!> These files are opened by @ref RM_OpenFiles.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param prefix           Prefix used when opening the output and log files.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_OpenFiles, @ref RM_CloseFiles.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetFilePrefix(id, "Advect_f90")
!> status = RM_OpenFiles(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_SetFilePrefix(id, prefix) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetFilePrefix(id, prefix) &
			BIND(C, NAME='RMF_SetFilePrefix') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: prefix(*)
        END FUNCTION RMF_SetFilePrefix  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: prefix
    RM_SetFilePrefix = RMF_SetFilePrefix(id, trim(prefix)//C_NULL_CHAR) 
END FUNCTION RM_SetFilePrefix  

!> MPI only. Defines a callback function that allows additional tasks to be done
!> by the workers. The method @ref RM_MpiWorker contains a loop,
!> where the workers receive a message (an integer),
!> run a function corresponding to that integer,
!> and then wait for another message.
!> RM_SetMpiWorkerCallback allows the developer to add another function
!> that responds to additional integer messages by calling developer-defined functions
!> corresponding to those integers.
!> @ref RM_MpiWorker calls the callback function when the message number
!> is not one of the PhreeqcRM message numbers.
!> Messages are unique integer numbers. PhreeqcRM uses integers in a range
!> beginning at 0. It is suggested that developers use message numbers starting
!> at 1000 or higher for their tasks.
!> The callback function calls a developer-defined function specified
!> by the message number and then returns to @ref RM_MpiWorker to wait for
!> another message.
!> @n@n
!> For Fortran, the functions that are called from the callback function
!> can use USE statements to find the data necessary to perform the tasks, and
!> the only argument to the callback function is an integer message argument.
!> @a RM_SetMpiWorkerCallback
!> must be called by each worker before @ref RM_MpiWorker is called.
!> @n@n
!> The motivation for this method is to allow the workers to perform other
!> tasks, for instance, parallel transport calculations, within the structure
!> of @ref RM_MpiWorker. The callback function
!> can be used to allow the workers to receive data, perform transport calculations,
!> and (or) send results, without leaving the loop of @ref RM_MpiWorker. Alternatively,
!> it is possible for the workers to return from @ref RM_MpiWorker
!> by a call to @ref RM_MpiWorkerBreak by root. The workers could then call
!> subroutines to receive data, calculate transport, and send data,
!> and then resume processing PhreeqcRM messages from root with another
!> call to @ref RM_MpiWorker.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param fcn              A function that returns an integer and has an integer argument.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_MpiWorker, @ref RM_MpiWorkerBreak.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> Code executed by root:
!> status = do_something()
!> 
!> Code executed by workers:
!> status = RM_SetMpiWorkerCallback(id, worker_tasks_f)
!> status = RM_MpiWorker(id)
!> 
!> Code executed by root and workers:    
!> integer function do_something
!>   implicit none
!>   INCLUDE 'mpif.h'
!>   integer status
!>   integer i, method_number, mpi_myself, mpi_task, mpi_tasks, worker_number;
!>   method_number = 1000
!>   call MPI_Comm_size(MPI_COMM_WORLD, mpi_tasks, status)
!>   call MPI_Comm_rank(MPI_COMM_WORLD, mpi_myself, status)
!>   if (mpi_myself .eq. 0) then     
!>     CALL MPI_Bcast(method_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status)
!>     write(*,*) "I am root."
!>     do i = 1, mpi_tasks-1
!>       CALL MPI_Recv(worker_number, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, status)
!>       write(*,*) "Recieved data from worker number ", worker_number, "."
!>     enddo
!>   else
!> 		CALL MPI_Send(mpi_myself, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status)
!>   endif
!>   do_something = 0
!> end function do_something
!> 
!> Code called by workers from method MpiWorker:
!> integer(kind=C_INT) function worker_tasks_f(method_number) BIND(C, NAME='worker_tasks_f')
!>   USE ISO_C_BINDING
!>   implicit none
!>   interface
!>     integer function do_something
!>     end function do_something
!>   end interface
!>   integer(kind=c_int), intent(in) :: method_number
!>   integer :: status
!>   if (method_number .eq. 1000) then
!>     status = do_something()
!>   endif
!>   worker_tasks_f = 0
!> end function worker_tasks_f
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by workers, before call to @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetMpiWorkerCallback(id, fcn)
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTERFACE
     INTEGER(KIND=C_INT) FUNCTION RMF_SetMpiWorkerCallback(id, fcn) &
          BIND(C, NAME='RMF_SetMpiWorkerCallback')
       USE ISO_C_BINDING
       INTEGER(KIND=C_INT), INTENT(in) :: id
!> \cond
       INTERFACE
          INTEGER(KIND=C_INT) FUNCTION fcn(method_number) BIND(C)
            USE ISO_C_BINDING
            INTEGER(KIND=C_INT), INTENT(in) :: method_number
          END FUNCTION fcn
       END INTERFACE
!> \endcond       
     END FUNCTION RMF_SetMpiWorkerCallback
  END INTERFACE
  INTEGER, INTENT(IN) :: id

  INTERFACE
     INTEGER(kind=c_int) FUNCTION fcn(method_number) BIND(C)
       USE ISO_C_BINDING
       !INTEGER, INTENT(in) :: method_number
       INTEGER(kind=c_int), INTENT(in) :: method_number
     END FUNCTION fcn
  END INTERFACE

  RM_SetMpiWorkerCallback = RMF_SetMpiWorkerCallback(id, fcn)
END FUNCTION RM_SetMpiWorkerCallback

!> Sets the property for partitioning solids between the saturated and unsaturated 
!> parts of a partially saturated cell. 
!> 
!> The option is intended to be used by saturated-only
!> flow codes that allow a variable water table.
!> The value has meaning only when saturations
!> less than 1.0 are encountered. The partially saturated cells
!> may have a small water-to-rock ratio that causes
!> reactions to proceed differently relative to fully saturated cells.
!> By setting  @a RM_SetPartitionUZSolids to true, the
!> amounts of solids and gases are partioned according to the saturation.
!> If a cell has a saturation of 0.5, then
!> the water interacts with only half of the solids and gases; the other half is unreactive
!> until the water table rises. As the saturation in a cell varies,
!> solids and gases are transferred between the
!> saturated and unsaturated (unreactive) reservoirs of the cell.
!> Unsaturated-zone flow and transport codes will probably use the default (false),
!> which assumes all gases and solids are reactive regardless of saturation.
!> 
!> @param id       The instance @a id returned from @ref RM_Create.
!> @param tf       @a True, the fraction of solids and gases available for 
!> reaction is equal to the saturation; 
!> @a False (default), all solids and gases are reactive regardless of saturation.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetPartitionUZSolids(id, 0)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetPartitionUZSolids(id, tf)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetPartitionUZSolids(id, tf) &
			BIND(C, NAME='RMF_SetPartitionUZSolids')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in)  :: tf
        END FUNCTION RMF_SetPartitionUZSolids 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in)  :: tf
    RM_SetPartitionUZSolids = RMF_SetPartitionUZSolids(id, tf)
END FUNCTION RM_SetPartitionUZSolids 

!> Set the porosity for each reaction cell. 
!> The volume of water in a reaction cell is the product of the porosity, the saturation
!> (@ref RM_SetSaturation), and the representative volume (@ref RM_SetRepresentativeVolume).
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param por              Array of porosities, unitless. Default is 0.1. Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetSaturation, @ref RM_SetRepresentativeVolume, @ref RM_SetSaturation.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(por(nxyz))
!> por = 0.2
!> status = RM_SetPorosity(id, por(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetPorosity(id, por)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetPorosity(id, por) &
			BIND(C, NAME='RMF_SetPorosity')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: por(*)
        END FUNCTION RMF_SetPorosity 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: por
	if (rmf_debug) call Chk_SetPorosity(id, por)
    RM_SetPorosity = RMF_SetPorosity(id, por)
END FUNCTION RM_SetPorosity 

SUBROUTINE Chk_SetPorosity(id, por)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: por
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, por, rmf_nxyz, "porosity", "RM_SetPorosity")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetPorosity")
    endif
END SUBROUTINE Chk_SetPorosity

!> Set the pressure for each reaction cell. Pressure effects are considered only in three of the
!> databases distributed with PhreeqcRM: phreeqc.dat, Amm.dat, and pitzer.dat.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param p                Array of pressures, in atm. Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetTemperature.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(pressure(nxyz))
!> pressure = 2.0
!> status = RM_SetPressure(id, pressure(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetPressure(id, p)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetPressure(id, p) &
			BIND(C, NAME='RMF_SetPressure')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: p(*)
        END FUNCTION RMF_SetPressure   
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: p
	if (rmf_debug) call Chk_SetPressure(id, p)
    RM_SetPressure = RMF_SetPressure(id, p)
END FUNCTION RM_SetPressure        

SUBROUTINE Chk_SetPressure(id, p) 
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: p
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, p, rmf_nxyz, "pressure", "RM_SetPressure")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetPressure")
    endif
END SUBROUTINE Chk_SetPressure

!> Enable or disable detailed output for each reaction cell. 
!> Printing for a cell will occur only when the
!> printing is enabled with @ref RM_SetPrintChemistryOn and the @a cell_mask value is 1.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param cell_mask        Array of integers. Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount). A value of 0 will
!> disable printing detailed output for the cell; a value of 1 will enable printing detailed output for a cell.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetPrintChemistryOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(print_chemistry_mask(nxyz))
!>   do i = 1, nxyz/2
!>   print_chemistry_mask(i) = 1
!>   print_chemistry_mask(i+nxyz/2) = 0
!> enddo
!> status = RM_SetPrintChemistryMask(id, print_chemistry_mask(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetPrintChemistryMask(id, cell_mask)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetPrintChemistryMask(id, cell_mask) &
			BIND(C, NAME='RMF_SetPrintChemistryMask') 
			USE ISO_C_BINDING  
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: cell_mask(*)
        END FUNCTION RMF_SetPrintChemistryMask 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, DIMENSION(:), INTENT(in) :: cell_mask
	if (rmf_debug) call Chk_SetPrintChemistryMask(id, cell_mask)
    RM_SetPrintChemistryMask = RMF_SetPrintChemistryMask(id, cell_mask)
END FUNCTION RM_SetPrintChemistryMask 

SUBROUTINE Chk_SetPrintChemistryMask(id, cell_mask)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:) :: cell_mask
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Integer1D(id, cell_mask, rmf_nxyz, "cell_mask", "RM_SetPrintChemistryMask")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetPrintChemistryMask")
    endif
END SUBROUTINE Chk_SetPrintChemistryMask

!> Setting to enable or disable printing detailed output from reaction calculations to the output file for a set of
!> cells defined by @ref RM_SetPrintChemistryMask. The detailed output prints all of the output
!> typical of a PHREEQC reaction calculation, which includes solution descriptions and the compositions of
!> all other reactants. The output can be several hundred lines per cell, which can lead to a very
!> large output file (prefix.chem.txt, @ref RM_OpenFiles). For the worker instances, the output can be limited to a set of cells
!> (@ref RM_SetPrintChemistryMask) and, in general, the
!> amount of information printed can be limited by use of options in the PRINT data block of PHREEQC 
!> (applied by using @ref RM_RunFile or @ref RM_RunString). 
!> Printing the detailed output for the workers is generally used only for debugging, and PhreeqcRM will run
!> significantly faster when printing detailed output for the workers is disabled.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param workers          0, disable detailed printing in the worker instances, 1, enable detailed printing
!> in the worker instances.
!> @param initial_phreeqc  0, disable detailed printing in the InitialPhreeqc instance, 1, enable detailed printing
!> in the InitialPhreeqc instances.
!> @param utility          0, disable detailed printing in the Utility instance, 1, enable detailed printing
!> in the Utility instance.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetPrintChemistryMask.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetPrintChemistryOn(id, 0, 1, 0)  ! workers, initial_phreeqc, utility
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetPrintChemistryOn(id, workers, initial_phreeqc, utility)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetPrintChemistryOn(id, workers, initial_phreeqc, utility) &
			BIND(C, NAME='RMF_SetPrintChemistryOn')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: workers, initial_phreeqc, utility
        END FUNCTION RMF_SetPrintChemistryOn 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: workers, initial_phreeqc, utility
    RM_SetPrintChemistryOn = RMF_SetPrintChemistryOn(id, workers, initial_phreeqc, utility)
END FUNCTION RM_SetPrintChemistryOn 

!> Set the load-balancing algorithm.
!> PhreeqcRM attempts to rebalance the load of each thread or process such that each
!> thread or process takes the same amount of time to run its part of a @ref RM_RunCells
!> calculation. Two algorithms are available; one uses individual times for each cell and
!> accounts for cells that were not run because
!> saturation was zero (default), and
!> the other assigns an average time to all cells.
!> The methods are similar, but limited testing indicates the default method performs better.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param method           0, indicates average times are used in rebalancing; 1 indicates individual
!> cell times are used in rebalancing (default).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetRebalanceFraction.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetRebalanceByCell(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetRebalanceByCell(id, method)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetRebalanceByCell(id, method) &
			BIND(C, NAME='RMF_SetRebalanceByCell')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in)  :: method
        END FUNCTION RMF_SetRebalanceByCell
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in)  :: method
    RM_SetRebalanceByCell = RMF_SetRebalanceByCell(id, method)
END FUNCTION RM_SetRebalanceByCell

!> Sets the fraction of cells that are transferred among threads or processes when rebalancing.
!> PhreeqcRM attempts to rebalance the load of each thread or process such that each
!> thread or process takes the same amount of time to run its part of a @ref RM_RunCells
!> calculation. The rebalancing transfers cell calculations among threads or processes to
!> try to achieve an optimum balance. @a RM_SetRebalanceFraction
!> adjusts the calculated optimum number of cell transfers by a fraction from 0 to 1.0 to
!> determine the actual number of cell transfers. A value of zero eliminates
!> load rebalancing. A value less than 1.0 is suggested to slow the approach to the optimum cell
!> distribution and avoid possible oscillations
!> when too many cells are transferred at one iteration, requiring reverse transfers at the next iteration.
!> Default is 0.5.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param f                Fraction from 0.0 to 1.0.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetRebalanceByCell.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetRebalanceFraction(id, 0.5d0)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetRebalanceFraction(id, f)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetRebalanceFraction(id, f) &
			BIND(C, NAME='RMF_SetRebalanceFraction')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in)  :: f
        END FUNCTION RMF_SetRebalanceFraction
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in)  :: f
    RM_SetRebalanceFraction = RMF_SetRebalanceFraction(id, f)
END FUNCTION RM_SetRebalanceFraction

!> Set the representative volume of each reaction cell.
!> By default the representative volume of each reaction cell is 1 liter.
!> The volume of water in a reaction cell is determined by the procuct of the representative volume,
!> the porosity (@ref RM_SetPorosity), and the saturation (@ref RM_SetSaturation).
!> The numerical method of PHREEQC is more robust if the water volume for a reaction cell is
!> within a couple orders of magnitude of 1.0.
!> Small water volumes caused by small porosities and (or) small saturations (and (or) small representative volumes)
!> may cause non-convergence of the numerical method.
!> In these cases, a larger representative volume may help. Note
!> that increasing the representative volume also increases
!> the number of moles of the reactants in the reaction cell (minerals, surfaces, exchangers,
!> and others), which are defined as moles per representative volume.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param rv              Vector of representative volumes, in liters. Default is 1.0 liter.
!> Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetPorosity, @ref RM_SetSaturation.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> double precision, dimension(:), allocatable   :: rv
!> allocate(rv(nxyz))
!> rv = 1.0
!> status = RM_SetRepresentativeVolume(id, rv(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetRepresentativeVolume(id, rv)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetRepresentativeVolume(id, rv) &
			BIND(C, NAME='RMF_SetRepresentativeVolume') 
			USE ISO_C_BINDING  
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: rv(*)
        END FUNCTION RMF_SetRepresentativeVolume 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: rv
    RM_SetRepresentativeVolume = RMF_SetRepresentativeVolume(id, rv)
END FUNCTION RM_SetRepresentativeVolume 

!> Set the saturation of each reaction cell. Saturation is a fraction ranging from 0 to 1.
!> The volume of water in a cell is the product of porosity (@ref RM_SetPorosity), saturation (@a RM_SetSaturation),
!> and representative volume (@ref RM_SetRepresentativeVolume). As a result of a reaction calculation,
!> solution properties (density and volume) will change;
!> the databases phreeqc.dat, Amm.dat, and pitzer.dat have the molar volume data to calculate these changes. 
!> The methods @ref RM_GetDensity, @ref RM_GetSolutionVolume, and @ref RM_GetSaturation 
!> can be used to account for these changes in the succeeding transport calculation.
!> @a RM_SetRepresentativeVolume should be called before initial conditions are defined for the reaction cells.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param sat              Array of saturations, unitless. Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetDensity, @ref RM_GetSaturation, @ref RM_GetSolutionVolume,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(sat(nxyz))
!> sat = 1.0
!> status = RM_SetSaturation(id, sat(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetSaturation(id, sat)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetSaturation(id, sat) &
			BIND(C, NAME='RMF_SetSaturation')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: sat(*)
        END FUNCTION RMF_SetSaturation 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: sat
	if (rmf_debug) call Chk_SetSaturation(id, sat)
    RM_SetSaturation = RMF_SetSaturation(id, sat)
END FUNCTION RM_SetSaturation 

SUBROUTINE Chk_SetSaturation(id, sat)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: sat
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, sat, rmf_nxyz, "sataturation", "RM_SetSaturation")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetSaturation")
    endif
END SUBROUTINE Chk_SetSaturation

!> Set the property that controls whether messages are written to the screen.
!> Messages include information about rebalancing during @ref RM_RunCells, and
!> any messages written with @ref RM_ScreenMessage.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param tf  @a 1, enable screen messages; @a 0, disable screen messages. Default is 1.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_RunCells, @ref RM_ScreenMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetScreenOn(rm_id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root.

INTEGER FUNCTION RM_SetScreenOn(id, tf)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetScreenOn(id, tf) &
			BIND(C, NAME='RMF_SetScreenOn')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: tf
        END FUNCTION RMF_SetScreenOn 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: tf
    RM_SetScreenOn = RMF_SetScreenOn(id, tf)
END FUNCTION RM_SetScreenOn 

!> Setting determines whether selected-output results are available to be retrieved
!> with @ref RM_GetSelectedOutput. @a 1 indicates that selected-output results
!> will be accumulated during @ref RM_RunCells and can be retrieved with @ref RM_GetSelectedOutput;
!> @a 0 indicates that selected-output results will not
!> be accumulated during @ref RM_RunCells.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param tf               0, disable selected output; 1, enable selected output.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_GetSelectedOutput, @ref RM_SetPrintChemistryOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetSelectedOutputOn(id, 1)        ! enable selected output
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetSelectedOutputOn(id, tf)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetSelectedOutputOn(id, tf) &
			BIND(C, NAME='RMF_SetSelectedOutputOn')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: tf
        END FUNCTION RMF_SetSelectedOutputOn  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: tf
    RM_SetSelectedOutputOn = RMF_SetSelectedOutputOn(id, tf)
END FUNCTION RM_SetSelectedOutputOn   

!> Sets the value of the species-save property.
!> This method enables use of PhreeqcRM with multicomponent-diffusion transport calculations.
!> By default, concentrations of aqueous species are not saved. Setting the species-save property to 1 allows
!> aqueous species concentrations to be retrieved
!> with @ref RM_GetSpeciesConcentrations, and solution compositions to be set with
!> @ref RM_SpeciesConcentrations2Module.
!> RM_SetSpeciesSaveOn must be called before calls to @ref RM_FindComponents.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param save_on          0, indicates species concentrations are not saved; 1, indicates species concentrations are
!> saved.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, @ref RM_GetSpeciesCount,
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesSaveOn, @ref RM_GetSpeciesZ,
!> @ref RM_GetSpeciesName, @ref RM_SpeciesConcentrations2Module.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> save_on = RM_SetSpeciesSaveOn(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers.

INTEGER FUNCTION RM_SetSpeciesSaveOn(id, save_on)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetSpeciesSaveOn(id, save_on) &
			BIND(C, NAME='RMF_SetSpeciesSaveOn')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: save_on
        END FUNCTION RMF_SetSpeciesSaveOn
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: save_on
    RM_SetSpeciesSaveOn = RMF_SetSpeciesSaveOn(id, save_on)
END FUNCTION RM_SetSpeciesSaveOn

!> Set the temperature for each reaction cell. If @a RM_SetTemperature is not called,
!> worker solutions will have temperatures as defined by initial conditions 
!> (@ref RM_InitialPhreeqc2Module and @ref RM_InitialPhreeqcCell2Module).
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param t                Array of temperatures, in degrees C. Size of array is @a nxyz, where @a nxyz is the number
!> of grid cells in the user's model (@ref RM_GetGridCellCount).
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_InitialPhreeqc2Module, 
!> @ref RM_InitialPhreeqcCell2Module, @ref RM_SetPressure.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> allocate(temperature(nxyz))
!> temperature = 20.0
!> status = RM_SetTemperature(id, temperature(1))
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetTemperature(id, t)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetTemperature(id, t) &
			BIND(C, NAME='RMF_SetTemperature')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: t(*)
        END FUNCTION RMF_SetTemperature 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: t
	if (rmf_debug) call Chk_SetTemperature(id, t)
    RM_SetTemperature = RMF_SetTemperature(id, t)
END FUNCTION RM_SetTemperature 

SUBROUTINE Chk_SetTemperature(id, t)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: t
    INTEGER :: errors
    errors = 0
    errors = errors + Chk_Double1D(id, t, rmf_nxyz, "temperature", "RM_SetTemperature")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SetTemperature")
    endif
END SUBROUTINE Chk_SetTemperature

!> Set current simulation time for the reaction module.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param time             Current simulation time, in seconds.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetTimeStep, @ref RM_SetTimeConversion.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetTime(id, time)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetTime(id, time)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetTime(id, time) &
			BIND(C, NAME='RMF_SetTime')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: time
        END FUNCTION RMF_SetTime 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in) :: time
    RM_SetTime = RMF_SetTime(id, time)
END FUNCTION RM_SetTime 

!> Set a factor to convert to user time units. Factor times seconds produces user time units.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param conv_factor      Factor to convert seconds to user time units.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetTime, @ref RM_SetTimeStep.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetTimeConversion(id, dble(1.0 / 86400.0)) ! days
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetTimeConversion(id, conv_factor)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetTimeConversion(id, conv_factor) &
			BIND(C, NAME='RMF_SetTimeConversion') 
			USE ISO_C_BINDING  
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: conv_factor
        END FUNCTION RMF_SetTimeConversion 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in) :: conv_factor
    RM_SetTimeConversion = RMF_SetTimeConversion(id, conv_factor)
END FUNCTION RM_SetTimeConversion 

!> Set current time step for the reaction module. This is the length
!> of time over which kinetic reactions are integrated.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param time_step        Current time step, in seconds.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetTime, @ref RM_SetTimeConversion.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetTimeStep(id, time_step)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetTimeStep(id, time_step)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetTimeStep(id, time_step) &
			BIND(C, NAME='RMF_SetTimeStep')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: time_step
        END FUNCTION RMF_SetTimeStep 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in) :: time_step
    RM_SetTimeStep = RMF_SetTimeStep(id, time_step)
END FUNCTION RM_SetTimeStep 

!> Sets input units for exchangers.
!> In PHREEQC input, exchangers are defined by moles of exchange sites (@a Mp).
!> @a RM_SetUnitsExchange specifies how the number of moles of exchange sites in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> 
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-P)*RV.
!> 
!> If a single EXCHANGE definition is used for cells with different initial porosity, 
!>    the three options scale quite differently. 
!> For option 0, the number of moles of exchangers will be the same regardless of porosity. 
!> For option 1, the number of moles of exchangers will be vary directly with porosity and inversely with rock volume. 
!> For option 2, the number of moles of exchangers will vary directly with rock volume and inversely with porosity.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for exchangers: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_InitialPhreeqc2Module, @ref RM_InitialPhreeqcCell2Module,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsExchange(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsExchange(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsExchange(id, option) &
			BIND(C, NAME='RMF_SetUnitsExchange')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsExchange 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsExchange = RMF_SetUnitsExchange(id, option)
END FUNCTION RM_SetUnitsExchange 

!> Set input units for gas phases.
!> In PHREEQC input, gas phases are defined by moles of component gases (@a Mp).
!> @a RM_SetUnitsGasPhase specifies how the number of moles of component gases in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> 
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> 
!> If a single GAS_PHASE definition is used for cells with different initial porosity, 
!>    the three options scale quite differently. 
!> For option 0, the number of moles of a gas component will be the same regardless of porosity. 
!> For option 1, the number of moles of a gas component will be vary directly with porosity and inversely with rock volume. 
!> For option 2, the number of moles of a gas component will vary directly with rock volume and inversely with porosity.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for gas phases: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_InitialPhreeqc2Module, @ref RM_InitialPhreeqcCell2Module,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsGasPhase(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsGasPhase(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsGasPhase(id, option) &
			BIND(C, NAME='RMF_SetUnitsGasPhase')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsGasPhase 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsGasPhase = RMF_SetUnitsGasPhase(id, option)
END FUNCTION RM_SetUnitsGasPhase 

!> Set input units for kinetic reactants.
!> 
!> In PHREEQC input, kinetics are defined by moles of kinetic reactants (@a Mp).
!> @a RM_SetUnitsKinetics specifies how the number of moles of kinetic reactants in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> 
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> 
!> If a single KINETICS definition is used for cells with different initial porosity, 
!>    the three options scale quite differently. 
!> For option 0, the number of moles of kinetic reactants will be the same regardless of porosity. 
!> For option 1, the number of moles of kinetic reactants will be vary directly with porosity and inversely with rock volume. 
!> For option 2, the number of moles of kinetic reactants will vary directly with rock volume and inversely with porosity.
!> 
!> Note that the volume of water in a cell in the reaction module is equal to the product of
!> porosity (@ref RM_SetPorosity), the saturation (@ref RM_SetSaturation), and representative volume (@ref
!> RM_SetRepresentativeVolume), which is usually less than 1 liter. It is important to write the RATES
!> definitions for homogeneous (aqueous) kinetic reactions to account for the current volume of
!> water, often by calculating the rate of reaction per liter of water and multiplying by the volume
!> of water (Basic function SOLN_VOL). 
!> 
!> Rates that depend on surface area of solids, are not dependent
!> on the volume of water. However, it is important to get the correct surface area for the kinetic
!> reaction. To scale the surface area with the number of moles, the specific area (m^2 per mole of reactant) 
!> can be defined as a parameter (KINETICS; -parm), which is multiplied by the number of moles of 
!> reactant (Basic function M) in RATES to obtain the surface area.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for kinetic reactants: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                     @ref RM_InitialPhreeqc2Module, @ref RM_InitialPhreeqcCell2Module,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume, @ref RM_SetSaturation.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsKinetics(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsKinetics(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsKinetics(id, option) &
			BIND(C, NAME='RMF_SetUnitsKinetics')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsKinetics 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsKinetics = RMF_SetUnitsKinetics(id, option)
END FUNCTION RM_SetUnitsKinetics 

!> Set input units for pure phase assemblages (equilibrium phases).
!> In PHREEQC input, equilibrium phases are defined by moles of each phase (@a Mp).
!> @a RM_SetUnitsPPassemblage specifies how the number of moles of phases in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> 
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> 
!> If a single EQUILIBRIUM_PHASES definition is used for cells with different initial porosity, 
!>    the three options scale quite differently. 
!> For option 0, the number of moles of a mineral will be the same regardless of porosity. 
!> For option 1, the number of moles of a mineral will be vary directly with porosity and inversely with rock volume. 
!> For option 2, the number of moles of a mineral will vary directly with rock volume and inversely with porosity.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for equilibrium phases: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_InitialPhreeqc2Module, @ref RM_InitialPhreeqcCell2Module,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsPPassemblage(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsPPassemblage(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsPPassemblage(id, option) &
			BIND(C, NAME='RMF_SetUnitsPPassemblage')  
			USE ISO_C_BINDING 
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsPPassemblage
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsPPassemblage = RMF_SetUnitsPPassemblage(id, option)
END FUNCTION RM_SetUnitsPPassemblage

!> Solution concentration units used by the transport model.
!> Options are 1, mg/L; 2 mol/L; or 3, mass fraction, kg/kgs.
!> PHREEQC defines solutions by the number of moles of each
!> element in the solution.
!> @n@n
!> To convert from mg/L to moles
!> of element in the representative volume of a reaction cell, mg/L is converted to mol/L and
!> multiplied by the solution volume,
!> which is the product of porosity (@ref RM_SetPorosity), saturation (@ref RM_SetSaturation),
!> and representative volume (@ref RM_SetRepresentativeVolume).
!> To convert from mol/L to moles
!> of element in the representative volume of a reaction cell, mol/L is
!> multiplied by the solution volume.
!> To convert from mass fraction to moles
!> of element in the representative volume of a reaction cell, kg/kgs is converted to mol/kgs, multiplied by density
!> (@ref RM_SetDensity) and
!> multiplied by the solution volume.
!> 
!> To convert from moles
!> of element in the representative volume of a reaction cell to mg/L, the number of moles of an element is divided by the
!> solution volume resulting in mol/L, and then converted to mg/L.
!> To convert from moles
!> of element in a cell to mol/L,  the number of moles of an element is divided by the
!> solution volume resulting in mol/L.
!> To convert from moles
!> of element in a cell to mass fraction, the number of moles of an element is converted to kg and divided
!> by the total mass of the solution.
!> Two options are available for the volume and mass of solution
!> that are used in converting to transport concentrations: (1) the volume and mass of solution are
!> calculated by PHREEQC, or (2) the volume of solution is the product of porosity (@ref RM_SetPorosity),
!> saturation (@ref RM_SetSaturation), and representative volume (@ref RM_SetRepresentativeVolume),
!> and the mass of solution is volume times density as defined by @ref RM_SetDensity.
!> Which option is used is determined by @ref RM_UseSolutionDensityVolume.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for solutions: 1, 2, or 3, default is 1, mg/L.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_SetDensity, @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume, @ref RM_SetSaturation,
!> @ref RM_UseSolutionDensityVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsSolution(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsSolution(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsSolution(id, option) &
			BIND(C, NAME='RMF_SetUnitsSolution')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsSolution  
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsSolution = RMF_SetUnitsSolution(id, option)
END FUNCTION RM_SetUnitsSolution  

!> Set input units for solid-solution assemblages.
!> In PHREEQC, solid solutions are defined by moles of each component (@a Mp).
!> @a RM_SetUnitsSSassemblage specifies how the number of moles of solid-solution components in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> 
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for solid solutions: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_InitialPhreeqc2Module, @ref RM_InitialPhreeqcCell2Module,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume.
!> 
!> If a single SOLID_SOLUTION definition is used for cells with different initial porosity, 
!>    the three options scale quite differently. 
!> For option 0, the number of moles of a solid-solution component will be the same regardless of porosity. 
!> For option 1, the number of moles of a solid-solution component will be vary directly with porosity and inversely with rock volume. 
!> For option 2, the number of moles of a solid-solution component will vary directly with rock volume and inversely with porosity.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsSSassemblage(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsSSassemblage(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsSSassemblage(id, option) &
			BIND(C, NAME='RMF_SetUnitsSSassemblage')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsSSassemblage 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsSSassemblage = RMF_SetUnitsSSassemblage(id, option)
END FUNCTION RM_SetUnitsSSassemblage  

!> Set input units for surfaces.
!> In PHREEQC input, surfaces are defined by moles of surface sites (@a Mp).
!> @a RM_SetUnitsSurface specifies how the number of moles of surface sites in a reaction cell (@a Mc)
!> is calculated from the input value (@a Mp).
!> 
!> Options are
!> 0, @a Mp is mol/L of RV (default),    @a Mc = @a Mp*RV, where RV is the representative volume (@ref RM_SetRepresentativeVolume);
!> 1, @a Mp is mol/L of water in the RV, @a Mc = @a Mp*P*RV, where @a P is porosity (@ref RM_SetPorosity); or
!> 2, @a Mp is mol/L of rock in the RV,  @a Mc = @a Mp*(1-@a P)*RV.
!> 
!> If a single SURFACE definition is used for cells with different initial porosity, 
!>    the three options scale quite differently. 
!> For option 0, the number of moles of surface sites will be the same regardless of porosity. 
!> For option 1, the number of moles of surface sites will be vary directly with porosity and inversely with rock volume. 
!> For option 2, the number of moles of surface sites will vary directly with rock volume and inversely with porosity.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param option           Units option for surfaces: 0, 1, or 2.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_InitialPhreeqc2Module, @ref RM_InitialPhreeqcCell2Module,
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetUnitsSurface(id, 1)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SetUnitsSurface(id, option)   
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SetUnitsSurface(id, option) &
			BIND(C, NAME='RMF_SetUnitsSurface')   
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: option
        END FUNCTION RMF_SetUnitsSurface 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: option
    RM_SetUnitsSurface = RMF_SetUnitsSurface(id, option)
END FUNCTION RM_SetUnitsSurface  

!> Set solution concentrations in the reaction cells
!> based on the vector of aqueous species concentrations (@a species_conc).
!> This method is intended for use with multicomponent-diffusion transport calculations,
!> and @ref RM_SetSpeciesSaveOn must be set to @a true.
!> The list of aqueous species is determined by @ref RM_FindComponents and includes all
!> aqueous species that can be made from the set of components.
!> The method determines the total concentration of a component
!> by summing the molarities of the individual species times the stoichiometric
!> coefficient of the element in each species.
!> Solution compositions in the reaction cells are updated with these component concentrations.
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param species_conc     Array of aqueous species concentrations. Dimension of the array is (@a nxyz, @a nspecies),
!> where @a nxyz is the number of user grid cells (@ref RM_GetGridCellCount), and @a nspecies is the number of aqueous species (@ref RM_GetSpeciesCount).
!> Concentrations are moles per liter.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_FindComponents, @ref RM_GetSpeciesConcentrations, @ref RM_GetSpeciesCount, 
!> @ref RM_GetSpeciesD25, @ref RM_GetSpeciesZ,
!> @ref RM_GetSpeciesName, @ref RM_GetSpeciesSaveOn, @ref RM_SetSpeciesSaveOn.
!> 
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_SetSpeciesSaveOn(id, 1)
!> ncomps = RM_FindComponents(id)
!> nspecies = RM_GetSpeciesCount(id)
!> nxyz = RM_GetGridCellCount(id)
!> allocate(species_c(nxyz, nspecies))
!> ...
!> status = RM_SpeciesConcentrations2Module(id, species_c(1,1))
!> status = RM_RunCells(id)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_SpeciesConcentrations2Module(id, species_conc)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_SpeciesConcentrations2Module(id, species_conc) &
			BIND(C, NAME='RMF_SpeciesConcentrations2Module')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            REAL(KIND=C_DOUBLE), INTENT(in) :: species_conc(*)
        END FUNCTION RMF_SpeciesConcentrations2Module
	END INTERFACE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: species_conc
	if (rmf_debug) call Chk_SpeciesConcentrations2Module(id, species_conc)
    RM_SpeciesConcentrations2Module = RMF_SpeciesConcentrations2Module(id, species_conc)
END FUNCTION RM_SpeciesConcentrations2Module  

SUBROUTINE Chk_SpeciesConcentrations2Module(id, species_conc)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: species_conc
    INTEGER :: errors, nspecies
    nspecies = RM_GetSpeciesCount(id)
    errors = 0
    errors = errors + Chk_Double2D(id, species_conc, rmf_nxyz, nspecies, "species_conc", "RM_SpeciesConcentrations2Module")
    if (errors .gt. 0) then
        errors = RM_Abort(id, -3, "Invalid argument in RM_SpeciesConcentrations2Module")
    endif
END SUBROUTINE Chk_SpeciesConcentrations2Module

!> Determines the volume and density to use when converting from the reaction-module concentrations
!> to transport concentrations (@ref RM_GetConcentrations). 
!> Two options are available to convert concentration units: 
!> (1) the density and solution volume calculated by PHREEQC are used, or 
!> (2) the specified density (@ref RM_SetDensity) 
!> and solution volume are defined by the product of 
!> saturation (@ref RM_SetSaturation), porosity (@ref RM_SetPorosity), 
!> and representative volume (@ref RM_SetRepresentativeVolume).
!> Transport models that consider density-dependent flow will probably use the 
!> PHREEQC-calculated density and solution volume (default), 
!> whereas transport models that assume constant-density flow will probably use
!> specified values of density and solution volume. 
!> Only the following databases distributed with PhreeqcRM have molar volume information 
!> needed to accurately calculate density and solution volume: phreeqc.dat, Amm.dat, and pitzer.dat.
!> Density is only used when converting to transport units of mass fraction. 
!> 
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param tf               @a True indicates that the solution density and volume as 
!> calculated by PHREEQC will be used to calculate concentrations. 
!> @a False indicates that the solution density set by @ref RM_SetDensity and the volume determined by the 
!> product of  @ref RM_SetSaturation, @ref RM_SetPorosity, and @ref RM_SetRepresentativeVolume,
!> will be used to calculate concentrations retrieved by @ref RM_GetConcentrations.
!> @see                    @ref RM_GetConcentrations, @ref RM_SetDensity, 
!> @ref RM_SetPorosity, @ref RM_SetRepresentativeVolume, @ref RM_SetSaturation.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_UseSolutionDensityVolume(id, 0)
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root, workers must be in the loop of @ref RM_MpiWorker.

INTEGER FUNCTION RM_UseSolutionDensityVolume(id, tf)
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_UseSolutionDensityVolume(id, tf) &
			BIND(C, NAME='RMF_UseSolutionDensityVolume')
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTEGER(KIND=C_INT), INTENT(in) :: tf
        END FUNCTION RMF_UseSolutionDensityVolume 
	END INTERFACE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in) :: tf
    RM_UseSolutionDensityVolume = RMF_UseSolutionDensityVolume(id, tf)
END FUNCTION RM_UseSolutionDensityVolume 

!> Print a warning message to the screen and the log file.
!> @param id               The instance @a id returned from @ref RM_Create.
!> @param warn_str         String to be printed.
!> @retval IRM_RESULT      0 is success, negative is failure (See @ref RM_DecodeError).
!> @see                    @ref RM_OpenFiles, @ref RM_LogMessage, @ref RM_OutputMessage, @ref RM_ScreenMessage, @ref RM_ErrorMessage.
!> @par Fortran Example:
!> @htmlonly
!> <CODE>
!> <PRE>
!> status = RM_WarningMessage(id, "Parameter is out of range, using default")
!> </PRE>
!> </CODE>
!> @endhtmlonly
!> @par MPI:
!> Called by root and (or) workers; only root writes to the log file.

INTEGER FUNCTION RM_WarningMessage(id, warn_str) 
	USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RMF_WarningMessage(id, warn_str) &
			BIND(C, NAME='RMF_WarningMessage') 
			USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: warn_str(*)
        END FUNCTION RMF_WarningMessage
	END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: warn_str
    RM_WarningMessage = RMF_WarningMessage(id, trim(warn_str)//C_NULL_CHAR)
END FUNCTION RM_WarningMessage

INTEGER FUNCTION Chk_Double1D(id, t, n1, var, func)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:) :: t
    INTEGER, INTENT(in) :: n1
    CHARACTER(len=*), INTENT(in) :: var, func
    CHARACTER(len=200) :: error_string
    INTEGER :: errors, status, t1
    t1 = size(t,1)
    errors = 0
    if (t1 .lt. n1)  then
        errors = errors + 1
        write(error_string, '(A,A,A,I8,A,A)') "Dimension of ", var, " is less than ", n1, " in ", func
        status = RM_ErrorMessage(id, trim(error_string)) 
    endif    
    Chk_Double1D = errors
END FUNCTION Chk_Double1D

INTEGER FUNCTION Chk_Double2D(id, t, n1, n2, var, func)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    DOUBLE PRECISION, INTENT(in), DIMENSION(:,:) :: t
    INTEGER, INTENT(in) :: n1, n2
    CHARACTER(len=*), INTENT(in) :: var, func
    CHARACTER(len=200) :: error_string
    INTEGER :: errors, status, t1, t2
    t1 = size(t,1)
    t2 = size(t,2)
    errors = 0
    if (t2 .ne. n2) then
        errors = errors + 1
        write(error_string, '(A,A,A,I8,A,A)') "Second dimension of ", var, " is not equal to ", n2, " in ", func
        status = RM_ErrorMessage(id, trim(error_string))  
    endif
    if (t1 .lt. n1)  then
        errors = errors + 1
        write(error_string, '(A,A,A,I8,A,A)') "First dimension of ", var, " is less than ", n1, " in ", func
        status = RM_ErrorMessage(id, trim(error_string)) 
    endif    
    Chk_Double2D = errors
END FUNCTION Chk_Double2D

INTEGER FUNCTION Chk_Integer1D(id, t, n1, var, func)
    IMPLICIT NONE
    INTEGER RMF_ErrorMessage
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:) :: t
    INTEGER, INTENT(in) :: n1
    CHARACTER(len=*), INTENT(in) :: var, func
    CHARACTER(len=200) :: error_string
    INTEGER :: errors, status, t1
    t1 = size(t,1)
    errors = 0
    if (t1 .lt. n1)  then
        errors = errors + 1
        write(error_string, '(A,A,A,I8,A,A)') "Dimension of ", var, " is less than ", n1, " in ", func
        status = RM_ErrorMessage(id, trim(error_string)) 
    endif    
    Chk_Integer1D = errors
END FUNCTION Chk_Integer1D

INTEGER FUNCTION Chk_Integer2D(id, t, n1, n2, var, func)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: id
    INTEGER, INTENT(in), DIMENSION(:,:) :: t
    INTEGER, INTENT(in) :: n1, n2
    CHARACTER(len=*), INTENT(in) :: var, func
    CHARACTER(len=200) :: error_string
    INTEGER :: errors, status, t1, t2
    t1 = size(t,1)
    t2 = size(t,2)
    errors = 0
    if (t2 .ne. n2) then
        errors = errors + 1
        write(error_string, '(A,A,A,I8,A,A)') "Second dimension of ", var, " is not equal to ", n2, " in ", func
        status = RM_ErrorMessage(id, trim(error_string))  
    endif
    if (t1 .lt. n1)  then
        errors = errors + 1
        write(error_string, '(A,A,A,I8,A,A)') "First dimension of ", var, " is less than ", n1, " in ", func
        status = RM_ErrorMessage(id, trim(error_string)) 
    endif    
    Chk_Integer2D = errors
END FUNCTION Chk_Integer2D

END MODULE PhreeqcRM

    

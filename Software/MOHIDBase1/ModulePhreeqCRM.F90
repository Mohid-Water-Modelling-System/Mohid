!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : PhreeqCRM
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2015
! REVISION      : Eduardo Jauch - v4.0
! DESCRIPTION   : Zero-dimensional model for chemistry equilibrium of solution, 
!                 pure phases, gas phase, solid phase, exchangers and surfaces
!                 
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

Module ModulePhreeqCRM

    use ModuleFunctions            
    use ModuleEnterData
    use ModuleGlobalData
    use ModuleTime
    use PhreeqcRM
    use IPhreeqc
    
    implicit none

    private 

    !Constructor---------------------------------------------------------------
    public  :: ConstructPhreeqCRM
    private ::      AllocateInstance
    private ::      ReadInputFile
    private ::          ReadConfiguration
    private ::              ReadInitialSolution
    private ::              ReadUnits
    private ::                  SelectUnitsID
    private ::              ConstructPhreeqcList
    private ::                  ConstructPhreeqc
    private ::                      ReadPhreeqCInitializationConfig
    private ::                      CreateMapping
    private ::                      LoadPhases
    private ::                          ConstructPhase
    private ::              ConstructPropertyList
    private ::                  ConstructProperty
    public  :: InitializePhreeqCRM
    private ::      StartPhreeqCRMInstance
    private ::      InititalizePhreeqCCells
    
    !Selector------------------------------------------------------------------    
    public  :: GetPhreeqCIndex
    public  :: GetPhreeqCDT
    public  :: GetPhreeqCSize
    !public  :: GetPhreeqCOutputSize
    public  :: GetPhreeqCPropertiesList
    public  :: UngetPhreeqC
    
    !Modifier------------------------------------------------------------------
    public  :: ModifyPhreeqCRM
    
    !Destructor----------------------------------------------------------------
    public  :: KillPhreeqCRM                                                     
    private ::      DeAllocateInstance

    !Management----------------------------------------------------------------
    private ::      Ready
    private ::          LocateObjPhreeqCRM 
    
    !Parameters----------------------------------------------------------------
    integer, parameter                  :: mg_L       = 1 !mg/L   => for solutions only
    integer, parameter                  :: mol_L      = 2 !mol/L  => for solutions only
    integer, parameter                  :: kg_kgs     = 3 !kg/kgs => for solutions only
    integer, parameter                  :: mol_Lcell  = 0 !mol/L cell
    integer, parameter                  :: mol_Lwater = 1 !mol/L water
    integer, parameter                  :: mol_kgrock = 2 !mol/kg rock
    
    character(*), parameter             :: ClearContentsStr = "DELETE; -all"
    
    !Data Types----------------------------------------------------------------
    
    private :: T_Property
    type T_Property

        type(T_PropertyID)              :: ID
        
        character(StringLength)         :: PhreeqCName = "*****"        
        integer                         :: PhreeqCIndex = -1
        
        logical                         :: IsOutputOnly = .false.
        real                            :: GFW = 0.0        
    end type T_Property
    
    private :: T_PositionTouple
    type T_PositionTouple
        character(StringLength)         :: Name
        integer                         :: Position
    end type T_PositionTouple
    
    public :: T_Point
    type T_Point
        integer                         :: I = -1
        integer                         :: J = -1
        integer                         :: K = -1
    end type T_Point
    
    private :: T_Line
    type T_Line
        character(len=:), allocatable   :: Text
    end type T_Line
    
    private :: T_PhreeqcInitialSolutions
    type T_PhreeqcInitialSolutions
        type(T_Line), allocatable       :: Lines(:)
        integer                         :: NumberOfLines = 0
        logical                         :: UseThis = .false.
    end type
    
    private :: T_PhaseLine
    type T_PhaseLine
        character(len=:), allocatable   :: LineStart
        character(len=:), allocatable   :: LineEnd
        integer                         :: PropertyIndex = -1
        logical                         :: CellIndex = .false.
    end type
    
    private :: T_Phase
    type T_Phase
        character(len=:), allocatable   :: Name
        type(T_PhaseLine), allocatable  :: Lines(:)
        integer                         :: NumberOfLines = 0      
    end type
    
    !T_Phreeqc stores the required information to run PhreeqC on a specific phreeqc in the model.
    !Right now, due the characteristics of INTERFACE module, only ONE phreeqc is allowed
    private :: T_Phreeqc
    type :: T_Phreeqc
        
        type(T_PropertyID)              :: ID
        
        type(T_Size3D)                  :: WorkSize
        logical                         :: FullArea = .true.
        
        integer                         :: NumberOfCells = 0
        integer                         :: NumberOfComponents = 0
        integer                         :: NumberOfOutputs = 0
        
        integer, allocatable            :: ComponentsToProperties(:)
        integer, allocatable            :: PropertiesToComponents(:)
        
        integer, allocatable            :: OutputsToProperties(:)
        integer, allocatable            :: PropertiesToOutputs(:)
        
        integer, pointer                :: Mapping (:) => null()
        
        !Index 1: 1 Solution, 2 Equilibrium, 3 Exchange, 4 Surface
        !         5 Gas Phase, 6 Solid Solution, 7 Kinectics
        !Index 2.1 C1: Initial solution/reactant, 
        !Index 2.2 C2: Solution/reactant to mix with C1,         
        integer                         :: Initialization(7, 2) = -1
        real                            :: MixFraction(7) = 1.0
        
        type(T_Phase), allocatable      :: Phases(:)
        integer                         :: NumberOfPhases
        
        real, pointer                   :: Concentrations(:,:) => null()
        real, pointer                   :: Outputs(:,:) => null()
        
        integer                         :: ErrorHandlerMode = 2 !Exit on error
        logical                         :: UseSolutionDensityVolume = .false.
        logical                         :: ExcessHO = .false.
        
        character(PathLength)           :: Database
        character(PathLength)           :: InputFile
        
        integer                         :: ObjPhreeqC = -1
        
        integer                         :: UserOutput
        
    end type T_Phreeqc
    
    private :: T_PhreeqCRM
    type :: T_PhreeqCRM
        
        integer                         :: InstanceID = 0
        
        integer                         :: ObjEnterData = 0
        
        integer                         :: NumberOfGridCells = 0
               
        integer                         :: NumberOfPhreeqcs = 0
        type(T_Phreeqc), pointer        :: Phreeqcs(:) => null()
        
        integer                         :: NumberOfProperties = 0
        integer                         :: NumberOfComponents = 0
        integer                         :: NumberOfOutputs = 0
                
        type(T_Property), pointer       :: Properties(:) => null()
        integer, pointer                :: PropertiesID(:) => null()
        
        integer                         :: UnitsSolution = -1
        integer                         :: UnitsGasPhase = -1
        integer                         :: UnitsExchange = -1
        integer                         :: UnitsSSassemblage = -1
        integer                         :: UnitsPPassemblage = -1
        integer                         :: UnitsKinetics = -1
        integer                         :: UnitsSurface = -1
        
        integer                         :: ArrayLB = 0
        integer                         :: ArrayUB = 0
        
        type(T_Time)                    :: BeginTime
        type(T_Time)                    :: EndTime
        type(T_Time)                    :: CurrentTime
        real                            :: ElapsedTime = 0.0
        real                            :: CurrentDT = 0.0
        
        integer                         :: Dimensions = 0 !1 for 1D, 2 for 2D and 3 for 3D
        
        character(PathLength)           :: DefaultDatabase = "*****"
        character(PathLength)           :: DefaultInputFile = "*****"
        
        type(T_Size3D)                  :: WorkSize
        
        logical                         :: LogOnError = .true.
        character(PathLength)           :: LogErrorFile = "log_error_phreeqc.txt"
        
        real                            :: PhreeqCDT
        
        type(T_PhreeqcInitialSolutions) :: InitialSolutions
        
        real, pointer                   :: SolutionMapping(:) => null()
        real, pointer                   :: EquilibriumMapping(:) => null()
        real, pointer                   :: ExchangeMapping(:) => null()
        real, pointer                   :: SurfaceMapping(:) => null()
        real, pointer                   :: GasPhaseMapping(:) => null()
        real, pointer                   :: SolidSolutionMapping(:) => null()
        real, pointer                   :: KineticsMapping(:) => null()
        
        type(T_PhreeqCRM), pointer      :: Next => null()
        type(T_PhreeqCRM), pointer      :: Prior => null()
        
    end type   T_PhreeqCRM
    
    !Global Module Variables---------------------------------------------------
    type (T_PhreeqCRM), pointer :: FirstObjPhreeqCRM
    type (T_PhreeqCRM), pointer :: Me
    
    !Subroutines---------------------------------------------------------------

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine ConstructPhreeqCRM (ObjID,                   &
                                   file_name,               &
                                   begin_time,              &
                                   end_time,                &
                                   array_lb,                &
                                   array_ub,                &
                                   stat)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjID  
        character(len=*)                                :: file_name
        type(T_Time)                                    :: begin_time
        type(T_Time)                                    :: end_time
        integer                                         :: array_lb
        integer                                         :: array_ub
        integer, optional, intent(OUT)                  :: stat

        !External----------------------------------------------------------------
        integer                                         :: ready_, stat_call

        !Local-------------------------------------------------------------------
        integer                                         :: stat_

        !------------------------------------------------------------------------

        stat_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mPHREEQC_)) then
            nullify (FirstObjPhreeqCRM)
            call RegisterModule (mPHREEQC_) 
        endif
        
        call Ready (ObjID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%NumberOfGridCells = array_ub - array_lb + 1
            !Me%WorkSize = worksize
            Me%BeginTime = begin_time
            Me%EndTime = end_time
            
            Me%ArrayLB = array_lb
            Me%ArrayUB = array_ub

            call ConstructEnterData(Me%ObjEnterData, file_name, STAT = stat_call) 
            if (stat_call .NE. SUCCESS_) then
                print *, "Was not possible to read the input file '", trim(file_name), "'."
                stop "ConstructPhreeqCRM - ModulePhreeqCRM - ERR-002"
            endif
            
            call ReadInputFile
       
            call KillEnterData(Me%ObjEnterData, STAT = stat_call) 
            if (stat_call .NE. SUCCESS_) then
                print *, "Was not possible to close the input file reading object."
                stop "ConstructPhreeqCRM - ModulePhreeqCRM - ERR-003"
            endif

            !Returns ID
            ObjID = Me%InstanceID

            stat_ = SUCCESS_
        else 
            
            print *, "Was not possible to construct the PhreeqCRM module"
            stop "ConstructPhreeqCRM - ModulePhreeqCRM - ERR-004"

        end if cd0


        if (present(STAT)) STAT = stat_

        !----------------------------------------------------------------------

    end subroutine ConstructPhreeqCRM
    
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_PhreeqCRM), pointer                             :: NewObject
        type (T_PhreeqCRM), pointer                             :: PreviousObj

        !----------------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewObject)
        nullify  (NewObject%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPhreeqCRM)) then
            
            FirstObjPhreeqCRM => NewObject
            Me                => NewObject
            
        else
            
            PreviousObj => FirstObjPhreeqCRM
            Me          => FirstObjPhreeqCRM%Next
            
            do while (associated(Me))
                
                PreviousObj => Me
                Me          => Me%Next
                
            enddo
            
            Me               => NewObject
            PreviousObj%Next => NewObject
            
        endif

        Me%InstanceID = RegisterNewInstance (mPhreeqC_)

    end subroutine AllocateInstance
    
    !--------------------------------------------------------------------------
      
    subroutine ReadInputFile
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
                
        call ReadConfiguration        
        call ConstructPropertyList
        call ConstructPhreeqcList
        
        !----------------------------------------------------------------------
    
    end subroutine ReadInputFile
    
    !--------------------------------------------------------------------------
    
    subroutine ReadConfiguration
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: iflag
        integer                                         :: stat_call
        integer                                         :: client_number
        integer                                         :: first
        integer                                         :: last
        logical                                         :: found
        
        !----------------------------------------------------------------------
        
        call ReadUnits
        
        call GetData (Me%DefaultDatabase,                           &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'DEFAULT_DATABASE',           &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .ne. SUCCESS_) then
            print *, "Was not possible to read the DEFAULT_DATABASE keyword."
            stop "ReadConfiguration - ModulePhreeqCRM - ERR-001"
        endif
        
        call ExtractBlockFromBuffer (Me%ObjEnterData,                   &
                                     ClientNumber    = client_number,   &
                                     block_begin     = '<begininit>',   &
                                     block_end       = '<endinit>',     &
                                     BlockFound      = found,           &
                                     FirstLine       = first,           &
                                     LastLine        = last,            &            
                                     STAT            = stat_call)
        if (stat_call /= SUCCESS_) then
            print *, "An unspecified error happened when retrieving the init block."
            stop "ReadConfiguration - ModulePhreeqCRM - ERR-002"
        endif        
        if (found) then
            Me%InitialSolutions%UseThis = .true.
            call ReadInitialSolution (first + 1, last - 1)
        else
            Me%InitialSolutions%UseThis = .false.            
        endif
        call Block_Unlock(Me%ObjEnterData, client_number, STAT = stat_call)
        call RewindBuffer(Me%ObjEnterData, STAT = stat_call)
        if (stat_call .NE. SUCCESS_) then
            print *, "An unspecified error happened when rewinding the input data file buffer."
            stop 'ReadConfiguration - ModulePhreeqCRM - ERR-003'
        endif
        
        if (.not. Me%InitialSolutions%UseThis) then
        
            call GetData (Me%DefaultInputFile,                          &
                          Me%ObjEnterData, iflag,                       &
                          Keyword       = 'DEFAULT_INPUT_FILE',         &
                          SearchType    = FromFile,                     &
                          ClientModule  = 'ModulePhreeqCRM',            &
                          STAT          = stat_call)
            if (stat_call .ne. SUCCESS_) then
                print *, "Was not possible to read the DEFAULT_INPUT_FILE keyword."
                stop "ReadConfiguration - ModulePhreeqCRM - ERR-004"
            endif

        endif
            
        call GetData (Me%PhreeqCDT,                                 &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'MODULE_DT',                  &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .ne. SUCCESS_ .or. iflag /= 1) then
            print *, "Was not possible to read the MODULE_DT keyword."
            stop "ReadConfiguration - ModulePhreeqCRM - ERR-005"
        endif                
        
        !----------------------------------------------------------------------
    
    end subroutine ReadConfiguration
    
    !--------------------------------------------------------------------------
    
    subroutine ReadInitialSolution (first, last)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: first
        integer                                         :: last
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: iflag
        integer                                         :: status
        type(T_PhreeqcInitialSolutions), pointer        :: is
        integer                                         :: index
        integer                                         :: l
        character(StringLength)                         :: item
        
        !----------------------------------------------------------------------    

        is => Me%InitialSolutions
        
        is%NumberOfLines = last - first + 1
        allocate (is%Lines (is%NumberOfLines))
        
        index = 0
        do l = first, last
            index = index + 1
            call GetData(item,                                      &
                         Me%ObjEnterData, iflag,                    &
                         Buffer_Line  = l,                          &
                         ClientModule = 'ModulePhreeqCRM',          &
                         STAT         = status)               
            if (status /= SUCCESS_) &
                stop 'ConstructPhase - ModuleWaterProperties - ERR-010'
            
            is%Lines(index)%Text = trim(item)
        enddo        
        
        !----------------------------------------------------------------------    
    
    end subroutine ReadInitialSolution
    
    !--------------------------------------------------------------------------
    
    subroutine ReadUnits
    
        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        character(StringLength)                         :: aux
        integer                                         :: iflag
        integer                                         :: stat_call
        
        !----------------------------------------------------------------------
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_SOLUTION',             &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_SOLUTION keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-001"
        endif
        Me%UnitsSolution = SelectUnitsID (aux, is_solution = .true.)
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_GASPHASE',             &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call) 
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_GASPHASE keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-002"
        endif
        Me%UnitsGasPhase = SelectUnitsID (aux, is_solution = .false.)
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_EXCHANGE',             &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call) 
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_EXCHANGE keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-003"
        endif
        Me%UnitsExchange = SelectUnitsID (aux, is_solution = .false.)
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_SSASSEMBLAGE',         &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call) 
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_SSASSEMBLAGE keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-004"
        endif
        Me%UnitsSSassemblage = SelectUnitsID (aux, is_solution = .false.)
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_PPASSEMBLAGE',         &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call) 
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_PPASSEMBLAGE keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-005"
        endif
        Me%UnitsPPassemblage = SelectUnitsID (aux, is_solution = .false.)
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_KINETICS',             &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call) 
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_KINECTICS keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-006"
        endif
        Me%UnitsKinetics = SelectUnitsID (aux, is_solution = .false.)
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'UNITS_SURFACE',              &
                      SearchType    = FromFile,                     &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call) 
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the UNITS_SURFACE keyword."
            stop "ReadUnits - ModulePhreeqCRM - ERR-007"
        endif
        Me%UnitsSurface = SelectUnitsID (aux, is_solution = .false.)
        !----------------------------------------------------------------------    
    
    end subroutine ReadUnits
    
    !--------------------------------------------------------------------------
    
    function SelectUnitsID (units, is_solution) result(res)
    
        !Arguments-------------------------------------------------------------
        character(*)                                    :: units
        logical, optional, intent(IN)                   :: is_solution
        integer                                         :: res
        
        !Local-----------------------------------------------------------------
        logical                                         :: is_solution_
        
        !----------------------------------------------------------------------
        
        if (present(is_solution)) then
            is_solution_ = is_solution
        else
            is_solution_ = .false.
        endif
        
        if (trim(units) .eq. "mg/L") then
            
            if (.not. is_solution) then
                print *, "The units 'mg/L' can be applied only for UNITS_SOLUTION"
                print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
                print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords." 
                stop "SelectUnitsID - ModulePhreeqCRM - ERR-001"
            endif
            
            res = mg_L
            
        else if (trim(units) .eq. "mol/L") then
            
            if (.not. is_solution) then
                print *, "The units 'mol/L' can be applied only for UNITS_SOLUTION"
                print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
                print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords." 
                stop "SelectUnitsID - ModulePhreeqCRM - ERR-002"
            endif            
            
            res = mol_L
            
        else if (trim(units) .eq. "kg/kgs") then
            
            if (.not. is_solution) then
                print *, "The units 'kg/kgs' can be applied only for UNITS_SOLUTION"
                print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
                print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords." 
                stop "SelectUnitsID - ModulePhreeqCRM - ERR-003"                
            endif     
            
            res = kg_kgs
            
        else if (trim(units) .eq. "mol/Lcell") then
            
            if (is_solution) then
                print *, "The units 'mol/Lcell' can not be applied for UNITS_SOLUTION"
                print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
                print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords." 
                stop "SelectUnitsID - ModulePhreeqCRM - ERR-004"  
            endif
            
            res = mol_Lcell
            
        else if (trim(units) .eq. "mol/Lwater") then
            
            if (is_solution) then
                print *, "The units 'mol/Lwater' can not be applied for UNITS_SOLUTION"
                print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
                print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords." 
                stop "SelectUnitsID - ModulePhreeqCRM - ERR-005" 
            endif   
            
            res = mol_Lwater
            
        else if (trim(units) .eq. "mol/kgrock") then
            
            if (is_solution) then
                print *, "The units 'mol/kgrock' can not be applied for UNITS_SOLUTION"
                print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
                print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords."                
                stop "SelectUnitsID - ModulePhreeqCRM - ERR-006" 
            endif
            
            res = mol_kgrock
            
        else
            
            print *, "Unknown units '", trim(units), "'."
            print *, "Use mg/L, mol/L and kg/kgs for UNITS_SOLUTION and"
            print *, "mol/Lcell, mol/Lwater and mol/kgrock for the other UNITS_* keywords."
            stop "SelectUnitsID - ModulePhreeqCRM - ERR-007" 
            
        endif
        
        !----------------------------------------------------------------------
    
    end function SelectUnitsID
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructPhreeqcList
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer                                         :: index
        logical                                         :: found
        integer                                         :: client_number
        type(T_Phreeqc), pointer                        :: phreeqc
        
        !----------------------------------------------------------------------
        call GetNumberOfBlocks (Me%ObjEnterData, '<beginphreeqc>', '<endphreeqc>',  &
                                FromFile_,                                          &
                                Me%NumberOfPhreeqcs,                                   &
                                stat = stat_call)        
        if (stat_call /= SUCCESS_) then
            print *, "An unspecified error happened when retrieving the number of phreeqc blocks"
            print *, "in the imput file."
            print *, "The error returned was: ", stat_call
            stop "ConstructPhreeqcsList - ModulePhreeqCRM - ERR-001"
        endif
    
        !if (Me%NumberOfPhreeqcs > 1) then
        !    print *, "Due to curent model characteristics, no more than ONE phreeqc can be setup."
        !    print *, "It was found the following number of 'phreeqc' blocks: ", Me%NumberOfPhreeqcs
        !    stop "ConstructPhreeqcsList - ModulePhreeqCRM - ERR-002"
        !endif
        
        if (Me%NumberOfPhreeqcs > 0) then
        
            allocate (Me%Phreeqcs(Me%NumberOfPhreeqcs))
            
            do index = 1, Me%NumberOfPhreeqcs
            
                phreeqc => Me%Phreeqcs(index)
                phreeqc%ID%IDNumber = index
                
                call ExtractBlockFromBuffer (Me%ObjEnterData,                                           &
                                             ClientNumber    = client_number,                           &
                                             block_begin     = '<beginphreeqc>',                        &
                                             block_end       = '<endphreeqc>',                          &
                                             BlockFound      = found,                                   &
                                             STAT            = stat_call)
                if (stat_call /= SUCCESS_) then
                    print *, "An unspecified error happened when retrieving the phreeqc block number '", index, "'."
                    stop "ConstructPhreeqcsList - ModulePhreeqCRM - ERR-003"
                elseif (.not. found) then
                    print *, "An phreeqc block was expected, but it was not found."
                    stop "ConstructPhreeqcsList - ModulePhreeqCRM - ERR-004"
                endif
                
                call ConstructPhreeqc (phreeqc, index, client_number, stat_call) 
                if (stat_call /= SUCCESS_) then
                    print *, "It was not possible to construct the property number '", index, "'."
                    stop "ConstructPhreeqcsList - ModulePhreeqCRM - ERR-005"                    
                endif
                
            enddo
            
            call Block_Unlock(Me%ObjEnterData, client_number, STAT = stat_call)
        
            call RewindBuffer(Me%ObjEnterData, STAT = stat_call)
            if (stat_call .NE. SUCCESS_) then
                print *, "An unspecified error happened when rewinding the input data file buffer."
                stop 'ConstructPhreeqcsList - ModulePhreeqCRM - ERR-006'
            endif
            
        else
            
            print *, "It is necessary to have at least ONE 'phreeqc' block defined in the MOHID PhreeqC input data file."
            stop 'ConstructPhreeqcsList - ModulePhreeqCRM - ERR-007'
                
        endif
    
    end subroutine ConstructPhreeqcList
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructPhreeqc (phreeqc, index, client_number, stat)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: client_number
        integer                                         :: index
        integer                                         :: stat
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer                                         :: limits(6)
        integer                                         :: number_of_limits
        integer                                         :: iflag
        
        !----------------------------------------------------------------------
        
        stat = SUCCESS_
        
        call GetData (phreeqc%ID%Name,                              &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'NAME',                       &
                      SearchType    = FromBlock,                    &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the NAME keyword."
            print *, "Phreeqc number: ", index
            stat = UNKNOWN_
        endif

        call GetData (phreeqc%ErrorHandlerMode,                     &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'ERROR_HANDLER_MODE',         &
                      SearchType    = FromBlock,                    &
                      Default       = 2,                            &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the ERROR_HANDLER_MODE keyword."
            print *, "Phreeqc number: ", index
            stat = UNKNOWN_
        endif
        
        call GetData (phreeqc%UseSolutionDensityVolume,             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'USE_SOLUTION_DENSITY_VOLUME',&
                      SearchType    = FromBlock,                    &
                      Default       = .false.,                      &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the USE_SOLUTION_DENSITY_VOLUME keyword."
            print *, "Phreeqc number: ", index
            stat = UNKNOWN_
        endif
        
        call GetData (phreeqc%ExcessHO,                             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'SEPARATE_H_O_EXCESS',        &
                      SearchType    = FromBlock,                    &
                      Default       = .false.,                      &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the SEPARATE_H_O_EXCESS keyword."
            print *, "Phreeqc number: ", index
            stat = UNKNOWN_
        endif
        
        call GetData (phreeqc%Database,                             &
                        Me%ObjEnterData, iflag,                       &
                        Keyword       = 'DATABASE',                   &
                        SearchType    = FromBlock,                    &
                        Default       = Me%DefaultDatabase,           &
                        ClientModule  = 'ModulePhreeqCRM',            &
                        STAT          = stat_call)
        if (stat_call /= SUCCESS_) then
            print *, "Was not possible to read the DATABASE keyword."
            print *, "Phreeqc number: ", index
            stat = UNKNOWN_
        endif
        
        if (.not. Me%InitialSolutions%UseThis) then
        
            call GetData (phreeqc%InputFile,                            &
                            Me%ObjEnterData, iflag,                       &
                            Keyword       = 'INPUT_FILE',                 &
                            SearchType    = FromBlock,                    &
                            Default       = Me%DefaultInputFile,          &
                            ClientModule  = 'ModulePhreeqCRM',            &
                            STAT          = stat_call)
            if (stat_call /= SUCCESS_) then
                print *, "Was not possible to read the INPUT_FILE keyword."
                print *, "Phreeqc number: ", index
                stat = UNKNOWN_
            endif

        endif
            
        call GetData (phreeqc%UserOutput,                           &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'USER_OUTPUT',                &
                      SearchType    = FromBlock,                    &
                      Default       = -1,                           &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call /= SUCCESS_) then
            print *, "Was not possible to read the USER_OUTPUT keyword."
            print *, "Phreeqc number: ", index
            stat = UNKNOWN_
        endif
        
        stat = ReadPhreeqCInitializationConfig ("SOLUTION", 1, phreeqc)
        if (stat /= SUCCESS_) return

        stat = ReadPhreeqCInitializationConfig ("EQUILIBRIUM", 2, phreeqc)
        if (stat /= SUCCESS_) return

        stat = ReadPhreeqCInitializationConfig ("EXCHANGE", 3, phreeqc)
        if (stat /= SUCCESS_) return

        stat = ReadPhreeqCInitializationConfig ("SURFACE", 4, phreeqc)
        if (stat /= SUCCESS_) return

        stat = ReadPhreeqCInitializationConfig ("GAS_PHASE", 5, phreeqc)
        if (stat /= SUCCESS_) return
        
        stat = ReadPhreeqCInitializationConfig ("SOLID_SOLUTION", 6, phreeqc)
        if (stat /= SUCCESS_) return
        
        stat = ReadPhreeqCInitializationConfig ("KINETICS", 7, phreeqc)
        if (stat /= SUCCESS_) return

        call GetData (limits,                                       &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'LIMITS',                     &
                      SearchType    = FromBlock,                    &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call /= SUCCESS_ .and. stat_call /= KEYWORD_NOT_FOUND_ERR_) then
            print *, "An unknown error happened when trying to read the LIMITS keyword."
            print *, "The error returned was: ", stat_call
            print *, "Phreeqc number: ", index
            print *, "Phreeqc name: ", trim(phreeqc%ID%Name)
            stat = UNKNOWN_
        endif
        
        if (iflag == 0) then
            
            phreeqc%FullArea = .true.
            
            phreeqc%NumberOfCells = Me%ArrayUB - Me%ArrayLB + 1
            print *, "Phreeqc :", phreeqc%ID%IDNumber
            print *, "Number of cells: ", phreeqc%NumberOfCells
            
        else
            
            if (iflag == 2) then 
                number_of_limits = 1
            elseif (iflag == 4) then 
                number_of_limits = 2
            elseif (iflag == 6) then 
                number_of_limits = 3
            else
                print *, "The number of limits must be 2 (ilb iub), 4(ilb iub jlb jub) or 6(ilb iub jlb jub klb kub)."
                print *, "Phreeqc number: ", index
                print *, "Phreeqc name: ", trim(phreeqc%ID%Name)
                stat = UNKNOWN_                                
            endif
            
            phreeqc%FullArea = .false.
            
        endif
                
        allocate (phreeqc%ComponentsToProperties (Me%NumberOfProperties))
        phreeqc%ComponentsToProperties = -1
        
        allocate (phreeqc%PropertiesToComponents (Me%NumberOfProperties))
        phreeqc%PropertiesToComponents = -1
        
        allocate (phreeqc%OutputsToProperties (Me%NumberOfProperties))
        phreeqc%OutputsToProperties = -1
        
        allocate (phreeqc%PropertiesToOutputs (Me%NumberOfProperties))
        phreeqc%PropertiesToOutputs = -1
        
        !call CreateMapping (phreeqc, limits, number_of_limits)
        call CreateMapping (phreeqc)
        
        call LoadPhases (phreeqc, client_number)
        
        !----------------------------------------------------------------------
        
    end subroutine ConstructPhreeqc
    
    !--------------------------------------------------------------------------
    
    function ReadPhreeqCInitializationConfig (keyword, index, phreeqc) result (stat)
    
        !Arguments-------------------------------------------------------------
        character(*)                                    :: keyword
        integer                                         :: index
        type(T_Phreeqc)                                    :: phreeqc
        integer                                         :: stat
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag
        integer                                         :: stat_call
        real                                            :: aux(3)
        
        !----------------------------------------------------------------------
        
        stat = SUCCESS_
        
        call GetData (aux,                                          &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = trim(keyword),                &
                      SearchType    = FromBlock,                    &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call /= SUCCESS_ .or. iflag /= 3) then
            print *, "Was not possible to read the keyword '", trim(keyword), "'."
            print *, "Phreeqc number: ", phreeqc%ID%IDNumber
            stat = UNKNOWN_
        endif
        
        if (stat==SUCCESS_) then
            if (iflag /= 1 .and. iflag /= 3) then
                print *, "Invalid keyword '", trim(keyword), "'."
                print *, "Phreeqc number: ", phreeqc%ID%IDNumber
                print *, "It must contain 1 ID or 2 ID's + 1 Fraction" 
                stat = UNKNOWN_
            else
                phreeqc%Initialization (index, 1) = int(aux (1))
                phreeqc%Initialization (index, 2) = int(aux (2))
                phreeqc%MixFraction (index) = aux (3)
            endif               
        endif
        
        !----------------------------------------------------------------------
    
    end function ReadPhreeqCInitializationConfig
    
    !--------------------------------------------------------------------------
    
    !subroutine CreateMapping (phreeqc, limits, number_of_limits)
    subroutine CreateMapping (phreeqc)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                           :: phreeqc
        !integer                                         :: limits(6)
        !integer                                         :: number_of_limits
    
        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        
        if (phreeqc%FullArea) then
            
            allocate (phreeqc%Mapping(Me%ArrayLB:Me%ArrayUB))
            phreeqc%Mapping = 1;
            
        else
            
            !ToDo: Adapt ModuleInterface so it is possible to create different phreeqcs for PhreeqC
            !ToDo: Prepare the mapping of different phreeqcs here
            
        endif
        
        !----------------------------------------------------------------------
    
    end subroutine CreateMapping
    
    !--------------------------------------------------------------------------
    
    subroutine LoadPhases (phreeqc, client_number)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: client_number
    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer                                         :: i
        integer                                         :: first
        integer                                         :: last
        logical                                         :: found

        !----------------------------------------------------------------------
        
        print *, "Constructing ModulePhreeqC Phases..."
        print *, ""
        
        call GetNumberOfBlocks (Me%ObjEnterData, '<beginphase>', '<endphase>',          &
                                FromBlock_,                                             &
                                phreeqc%NumberOfPhases,                                 &
                                ClientNumber = client_number,                           &
                                stat = stat_call)        
        if (stat_call /= SUCCESS_) then
            print *, "An unspecified error happened when retrieving the number of phases"
            print *, "in the input file."
            print *, "The error returned was: ", stat_call
            stop "LoadPhases - ModulePhreeqCRM - ERR-001"
        endif
    
        if (phreeqc%NumberOfPhases > 0) then
        
            allocate (phreeqc%Phases (phreeqc%NumberOfPhases))
        
            do i = 1, phreeqc%NumberOfPhases
            
                call ExtractBlockFromBlock (Me%ObjEnterData,                                            &
                                            ClientNumber        = client_number,                        &
                                            block_begin         = '<beginphase>',                       &
                                            block_end           = '<endphase>',                         &
                                            BlockInBlockFound   = found,                                &
                                            FirstLine           = first,                                &
                                            LastLine            = last,                                 &
                                            STAT                = stat_call)
                if (stat_call /= SUCCESS_) then
                    print *, "An unspecified error happened when retrieving the phase block number '", i, "'."
                    stop "LoadPhases - ModulePhreeqCRM - ERR-002"
                elseif (.not. found) then
                    print *, "A phase block was expected, but it was not found."
                    stop "LoadPhases - ModulePhreeqCRM - ERR-003"
                endif
                
                call ConstructPhase (phreeqc, i, first, last, stat_call) 
                if (stat_call /= SUCCESS_) then
                    print *, "It was not possible to construct the phase number '", i, "'."
                    stop "LoadPhases - ModulePhreeqCRM - ERR-004"                    
                endif
                
            enddo
            
            call Block_Unlock(Me%ObjEnterData, client_number, STAT = stat_call)
        
            !call RewindBuffer(Me%ObjEnterData, STAT = stat_call)
            !if (stat_call .NE. SUCCESS_) then
            !    print *, "An unspecified error happened when rewinding the input data file buffer."
            !    stop 'LoadPhases - ModulePhreeqCRM - ERR-005'
            !endif
        
        else
            
            print *, "No Phases are present."
            
        endif
            
        !----------------------------------------------------------------------
        
    end subroutine
    
    !--------------------------------------------------------------------------

    subroutine ConstructPhase (phreeqc, phase_index, first, last, stat)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: phase_index
        integer                                         :: first
        integer                                         :: last
        integer                                         :: stat
                                                    
        !Local-----------------------------------------------------------------
        type(T_Phase), pointer                          :: phase
        character(StringLength)                         :: item
        integer                                         :: l
        integer                                         :: status
        integer                                         :: index
        integer                                         :: iflag
        
        !----------------------------------------------------------------------
        
        stat = SUCCESS_
        
        phase => phreeqc%Phases(phase_index)
        
        phase%NumberOfLines = last - first - 2
        allocate (phase%Lines (phase%NumberOfLines))
        
        index = 0
        do l = first + 1, last - 1
            index = index + 1
            call GetData(item,                                      &
                         Me%ObjEnterData, iflag,                    &
                         Buffer_Line  = l,                          &
                         ClientModule = 'ModulePhreeqCRM',          &
                         STAT         = status)               
            if (status /= SUCCESS_) &
                stop 'ConstructPhase - ModuleWaterProperties - ERR-010'
            
            if (index == 1) then
                phase%Name = trim(item)
            else
                call ProcessLine(phase, index - 1, item)
            endif
        enddo
        
        !----------------------------------------------------------------------
        
    end subroutine ConstructPhase
        
    !--------------------------------------------------------------------------
    
    subroutine ProcessLine(phase, line_index, text)
    
        !Arguments-------------------------------------------------------------
        type(T_Phase), pointer                          :: phase
        integer                                         :: line_index
        character(len=*)                                :: text
                                                    
        !Local-----------------------------------------------------------------
        character(len=:), allocatable                   :: trimmed_line
        integer                                         :: i_start
        integer                                         :: i_end
        integer                                         :: line_size
        logical                                         :: process_line
        integer                                         :: i
        
        !----------------------------------------------------------------------
        trimmed_line = trim(text)
        
        i_start = -1
        i_end = -1
        process_line = .false.
        
        line_size = len(trimmed_line)
        
        if (line_size <= 4) then
            phase%Lines(line_index)%LineStart = trimmed_line
            phase%Lines(line_index)%LineEnd = ""
            return
        endif
        
        
do1:    do i = 1, line_size - 1
            if (trimmed_line(i:i) == '{' .and. trimmed_line(i+1:i+1)=='{') then
                i_start = i + 2
                exit do1
            endif
        enddo do1
    
        if (i_start > 0) then
do2:        do i = i_start, len(trimmed_line) - 1
                if (trimmed_line(i:i) == '}' .and. trimmed_line(i+1:i+1) == '}') then
                    i_end = i - 1
                    exit do2
                endif
            enddo do2

            if (i_end <= i_start) then
                print *, "While processing lines on phase '"//trim(phase%Name)//"',"
                print *, "an error was found in a line. It is ill formed."
                stop 'ProcessLine - ModuleWaterProperties - ERR-010'
            endif

            process_line = .true.
        endif
        
        if (process_line) then            
            if (i_start > 3) then
                phase%Lines(line_index)%LineStart = trim(trimmed_line(1:i_start-3))
            else
                phase%Lines(line_index)%LineStart = ""
            endif
            if (i_end < line_size - 3) then
                phase%Lines(line_index)%LineEnd = trim(trimmed_line(i_end+3:line_size))
            else
                phase%Lines(line_index)%LineEnd = ""
            endif
            
            if (trimmed_line(i_start:i_end) == "cell_index") then
                phase%Lines(line_index)%CellIndex = .true.
            else
                phase%Lines(line_index)%PropertyIndex = FindPhreeqCPropertyIndex(trimmed_line(i_start:i_end))
            endif
        else
            phase%Lines(line_index)%LineStart = trimmed_line
            phase%Lines(line_index)%LineEnd = ""
        endif
        
        !----------------------------------------------------------------------
        
    end subroutine ProcessLine
  
    !--------------------------------------------------------------------------
    
    subroutine InitializePhreeqCRM (ObjID,                   &
                                    SolutionMapping,         &
                                    EquilibriumMapping,      &
                                    ExchangeMapping,         &
                                    SurfaceMapping,          &
                                    GasPhaseMapping,         &
                                    SolidSolutionMapping,    &
                                    KineticsMapping,         &
                                    stat)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjID  
        real, pointer, optional                         :: SolutionMapping(:)
        real, pointer, optional                         :: EquilibriumMapping(:)
        real, pointer, optional                         :: ExchangeMapping(:)
        real, pointer, optional                         :: SurfaceMapping(:)
        real, pointer, optional                         :: GasPhaseMapping(:)
        real, pointer, optional                         :: SolidSolutionMapping(:)
        real, pointer, optional                         :: KineticsMapping(:)
        integer, optional, intent(OUT)                  :: stat
        
        !Local-------------------------------------------------------------------
        integer                                         :: ready_
        integer                                         :: stat_
        integer                                         :: index
        type(T_Phreeqc), pointer                        :: phreeqc

        !------------------------------------------------------------------------

        stat_ = UNKNOWN_
        
        call Ready (ObjID, ready_)    

cd0 :   if (ready_ .EQ. IDLE_ERR_) then
    
            if (present(SolutionMapping)) Me%SolutionMapping => SolutionMapping
            if (present(EquilibriumMapping)) Me%EquilibriumMapping => EquilibriumMapping
            if (present(ExchangeMapping)) Me%ExchangeMapping => ExchangeMapping
            if (present(SurfaceMapping)) Me%SurfaceMapping => SurfaceMapping
            if (present(GasPhaseMapping)) Me%GasPhaseMapping => GasPhaseMapping
            if (present(SolidSolutionMapping)) Me%SolidSolutionMapping => SolidSolutionMapping
            if (present(KineticsMapping)) Me%KineticsMapping => KineticsMapping
            
            do index = 1, Me%NumberOfPhreeqcs
            
                phreeqc => Me%Phreeqcs(index)
            
                call StartPhreeqCRMInstance (phreeqc)

            enddo
                
            stat_ = SUCCESS_
        else 
            
            print *, "Was not possible to initialize the PhreeqCRM module"
            stop "InitializePhreeqCRM - ModulePhreeqCRM - ERR-001"

        end if cd0


        if (present(STAT)) STAT = stat_

        !----------------------------------------------------------------------

    end subroutine InitializePhreeqCRM
                                    
    !--------------------------------------------------------------------------
    
    subroutine StartPhreeqCRMInstance (phreeqc)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                           :: phreeqc
    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer, allocatable                            :: grid2chem(:)
        integer                                         :: i
        integer                                         :: aux_logical
        real, allocatable                               :: aux_real_array(:)

        !----------------------------------------------------------------------
    
        phreeqc%ObjPhreeqC = RM_Create (Me%NumberOfGridCells, openmp_num_threads)
        if (phreeqc%ObjPhreeqC < 0) then
            print *, "It was not possible to create a PhreeqCRM instance."
            print *, "The error returned was ", phreeqc%ObjPhreeqC
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-001"
        endif
        
        stat_call = RM_SetScreenOn (phreeqc%ObjPhreeqC, 0)
        if (stat_call < 0) then
            print *, "Was not possible to set to off the messages on screen."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-002"
        endif
        
        stat_call = RM_SetErrorHandlerMode(phreeqc%ObjPhreeqC, phreeqc%ErrorHandlerMode)
        if (stat_call < 0) then
            print *, "Was not possible to setup the error handler mode."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-003"
        endif
        
        if (phreeqc%ExcessHO) then
            aux_logical = 1
        else
            aux_logical = 0
        endif
        stat_call = RM_SetComponentH2O(phreeqc%ObjPhreeqC, aux_logical)
        if (stat_call < 0) then
            print *, "Was not possible to select the excess H and O mode."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-004"
        endif
        
        stat_call = RM_SetRebalanceFraction(phreeqc%ObjPhreeqC, 0.0d0)
        if (stat_call < 0) then
            print *, "Was not possible to set the rebalance fraction."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-005"
        endif
        
        stat_call = RM_SetRebalanceByCell(phreeqc%ObjPhreeqC, 1)
        if (stat_call < 0) then
            print *, "Was not possible to activate the rebalance by cell."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-006"
        endif
        
        if (phreeqc%UseSolutionDensityVolume) then
            aux_logical = 1
        else
            aux_logical = 0
        endif
        stat_call = RM_UseSolutionDensityVolume(phreeqc%ObjPhreeqC, aux_logical)
        if (stat_call < 0) then
            print *, "Was not possible to set the use of the solution density volume."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-007"
        endif
        
        stat_call = RM_SetPartitionUZSolids(phreeqc%ObjPhreeqC, 0)        
        if (stat_call < 0) then
            print *, "Was not possible to set the partition UZ solids."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-008"
        endif
        
        allocate (grid2chem(phreeqc%NumberOfCells))
        do i = 1, phreeqc%NumberOfCells
            grid2chem (i) = i - 1
        enddo        
        stat_call = RM_Createmapping (phreeqc%ObjPhreeqC, grid2chem)
        if (stat_call < 0) then
            print *, "Was not possible to screate the grid 2 chem mapping."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-009"
        endif
        
        stat_call = RM_SetUnitsSolution (phreeqc%ObjPhreeqC, Me%UnitsSolution)        
        if (stat_call < 0) then
            print *, "Was not possible to set the solution units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-010"
        endif
        
        stat_call = RM_SetUnitsExchange (phreeqc%ObjPhreeqC, Me%UnitsExchange)
        if (stat_call < 0) then
            print *, "Was not possible to set the solution units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-011"
        endif
                
        stat_call = RM_SetUnitsGasPhase (phreeqc%ObjPhreeqC, Me%UnitsGasPhase)
        if (stat_call < 0) then
            print *, "Was not possible to set the gas phase units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-012"
        endif
        
        stat_call = RM_SetUnitsKinetics (phreeqc%ObjPhreeqC, Me%UnitsKinetics)
        if (stat_call < 0) then
            print *, "Was not possible to set the knetics units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-013"
        endif
        
        stat_call = RM_SetUnitsPPassemblage (phreeqc%ObjPhreeqC, Me%UnitsPPassemblage)
        if (stat_call < 0) then
            print *, "Was not possible to set the pp assemblage units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-014"
        endif
        
        stat_call = RM_SetUnitsSSassemblage (phreeqc%ObjPhreeqC, Me%UnitsSSassemblage)
        if (stat_call < 0) then
            print *, "Was not possible to set the ss assemblage units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-015"
        endif
        
        stat_call = RM_SetUnitsSurface (phreeqc%ObjPhreeqC, Me%UnitsSurface)
        if (stat_call < 0) then
            print *, "Was not possible to set the surface units."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-016"
        endif
        
        stat_call = RM_SetTimeconversion(phreeqc%ObjPhreeqC, 1.0d0)
        if (stat_call < 0) then
            print *, "An error happened when trying to set the time conversion factor."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-017"
        endif        
        
        stat_call = RM_LoadDatabase (phreeqc%ObjPhreeqC, phreeqc%Database)
        if (stat_call < 0) then
            print *, "An error happened when loading the PhreeqC database '", trim(phreeqc%Database), "'."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-018"
        endif     
        
        allocate (aux_real_array(Me%NumberOfGridCells))
        aux_real_array = 1.0d0
        stat_call = RM_SetRepresentativeVolume (phreeqc%ObjPhreeqC, aux_real_array)
        if (stat_call < 0) then
            print *, "An error happened when trying to set the representative volume to 1.0."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-019"
        endif
        
        !ToDo: Both Porosity and Saturation must change, so MOHID Land can be used also
        stat_call = RM_SetPorosity (phreeqc%ObjPhreeqC, aux_real_array)
        if (stat_call < 0) then
            print *, "An error happened when trying to set the porosity to 1.0."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-020"
        endif
        
        stat_call = RM_SetSaturation (phreeqc%ObjPhreeqC, aux_real_array)
        if (stat_call < 0) then
            print *, "An error happened when trying to set the saturation to 1.0."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-021"
        endif
        
        stat_call = RM_SetSelectedOutputOn(phreeqc%ObjPhreeqC, 1)
        if (stat_call < 0) then
            print *, "An error happened when trying to enable user outputs."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-022"
        endif
        
        stat_call = RM_SetPrintChemistryOn(phreeqc%ObjPhreeqC, 0, 0, 0)
        if (stat_call < 0) then
            print *, "An error happened when trying to enable printing."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-022"
        endif   
        
        if (.not. Me%InitialSolutions%UseThis) then
        
            stat_call = RM_RunFile(phreeqc%ObjPhreeqC, 1, 1, 1, phreeqc%InputFile)
            if (stat_call < 0) then
                print *, "An error happened when trying to run the file '", trim(phreeqc%InputFile), "'."
                stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-023"
            endif

        else
            
            stat_call = RunInitialSolutions(phreeqc)
            if (stat_call < 0) then
                print *, "An error happened when trying to run the initial solutions."
                stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-024"
            endif
            
        endif
        
        stat_call = RM_RunString(phreeqc%ObjPhreeqC, 1, 0, 1, ClearContentsStr)  ! workers, initial_phreeqc, utility
        if (stat_call < 0) then
            print *, "An error happened when trying to clean the contents."
            stop "StartPhreeqCRMInstance - ModulePhreeqCRM - ERR-025"
        endif       
        
        call InititalizePhreeqCCells (phreeqc)
        
        deallocate (aux_real_array)
        
        !----------------------------------------------------------------------
    
    end subroutine StartPhreeqCRMInstance
    
    !--------------------------------------------------------------------------
    
    function RunInitialSolutions (phreeqc) result(res)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: res
        
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer                                         :: l
        character(len=:), allocatable                   :: text
        integer                                         :: n_characters
        !----------------------------------------------------------------------
        
        res = 0
        n_characters = 0
    
        do l = 1, Me%InitialSolutions%NumberOfLines
            n_characters = n_characters + len(trim(Me%InitialSolutions%Lines(l)%Text))
        enddo
        
        n_characters = n_characters + Me%InitialSolutions%NumberOfLines        
        allocate (character(n_characters)::text)

        do l = 1, Me%InitialSolutions%NumberOfLines
            text = trim(text)//trim(Me%InitialSolutions%Lines(l)%Text)//';'
        enddo
        
        stat_call = RM_RunString(phreeqc%ObjPhreeqC, 1, 1, 1, text)
        if (stat_call < 0) then
            res = stat_call
            return
        endif        
        
        !----------------------------------------------------------------------
    
    end function RunInitialSolutions
    
    !--------------------------------------------------------------------------
    
    subroutine InititalizePhreeqCCells (phreeqc)
        
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer, allocatable                            :: ic1(:,:)
        integer, allocatable                            :: ic2(:,:)
        real, allocatable                               :: f1(:,:)
        integer                                         :: i
        
        !----------------------------------------------------------------------
        allocate (ic1(phreeqc%NumberOfCells,7))
        allocate (ic2(phreeqc%NumberOfCells,7))
        allocate (f1(phreeqc%NumberOfCells,7))
        
        ic1 = -1        
        ic2 = -1          
        f1 = 1.0
          
        do i = 1, phreeqc%NumberOfCells
            if(associated(Me%SolutionMapping)) then
                ic1(i,1) = Me%SolutionMapping(i)
            else
                ic1(i,1) = phreeqc%Initialization(1, 1) ! Solution
            endif
            if(associated(Me%EquilibriumMapping)) then
                ic1(i,1) = Me%EquilibriumMapping(i)
            else
                ic1(i,2) = phreeqc%Initialization(2, 1) ! Equilibrium phases
            endif
            if(associated(Me%ExchangeMapping)) then
                ic1(i,1) = Me%ExchangeMapping(i)
            else
                ic1(i,3) = phreeqc%Initialization(3, 1) ! Exchange
            endif
            if(associated(Me%SurfaceMapping)) then
                ic1(i,1) = Me%SurfaceMapping(i)
            else
                ic1(i,4) = phreeqc%Initialization(4, 1) ! Surface
            endif
            if(associated(Me%GasPhaseMapping)) then
                ic1(i,1) = Me%GasPhaseMapping(i)
            else
                ic1(i,5) = phreeqc%Initialization(5, 1) ! Gas phase
            endif
            if(associated(Me%SolidSolutionMapping)) then
                ic1(i,1) = Me%SolidSolutionMapping(i)
            else            
                ic1(i,6) = phreeqc%Initialization(6, 1) ! Solid solutions
            endif
            if(associated(Me%KineticsMapping)) then
                ic1(i,7) = Me%KineticsMapping(i)
            else
                ic1(i,7) = phreeqc%Initialization(7, 1) ! Kinetics
            endif
    
            ic2(i,1) = phreeqc%Initialization(1, 2) ! Solution
            ic2(i,2) = phreeqc%Initialization(2, 2) ! Equilibrium phases
            ic2(i,3) = phreeqc%Initialization(3, 2) ! Exchange
            ic2(i,4) = phreeqc%Initialization(4, 2) ! Surface
            ic2(i,5) = phreeqc%Initialization(5, 2) ! Gas phase
            ic2(i,6) = phreeqc%Initialization(6, 2) ! Solid solutions
            ic2(i,7) = phreeqc%Initialization(7, 2) ! Kinetics
            
            f1(i,1) = phreeqc%MixFraction(1) ! Solution
            f1(i,2) = phreeqc%MixFraction(2) ! Equilibrium phases
            f1(i,3) = phreeqc%MixFraction(3) ! Exchange
            f1(i,4) = phreeqc%MixFraction(4) ! Surface
            f1(i,5) = phreeqc%MixFraction(5) ! Gas phase
            f1(i,6) = phreeqc%MixFraction(6) ! Solid solutions
            f1(i,7) = phreeqc%MixFraction(7) ! Kinetics
        enddo
        
        stat_call = RM_InitialPhreeqc2Module(phreeqc%ObjPhreeqC, ic1, ic2, f1)
        if (stat_call < 0) then
            print *, "An error happened when trying to initialize PhreeqC cells."
            stop "InititalizePhreeqCCells - ModulePhreeqCRM - ERR-001"
        endif
        
        !----------------------------------------------------------------------
    
    end subroutine InititalizePhreeqCCells
    
    !--------------------------------------------------------------------------
    
    function FindPhreeqCProperty (phreeqc_name, index, stat) result (property)
    
        !Arguments-------------------------------------------------------------
        character(*)                                    :: phreeqc_name
        integer                                         :: index
        integer                                         :: stat
        type(T_Property), pointer                       :: property
        
        !Local-----------------------------------------------------------------
        integer                                         :: i
        
        !----------------------------------------------------------------------
        
        stat = UNKNOWN_
        index = -1
        
do1:    do i = 1, Me%NumberOfProperties
            
            property => Me%Properties(i)
            
            if (trim(property%PhreeqCName) == trim(phreeqc_name)) then
                stat = SUCCESS_
                index = i
                exit do1
            endif
            
        enddo do1
        
        if (stat == UNKNOWN_) property => null()
        
        !----------------------------------------------------------------------
        
    end function FindPhreeqCProperty
    
    !--------------------------------------------------------------------------
    
    function FindPhreeqCPropertyIndex (phreeqc_name) result (index)
    
        !Arguments-------------------------------------------------------------
        character(*)                                    :: phreeqc_name
        integer                                         :: index
        
        !Local-----------------------------------------------------------------
        integer                                         :: i
        type(T_Property), pointer                       :: property
        
        !----------------------------------------------------------------------
        
        index = -1
        
do1:    do i = 1, Me%NumberOfProperties
            
            property => Me%Properties(i)
            
            if (trim(property%PhreeqCName) == trim(phreeqc_name)) then
                index = i
                exit do1
            endif
            
        enddo do1
        
        !----------------------------------------------------------------------
        
    end function FindPhreeqCPropertyIndex
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructPropertyList
    
        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer                                         :: index
        logical                                         :: found
        integer                                         :: client_number
        integer                                         :: number_of_outputs
        
        !----------------------------------------------------------------------
        
        print *, "Constructing ModulePhreeqC properties..."
        print *, ""
        
        number_of_outputs = 0
        
        call GetNumberOfBlocks (Me%ObjEnterData, '<beginproperty>', '<endproperty>',    &
                                FromFile_,                                              &
                                Me%NumberOfProperties,                                  &
                                stat = stat_call)        
        if (stat_call /= SUCCESS_) then
            print *, "An unspecified error happened when retrieving the number of properties"
            print *, "in the input file."
            print *, "The error returned was: ", stat_call
            stop "ConstructPropertyList - ModulePhreeqCRM - ERR-001"
        endif
    
        if (Me%NumberOfProperties > 0) then
        
            allocate (Me%Properties (Me%NumberOfProperties))
            allocate (Me%PropertiesID (Me%NumberOfProperties))
           
            do index = 1, Me%NumberOfProperties
            
                call ExtractBlockFromBuffer (Me%ObjEnterData,                                           &
                                             ClientNumber    = client_number,                           &
                                             block_begin     = '<beginproperty>',                       &
                                             block_end       = '<endproperty>',                         &
                                             BlockFound      = found,                                   &
                                             STAT            = stat_call)
                if (stat_call /= SUCCESS_) then
                    print *, "An unspecified error happened when retrieving the property block number '", index, "'."
                    stop "ConstructPropertyList - ModulePhreeqCRM - ERR-002"
                elseif (.not. found) then
                    print *, "A property block was expected, but it was not found."
                    stop "ConstructPropertyList - ModulePhreeqCRM - ERR-003"
                endif
                
                call ConstructProperty (index, number_of_outputs, stat_call) 
                if (stat_call /= SUCCESS_) then
                    print *, "It was not possible to construct the property number '", index, "'."
                    stop "ConstructPhreeqcsList - ModulePhreeqCRM - ERR-004"                    
                endif
                
            enddo
            
            call Block_Unlock(Me%ObjEnterData, client_number, STAT = stat_call)
        
            call RewindBuffer(Me%ObjEnterData, STAT = stat_call)
            if (stat_call .NE. SUCCESS_) then
                print *, "An unspecified error happened when rewinding the input data file buffer."
                stop 'ConstructPropertyList - ModulePhreeqCRM - ERR-005'
            endif
            
        else
            
            print *, "It is necessary to have at least ONE 'property' block defined in the MOHID PhreeqC input data file."
            stop 'ConstructPropertyList - ModulePhreeqCRM - ERR-006'
                
        endif   
        
        Me%NumberOfComponents = Me%NumberOfProperties - number_of_outputs
        
        !----------------------------------------------------------------------
    
    end subroutine ConstructPropertyList
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructProperty (index, number_of_outputs, stat)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: index
        integer                                         :: number_of_outputs
        integer                                         :: stat
                                                    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_call
        integer                                         :: iflag
        
        !----------------------------------------------------------------------
        
        stat = SUCCESS_
        
        call ConstructPropertyID (Me%Properties(index)%ID, Me%ObjEnterData, FromBlock)
        
        print *, "Index: ", index
        print *, "Property '", trim(Me%Properties(index)%ID%Name), "'"
        print *, "ID: ", Me%Properties(index)%ID%IDNumber
        print *, ""
        
        Me%PropertiesID (index) = Me%Properties(index)%ID%IDNumber
        
        call GetData (Me%Properties(index)%PhreeqCName,             &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'PHREEQC_NAME',               &
                      SearchType    = FromBlock,                    &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_ .or. iflag <= 0) then
            print *, "Was not possible to read the PHREEQC_NAME keyword."
            stat = UNKNOWN_
        endif
        
        call GetData (Me%Properties(index)%IsOutputOnly,            &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'FOR_OUTPUT_ONLY',            &
                      SearchType    = FromBlock,                    &
                      Default       = .false.,                      &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_) then
            print *, "Was not possible to read the IS_OUTPUT_ONLY keyword."
            stat = UNKNOWN_
        endif
        
        call GetData (Me%Properties(index)%GFW,                     &
                      Me%ObjEnterData, iflag,                       &
                      Keyword       = 'GFW',                        &
                      SearchType    = FromBlock,                    &
                      Default       = 0.0,                          &
                      ClientModule  = 'ModulePhreeqCRM',            &
                      STAT          = stat_call)
        if (stat_call .NE. SUCCESS_) then
            print *, "Was not possible to read the GFW keyword."
            stat = UNKNOWN_
        endif
        
        if (Me%Properties(index)%IsOutputOnly) number_of_outputs = number_of_outputs + 1
        
        !----------------------------------------------------------------------
        
    end subroutine ConstructProperty
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine GetPhreeqCIndex (obj_id, property_id, index, stat)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: obj_id
        integer                                         :: property_id
        integer, intent(out)                            :: index
        integer                                         :: stat

        !Local-----------------------------------------------------------------
        integer                                         :: ready_
        integer                                         :: i

        !----------------------------------------------------------------------

        stat = UNKNOWN_

        call Ready(obj_id, ready_)    
        
        if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

do1:        do i = 1, Me%NumberOfProperties
            
                if (property_id == Me%Properties(i)%ID%IDNumber) then
                    
                    stat = SUCCESS_
                    
                    index = i
                    
                endif
            
            enddo do1

            stat = SUCCESS_
            
        else 
        
            stat = ready_
        
        end if

        !----------------------------------------------------------------------
    
    end subroutine GetPhreeqCIndex
    
    !--------------------------------------------------------------------------
    
    subroutine GetPhreeqCDT (obj_id, dt, stat)

        !Arguments-------------------------------------------------------------
        integer                                         :: obj_id
        real                                            :: dt
        integer                                         :: stat
        
        !Local-----------------------------------------------------------------
        integer :: ready_              
        
        !----------------------------------------------------------------------

        stat = UNKNOWN_

        call Ready(obj_id, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            dt = Me%PhreeqCDT

            stat = SUCCESS_
            
        else 
    
            stat = ready_
            
        end if cd1

        !----------------------------------------------------------------------

    end subroutine GetPhreeqCDT 

    !--------------------------------------------------------------------------
    
    subroutine GetPhreeqCSize (obj_id, lb, ub, stat)

        !Arguments-------------------------------------------------------------
        integer                                         :: obj_id
        integer                                         :: lb
        integer                                         :: ub
        integer                                         :: stat
        
        !Local-----------------------------------------------------------------
        integer :: ready_              
        
        !----------------------------------------------------------------------

        stat = UNKNOWN_

        call Ready(obj_id, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            lb = 1
            ub = Me%NumberOfProperties

            stat = SUCCESS_
            
        else 
    
            stat = ready_
            
        end if cd1

        !----------------------------------------------------------------------

    end subroutine GetPhreeqCSize 
    
    !--------------------------------------------------------------------------
    
    subroutine GetPhreeqCPropertiesList (obj_id, list, stat)

        !Arguments-------------------------------------------------------------
        integer                                         :: obj_id
        integer, dimension(:), pointer                  :: list
        integer, optional, intent(OUT)                  :: stat

        !Local-----------------------------------------------------------------
        integer                                         :: ready_

        !----------------------------------------------------------------------

        stat = UNKNOWN_

        call Ready(obj_id, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mPhreeqC_, Me%InstanceID)

            list =>  Me%PropertiesID

            stat = SUCCESS_
            
        else 
    
            stat = ready_
            
        end if if1

        !----------------------------------------------------------------------

    end subroutine GetPhreeqCPropertiesList

    !--------------------------------------------------------------------------
          
    subroutine UnGetPhreeqC (obj_id, array, stat)

        !Arguments-------------------------------------------------------------
        integer                                         :: obj_id
        integer, dimension(:), pointer                  :: array
        integer, intent(OUT)                            :: stat

        !Local-----------------------------------------------------------------
        integer                                         :: ready_

        !----------------------------------------------------------------------

        stat = UNKNOWN_

        call Ready(obj_id, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(array)
            call Read_Unlock(mPhreeqC_, Me%InstanceID, "UnGetPhreeqC")

            stat = SUCCESS_
            
        else
            
            stat = ready_
            
        end if

        !----------------------------------------------------------------------
            
    end subroutine UnGetPhreeqC
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    subroutine ModifyPhreeqCRM (obj_id,                 &
                                temperature,            &
                                pressure,               &
                                porosity,               &
                                saturation,             &
                                mass,                   &                                
                                openpoints,             &
                                density,                &
                                points,                 &
                                is_starting,            &
                                stat)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: obj_id
        real, pointer                                   :: temperature(:)
        real, pointer                                   :: pressure(:)
        real, pointer                                   :: porosity(:)
        real, pointer                                   :: saturation(:)
        real, pointer                                   :: mass(:,:)
        integer, pointer                                :: openpoints(:)
        real, pointer, optional                         :: density(:)
        type(T_Point), pointer                          :: points(:)
        logical                                         :: is_starting
        integer, intent(OUT), optional                  :: stat
    
        !Local-----------------------------------------------------------------
        integer                                         :: stat_, ready_
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: phreeqc_i        
        integer                                         :: i
        integer                                         :: ncomps
        integer                                         :: property_i
        integer                                         :: first_point
        integer                                         :: last_point
        !$ integer                                      :: CHUNK
        
        !----------------------------------------------------------------------
       
        stat_ = UNKNOWN_
                      
        call Ready(obj_id, ready_)
        if (ready_ .EQ. IDLE_ERR_) then

            first_point = lbound(openpoints, 1)
            last_point = ubound(openpoints, 1)
            
            do phreeqc_i = 1, Me%NumberOfPhreeqcs
                
                phreeqc => Me%Phreeqcs(phreeqc_i)                
                
                ncomps = RM_FindComponents(phreeqc%ObjPhreeqC) 
                if (ncomps /= phreeqc%NumberOfComponents) then                    
                    call UpdatePhreeqcComponentsLookupTable (phreeqc, ncomps)
                endif

                do i = first_point, last_point
                if (openpoints(i) == 0) then                    
                    !This will ensure that PhreeqC will not run for this cell
                    saturation (i) = 0.0
                endif    
                enddo
                
                if (.not. is_starting) then
                    
                    !$ CHUNK = CHUNK_I(1, property_i)
                    !$OMP PARALLEL PRIVATE(i,property_i)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do property_i = 1, (phreeqc%NumberOfComponents)
                    do i = first_point, last_point
                        phreeqc%Concentrations(i,property_i) = mass(phreeqc%ComponentsToProperties(property_i),i)
                    enddo
                    enddo
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL
                    
                endif    
                
                if (is_starting) then
                    stat_ = RM_SetTime(phreeqc%ObjPhreeqC, 0.0)
                    if (stat_ < 0) then
                        print *, "An error happened when trying to set the time."
                        stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-001"
                    endif
                
                    stat_ = RM_SetTimeStep(phreeqc%ObjPhreeqC, 0.0)
                    if (stat_ < 0) then
                        print *, "An error happened when trying to set the time step."
                        stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-010"
                    endif
                else
                    stat_ = RM_SetTime(phreeqc%ObjPhreeqC, Me%ElapsedTime)
                    if (stat_ < 0) then
                        print *, "An error happened when trying to set the time."
                        stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-020"
                    endif
                
                    stat_ = RM_SetTimeStep(phreeqc%ObjPhreeqC, Me%PhreeqCDT)
                    if (stat_ < 0) then
                        print *, "An error happened when trying to set the time step."
                        stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-030"
                    endif
                endif
                       
                stat_ = RM_SetTemperature(phreeqc%ObjPhreeqC, temperature)
                if (stat_ < 0) then
                    print *, "An error happened when trying to set the temperature."
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-040"
                endif 
                
                stat_ = RM_SetSaturation(phreeqc%ObjPhreeqC, saturation)
                if (stat_ < 0) then
                    print *, "An error happened when trying to set the saturation."
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-050"
                endif
                
                stat_ = RM_SetPorosity(phreeqc%ObjPhreeqC, porosity)
                if (stat_ < 0) then
                    print *, "An error happened when trying to set the porosity."
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-060"
                endif
                
                stat_ = RM_SetPressure(phreeqc%ObjPhreeqC, pressure)
                if (stat_ < 0) then
                    print *, "An error happened when trying to set the pressure."
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-070"
                endif
                
                if (phreeqc%UseSolutionDensityVolume) then
                stat_ = RM_SetDensity(phreeqc%ObjPhreeqC, density)
                if (stat_ < 0) then
                    print *, "An error happened when trying to set the density."
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-080"
                endif
                endif

                if (.not. is_starting) then
                    stat_ = RM_SetConcentrations(phreeqc%ObjPhreeqC, phreeqc%Concentrations)
                    if (stat_ < 0) then
                        print *, "An error happened when trying to set the concentrations."
                        stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-090"
                    endif
                    
                    do i = 1, phreeqc%NumberOfPhases
                        call UpdatePhase(phreeqc, i, mass, openpoints)
                    enddo
                endif
                
                stat_ = RM_RunCells(phreeqc%ObjPhreeqC)
                if (stat_ < 0) then
                    if (Me%LogOnError) then
                    if (.not. is_starting) then
                        call WriteLogFile(phreeqc, first_point, last_point, points)
                    endif
                    endif
                    
                    print *, "An error happened when trying to run PhreeqC."                    
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-100"
                endif                
                
                stat_ = RM_GetConcentrations(phreeqc%ObjPhreeqC, phreeqc%Concentrations)
                if (stat_ < 0) then
                    print *, "An error happened when trying to get concentrations."
                    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-110"
                endif
                
                !stat_ = RM_GetDensity(phreeqc%ObjPhreeqC, density)
                !if (stat_ < 0) then
                !    print *, "An error happened when trying to get the density."
                !    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-011"
                !endif
            
                !stat_ = RM_GetSolutionVolume(phreeqc%ObjPhreeqC, solution_volume)
                !if (stat_ < 0) then
                !    print *, "An error happened when trying to get the solution volume."
                !    stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-012"
                !endif 
                
                !$ CHUNK = CHUNK_I(1, property_i)
                !$OMP PARALLEL PRIVATE(i,property_i)
                !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do property_i = 1, phreeqc%NumberOfComponents
                do i = first_point, last_point
                    if (openpoints(i) == 0) then
                        mass(phreeqc%ComponentsToProperties(property_i),i) = 0.0d0
                    else
                        mass(phreeqc%ComponentsToProperties(property_i),i) = phreeqc%Concentrations(i,property_i)
                    endif
                enddo
                enddo
                !$OMP END DO NOWAIT
                !$OMP END PARALLEL
                
                if (phreeqc%UserOutput >= 0) then 
                    if (phreeqc%NumberOfOutputs <= 0) then
                        call UpdatePhreeqcOututsLookupTable(phreeqc)
                    endif
    
                    stat_ = RM_GetSelectedOutput (phreeqc%ObjPhreeqC, phreeqc%Outputs)
                    if (stat_ < 0) then
                        print *, "An error happened when trying to get the user output."
                        stop "ModifyPhreeqCRM - ModulePhreeqCRM - ERR-120"
                    endif
                    
                    !$ CHUNK = CHUNK_I(1, property_i)
                    !$OMP PARALLEL PRIVATE(i,property_i)
                    !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                    do property_i = 1, phreeqc%NumberOfOutputs
                    do i = 1, phreeqc%NumberOfCells
                        if (openpoints(i) == 0) then
                            mass (phreeqc%OutputsToProperties(property_i),i) = 0.0
                        else                            
                            mass (phreeqc%OutputsToProperties(property_i),i) = phreeqc%Outputs(i,property_i)
                        endif
                    enddo
                    enddo
                    !$OMP END DO NOWAIT
                    !$OMP END PARALLEL
                    
                endif

            enddo
            
            stat_ = SUCCESS_
            
        else               
            
            stat_ = ready_
            
        end if

        Me%ElapsedTime = Me%ElapsedTime + Me%PhreeqCDT
        
        if (present(stat)) stat = stat_
        
        !----------------------------------------------------------------------        
    
    end subroutine ModifyPhreeqCRM
    
    !--------------------------------------------------------------------------
                                
    subroutine UpdatePhase (phreeqc, index, mass, openpoints)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: index
        real, pointer                                   :: mass(:,:)
        integer, pointer                                :: openpoints(:)
        
        !Local-----------------------------------------------------------------
        type(T_Phase), pointer                          :: phase
        type(T_PhaseLine), pointer                      :: line
        integer                                         :: l
        integer                                         :: w
        integer                                         :: n_workers
        integer                                         :: id
        integer, allocatable                            :: start_cell(:)
        integer, allocatable                            :: end_cell(:)
        integer                                         :: sc
        integer                                         :: ec
        integer                                         :: c
        integer                                         :: stat_call
        real                                            :: mass_mol
        real                                            :: gfw
        character(len=30)                               :: text
        integer                                         :: cell
        !----------------------------------------------------------------------
        phase => phreeqc%Phases(index)
                
        n_workers = RM_GetThreadCount(phreeqc%ObjPhreeqC)        
        
        allocate (start_cell(n_workers))
        allocate (end_cell(n_workers))        
        stat_call = RM_GetStartCell(phreeqc%ObjPhreeqC, start_cell)
        stat_call = RM_GetEndCell(phreeqc%ObjPhreeqC, end_cell)
        
        do w = 1, n_workers
            
            sc = start_cell(w)
            ec = end_cell(w)
            
            id = RM_GetIPhreeqcId(phreeqc%ObjPhreeqC, w - 1)
            
            do c = sc, ec
                cell = c
                if (openpoints(c+1) == 1) then
                if (phase%NumberOfLines > 0) then                    
                    do l = 1, phase%NumberOfLines        
                        line => phase%Lines(l)
                
                        if (line%PropertyIndex > -1) then
                            gfw = Me%Properties(line%PropertyIndex)%GFW
                            if (gfw <= 0) then
                                print *, 'Property '//trim(Me%Properties(line%PropertyIndex)%ID%Name)//' has a GFW <= 0.0'
                                print *, 'Set the correct GFW or deactivate the phase where it is used, please.'
                                stop 'UpdatePhase - ModulePhreeqCRM - ERR-010'
                            endif
                            !mol = (mg/L) * mol/g * g/mg * L
                            mass_mol = mass(line%PropertyIndex, c+1) / gfw / 1000.0 * 1.0
                            write (text, '(F25.15)') mass_mol 
                            stat_call = AccumulateLine (id, trim(line%LineStart)//" "//trim(text)//" "//trim(line%LineEnd))
                        elseif (line%CellIndex) then
                            write (text, '(I10)') cell
                            !print *, trim(line%LineStart)//" "//trim(text)//" "//trim(line%LineEnd)
                            stat_call = AccumulateLine (id, trim(line%LineStart)//" "//trim(text)//" "//trim(line%LineEnd))
                        else
                            !print *, line%LineStart
                            stat_call = AccumulateLine (id, line%LineStart)
                        endif
                    enddo
                    
                    stat_call = RunAccumulated (id)
                endif
                endif
            enddo
        enddo
        
        !----------------------------------------------------------------------
    
    end subroutine UpdatePhase
                                
    !--------------------------------------------------------------------------
                       
    subroutine UpdatePhreeqcComponentsLookupTable (phreeqc, ncomps)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: ncomps
        
        !Local-----------------------------------------------------------------
        character(StringLength)                         :: component_name
        integer                                         :: i
        integer                                         :: index
        type(T_Property), pointer                       :: property
        integer                                         :: stat_call
        
        !----------------------------------------------------------------------
        if (associated(phreeqc%Concentrations)) then
            deallocate (phreeqc%Concentrations)
        endif
        
        do i = 1, ncomps
            stat_call = RM_GetComponent(phreeqc%ObjPhreeqC, i, component_name)
            if (stat_call < 0) then
                print *, "An error happened when trying to get Phreeqc component names."
                stop "UpdatePhreeqcComponentsLookupTable - ModulePhreeqCRM - ERR-010"
            endif
            
            property => FindPhreeqCProperty(component_name, index, stat_call)
            if (stat_call /= SUCCESS_) then
                print *, "The Phreeqc component '"//trim(component_name)//"' does not have"
                print *, "an associated MOHID property."
                stop "UpdatePhreeqcComponentsLookupTable - ModulePhreeqCRM - ERR-010"
            else
                phreeqc%ComponentsToProperties(i) = index
            endif            
        enddo
        
        phreeqc%NumberOfComponents = ncomps
                
        allocate (phreeqc%Concentrations(phreeqc%NumberOfCells, phreeqc%NumberOfComponents))
        
        !----------------------------------------------------------------------
    
    end subroutine UpdatePhreeqcComponentsLookupTable
    
    !--------------------------------------------------------------------------
                       
    subroutine UpdatePhreeqcOututsLookupTable (phreeqc)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        
        !Local-----------------------------------------------------------------
        integer                                         :: i
        integer                                         :: index
        type(T_Property), pointer                       :: property
        integer                                         :: stat_call
        character(StringLength)                         :: output_header        
        
        !----------------------------------------------------------------------
        if (associated(phreeqc%Outputs)) then
            deallocate (phreeqc%Outputs)
        endif
        
        phreeqc%NumberOfOutputs = RM_GetSelectedOutputcolumncount(phreeqc%ObjPhreeqC)
        if (phreeqc%NumberOfOutputs <= 0) then
            print *, "The number of user output columns for phreeqc should be > 0"
            stop "UpdatePhreeqcOututsLookupTable - ModulePhreeqCRM - ERR-001"
        endif
        
        do i = 1, phreeqc%NumberOfOutputs
            stat_call = RM_GetSelectedOutputheading(phreeqc%ObjPhreeqC, i, output_header)
            if (stat_call < 0) then
                print *, "An error happened when trying to get the user output header."
                stop "UpdatePhreeqcOututsLookupTable - ModulePhreeqCRM - ERR-005"
            endif 
            
            property => FindPhreeqCProperty(output_header, index, stat_call)
            if (stat_call /= SUCCESS_) then
                print *, "The Phreeqc component '"//trim(output_header)//"' does not have"
                print *, "an associated MOHID property."
                stop "UpdatePhreeqcOututsLookupTable - ModulePhreeqCRM - ERR-010"
            else
                phreeqc%OutputsToProperties(i) = index
            endif            
        enddo
                
        allocate (phreeqc%Outputs(phreeqc%NumberOfCells, phreeqc%NumberOfOutputs))
        
        !----------------------------------------------------------------------
    
    end subroutine UpdatePhreeqcOututsLookupTable
    
    !--------------------------------------------------------------------------
    
    subroutine WriteLogFile (phreeqc, first_point, last_point, points)
    
        !Arguments-------------------------------------------------------------
        type(T_Phreeqc), pointer                        :: phreeqc
        integer                                         :: first_point
        integer                                         :: last_point
        type(T_Point), pointer                          :: points(:)
        
        !Local-----------------------------------------------------------------
        integer                                         :: status
        integer                                         :: UnitOutput
        integer                                         :: ncomps
        integer                                         :: i
        integer                                         :: property_i
        character(StringLength)                         :: component_name
        character(StringLength)                         :: token
        character(2048)                                 :: text, aux
        
        !----------------------------------------------------------------------
        call UnitsManager (UnitOutput, OPEN_FILE)      
        open(UNIT = UnitOutput, FILE = Me%LogErrorFile, STATUS  = "UNKNOWN", IOSTAT  = status)
        if (status == SUCCESS_) then
        
            write(UnitOutput, *)"MOHID PhreeqCRM Error Log File"
            write(UnitOutput, *)"This file contains the concentrations sent to PhreeqCRM."
            write(UnitOutput, *)"Phreeqc: "//trim(phreeqc%ID%Name)
            
                        
            ncomps = RM_FindComponents(phreeqc%ObjPhreeqC)            
            
            write(UnitOutput, *)"Number of Components: ", ncomps
            write(UnitOutput, *)"--------------------------------------------------------"
            
            text = "Cell I J K"
            do i = 1, ncomps
              status = RM_GetComponent(phreeqc%ObjPhreeqC, i, component_name)
              aux = trim(text)//" "//trim(component_name)
              text = aux
            enddo               
                       
            write(UnitOutput, '(A)')trim(text)
            
            do i = first_point, last_point
                write(text, '(I8,3I5)') i, points(i)%I, points(i)%J, points(i)%K 
                do property_i = 1, (phreeqc%NumberOfComponents)
                    write(token, '(E15.7)') phreeqc%Concentrations(i-1,property_i-1)
                    aux = trim(text)//" "//trim(token)
                    text = aux
                enddo
                write(UnitOutput, '(A)') trim(text)
            enddo
               
            call UnitsManager (UnitOutput, CLOSE_FILE)
            
        endif        
        
        !----------------------------------------------------------------------
        
    end subroutine WriteLogFile
                                
    !--------------------------------------------------------------------------
    
    !subroutine InitializeEquilibrium (phreeqc, mass)
    !
    !    !Arguments-------------------------------------------------------------
    !    type(T_Phreeqc), pointer                           :: phreeqc
    !    real, pointer                                   :: mass(:,:)
    !    
    !    !Local-----------------------------------------------------------------
    !    integer                                         :: stat_call
    !    real, allocatable                               :: concentrations(:,:)
    !    integer                                         :: i, j
    !    
    !    !----------------------------------------------------------------------
    !            
    !    allocate (concentrations(lbound(mass,2)-1:ubound(mass,2)-1,lbound(mass,1)-1:ubound(mass,1)-1))
    !    
    !    stat_call = RM_SetTime (phreeqc%ObjPhreeqC, 0.0)
    !    if (stat_call < 0) then
    !        print *, "An error happened when trying to initialize PhreeqC cells."
    !        stop "InititalizePhreeqCCells - ModulePhreeqCRM - ERR-001"
    !    endif
    !    
    !    stat_call = RM_SetTimeStep (phreeqc%ObjPhreeqC, 0.0)
    !    if (stat_call < 0) then
    !        print *, "An error happened when trying to initialize PhreeqC cells."
    !        stop "InititalizePhreeqCCells - ModulePhreeqCRM - ERR-001"
    !    endif
    !    
    !    stat_call = RM_RunCells (phreeqc%ObjPhreeqC) 
    !    if (stat_call < 0) then
    !        print *, "An error happened when trying to initialize PhreeqC cells."
    !        stop "InititalizePhreeqCCells - ModulePhreeqCRM - ERR-001"
    !    endif
    !    
    !    stat_call = RM_GetConcentrations (phreeqc%ObjPhreeqC, concentrations)
    !    if (stat_call < 0) then
    !        print *, "An error happened when trying to initialize PhreeqC cells."
    !        stop "InititalizePhreeqCCells - ModulePhreeqCRM - ERR-001"
    !    endif
    !    
    !    do i = 0, ubound(concentrations, 2)
    !    do j = 0, ubound(concentrations, 1)
    !        mass(i+1,j+1) = concentrations(j,i)
    !    enddo
    !    enddo
    !    
    !    !----------------------------------------------------------------------    
    !
    !end subroutine InitializeEquilibrium                                
                                
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillPhreeqCRM (ObjID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: stat_
        integer                             :: nUsers
        integer                             :: i
        integer                             :: stat_call

        !------------------------------------------------------------------------

        stat_ = UNKNOWN_

        call Ready(ObjID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPhreeqC_,  Me%InstanceID)

            if (nUsers == 0) then

                do i = 1, Me%NumberOfPhreeqcs
                    
                    stat_call = RM_Destroy (Me%Phreeqcs(i)%ObjPhreeqC)
                    if (stat_call < 0) then
                        print *, "WARNING: Was not possible to destroy instance"
                        print *, "of PhreeqC for Phreeqc '", trim(Me%Phreeqcs(i)%ID%Name), "'."
                    endif
                    
                enddo
                
                deallocate(Me%Properties)                                
                deallocate(Me%Phreeqcs)
                
                call DeallocateInstance ()

                ObjID = 0
                stat_ = SUCCESS_

            end if
            
        else 
    
            stat_ = ready_
            
        end if cd1

        if (present(STAT)) STAT = stat_

        !----------------------------------------------------------------------

    end subroutine KillPhreeqCRM

    !--------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_PhreeqCRM), pointer                     :: AuxObj
        type (T_PhreeqCRM), pointer                     :: PreviousObj
        
        !----------------------------------------------------------------------

        if (Me%InstanceID == FirstObjPhreeqCRM%InstanceID) then
            
            FirstObjPhreeqCRM => FirstObjPhreeqCRM%Next
            
        else
            
            PreviousObj => FirstObjPhreeqCRM
            AuxObj      => FirstObjPhreeqCRM%Next
            
            do while (AuxObj%InstanceID /= Me%InstanceID)
                
                PreviousObj => AuxObj
                AuxObj      => AuxObj%Next
                
            enddo

            PreviousObj%Next => AuxObj%Next

        endif

        deallocate (Me)
        nullify    (Me) 
        
        !----------------------------------------------------------------------
            
    end subroutine DeallocateInstance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Ready (ObjID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjID > 0) then
    
            call LocateObjPhreeqCRM (ObjID)
            ready_ = VerifyReadLock (mPhreeqC_, Me%InstanceID)
            
        else
    
            ready_ = OFF_ERR_
            
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjPhreeqCRM (ObjId)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjId

        !Local-----------------------------------------------------------------

        !----------------------------------------------------------------------
        
        Me => FirstObjPhreeqCRM        
        do while (associated (Me))
            
            if (Me%InstanceID == ObjId) exit
            Me => Me%Next
            
        enddo

        if (.not. associated(Me)) then
            
            print *, "The PhreeqCRM object with ID '", ObjId, "' does not exist."
            stop "LocateObjPhreeqCRM - ModulePhreeqCRM - ERR-001"
            
        endif
        
        !----------------------------------------------------------------------
        
    end subroutine LocateObjPhreeqCRM
    
    !--------------------------------------------------------------------------
    
End Module ModulePhreeqCRM
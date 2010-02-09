!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : PhreeqC
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2009
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

Module ModulePhreeqC
    use ModulePhreeqCData
    use ModuleFunctions            
    use ModuleEnterData
    use ModuleGlobalData             

    implicit none

    private 

    !Subroutines---------------------------------------------------------------
    !Constructor
    public  :: StartPhreeqC
    private ::      AllocateInstance
    private ::      InitializeInstance
    private ::      ReadInputFile
    private ::          ReadPhreeqCOptions
    private ::              ReadPhreeqCDatabase
                
    !Selector
    public  :: SetPhreeqCProperty
    private ::      SetSolutionProperty
    private ::      SetPhaseProperty
    private ::      SetExchangeProperty 
    public  :: GetPhreeqCDT 
    public  :: GetPhreeqCPropIndex    
    public  :: UngetPhreeqC           
                        
    !Modifier
    public  :: ModifyPhreeqC 
    private ::      MakeCalculations
    private ::          CalculateSolutionParameters
    private ::          ConvertInputs
    private ::          ConvertResults        
        
    !Destructor
    public  :: KillPhreeqC                                                     
    private ::      DeallocateInstance

    !Management
    private :: Ready
    private ::      LocateObjPhreeqC
        
    !Constants-----------------------------------------------------------------
    integer, parameter :: MaxStrLength  = StringLength !StringLength is defined in ModuleGlobalData 
    integer, parameter :: maxProperties = 25    
    
    integer, parameter :: SOLUTION    = 1
    integer, parameter :: PHASE       = 2
    integer, parameter :: SOLID_PHASE = 3
    integer, parameter :: GAS_PHASE   = 4
    integer, parameter :: SURFACE     = 5
    integer, parameter :: SPECIES     = 6
    integer, parameter :: EXCHANGE    = 7
    
    !Solution units (for concentration properties)
    integer, parameter :: mol_l     = 1
    integer, parameter :: mmol_l    = 2
    integer, parameter :: umol_l    = 3
    integer, parameter :: g_l       = 4
    integer, parameter :: mg_l      = 5 !default for user input data
    integer, parameter :: ug_l      = 6
    integer, parameter :: eq_l      = 7
    integer, parameter :: meq_l     = 8
    integer, parameter :: ueq_l     = 9
    integer, parameter :: mol_kgs   = 10
    integer, parameter :: mmol_kgs  = 11
    integer, parameter :: umol_kgs  = 12
    integer, parameter :: g_kgs     = 13
    integer, parameter :: mg_kgs    = 14
    integer, parameter :: ug_kgs    = 15
    integer, parameter :: eq_kgs    = 16
    integer, parameter :: meq_kgs   = 17
    integer, parameter :: ueq_kgs   = 18
    integer, parameter :: mol_kgw   = 19 !default for use with phreeqc
    integer, parameter :: mmol_kgw  = 20
    integer, parameter :: umol_kgw  = 21
    integer, parameter :: g_kgw     = 22
    integer, parameter :: mg_kgw    = 23
    integer, parameter :: ug_kgw    = 24
    integer, parameter :: eq_kgw    = 25
    integer, parameter :: meq_kgw   = 26
    integer, parameter :: ueq_kgw   = 27
    
    !Phase/Exchange units
    integer, parameter :: mol       = 29 !default for use with phreeqc
    integer, parameter :: mmol      = 30 
    integer, parameter :: umol      = 31
    
    !Phase units
    integer, parameter :: g_kgsoil  = 32
    integer, parameter :: mg_kgsoil = 33 !default for user solid phases
    integer, parameter :: ug_kgsoil = 34
            
    !Types---------------------------------------------------------------------            
    type T_PhreeqCUnits
        integer             :: MasterSpecies = mg_l
        integer             :: Species       = mg_l
        integer             :: Alkalinity    = mg_l
        integer             :: SolidPhases   = ug_kgsoil !Check to see if this is reasonable
        integer             :: Exchangers    = ug_kgsoil
    end type T_PhreeqCUnits
        
    type T_PhreeqCProperty
        integer                     :: PropertyID
        integer                     :: PhreeqCInputID
        integer                     :: PhreeqCResultID
        type(T_ChemistryParameters) :: Params
        real                        :: Mass
        real                        :: Volume
    end type T_PhreeqCProperty
        
        
    type T_PhreeqCOptions
        character(len=2048)  :: Database                  !Path for database
        character(len=2048)  :: DatabaseAux  = ''         !Path for auxiliary database
        logical              :: PrintAlways  = .false.    !For DEBUG    
        real                 :: DTSeconds    = null_real 
        real                 :: DTDay        = null_real
        real                 :: HPlusDensity              !g/L
        real                 :: WaterDensity              !g/L
        type(T_RedoxPair)    :: Redox
        type(T_PhreeqCUnits) :: Units
        integer              :: pHCharge
        integer              :: pECharge
    end type T_PhreeqCOptions
    
    type T_External
        real, pointer, dimension(:,:) :: PropertiesValues
        real, pointer, dimension(:  ) :: SolutionVolume
        real, pointer, dimension(:  ) :: SolutionTemperature
        real, pointer, dimension(:  ) :: SolutionpH
        real, pointer, dimension(:  ) :: SolutionpE
    end type T_External
    
    type T_Calculations
        real :: MassOfAllSolutes
        real :: VolumeOfAllSolutes
        real :: MassOfWater
        real :: VolumeOfWater
        real :: MassOfSolution
        real :: VolumeOfSolution
        real :: DensityOfSolution
        real :: Factor
    end type T_Calculations
    
    type T_PhreeqC
        private
        integer                                           :: InstanceID           !ID of the ModulePhreeqC instance 
        integer                                           :: PhreeqCInstanceID    !ID of the PhreeqC Object instance linked to the InstanceID instance
        integer                                           :: ObjEnterData = 0     !Instance of ModuleEnterData
        type(T_PhreeqC) , pointer                         :: Next                 !Collection of instances of ModulePhreeqC
        type(T_PhreeqCOptions)                            :: PhreeqCOptions       !Global options read from the PhreeqC Input File
        type(T_External)                                  :: Ext                  !Pointers to Water Mass, Properties Values and other required data 
        type(T_PhreeqCProperty), dimension(MaxProperties) :: Properties           !Info about each property
        integer                                           :: PropertyCount        !Number of properties
        type(T_Calculations)                              :: CalcData             !Temporary data for calculations
        
    end type T_PhreeqC


    !Global Module Variables---------------------------------------------------
    type (T_PhreeqC), pointer :: FirstObjPhreeqC
    type (T_PhreeqC), pointer :: Me      
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !INTERFACE TO PhreeqCLIB

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    interface to subroutine pm_start [ALIAS:'?pm_start@@YAXPAH0@Z'](a, b)
	    integer*4 a [REFERENCE]
	    integer*4 b [REFERENCE]
    end
         
    interface to subroutine pm_solution_redox [ALIAS: '?pm_solution_redox@@YAXPAHPADPAN120HH@Z'](a, b, c, d, e, f)
      integer*4        a [REFERENCE]
      character(Len=*) b [REFERENCE]
      real*8           c [REFERENCE]
      character(Len=*) d [REFERENCE]
      real*8           e [REFERENCE]
      integer*4        f [REFERENCE]
    end
    
    interface to subroutine pm_conc_add [ALIAS:'?pm_conc_add@@YAXPAHPADPAN0H@Z'] (a, b, c, d)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE]
      real*8           c [REFERENCE]
      integer*4        d [REFERENCE]
    end 
    
    interface to subroutine pm_conc_use [ALIAS:'?pm_conc_use@@YAXPAH0000000@Z'] (a, b, c, d, e, f, g, h)
      integer*4 a [REFERENCE] 
      integer*4 b [REFERENCE] 
      integer*4 c [REFERENCE]
      integer*4 d [REFERENCE]
      integer*4 e [REFERENCE] 
      integer*4 f [REFERENCE] 
      integer*4 g [REFERENCE]
      integer*4 h [REFERENCE]
    end 
    
    interface to subroutine pm_conc_as [ALIAS:'?pm_conc_as@@YAXPAHPAD0H@Z'] (a, b, c)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      integer*4        c [REFERENCE]
    end 

    interface to subroutine pm_conc_gfw [ALIAS:'?pm_conc_gfw@@YAXPAHPAN0@Z'] (a, b, c)
      integer*4 a [REFERENCE] 
      real*8    b [REFERENCE] 
      integer*4 c [REFERENCE]
    end 

    interface to subroutine pm_conc_phase [ALIAS:'?pm_conc_phase@@YAXPAHPADPAN0@Z'] (a, b, c, d)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      real*8           c [REFERENCE]
      integer*4        d [REFERENCE]
    end 

    interface to subroutine pm_conc_redox [ALIAS:'?pm_conc_redox@@YAXPAHPADPAN120HH@Z'] (a, b, c, d, e, f)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      real*8           c [REFERENCE]
      character(Len=*) d [REFERENCE] 
      real*8           e [REFERENCE]
      integer*4        f [REFERENCE]
    end 

    interface to subroutine pm_conc_save [ALIAS:'?pm_conc_save@@YAXPAH00PAN0@Z'] (a, b, c, d, e)
      integer*4 a [REFERENCE] 
      integer*4 b [REFERENCE] 
      integer*4 c [REFERENCE]
      real*8    d [REFERENCE]
      integer*4 e [REFERENCE]
    end 

    interface to subroutine pm_use_ppa [ALIAS:'?pm_use_ppa@@YAXPAH00@Z'] (a, b, c)
      integer*4 a [REFERENCE] 
      integer*4 b [REFERENCE] 
      integer*4 c [REFERENCE]
    end 

    interface to subroutine pm_use_exa [ALIAS:'?pm_use_exa@@YAXPAH00@Z'] (a, b, c)
      integer*4 a [REFERENCE] 
      integer*4 b [REFERENCE] 
      integer*4 c [REFERENCE]
    end 

    interface to subroutine pm_use_sa [ALIAS:'?pm_use_sa@@YAXPAH00@Z'] (a, b, c)
      integer*4 a [REFERENCE] 
      integer*4 b [REFERENCE] 
      integer*4 c [REFERENCE]
    end 

    interface to subroutine pm_ppa_pp [ALIAS:'?pm_ppa_pp@@YAXPAHPAD1PAN20000HH@Z'] (a, b, c, d, e, f, g, h, i)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      character(Len=*) c [REFERENCE]
      real*8           d [REFERENCE]
      real*8           e [REFERENCE] 
      integer*4        f [REFERENCE] 
      integer*4        g [REFERENCE]
      integer*4        h [REFERENCE]
      integer*4        i [REFERENCE]
    end 

    interface to subroutine pm_exa_exchanger [ALIAS:'?pm_exa_exchanger@@YAXPAHPAD01PAN00HH@Z'] (a, b, c, d, e, f, g)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      integer*4        c [REFERENCE]
      character(Len=*) d [REFERENCE] 
      real*8           e [REFERENCE]
      integer*4        f [REFERENCE]
      integer*4        g [REFERENCE]
    end 

    interface to subroutine pm_sa_options [ALIAS:'?pm_sa_options@@YAXPAH00PAN00@Z'] (a, b, c, d, e, f)
      integer*4        a [REFERENCE] 
      integer*4        b [REFERENCE] 
      integer*4        c [REFERENCE]
      real*8           d [REFERENCE] 
      integer*4        e [REFERENCE]
      integer*4        f [REFERENCE]
    end 

    interface to subroutine pm_sa_surface [ALIAS:'?pm_sa_surface@@YAXPAH0PAD1PAN22000HH@Z'] (a, b, c, d, e, f, g, h)
      integer*4        a [REFERENCE] 
      integer*4        b [REFERENCE] 
      character(Len=*) c [REFERENCE]
      character(Len=*) d [REFERENCE]
      real*8           e [REFERENCE] 
      real*8           f [REFERENCE] 
      real*8           g [REFERENCE]
      integer*4        h [REFERENCE]
    end 

    interface to subroutine pm_get_species_index [ALIAS:'?pm_get_species_index@@YAXPAHPAD00H@Z'] (a, b, c, d)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      integer*4        c [REFERENCE]
      integer*4        d [REFERENCE]
    end 
    
    interface to subroutine pm_run_model [ALIAS:'?pm_run_model@@YAXPAH0@Z'](a, b)
	    integer*4 a [REFERENCE]
	    integer*4 b [REFERENCE]
    end

    interface to subroutine pm_setup_model [ALIAS:'?pm_setup_model@@YAXPAH0@Z'](a, b)
	    integer*4 a [REFERENCE]
	    integer*4 b [REFERENCE]
    end

    interface to subroutine pm_read_database [ALIAS:'?pm_read_database@@YAXPAHPAD0@Z'] (a, b, c)
      integer*4        a [REFERENCE] 
      character(Len=*) b [REFERENCE] 
      integer*4        c [REFERENCE]
    end 

    interface to subroutine pm_kill [ALIAS:'?pm_kill@@YAXPAH0@Z'](a, b)
	    integer*4 a [REFERENCE]
	    integer*4 b [REFERENCE]
    end       
    
    interface to subroutine pm_set_required_data [ALIAS:'?pm_set_required_data@@YAXPAHPAN110@Z'] (a, b, c, d, e)
      integer*4        a [REFERENCE] 
      real*8           b [REFERENCE] 
      real*8           c [REFERENCE] 
      real*8           d [REFERENCE] 
      integer*4        e [REFERENCE]
    end 

    interface to subroutine pm_set_input_value [ALIAS:'?pm_set_input_value@@YAXPAH00PAN0@Z'] (a, b, c, d, e)
      integer*4        a [REFERENCE] 
      integer*4        b [REFERENCE] 
      integer*4        c [REFERENCE] 
      real*8           d [REFERENCE] 
      integer*4        e [REFERENCE]
    end 

    interface to subroutine pm_get_result_value [ALIAS:'?pm_get_result_value@@YAXPAH00PAN0@Z'] (a, b, c, d, e)
      integer*4        a [REFERENCE] 
      integer*4        b [REFERENCE] 
      integer*4        c [REFERENCE] 
      real*8           d [REFERENCE] 
      integer*4        e [REFERENCE]
    end 

    interface to subroutine pm_get_data [ALIAS:'?pm_get_data@@YAXPAHPAN110@Z'] (a, b, c, d, e)
      integer*4        a [REFERENCE] 
      real*8           b [REFERENCE] 
      real*8           c [REFERENCE] 
      real*8           d [REFERENCE] 
      integer*4        e [REFERENCE]
    end 

    interface to subroutine pm_set_ph [ALIAS: '?pm_set_ph@@YAXPAH0PAN0@Z'] (a, b, c, d)
      integer*4        a [REFERENCE]
      integer*4        b [REFERENCE]
      real*8           c [REFERENCE]
      integer*4        d [REFERENCE]
    end

    interface to subroutine pm_set_pe [ALIAS: '?pm_set_pe@@YAXPAH000@Z'] (a, b, c, d)
      integer*4        a [REFERENCE]
      integer*4        b [REFERENCE]
      real*8           c [REFERENCE]
      integer*4        d [REFERENCE]
    end

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartPhreeqC(PhreeqCID, Filename, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                        :: PhreeqCID
        character(LEN = *)             :: FileName 
        integer, optional, intent(OUT) :: STAT     

        !Local-----------------------------------------------------------------
        integer :: STAT_
        integer :: ready_         

        !----------------------------------------------------------------------
        
        STAT_ = UNKNOWN_        

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mPHREEQC_)) then
        
            nullify (FirstObjPhreeqC)
            call RegisterModule (mPHREEQC_) 
            
        endif
        
        call Ready(PhreeqCID, ready_)

cd0 :   if (ready_ .EQ. OFF_ERR_) then
            
            !Allocate a new instance of ModulePhreeqC (this module...)
            call AllocateInstance 
            
            !Initialize this instance variables
            call InitializeInstance
            
            !Create PhreeqC C++ object
            call pm_start (Me%PhreeqCInstanceID, STAT_) 
            if (STAT_ .EQ. 0) &
                 stop 'Subroutine StartPhreeqC; module ModulePhreeqC. ERR001.'
                       
            !Read ModulePhreeqC Input File  
            call ReadInputFile (FileName)
            
            !Returns ID
            PhreeqCID = Me%InstanceID            

            STAT_ = SUCCESS_
            
        else  cd0
            
            stop 'Subroutine StartPhreeqC; module ModulePhreeqC. ERR002'

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartPhreeqC
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine AllocateInstance

        !Local-----------------------------------------------------------------
        type (T_PhreeqC), pointer           :: NewObjPhreeqC
        type (T_PhreeqC), pointer           :: PreviousObjPhreeqC

        !----------------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewObjPhreeqC)
        nullify  (NewObjPhreeqC%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjPhreeqC)) then
        
            FirstObjPhreeqC    => NewObjPhreeqC
            Me                 => NewObjPhreeqC
            
        else
        
            PreviousObjPhreeqC => FirstObjPhreeqC
            Me                 => FirstObjPhreeqC%Next
            
            do while (associated(Me))
            
                PreviousObjPhreeqC  => Me
                Me                  => Me%Next
                
            enddo
            
            Me                      => NewObjPhreeqC
            PreviousObjPhreeqC%Next => NewObjPhreeqC
            
        endif

        Me%InstanceID = RegisterNewInstance (mPHREEQC_)
        
        !----------------------------------------------------------------------

    end subroutine AllocateInstance    
    !----------------------------------------------------------------------------
       
   
    !----------------------------------------------------------------------------    
    subroutine InitializeInstance
    
        !----------------------------------------------------------------------
        
        Me%PropertyCount = 0  
              
        Me%CalcData%MassOfAllSolutes   = 0
        Me%CalcData%VolumeOfAllSolutes = 0
        Me%CalcData%MassOfWater        = 0
        Me%CalcData%VolumeOfWater      = 0
        Me%CalcData%MassOfSolution     = 0
        Me%CalcData%factor             = 0
        
        !----------------------------------------------------------------------

    end subroutine InitializeInstance    
    !----------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ReadInputFile (FileName)
     
        !Arguments-------------------------------------------------------------
        character(LEN = *)             :: FileName      
        
        !Local-----------------------------------------------------------------
        integer :: zero = 0
        integer :: STAT_
        
        !Begin-----------------------------------------------------------------
        !Associate EnterData
        call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_)
        if (STAT_ .NE. SUCCESS_) &
            stop 'Subroutine ReadInputFile; module ModulePhreeqC. ERR01.'      
                
        call ReadPhreeqCOptions                
        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_) 
        if (STAT_ .NE. SUCCESS_) &
            stop 'Subroutine ReadInputFile; ModuleSedimentQuality. ERR02.'

        !----------------------------------------------------------------------

    end subroutine ReadInputFile
    !--------------------------------------------------------------------------
    
   
    !--------------------------------------------------------------------------
    subroutine ReadPhreeqCOptions

        !Local-----------------------------------------------------------------
        integer :: FromFile
        integer :: STAT_CALL
        integer :: flag
        
        !Begin-----------------------------------------------------------------       
        call GetExtractType (FromFile = FromFile)

        call GetData(Me%PhreeqCOptions%Database     , &
                     Me%ObjEnterData, flag          , &
                     SearchType   = FromFile        , &
                     keyword      = 'DATABASE'      , &
                     ClientModule = 'ModulePhreeqC' , &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR001.'             
        if (flag .EQ. 0) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR002.'
            
        call GetData(Me%PhreeqCOptions%DatabaseAux  , &
                     Me%ObjEnterData, flag          , &
                     SearchType   = FromFile        , &
                     keyword      = 'DATABASE_AUX'  , &
                     default      = ''              , &
                     ClientModule = 'ModulePhreeqC' , &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR003.' 

!        call GetData(Me%PhreeqCOptions%PrintAlways  , &
!                     Me%ObjEnterData, flag          , &
!                     SearchType   = FromFile        , &
!                     keyword      = 'PRINT_ALWAYS'  , &
!                     default      = 0               , &
!                     ClientModule = 'ModulePhreeqC' , &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR004.' 
        
        call GetData(Me%PhreeqCOptions%HPlusDensity , &
                     Me%ObjEnterData, flag          , &
                     SearchType   = FromFile        , &
                     keyword      = 'HPLUS_DENSITY' , &
                     default      = 0.09            , &
                     ClientModule = 'ModulePhreeqC' , &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR005.' 
        
        call GetData(Me%PhreeqCOptions%WaterDensity , &
                     Me%ObjEnterData, flag          , &
                     SearchType   = FromFile        , &
                     keyword      = 'WATER_DENSITY' , &
                     default      = 998.0           , &
                     ClientModule = 'ModulePhreeqC' , &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR006.' 
               
        call GetData(Me%PhreeqCOptions%Redox%Element1,      & 
                     Me%ObjEnterData, flag,                 &
                     SearchType   = FromFile,               &
                     keyword      = 'REDOX_PAIR_ELEMENT_1', & 
                     ClientModule = 'ModulePhreeqC',        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR007.'  
        
        if (flag .NE. 0) then
        
            call GetData(Me%PhreeqCOptions%Redox%Valence1,      & 
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'REDOX_PAIR_VALENCE_1', & 
                         ClientModule = 'ModulePhreeqC',        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. flag .EQ. 0) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR008.'                      

            call GetData(Me%PhreeqCOptions%Redox%Element2,      & 
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'REDOX_PAIR_ELEMENT_2', & 
                         ClientModule = 'ModulePhreeqC'  , &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. flag .EQ. 0) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR009.'                      

            ! Valence 2 ------------------------------------------------------
            call GetData(Me%PhreeqCOptions%Redox%Valence2,      & 
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'REDOX_PAIR_VALENCE_2', & 
                         ClientModule = 'ModulePhreeqC',        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. flag .EQ. 0) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR010.'                      

            call pm_solution_redox(Me%PhreeqCInstanceID,                            &
                                   trim(Me%PhreeqCOptions%Redox%Element1)//char(0), &
                                   Me%PhreeqCOptions%Redox%Valence1,                &
                                   trim(Me%PhreeqCOptions%Redox%Element2)//char(0), &
                                   Me%PhreeqCOptions%Redox%Valence2,                &
                                   STAT_CALL)
            if (STAT_CALL .EQ. 0) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR012.'                                     

        end if        
       
        call GetData(Me%PhreeqCOptions%pHCharge,     &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'PH_CHARGE',     &
                     default      = 0,               &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR013.'                              
        
        call GetData(Me%PhreeqCOptions%pECharge,     &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'PE_CHARGE',     &
                     default      = 0,               &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR014.'        

        if (Me%PhreeqCOptions%pECharge .AND. Me%PhreeqCOptions%pHCharge) &
            stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR015.'
        
        call GetData(Me%PhreeqCOptions%DTSeconds    , &
                     Me%ObjEnterData, flag          , &
                     SearchType   = FromFile        , &
                     keyword      ='DTSECONDS'      , & 
                     default      = 3600.           , & 
                     ClientModule = 'ModulePhreeqC' , &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. ERR016.' 

cd1:   if (flag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Keyword DTSECONDS not found in PhreeqC data file.'
            write(*,*) 'Subroutine ReadPhreeqCOptions; Module ModulePhreeqC. WRN001.'
            write(*,*) 'Assumed ', Me%PhreeqCOptions%DTSeconds, 'seconds (',  Me%PhreeqCOptions%DTSeconds / 60.0, 'hour).'
            write(*,*) 
        end if cd1
        
        !For compatibility with the rest of the program,  
        Me%PhreeqCOptions%DTDay = Me%PhreeqCOptions%DTSeconds / 24.0 / 60.0 / 60.0
        
        call ReadPhreeqCDatabase 

        !----------------------------------------------------------------------

    end subroutine ReadPhreeqCOptions         
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine ReadPhreeqCDatabase

        !Local ----------------------------------------------------------------
        integer :: status

        !----------------------------------------------------------------------

        call pm_read_database (Me%PhreeqCInstanceID, trim(Me%PhreeqCOptions%Database)//char(0), status)      
        if (status .EQ. 0) stop 'Subroutine PhreeqCReadDatabase; Module ModulePhreeqC. ERR001.'
            
        if (Me%PhreeqCOptions%DatabaseAux .NE. '') then
        
            call pm_read_database (Me%PhreeqCInstanceID, trim(Me%PhreeqCOptions%DatabaseAux)//char(0), status)
            if (status .EQ. 0) stop 'Subroutine PhreeqCReadDatabase; Module ModulePhreeqC. ERR002.'
        
        end if
            
        call pm_setup_model (Me%PhreeqCInstanceID, status)      
        if (status .EQ. 0) stop 'Subroutine PhreeqCReadDatabase; Module ModulePhreeqC. ERR003.'
        !----------------------------------------------------------------------

    end subroutine ReadPhreeqCDatabase
    !--------------------------------------------------------------------------
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine SetPhreeqCProperty (PhreeqCID, PropertyID, Property, STAT)
        
        !Arguments-------------------------------------------------------------
        integer                        :: PhreeqCID
        integer                        :: PropertyID
        type(T_ChemistryParameters)    :: Property
        integer, optional, intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer :: ready_              
        integer :: STAT_ !Auxiliar local variable
        real    :: zero = 0.0
        integer :: status
        
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            Me%PropertyCount = Me%PropertyCount + 1
                       
            Me%Properties(Me%PropertyCount)%Params = Property   
            Me%Properties(Me%PropertyCount)%PropertyID = PropertyID        

            select case (Property%Group)
                case (0) !No group
                case (SOLUTION)
                    call SetSolutionProperty(Me%Properties(Me%PropertyCount))
                case (PHASE)
                    call SetPhaseProperty(Me%Properties(Me%PropertyCount))
                case (EXCHANGE)
                    call SetExchangeProperty(Me%Properties(Me%PropertyCount))
                case default
                    stop 'Subroutine SetPhreeqCProperty; Module ModulePhreeqC. ERR001.' !Group not recognized
            end select
                
            STAT_ = SUCCESS_
            
        else
         
            STAT_ = ready_
            
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------
        
    end subroutine SetPhreeqCProperty   
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    subroutine SetSolutionProperty (Property)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty) :: Property

        !Local-----------------------------------------------------------------
        real    :: zero = 0.0
        real    :: database_gfw
        integer :: status
        
        !----------------------------------------------------------------------

        !Pass to PhreeqC Object       
        call pm_conc_add(Me%PhreeqCInstanceID, trim(Property%Params%PhreeqCName)//char(0), zero, status)
        if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR001.'                      
                    
        !Units are required to be mg/L for solution concentrations
        Property%Params%UseUnits = 0
        call pm_conc_use(Me%PhreeqCInstanceID, Property%Params%Charge, Property%Params%UsePhase, Property%Params%UseAs, &
                         Property%Params%UseGFW, Property%Params%UseUnits, Property%Params%UseRedox, status)
        if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR002.'                      
                    
        if (Property%Params%UseAs .EQ. 1) call pm_conc_as(Me%PhreeqCInstanceID, trim(Property%Params%As)//char(0), status)
        if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR003.'                      
                    
        if (Property%Params%UseGFW .EQ. 1) call pm_conc_gfw(Me%PhreeqCInstanceID, Property%Params%GFW, status)
        if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR004.'                      
                           
        if (Property%Params%UsePhase .EQ. 1) call pm_conc_phase(Me%PhreeqCInstanceID, trim(Property%Params%PhaseName)//char(0), & 
                                                         Property%Params%SI, status)
        if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR006.'                      
                    
        if (Property%Params%UseRedox .EQ. 1) call pm_conc_redox(Me%PhreeqCInstanceID,                       &
                                                         trim(Property%Params%RedoxPair%Element1)//char(0), &
                                                         Property%Params%RedoxPair%Valence1,                &
                                                         trim(Property%Params%RedoxPair%Element2)//char(0), &
                                                         Property%Params%RedoxPair%Valence2,                &   
                                                         status)
    if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR007.'                      
                                                                                  
        !Now, save the solution property in the definitive structure
        call pm_conc_save (Me%PhreeqCInstanceID, Property%PhreeqCInputID, Property%PhreeqCResultID, database_gfw, status)
        if (status .EQ. 0) stop 'Subroutine SetSolutionProperty; Module ModulePhreeqC. ERR008.'  
        
        if (Property%Params%UseGFW .EQ. 0) Property%Params%GFW = database_gfw                   
         
        !----------------------------------------------------------------------                    

    end subroutine SetSolutionProperty
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine SetPhaseProperty (Property)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty) :: Property

        !Local-----------------------------------------------------------------
        real                    :: zero = 0.0
        integer                 :: STAT_
        character(StringLength) :: Alternative
        
        !Begin-----------------------------------------------------------------
       
        if (Property%Params%UseAlternativeFormula .EQ. 1) then
            Alternative = Property%Params%AlternativeFormula
        else if (Property%Params%UseAlternativePhase .EQ. 1) then
            Alternative = Property%Params%AlternativePhase
        else
            Alternative = ''
        endif
       
        !Pass to PhreeqC Object
        call pm_ppa_pp(Me%PhreeqCInstanceID, trim(Property%Params%PhreeqCName)//char(0), trim(Alternative)//char(0), &
                       Property%Params%SI, zero, Property%Params%ForceEquality, Property%Params%DissolveOnly, Property%PhreeqCInputID, STAT_)                       
        if (STAT_ .EQ. 0) stop 'Subroutine SetPhaseProperty; Module ModulePhreeqC. ERR001.' 
            
        Property%Params%PhreeqCResultID = Property%Params%PhreeqCInputID                   
        !----------------------------------------------------------------------
    
    end subroutine SetPhaseProperty
    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------
    subroutine SetExchangeProperty (Property)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty) :: Property

        !Local-----------------------------------------------------------------
        real                    :: zero = 0.0
        integer                 :: STAT_
        integer                 :: PhreeqCIndex
        character(StringLength) :: Formula
        
        !Begin-----------------------------------------------------------------
        
        if (Property%Params%UsePhase .EQ. 1) then
            Property%Params%ExType = 1
            Formula = Property%Params%PhaseName
        else if (Property%Params%UseKinetic .EQ. 1) then
            Property%Params%ExType = 2
            Formula = Property%Params%KineticName
        else
            Property%Params%ExType = 0
            Formula = ''
        endif
        
        ! Pass to PhreeqC Object    
        call pm_exa_exchanger(Me%PhreeqCInstanceID, trim(Formula)//char(0), Property%Params%ExType, &
                              trim(Property%Params%Has)//char(0), zero, Property%Params%PhreeqCInputID, STAT_)
        if (STAT_ .EQ. 0) stop 'Subroutine SetExchangeProperty; Module ModulePhreeqC. ERR001.'
         
        call pm_get_species_index(Me%PhreeqCInstanceID, trim(Formula)//char(0), Property%PhreeqCResultID, STAT_) 
        if (STAT_ .EQ. 0) stop 'Subroutine SetExchangeProperty; Module ModulePhreeqC. ERR002.'        
        !----------------------------------------------------------------------
    
    end subroutine SetExchangeProperty
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine GetPhreeqCDT (PhreeqCID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: PhreeqCID
        real,    optional, intent(OUT) :: DTDay
        real,    optional, intent(OUT) :: DTSecond
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_ !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%PhreeqCOptions%DTDay
            if (present(DTSecond)) DTSecond = Me%PhreeqCOptions%DTSeconds

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetPhreeqCDT 
    !--------------------------------------------------------------------------    
    
    
    !--------------------------------------------------------------------------    
    subroutine GetPhreeqCPropIndex (PhreeqCID, PropertyIDNumber, PropertyIndex, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                        :: PhreeqCID
        integer,           intent(IN ) :: PropertyIDNumber
        integer,           intent(OUT) :: PropertyIndex
        integer, optional, intent(OUT) :: STAT
 
        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_              !Auxiliar local variable
        logical :: found
        integer :: CurrentIndex
        
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)    
        
cd1:    if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            found = .false.
            
            do CurrentIndex = 1, Me%PropertyCount

                if (PropertyIDNumber .EQ. Me%Properties(CurrentIndex)%PropertyID)then
                
                    PropertyIndex = CurrentIndex
                    found         = .true.
                    exit
                    
                end if
            
            end do

            if(.not. found) then
                STAT_ = NOT_FOUND_ERR_
            else
                STAT_ = SUCCESS_
            endif

        else 
        
            STAT_ = ready_
            
        end if cd1

        if (present(STAT)) STAT = STAT_
        !----------------------------------------------------------------------
    
    end subroutine GetPhreeqCPropIndex
    !--------------------------------------------------------------------------    
    
    !--------------------------------------------------------------------------    
    subroutine UnGetPhreeqC(PhreeqCID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: PhreeqCID
        integer, dimension(:), pointer :: Array
        integer, intent(OUT), optional :: STAT

        !Local-----------------------------------------------------------------
        integer :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mPHREEQC_, Me%InstanceID, "UnGetPhreeqC")

            STAT_ = SUCCESS_
            
        else
                       
            STAT_ = ready_
            
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetPhreeqC
    !--------------------------------------------------------------------------
    
        
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyPhreeqC (PhreeqCID,           &
                              PropertiesValues,    &
                              SolutionVolume,      &
                              SolutionTemperature, &
                              SolutionpH,          &
                              SolutionpE,          &
                              CellsArrayLB,        &
                              CellsArrayUB,        &
                              OpenPoints,          &
                              STAT)  

        !Arguments---------------------------------------------------------------
        integer                                     :: PhreeqCID
        real,               pointer, dimension(:,:) :: PropertiesValues
        real,               pointer, dimension(:  ) :: SolutionVolume
        real,               pointer, dimension(:  ) :: SolutionTemperature
        real,               pointer, dimension(:  ) :: SolutionpH
        real,               pointer, dimension(:  ) :: SolutionpE
        integer,            intent(IN)              :: CellsArrayLB, CellsArrayUB
        integer, optional,  pointer, dimension(:  ) :: OpenPoints
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_  
        logical                                     :: CalcPoint        
        integer                                     :: CellIndex
        integer                                     :: ready_   

        !------------------------------------------------------------------------                         
            
        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            Me%Ext%PropertiesValues => PropertiesValues
            if (.NOT. associated(Me%Ext%PropertiesValues)) stop 'Subroutine ModifyPhreeqC; Module ModulePhreeqC. ERR001.'

            Me%Ext%SolutionVolume => SolutionVolume
            if (.NOT. associated(Me%Ext%SolutionVolume)) stop 'Subroutine ModifyPhreeqC; Module ModulePhreeqC. ERR002.'
                
            Me%Ext%SolutionTemperature => SolutionTemperature
            if (.NOT. associated(Me%Ext%SolutionTemperature)) stop 'Subroutine ModifyPhreeqC; Module ModulePhreeqC. ERR003.'

            Me%Ext%SolutionpH => SolutionpH
            if (.NOT. associated(Me%Ext%SolutionpH)) stop 'Subroutine ModifyPhreeqC; Module ModulePhreeqC. ERR004.'

            Me%Ext%SolutionpE => SolutionpE
            if (.NOT. associated(Me%Ext%SolutionpE)) stop 'Subroutine ModifyPhreeqC; Module ModulePhreeqC. ERR005.'            
            
do1 :       do CellIndex = CellsArrayLB, CellsArrayUB
            
                !If this module is called from the Interface module, OpenPoint is present
                !and the  module runs for all Openpoints
                !If this module is called from the Lagrangian module, OpenPoint is not present
                !and the  module runs for all volumes
                if (present(OpenPoints)) then
                    if (OpenPoints(CellIndex) == OpenPoint) then !Question: Where came OpenPoint?
                        CalcPoint = .true.
                    else
                        CalcPoint = .false.
                    endif
                else
                    CalcPoint = .true.
                endif

                if (CalcPoint) call MakeCalculations(CellIndex)
                              
            end do do1
                     
            STAT_ = SUCCESS_
            
        else              
         
            STAT_ = ready_
            
        end if cd1

        nullify(Me%Ext%PropertiesValues)
        nullify(Me%Ext%SolutionVolume)
        nullify(Me%Ext%SolutionTemperature)
        nullify(Me%Ext%SolutionpH)
        nullify(Me%Ext%SolutionpE)

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ModifyPhreeqC
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine MakeCalculations(CellIndex)

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex 

        !Local-------------------------------------------------------------------
        integer :: PropertyIndex
        integer :: STAT_
        real    :: ph, pe

        !Begin-------------------------------------------------------------------  
        call CalculateSolutionParameters(CellIndex) !Calculates mass of water and solution density
        
        call pm_set_ph(Me%PhreeqCInstanceID, Me%PhreeqCOptions%pHCharge, Me%Ext%SolutionpH(CellIndex), STAT_)
        if (STAT_ .EQ. 0) stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR001.'

        call pm_set_pe(Me%PhreeqCInstanceID, Me%PhreeqCOptions%pECharge, Me%Ext%SolutionpE(CellIndex), STAT_)
        if (STAT_ .EQ. 0) stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR002.'

        call pm_set_required_data(Me%PhreeqCInstanceID,                  &
                                  Me%CalcData%MassOfWater,               & 
                                  Me%Ext%SolutionTemperature(CellIndex), &
                                  Me%CalcData%DensityOfSolution,         & 
                                  STAT_)                                
        if (STAT_ .EQ. 0) stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR003.'                                  

        !Convert input units from MOHID format to PhreeqC format
        call ConvertInputs (CellIndex)
             
        !Pass to PhreeqCObject all the properties
do1:    do PropertyIndex = 1, Me%PropertyCount !'i' is the index for the current property
          
            select case (Me%Properties(PropertyIndex)%Params%Group)
                case (1, 2, 3, 4, 5, 7)
                    call pm_set_input_value (Me%PhreeqCInstanceID,                              & 
                                             Me%Properties(PropertyIndex)%PhreeqCInputID,       &
                                             Me%Properties(PropertyIndex)%Params%Group,         &
                                             Me%Ext%PropertiesValues(PropertyIndex, CellIndex), &
                                             STAT_)
                    if (STAT_ .EQ. 0) stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR004.'
                case default 
                    stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR005.'
                end select
        end do do1
        
        !Run model       
        call pm_run_model (Me%PhreeqCInstanceID, STAT_) 
        if (STAT_ .EQ. 0) &
            stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR006.'
        
do2:    do PropertyIndex = 1, Me%PropertyCount

            call pm_get_result_value (Me%PhreeqCInstanceID,                               & 
                                      Me%Properties(PropertyIndex)%PhreeqCResultID,       &
                                      Me%Properties(PropertyIndex)%Params%Group,          &
                                      Me%Ext%PropertiesValues(PropertyIndex, CellIndex), &
                                      STAT_)
            if (STAT_ .EQ. 0) stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR007.'            
        
        end do do2
        
        call pm_get_data(Me%PhreeqCInstanceID, Me%CalcData%MassOfWater, ph, pe, STAT_)
        if (STAT_ .EQ. 0) stop 'Subroutine MakeCalculations; module ModulePhreeqC. ERR008.'            
        
        call ConvertResults (CellIndex, ph)
        
        !------------------------------------------------------------------------
       
    end subroutine MakeCalculations
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------
    subroutine CalculateSolutionParameters (CellIndex)

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex

        !Local-------------------------------------------------------------------
        integer :: STAT_
        integer :: PropertyIndex
        real    :: SolutionVolume_L
        
        !Begin-------------------------------------------------------------------       
        !      L         =                m3                * 1000
        SolutionVolume_L = Me%Ext%SolutionVolume(CellIndex) * 1000

        !Find the mass of Solution H+ using pH
        !Considering the GFW of H+ as aprox. 1, and the activity of H+ same as concentration of H+, the result of the expression below gives the mass (in grams) of H+ directly 
        Me%CalcData%MassOfAllSolutes = (10 ** (-Me%Ext%SolutionpH(CellIndex))) * SolutionVolume_L
        
        !With the mass of H+, find the volume of H+ (uses H+ density given by the user or the 0.09 g/L default)
        !             L                =               g              /               g/L
        Me%CalcData%VolumeOfAllSolutes = Me%CalcData%MassOfAllSolutes / Me%PhreeqCOptions%HPlusDensity

        !Second, find the mass and volume of all other solutes (solution concentration properties)
do1:    do PropertyIndex = 1, Me%PropertyCount

            if (Me%Properties(PropertyIndex)%Params%Group .EQ. SOLUTION) then !The units of solution concentrations MUST be in the mg/L format
                
                !                g                = 0.001 *                      mg/L                         *        L        
                Me%Properties(PropertyIndex)%Mass = 0.001 * Me%Ext%PropertiesValues(PropertyIndex, CellIndex) * SolutionVolume_L
                !            g               =              g               +                 g
                Me%CalcData%MassOfAllSolutes = Me%CalcData%MassOfAllSolutes + Me%Properties(PropertyIndex)%Mass 
                
                !                L                  =                   g               /                      g/L
                Me%Properties(PropertyIndex)%Volume = Me%Properties(PropertyIndex)%Mass / Me%Properties(PropertyIndex)%Params%Density
                !             L                =               L                +       L
                Me%CalcData%VolumeOfAllSolutes = Me%CalcData%VolumeOfAllSolutes + Me%Properties(PropertyIndex)%Volume
                
            end if

        end do do1 
        
        !Find the volume of water
        !          L              =        L         -               L
        Me%CalcData%VolumeOfWater = SolutionVolume_L - Me%CalcData%VolumeOfAllSolutes
        
        !Find the mass of water
        !          kg           =             L             *                g/L              * 0.001
        Me%CalcData%MassOfWater = Me%CalcData%VolumeOfWater * (Me%PhreeqCOptions%WaterDensity * 0.001)
        
        !Find Mass of solution
        !           kg             =           kg            + (0.001 *             g               )
        Me%CalcData%MassOfSolution = Me%CalcData%MassOfWater + (0.001 * Me%CalcData%MassOfAllSolutes)
        
        !Find Solution Density
        !            kg/L             =             kg             /       L
        Me%CalcData%DensityOfSolution = Me%CalcData%MassOfSolution / SolutionVolume_L
        
        !Find Factor
        !      kg/kg       =            kg           /            kg
        Me%CalcData%Factor = Me%CalcData%MassOfWater / Me%CalcData%MassOfSolution

        !------------------------------------------------------------------------       
        
    end subroutine CalculateSolutionParameters
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------
    subroutine ConvertInputs (CellIndex)
        !For now, the input units for each group are fixed.
        !In the future, if possible, this will be made more flexible

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex

        !Local-------------------------------------------------------------------
        integer                     :: STAT_ 
        integer                     :: PropertyIndex

        !------------------------------------------------------------------------       
                                       
do1:    do PropertyIndex = 1, Me%PropertyCount
            
            select case (Me%Properties(PropertyIndex)%Params%Group)
                case (SOLUTION) !The input units for SOLUTION concentration properties MUST be mg/L
                                                        
                    if (.NOT.((Me%Properties(PropertyIndex)%Params%PhreeqCName .EQ. 'H(1)') .OR. (Me%Properties(PropertyIndex)%Params%PhreeqCName .EQ. 'E'))) then   
                                                                                         
                        !The factor turns the molality of a solution as if the solution has 1kg of water
                        !                    mol/kgw                      = (                   g              /                   g/mol                ) / (          kg            *       factor      )
                        Me%Ext%PropertiesValues(PropertyIndex, CellIndex) = (Me%Properties(PropertyIndex)%Mass / Me%Properties(PropertyIndex)%Params%GFW) / (Me%CalcData%MassOfWater * Me%CalcData%Factor)
                        
                    end if  
                    
                case (PHASE)
                
                case (EXCHANGE)
                               
                case default                                                      

            end select
                    
        end do do1               
                
        !------------------------------------------------------------------------       
        
    end subroutine ConvertInputs
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------
    subroutine ConvertResults (CellIndex, pH)
        !For now, the input units for each group are fixed.
        !In the future, if possible, this will be made more flexible

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex
        real,    intent(IN) :: pH

        !Local-------------------------------------------------------------------
        integer                     :: PropertyIndex

        !------------------------------------------------------------------------       
                
        !          L              = kg / (0.001 *              g/L              )
        Me%CalcData%VolumeOfWater = 1  / (0.001 * Me%PhreeqCOptions%WaterDensity)
                
        !Find the mass of Solution H+ using pH
        !Considering the GFW of H+ as aprox. 1, and the activity of H+ same as concentration of H+, the result of the expression below gives the mass (in grams) of H+ directly
        Me%CalcData%MassOfAllSolutes = 10 ** (-Me%Ext%SolutionpH(CellIndex))
        
        !With the mass of H+, find the volume of H+ (uses H+ density given by the user or the 0.09 g/L default)
        !             L                =               g              /               g/L
        Me%CalcData%VolumeOfAllSolutes = Me%CalcData%MassOfAllSolutes / Me%PhreeqCOptions%HPlusDensity
                
                
do1:    do PropertyIndex = 1, Me%PropertyCount

            if (Me%Properties(PropertyIndex)%Params%Group .EQ. SOLUTION) then 
                            
                !                g                = 1 kgw *                     mol/kgw                       *                   g/mol                 
                Me%Properties(PropertyIndex)%Mass = 1     * Me%Ext%PropertiesValues(PropertyIndex, CellIndex) * Me%Properties(PropertyIndex)%Params%GFW 
                !            g               =              g               +                 g
                Me%CalcData%MassOfAllSolutes = Me%CalcData%MassOfAllSolutes + Me%Properties(PropertyIndex)%Mass 
                
                !                L                  =                   g               /                      g/L
                Me%Properties(PropertyIndex)%Volume = Me%Properties(PropertyIndex)%Mass / Me%Properties(PropertyIndex)%Params%Density
                !             L                =               L                +       L
                Me%CalcData%VolumeOfAllSolutes = Me%CalcData%VolumeOfAllSolutes + Me%Properties(PropertyIndex)%Volume
                
            end if

        end do do1                 
                    
        !           g              = 1000 g  +              g
        Me%CalcData%MassOfSolution = 1000    + Me%CalcData%MassOfAllSolutes
        
        !            L               =             L             +               L
        Me%CalcData%VolumeOfSolution = Me%CalcData%VolumeOfWater + Me%CalcData%VolumeOfAllSolutes
         
        !           g/L               =              g             /              L
        Me%CalcData%DensityOfSolution = Me%CalcData%MassOfSolution / Me%CalcData%VolumeOfSolution 
               
do2:    do PropertyIndex = 1, Me%PropertyCount

            if (Me%Properties(PropertyIndex)%Params%Group .EQ. SOLUTION) then 
            
                !                   mg/L                          = 1000 *               g                  /              L  
                Me%Ext%PropertiesValues(PropertyIndex, CellIndex) = 1000 * Me%Properties(PropertyIndex)%Mass / Me%CalcData%VolumeOfSolution
                                            
            end if

        end do do2                 
        
        !------------------------------------------------------------------------       
                
    end subroutine ConvertResults
    !----------------------------------------------------------------------------
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !----------------------------------------------------------------------------
    subroutine KillPhreeqC(PhreeqCID, STAT)

        !Arguments---------------------------------------------------------------
        integer                        :: PhreeqCID
        integer, optional, intent(OUT) :: STAT

        !External----------------------------------------------------------------
        integer :: ready_ 
        integer :: status             

        !Local-------------------------------------------------------------------
        integer :: STAT_
        integer :: nUsers 

        !------------------------------------------------------------------------                      

        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mPHREEQC_,  Me%InstanceID)

cd2 :       if (nUsers == 0) then

                !Question: If I create a instance and after that kill this instance, this piece of code will run?
                call pm_kill(Me%PhreeqCInstanceID, status)
                if (status .EQ. 0) &
                    stop 'Subroutine KillPhreeqC; module ModulePhreeqC. ERR01.'
                
                call DeallocateInstance 

                PhreeqCID = 0
                STAT_ = SUCCESS_

            end if cd2
        
        else cd1
        
            STAT_ = ready_
        
        end if cd1


        if (present(STAT)) STAT = STAT_
    !------------------------------------------------------------------------

    end subroutine KillPhreeqC
    !------------------------------------------------------------------------


    !------------------------------------------------------------------------
    subroutine DeallocateInstance 

        !Local-----------------------------------------------------------------
        type (T_PhreeqC), pointer :: AuxObjPhreeqC
        type (T_PhreeqC), pointer :: PreviousObjPhreeqC

        !Updates pointers
        if (Me%InstanceID == FirstObjPhreeqC%InstanceID) then
        
            FirstObjPhreeqC => FirstObjPhreeqC%Next
            
        else
        
            PreviousObjPhreeqC => FirstObjPhreeqC
            AuxObjPhreeqC      => FirstObjPhreeqC%Next
            
            do while (AuxObjPhreeqC%InstanceID /= Me%InstanceID)
            
                PreviousObjPhreeqC => AuxObjPhreeqC
                AuxObjPhreeqC      => AuxObjPhreeqC%Next
                
            enddo

            !Now update linked list
            PreviousObjPhreeqC%Next => AuxObjPhreeqC%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 
            
    end subroutine DeallocateInstance
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !----------------------------------------------------------------------------
    subroutine Ready (PhreeqCID, ready_) 

        !Arguments-------------------------------------------------------------
        integer :: PhreeqCID
        integer :: ready_
        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (PhreeqCID > 0) then
            
            call LocateObjPhreeqC (PhreeqCID)
            ready_ = VerifyReadLock (mPHREEQC_, Me%InstanceID)

        else

            ready_ = OFF_ERR_

        end if cd1
        !----------------------------------------------------------------------

    end subroutine Ready
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine LocateObjPhreeqC (PhreeqCID)

    !Arguments-------------------------------------------------------------
        integer :: PhreeqCID
    !--------------------------------------------------------------------------

        Me => FirstObjPhreeqC
        do while (associated (Me))
            if (Me%InstanceID == PhreeqCID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))   &
            stop 'Subroutine LocateObjPhreeqC; Module ModulePhreeqC. ERR001.'
    !--------------------------------------------------------------------------
    
    end subroutine LocateObjPhreeqC
    !--------------------------------------------------------------------------


end module ModulePhreeqC    
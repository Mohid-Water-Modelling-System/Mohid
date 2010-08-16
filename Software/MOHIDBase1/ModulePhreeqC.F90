!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
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

#ifdef _PHREEQC_

Module ModulePhreeqC

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
    private ::          ReadPhreeqCProperties
    private ::              ConstructProperty
    private ::                  ReadChemistryParameters
    private ::                      ReadChemistryConcGroupParam
    private ::                      ReadChemistryPhasesGroupParam
!    private ::                      ReadChemistrySolidGroupParam
!    private ::                      ReadChemistryGasGroupParam
!    private ::                      ReadChemistrySurfGroupParam
    private ::                      ReadChemistrySpeciesGroupParam 
    private ::                      ReadChemistryExcGroupParam       
    private ::              AddProperty
                
    private :: SetPhreeqCProperty
    private ::      SetSpeciesProperty
    private ::      SetConcentrationProperty
    private ::      SetPhaseProperty
    private ::      SetExchangeProperty 
    
    
    !Selector
    public  :: GetPhreeqCDT 
    public  :: GetPhreeqCNeedSoilDryDensity
    public  :: GetPhreeqCPropIndex      
    public  :: UngetPhreeqC           
                        
    !Modifier
    public  :: ModifyPhreeqC 
    private ::      MakeCalculations
    private ::          CalculateSolutionDensity
    private ::          ConvertInput
    private ::          ConvertResult
!    private ::          CalculateSolutionParameters
!    private ::          ConvertInputs
!    private ::          ConvertResults 
    private ::          PrintDataInput
    private ::          PrintDataOutput       
        
    !Destructor
    public  :: KillPhreeqC                                                     
    private ::      DeallocateInstance

    !Management
    private :: EndWithError
    private :: OpenDebugFile
    private :: CloseDebugFile
    private :: Ready
    private ::      LocateObjPhreeqC
        
    !Constants-----------------------------------------------------------------
    integer,                       parameter :: MaxStrLength     = StringLength !StringLength is defined in ModuleGlobalData 
    integer,                       parameter :: MaxProperties    = 25    
    character(LEN = StringLength), parameter :: prop_block_begin = '<beginproperty>'
    character(LEN = StringLength), parameter :: prop_block_end   = '<endproperty>'       
    
    real,                          parameter :: DefaultWaterDensity = .998 !kg/L at 20ºC

    !Constants-----------------------------------------------------------------
    integer, parameter :: OTHER          = 0
    integer, parameter :: CONCENTRATION  = 1
    integer, parameter :: PHASE          = 2
    integer, parameter :: GAS            = 3
    integer, parameter :: SURFACE        = 4
    integer, parameter :: SPECIES        = 5
    integer, parameter :: EXCHANGE       = 6
    integer, parameter :: PUREPHASESOLID = 7
    integer, parameter :: PUREPHASEGAS   = 8 
        
    integer, parameter :: SOLID_PHASE = 1
    integer, parameter :: GAS_PHASE   = 2
    
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
    integer, parameter :: mol_kgw   = 19 
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
    
    type T_RedoxPair
        character(20) :: Element1
        character(20) :: Element2
        real          :: Valence1
        real          :: Valence2        
    end type T_RedoxPair
    
    type T_ChemistryParameters
        character(StringLength) :: PhreeqCName
        character(StringLength) :: Units
        character(StringLength) :: PhaseName
        character(StringLength) :: KineticName
        character(StringLength) :: AlternativePhase
        character(StringLength) :: AlternativeFormula
        character(StringLength) :: As
        character(StringLength) :: Has
        character(StringLength) :: Formula
        type(T_RedoxPair)       :: RedoxPair
        integer                 :: Group
        integer                 :: ExType
        integer                 :: Charge
        integer                 :: DoNotChange           = 0
        real                    :: SI                        !Saturation Index
        real                    :: GFW                       !Gram Formula Weight        
        real                    :: Density                   !Density for solution concentrations.
        integer                 :: UseGFW                = 0     
        integer                 :: UseAs                 = 0
        integer                 :: UseAlternativePhase   = 0
        integer                 :: UseAlternativeFormula = 0
        integer                 :: UsePhase              = 0
        integer                 :: UseRedox              = 0
        integer                 :: UseUnits              = 0
        integer                 :: UseKinetic            = 0
        integer                 :: ForceEquality         = 0
        integer                 :: DissolveOnly          = 0
        integer                 :: PhaseType             = 1 ! 1 - Solid, 2 - Gas
        logical                 :: Debug                 = .false.
    end type T_ChemistryParameters    
    
    type T_PhreeqCUnits
        integer             :: MasterSpecies 
        integer             :: Species       
        integer             :: Alkalinity    
        integer             :: SolidPhases   
        integer             :: Exchangers    
        integer             :: Solution      = mol_kgw
    end type T_PhreeqCUnits      
        
    type T_PhreeqCOptions
        character(len=2048)  :: Database                  !Path for database
        character(len=2048)  :: DatabaseAux  = ''         !Path for auxiliary database
        logical              :: PrintInput   = .false.    !For DEBUG 
        logical              :: PrintOutput  = .false.   
        real                 :: DTSeconds    = null_real 
        real                 :: DTDay        = null_real
        real                 :: HPlusDensity              !g/L
        real                 :: WaterDensity              !g/L
        type(T_RedoxPair)    :: Redox
        type(T_PhreeqCUnits) :: Units
        integer              :: pHCharge
        integer              :: pECharge
        integer              :: UseFixedTemperature = 0
        integer              :: UseFixedpH          = 0
        integer              :: UseFixedpE          = 0
        integer              :: UseFixedSoilDensity = 0
        real                 :: FixedTemperature
        real                 :: FixedpH
        real                 :: FixedpE
        real                 :: FixedSoilDensity
        integer              :: UseExchanger     = 0
        integer              :: UseSolidPhase    = 0
        integer              :: UseGasPhase      = 0
        integer              :: UseGas           = 0
        integer              :: UseSolidSolution = 0
        integer              :: UseSurface       = 0
        logical              :: Debug            = .false.
        logical              :: PrintAllways     = .false.        
    end type T_PhreeqCOptions


    !Types---------------------------------------------------------------------            
       
    type T_PhreeqCProperty
        type (T_PropertyID)              :: ID
        integer                          :: Index
        integer                          :: PropertyID 
        integer                          :: PhreeqCInputID
        integer                          :: PhreeqCResultID
        type(T_ChemistryParameters)      :: Params
        real                             :: PropertyValue !Used to store temporary cell value for the property
        real                             :: Volume
        type(T_PhreeqCProperty), pointer :: Next => null()
        type(T_PhreeqCProperty), pointer :: Prev => null()
    end type T_PhreeqCProperty
    
    type T_External
        real, pointer, dimension(:,:) :: PropertiesValues
        real, pointer, dimension(:  ) :: WaterVolume
        real, pointer, dimension(:  ) :: WaterMass
        real, pointer, dimension(:  ) :: Temperature
        real, pointer, dimension(:  ) :: pH
        real, pointer, dimension(:  ) :: pE
        real, pointer, dimension(:  ) :: SolidMass
    end type T_External
    
    type T_Calculations
        real :: MassOfAllSolutes
        real :: VolumeOfAllSolutes
        real :: MassOfWater
        real :: VolumeOfWater
        real :: MassOfSolution
        real :: VolumeOfSolution
        real :: DensityOfSolution
    end type T_Calculations
    
    type T_PhreeqC
        private
        integer                                           :: InstanceID           !ID of the ModulePhreeqC instance 
        integer                                           :: PhreeqCInstanceID    !ID of the PhreeqC Object instance linked to the InstanceID instance
        integer                                           :: ObjEnterData = 0     !Instance of ModuleEnterData
        type(T_PhreeqC) , pointer                         :: Next                 !Collection of instances of ModulePhreeqC
        type(T_PhreeqCOptions)                            :: MOptions             !Global options read from the PhreeqC Input File
        type(T_External)                                  :: Ext                  !Pointers to Water Mass, Properties Values and other required data 
        !type(T_PhreeqCProperty), dimension(MaxProperties) :: Properties           !Info about each property. Use or delete?
        integer                                           :: PropertyCount        !Number of properties
        integer                                           :: PropertiesNumber = 0 !Use this or PropertyCount or both?
        real, dimension(:), pointer                       :: PropertyValues => null()      
        type(T_Calculations)                              :: CalcData             !Temporary data for calculations   
        type(T_PhreeqCProperty), pointer                  :: FirstProperty => null()
        type(T_PhreeqCProperty), pointer                  :: LastProperty  => null()  
        integer                                           :: DebugFileUnit = -1
    end type T_PhreeqC


    !Global Module Variables---------------------------------------------------
    type (T_PhreeqC), pointer :: FirstObjPhreeqC
    type (T_PhreeqC), pointer :: Me      
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !INTERFACE TO PhreeqCLIB

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    interface 
        subroutine pmStart [ALIAS:'?pm_start@@YAXPAH0@Z'](a, b)
            integer*4 a [REFERENCE]
            integer*4 b [REFERENCE]
        end
	    
        subroutine pm_conc_use [ALIAS:'?pm_conc_use@@YAXPAH0000000@Z'] (a, b, c, d, e, f, g, h)
            integer*4 a [REFERENCE] 
            integer*4 b [REFERENCE] 
            integer*4 c [REFERENCE]
            integer*4 d [REFERENCE]
            integer*4 e [REFERENCE] 
            integer*4 f [REFERENCE] 
            integer*4 g [REFERENCE]
            integer*4 h [REFERENCE]
        end
        
        subroutine pm_conc_as [ALIAS:'?pm_conc_as@@YAXPAHPAD0H@Z'] (a, b, c)
            integer*4        a [REFERENCE] 
            character(Len=*) b [REFERENCE] 
            integer*4        c [REFERENCE]
        end
        
        subroutine pm_use_ppa [ALIAS:'?pm_use_ppa@@YAXPAH00@Z'] (a, b, c)
            integer*4 a [REFERENCE] 
            integer*4 b [REFERENCE] 
            integer*4 c [REFERENCE]
        end          

        subroutine pm_use_exa [ALIAS:'?pm_use_exa@@YAXPAH00@Z'] (a, b, c)
            integer*4 a [REFERENCE] 
            integer*4 b [REFERENCE] 
            integer*4 c [REFERENCE]
        end

        subroutine pm_use_sa [ALIAS:'?pm_use_sa@@YAXPAH00@Z'] (a, b, c)
            integer*4 a [REFERENCE] 
            integer*4 b [REFERENCE] 
            integer*4 c [REFERENCE]
        end        
        
        subroutine pmSetUse [ALIAS:'?pm_set_use@@YAXPAH000000@Z'] (a, b, c, d, e, f, g)
            integer*4 a [REFERENCE] 
            integer*4 b [REFERENCE] 
            integer*4 c [REFERENCE]
            integer*4 d [REFERENCE] 
            integer*4 e [REFERENCE]
            integer*4 f [REFERENCE] 
            integer*4 g [REFERENCE] 
        end        
                
        subroutine pm_get_species_index [ALIAS:'?pm_get_species_index@@YAXPAHPAD00H@Z'] (a, b, c, d)
            integer*4        a [REFERENCE] 
            character(Len=*) b [REFERENCE] 
            integer*4        c [REFERENCE]
            integer*4        d [REFERENCE]
        end
    
        subroutine pmRunPhreeqC [ALIAS:'?pm_run_model@@YAXPAH0@Z'](a, b)
	        integer*4 a [REFERENCE]
	        integer*4 b [REFERENCE]
	    end

        subroutine pmSetupModel [ALIAS:'?pm_setup_model@@YAXPAH0@Z'](a, b)
	        integer*4 a [REFERENCE]
	        integer*4 b [REFERENCE]
	    end

        subroutine pmLoadDatabase [ALIAS:'?pm_read_database@@YAXPAHPAD0@Z'] (a, b, c)
            integer*4        a [REFERENCE] 
            character(Len=*) b [REFERENCE] 
            integer*4        c [REFERENCE]
        end

        subroutine pmKill [ALIAS:'?pm_kill@@YAXPAH0@Z'](a, b)
	        integer*4 a [REFERENCE]
	        integer*4 b [REFERENCE]
	    end        
    end interface
       
    interface pmSetSolutionData
	    subroutine PMSetSolutionDataA [ALIAS: '?pm_solution_data@@YAXPAHPAN0110@Z'] (a, b, c, d, e, f)
	        integer*4 a [REFERENCE]
	        real*8    b [REFERENCE]
	        integer*4 c [REFERENCE]
	        real*8    d [REFERENCE]
	        real*8    e [REFERENCE]
	        integer*4 f [REFERENCE]
	    end
	    
!	    subroutine PMSetSolutionDataB [ALIAS: '?pm_solution_data@@YAXPAHPAN0110@Z'] (a, b, c, d, e, f)
!	        integer*4 a [REFERENCE]
!	        real*16   b [REFERENCE]
!	        integer*4 c [REFERENCE]
!	        real*16   d [REFERENCE]
!	        real*16   e [REFERENCE]
!	        integer*4 f [REFERENCE]
!	    end
!	    
!	    subroutine PMSetSolutionDataC [ALIAS: '?pm_solution_data@@YAXPAHPAN0110@Z'] (a, b, c, d, e, f)
!	        integer*4 a [REFERENCE]
!	        real*4    b [REFERENCE]
!	        integer*4 c [REFERENCE]
!	        real*4    d [REFERENCE]
!	        real*4    e [REFERENCE]
!	        integer*4 f [REFERENCE]
!	    end
	end interface
           
    interface pm_solution_redox
        subroutine pm_solution_redoxA [ALIAS: '?pm_solution_redox@@YAXPAHPADPAM120HH@Z'](a1, a2, a3, a4, a5, a6)
            integer*4        a1 [REFERENCE]
            character(Len=*) a2 [REFERENCE]
            real*4           a3 [REFERENCE]
            character(Len=*) a4 [REFERENCE]
            real*4           a5 [REFERENCE]
            integer*4        a6 [REFERENCE]
        end
        
        subroutine pm_solution_redoxB [ALIAS: '?pm_solution_redox@@YAXPAHPADPAN120HH@Z'](b1, b2, b3, b4, b5, b6)
            integer*4        b1 [REFERENCE]
            character(Len=*) b2 [REFERENCE]
            real*8           b3 [REFERENCE]
            character(Len=*) b4 [REFERENCE]
            real*8           b5 [REFERENCE]
            integer*4        b6 [REFERENCE]
        end
        
        subroutine pm_solution_redoxC [ALIAS: '?pm_solution_redox@@YAXPAHPADPAO120HH@Z'](c1, c2, c3, c4, c5, c6)
            integer*4        c1 [REFERENCE]
            character(Len=*) c2 [REFERENCE]
            real*16          c3 [REFERENCE]
            character(Len=*) c4 [REFERENCE]
            real*16          c5 [REFERENCE]
            integer*4        c6 [REFERENCE]
        end        
    end interface       
    
    interface pm_conc_add
        subroutine pm_conc_addA [ALIAS:'?pm_conc_add@@YAXPAHPADPAM0H@Z'] (a1, a2, a3, a4)
            integer*4        a1 [REFERENCE] 
            character(Len=*) a2 [REFERENCE]
            real*4           a3 [REFERENCE]
            integer*4        a4 [REFERENCE]
        end
        
        subroutine pm_conc_addB [ALIAS:'?pm_conc_add@@YAXPAHPADPAN0H@Z'] (b1, b2, b3, b4)
            integer*4        b1 [REFERENCE] 
            character(Len=*) b2 [REFERENCE]
            real*8           b3 [REFERENCE]
            integer*4        b4 [REFERENCE]
        end
        
        subroutine pm_conc_addC [ALIAS:'?pm_conc_add@@YAXPAHPADPAO0H@Z'] (c1, c2, c3, c4)
            integer*4        c1 [REFERENCE] 
            character(Len=*) c2 [REFERENCE]
            real*16          c3 [REFERENCE]
            integer*4        c4 [REFERENCE]
        end
    end interface
    
    interface pm_conc_gfw
        subroutine pm_conc_gfwA [ALIAS:'?pm_conc_gfw@@YAXPAHPAM0@Z'] (a1, a2, a3)
            integer*4 a1 [REFERENCE] 
            real*4    a2 [REFERENCE] 
            integer*4 a3 [REFERENCE]
        end
        
        subroutine pm_conc_gfwB [ALIAS:'?pm_conc_gfw@@YAXPAHPAN0@Z'] (b1, b2, b3)
            integer*4 b1 [REFERENCE] 
            real*8    b2 [REFERENCE] 
            integer*4 b3 [REFERENCE]
        end
                
        subroutine pm_conc_gfwC [ALIAS:'?pm_conc_gfw@@YAXPAHPAO0@Z'] (c1, c2, c3)
            integer*4 c1 [REFERENCE] 
            real*16   c2 [REFERENCE] 
            integer*4 c3 [REFERENCE]
        end    
    end interface

    interface pm_conc_phase
        subroutine pm_conc_phaseA [ALIAS:'?pm_conc_phase@@YAXPAHPADPAM0@Z'] (a1, a2, a3, a4)
            integer*4        a1 [REFERENCE] 
            character(Len=*) a2 [REFERENCE] 
            real*4           a3 [REFERENCE]
            integer*4        a4 [REFERENCE]
        end

        subroutine pm_conc_phaseB [ALIAS:'?pm_conc_phase@@YAXPAHPADPAN0@Z'] (b1, b2, b3, b4)
            integer*4        b1 [REFERENCE] 
            character(Len=*) b2 [REFERENCE] 
            real*8           b3 [REFERENCE]
            integer*4        b4 [REFERENCE]
        end

        subroutine pm_conc_phaseC [ALIAS:'?pm_conc_phase@@YAXPAHPADPAO0@Z'] (c1, c2, c3, c4)
            integer*4        c1 [REFERENCE] 
            character(Len=*) c2 [REFERENCE] 
            real*16          c3 [REFERENCE]
            integer*4        c4 [REFERENCE]
        end
    end interface

    interface pm_conc_redox
        subroutine pm_conc_redoxA [ALIAS:'?pm_conc_redox@@YAXPAHPADPAM120HH@Z'] (a1, a2, a3, a4, a5, a6)
            integer*4        a1 [REFERENCE] 
            character(Len=*) a2 [REFERENCE] 
            real*4           a3 [REFERENCE]
            character(Len=*) a4 [REFERENCE] 
            real*4           a5 [REFERENCE]
            integer*4        a6 [REFERENCE]
        end
        
        subroutine pm_conc_redoxB [ALIAS:'?pm_conc_redox@@YAXPAHPADPAN120HH@Z'] (b1, b2, b3, b4, b5, b6)
            integer*4        b1 [REFERENCE] 
            character(Len=*) b2 [REFERENCE] 
            real*8           b3 [REFERENCE]
            character(Len=*) b4 [REFERENCE] 
            real*8           b5 [REFERENCE]
            integer*4        b6 [REFERENCE]
        end

        subroutine pm_conc_redoxC [ALIAS:'?pm_conc_redox@@YAXPAHPADPAO120HH@Z'] (c1, c2, c3, c4, c5, c6)
            integer*4        c1 [REFERENCE] 
            character(Len=*) c2 [REFERENCE] 
            real*16          c3 [REFERENCE]
            character(Len=*) c4 [REFERENCE] 
            real*16          c5 [REFERENCE]
            integer*4        c6 [REFERENCE]
        end
    end interface

    interface pm_conc_save
        subroutine pm_conc_saveA [ALIAS:'?pm_conc_save@@YAXPAH00PAM0@Z'] (a1, a2, a3, a4, a5)
            integer*4 a1 [REFERENCE] 
            integer*4 a2 [REFERENCE] 
            integer*4 a3 [REFERENCE]
            real*4    a4 [REFERENCE]
            integer*4 a5 [REFERENCE]
        end

        subroutine pm_conc_saveB [ALIAS:'?pm_conc_save@@YAXPAH00PAN0@Z'] (b1, b2, b3, b4, b5)
            integer*4 b1 [REFERENCE] 
            integer*4 b2 [REFERENCE] 
            integer*4 b3 [REFERENCE]
            real*8    b4 [REFERENCE]
            integer*4 b5 [REFERENCE]
        end

        subroutine pm_conc_saveC [ALIAS:'?pm_conc_save@@YAXPAH00PAO0@Z'] (c1, c2, c3, c4, c5)
            integer*4 c1 [REFERENCE] 
            integer*4 c2 [REFERENCE] 
            integer*4 c3 [REFERENCE]
            real*16   c4 [REFERENCE]
            integer*4 c5 [REFERENCE]
        end
    end interface

    interface pm_ppa_pp 
        subroutine pm_ppa_ppA [ALIAS:'?pm_ppa_pp@@YAXPAHPAD1PAM20000HH@Z'] (a1, a2, a3, a4, a5, a6, a7, a8, a9)
            integer*4        a1 [REFERENCE] 
            character(Len=*) a2 [REFERENCE] 
            character(Len=*) a3 [REFERENCE]
            real*4           a4 [REFERENCE]
            real*4           a5 [REFERENCE] 
            integer*4        a6 [REFERENCE] 
            integer*4        a7 [REFERENCE]
            integer*4        a8 [REFERENCE]
            integer*4        a9 [REFERENCE]
        end
    
        subroutine pm_ppa_ppB [ALIAS:'?pm_ppa_pp@@YAXPAHPAD1PAN20000HH@Z'] (b1, b2, b3, b4, b5, b6, b7, b8, b9)
            integer*4        b1 [REFERENCE] 
            character(Len=*) b2 [REFERENCE] 
            character(Len=*) b3 [REFERENCE]
            real*8           b4 [REFERENCE]
            real*8           b5 [REFERENCE] 
            integer*4        b6 [REFERENCE] 
            integer*4        b7 [REFERENCE]
            integer*4        b8 [REFERENCE]
            integer*4        b9 [REFERENCE]
        end

        subroutine pm_ppa_ppC [ALIAS:'?pm_ppa_pp@@YAXPAHPAD1PAO20000HH@Z'] (c1, c2, c3, c4, c5, c6, c7, c8, c9)
            integer*4        c1 [REFERENCE] 
            character(Len=*) c2 [REFERENCE] 
            character(Len=*) c3 [REFERENCE]
            real*16          c4 [REFERENCE]
            real*16          c5 [REFERENCE] 
            integer*4        c6 [REFERENCE] 
            integer*4        c7 [REFERENCE]
            integer*4        c8 [REFERENCE]
            integer*4        c9 [REFERENCE]
        end
    end interface

    interface pm_exa_exchanger
        subroutine pm_exa_exchangerA [ALIAS:'?pm_exa_exchanger@@YAXPAHPAD01PAM00HH@Z'] (a1, a2, a3, a4, a5, a6, a7)
            integer*4        a1 [REFERENCE] 
            character(Len=*) a2 [REFERENCE] 
            integer*4        a3 [REFERENCE]
            character(Len=*) a4 [REFERENCE] 
            real*4           a5 [REFERENCE]
            integer*4        a6 [REFERENCE]
            integer*4        a7 [REFERENCE]
        end

        subroutine pm_exa_exchangerB [ALIAS:'?pm_exa_exchanger@@YAXPAHPAD01PAN00HH@Z'] (b1, b2, b3, b4, b5, b6, b7)
            integer*4        b1 [REFERENCE] 
            character(Len=*) b2 [REFERENCE] 
            integer*4        b3 [REFERENCE]
            character(Len=*) b4 [REFERENCE] 
            real*8           b5 [REFERENCE]
            integer*4        b6 [REFERENCE]
            integer*4        b7 [REFERENCE]
        end

        subroutine pm_exa_exchangerC [ALIAS:'?pm_exa_exchanger@@YAXPAHPAD01PAO00HH@Z'] (c1, c2, c3, c4, c5, c6, c7)
            integer*4        c1 [REFERENCE] 
            character(Len=*) c2 [REFERENCE] 
            integer*4        c3 [REFERENCE]
            character(Len=*) c4 [REFERENCE] 
            real*16          c5 [REFERENCE]
            integer*4        c6 [REFERENCE]
            integer*4        c7 [REFERENCE]
        end
    end interface

    interface pm_sa_options
        subroutine pm_sa_optionsA [ALIAS:'?pm_sa_options@@YAXPAH00PAM00@Z'] (a1, a2, a3, a4, a5, a6)
            integer*4        a1 [REFERENCE] 
            integer*4        a2 [REFERENCE] 
            integer*4        a3 [REFERENCE]
            real*4           a4 [REFERENCE] 
            integer*4        a5 [REFERENCE]
            integer*4        a6 [REFERENCE]
        end

        subroutine pm_sa_optionsB [ALIAS:'?pm_sa_options@@YAXPAH00PAN00@Z'] (b1, b2, b3, b4, b5, b6)
            integer*4        b1 [REFERENCE] 
            integer*4        b2 [REFERENCE] 
            integer*4        b3 [REFERENCE]
            real*8           b4 [REFERENCE] 
            integer*4        b5 [REFERENCE]
            integer*4        b6 [REFERENCE]
        end

        subroutine pm_sa_optionsC [ALIAS:'?pm_sa_options@@YAXPAH00PAO00@Z'] (c1, c2, c3, c4, c5, c6)
            integer*4        c1 [REFERENCE] 
            integer*4        c2 [REFERENCE] 
            integer*4        c3 [REFERENCE]
            real*16          c4 [REFERENCE] 
            integer*4        c5 [REFERENCE]
            integer*4        c6 [REFERENCE]
        end
    end interface

    interface pm_sa_surface
        subroutine pm_sa_surfaceA [ALIAS:'?pm_sa_surface@@YAXPAH0PAD1PAM22000HH@Z'] (a1, a2, a3, a4, a5, a6, a7, a8)
            integer*4        a1 [REFERENCE] 
            integer*4        a2 [REFERENCE] 
            character(Len=*) a3 [REFERENCE]
            character(Len=*) a4 [REFERENCE]
            real*4           a5 [REFERENCE] 
            real*4           a6 [REFERENCE] 
            real*4           a7 [REFERENCE]
            integer*4        a8 [REFERENCE]
        end

        subroutine pm_sa_surfaceB [ALIAS:'?pm_sa_surface@@YAXPAH0PAD1PAN22000HH@Z'] (b1, b2, b3, b4, b5, b6, b7, b8)
            integer*4        b1 [REFERENCE] 
            integer*4        b2 [REFERENCE] 
            character(Len=*) b3 [REFERENCE]
            character(Len=*) b4 [REFERENCE]
            real*8           b5 [REFERENCE] 
            real*8           b6 [REFERENCE] 
            real*8           b7 [REFERENCE]
            integer*4        b8 [REFERENCE]
        end

        subroutine pm_sa_surfaceC [ALIAS:'?pm_sa_surface@@YAXPAH0PAD1PAO22000HH@Z'] (c1, c2, c3, c4, c5, c6, c7, c8)
            integer*4        c1 [REFERENCE] 
            integer*4        c2 [REFERENCE] 
            character(Len=*) c3 [REFERENCE]
            character(Len=*) c4 [REFERENCE]
            real*16          c5 [REFERENCE] 
            real*16          c6 [REFERENCE] 
            real*16          c7 [REFERENCE]
            integer*4        c8 [REFERENCE]
        end
    end interface

    interface pm_set_required_data
        subroutine pm_set_required_dataA [ALIAS:'?pm_set_required_data@@YAXPAHPAM110@Z'] (a1, a2, a3, a4, a5)
            integer*4        a1 [REFERENCE] 
            real*4           a2 [REFERENCE] 
            real*4           a3 [REFERENCE] 
            real*4           a4 [REFERENCE] 
            integer*4        a5 [REFERENCE]
        end

        subroutine pm_set_required_dataB [ALIAS:'?pm_set_required_data@@YAXPAHPAN110@Z'] (b1, b2, b3, b4, b5)
            integer*4        b1 [REFERENCE] 
            real*8           b2 [REFERENCE] 
            real*8           b3 [REFERENCE] 
            real*8           b4 [REFERENCE] 
            integer*4        b5 [REFERENCE]
        end

        subroutine pm_set_required_dataC [ALIAS:'?pm_set_required_data@@YAXPAHPAO110@Z'] (c1, c2, c3, c4, c5)
            integer*4        c1 [REFERENCE] 
            real*16          c2 [REFERENCE] 
            real*16          c3 [REFERENCE] 
            real*16          c4 [REFERENCE] 
            integer*4        c5 [REFERENCE]
        end
    end interface

    interface pmSetInputValue
        subroutine pm_set_input_valueA [ALIAS:'?pm_set_input_value@@YAXPAH00PAM0@Z'] (a1, a2, a3, a4, a5)
            integer*4        a1 [REFERENCE] 
            integer*4        a2 [REFERENCE] 
            integer*4        a3 [REFERENCE] 
            real*4           a4 [REFERENCE] 
            integer*4        a5 [REFERENCE]
        end           

        subroutine pm_set_input_valueB [ALIAS:'?pm_set_input_value@@YAXPAH00PAN0@Z'] (b1, b2, b3, b4, b5)
            integer*4        b1 [REFERENCE] 
            integer*4        b2 [REFERENCE] 
            integer*4        b3 [REFERENCE] 
            real*8           b4 [REFERENCE] 
            integer*4        b5 [REFERENCE]
        end           

        subroutine pm_set_input_valueC [ALIAS:'?pm_set_input_value@@YAXPAH00PAO0@Z'] (c1, c2, c3, c4, c5)
            integer*4        c1 [REFERENCE] 
            integer*4        c2 [REFERENCE] 
            integer*4        c3 [REFERENCE] 
            real*16          c4 [REFERENCE] 
            integer*4        c5 [REFERENCE]
        end           
    end interface

    interface pmGetResultValue
        subroutine pm_get_result_valueA [ALIAS:'?pm_get_result_value@@YAXPAH00PAM0@Z'] (a1, a2, a3, a4, a5)
            integer*4        a1 [REFERENCE] 
            integer*4        a2 [REFERENCE] 
            integer*4        a3 [REFERENCE] 
            real*4           a4 [REFERENCE] 
            integer*4        a5 [REFERENCE]
        end
        
        subroutine pm_get_result_valueB [ALIAS:'?pm_get_result_value@@YAXPAH00PAN0@Z'] (b1, b2, b3, b4, b5)
            integer*4        b1 [REFERENCE] 
            integer*4        b2 [REFERENCE] 
            integer*4        b3 [REFERENCE] 
            real*8           b4 [REFERENCE] 
            integer*4        b5 [REFERENCE]
        end
        
        subroutine pm_get_result_valueC [ALIAS:'?pm_get_result_value@@YAXPAH00PAO0@Z'] (c1, c2, c3, c4, c5)
            integer*4        c1 [REFERENCE] 
            integer*4        c2 [REFERENCE] 
            integer*4        c3 [REFERENCE] 
            real*16          c4 [REFERENCE] 
            integer*4        c5 [REFERENCE]
        end
    end interface

    interface pmGetSolutionData
        subroutine pm_get_dataA [ALIAS:'?pm_get_data@@YAXPAHPAM110@Z'] (a1, a2, a3, a4, a5)
            integer*4        a1 [REFERENCE] 
            real*4           a2 [REFERENCE] 
            real*4           a3 [REFERENCE] 
            real*4           a4 [REFERENCE] 
            integer*4        a5 [REFERENCE]
        end

        subroutine pm_get_dataB [ALIAS:'?pm_get_data@@YAXPAHPAN110@Z'] (b1, b2, b3, b4, b5)
            integer*4        b1 [REFERENCE] 
            real*8           b2 [REFERENCE] 
            real*8           b3 [REFERENCE] 
            real*8           b4 [REFERENCE] 
            integer*4        b5 [REFERENCE]
        end

        subroutine pm_get_dataC [ALIAS:'?pm_get_data@@YAXPAHPAO110@Z'] (c1, c2, c3, c4, c5)
            integer*4        c1 [REFERENCE] 
            real*16          c2 [REFERENCE] 
            real*16          c3 [REFERENCE] 
            real*16          c4 [REFERENCE] 
            integer*4        c5 [REFERENCE]
        end
    end interface

    interface pmSetPH
        subroutine pm_set_phA [ALIAS: '?pm_set_ph@@YAXPAH0PAM0@Z'] (a1, a2, a3, a4)
            integer*4        a1 [REFERENCE]
            integer*4        a2 [REFERENCE]
            real*4           a3 [REFERENCE]
            integer*4        a4 [REFERENCE]
        end

        subroutine pm_set_phB [ALIAS: '?pm_set_ph@@YAXPAH0PAN0@Z'] (b1, b2, b3, b4)
            integer*4        b1 [REFERENCE]
            integer*4        b2 [REFERENCE]
            real*8           b3 [REFERENCE]
            integer*4        b4 [REFERENCE]
        end

        subroutine pm_set_phC [ALIAS: '?pm_set_ph@@YAXPAH0PAO0@Z'] (c1, c2, c3, c4)
            integer*4        c1 [REFERENCE]
            integer*4        c2 [REFERENCE]
            real*16          c3 [REFERENCE]
            integer*4        c4 [REFERENCE]
        end
    end interface

    interface pmSetPE
        subroutine pm_set_peA [ALIAS: '?pm_set_pe@@YAXPAH0PAM0@Z'] (a1, a2, a3, a4)
            integer*4        a1 [REFERENCE]
            integer*4        a2 [REFERENCE]
            real*4           a3 [REFERENCE]
            integer*4        a4 [REFERENCE]
        end
        
        subroutine pm_set_peB [ALIAS: '?pm_set_pe@@YAXPAH0PAN0@Z'] (b1, b2, b3, b4)
            integer*4        b1 [REFERENCE]
            integer*4        b2 [REFERENCE]
            real*8           b3 [REFERENCE]
            integer*4        b4 [REFERENCE]
        end

        subroutine pm_set_peC [ALIAS: '?pm_set_pe@@YAXPAH0PAO0@Z'] (c1, c2, c3, c4)
            integer*4        c1 [REFERENCE]
            integer*4        c2 [REFERENCE]
            real*16          c3 [REFERENCE]
            integer*4        c4 [REFERENCE]
        end
    end interface

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
            call pmStart (Me%PhreeqCInstanceID, STAT_) 
            if (STAT_ .EQ. 0) call EndWithError('Subroutine StartPhreeqC; ModulePhreeqC. ERR010.')
                       
            !Read ModulePhreeqC Input File  
            call ReadInputFile (FileName)
            
            !Returns ID
            PhreeqCID = Me%InstanceID            

            STAT_ = SUCCESS_
            
        else  cd0
            
            call EndWithError ('Subroutine StartPhreeqC; ModulePhreeqC. ERR020.')

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
         
        !ToDo: Check if this is necessary
        Me%CalcData%MassOfAllSolutes   = 0
        Me%CalcData%VolumeOfAllSolutes = 0
        Me%CalcData%MassOfWater        = 0
        Me%CalcData%VolumeOfWater      = 0
        Me%CalcData%MassOfSolution     = 0
        
        !----------------------------------------------------------------------

    end subroutine InitializeInstance    
    !----------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    subroutine ReadInputFile (FileName)
     
        !Arguments-------------------------------------------------------------
        character(LEN = *)             :: FileName      
        
        !Local-----------------------------------------------------------------
        integer :: STAT_
        
        !Begin-----------------------------------------------------------------
        !Associate EnterData
        call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_)
        if (STAT_ .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadInputFile; ModulePhreeqC. ERR010.') 
                
        call ReadPhreeqCOptions   
        call ReadPhreeqCProperties             
        
        call KillEnterData(Me%ObjEnterData, STAT = STAT_) 
        if (STAT_ .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadInputFile; ModulePhreeqC. ERR020.')

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

        call GetData(Me%MOptions%Database,           &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'DATABASE',      &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR010.')
        if (flag .EQ. 0) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR020.')
        
            
        call GetData(Me%MOptions%DatabaseAux,        &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'DATABASE_AUX',  &
                     default      = '',              &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR030.') 

        call GetData(Me%MOptions%PrintInput,         &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'PRINT_INPUT',   &
                     default      = .false.,         &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR040.') 
        
        call GetData(Me%MOptions%PrintOutput,        &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'PRINT_OUTPUT',  &
                     default      = .false.,         &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR050.') 

        call GetData(Me%MOptions%HPlusDensity,       &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'HPLUS_DENSITY', &
                     default      = 0.09,            &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR060.')
        
        call GetData(Me%MOptions%WaterDensity,       &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'WATER_DENSITY', &
                     default      = 998.0,           &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR070.') 
               
        call GetData(Me%MOptions%Redox%Element1,            & 
                     Me%ObjEnterData, flag,                 &
                     SearchType   = FromFile,               &
                     keyword      = 'REDOX_PAIR_ELEMENT_1', & 
                     ClientModule = 'ModulePhreeqC',        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR080.')  
        
        if (flag .NE. 0) then
        
            call GetData(Me%MOptions%Redox%Valence1,            & 
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'REDOX_PAIR_VALENCE_1', & 
                         ClientModule = 'ModulePhreeqC',        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. flag .EQ. 0) &
                call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR090.')                      

            call GetData(Me%MOptions%Redox%Element2,            & 
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'REDOX_PAIR_ELEMENT_2', & 
                         ClientModule = 'ModulePhreeqC',        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. flag .EQ. 0) &
                call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR100.')                      

            ! Valence 2 ------------------------------------------------------
            call GetData(Me%MOptions%Redox%Valence2,            & 
                         Me%ObjEnterData, flag,                 &
                         SearchType   = FromFile,               &
                         keyword      = 'REDOX_PAIR_VALENCE_2', & 
                         ClientModule = 'ModulePhreeqC',        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. flag .EQ. 0) &
                call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR110.')

            call pm_solution_redox(Me%PhreeqCInstanceID,                      &
                                   trim(Me%MOptions%Redox%Element1)//char(0), &
                                   Me%MOptions%Redox%Valence1,                &
                                   trim(Me%MOptions%Redox%Element2)//char(0), &
                                   Me%MOptions%Redox%Valence2,                &
                                   STAT_CALL)
            if (STAT_CALL .EQ. 0) &
                call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR120.')

        end if        
       
        call GetData(Me%MOptions%pHCharge,           &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'PH_CHARGE',     &
                     default      = 0,               &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR130.')
        
        call GetData(Me%MOptions%pECharge,           &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      = 'PE_CHARGE',     &
                     default      = 0,               &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR140.')

        if ((Me%MOptions%pECharge .EQ. 1) .AND. (Me%MOptions%pHCharge .EQ. 1)) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR150.')
        
        call GetData(Me%MOptions%DTSeconds,          &
                     Me%ObjEnterData, flag,          &
                     SearchType   = FromFile,        &
                     keyword      ='DTSECONDS',      & 
                     default      = 3600.,           & 
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR160.')

        call GetData(Me%MOptions%Debug,               &
                     Me%ObjEnterData, flag,           &
                     SearchType   = FromFile,         &
                     keyword      = 'DEBUG',          &
                     default      = .false.,          &
                     ClientModule = 'ModulePhreeqC',  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR170.') 

        call GetData(Me%MOptions%PrintAllways,        &
                     Me%ObjEnterData, flag,           &
                     SearchType   = FromFile,         &
                     keyword      = 'PRINT_ALLWAYS',  &
                     default      = .false.,          &
                     ClientModule = 'ModulePhreeqC',  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) &
            call EndWithError ('Subroutine ReadPhreeqCOptions; ModulePhreeqC. ERR180.') 

cd1:   if (flag .EQ. 0) then
            write(*,*) 
            write(*,*) 'Keyword DTSECONDS not found in PhreeqC data file.'
            write(*,*) 'Subroutine ReadPhreeqCOptions; ModulePhreeqC. WRN010.'
            write(*,*) 'Assumed ', Me%MOptions%DTSeconds, 'seconds (',  Me%MOptions%DTSeconds / 60.0, 'hour).'
            write(*,*) 
        end if cd1
        
        !For compatibility with the rest of the program,  
        Me%MOptions%DTDay = Me%MOptions%DTSeconds / 24.0 / 60.0 / 60.0
                
        call OpenDebugFile
        
        call ReadPhreeqCDatabase 

        !----------------------------------------------------------------------

    end subroutine ReadPhreeqCOptions         
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine ReadPhreeqCDatabase

        !Local ----------------------------------------------------------------
        integer :: status

        !----------------------------------------------------------------------

        call pmLoadDatabase (Me%PhreeqCInstanceID, trim(Me%MOptions%Database)//char(0), status)      
        if (status .EQ. 0) call EndWithError ('Subroutine PhreeqCReadDatabase; ModulePhreeqC. ERR010.')
            
        if (Me%MOptions%DatabaseAux .NE. '') then
        
            call pmLoadDatabase (Me%PhreeqCInstanceID, trim(Me%MOptions%DatabaseAux)//char(0), status)
            if (status .EQ. 0) call EndWithError ('Subroutine PhreeqCReadDatabase; ModulePhreeqC. ERR020.')
        
        end if
            
        call pmSetupModel (Me%PhreeqCInstanceID, status)      
        if (status .EQ. 0) call EndWithError ('Subroutine PhreeqCReadDatabase; ModulePhreeqC. ERR030.')
        !----------------------------------------------------------------------

    end subroutine ReadPhreeqCDatabase
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    subroutine ReadPhreeqCProperties

        !External--------------------------------------------------------------
        integer                             :: ClientNumber
        integer                             :: STAT_CALL
        logical                             :: BlockFound

        !Local-----------------------------------------------------------------
        type (T_PhreeqCProperty), pointer   :: NewProperty
        type (T_PhreeqCProperty), pointer   :: PropertyX 
        integer                             :: Index       

        !----------------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                 &
                                        ClientNumber = ClientNumber,     &
                                        block_begin  = prop_block_begin, &
                                        block_end    = prop_block_end,   &
                                        BlockFound   = BlockFound,       &
                                        STAT         = STAT_CALL)
cd1 :       if (STAT_CALL .EQ. SUCCESS_) then    

cd2 :           if (BlockFound) then                                                  
                    
                    !Construct a New Property 
                    Call ConstructProperty (NewProperty)

                    !Add new Property to the Properties List 
                    Call AddProperty (NewProperty)

                else cd2

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        call EndWithError ('ReadPhreeqCProperties - ModulePhreeqC - ERR010')
                    exit do1    !No more blocks
                
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                call EndWithError ('ReadPhreeqCProperties - ModulePhreeqC - ERR020')
            
            else cd1
                
                call EndWithError ('ReadPhreeqCProperties - ModulePhreeqC - ERR030')
            
            end if cd1
        
        end do do1

        nullify (Me%PropertyValues)
        allocate (Me%PropertyValues(Me%PropertiesNumber))
        
        Index = 0
        PropertyX => Me%FirstProperty
        do while (associated(PropertyX))
            PropertyX%PropertyID = Index
            Index = Index + 1
            PropertyX => PropertyX%Next
        end do


        !----------------------------------------------------------------------
    
    end subroutine ReadPhreeqCProperties
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine ConstructProperty (NewProperty)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer    :: NewProperty

        !External--------------------------------------------------------------
        integer                             :: STAT_CALL

        !----------------------------------------------------------------------
             
        allocate (NewProperty, STAT = STAT_CALL)            
        if(STAT_CALL .NE. SUCCESS_) call EndWithError ('ConstructProperty - ModulePhreeqC - ERR010')
        
        nullify(NewProperty%Prev, NewProperty%Next)

        call ConstructPropertyID (NewProperty%ID, Me%ObjEnterData, FromBlock)
        call ReadChemistryParameters (NewProperty)
        call SetPhreeqCProperty (NewProperty)
        !----------------------------------------------------------------------
    
    end subroutine ConstructProperty
    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------
    subroutine ReadChemistryParameters (NewProperty)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: NewProperty
        
        !Local-----------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: iflag
    
        !----------------------------------------------------------------------
       
        call GetData(NewProperty%Params%Group,       &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_GROUP', &
                     Default      = 0,               &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryParameters - ModulePhreeqC - ERR010')
        
        call GetData(NewProperty%Params%PhreeqCName, &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_NAME',  &
                     Default      = '',              &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_ .OR. (iflag .EQ. 0 .AND. NewProperty%Params%Group .NE. 0)) &
            call EndWithError ('ReadChemistryParameters - ModulePhreeqC - ERR020')
        
        call GetData(NewProperty%Params%DoNotChange,         &
                     Me%ObjEnterData, iflag,                 &
                     SearchType   = FromBlock,               &
                     keyword      = 'PHREEQC_DO_NOT_CHANGE', &
                     Default      = 0,                       &
                     ClientModule = 'ModulePhreeqC',         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryParameters - ModulePhreeqC - ERR030')        
        
        call GetData(NewProperty%Params%Debug,       &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_DEBUG', &
                     Default      = .false.,         &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryParameters - ModulePhreeqC - ERR040')
        
        
!        call GetData(NewProperty%Params%Charge,                    &
!                     Me%ObjEnterData, iflag,                       &
!                     SearchType   = FromBlock,                     &
!                     keyword      = 'PHREEQC_CHARGE',              &
!                     Default      = 0,                             &
!                     ClientModule = 'ModulePorousMediaProperties', &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryParameters - ModulePhreeqC - ERR004')
               
        select case (NewProperty%Params%Group)
        case (OTHER) !NO GROUP - pH, pE, Temperature, Density, etc
            !Do nothing
        case (CONCENTRATION) !CONCENTRATION - Only Concentration properties
            call ReadChemistryConcGroupParam (NewProperty)
        case (PHASE) !PHASES
            call ReadChemistryPhasesGroupParam (NewProperty)
!        case (3) !SOLID PHASE
!            call ReadChemistrySolidGroupParam (NewProperty)
!        case (4) !GAS PHASE
!            call ReadChemistryGasGroupParam (NewProperty)
!        case (SURFACE) !SURFACE
!            call ReadChemistrySurfGroupParam (NewProperty)
        case (SPECIES) !SPECIES
            call ReadChemistrySpeciesGroupParam (NewProperty)
        case (EXCHANGE) !EXCHANGE
            call ReadChemistryExcGroupParam (NewProperty)
        case default
            !ToDo: If the group does not exist, must raise an exception and warn the user
        end select
        
    end subroutine ReadChemistryParameters
    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------
    subroutine ReadChemistryConcGroupParam (NewProperty)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: NewProperty
        
        !Local-----------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: iflag
    
        !----------------------------------------------------------------------
        call GetData(NewProperty%Params%Charge,       &
                     Me%ObjEnterData, iflag,          &
                     SearchType   = FromBlock,        &
                     keyword      = 'PHREEQC_CHARGE', &
                     Default      = 0,                &
                     ClientModule = 'ModulePhreeqC',  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryConcGroupParam - ModulePhreeqC - ERR010')

        call GetData(NewProperty%Params%GFW,         &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_GFW',   &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryConcGroupParam - ModulePhreeqC - ERR020')
        if (iflag .EQ. 0) then
            NewProperty%Params%UseGFW = 0
        else              
            NewProperty%Params%UseGFW = 1
        endif
        
        if (NewProperty%Params%UseGFW .EQ. 0) then
            call GetData(NewProperty%Params%As,          &
                         Me%ObjEnterData, iflag,         &
                         SearchType   = FromBlock,       &
                         keyword      = 'PHREEQC_AS',    &
                         ClientModule = 'ModulePhreeqC', &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryConcGroupParam - ModulePhreeqC - ERR030')
            if (iflag .EQ. 0) then
                NewProperty%Params%UseAs = 0
            else              
                NewProperty%Params%UseAs = 1
            end if
        end if
        
        call GetData(NewProperty%Params%RedoxPair%Element1,    & 
                     Me%ObjEnterData, iflag,                   &
                     SearchType   = FromBlock,                 &
                     keyword      = 'PHREEQC_REDOX_ELEMENT_1', & 
                     ClientModule = 'ModulePhreeqC',           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryConcGroupParam; ModulePhreeqC. ERR040.')          
        if (iflag .NE. 0) then        
            call GetData(NewProperty%Params%RedoxPair%Valence1,    & 
                         Me%ObjEnterData, iflag,                   &
                         SearchType   = FromBlock,                 &
                         keyword      = 'PHREEQC_REDOX_VALENCE_1', & 
                         ClientModule = 'ModulePhreeqC',           &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0) call EndWithError ('ReadChemistryConcGroupParam; ModulePhreeqC. ERR050.')

            call GetData(NewProperty%Params%RedoxPair%Element2,    & 
                         Me%ObjEnterData, iflag,                   &
                         SearchType   = FromBlock,                 &
                         keyword      = 'PHREEQC_REDOX_ELEMENT_2', & 
                         ClientModule = 'ModulePhreeqC',           &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0) call EndWithError ('ReadChemistryConcGroupParam; ModulePhreeqC. ERR060.')

            ! Valence 2 ------------------------------------------------------
            call GetData(NewProperty%Params%RedoxPair%Valence2,    & 
                         Me%ObjEnterData, iflag,                   &
                         SearchType   = FromBlock,                 &
                         keyword      = 'PHREEQC_REDOX_VALENCE_2', & 
                         ClientModule = 'ModulePhreeqC',           &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0) call EndWithError ('ReadChemistryConcGroupParam; ModulePhreeqC. ERR070.')
            NewProperty%Params%UseRedox = 1
        else
            NewProperty%Params%UseRedox = 0        
        end if        
        
        call GetData(NewProperty%Params%PhaseName,        &
                     Me%ObjEnterData, iflag,              &
                     SearchType   = FromBlock,            &
                     keyword      = 'PHREEQC_PHASE_NAME', &
                     ClientModule = 'ModulePhreeqC',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryConcGroupParam - ModulePhreeqC- ERR080')
        if (iflag .NE. 0) then
            NewProperty%Params%usePhase = 1
            
            call GetData(NewProperty%Params%SI,          &
                         Me%ObjEnterData, iflag,         &
                         SearchType   = FromBlock,       &
                         keyword      = 'PHREEQC_SI',    &
                         default      = 0.0,             &
                         ClientModule = 'ModulePhreeqC', &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryConcGroupParam - ModulePhreeqC - ERR090')
 
        else
            NewProperty%Params%UsePhase = 0
        end if
       
!        call GetData(NewProperty%Params%Density,       &
!                     Me%ObjEnterData, iflag,           &
!                     SearchType   = FromBlock,         &
!                     keyword      = 'PHREEQC_DENSITY', &
!                     ClientModule = 'ModulePhreeqC',   &
!                     STAT         = STAT_CALL)
!        if (STAT_CALL .NE. SUCCESS_ .OR. iflag .EQ. 0) call EndWithError ('ReadChemistryConcGroupParam - ModulePhreeqC - ERR100')
         
        !----------------------------------------------------------------------

    end subroutine ReadChemistryConcGroupParam   
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    subroutine ReadChemistryPhasesGroupParam (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: NewProperty
        
        !Local-----------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: iflag
    
        !----------------------------------------------------------------------
        call GetData(NewProperty%Params%PhaseType,        &
                     Me%ObjEnterData, iflag,              &
                     SearchType   = FromBlock,            &
                     keyword      = 'PHREEQC_PHASE_TYPE', &
                     Default      = 1,                    &
                     ClientModule = 'ModulePhreeqC',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR010')
                
        call GetData(NewProperty%Params%SI,          &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_SI',    &
                     Default      = 0.0,             &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR020')
        if (NewProperty%Params%SI .NE. 0.0) then            
            call GetData(NewProperty%Params%AlternativeFormula,        &
                         Me%ObjEnterData, iflag,                       &
                         SearchType   = FromBlock,                     &
                         keyword      = 'PHREEQC_ALTERNATIVE_FORMULA', &
                         ClientModule = 'ModulePhreeqC',               &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR030')
            if (iflag .NE. 0) then
                NewProperty%Params%UseAlternativeFormula = 1 
                NewProperty%Params%UseAlternativePhase   = 0
            else
                NewProperty%Params%UseAlternativeFormula = 0
                call GetData(NewProperty%Params%AlternativePhase,        &
                             Me%ObjEnterData, iflag,                     &
                             SearchType   = FromBlock,                   &
                             keyword      = 'PHREEQC_ALTERNATIVE_PHASE', &
                             ClientModule = 'ModulePhreeqC',             &
                             STAT         = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR040')
                if (iflag .NE. 0) then
                    NewProperty%Params%UseAlternativePhase = 1
                else
                    NewProperty%Params%UseAlternativePhase = 0
                endif
            endif
        endif
        
        call GetData(NewProperty%Params%GFW,         &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_GFW',   &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR050')
        if (iflag .EQ. 0) then
            NewProperty%Params%UseGFW = 0
        else              
            NewProperty%Params%UseGFW = 1
        endif
        
        call GetData(NewProperty%Params%ForceEquality,        &
                     Me%ObjEnterData, iflag,                  &
                     SearchType   = FromBlock,                &
                     keyword      = 'PHREEQC_FORCE_EQUALITY', &
                     Default      = 0,                        &
                     ClientModule = 'ModulePhreeqC',          &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR060')
        
        call GetData(NewProperty%Params%DissolveOnly,        &
                     Me%ObjEnterData, iflag,                 &
                     SearchType   = FromBlock,               &
                     keyword      = 'PHREEQC_DISSOLVE_ONLY', &
                     Default      = 0,                       &
                     ClientModule = 'ModulePhreeqC',         &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryPhasesGroupParam - ModulePhreeqC - ERR070')

        !----------------------------------------------------------------------
        
    end subroutine ReadChemistryPhasesGroupParam
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
!    subroutine ReadChemistrySolidGroupParam (NewProperty)
!    
!        !Arguments-------------------------------------------------------------
!        type(T_PhreeqCProperty), pointer :: NewProperty
!        
!        !Local-----------------------------------------------------------------
!        integer                   :: STAT_CALL
!    
!        !----------------------------------------------------------------------
!
!    end subroutine ReadChemistrySolidGroupParam    
!    !--------------------------------------------------------------------------
!
!
!    !--------------------------------------------------------------------------
!    subroutine ReadChemistryGasGroupParam (NewProperty)
!
!        !Arguments-------------------------------------------------------------
!        type(T_PhreeqCProperty), pointer :: NewProperty
!        
!        !Local-----------------------------------------------------------------
!        integer                   :: STAT_CALL
!    
!        !----------------------------------------------------------------------
!
!    end subroutine ReadChemistryGasGroupParam    
!    !--------------------------------------------------------------------------
!
!
!    !--------------------------------------------------------------------------
!    subroutine ReadChemistrySurfGroupParam (NewProperty)
!
!        !Arguments-------------------------------------------------------------
!        type(T_PhreeqCProperty), pointer :: NewProperty
!        
!        !Local-----------------------------------------------------------------
!        integer                   :: STAT_CALL
!    
!        !----------------------------------------------------------------------
!
!    endsubroutine ReadChemistrySurfGroupParam
!    !--------------------------------------------------------------------------
!
!
    !--------------------------------------------------------------------------
    subroutine ReadChemistrySpeciesGroupParam (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: NewProperty
        
        !Local-----------------------------------------------------------------
        integer                   :: STAT_CALL
        integer                   :: iflag
    
        !----------------------------------------------------------------------

        call GetData(NewProperty%Params%GFW,         &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_GFW',   &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistrySpeciesGroupParam - ModulePhreeqC - ERR010')
        if (iflag .EQ. 0) then
            NewProperty%Params%UseGFW = 0
        else              
            NewProperty%Params%UseGFW = 1
        endif

        !----------------------------------------------------------------------

    end subroutine
    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------
    subroutine ReadChemistryExcGroupParam (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: NewProperty
        
        !Local-----------------------------------------------------------------
        integer :: STAT_CALL
        integer :: iflag
    
        !----------------------------------------------------------------------
        call GetData(NewProperty%Params%PhaseName,        &
                     Me%ObjEnterData, iflag,              &
                     SearchType   = FromBlock,            &
                     keyword      = 'PHREEQC_PHASE_NAME', &
                     ClientModule = 'ModulePhreeqC',      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryExcGroupParam - ModulePhreeqC - ERR001')
        if (iflag .NE. 0) then
            NewProperty%Params%UsePhase   = 1
            NewProperty%Params%UseKinetic = 0
        else
             NewProperty%Params%UsePhase = 0
           
            call GetData(NewProperty%Params%KineticName,        &
                         Me%ObjEnterData, iflag,                &
                         SearchType   = FromBlock,              &
                         keyword      = 'PHREEQC_KINETIC_NAME', &
                         ClientModule = 'ModulePhreeqC',        &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryExcGroupParam - ModulePhreeqC - ERR020')
            if (iflag .NE. 0) then
                NewProperty%Params%UseKinetic = 1 
            else
                NewProperty%Params%UseKinetic = 0
            endif
        endif
             
        call GetData(NewProperty%Params%GFW,         &
                     Me%ObjEnterData, iflag,         &
                     SearchType   = FromBlock,       &
                     keyword      = 'PHREEQC_GFW',   &
                     ClientModule = 'ModulePhreeqC', &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) call EndWithError ('ReadChemistryExcGroupParam - ModulePhreeqC - ERR030')
        if (iflag .EQ. 0) then
            NewProperty%Params%UseGFW = 0
        else              
            NewProperty%Params%UseGFW = 1
        endif                         
        !----------------------------------------------------------------------
        
    end subroutine ReadChemistryExcGroupParam    
    !--------------------------------------------------------------------------
    
    
    subroutine SetPhreeqCProperty (Property)
        
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty)        :: Property
 
        !----------------------------------------------------------------------  
        select case (Property%Params%Group)
            case (0) 
                !Do nothing. 
            case (SPECIES)
                call SetSpeciesProperty(Property)                    
            case (CONCENTRATION)
                call SetConcentrationProperty(Property)
            case (PHASE)
                call SetPhaseProperty(Property)
            case (EXCHANGE)
                call SetExchangeProperty(Property)
            case default
                call EndWithError ('Subroutine SetPhreeqCProperty; ModulePhreeqC. ERR010.') !Group not recognized
        end select
        !----------------------------------------------------------------------
        
    end subroutine SetPhreeqCProperty   
    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------
    subroutine SetSpeciesProperty (Property)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty) :: Property

        !Local-----------------------------------------------------------------
        integer :: status
        
        !----------------------------------------------------------------------

        Property%PhreeqCInputID = -1
        call pm_get_species_index(Me%PhreeqCInstanceID,                       &
                                  trim(Property%Params%PhreeqCName)//char(0), &
                                  Property%PhreeqCResultID,                   &
                                  status)
        if (status .EQ. 0) &
            call EndWithError ('Subroutine SetSpeciesProperty; ModulePhreeqC. ERR010.')                      
                             
        !----------------------------------------------------------------------                    
    
    end subroutine SetSpeciesProperty
    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------
    subroutine SetConcentrationProperty (Property)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty) :: Property

        !Local-----------------------------------------------------------------
        real    :: zero = 0.0
        real    :: database_gfw
        integer :: status
        
        !----------------------------------------------------------------------

        !Pass to PhreeqC Object       
        call pm_conc_add(Me%PhreeqCInstanceID, trim(Property%Params%PhreeqCName)//char(0), zero, status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR010.')                      
                    
        !Units are required to be mg/L for solution concentrations
        Property%Params%UseUnits = 0
        call pm_conc_use(Me%PhreeqCInstanceID, Property%Params%Charge, Property%Params%UsePhase, Property%Params%UseAs, &
                         Property%Params%UseGFW, Property%Params%UseUnits, Property%Params%UseRedox, status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR020.')                      
                    
        if (Property%Params%UseAs .EQ. 1) call pm_conc_as(Me%PhreeqCInstanceID, trim(Property%Params%As)//char(0), status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR030.')                      
                    
        if (Property%Params%UseGFW .EQ. 1) call pm_conc_gfw(Me%PhreeqCInstanceID, Property%Params%GFW, status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR040.')                      
                           
        if (Property%Params%UsePhase .EQ. 1) call pm_conc_phase(Me%PhreeqCInstanceID, trim(Property%Params%PhaseName)//char(0), & 
                                                         Property%Params%SI, status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR050.')                      
                    
        if (Property%Params%UseRedox .EQ. 1) call pm_conc_redox(Me%PhreeqCInstanceID,                       &
                                                         trim(Property%Params%RedoxPair%Element1)//char(0), &
                                                         Property%Params%RedoxPair%Valence1,                &
                                                         trim(Property%Params%RedoxPair%Element2)//char(0), &
                                                         Property%Params%RedoxPair%Valence2,                &   
                                                         status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR060.')                      
                                                                                  
        !Now, save the solution property in the definitive structure
        call pm_conc_save (Me%PhreeqCInstanceID, Property%PhreeqCInputID, Property%PhreeqCResultID, database_gfw, status)
        if (status .EQ. 0) call EndWithError ('Subroutine SetSolutionProperty; ModulePhreeqC. ERR070.')  
        
        if (Property%Params%UseGFW .EQ. 0) Property%Params%GFW = database_gfw                   
         
        !----------------------------------------------------------------------                    

    end subroutine SetConcentrationProperty
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine SetPhaseProperty (Property)
    
        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty) :: Property

        !Local-----------------------------------------------------------------
        real                    :: zero = 0.0
        real                    :: SI
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
       
        if (Property%Params%PhaseType .EQ. SOLID_PHASE) then
            Me%MOptions%UseSolidPhase = 1
            SI = Property%Params%SI
        else
            Me%MOptions%UseGasPhase = 1
            SI = log10(Property%Params%SI)
        end if
       
        !Pass to PhreeqC Object
        call pm_ppa_pp(Me%PhreeqCInstanceID,                       &
                       trim(Property%Params%PhreeqCName)//char(0), &
                       trim(Alternative)//char(0),                 &
                       SI,                                         &
                       zero,                                       &
                       Property%Params%ForceEquality,              &
                       Property%Params%DissolveOnly,               &
                       Property%PhreeqCInputID,                    &
                       STAT_)                       
        if (STAT_ .EQ. 0) call EndWithError ('Subroutine SetPhaseProperty; ModulePhreeqC. ERR010.') 
            
        Property%PhreeqCResultID = Property%PhreeqCInputID   
        
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
        call pm_exa_exchanger(Me%PhreeqCInstanceID,                       &
                              trim(Property%Params%PhreeqCName)//char(0), &
                              Property%Params%ExType,                     &
                              trim(Formula)//char(0),                     &
                              zero,                                       &
                              Property%PhreeqCInputID,                    &
                              STAT_)
        if (STAT_ .EQ. 0) &
            call EndWithError ('Subroutine SetExchangeProperty; ModulePhreeqC. ERR010.')
         
        call pm_get_species_index(Me%PhreeqCInstanceID,                       &
                                  trim(Property%Params%PhreeqCName)//char(0), &
                                  Property%PhreeqCResultID,                   &
                                  STAT_) 
        if (STAT_ .EQ. 0) &
            call EndWithError ('Subroutine SetExchangeProperty; ModulePhreeqC. ERR020.') 
        
        Me%MOptions%UseExchanger = 1       
        !----------------------------------------------------------------------
    
    end subroutine SetExchangeProperty
    !--------------------------------------------------------------------------
    
        
    !--------------------------------------------------------------------------
    subroutine AddProperty (NewProperty)

        !Arguments-------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: NewProperty

        !----------------------------------------------------------------------

        ! Add to the Property List a new property
        if (.NOT. associated(Me%FirstProperty)) then
            Me%PropertiesNumber     = 1
            Me%FirstProperty        => NewProperty
            Me%LastProperty         => NewProperty
        else
            NewProperty%Prev        => Me%LastProperty
            Me%LastProperty%Next    => NewProperty
            Me%LastProperty         => NewProperty
            Me%PropertiesNumber     = Me%PropertiesNumber + 1
        end if 
        
        Me%LastProperty%Index = Me%PropertiesNumber 
        !----------------------------------------------------------------------

    end subroutine AddProperty 
    !--------------------------------------------------------------------------
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
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

            if (present(DTDay   )) DTDay    = Me%MOptions%DTDay
            if (present(DTSecond)) DTSecond = Me%MOptions%DTSeconds

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetPhreeqCDT 
    !--------------------------------------------------------------------------    
    
    
    !--------------------------------------------------------------------------    
    subroutine GetPhreeqCNeedSoilDryDensity (PhreeqCID, NeedSoilDryDensity, STAT)
        !Arguments-------------------------------------------------------------
        integer                        :: PhreeqCID
        logical,           intent(OUT) :: NeedSoilDryDensity
        integer, optional, intent(OUT) :: STAT

        !External--------------------------------------------------------------
        integer :: ready_              

        !Local-----------------------------------------------------------------
        integer :: STAT_ !Auxiliar local variable
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            if (Me%MOptions%UseExchanger .EQ. 1 .OR. Me%MOptions%UseSolidPhase .EQ. 1) then
                NeedSoilDryDensity = .true.
            else
                NeedSoilDryDensity = .false.
            endif            

            STAT_ = SUCCESS_
        else 
                
            STAT_ = ready_
            
        end if cd1

        if (present(STAT))  STAT = STAT_

        !----------------------------------------------------------------------
    
    end subroutine GetPhreeqCNeedSoilDryDensity
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
       
        type(T_PhreeqCProperty), pointer :: PropertyX
        
        !----------------------------------------------------------------------
        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)    
        
cd1:    if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then

            found = .false.
            
            PropertyX => Me%FirstProperty
            do while (associated(PropertyX))

                if (PropertyIDNumber .EQ. PropertyX%ID%IDNumber)then
                
                    PropertyIndex = PropertyX%Index
                    found         = .true.
                    exit
                    
                end if
            
                PropertyX => PropertyX%Next
            
            end do

            if(.NOT. found) then
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
    
    
!    !--------------------------------------------------------------------------    
!    subroutine GetPhreeqCOptions (PhreeqCID, MOptions, STAT)
!
!        !Arguments-------------------------------------------------------------
!        integer                                       :: PhreeqCID
!        type(T_PhreeqCOptions),           intent(OUT) :: MOptions
!        integer,                optional, intent(OUT) :: STAT
! 
!        !External--------------------------------------------------------------
!        integer :: ready_              
!
!        !Local-----------------------------------------------------------------
!        integer :: STAT_              !Auxiliar local variable
!        
!        !----------------------------------------------------------------------
!        STAT_ = UNKNOWN_
!
!        call Ready(PhreeqCID, ready_)    
!        
!cd1:    if ((ready_ .EQ. IDLE_ERR_) .OR. (ready_ .EQ. READ_LOCK_ERR_)) then
!
!            MOptions = Me%MOptions
!                
!            STAT_ = SUCCESS_            
!
!        else 
!        
!            STAT_ = ready_
!            
!        end if cd1
!
!        if (present(STAT)) STAT = STAT_
!        !----------------------------------------------------------------------
!        
!    end subroutine GetPhreeqCOptions
!    !--------------------------------------------------------------------------    

    
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
                              WaterVolume,         &
                              WaterMass,           &
                              Temperature,         &
                              pH,                  &
                              pE,                  &
                              SolidMass,           &
                              CellsArrayLB,        &
                              CellsArrayUB,        &
                              OpenPoints,          &
                              STAT)  

        !Arguments---------------------------------------------------------------
        integer                                     :: PhreeqCID
        real,               pointer, dimension(:,:) :: PropertiesValues
        real,               pointer, dimension(:  ) :: WaterVolume
        real,               pointer, dimension(:  ) :: WaterMass
        real,               pointer, dimension(:  ) :: Temperature
        real,               pointer, dimension(:  ) :: pH
        real,               pointer, dimension(:  ) :: pE
        real,    optional,  pointer, dimension(:  ) :: SolidMass
        integer,            intent(IN)              :: CellsArrayLB, CellsArrayUB        
        integer, optional,  pointer, dimension(:  ) :: OpenPoints
        integer, optional,  intent(OUT)             :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_  
        logical                                     :: CalcPoint        
        integer                                     :: CellIndex
        integer                                     :: ready_ 
        integer                                     :: UsePhase 

        !------------------------------------------------------------------------                         
            
        STAT_ = UNKNOWN_

        call Ready(PhreeqCID, ready_)

cd1 :   if (ready_ .EQ. IDLE_ERR_) then
            
            if ((Me%MOptions%UseSolidPhase .eq. 1) .or. (Me%MOptions%UseGasPhase .eq. 1)) then
                UsePhase = 1
            else
                UsePhase = 0
            endif            
            call pmSetUse(Me%PhreeqCInstanceID,         &
                          UsePhase,                     &
                          Me%MOptions%UseGas,           & 
                          Me%MOptions%UseSolidSolution, &
                          Me%MOptions%UseSurface,       &
                          Me%MOptions%UseExchanger,     &
                          STAT_)
            if (STAT_ .EQ. 0) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR010.')

            Me%Ext%PropertiesValues => PropertiesValues
            if (.NOT. associated(Me%Ext%PropertiesValues)) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR020.')

            Me%Ext%WaterVolume => WaterVolume
            if (.NOT. associated(Me%Ext%Watervolume)) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR030.')
                
            Me%Ext%WaterMass => WaterMass
            if (.NOT. associated(Me%Ext%WaterMass)) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR030.')
                            
            Me%Ext%Temperature => Temperature
            if (.NOT. associated(Me%Ext%Temperature)) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR040.')

            Me%Ext%pH => pH
            if (.NOT. associated(Me%Ext%pH)) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR050.')

            Me%Ext%pE => pE
            if (.NOT. associated(Me%Ext%pE)) &
                call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR060.')            
            
            if (present(SolidMass)) then
                Me%Ext%SolidMass => SolidMass
                if (.NOT. associated(Me%Ext%SolidMass)) &
                    call EndWithError ('Subroutine ModifyPhreeqC; ModulePhreeqC. ERR070.')            
            else
                nullify(Me%Ext%SolidMass)
            end if
            
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
        nullify(Me%Ext%WaterVolume)
        nullify(Me%Ext%WaterMass)
        nullify(Me%Ext%Temperature)
        nullify(Me%Ext%pH)
        nullify(Me%Ext%pE)
        nullify(Me%Ext%SolidMass)

        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine ModifyPhreeqC
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine MakeCalculations(CellIndex)

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex 

        !Local-------------------------------------------------------------------
        integer :: STAT_
        real    :: ph, pe, mass_of_water

        type(T_PhreeqCProperty), pointer :: PropertyX

        !Begin-------------------------------------------------------------------  
        call CalculateSolutionDensity(CellIndex, Me%CalcData%DensityOfSolution)
               
        call pmSetPH(Me%PhreeqCInstanceID, &
                     Me%MOptions%pHCharge, &
                     Me%Ext%pH(CellIndex), &
                     STAT_)
        if (STAT_ .EQ. 0) &
            call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR010.')

        call pmSetPE(Me%PhreeqCInstanceID, &
                     Me%MOptions%pECharge, &
                     Me%Ext%pE(CellIndex), &
                     STAT_)
        if (STAT_ .EQ. 0) &
            call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR020.')

        call pmSetSolutionData(Me%PhreeqCInstanceID,          &
                               Me%Ext%Temperature(CellIndex), &
                               Me%MOptions%Units%Solution,    & 
                               Me%CalcData%DensityOfSolution, & 
                               Me%Ext%WaterMass(CellIndex),   &                                                               
                               STAT_)                                
        if (STAT_ .EQ. 0) &
            call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR030.')                                  
                       
        !Pass to PhreeqCObject all the properties
        PropertyX => Me%FirstProperty
do1:    do while (associated(PropertyX))
          
            select case (PropertyX%Params%Group)
                case (CONCENTRATION, PHASE, GAS, SURFACE, EXCHANGE) 
                    call ConvertInput (CellIndex, PropertyX)
                              
                    call pmSetInputValue (Me%PhreeqCInstanceID,     & 
                                          PropertyX%PhreeqCInputID, &
                                          PropertyX%Params%Group,   &
                                          PropertyX%PropertyValue,  &
                                          STAT_)
                    if (STAT_ .EQ. 0) &
                        call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR040.')
                    
                case (OTHER, SPECIES)
                    !Do nothing
                    
                case default 
                    call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR050.')
                    
            end select
                
            PropertyX => PropertyX%Next    
                
        end do do1
        
        !This will print the input to a file if user turn on the option
        if (Me%MOptions%PrintAllways) then
            call PrintDataInput (CellIndex)        
        endif
        
        call pmRunPhreeqC (Me%PhreeqCInstanceID, STAT_) 
        if (STAT_ .EQ. 0) then
            if (Me%MOptions%Debug) then
                call PrintDataInput (CellIndex)        
            endif
            call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR060.')
        endif
        
        call pmGetSolutionData(Me%PhreeqCInstanceID, mass_of_water, ph, pe, STAT_)
        if (STAT_ .EQ. 0) &
            call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR070.')                    
        
        PropertyX => Me%FirstProperty
do2:    do while (associated(PropertyX))

            if (PropertyX%Params%DoNotChange .NE. 1) then

                call pmGetResultValue (Me%PhreeqCInstanceID,      & 
                                       PropertyX%PhreeqCResultID, &
                                       PropertyX%Params%Group,    &
                                       PropertyX%PropertyValue,   &
                                       STAT_)
                if (STAT_ .EQ. 0) &
                    call EndWithError ('Subroutine MakeCalculations; ModulePhreeqC. ERR080.')

                !PhreeqC return results in mol's
                call ConvertResult (CellIndex, PropertyX)

            end if
            
            PropertyX => PropertyX%Next
                    
        end do do2
                     
        if (Me%MOptions%PrintAllways) then
            call PrintDataOutput (CellIndex)
        endif
                                
        !------------------------------------------------------------------------
       
    end subroutine MakeCalculations
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine CalculateSolutionDensity (CellIndex, SolutionDensity)

        !Argument----------------------------------------------------------------
        integer, intent(IN)  :: CellIndex
        real,    intent(OUT) :: SolutionDensity

        !Local-------------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: PropertyX
        real                             :: Volume
        
        !Begin-------------------------------------------------------------------   
        Me%CalcData%MassOfSolution = Me%Ext%WaterMass(CellIndex)       
        !  L   =             m3                * 1E3
        Volume = Me%Ext%WaterVolume(CellIndex) * 1E3
        
        PropertyX => Me%FirstProperty        
do1:    do while (associated(PropertyX))
            !kg  =  kg  + (                        mg                          * 1E-6)
            Me%CalcData%MassOfSolution = Me%CalcData%MassOfSolution + (Me%Ext%PropertiesValues(PropertyX%Index, CellIndex) * 1E-6)
            PropertyX => PropertyX%Next            
        end do do1             
        
        !     kg/L      =  kg  /   L
        SolutionDensity = Me%CalcData%MassOfSolution / Volume
        !------------------------------------------------------------------------

    end subroutine CalculateSolutionDensity
    !----------------------------------------------------------------------------
    
    
    !----------------------------------------------------------------------------
    subroutine ConvertInput (CellIndex, Property)
        !For now, the input units for each group are fixed.
        !In the future, if possible, this will be made more flexible

        !Argument----------------------------------------------------------------
        integer, intent(IN)              :: CellIndex
        type(T_PhreeqCProperty), pointer :: Property  
        
        !Local------------------------------------------------------------------- 
        character(StringLength) :: Name    

        !------------------------------------------------------------------------       
                                       
        select case (Property%Params%Group)
            case (CONCENTRATION) !The input units for CONCENTRATION concentration properties MUST be mg/L
                Name = Property%Params%PhreeqCName                                           
                if (.NOT. ((Name .EQ. 'H(1)') .OR. (Name .EQ. 'E'))) then                       
                    !mol/kgw = ((mg/L) * 1E-3 / (kgw/L) / (g/mol))
                    Property%PropertyValue = (Me%Ext%PropertiesValues(Property%Index, CellIndex) * 1E-3 / &
                                              DefaultWaterDensity / Property%Params%GFW)                    
                end if  
                
            case (PHASE)
                
                if (Property%Params%PhaseType .EQ. SOLID_PHASE) then !It's a solid pure phase
                    if (.NOT. associated(Me%Ext%SolidMass)) &
                        call EndWithError ('Subroutine ConvertInputs; ModulePhreeqC. ERR010.')
                    
                    Property%PropertyValue = Me%Ext%PropertiesValues(Property%Index, CellIndex)
                else !it's a gas pure phase                 
                    !For now, the "moles" of the gas at disposition are always 10 moles and the volume is "ignored"
                    !Basically, at the partial pressure given there is an "infinite supply" of the gas.
                    !This must be changed in the future...   
                    Property%PropertyValue = Me%Ext%PropertiesValues(Property%Index, CellIndex)
                end if
            
            case (EXCHANGE)

                if (.NOT. associated(Me%Ext%SolidMass)) &
                    call EndWithError ('Subroutine ConvertInputs; ModulePhreeqC. ERR020.')
                !mols = (mg/kgs * kgs / 1000 / (g/mol))            
                Property%PropertyValue = (Me%Ext%PropertiesValues(Property%Index, CellIndex) * &
                                          Me%Ext%SolidMass(CellIndex) / 1000 / Property%Params%GFW)

            case default                                                      

                call EndWithError ('Subroutine ConvertInputs; ModulePhreeqC. ERR030.')
                
        end select
                    
        !------------------------------------------------------------------------       
        
    end subroutine ConvertInput
    !----------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------
    subroutine ConvertResult (CellIndex, Property)
        !For now, the input units for each group are fixed.
        !In the future, if possible, this will be made more flexible

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex                
        type(T_PhreeqCProperty), pointer :: Property  

        !------------------------------------------------------------------------  
                     
        select case (Property%Params%Group)                
            case (CONCENTRATION, SPECIES) 
                      
                !mg/L = mol * (g/mol) * 1000 / (m3 * 1000)
                Me%Ext%PropertiesValues(Property%Index, CellIndex) = Property%PropertyValue * &
                                                                     Property%Params%GFW * 1E3 / &
                                                                     (Me%Ext%WaterVolume(CellIndex) * 1E3)
           
            case (PHASE)
                            
                if (Property%Params%PhaseType .EQ. SOLID_PHASE) then

                    if (.NOT. associated(Me%Ext%SolidMass)) &
                        call EndWithError ('Subroutine ConvertInputs; ModulePhreeqC. ERR010.')                
                    Me%Ext%PropertiesValues(Property%Index, CellIndex) = Property%PropertyValue
                else
                                       
                    !For now, partial pressure is a FIXED parameter and do not change beteen DT's.
                    !Also the number of "moles" and volume do not change (the first is 10 moles and the second is ignored)
                    !Me%Ext%PropertiesValues(Property%Index, CellIndex) = Property%PropertyValue
                                     
                end if
            
            case (EXCHANGE)
            
                    !mg/kgs = (mols * g/mol * 1000 / kgs)
                    Me%Ext%PropertiesValues(Property%Index, CellIndex) = (Property%PropertyValue * &
                                                                          Property%Params%GFW * 1000 / &
                                                                          Me%Ext%SolidMass(CellIndex)) 
                    
            case (OTHER)
                !Do nothing
                
            case default
                
                !ToDo: Put an error message here
                
        end select                             
        
        !------------------------------------------------------------------------       
                
    end subroutine ConvertResult
    !----------------------------------------------------------------------------      
    
    
    !----------------------------------------------------------------------------
    subroutine PrintDataInput(CellIndex)
    
        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex

        !Local-------------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: PropertyX

        !------------------------------------------------------------------------
        if (Me%MOptions%PrintInput) then

            write (Me%DebugFileUnit, *) 'INPUT FOR CELL: ', CellIndex
            write (Me%DebugFileUnit, *) 'Volume of Water: ', Me%Ext%WaterVolume(CellIndex)
            write (Me%DebugFileUnit, *) 'Mass of Water: ', Me%Ext%WaterMass(CellIndex)
            write (Me%DebugFileUnit, *) 'Volume of Solution: ', Me%Ext%WaterVolume(CellIndex)
            write (Me%DebugFileUnit, *) 'Mass of Solution: ', Me%CalcData%MassOfSolution            
            write (Me%DebugFileUnit, *) 'Solution Density: ', Me%CalcData%DensityOfSolution
            write (Me%DebugFileUnit, *) 'Mass of Soil: ', Me%Ext%SolidMass(CellIndex)
            write (Me%DebugFileUnit, *) 'pH: ', Me%Ext%pH(CellIndex)
            write (Me%DebugFileUnit, *) 'pE: ', Me%Ext%pE(CellIndex)
            
            PropertyX => Me%FirstProperty
do1:        do while (associated(PropertyX))

                if (PropertyX%Params%Debug) then
                    
                    write (Me%DebugFileUnit, *) trim(PropertyX%Params%PhreeqCName), ': ', &
                                                Me%Ext%PropertiesValues(PropertyX%Index, CellIndex), &
                                                " (", PropertyX%PropertyValue , ")"
                
                endif
                                                  
                PropertyX => PropertyX%Next    
                    
            end do do1

        endif
        !------------------------------------------------------------------------
    
    end subroutine PrintDataInput
    !----------------------------------------------------------------------------
    

    !----------------------------------------------------------------------------
    subroutine PrintDataOutput (CellIndex)

        !Argument----------------------------------------------------------------
        integer, intent(IN) :: CellIndex

        !Local-------------------------------------------------------------------
        type(T_PhreeqCProperty), pointer :: PropertyX

        !------------------------------------------------------------------------
        if (Me%MOptions%PrintOutput) then


            write (Me%DebugFileUnit, *) 'OUTPUT FOR CELL: ', CellIndex

            PropertyX => Me%FirstProperty
do1:        do while (associated(PropertyX))

                if (PropertyX%Params%Debug) then
                    
                    write (Me%DebugFileUnit, *) trim(PropertyX%Params%PhreeqCName), ': ', &
                                                Me%Ext%PropertiesValues(PropertyX%Index, CellIndex), &
                                                "(", PropertyX%PropertyValue, ")"
                
                endif
                                                  
                PropertyX => PropertyX%Next    
                    
            end do do1

            write (Me%DebugFileUnit, *) '--------------------------------'

        endif
        !------------------------------------------------------------------------

    end subroutine PrintDataOutput
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
                call pmKill(Me%PhreeqCInstanceID, status)
                if (status .EQ. 0) &
                    call EndWithError ('Subroutine KillPhreeqC; ModulePhreeqC. ERR010.')
                
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
    subroutine EndWithError (ErrorMessage)

        !Arguments---------------------------------------------------------------
        character(LEN=*) :: ErrorMessage
    
        !------------------------------------------------------------------------
        call CloseDebugFile
        write (*,*) ErrorMessage
        stop
        
        !------------------------------------------------------------------------

    end subroutine EndWithError
    !----------------------------------------------------------------------------
    
    
    !----------------------------------------------------------------------------
    subroutine OpenDebugFile 
    
        !Local-------------------------------------------------------------------
        integer :: STAT
    
        !------------------------------------------------------------------------
        if (Me%MOptions%PrintInput .OR. Me%MOptions%PrintOutput) then
        
            call UnitsManager (Me%DebugFileUnit, OPEN_FILE, STAT)
            
            if (STAT .NE. SUCCESS_) &
                stop 'Subroutine OpenDebugFile; ModulePhreeqC. ERR010.'            

            Open (UNIT=Me%DebugFileUnit,   &
                  FILE='PhreeqCDebug.txt', &
                  STATUS='REPLACE',        &
                  ACTION='WRITE',          &
                  IOSTAT=STAT)
            
            if (STAT .NE. SUCCESS_) &
                stop 'Subroutine OpenDebugFile; ModulePhreeqC. ERR020.'            

        else
        
            Me%DebugFileUnit = -1
        
        endif
        !------------------------------------------------------------------------
    
    end subroutine OpenDebugFile 
    !----------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    subroutine CloseDebugFile
    
        !Local-------------------------------------------------------------------
        integer :: STAT
    
        !------------------------------------------------------------------------
        if (Me%DebugFileUnit .NE. -1) then

            call UnitsManager (Me%DebugFileUnit, CLOSE_FILE, STAT)
            
            if (STAT .NE. SUCCESS_) &
                stop 'Subroutine CloseDebugFile; ModulePhreeqC. ERR010.'  

        endif        
        !------------------------------------------------------------------------
        
    end subroutine CloseDebugFile
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine Ready (PhreeqCID, ready_) 

        !Arguments---------------------------------------------------------------
        integer :: PhreeqCID
        integer :: ready_
        !------------------------------------------------------------------------

        nullify (Me)

cd1:    if (PhreeqCID > 0) then
            
            call LocateObjPhreeqC (PhreeqCID)
            ready_ = VerifyReadLock (mPHREEQC_, Me%InstanceID)

        else

            ready_ = OFF_ERR_

        end if cd1
        !------------------------------------------------------------------------

    end subroutine Ready
    !----------------------------------------------------------------------------


    !----------------------------------------------------------------------------
    subroutine LocateObjPhreeqC (PhreeqCID)

    !Arguments-------------------------------------------------------------------
        integer :: PhreeqCID
    !----------------------------------------------------------------------------

        Me => FirstObjPhreeqC
        do while (associated (Me))
            if (Me%InstanceID == PhreeqCID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))   &
            call EndWithError ('Subroutine LocateObjPhreeqC; ModulePhreeqC. ERR010.')
    !----------------------------------------------------------------------------
    
    end subroutine LocateObjPhreeqC
    !----------------------------------------------------------------------------


end module ModulePhreeqC    

#endif
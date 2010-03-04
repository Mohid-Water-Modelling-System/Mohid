!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : PhreeqCData
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : February 2010
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

module ModulePhreeqCData
    
    use ModuleGlobalData
    
    implicit none

    public

    !Constants-----------------------------------------------------------------
    integer, parameter :: OTHER         = 0
    integer, parameter :: CONCENTRATION = 1
    integer, parameter :: PHASE         = 2
    integer, parameter :: GAS           = 3
    integer, parameter :: SURFACE       = 4
    integer, parameter :: SPECIES       = 5
    integer, parameter :: EXCHANGE      = 6
    
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
    end type T_ChemistryParameters    
    
    type T_PhreeqCUnits
        integer             :: MasterSpecies = mg_l
        integer             :: Species       = mg_l
        integer             :: Alkalinity    = mg_l
        integer             :: SolidPhases   = g_kgsoil !Check to see if this is reasonable
        integer             :: Exchangers    = ug_kgsoil
    end type T_PhreeqCUnits      
        
    type T_PhreeqCOptions
        character(len=2048)  :: Database                  !Path for database
        character(len=2048)  :: DatabaseAux  = ''         !Path for auxiliary database
        logical              :: PrintInput   = .false.    !For DEBUG    
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
    end type T_PhreeqCOptions



end module ModulePhreeqCData

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006, 2010. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 
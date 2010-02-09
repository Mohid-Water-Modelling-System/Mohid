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
        
        integer                 :: PhreeqCInputID
        integer                 :: PhreeqCResultID
    end type T_ChemistryParameters    

end module ModulePhreeqCData

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. MARETEC, Instituto Superior Técnico, Technical University of Lisbon. 
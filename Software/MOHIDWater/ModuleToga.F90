!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Water
! MODULE        : Toga
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : This module generates returns a water level estimate from a set of tidal components.
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
! This module was adapted from The University of Hawaii (Pat Caldwell) Tidal analysis Package.
! The Hawaii was also ana adaptation of the package develop by Michael Foreman from the 
! Institute of Ocean Sciences in Canada. 
! The original Foreman tidal analysis software is available via:
! http://www.pac.dfo-mpo.gc.ca/sci/osap/projects/tidpack/tidpack_e.htm 
!
!
!------------------------------------------------------------------------------

Module ModuleToga

    use ModuleGlobalData
    use ModuleTime                

    Implicit none

    Private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public :: ConstructToga

    !Selector


    !Modifier
    Public  :: TogaLevel
    private ::      ASTRO
    private ::          TideFreq


    !Destructor
    public  :: KillToga


    !Management
    private ::  Ready
    private ::      LocateObjToga

    !Parameter
    integer, parameter          :: NTidal  =  45
    integer, parameter          :: NCompt  = 146
    integer, parameter          :: NConsH1 = 170
    integer, parameter          :: NConsH2 = 251

    real, parameter             :: DTHora = 1.0 / 3600.0
    real, parameter             :: TwoPI  = 2.0 * PI




    !Type
    type      T_Toga
        integer                                 :: InstanceID = null_int !inicialization: Carina
        real(8)                                 :: DS    = null_real
        real(8)                                 :: DH    = null_real
        real(8)                                 :: DP    = null_real
        real(8)                                 :: DNP   = null_real
        real(8)                                 :: DPP   = null_real
        real(8)                                 :: D1    = null_real
        real(8)                                 :: HH    = null_real
        real(8)                                 :: Tau   = null_real
        real(8)                                 :: Hour  = null_real
        real(8)                                 :: H     = null_real
        real(8)                                 :: S     = null_real
        real(8)                                 :: P     = null_real
        real(8)                                 :: PP    = null_real
        real(8)                                 :: ENP   = null_real
        character(len=5), pointer, dimension(:) :: kon   => null()  !inicialization: Carina
        integer, pointer, dimension(:)          :: indx  => null()  !inicialization: Carina
        real, pointer, dimension(:)             :: Sig, twoc  => null()  !inicialization: Carina
        real, pointer, dimension(:)             :: CH, CHP    => null()  !inicialization: Carina  
        real, pointer, dimension(:)             :: VMare, UMare, FMare, Freq  => null()  !inicialization: Carina
        type (T_Toga), pointer                  :: Next  => null()  !inicialization: Carina
    end type T_Toga

    !Global Module Variables
    type (T_Toga), pointer                      :: FirstToga   => null()  !inicialization: Carina
    type (T_Toga), pointer                      :: Me          => null()  !inicialization: Carina

    !Data----------------------------------------------------------------------
    character(LEN = 5), dimension(NCompt ), parameter  :: KonTB = (/                             &
               'Z0   ','SA   ','SSA  ','MSM  ','MM   ','MSF  ','MF   ','ALP1 ','2Q1  ','SIG1 ',  &
               'Q1   ','RHO1 ','O1   ','TAU1 ','BET1 ','NO1  ','CHI1 ','PI1  ','P1   ','S1   ',  &
               'K1   ','PSI1 ','PHI1 ','THE1 ','J1   ','OO1  ','UPS1 ','OQ2  ','EPS2 ','2N2  ',  &
               'MU2  ','N2   ','NU2  ','GAM2 ','H1   ','M2   ','H2   ','LDA2 ','L2   ','T2   ',  &
               'S2   ','R2   ','K2   ','ETA2 ','M3   ','2PO1 ','SO1  ','ST36 ','2NS2 ','ST37 ',  &
               'ST1  ','ST2  ','ST3  ','O2   ','ST4  ','SNK2 ','OP2  ','MKS2 ','ST5  ','ST6  ',  &
               '2SK2 ','MSN2 ','ST7  ','2SM2 ','ST38 ','SKM2 ','2SN2 ','NO3  ','MO3  ','NK3  ',  &
               'SO3  ','MK3  ','SP3  ','SK3  ','ST8  ','N4   ','3MS4 ','ST39 ','MN4  ','ST40 ',  &
               'ST9  ','M4   ','ST10 ','SN4  ','KN4  ','MS4  ','MK4  ','SL4  ','S4   ','SK4  ',  &
               'MNO5 ','2MO5 ','3MP5 ','MNK5 ','2MP5 ','2MK5 ','MSK5 ','3KM5 ','2SK5 ','ST11 ',  &
               '2NM6 ','ST12 ','ST41 ','2MN6 ','ST13 ','M6   ','MSN6 ','MKN6 ','2MS6 ','2MK6 ',  &
               'NSK6 ','2SM6 ','MSK6 ','ST42 ','S6   ','ST14 ','ST15 ','M7   ','ST16 ','3MK7 ',  &
               'ST17 ','ST18 ','3MN8 ','ST19 ','M8   ','ST20 ','ST21 ','3MS8 ','3MK8 ','ST22 ',  &
               'ST23 ','ST24 ','ST25 ','ST26 ','4MK9 ','ST27 ','ST28 ','M10  ','ST29 ','ST30 ',  &
               'ST31 ','ST32 ','ST33 ','M12  ','ST34 ','ST35 '/)

    character(LEN = 5), dimension(NConsH2), parameter  :: KONCO = (/                 &
                'P1   ','O1   ','S2   ','O1   ','M2   ','N2   ','S2   ','N2   ','S2   ','M2   ', & 
                'S2   ','N2   ','K2   ','S2   ','M2   ','N2   ','K2   ','S2   ','M2   ','S2   ', &
                'K2   ','O1   ','K2   ','N2   ','S2   ','S2   ','N2   ','K2   ','O1   ','P1   ', &
                'M2   ','K2   ','S2   ','M2   ','K2   ','S2   ','S2   ','N2   ','M2   ','K2   ', &
                'S2   ','K2   ','M2   ','S2   ','N2   ','K2   ','M2   ','S2   ','N2   ','S2   ', &
                'M2   ','M2   ','S2   ','N2   ','S2   ','K2   ','M2   ','S2   ','N2   ','N2   ', &
                'O1   ','M2   ','O1   ','N2   ','K1   ','S2   ','O1   ','M2   ','K1   ','S2   ', &
                'P1   ','S2   ','K1   ','M2   ','N2   ','S2   ','N2   ','M2   ','S2   ','M2   ', & 
                'S2   ','N2   ','K2   ','M2   ','N2   ','M2   ','S2   ','K2   ','M2   ','N2   ', &
                'K2   ','S2   ','M2   ','M2   ','K2   ','S2   ','S2   ','N2   ','K2   ','N2   ', &
                'M2   ','S2   ','M2   ','K2   ','S2   ','L2   ','S2   ','S2   ','K2   ','M2   ', &
                'N2   ','O1   ','M2   ','O1   ','M2   ','P1   ','M2   ','N2   ','K1   ','M2   ', &
                'P1   ','M2   ','K1   ','M2   ','S2   ','K1   ','K2   ','K1   ','M2   ','S2   ', &
                'K1   ','N2   ','K2   ','S2   ','N2   ','M2   ','N2   ','M2   ','K2   ','S2   ', &
                'M2   ','S2   ','K2   ','M2   ','N2   ','M2   ','N2   ','K2   ','S2   ','M2   ', &
                'M2   ','S2   ','N2   ','M2   ','K2   ','N2   ','M2   ','S2   ','M2   ','K2   ', &
                'N2   ','S2   ','K2   ','S2   ','M2   ','M2   ','S2   ','K2   ','M2   ','S2   ', & 
                'K2   ','S2   ','M2   ','N2   ','O1   ','N2   ','M2   ','K1   ','M2   ','M2   ', &
                'S2   ','O1   ','M2   ','K1   ','M2   ','S2   ','K2   ','O1   ','M2   ','N2   ', &
                'M2   ','N2   ','M2   ','N2   ','K2   ','S2   ','M2   ','M2   ','S2   ','N2   ', &
                'M2   ','N2   ','K2   ','M2   ','S2   ','M2   ','K2   ','M2   ','S2   ','N2   ', &
                'K2   ','M2   ','S2   ','M2   ','S2   ','K2   ','M2   ','N2   ','K1   ','M2   ', &
                'N2   ','K1   ','M2   ','K1   ','M2   ','S2   ','K1   ','M2   ','N2   ','M2   ', &
                'M2   ','N2   ','S2   ','M2   ','S2   ','M2   ','N2   ','S2   ','K2   ','M2   ', &
                'S2   ','M2   ','S2   ','K1   ','M2   ','M2   ','S2   ','M2   ','N2   ','K2   ', &
                'S2   '/)
                                    
    integer, dimension(NTidal), parameter :: II = (/                         &
                  0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,  &
                  1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  &
                  1,   1,   1,   2,   2,   2,   2,   2,   2,   2,   2,   2,  &
                  2,   2,   2,   2,   2,   2,   2,   2,   3/)
    
    integer, dimension(NTidal), parameter :: JJ = (/                         &
                  0,   0,   0,   1,   1,   2,   2,  -4,  -3,  -3,  -2,  -2,  &
                 -1,  -1,   0,   0,   0,   1,   1,   1,   1,   1,   1,   2,  &
                  2,   3,   4,  -3,  -3,  -2,  -2,  -1,  -1,   0,   0,   0,  &
                  0,   1,   1,   2,   2,   2,   2,   3,   0/)

    integer, dimension(NTidal), parameter :: KK = (/                         &
                  0,   1,   2,  -2,   0,  -2,   0,   2,   0,   2,   0,   2,  &
                  0,   2,  -2,   0,   2,  -3,  -2,  -1,   0,   1,   2,  -2,  &
                  0,   0,   0,   0,   2,   0,   2,   0,   2,  -2,  -1,   0,  &
                  1,  -2,   0,  -3,  -2,  -1,   0,   0,   0/)

    integer, dimension(NTidal), parameter :: LL = (/                         &
                  0,   0,   0,   1,  -1,   0,   0,   1,   2,   0,   1,  -1,  &
                  0,   0,   1,   1,  -1,   0,   0,   0,   0,   0,   0,   1,  &
                 -1,   0,  -1,   3,   1,   2,   0,   1,  -1,   2,   0,   0,  &
                  0,   1,  -1,   0,   0,   0,   0,  -1,   0/)
                  
    integer, dimension(NTidal), parameter :: MM = (/                         &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0/)
                  
    integer, dimension(NTidal), parameter :: NN = (/                         &
                  0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   1,   0,   1,   0,  -1,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,  &
                 -1,   0,   0,   1,   0,  -1,   0,   0,   0/)

    integer, dimension(NCompt), parameter :: NJ = (/                         &                 
                  1,   1,   1,   1,   1,   1,   1,   2,   5,   4,  10,   5,  &
                  8,   5,   1,   9,   2,   1,   6,   2,  10,   1,   5,   4,  &
                 10,   8,   5,   2,   3,   4,   3,   4,   4,   3,   2,   9,  &
                  1,   1,   5,   1,   3,   2,   5,   7,   1,   2,   2,   3,  &
                  2,   2,   3,   4,   3,   1,   3,   3,   2,   3,   3,   4,  &
                  2,   3,   4,   2,   3,   3,   2,   2,   2,   2,   2,   2,  &
                  2,   2,   3,   1,   2,   4,   2,   3,   4,   1,   3,   2,  &
                  2,   2,   2,   2,   1,   2,   3,   2,   2,   3,   2,   2,  &
                  3,   3,   2,   3,   2,   4,   3,   2,   4,   1,   3,   3,  &
                  2,   2,   3,   2,   3,   3,   1,   3,   3,   1,   3,   2,  &
                  4,   2,   2,   4,   1,   3,   3,   2,   2,   4,   2,   3,  &
                  3,   3,   2,   3,   2,   1,   3,   2,   4,   2,   3,   1,  &
                  2,   4/)

    real, dimension(NTidal ), parameter :: Semi = (/                         &
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,-0.25,-0.25,-0.25,  &
               -0.25,-0.25,-0.25,-0.75,-0.75,-0.75,-0.75,-0.25,-0.25,-0.75,  &
               -0.75,-0.75,-0.75,-0.75,-0.75,-0.75,-0.75, 0.00, 0.00, 0.00,  &
                0.00, 0.00, 0.00,-0.50,-0.50, 0.00, 0.00,-0.50,-0.50, 0.00,  &
                0.00,-0.50, 0.00, 0.00,-0.50/)

    integer, dimension(NConsH1), parameter :: LDEL = (/                      &
                  0,   0,   0,   0,   0,   0,   0,  -1,   0,  -2,  -1,  -1,  &
                  0,   0,  -1,   0,   0,   2,  -2,  -2,  -1,  -1,  -1,   0,  &
                 -1,   0,   1,   2,   0,   0,   1,   2,   2,  -1,   0,   0,  &
                  1,   1,   1,   2,   2,  -2,  -1,   0,   0,   0,   0,  -2,  &
                 -2,  -2,  -1,  -1,  -1,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   1,   2,   2,   0,   0,  -2,  -1,  -1,  -1,   0,  &
                  0,   0,   0,   1,   1,   0,  -2,  -2,   0,   0,   0,  -2,  &
                 -1,   0,   0,   0,   0,   0,   1,   1,   1,   1,   2,   2,  &
                  2,  -2,  -2,  -2,  -1,  -1,   0,   0,   0,  -2,   0,   0,  &
                  1,   1,  -1,   0,  -1,  -1,   0,  -2,  -1,  -1,   0,  -1,  &
                 -1,   0,  -2,  -1,   0,   0,   0,   1,   2,   2,  -2,  -1,  &
                  0,   0,   1,  -1,  -1,   0,   0,   1,   1,   1,   2,   2,  &
                  0,   0,   0,   2,   2,   2,   2,   0,   0,   1,   2,   0,  &
                  0,  -1,  -1,   0,   0,   0,   0,   0,   0,   1,   1,   1,  &
                  2,   0/)
                  
    integer, dimension(NConsH1), parameter :: MDEL = (/                      &
                  0,   0,   0,   0,   0,   0,   0,   0,  -1,  -2,  -1,   0,  &
                 -2,  -1,   0,  -2,  -1,   0,  -3,  -2,  -2,  -1,   0,  -2,  &
                  0,  -1,   0,   0,  -2,  -1,   0,   0,   1,   0,  -2,  -1,  &
                 -1,   0,   1,   0,   1,   0,   0,  -1,   1,   2,  -1,  -2,  &
                 -1,   0,  -1,   0,   1,  -1,   1,   2,  -1,   1,  -1,  -2,  &
                 -1,   0,   0,   0,   1,   0,   1,  -1,  -1,   0,   1,  -2,  &
                 -1,   1,   2,   0,   1,   1,   0,   1,   0,   1,   2,  -1,  &
                  0,  -1,   1,  -1,   1,   2,  -1,   0,   1,   2,   0,   1,  &
                  2,  -1,   0,   1,   0,   1,   1,   2,   3,   0,   1,   2,  &
                  0,   1,   0,  -1,  -1,   0,  -1,  -2,  -1,   0,  -1,  -1,  &
                  0,  -1,  -2,   0,  -2,  -1,  -1,   0,   0,   1,  -2,   0,  &
                 -1,  -1,   0,  -1,   0,  -2,  -1,  -1,   0,   1,   0,   1,  &
                 -1,  -1,  -1,  -1,   0,   1,   2,   0,  -1,   0,   0,   0,  &
                  1,   0,   1,  -1,   1,   2,  -1,   1,   2,   0,   1,   2,  &
                  0,  -1/)
                  
    integer, dimension(NConsH1), parameter :: NDEL = (/                      &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   2,   0,   0,   0,  -2,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,  -2,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,  &
                  2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0/)
                  
    integer, dimension(NConsH1), parameter :: IR = (/                        &
                  0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   1,   1,  &
                  0,   0,   1,   0,   0,   0,   0,   0,   1,   1,   1,   0,  &
                  0,   0,   1,   0,   0,   0,   1,   0,   0,   1,   0,   0,  &
                  1,   1,   1,   0,   0,   0,   1,   0,   0,   0,   0,   0,  &
                  0,   0,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0,  &
                  0,   0,   1,   0,   0,   0,   0,   0,   1,   1,   1,   0,  &
                  0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,   0,  &
                  1,   0,   0,   0,   0,   0,   1,   1,   1,   1,   0,   0,  &
                  0,   0,   0,   0,   1,   1,   0,   0,   0,   0,   0,   0,  &
                  1,   1,   2,   0,   2,   2,   0,   0,   2,   2,   0,   2,  &
                  2,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   2,  &
                  0,   0,   0,   2,   2,   0,   0,   2,   2,   2,   0,   0,  &
                  0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,  &
                  0,   2,   2,   0,   0,   0,   0,   0,   0,   2,   2,   2,  &
                  0,   0/)
                  
                  
    real, dimension(NConsH1), parameter :: PH = (/                           &
                0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.75, 0.00, 0.50,  &
                0.75, 0.75, 0.50, 0.00, 0.75, 0.50, 0.00, 0.50, 0.50, 0.50,  &
                0.75, 0.75, 0.75, 0.50, 0.00, 0.00, 0.75, 0.50, 0.50, 0.00,  &
                0.75, 0.50, 0.00, 0.25, 0.50, 0.00, 0.25, 0.75, 0.25, 0.50,  &
                0.50, 0.00, 0.25, 0.50, 0.50, 0.50, 0.00, 0.50, 0.00, 0.00,  &
                0.75, 0.25, 0.75, 0.50, 0.00, 0.50, 0.50, 0.00, 0.50, 0.00,  &
                0.50, 0.50, 0.75, 0.50, 0.50, 0.00, 0.50, 0.00, 0.75, 0.25,  &
                0.75, 0.00, 0.50, 0.00, 0.50, 0.25, 0.25, 0.00, 0.00, 0.00,  &
                0.00, 0.50, 0.50, 0.00, 0.25, 0.50, 0.00, 0.50, 0.00, 0.50,  &
                0.75, 0.25, 0.25, 0.25, 0.50, 0.50, 0.50, 0.50, 0.00, 0.00,  &
                0.25, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.25,  &
                0.25, 0.50, 0.25, 0.25, 0.50, 0.50, 0.25, 0.25, 0.50, 0.25,  &
                0.25, 0.50, 0.50, 0.00, 0.00, 0.50, 0.50, 0.75, 0.00, 0.50,  &
                0.00, 0.25, 0.50, 0.50, 0.50, 0.75, 0.75, 0.00, 0.50, 0.25,  &
                0.75, 0.75, 0.00, 0.00, 0.50, 0.50, 0.50, 0.00, 0.50, 0.50,  &
                0.50, 0.00, 0.00, 0.75, 0.00, 0.50, 0.00, 0.75, 0.75, 0.50,  &
                0.00, 0.00, 0.50, 0.00, 0.00, 0.75, 0.75, 0.75, 0.50, 0.50/)

    real, dimension(NConsH1), parameter :: EE = (/                              &
                 0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  &
                 0.0360,  0.1906,  0.0063,  0.0241,  0.0607,  0.0063,  0.1885,  &
                 0.0095,  0.0061,  0.1884,  0.0087,  0.0007,  0.0039,  0.0010,  &
                 0.0115,  0.0292,  0.0057,  0.0008,  0.1884,  0.0018,  0.0028,  &
                 0.0058,  0.1882,  0.0131,  0.0576,  0.0175,  0.0003,  0.0058,  &
                 0.1885,  0.0004,  0.0029,  0.0004,  0.0064,  0.0010,  0.0446,  &
                 0.0426,  0.0284,  0.2170,  0.0142,  0.2266,  0.0057,  0.0665,  &
                 0.3596,  0.0331,  0.2227,  0.0290,  0.0290,  0.2004,  0.0054,  &
                 0.0282,  0.2187,  0.0078,  0.0008,  0.0112,  0.0004,  0.0004,  &
                 0.0015,  0.0003,  0.3534,  0.0264,  0.0002,  0.0001,  0.0007,  &
                 0.0001,  0.0001,  0.0198,  0.1356,  0.0029,  0.0002,  0.0001,  &
                 0.0190,  0.0344,  0.0106,  0.0132,  0.0384,  0.0185,  0.0300,  &
                 0.0141,  0.0317,  0.1993,  0.0294,  0.1980,  0.0047,  0.0027,  &
                 0.0816,  0.0331,  0.0027,  0.0152,  0.0098,  0.0057,  0.0037,  &
                 0.1496,  0.0296,  0.0240,  0.0099,  0.6398,  0.1342,  0.0086,  &
                 0.0611,  0.6399,  0.1318,  0.0289,  0.0257,  0.1042,  0.0386,  &
                 0.0075,  0.0402,  0.0373,  0.0061,  0.0117,  0.0678,  0.0374,  &
                 0.0018,  0.0104,  0.0375,  0.0039,  0.0008,  0.0005,  0.0373,  &
                 0.0373,  0.0042,  0.0042,  0.0036,  0.1429,  0.0293,  0.0330,  &
                 0.0224,  0.0447,  0.0001,  0.0004,  0.0005,  0.0373,  0.0001,  &
                 0.0009,  0.0002,  0.0006,  0.0002,  0.0217,  0.0448,  0.0366,  &
                 0.0047,  0.2505,  0.1102,  0.0156,  0.0000,  0.0022,  0.0001,  &
                 0.0001,  0.2535,  0.0141,  0.0024,  0.0004,  0.0128,  0.2980,  &
                 0.0324,  0.0187,  0.4355,  0.0467,  0.0747,  0.0482,  0.0093,  &
                 0.0078,  0.0564/)
             
    real, dimension(NConsH2), parameter :: Coef = (/                         &
                2.00,-1.00, 1.00,-1.00, 2.00, 1.00,-2.00, 2.00,-1.00, 3.00,  &
               -2.00, 2.00, 1.00,-2.00, 1.00, 1.00, 1.00,-2.00, 2.00, 1.00,  &
               -2.00, 2.00, 2.00, 1.00,-2.00, 1.00, 1.00,-1.00, 1.00, 1.00,  &
                1.00, 1.00,-1.00, 1.00, 2.00,-2.00, 2.00, 1.00,-1.00,-1.00,  &
                2.00,-1.00, 1.00, 1.00,-1.00, 2.00, 1.00,-1.00,-1.00, 2.00,  &
               -1.00, 2.00, 1.00,-2.00, 1.00, 1.00,-1.00, 2.00,-1.00, 1.00,  &
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,  &
                1.00, 1.00, 1.00, 2.00, 1.00,-1.00, 2.00, 3.00,-1.00, 1.00,  &
                1.00, 1.00,-1.00, 1.00, 1.00, 2.00, 1.00,-1.00, 1.00, 1.00,  &
                1.00,-1.00, 2.00, 2.00, 1.00,-1.00, 1.00, 1.00, 1.00, 1.00,  &
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 1.00, 1.00, 1.00,  &
                1.00, 1.00, 2.00, 1.00, 3.00,-1.00, 1.00, 1.00, 1.00, 2.00,  &
                1.00, 2.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 2.00,  &
                1.00, 3.00, 1.00,-1.00, 2.00, 1.00, 2.00, 1.00, 1.00,-1.00,  &
                3.00, 1.00,-1.00, 2.00, 1.00, 2.00, 1.00, 1.00,-1.00, 3.00,  &
                1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 1.00, 2.00, 1.00,  &
                1.00, 1.00, 1.00, 2.00, 1.00, 1.00, 1.00, 1.00, 2.00, 2.00,  &
               -1.00, 3.00, 2.00, 1.00, 1.00, 2.00, 1.00, 1.00, 3.50, 2.00,  &
                1.00, 1.00, 3.00, 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 2.00,  &
                3.00, 1.00, 3.00, 1.00, 1.00,-1.00, 4.00, 2.00, 1.00, 1.00,  &
                2.00, 1.00, 1.00, 3.00, 1.00, 3.00, 1.00, 1.00, 1.00, 1.00,  &
                1.00, 2.00, 2.00, 2.00, 1.00, 1.00, 2.00, 2.00, 1.00, 3.00,  &
                1.00, 1.00, 4.00, 1.00, 3.00, 1.00, 1.00, 4.00, 1.00, 5.00,  &
                3.00, 1.00, 1.00, 4.00, 1.00, 2.00, 1.00, 1.00, 1.00, 3.00,  &
                2.00, 4.00, 1.00, 1.00, 6.00, 5.00, 1.00, 3.00, 1.00, 1.00,  &
                1.00/)

    !--------------------------------------------------------------------------


    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructToga(TogaID, NWaves, WaveName, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TogaID
        integer,           intent(IN )              :: NWaves
        character(LEN = 5), pointer, dimension(:)   :: WaveName
        integer, optional, intent(OUT)              :: STAT     

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        logical                                     :: validwave = .FALSE.
        integer                                     :: KT, I
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mToga_)) then
            nullify (FirstToga)
            call RegisterModule (mToga_) 
        endif

        call Ready(TogaID, ready_)    

        if (ready_ .EQ. OFF_ERR_) then
  
            call AllocateInstance

            !Allocates memory for Toga data
            allocate(Me%indx (NCompt))
            allocate(Me%Sig  (NCompt))
            allocate(Me%kon  (NCompt))
            allocate(Me%twoc (NCompt))
            allocate(Me%CH   (NCompt))     
            allocate(Me%CHP  (NCompt))     
            allocate(Me%VMare(NCompt))  
            allocate(Me%UMare(NCompt))  
            allocate(Me%FMare(NCompt))  
            allocate(Me%Freq (NCompt)) 

            Me%Freq = 0.0

do2 :       do I = 1, NWaves
                validwave = .FALSE.
do1 :           do KT=1,NCompt
                    if (KonTB(KT) == WaveName(I)) then 
                        validwave = .TRUE.
                        exit do1
                    endif
                end do do1

if1 :           if (validwave) then
                    Me%indx(I)  = KT
                    Me%Sig (I)  = Me%Freq(KT)
                    Me%kon (I)  = WaveName(I) 

                    Me%twoc(I) = 2.D0*COS(TwoPI*Me%Sig(I)*DTHora)
                else
                    write (*,90) WaveName(I)
90                  format(//////T25,'W A R N I N G !!!',//T25, 'Routine ASTRO cannot find constituent ',A5 /////)
                    stop 'Subroutine ConstructToga; Module ModuleToga. ERR02'
                end if if1
            end do do2


            !Returns ID
            TogaID = Me%InstanceID
            STAT_  = SUCCESS_

        else 
            stop 'ModuleToga - ConstructToga - ERR99' 
        end if 


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructToga

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
    
        !Local-----------------------------------------------------------------
        type (T_Toga), pointer                      :: NewObjToga
        type (T_Toga), pointer                      :: PreviousObjToga


        !Allocates new instance
        allocate (NewObjToga)
        nullify  (NewObjToga%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstToga)) then
            FirstToga            => NewObjToga
            Me                   => NewObjToga
        else
            PreviousObjToga      => FirstToga
            Me                   => FirstToga%Next
            do while (associated(Me))
                PreviousObjToga  => Me
                Me               => Me%Next
            enddo
            Me                   => NewObjToga
            PreviousObjToga%Next => NewObjToga
        endif

        Me%InstanceID = RegisterNewInstance (mTOGA_)

    end subroutine AllocateInstance

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine TogaLevel(TogaID, WaterLevel, Decimal_Latitude, TimeReference,            &
                         ReferenceLevel, NWaves, WaveAmplitude, WavePhase, time_,        &
                         STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TogaID
        real,              intent(OUT)              :: WaterLevel
        real,              intent(IN )              :: Decimal_Latitude
        real,              intent(IN )              :: TimeReference
        real,              intent(IN )              :: ReferenceLevel
        integer,           intent(IN )              :: NWaves
        real, pointer, dimension(:)                 :: WaveAmplitude
        real, pointer, dimension(:)                 :: WavePhase
        type(T_Time)                                :: time_
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: K, INDK
        integer                                     :: GregDay
        integer                                     :: ready_
        real                                        :: TempHora
        real                                        :: ANG, DTPY, HGT
        integer                                     :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TogaID, ready_)

if1 :   if (ready_ .NE. OFF_ERR_) then

            TempHora = TimeHours(time_)

            Me%DS   = null_real
            Me%DH   = null_real
            Me%DP   = null_real
            Me%DNP  = null_real
            Me%DPP  = null_real
            Me%D1   = null_real
            Me%HH   = null_real
            Me%Tau  = null_real
            Me%Hour = null_real
            Me%H    = null_real
            Me%S    = null_real
            Me%ENP  = null_real
            Me%P    = null_real
            Me%PP   = null_real

            ! calculates the Gregorian day and the number of hours correspnding to it

            call DateToGregorianDay(time_, GregDay)

            Me%Hour = GregDay * 24d0

            ! calculates the astronomic factors

            call ASTRO (Decimal_Latitude)

            Do K =1, NWaves
                Me%Sig(K)=Me%Freq(Me%indx(K))
            EndDo

            ! Calculates the level

            DO K = 1, NWaves
                INDK   = Me%indx(K)
                DTPY   = Me%FMare(INDK)*WaveAmplitude(K)
                ANG    = Me%VMare(INDK) + (TempHora - TimeReference)*Me%Sig(K)  &
                                  + Me%UMare(INDK)-WavePhase(K)
                ANG    = ANG - INT(ANG)
                Me%CHP(K) = DTPY*COS((ANG-DTHora*Me%Sig(K))*TwoPI)
                Me%CH(K)  = DTPY*COS(ANG*TwoPI)
            ENDDO

            HGT = 0.0
            DO K=1,NWaves
                DTPY   = Me%CH(K)
                Me%CH(K)  = Me%twoc(K)*DTPY-Me%CHP(K)
                Me%CHP(K) = DTPY
                HGT    = HGT+Me%CH(K)
            ENDDO

            WaterLevel = HGT + ReferenceLevel

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine TogaLevel

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    !                        SUBROUTINE ASTRO                               
    !--------------------------------------------------------------------------
    !     THIS ROUTINE CALCULATES THE FREQUENCY Sig IN CPH, THE  ASTRO  ARG 
    !     V IN CYCLES , THE  NODAL  CORRECTION  PHASE U AND  AMP F  FOR THE 
    !     CONSTITUENTS KONTAB(1 ... NCOMPT).                                
    !     THE INPUT parameter HOUR GIVES THE TIME AT WHICH THE  RESULTS ARE 
    !     TO BE CALCULATED. HOUR HAS ITS ORIGIN AROUND 1900 AND IS OBTAINED 
    !     USING ROUTINE CDAY.                                               
    !                                                                       
    ! DATA INPUT                                                            
    ! ==========                                                            
    !       KON = CONSTITUENT NAME                                          
    !       II,JJ,KK,LL,MM,NN = THE SIX DOODSON NUMBERS                     
    !       Semi= PHASE CORRECTION                                          
    !       NJ  = THE NUMBER OF SATELLITES FOR THIS CONSTITUENT.            
    !       LDEL,MDEL,NDEL= THE CHANGES IN THE LAST THREE DOODSON NUMBERS   
    !                       FROM THOSE OF THE MAIN CONSTITUENT.             
    !       PH  = THE PHASE CORRECTION                                      
    !       EE  = THE AMPLITUDE RATIO OF THE SATELLITE TIDAL POTENTIAL TO   
    !             THAT OF THE MAIN CONSTITUENT.                             
    !       IR  = 1 IF THE AMPLITUDE  RATIO  HAS  TO BE MULTIPLIED BY THE   
    !               LATITUDE CORRECTION  FACTOR  FOR DIURNAL CONSTITUENTS   
    !             2 IF THE AMPLITUDE RATIO  HAS  TO  BE MULTIPLIED BY THE   
    !               LATITUDE CORRECTION  FACTOR FOR  Semi-DIURNAL CONSTI-   
    !               TUENTS. OTHERWISE IF NO CORRECTION IS REQUIRED TO THE   
    !               AMPLITUDE RATIO.                                        
    !      KON  = NAME  OF THE SHALLOW WATER CONSTITUENT                    
    !      NJ   = NUMBER OF MAIN CONSTITUENTS FROM WHICH IT IS DERIVED.     
    !      Coef,KONCO = COMBINATION NUMBER AND NAME OF THESE MAIN           
    !                   CONSTITUENTS.                                      

    subroutine ASTRO(Decimal_Latitude)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: Decimal_Latitude

        !Local-----------------------------------------------------------------
        type (T_Time)                               :: Time0
        integer                                     :: KD0, K, K1, IUU, L, JL, J1, J
        integer                                     :: IV, JBASE, INTDYS
        real(8)                                     :: DTau, RR, UUDBL, UU
        real                                        :: SLAT, VDBL, SUMC, SUMS

        !----------------------------------------------------------------------

        SLAT=SIN(PI * Decimal_Latitude / 180.0)

        Me%D1 = Me%Hour / 24.d0


        call SetDate(Time0, 1899, 12, 31, 0, 0, 0)

        call DateToGregorianDay(Time0,KD0)

        Me%D1 = Me%D1 - float(KD0) - 0.5d0

        Call TideFreq

        INTDYS      = idint(Me%Hour / 24.d0)
        Me%HH  = Me%Hour - float(INTDYS*24)
        Me%Tau = Me%HH / 24.D0 + Me%H - Me%S
        DTau   = 365.d0 + Me%DH - Me%DS
        JBASE       = 0



do4 :   DO K = 1, NTidal
            Me%Freq(K)=(II(K)*DTau+JJ(K)* Me%DS +KK(K)* Me%DH                      &
                          + LL(K)*Me%DP+MM(K) * Me%DNP +NN(K)* Me%DPP) / (365.0 * 24.0)

            VDBL = II(K)* Me%Tau +JJ(K)* Me%S +KK(K)* Me%H                         &
                          + LL(K)* Me%P +MM(K)* Me%ENP +NN(K)* Me%PP +Semi(K)
            IV   = VDBL
            IV   = (IV/2)*2
            Me%VMare(K) = VDBL-IV
            J1   = JBASE+1
            JL   = JBASE+NJ(K)
            SUMC = 1.
            SUMS = 0.

do5 :       DO J = J1, JL
                ! Here the satellite amplitude ratio adjustment for Latitude is made
                RR = EE(J)
                L  = IR(J)+1
if4 :           if     (L .EQ. 0) then
                    RR=EE(J)*0.36309*(1.-5.*SLAT*SLAT)/SLAT
                else if (L .GT. 0) then
                    RR=EE(J)*2.59808*SLAT
                end if if4

                UUDBL = LDEL(J)* Me%P +MDEL(J)* Me%ENP +NDEL(J)* Me%PP +PH(J)
                IUU   = UUDBL
                UU    = UUDBL-IUU
                SUMC  = SUMC+RR*COS(UU*TwoPI)
                SUMS  = SUMS+RR*SIN(UU*TwoPI)
            end do do5

            Me%FMare(K) = SQRT(SUMC*SUMC+SUMS*SUMS)
            Me%UMare(K) = ATAN2(SUMS,SUMC)/TwoPI
            JBASE    = JL
        end do do4


        JBASE = 0
        K1    = NTidal + 1
if2 :   if (.NOT. (K1 .GT. NCompt)) then
do1 :       DO K = K1, Ncompt
                Me%FMare(K) = 1.0
                Me%VMare(K) = 0.0
                Me%UMare(K) = 0.0
                J1               = JBASE + 1
                JL               = JBASE + NJ(K)
                Me%Freq(K)  = 0.0

do3 :           do J = J1, JL
do2 :               do L = 1,  NTidal
if1 :               if (KonTB(L) .EQ. KONCO(J)) Then
                        Me%FMare(K) = Me%FMare(K)*Me%FMare(L)**ABS(Coef(J))
                        Me%VMare(K) = Me%VMare(K)+Coef(J)*Me%VMare(L)
                        Me%UMare(K) = Me%UMare(K)+Coef(J)*Me%UMare(L)
                        Me%Freq (K) = Me%Freq (K)+Coef(J)*Me%Freq(L)    
                                
                        cycle do3
                    end if if1
                    end do do2

                    Write (*,241) KONCO(J)
                    STOP 'Subroutine ASTRO; Module ModuleToga. ERR01' 
                end do do3
                JBASE = JL
            end do do1
        end if if2

        !Format----------------------------------------

241     format(///T10,'A T T E N T I O N !!!',//T10, &
               'ROUTINE ASTRO UNABLE TO FIND CONSTITUENT ' ,A5 //////)

        !----------------------------------------------------------------------

    end subroutine ASTRO

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    !                                     SUBROUTINE TideFreq
    !   this subroutine calculates the following five ephermides of the sun and moon:
    !   
    !   h  = mean longitude of the sum
    !   pp = mean longitude of the solar perigee
    !   s  = mean longitude of the moon
    !   p  = mean longitude of the lunar perigee
    !   ENP = negative of the longitude of the mean ascending node and their rates of change.
    !   
    !   Units for the  ephermides  are  cycles and for their derivatives are cycles/365 days
    !   
    !   The formulae for calculating this ephermides were taken from pages 98 and 107 of the
    !   Explanatory Supplement  to  the  Astronomical Ephermeris and the American Ephermeris
    !   and Nautical Almanac (1961)

    subroutine TideFreq

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8)                                     :: D2, F, F2

        !----------------------------------------------------------------------

        D2 = Me%D1 * 1.d-4
        F  = 360.d0
        F2 = F/365.d0

        Me%H   = 279.696678d0+.9856473354d0* Me%D1 +.00002267d0*d2*D2
        Me%S   = 270.434164d0+13.1763965268d0* Me%D1 -.000085d0*d2*D2+.000000039d0*D2**3
        Me%P   = 334.329556d0+.1114040803d0* Me%D1 -.0007739d0*d2*D2-.00000026d0*D2**3
        Me%ENP =-259.183275d0+.0529539222d0* Me%D1 -.0001557d0*d2*D2-.00000005d0*D2**3

        Me%H   = Me%H   / F
!        Me%PP  = Me%PP  / F  !RCM
        Me%S   = Me%S   / F
        Me%P   = Me%P   / F
        Me%ENP = Me%ENP / F

        Me%H   = Me%H   - dint(Me%H  )
!        Me%PP  = Me%PP  - dint(Me%PP )  !RCM
        Me%PP  = 0.0
        Me%S   = Me%S   - dint(Me%S  )
        Me%P   = Me%P   - dint(Me%P  )
        Me%ENP = Me%ENP - dint(Me%ENP)

        Me%DH  = .9856473354d0+2.d-8*.00002267d0* Me%D1 
        Me%DPP = .0000470684d0+2.d-8*.0000339d0* Me%D1 +3.d-12*.00000007d0* Me%D1**2
        Me%DS  = 13.1763965268d0-2.d-8*.000085d0* Me%D1 +3.d-12*.000000039d0* Me%D1**2
        Me%DP  = .1114040803d0-2.d-8*.0007739d0* Me%D1 -3.d-12*.00000026d0* Me%D1**2
        Me%DNP = .0529539222d0-2.d-8*.0001557d0* Me%D1 -3.d-12*.00000005d0* Me%D1**2

        Me%DH  = Me%DH  / F2
        Me%DPP = Me%DPP / F2
        Me%DS  = Me%DS  / F2
        Me%DP  = Me%DP  / F2
        Me%DNP = Me%DNP / F2

        !----------------------------------------------------------------------

    end subroutine TideFreq


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillToga(TogaID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: TogaID
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_              
        integer                                     :: STAT_
        integer                                     :: nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TogaID, ready_)

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTOGA_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%indx )
                deallocate(Me%Sig  )
                deallocate(Me%kon  )
                deallocate(Me%twoc )
                deallocate(Me%CH   )
                deallocate(Me%CHP  )
                deallocate(Me%VMare)
                deallocate(Me%UMare)
                deallocate(Me%FMare)
                deallocate(Me%Freq )

                call DeallocateInstance

                TogaID = 0
                STAT_  = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_
              
        !----------------------------------------------------------------------

    end subroutine KillToga

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance 

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Toga), pointer                      :: AuxObjToga
        type (T_Toga), pointer                      :: PreviousObjToga

        !Updates pointers
        if (Me%InstanceID == FirstToga%InstanceID) then
            FirstToga => FirstToga%Next
        else
            PreviousObjToga => FirstToga
            AuxObjToga      => FirstToga%Next
            do while (AuxObjToga%InstanceID /= Me%InstanceID)
                PreviousObjToga => AuxObjToga
                AuxObjToga      => AuxObjToga%Next
            enddo

            !Now update linked list
            PreviousObjToga%Next => AuxObjToga%Next

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

    subroutine Ready (ObjToga_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjToga_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjToga_ID > 0) then
            call LocateObjToga(ObjToga_ID)
            ready_ = VerifyReadLock (mTOGA_,  Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjToga (ObjTogaID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTogaID

        !Local-----------------------------------------------------------------

        Me => FirstToga
        do while (associated (Me))
            if (Me%InstanceID == ObjTogaID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleToga - LocateObjToga - ERR01'

    end subroutine LocateObjToga

    !--------------------------------------------------------------------------

end module ModuleToga

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

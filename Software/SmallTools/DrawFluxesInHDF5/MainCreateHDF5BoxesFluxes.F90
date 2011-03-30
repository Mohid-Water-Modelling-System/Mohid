!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Shell
! PROGRAM       : MainShell
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig /Luis Fernandes - v4.0
! DESCRIPTION   : Shell to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program CreateHDF5BoxesFluxes

    use ModuleDrawFluxesInHDF5

    implicit none

    call ConstructCreateHDF5BoxesFluxes
    call ModifyCreateHDF5BoxesFluxes
    call KillCreateHDF5BoxesFluxes

end program CreateHDF5BoxesFluxes

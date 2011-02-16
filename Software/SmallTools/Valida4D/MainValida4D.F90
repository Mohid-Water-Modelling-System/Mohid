!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Valida4D
! PROGRAM       : MainValida4D
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : June 2010
! REVISION      : Paulo Leitão - v4.0
! DESCRIPTION   : Valida4D to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MainValida4D

    use ModuleGlobalData
    use ModuleValida4D

    implicit none
    
    
    !Begin---------------------------------------------------------------------
    
    call ConstructValida4D()

    call ModifyValida4D   () 
    
    call KillValida4D     () 
    

   
end program MainValida4D

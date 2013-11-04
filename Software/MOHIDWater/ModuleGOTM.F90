!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License, 
!as published by the Free Software Foundation.
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

!BOP
!
! !MODULE: turbulence: 
!
! !INTERFACE:
   module ModuleGOTM 
!
! !DESCRIPTION: 
! In this module all turbulence calculations are integrated. Further 
! descriptions are given in the single subroutines. 
!
! !USES:
!   use mtridiagonal
   !$ use omp_lib
   IMPLICIT NONE
   
!
!  Default all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
  
    public :: init_turbulence !Constructor
    public :: init_turbulence_parameters
    public :: clear_turbulence_parameters !griflet: new method
    public :: copy_turbulence_parameters !griflet: new method
    !public :: AssociateExportGOTM 
    !public :: Associate_to_ExportGOTM
    public :: do_turbulence    !Modifier
    ! Others are private
        private :: stabilityfunctions,do_tke,lengthscale
        private :: kolpran
        private :: internal_wave
        private :: algebraiclength, dissipationeq,ispralength,lengthscaleeq,potentialml
        private :: tkeeq,tkealgebraic
        private :: cmue_kc,cmue_ca,cmue_cb
    !    private :: cmue_bb,cmue_bbqe
        private :: cmue_kcqe,cmue_caqe,cmue_cbqe
        private :: cmue_ma,cmue_sg,cmue_rf

!   public stabilityfunctions,do_tke,lengthscale
!   public kolpran,do_turbulence,internal_wave 
    
    public :: kill_turbulence !Destructor ... griflet: The destructor was missing!

! !PUBLIC DATA TYPES
    public  :: T_Gotm
    public  :: T_GotmParameters
    public  :: T_GOTMexport
    public  :: T_GOTMimport !griflet
    public  :: T_GOTMTrid

! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  $Log$
!\
! !BUGS
!
!EOP
    
    !STAT
    integer, parameter :: UNKNOWN_        =-99
    integer, parameter :: SUCCESS_        = 0
    integer, parameter :: OFF_ERR_        = 1 
    integer, parameter :: ON_ERR_         = 2
    integer, parameter :: READ_LOCK_ERR_  = 3
    integer, parameter :: WRITE_LOCK_ERR_ = 4
    integer, parameter :: IDLE_ERR_       = 5
    integer, parameter :: NBUSERS_ERR_    = 6
    integer, parameter :: Null_int        = -9999
    real, parameter    :: Null_real       = -9.999E+15

    type T_GOTMexport
        double precision,  dimension(:), pointer    :: tke,L,eps 
        double precision,  dimension(:), pointer    :: num,nuh
    end type T_GOTMexport   

    type T_GOTMParameters
        !The data are written aligned
        double precision    :: cde,sig_e0,sig_e1
        double precision    :: craig_m
        
        !  These variables are read from the 'turbulence' namelist
        double precision    :: const_num= 5.e-4 
        double precision    :: const_nuh= 5.e-4 
        double precision    :: k_min=3e-6 
        double precision    :: eps_min=5e-10
        double precision    :: L_min=0.01 
        !  Turbulence parameters - 'turb_parameters' namelist.
        double precision    :: kappa=0.4
        double precision    :: Prandtl0=0.714 
        double precision    :: cm0=0.527
        double precision    :: cm_craig=0.73
        double precision    :: cw=100.
        double precision    :: galp=0.53
        !  The k-eps model - 'keps' namelist.
        double precision    :: ce1=1.44
        double precision    :: ce2=1.92
        double precision    :: ce3plus=1.
        double precision    :: ce3minus=-0.629
        double precision    :: sig_k=1.
        !
        !  The MY model - 'my' namelist.
        double precision    :: sl=0.2
        double precision    :: e1=1.8
        double precision    :: e2=1.33
        double precision    :: e3=1.8

        !  Stability functions - 'stabfunc' namelist.
        double precision    :: a1=0.92
        double precision    :: a2=0.74
        double precision    :: b1=16.6
        double precision    :: b2=10.1
        double precision    :: c1=0.08
        double precision    :: c2=0.7
        double precision    :: c3=0.2
        double precision    :: qeghmax=0.0233
        double precision    :: qeghmin=-0.28
        double precision    :: qeghcrit=0.02

        !  Internal wave model - the 'iw' namelist.
        double precision        :: alpha    = 0.0
        double precision        :: klimiw   = 1e-6
        double precision        :: rich_cr  = 0.7
        double precision        :: numiw    = 1e-4
        double precision        :: nuhiw    = 5e-5
        double precision        :: numshear = 5e-3
        integer                 :: iw_model = 0

        !MY namelist
        integer     :: MY_length=1

        !Namelist turbulence
        integer     :: turb_method=2 
        integer     :: tke_method=2
        integer     :: len_scale_method=8
        integer     :: stab_method=3
        logical     :: length_lim=.true.
        logical     :: craig_banner=.false.

        !Namelist stabilityfunctions
        logical     :: qesmooth=.true.

        !keps namelist
        logical     :: flux_bdy

    end type T_GOTMParameters

    type T_GOTMTrid
        double precision,  dimension(:), pointer    :: au,bu,cu,du,ru,qu 
    end type T_GOTMTrid

    !griflet start
    type T_GOTMimport
        double precision,  dimension(:), pointer    :: nn,ss,P,B,h
    end type T_GOTMimport
    !riflet end
    
    type T_Gotm
        double precision,  dimension(:), pointer    :: cmue1,cmue2
        double precision,  dimension(:), pointer    :: tkeo
        double precision,  dimension(:), pointer    :: as,an
        double precision,  dimension(:), pointer    :: xRf 

        !Instance of ModuleGOTM
        type(T_GOTMParameters        ), pointer :: Parameters
        type(T_GOTMexport            ), pointer :: Export
        !griflet start
        type(T_GOTMimport           ), pointer :: Impor
        !griflet end
        type(T_GOTMTrid                  ), pointer :: Trid
        integer                                     :: TID
         
    end type T_Gotm
    
    
    real   , parameter :: FillValueReal=-9.9e15
    integer, parameter :: FillValueInt =-9999999


!  Finished with all 'turbulence' related namelists.

!  Turbulent Kinetic Energy.
   integer, parameter   :: tke_local_eq=1
   integer, parameter   :: tke_keps=2 
   integer, parameter   :: tke_MY=3

!  Stability functions.
   integer, parameter   :: KanClay=1
   integer, parameter   :: BurBaum=2
   integer, parameter   :: CanutoA=3
   integer, parameter   :: CanutoB=4
   integer, parameter   :: KanClayQe=5 
   integer, parameter   :: BurBaumQe=6
   integer, parameter   :: CanutoAQe=7
   integer, parameter   :: CanutoBQe=8
   integer, parameter   :: Constant=9
   integer, parameter   :: MunkAnderson=10
   integer, parameter   :: SchumGerz=11
   integer, parameter   :: FluxRich=12

!  Length scale calculation.
   integer, parameter   :: diss_eq=8
   integer, parameter   :: length_eq=9
   integer, parameter   :: BougeaultAndre=6
   integer, parameter   :: ispra_length=6

!-----------------------------------------------------------------------

    contains

!griflet-----------------------------------------------------------------------
    subroutine clear_turbulence_parameters(ObjGotmparameters)
   
    !Arguments
    type(T_Gotmparameters   ), pointer :: ObjGOTMparameters

!    !Local variables
!    logical craig_banner,length_lim
!    logical qesmooth,flux_bdy
!    integer turb_method,tke_method,len_scale_method,stab_method,MY_length
!    integer iw_model   
!    double precision  const_num,const_nuh,k_min,L_min,eps_min
!    double precision kappa,Prandtl0,cm0,cm_craig,cw,galp
!    double precision ce1,ce2,ce3minus,ce3plus,sig_k
!    double precision sl,e1,e2,e3
!    double precision a1,a2,b1,b2,c2,c3,qeghmax,qeghmin,qeghcrit   
!    double precision alpha,klimiw,rich_cr,numiw,nuhiw,numshear
!    double precision c1,cde,craig_m,sig_e0,sig_e1
   
    ObjGOTMParameters%turb_method           = null_int
    ObjGOTMParameters%tke_method            = null_int
    ObjGOTMParameters%len_scale_method      = null_int
    ObjGOTMParameters%stab_method           = null_int
    ObjGOTMParameters%MY_length             = null_int
    ObjGOTMParameters%iw_model              = null_int

    ObjGOTMParameters%craig_banner = .false.
    ObjGOTMParameters%length_lim   = .false.
    ObjGOTMParameters%flux_bdy     = .false.
    ObjGOTMParameters%qesmooth     = .false.

    ObjGOTMParameters%c1           = null_real
    ObjGOTMParameters%cde          = null_real
    ObjGOTMParameters%craig_m      = null_real
    ObjGOTMParameters%sig_e0       = null_real
    ObjGOTMParameters%sig_e1       = null_real
    ObjGOTMParameters%const_num    = null_real
    ObjGOTMParameters%const_nuh    = null_real
    ObjGOTMParameters%k_min        = null_real
    ObjGOTMParameters%L_min        = null_real
    ObjGOTMParameters%eps_min      = null_real
    ObjGOTMParameters%kappa        = null_real
    ObjGOTMParameters%Prandtl0     = null_real
    ObjGOTMParameters%cm0          = null_real
    ObjGOTMParameters%cm_craig     = null_real
    ObjGOTMParameters%cw           = null_real
    ObjGOTMParameters%galp         = null_real
    ObjGOTMParameters%ce1          = null_real
    ObjGOTMParameters%ce2          = null_real
    ObjGOTMParameters%ce3minus     = null_real
    ObjGOTMParameters%ce3plus      = null_real
    ObjGOTMParameters%sig_k        = null_real
    ObjGOTMParameters%sl           = null_real
    ObjGOTMParameters%e1           = null_real
    ObjGOTMParameters%e2           = null_real
    ObjGOTMParameters%e3           = null_real
    ObjGOTMParameters%a1           = null_real
    ObjGOTMParameters%a2           = null_real
    ObjGOTMParameters%b1           = null_real
    ObjGOTMParameters%b2           = null_real
    ObjGOTMParameters%c2           = null_real
    ObjGOTMParameters%c3           = null_real
    ObjGOTMParameters%qeghmax      = null_real
    ObjGOTMParameters%qeghmin      = null_real
    ObjGOTMParameters%qeghcrit     = null_real
    ObjGOTMParameters%alpha        = null_real
    ObjGOTMParameters%klimiw       = null_real
    ObjGOTMParameters%rich_cr      = null_real
    ObjGOTMParameters%numiw        = null_real
    ObjGOTMParameters%nuhiw        = null_real
    ObjGOTMParameters%numshear     = null_real
    
    end subroutine clear_turbulence_parameters 
!griflet-----------------------------------------------------------------------

!griflet-----------------------------------------------------------------------
    subroutine copy_turbulence_parameters(targetP, sourceP)
   
    !Arguments (already allocated)
    type(T_Gotmparameters   ), pointer :: targetP, sourceP
   
    targetP%turb_method           = sourceP%turb_method
    targetP%tke_method            = sourceP%tke_method
    targetP%len_scale_method      = sourceP%len_scale_method
    targetP%stab_method           = sourceP%stab_method
    targetP%MY_length             = sourceP%MY_length
    targetP%iw_model              = sourceP%iw_model

    targetP%c1           = sourceP%c1
    targetP%cde          = sourceP%cde
    targetP%craig_m      = sourceP%craig_m
    targetP%sig_e0       = sourceP%sig_e0
    targetP%sig_e1       = sourceP%sig_e1
    targetP%craig_banner = sourceP%craig_banner
    targetP%length_lim   = sourceP%length_lim
    targetP%const_num    = sourceP%const_num
    targetP%const_nuh    = sourceP%const_nuh
    targetP%k_min        = sourceP%k_min
    targetP%L_min        = sourceP%L_min
    targetP%eps_min      = sourceP%eps_min
    targetP%kappa        = sourceP%kappa
    targetP%Prandtl0     = sourceP%Prandtl0
    targetP%cm0          = sourceP%cm0
    targetP%cm_craig     = sourceP%cm_craig
    targetP%cw           = sourceP%cw
    targetP%galp         = sourceP%galp
    targetP%ce1          = sourceP%ce1
    targetP%ce2          = sourceP%ce2
    targetP%ce3minus     = sourceP%ce3minus
    targetP%ce3plus      = sourceP%ce3plus
    targetP%sig_k        = sourceP%sig_k
    targetP%flux_bdy     = sourceP%flux_bdy
    targetP%sl           = sourceP%sl
    targetP%e1           = sourceP%e1
    targetP%e2           = sourceP%e2
    targetP%e3           = sourceP%e3
    targetP%a1           = sourceP%a1
    targetP%a2           = sourceP%a2
    targetP%b1           = sourceP%b1
    targetP%b2           = sourceP%b2
    targetP%c2           = sourceP%c2
    targetP%c3           = sourceP%c3
    targetP%qesmooth     = sourceP%qesmooth
    targetP%qeghmax      = sourceP%qeghmax
    targetP%qeghmin      = sourceP%qeghmin
    targetP%qeghcrit     = sourceP%qeghcrit
    targetP%alpha        = sourceP%alpha
    targetP%klimiw       = sourceP%klimiw
    targetP%rich_cr      = sourceP%rich_cr
    targetP%numiw        = sourceP%numiw
    targetP%nuhiw        = sourceP%nuhiw
    targetP%numshear     = sourceP%numshear

    end subroutine copy_turbulence_parameters
!griflet-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize the Turbulence Module parameters
!
! !INTERFACE:
   subroutine init_turbulence_parameters(ObjGotmparameters,namlst,fn,STAT)
!
! !DESCRIPTION:
!  Initialize all turbulence related stuff - reads a number of namelists
!  and allocates memory for turbulence related vectors.
!
! !USES:
!
! !INPUT PARAMETERS:
!Arguments-------------------------------------------------------------
   integer, optional, intent(OUT) :: STAT    
   integer, intent(in)  :: namlst
   character(len=*), intent(in) :: fn
!
! !INPUT/OUTPUT PARAMETERS:


! !OUTPUT PARAMETERS:
   type(T_Gotmparameters   ), pointer :: ObjGOTMparameters

! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
   logical craig_banner,length_lim
   logical qesmooth,flux_bdy
   integer turb_method,tke_method,len_scale_method,stab_method,MY_length
   integer iw_model   
   double precision  const_num,const_nuh,k_min,L_min,eps_min
   double precision kappa,Prandtl0,cm0,cm_craig,cw,galp
   double precision ce1,ce2,ce3minus,ce3plus,sig_k
   double precision sl,e1,e2,e3
   double precision a1,a2,b1,b2,c2,c3,qeghmax,qeghmin,qeghcrit   
   double precision alpha,klimiw,rich_cr,numiw,nuhiw,numshear
   !double precision c1,cde,craig_m,sig_e0,sig_e1
   
   namelist /turbulence/ turb_method,tke_method,len_scale_method,stab_method, &
                         craig_banner,length_lim,const_num,const_nuh,k_min,   &
                         L_min,eps_min
   namelist /turb_parameters/ kappa,Prandtl0,cm0,cm_craig,cw,galp
   namelist /keps/ ce1,ce2,ce3minus,ce3plus,sig_k,flux_bdy
   namelist /my/ sl,e1,e2,e3,MY_length
   namelist /stabfunc/ a1,a2,b1,b2,c2,c3,qesmooth,qeghmax,qeghmin,qeghcrit
   namelist /iw/ iw_model,alpha,klimiw,rich_cr,numiw,nuhiw,numshear
   integer STAT_,STAT_CALL

!
!EOP
!-----------------------------------------------------------------------
!BOC
     
     allocate (ObjGOTMparameters, STAT = STAT_CALL)
     if (STAT_CALL .NE. SUCCESS_)                                       &
                stop 'Subroutine init_turbulence_parameters - ModuleGOTM. Allocation of OBJGOTMParameters.'
     
     write(0,*) '   ', 'init_turbulence'

     open(namlst,file=fn,status='old',action='read',err=80)
     write(0,*) '       ', 'reading turbulence namelists..'
     read(namlst,nml=turbulence,err=81)
     read(namlst,nml=turb_parameters,err=82)
     read(namlst,nml=keps,err=83)
     read(namlst,nml=my,err=84)
     read(namlst,nml=stabfunc,err=85)
     read(namlst,nml=iw,err=86)
     close (namlst)

     ObjGOTMParameters%c1 =(1.-B1**(-1./3.)/A1-6*A1/B1)/3. !See Kantha & Clayson 1994, eq. (23)
     ObjGOTMParameters%cde=cm0*cm0*cm0
     ObjGOTMParameters%craig_m=sqrt(1.5*cm_craig**2*sig_k/kappa**2)
     ObjGOTMParameters%sig_e0=(4./3.*ObjGOTMParameters%craig_m+1.)*(ObjGOTMParameters%craig_m+1.)*kappa**2/(ce2*cm_craig**2)
     ObjGOTMParameters%sig_e1= kappa*kappa*cm0/(ce2-ce1)/ObjGOTMParameters%cde

     ObjGOTMParameters%turb_method  = turb_method
     ObjGOTMParameters%tke_method   = tke_method
     ObjGOTMParameters%len_scale_method    = len_scale_method
     ObjGOTMParameters%stab_method  = stab_method
     ObjGOTMParameters%craig_banner = craig_banner
     ObjGOTMParameters%length_lim   = length_lim
     ObjGOTMParameters%const_num    = const_num
     ObjGOTMParameters%const_nuh    = const_nuh
     ObjGOTMParameters%k_min        = k_min
     ObjGOTMParameters%L_min        = L_min
     ObjGOTMParameters%eps_min      = eps_min
     ObjGOTMParameters%kappa        = kappa
     ObjGOTMParameters%Prandtl0     = Prandtl0
     ObjGOTMParameters%cm0          = cm0
     ObjGOTMParameters%cm_craig     = cm_craig
     ObjGOTMParameters%cw           = cw
     ObjGOTMParameters%galp         = galp
     ObjGOTMParameters%ce1          = ce1
     ObjGOTMParameters%ce2          = ce2
     ObjGOTMParameters%ce3minus     = ce3minus
     ObjGOTMParameters%ce3plus      = ce3plus
     ObjGOTMParameters%sig_k        = sig_k
     ObjGOTMParameters%flux_bdy     = flux_bdy
     ObjGOTMParameters%sl           = sl
     ObjGOTMParameters%e1           = e1
     ObjGOTMParameters%e2           = e2
     ObjGOTMParameters%e3           = e3
     ObjGOTMParameters%MY_length    = MY_length
     ObjGOTMParameters%a1           = a1
     ObjGOTMParameters%a2           = a2
     ObjGOTMParameters%b1           = b1
     ObjGOTMParameters%b2           = b2
     ObjGOTMParameters%c2           = c2
     ObjGOTMParameters%c3           = c3
     ObjGOTMParameters%qesmooth     = qesmooth
     ObjGOTMParameters%qeghmax      = qeghmax
     ObjGOTMParameters%qeghmin      = qeghmin
     ObjGOTMParameters%qeghcrit     = qeghcrit
     ObjGOTMParameters%iw_model     = iw_model
     ObjGOTMParameters%alpha        = alpha
     ObjGOTMParameters%klimiw       = klimiw
     ObjGOTMParameters%rich_cr      = rich_cr
     ObjGOTMParameters%numiw        = numiw
     ObjGOTMParameters%nuhiw        = nuhiw
     ObjGOTMParameters%numshear     = numshear
     write(0,*) '       ', 'done.'

     STAT_=SUCCESS_

     if (present(STAT))                                                    &
            STAT = STAT_
     

   return
80 write(0,*) 'FATAL ERROR: ', 'I could not open "gotmturb.inp"'
   stop 'init_turbulence'
81 write(0,*) 'FATAL ERROR: ', 'I could not read "turbulence" namelist'
   stop 'init_turbulence'
82 write(0,*) 'FATAL ERROR: ', 'I could not read "turb_parameters" namelist'
   stop 'init_turbulence'
83 write(0,*) 'FATAL ERROR: ', 'I could not read "keps" namelist'
   stop 'init_turbulence'
84 write(0,*) 'FATAL ERROR: ', 'I could not read "my" namelist'
   stop 'init_turbulence'
85 write(0,*) 'FATAL ERROR: ', 'I could not read "stabfunc" namelist'
   stop 'init_turbulence'
86 write(0,*) 'FATAL ERROR: ', 'I could not read "iw" namelist'
   stop 'init_turbulence'
   

   end subroutine init_turbulence_parameters 
!EOC

!BOP
!
! !IROUTINE: Initialize the Turbulence Module variables
!
! !INTERFACE:
   subroutine init_turbulence(ObjGotm, ObjGOTMParameters, KLB, KUB, STAT)
!
! !DESCRIPTION:
!  Allocates memory for turbulence related vectors.
!
!  Arguments
   type(T_GOTM), pointer            :: ObjGOTM
   type(T_GOTMParameters), pointer  :: ObjGOTMParameters
   integer                          :: KLB,KUB
   integer, optional, intent(OUT)   :: STAT   
!  LOCAL VARIABLES:
   integer STAT_

   STAT_=UNKNOWN_

   allocate(ObjGOTM%tkeo(KLB:KUB))
   allocate(ObjGOTM%cmue1(KLB:KUB))
   allocate(ObjGOTM%cmue2(KLB:KUB))
   allocate(ObjGOTM%xRF(KLB:KUB))
   allocate(ObjGOTM%an(KLB:KUB))
   allocate(ObjGOTM%as(KLB:KUB))

   ObjGOTM%tkeo     = ObjGOTMParameters%k_min
   ObjGOTM%cmue1    = FillValueReal
   ObjGOTM%cmue2    = FillValueReal
   ObjGOTM%xRF      = FillValueReal
   ObjGOTM%an       = FillValueReal
   ObjGOTM%as       = FillValueReal

   allocate(ObjGOTM%Export)  
   allocate(ObjGOTM%Export%tke(KLB:KUB))
   allocate(ObjGOTM%Export%eps(KLB:KUB))
   allocate(ObjGOTM%Export%L(KLB:KUB))
   allocate(ObjGOTM%Export%num(KLB:KUB))
   allocate(ObjGOTM%Export%nuh(KLB:KUB))

   ObjGOTM%Export%tke   = ObjGOTMParameters%k_min 
   ObjGOTM%Export%eps   = ObjGOTMParameters%eps_min
   ObjGOTM%Export%L     = ObjGOTMParameters%L_min
   ObjGOTM%Export%num   = 0.
   ObjGOTM%Export%nuh   = 0.

   allocate(ObjGOTM%Trid)
   allocate(ObjGOTM%Trid%au(KLB:KUB))
   allocate(ObjGOTM%Trid%bu(KLB:KUB))
   allocate(ObjGOTM%Trid%cu(KLB:KUB))
   allocate(ObjGOTM%Trid%du(KLB:KUB))
   allocate(ObjGOTM%Trid%ru(KLB:KUB))
   allocate(ObjGOTM%Trid%qu(KLB:KUB))

   ObjGOTM%Trid%au = FillValueReal
   ObjGOTM%Trid%bu = FillValueReal
   ObjGOTM%Trid%cu = FillValueReal
   ObjGOTM%Trid%du = FillValueReal
   ObjGOTM%Trid%ru = FillValueReal
   ObjGOTM%Trid%qu = FillValueReal

   !griflet start
   allocate(ObjGOTM%Impor)
   allocate(ObjGOTM%Impor%nn(KLB:KUB))    
   allocate(ObjGOTM%Impor%ss(KLB:KUB))
   allocate(ObjGOTM%Impor%P(KLB:KUB))
   allocate(ObjGOTM%Impor%B(KLB:KUB))
   allocate(ObjGOTM%Impor%h(KLB:KUB))
   
   ObjGOTM%Impor%nn = FillValueReal
   ObjGOTM%Impor%ss = FillValueReal
   ObjGOTM%Impor%P  = FillValueReal
   ObjGOTM%Impor%B  = FillValueReal
   ObjGOTM%Impor%h  = FillValueReal

   STAT_=SUCCESS_

   if (present(STAT))                                                    &
            STAT = STAT_

   end subroutine init_turbulence
   !-------------------------------------------------------------------------

   !-------------------------------------------------------------------------
   !griflet: This subroutine was missing and contributing to a slight memory leak!!!
   subroutine kill_turbulence(ObjGotm, STAT)

    integer, optional, intent(OUT)  :: STAT
    type(T_GOTM), pointer           :: ObjGOTM
    
    integer                         :: STAT_

    deallocate(ObjGOTM%Export%tke)
    deallocate(ObjGOTM%Export%eps)
    deallocate(ObjGOTM%Export%L)
    deallocate(ObjGOTM%Export%num)
    deallocate(ObjGOTM%Export%nuh)
    deallocate(ObjGOTM%Export)  

    deallocate(ObjGOTM%Trid%au)
    deallocate(ObjGOTM%Trid%bu)
    deallocate(ObjGOTM%Trid%cu)
    deallocate(ObjGOTM%Trid%du)
    deallocate(ObjGOTM%Trid%ru)
    deallocate(ObjGOTM%Trid%qu)    
    deallocate(ObjGOTM%Trid)

    deallocate(ObjGOTM%Impor%nn)
    deallocate(ObjGOTM%Impor%ss)
    deallocate(ObjGOTM%Impor%P)
    deallocate(ObjGOTM%Impor%B)
    deallocate(ObjGOTM%Impor%h)    
    deallocate(ObjGOTM%Impor)    

    deallocate(ObjGOTM%tkeo)
    deallocate(ObjGOTM%cmue1)
    deallocate(ObjGOTM%cmue2)
    deallocate(ObjGOTM%xRF)
    deallocate(ObjGOTM%an)
    deallocate(ObjGOTM%as)
              
    STAT_=SUCCESS_

    if (present(STAT))                                                    &
        STAT = STAT_
            
    end subroutine kill_turbulence
   !-------------------------------------------------------------------------

   !-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_turbulence() - do the vertical 1D turbulence
!
! !INTERFACE:
   subroutine do_turbulence(ObjGOTM,nlev,dt,depth,u_taus,u_taub,z0s,z0b,  &
                            h,NN,SS,P,B)

!
! !DESCRIPTION:
!
! !USES:
!
! !INPUT PARAMETERS:
   type(T_GOTM), pointer     :: ObjGOTM
   integer, intent(in)  :: nlev
   double precision, intent(in) :: dt,depth,u_taus,u_taub,z0s,z0b
   double precision, intent(in) :: h(0:nlev)
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev),P(0:nlev),B(0:nlev)
! !OUTPUT PARAMETERS:
!
! !BUGS:
!
! !SEE ALSO: 
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
!
!EOP
!-------------------------------------------------------------------------
!BOC
        
!        select case (turb_method)
!         case (2)
!            call production(nlev,alpha,num,nuh,P,B)
            call stabilityfunctions(ObjGOTM,nlev,NN,SS)
            call do_tke(ObjGOTM,nlev,dt,u_taus,u_taub,z0s,h,NN,SS,P,B)
            call lengthscale(ObjGOTM,nlev,dt,z0b,z0s,u_taus,u_taub,depth,h,NN,P,B)
            call kolpran(ObjGOTM,nlev,u_taub,u_taus,z0b,z0s)
!         case default
!        end select 
        call internal_wave(ObjGOTM,nlev,NN,SS)

   end subroutine do_turbulence
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Stability Functions: 
!
! !INTERFACE:
   subroutine stabilityfunctions(ObjGOTM,nlev,NN,SS)
!
! !DESCRIPTION:
!  Based on user input - this routine calls the appropriate stability
!  method. 
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
   integer      :: i 
   double precision     :: LLk

   include 'GOTMVariables_in.F90'
!
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=0,nlev
      LLk=L(i)*L(i)/tke(i)
      as(i)=LLk*SS(i)
      an(i)=LLk*NN(i)
   end do
   
   select case(stab_method)
      case(KanClay)
         call cmue_kc(ObjGOTM,nlev)
      case(BurBaum)
!         call cmue_bb(ObjGOTM,nlev)
          stop 'This option not available'
      case(CanutoA)
         call cmue_ca(ObjGOTM,nlev)
      case(CanutoB)
         call cmue_cb(ObjGOTM,nlev)
      case(KanClayQe)
         call cmue_kcqe(ObjGOTM,nlev)
      case(BurBaumQe)
!         call cmue_bbqe(ObjGOTM,nlev)
          stop 'This option not available'
      case(CanutoAQe)
         call cmue_caqe(ObjGOTM,nlev)
      case(CanutoBQe)
         call cmue_cbqe(ObjGOTM,nlev)
      case(Constant)
         cmue1=cm0
         cmue2=cm0/Prandtl0
      case(MunkAnderson)
         call cmue_ma(ObjGOTM,nlev,NN,SS)
      case(SchumGerz)
         call cmue_sg(ObjGOTM,nlev,NN,SS)
      case(FluxRich)
         call cmue_rf(ObjGOTM,nlev,NN,SS)
      case default
   end select

   cmue1(0)=cmue1(1)
   cmue1(nlev)=cmue1(nlev-1)
   cmue2(0)=cmue2(1)
   cmue2(nlev)=cmue2(nlev-1)

   !include 'GOTMVariables_out.F90'
 
   end subroutine stabilityfunctions 
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Turbulent Kinetic Energy Calculation.
!
! !INTERFACE:
   subroutine do_tke(ObjGOTM,nlev,dt,u_taus,u_taub,z0s,h,NN,SS,P,B)
!
! !DESCRIPTION:
!  Based on user input - this routines calls the appropriate routine for
!  calculating the turbulent kinetic energy.
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: dt
   double precision, intent(in) :: u_taus,u_taub,z0s
   double precision, intent(in) :: h(0:nlev)
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
   double precision, intent(in) :: P(0:nlev),B(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision     :: numtke(0:nlev)

   include 'GOTMVariables_in.F90'
!
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (tke_method)
      case(tke_local_eq)
         call tkealgebraic(ObjGOTM,nlev,u_taus,u_taub,NN,SS)
      case(tke_keps)
         do i=1,nlev-1
            numtke(i)=num(i)/sig_k
         end do
         call tkeeq(ObjGOTM,nlev,dt,u_taus,u_taub,z0s,h,P,B,numtke)
      case(tke_MY)
         do i=1,nlev-1
            numtke(i)=Sl*sqrt(2.*tke(i))*L(i)
         end do
         call tkeeq(ObjGOTM,nlev,dt,u_taus,u_taub,z0s,h,P,B,numtke)
      case default
   end select
   
   !include 'GOTMVariables_out.F90'
      
   end subroutine do_tke 
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Length Scales: 
!
! !INTERFACE:
   subroutine lengthscale(ObjGOTM,nlev,dt,z0b,z0s,u_taus,u_taub,depth,h,NN,P,B)
!
! !DESCRIPTION:
!  Calls different subroutines that calculate the lengthscale $L$  
!  and the dissipation rate $\epsilon$ with different methods.
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: dt
   double precision, intent(in) :: z0b,z0s
   double precision, intent(in) :: u_taus,u_taub
   double precision, intent(in) :: depth
   double precision, intent(in) :: h(0:nlev),NN(0:nlev)
   double precision, intent(in) :: P(0:nlev),B(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision     :: numtke(0:nlev)

   include 'GOTMVariables_in.F90'
!
!EOP
!-----------------------------------------------------------------------
!BOC
   select case(len_scale_method)
      case(diss_eq)
         call dissipationeq(ObjGOTM,nlev,dt,z0b,z0s,u_taus,u_taub,h,NN,P,B)
      case(length_eq)
         do i=1,nlev-1
            numtke(i)=Sl*sqrt(2.*tke(i))*L(i)
         end do
         call lengthscaleeq(ObjGOTM,nlev,dt,z0b,z0s,depth,h,NN,P,B,numtke)
      case(BougeaultAndre)      ! Bougeault and Andre (1986)
         call potentialml(ObjGOTM,nlev,z0b,z0s,h,depth,NN)
      case default
         call algebraiclength(ObjGOTM,len_scale_method,nlev,z0b,z0s,depth,h,NN)
   end select

   !include 'GOTMVariables_out.F90'
 
   end subroutine lengthscale
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Length scale from dissipation equation. 
! 
! !INTERFACE:
   subroutine dissipationeq(ObjGOTM,nlev,dt,z0b,z0s,u_taus,u_taub,h,NN,P,B)
!
! !DESCRIPTION:
!  This subroutine calculates the dissipation rate in the framework
!  of the k-epsilon model:
!
!  \begin{equation}\label{eps_eq}
!  \partial_t \varepsilon -
!  \partial_z(\nu_{\varepsilon}\partial_z \varepsilon) =
!  \frac{\varepsilon}{k} \left(c_{\varepsilon 1}P + c_{\varepsilon 3}B - c_{\varepsilon 2}\varepsilon \right).
!  \end{equation}
!
!  As boundary condtions a choice between Dirichlet (flux-bdy=.false.)
!  and Neumann flux conditions (flux-bdy=.true.) has to be made.   
!
!  Dirichlet conditions:
!
!  \begin{equation}\label{Standardeps}
!  \varepsilon =
!  \left( c_{\mu}^0 \right)^3 \frac{k^{3/2}}{\kappa (\tilde z + z_0)}.
!  \end{equation}
!
!  Neumann flux conditions:
!
!  \begin{equation}\label{Fluxeps}
!  \frac{\nu_t}{\sigma_{\varepsilon}} \partial_{\tilde z} \varepsilon =
!  -\left( c_{\mu}^0 \right)^3
!  \frac{\nu_t}{\sigma_{\varepsilon}}
!  \frac{k^{3/2}}{\kappa (\tilde z + z_0)^2}. 
!  \end{equation}
!
!  If the {\it Craig and Banner} [1994] and the {\it Craig} [1996]
! surface wave breaking theory is applied, then the following extended
! surface boundary condition for $\eps$ is used:
!  
!  
!  \begin{equation}\label{NewBC}
!  -\frac{\nu_t}{\sigma_{\eps}}\partial_z \eps=
!  -\frac{\nu_t}{\sigma_{\eps}}
!  (c_{\mu}^0)^{3/4}
!  \frac{\frac{3}{2}\frac{\sigma_k(c_{\mu}^0)^{3/4}}{c_{\mu}}c_w(u_*^s)^3
!  +\kappa k^{3/2}}
!  {\kappa^2(z'+z_0^s)^2}.
!  \end{equation}
!
!  This has been constructed with the aid of the analytical solution
!  for the dissipation rate $\eps$ as suggested by {\it Craig} [1996]:
!
!  \begin{equation}\label{epsanalyt}
!  \eps=\frac{(u_s^*)^3}{\kappa (z'+z_0^s)}
!  \left[a+\left(\frac{3\sigma_k}{2}\right)^{1/2}
!  c_{\mu}^{1/4}c_w\left(\frac{z'+z_0^s}{z_0^s}\right)^{-m}\right]. 
!  \end{equation}
!
!  The turbulent Schmidt number $\sigma_{\eps}$ for the dissipation
!  rate $\eps$ is a linear interpolation between
!
!  \begin{equation}\label{sigma0}
!  \sigma_{\eps}=
!  \sigma_{\eps 0}
!  =\left(\frac43m+1\right)(m+1)\frac{\kappa^2}{c_2c_{\mu}^{1/2}}
!  \approx 2.4
!  \end {equation}
!  
!  for $(P+B)/\eps=0$ and 
!  
!  \begin{equation}\label{sigma1}
!  \sigma_{\eps}=\sigma_{\eps 1}=\frac{\kappa^2}{c_{\mu}^{1/2}(c_2-c_1)}
!  \approx 1.111 
!  \end{equation}
!
!  for $(P+B)/\eps=1$. For more details, see {\it Burchard} [2000]. 
!
!  At the end of the subroutine, the Galperin et al. [1988] limitation
!  and the calculation of the macro length scale is carried out. 
!
! !USES:
!   use mtridiagonal
!   use turbulence, ONLY: ce1,ce2,ce3plus,ce3minus,cde,kappa
!   use turbulence, ONLY: sig_e0,sig_e1,sig_k
!   use turbulence, ONLY: cm0,galp,flux_bdy,length_lim
!   use turbulence, ONLY: num,eps,L,tkeo,tke,eps_min,L_min
!   use turbulence, ONLY: cmue1,craig_banner,cw
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: dt
   double precision, intent(in) :: z0b,z0s
   double precision, intent(in) :: u_taus,u_taub
   double precision, intent(in) :: h(0:nlev)
   double precision, intent(in) :: P(0:nlev),B(0:nlev)
   double precision, intent(in) :: NN(0:nlev)

   
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   double precision         :: avh(0:nlev),flux(0:nlev)
   double precision         :: pminus(0:nlev),pplus(0:nlev)
   double precision         :: prod,buoyan,diss
   double precision         :: cee3,epslim 
   double precision         :: peps,sig_e(0:nlev) 
   double precision         :: kk
   integer      :: i

   include 'GOTMVariables_in.F90'

! 
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1 
      flux(i)=0. 
   end do

!  Determination of the turbulent Schmidt number for the dissipation rate:

   if (craig_banner) then     ! With wave breaking
      sig_e(nlev)=sig_e0
      do i=1,nlev-1
         peps=(P(i)+B(i))/eps(i)
         if (peps .gt. 1.) peps=1.
         sig_e(i)=peps*sig_e1+(1.-peps)*sig_e0
      end do
      sig_e(0)=sig_e1
   else                        ! No wave breaking
      do i=0,nlev
         sig_e(i)=sig_e1
      end do
   end if 

   do i=1,nlev
      avh(i)=0.5*(num(i-1)/sig_e(i-1)+num(i)/sig_e(i))
   end do

   if (flux_bdy) then
      flux(1  )=avh(1)*cde*(tkeo(1  )**1.5)/(kappa*(z0b+0.5*h(1))**2.)
      kk=tkeo(nlev-1)
      if (craig_banner) then
         flux(nlev-1)=cmue1(nlev-1)*sqrt(kk)*kappa*(0.5*h(nlev)+z0s)   &
                /sig_e(nlev-1)*cde*(kk**1.5+                           &
                1.5/kappa/cmue1(nlev-1)*sig_k*cw*u_taus**3)/           &
                      (kappa*(z0s+0.5*h(nlev))**2.)
      else
         flux(nlev-1)=cmue1(nlev-1)*sqrt(kk)*kappa*(0.5*h(nlev)+z0s)   &
                /sig_e(nlev-1)*cde*kk**1.5/                            &
                      (kappa*(z0s+0.5*h(nlev))**2.)
! A bug in the previous two lines has been found 
! by Patrick Luyten, MUMM, Belgium. kappa had been squared as well before.
! See the GOTM report, 1999 for the correct mathematical formulation. 
      end if 
      avh(1)=0.
      avh(nlev)=0.
   else       ! Not needed any more for Dirichlet conditions, only 
! kept for "historical" reasons, see Burchard et al. [1998]. 
      avh(1)=u_taub**4*2/sig_e1/(eps(0)+eps(1))
      avh(nlev)=u_taus**4*2/sig_e1/(eps(nlev)+eps(nlev-1))
   end if 

   do i=1,nlev-1
      if (B(i).gt.0) then
         cee3=ce3plus 
      else
         cee3=ce3minus 
      end if
      prod=ce1*eps(i)/tkeo(i)*P(i)
      buoyan=cee3*eps(i)/tkeo(i)*B(i)
      diss=ce2*eps(i)*eps(i)/tkeo(i)
      if (prod+buoyan.gt.0) then
         pplus(i)=prod+buoyan
         pminus(i)=diss
      else
         pplus(i)=prod
         pminus(i)=diss-buoyan
      end if
   end do
   do i=1,nlev-1
      au(i)=-2.*dt*avh(i)/(h(i)+h(i+1))/h(i)
      cu(i)=-2.*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
      bu(i)=1.-au(i)-cu(i)+pminus(i)*dt/eps(i)
      du(i)=(1+pplus(i)*dt/eps(i))*eps(i)+flux(i)*dt/(0.5*(h(i)+h(i+1)))
   end do

   if (flux_bdy) then
      call Tridiagonal(ObjGOTM,nlev,1,nlev-1,eps)
      eps(0) = cde*sqrt(tke(0)*tke(0)*tke(0))/kappa/z0b 
      eps(nlev) = cde*sqrt(tke(nlev)*tke(nlev)*tke(nlev))/kappa/z0s 
   else
      cu(1)=0.          ! lower boundary, one grid-box from bed
      bu(1)=1.
      du(1)=cde*tke(1)**1.5/kappa/(z0b+h(1)) 

      bu(nlev-1)=1.    ! upper boundary, one grid-box from surface 
      au(nlev-1)=0.
      du(nlev-1)=cde*tke(nlev-1)**1.5/kappa/(z0s+h(nlev)) 

      call Tridiagonal(ObjGOTM,nlev,1,nlev-1,eps)
   end if 
 
   do i=0,nlev
      if ((NN(i).gt.0).and.(length_lim)) then
         epslim=cm0*cm0*cm0/sqrt(2.)/galp*tke(i)*sqrt(NN(i)) 
      else
         epslim=eps_min 
      end if
      if (eps(i).lt.epslim) eps(i)=epslim
      L(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/eps(i)
      if (L(i).lt.L_min) L(i)=L_min
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine dissipationeq
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 'Simple' algebraic length scales.
!
! !INTERFACE:
   subroutine algebraiclength(ObjGOTM,method,nlev,z0b,z0s,depth,h,NN) 
!
! !DESCRIPTION:
!  These subroutine computes the vertical profile of the mixing lengthscale
!  from an algebraic relation.
!
!  1) Parabola
!  \begin{equation}
!  l=\kappa \frac {l_b l_s} {l_b+l_s}
!  \end{equation}
!  2) Triangle
!  \begin{equation}
!  l=\kappa max(l_b,l_s)
!  \end{equation}
!  3) Distorted Parabola. Robert-Ouellet
!  \begin{equation}
!  l=\kappa l_b (1-\frac {l_b} {D})
!  \end{equation}
!  4) Xing, parabolic with exponential $d_b$
!  \begin{equation}
!  l_b=\kappa l_b e^{-\beta l_b}
!  \end{equation}
!  5) Blackadar (in the presence of two boundaries)
!  \begin{equation}
!  l=\kappa \frac {1} {\frac {1} {l_b}+\frac {1} {l_s}+\frac {1} {l_a}} \\
!  l_a=\gamma_0 \frac {\int_{0}^{D} k^{1/2} zdz} {\int_{0}^{D} k^{1/2} dz}
!  \end{equation}
!  6) see ispralength.f
!
!  At the end of the subroutine, the dissipation rate is calculated using:
!
!  \begin{equation}\label{DefDissip}
!  \varepsilon = (c_{\mu}^0)^3 \frac{k^{3/2}}{L}. 
!  \end{equation}
!
! !USES:
!   use turbulence, ONLY: L,eps,tkeo,tke,L_min,eps_min,cde,galp,kappa,length_lim

!
! !INPUT PARAMETERS:
   integer, intent(in)  :: method,nlev
   double precision, intent(in) :: z0b,z0s
   double precision, intent(in) :: depth,h(0:nlev),NN(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision         :: ds,db,dbxing
   double precision         :: beta,gamma,La,La_up,La_down
   double precision         :: Lcrit

   integer, parameter   :: Parabola=1
   integer, parameter   :: Triangle=2
   integer, parameter   :: Xing=3
   integer, parameter   :: RobertOuellet=4
   integer, parameter   :: Blackadar=5
   integer, parameter   :: ispra_length=7

   include 'GOTMVariables_in.F90'
!
!EOP
!-----------------------------------------------------------------------
!BOC
   db=0.0
!
! Parabola shape
!
   select case (method)
      case(parabola)
         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
            L(i)=kappa*(ds+z0s)*(db+z0b)/(ds+db+z0b+z0s)
         end do
      case(triangle)
         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
            L(i)=kappa*min(ds+z0s,db+z0b)
         end do
!  
! Xing and Davies (1995). 
! Modification of parabolic mixing length. db changes:
!
      case(Xing)
         beta=-2.
         do i=1,nlev-1
           db=db+h(i)
           ds=depth-db
           dbxing=db*exp(-beta*db)
           L(i)=kappa*(ds+z0s)*(dbxing+z0b)/(ds+dbxing+z0s+z0b)
         end do
!
! Robert and Ouellet(1987). Similar to parabolic
!
      case(RobertOuellet)
         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
            L(i)= kappa*(db+z0b)*sqrt((ds+z0s)/(ds+db+z0b+z0s)) 
            L(i)= kappa*(db+z0b)*sqrt(1-db/(ds+db)) 
         end do
!
! Blackadar (1962). 
! In the form suggested by Luyten et al. (1996) for two boundary layers.
!
      case(Blackadar)
         La_up=0.
         La_down=0.
         do i=1,nlev-1
            db=db+h(i) 
            La_up=La_up+sqrt(tkeo(i))*(db+z0b)*h(i)
            La_down=La_down+sqrt(tkeo(i))*h(i) 
         end do
         gamma=0.2
         La=gamma*La_up/La_down
         db=0.0
         do i=1,nlev-1
            db=db+h(i)
            ds=depth-db
             L(i)=1/(1/(kappa*(ds+z0s))+1/(kappa*(db+z0b))+1/La)
         end do
!
!  Ispramix
!
      case(ispra_length)
         call ispralength(ObjGOTM,nlev,NN,h,depth) 
      case default
   end select

! Boundary conditions for L
!
   L(0)=kappa*z0b
   L(nlev)=kappa*z0s 
 
   do i=0,nlev
      if ((NN(i).gt.0).and.(length_lim)) then
         Lcrit=sqrt(2*galp*galp*tke(i)/NN(i))
         if (L(i).gt.Lcrit) L(i)=Lcrit
      end if
      if (L(i).lt.L_min) L(i)=L_min
      eps(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/L(i)
      if (eps(i).lt.eps_min) eps(i)=eps_min
   end do  

   !include 'GOTMVariables_out.F90'

   end subroutine algebraiclength 
!EOC

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ispralength() - length scale from ISPRAMIX
!
! !INTERFACE:
      subroutine Ispralength(ObjGOTM,nlev,NN,h,depth)
!
! !DESCRIPTION:
!  This subroutine calculates the 
!  lengthscale used in the ISPRAMIX model of JRC, Ispra, Italy
!  (Eifler and Schrimpf [1992] and Demirov et al. [1998]).
!  $L$ in both mixed layers is obtained from a {\it Blackadar} [1962]
!  type formula:
!  \begin{equation}\label {Lmixed}
!  L=\frac {\kappa \tilde z} {1+\frac {\kappa \tilde z} {c_2 \cdot h_m}}
!  (1-R_f)^e
!  \end{equation}
!  where $\tilde z$
!  is the distance from the interface (surface or bottom). The
!  fraction in (\ref{Lmixed})
!  predicts an approximation to a linear behavior of $L$ near boundaries 
!  and a value proportional to the thickness of the mixed
!  layer far from the interface, $L=c_2 h_m$, where $c_2=0.065$
!  is estimated from experimental data as discussed in
!  {\it Eifler and Schrimpf} [1992].
!  The factor $(1-R_f)$, with the flux Richardson
!  number $R_f=-B/P$, accounts for the effect
!  of stratification on the length scale.
!  The parameter $e$ is here a tuning parameter
!  (pers.\ comm.\ Walter Eifler, JRC, Ispra, Italy)
!  which is usually set to $e=1$.
!
! !USES:
!   use turbulence, ONLY: L,tke,L_min,xRF,kappa
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: NN(0:nlev)
   double precision, intent(in) :: h(0:nlev),depth
!
! !OUTPUT PARAMETERS:
!
! !BUGS:
!
! !SEE ALSO: 
!  potentialml(), algebraiclength()   
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
!  Original author(s): GOTM Team
!
!  $Log$
!
! !LOCAL VARIABLES:
  integer       :: i,SLind,BLind,Index,Index2
  double precision      :: hms,hmb,db,ds
  double precision      :: kml,c2_i,c3_i

  include 'GOTMVariables_in.F90'
!
!EOP
!-------------------------------------------------------------------------
!BOC
   kml   = 1.e-5
   c2_i  = 0.065

!  Calculation of surface mixed layer depth 
   hms=0.
   SLind=1
   do i=nlev,1,-1
      hms=hms+h(i)
      if (tke(i).le.kml) then 
         SLind=i
         goto 500
      end if   
   end do
500  continue
!  Calculation of bottom mixed layer depth
   hmb=0.
   BLind=nlev
   do i=1,nlev
      hmb=hmb+h(i)
      if (tke(i).le.kml) then
         BLind=i
         goto 501 
      end if
   end do
501  Continue

! If there is no point where k < kml, the water column is assumed to be mixed.
   if (BLind.gt.SLind) then 
      hms=0.5*depth 
      hmb=0.5*depth
      BLind=int(nlev/2)
      SLind=int(nlev/2)+1 
   endif

! Calculation of mixing length in bottom layer 
   db=0.
   do i=1,BLind 
      db=db+h(i)
      L(i)=kappa*db/(1.+kappa*db/(c2_i*hmb+L_min))*xRf(i)**3 
      if (L(i).lt.L_min) L(i)=L_min
   end do 

! Calculation of mixing length in surface layer
   ds=h(nlev)
   do i=nlev-1,SLind,-1
      ds=ds+h(i)
      L(i)=kappa*ds/(1.+kappa*ds/(c2_i*hms+L_min))*xRf(i)**3
      if (L(i).lt.L_min) L(i)=L_min
   end do

! Calculation of mixing length in the intermediate region
       
   c3_i=L(SLind)*sqrt(NN(SLind)/tke(SLind))
   if (c3_i.lt.1e-10) c3_i=0.
   Index=Slind-1
   do i=SLind-1,BLind+1,-1
      if (NN(i).le.0.) then
         L(i)=L_min
      else
         L(i)=max(c3_i*sqrt(tke(i)/NN(i)),L_min)
         if (L(i).gt.L(SLind)) L(i)=L(SLind) 
      endif
      if (L(i).eq.L_min) then
         Index=i 
         goto 503
      end if
   end do
503  continue
   c3_i=L(BLind)*sqrt(NN(BLind)/tke(BLind))
   if (c3_i.lt.1e-10) c3_i=0.
   Index2=BLind+1  
   do i=BLind+1,Index
      if (NN(i).le.0.) then
         L(i)=L_min
      else
         L(i)=max(c3_i*sqrt(tke(i)/NN(i)),L_min)
         if(L(i).gt.L(BLind)) L(i)=L(BLind) 
      endif
      if (L(i).eq.L_min) then
         Index2=i 
         goto 504
      end if
   end do 
504  continue 
   do i=Index2+1,Index-1
      L(i)=L_min
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine Ispralength
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Internal Waves. 
!
! !INTERFACE:
   subroutine internal_wave(ObjGOTM,nlev,NN,SS) 
!
! !DESCRIPTION:
!  Imposes eddy viscosity and diffisivity characteristic 
!  of internal wave activity and shear instability when there is extinction 
!  of turbulence as suggested by Kantha and Clayson [1994]. 
!  In this case, these new num and nuh 
!  are used instead of those computed with the model.
!
!  When k is small (extinction of turbulence, diagnosed by $k<klimiw$), 
!  $\nu_t$ and $\nu'_t$ are set to empirical values typical 
!  in the presence of internal wave activity (IW) and shear 
!  instability (SI). 
!  {\large
!  \begin{equation}
!  \nu_t=(\nu_t)^{IW}+(\nu_t)^{SI}, \quad
!  \nu_t'=(\nu_t')^{IW}+(\nu'_t)^{SI}
!  \end{equation}
!  \begin{equation}
!  (\nu_t)^{IW}=10^{-4}, \quad         
!  (\nu'_t)^{IW}=5 10^{-5}
!  \end{equation} 
!  \begin{eqnarray}
!  (\nu_t)^{SI}=(\nu_t')^{SI}=0, & R_i>0.7 \\
!  (\nu_t)^{SI}=(\nu_t')^{SI}=5 10^{-3} \left[1-\left(\frac {R_i} 
!  {0.7}\right)^2\right]^3, & 0<R_i<0.7 \\
!  (\nu_t)^{SI}= (\nu_t')^{SI}=5 10^{-3}, & R_i < 0
!  \end{eqnarray}
! 
!  The unit of all diffusivities is $m^2 s^{-1}$
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
   double precision         :: rich(0:nlev)  
   double precision         :: rich2,pot,x      
   integer      :: i

   include 'GOTMVariables_in.F90'

!
!EOP
!-----------------------------------------------------------------------
!BOC
   select case (iw_model)
      case (0)
      case (1)
      case (2)
         rich2 = rich_cr*rich_cr
         do i=1,nlev-1 
            if (tke(i).le.klimiw) then
               rich(i)=NN(i)/(SS(i)+1.e-10)
               if (rich(i).lt.rich_cr) then
                  if (rich(i).gt.0) then
                     pot=1-rich(i)*rich(i)/rich2 
                     x=numshear*pot*pot*pot
                     num(i)=numiw+x 
                     nuh(i)=nuhiw+x  
                  else
                     num(i)=numiw+numshear
                     nuh(i)=nuhiw+numshear
                  end if          
               else
                  num(i)=numiw
                  nuh(i)=nuhiw
               end if
            end if   
         end do
      case default
   end select

   !include 'GOTMVariables_out.F90'

   end subroutine internal_wave 
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: The Kolmogorov-Prandtl relation. 
!
! !INTERFACE:
   subroutine kolpran(ObjGOTM,nlev,u_taub,u_taus,z0b,z0s)
!
! !DESCRIPTION:
!  Eddy viscosity/diffusivity are calculated by means of the relation of 
!  Kolmogorov and Prandtl from the computed values of k, L and 
!  stability functions. 
!  \begin{equation}
!  \nu_t = c_{\mu} \sqrt{k}L;\quad \nu'_t = c'_{\mu} \sqrt{k}L,
!  \end{equation}
!
! !USES:
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: u_taub,u_taus
   double precision, intent(in) :: z0b,z0s
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  See turbulence module
!
! !LOCAL VARIABLES:
   double precision     :: x
   integer      :: i

   include 'GOTMVariables_in.F90'
!
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev
      x=sqrt(tke(i))*L(i)
      num(i)=cmue1(i)*x
      nuh(i)=cmue2(i)*x
   end do

   num(0)=kappa*u_taub*z0b
   num(nlev)=kappa*u_taus*z0s 
   nuh(0)=kappa*u_taub*z0b
   nuh(nlev)=kappa*u_taus*z0s 

   !include 'GOTMVariables_out.F90'
 
   end subroutine kolpran 
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: An equation for the length scale. 
! 
! !INTERFACE:
   subroutine lengthscaleeq(ObjGOTM,nlev,dt,z0b,z0s,depth,h,NN,P,B,numtke) 
!
! !DESCRIPTION:
!  This subroutine calculates the lengthscale equation according to
!  Mellor and Yamada [1982]:
!
!  \begin{equation}\label{kL_eq}
!  \partial_t (kL) - \partial_z\left(\nu_L\partial_z (kL)\right) =
!  L \left(c_{L1}P + c_{L3}B
!  -  \left(1 + E_2\left(\frac{L}{L_z}\right)^2\right)\varepsilon \right).
!  \end{equation}
!
!  At the end of the subroutine, the Galperin et al. [1988] length
!  limitation is applied and the disispation rate calculated. 
!
! !USES:
!   use mtridiagonal
!   use turbulence, ONLY: kappa,cde,galp,length_lim
!   use turbulence, ONLY: MY_length,e1,e2,e3,b1
!   use turbulence, ONLY: num,L,L_min,eps,eps_min,tkeo,tke
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: dt
   double precision, intent(in) :: z0b,z0s
   double precision, intent(in) :: depth,h(0:nlev)
   double precision, intent(in) :: NN(0:nlev)
   double precision, intent(in) :: P(0:nlev),B(0:nlev)
   double precision, intent(in) :: numtke(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   double precision         :: avh(0:nlev),q2l(0:nlev),q3(0:nlev)
   double precision         :: phi_minus(0:nlev),phi_plus(0:nlev)
   double precision         :: Lz(0:nlev)
   double precision             :: ds,db,prod,buoyan,diss,Lcrit 
   integer      :: i 

   include 'GOTMVariables_in.F90'

! 
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      q2l(i)=2.*tkeo(i)*L(i)
      q3 (i)=sqrt(8.*tke(i)*tke(i)*tke(i))
   end do
 
   do i=1,nlev
      avh(i)=0.5*(numtke(i-1)+numtke(i))
   end do
     
   db=0.0    ! Diagnostic Length Scale  
   ds=0.0 
   do i=1,nlev-1   
      db=db+h(i) 
      ds=depth-db
! Parabola shape
      if (MY_length.eq.1) Lz(i)=kappa*(ds+z0s)*(db+z0b)/(ds+z0s+db+z0b)
! Triangle shape
      if (MY_length.eq.2) Lz(i)=kappa*min(ds+z0s,db+z0b)
! For infinite depth
      if (MY_length.eq.3) Lz(i)=kappa*(ds+z0s)
   end do
       
   do i=1,nlev-1
      prod=e1*L(i)*P(i)
      buoyan=e3*L(i)*B(i)
      diss=-q3(i)/b1*(1.+e2*(L(i)/Lz(i))*(L(i)/Lz(i)))
      if (prod+buoyan .gt. 0) then
         phi_plus(i)=prod+buoyan
         phi_minus(i)=-diss
      else
         phi_plus(i)=prod
         phi_minus(i)=-buoyan-diss
      end if
   end do
 
   do i=2,nlev-2
      au(i)=-2.*dt*avh(i)/(h(i)+h(i+1))/h(i)
      cu(i)=-2.*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
      bu(i)=1.-au(i)-cu(i)+phi_minus(i)*dt/q2l(i)
      du(i)=(1+phi_plus(i)*dt/q2l(i))*q2l(i)
   end do
 
   cu(1)=0.
   bu(1)=1.
   du(1)=2.*tke(1)*kappa*(z0b+h(1)) 

   bu(nlev-1)=1.
   au(nlev-1)=0.
   du(nlev-1)=2.*tke(nlev-1)*kappa*(z0s+h(nlev)) 

   call Tridiagonal(ObjGOTM,nlev,1,nlev-1,q2l)
 
   do i=1,nlev-1
     L(i)=q2l(i)/(2.*tke(i))
     if ((NN(i).gt.0).and.(length_lim)) then
       Lcrit=sqrt(2*galp*galp*tke(i)/NN(i)) 
       if (L(i).gt.Lcrit) L(i)=Lcrit  
     end if
     if (L(i).lt.L_min) L(i)=L_min
     eps(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/L(i)
     if (eps(i).lt.eps_min) eps(i)=eps_min
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine lengthscaleeq
!EOC

!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: PotentialML - length sacle using 2 master lenghth scales
!
! !INTERFACE:
   subroutine potentialml(ObjGOTM,nlev,z0b,z0s,h,depth,NN)

! !DESCRIPTION:
!  Computes the length scale by defining two master 
!  length scales $l_u$ and $l_d$
!  \begin{equation}
!  \begin{array}{l}
!  \int_{z_0}^{z_0+l_u(z_0)} [b(z_0)-b(z)] dz =k(z_0) \\
!  \int_{z_0-l_d(z_0)}^{z_0} [b(z)-b(z_0)] dz =k(z_0)
!  \end{array}
!  \end{equation}
! 
!   From $l_u$ and $l_d$ two length scales are defined $l_k$ 
!   (characteristic mixing length)
!   and $l_\epsilon$ (characteristic dissipation length):
!   \begin{equation}
!   \begin{array}{l}
!   l_k(z_0)= min[l_d(z_0),l_u(z_0)] \\
!   l_{\epsilon}(z_0)={[l_d(z_0)l_u(z_0)]}^{1/2}
!   \end{array}
!   \end{equation}
! 
!   $l_k$ is used in kolpran() to compute eddy viscosity/difussivity  
!   (is transported as L()). $l_{\epsilon}$ is ed to compute $\epsilon$:
!   \begin{equation}
!   \epsilon=C_{\epsilon}k^{3/2}l_{\epsilon}^{-1}, with C_{\epsilon}=0.7
!   \end{equation}
!
! !USES:
!   use turbulence, ONLY: L,eps,tkeo,tke,L_min,eps_min,cde,galp,kappa,length_lim
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev 
   double precision, intent(in) :: h(0:nlev),depth,NN(0:nlev)
   double precision, intent(in) :: z0b,z0s
!
! !OUTPUT PARAMETERS:
!
! !BUGS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM Team
!
!  $Log$
!
! !LOCAL VARIABLES:
   double precision     :: ds(0:nlev),db(0:nlev)
   double precision     :: lu(0:nlev),ld(0:nlev)
   double precision     :: lk(0:nlev),leps(0:nlev)
   double precision, parameter  :: NNmin=1.e-8
   double precision     :: Lcrit,buoydiff,integral,ceps
   integer      :: i,j

   include 'GOTMVariables_in.F90'
!EOP
!-------------------------------------------------------------------------
!BOC
   db(0)=0.
   ds(nlev)=0.
   do i=1,nlev-1
      db(i)=db(i-1)+h(i)      ! distance of intercace i from bottom 
      ds(i)=depth-db(i)       ! distance of intercace i from surface 
   end do
!
!  Calculation of lu and ld by solving the integral equation following 
!  Gaspar (1990). Some other approximations of the integral equation 
!  are possible.
!
! Computation of lupward
!
   do i=1,nlev-1
      lu(i)=0.
      integral=0.
      buoydiff=0.
      do j=i+1,nlev
         buoydiff=buoydiff+NN(j-1)*0.5*(h(j)+h(j-1))
         integral=integral+buoydiff*h(j)
         if (integral.ge.tkeo(i)) then
            if(j.ne.nlev) then
               if(j.ne.i+1) then
                  lu(i)=lu(i)-(integral-tkeo(i))/buoydiff
               else 
!           To avoid lu(i) from becoming too large if NN(i) is too small
              if(NN(i).gt.NNmin) then
                 lu(i)=sqrt(2.)*sqrt(tkeo(i))/sqrt(NN(i))
                  else
                 lu(i)=h(i)
                  end if
               end if 
               goto 600
            end if
         end if
         lu(i)=lu(i)+h(j)
      end do 
600   continue
!     Implicitely done in the do loop: if (lu(i).gt.ds(i)) lu(i)=ds(i) 
!     lu limited by distance to surface 
   end do

!  Computation of ldownward
   do i=nlev-1,1,-1
      ld(i)=0.
      integral=0.
      buoydiff=0.
      do j=i-1,1,-1 
         buoydiff=buoydiff+NN(j)*0.5*(h(j+1)+h(j))
         integral=integral-buoydiff*h(j)
         if (integral.ge.tkeo(i)) then
            if(j.ne.0) then
               if(j.ne.i-1) then
                  ld(i)=ld(i)-(integral-tkeo(i))/buoydiff
               else
!              To avoid ld(i) from becoming too large if NN(i) is too small
                  if(NN(i).gt.NNmin) then
                     ld(i)=sqrt(2.)*sqrt(tkeo(i))/sqrt(NN(i))
                  else
                     ld(i)=h(i)
                  end if
               end if 
               goto 610
            end if
         end if
         ld(i)=ld(i)+h(j)
      end do
610   continue
!     if (ld(i).gt.db(i)) ld(i)=db(i) !ld limited by distance to bottom
   end do         

!   Calculation of lk and leps, mixing and dissipation lengths
   do i=nlev-1,1,-1 
!  Suggested by Gaspar:        lk(i)   = min(lu(i),ld(i))
      lk(i)=sqrt(lu(i)*ld(i))
      leps(i) = sqrt(lu(i)*ld(i)) 
   end do

!  We set L=lk because it is the one we use to calculate num and nuh
   ceps=0.7
   do i=1,nlev-1
      L(i)=lk(i)
   end do      
     
   L(0)=kappa*z0b
   L(nlev)=kappa*z0s
!  Gaspar uses null gradient
   do i=0,nlev
      if ((NN(i).gt.0).and.(length_lim)) then
         Lcrit=sqrt(2*galp*galp*tke(i)/NN(i))
         if (L(i).gt.Lcrit) L(i)=Lcrit
      end if
      if (L(i).lt.L_min) L(i)=L_min
      eps(i)=cde*sqrt(tke(i)*tke(i)*tke(i))/L(i)
      if(eps(i).lt.eps_min) eps(i)=eps_min
   end do  

   !include 'GOTMVariables_out.F90'

   end subroutine potentialml
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Algebraic equation for TKE (local equalibrium).
! 
! !INTERFACE:
   subroutine tkealgebraic(ObjGOTM,nlev,u_taus,u_taub,NN,SS)
!
! !DESCRIPTION:
!  This subroutine calculates the turbulent kinetic energy based
!  on the local equilibrium assumption
!
!  \begin{equation}
!  P+B-\varepsilon=0.
!  \end{equation}  
!
! !USES:
!   use turbulence, only: tkeo,tke,L,k_min,cmue2,cde,cmue1,cm0
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: u_taus,u_taub
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      tkeo(i)=tke(i)
      tke(i)=L(i)*L(i)/cde*(cmue1(i)*SS(i)-cmue2(i)*NN(i))
   end do 
   tke(0)=u_taub*u_taub/sqrt(cm0*cde) 
   tke(nlev)=u_taus*u_taus/sqrt(cm0*cde)
   do i=0,nlev
      if (tke(i).lt.k_min) tke(i)=k_min 
   end do 
 
   !include 'GOTMVariables_out.F90'
  
   end subroutine tkealgebraic
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Dynamic equation for TKE.
! 
! !INTERFACE:
   subroutine tkeeq(ObjGOTM,N,dt,u_taus,u_taub,z0s,h,P,B,numtke)
!
! !DESCRIPTION:
!  This subroutine calculates the turbulent kinetic energy as
!  needed for one- or two-equation models:
!
!  \begin{equation}\label{k_eq}
!  \partial_t k - \partial_z(\nu_k\partial_z k) =  P + B -\varepsilon,
!  \end{equation}
!
!  The diffusion coefficient depends on the type of model (k-epsilon
!  or Mellor-Yamada). 
!
!  As boundary conditions a choice between Dirichlet (flux-bdy=.false.)
!  and Neumann no-flux conditions (flux-bdy=.true.) has to be made.
!
!  Dirichlet condition:
!
!  \begin{equation}
!  k=\left(\frac{u_*}{c_{\mu}^0}\right)^2. 
!  \end{equation}
!
!  If flux conditions are chose, the {\it Craig and Banner} [1994] and
!  the {\it Craig} [1996] surface wave breaking theory can
!  be used. The boundarz condition is then:
!
!   \begin{equation}
!   -\nu_t \partial_zk =c_w (u_s^*)^3 \hfill z= 0. \qquad
!   \end{equation}
!
!   Since the flux is applied half a grid-box away from the boundary,
!   the {\it Craig} [1996] analytical solution is used for the
!   construction of this boundary condition:
!
!  \begin{equation}\label{tkeanalyt}
!  k=\frac{(u_s^*)^2}{c_{\mu}^{1/2}}
!  \left[a+
!  \left(\frac{3\sigma_k}{2}\right)^{1/2}
!  c_{\mu}^{1/4}c_w\left(\frac{z'+z_0^s}{z_0^s}\right)^{-m}
!  \right]^{2/3}. 
!  \end{equation}
! 
!   
!  The sink terms are treated quasi-implicitely in oder to guarantee
!  positivity. 
!
! !USES:
!   use mTridiagonal
!   use turbulence, ONLY: tkeo,tke,k_min,eps,num,cm0,cde,flux_bdy
!   use turbulence, ONLY: tke_method
!   use turbulence, ONLY: craig_banner,craig_m,cw
!   IMPLICIT NONE
!
! !INPUT PARAMETERS: 
   integer, intent(in)  :: N
   double precision, intent(in) :: dt 
   double precision, intent(in) :: u_taus,u_taub,z0s
   double precision, intent(in) :: h(0:N)
   double precision, intent(in) :: P(0:N),B(0:N)
   double precision, intent(in) :: numtke(0:N)
!
! !INPUT/OUTPUT PARAMETERS: 
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
! 
! !LOCAL VARIABLES:
   double precision         :: avh(0:N)
   double precision         :: pminus(0:N),pplus(0:N)
   double precision         :: prod,buoyan,diss
   double precision         :: zz 
   logical, save    :: warning=.true. 
   integer      :: i

   include 'GOTMVariables_in.F90'

!EOP
!-----------------------------------------------------------------------
!BOC


   tkeo=tke
 
   do i=2,N-1
      avh(i)=0.5*(numtke(i-1)+numtke(i))
   end do

   if (flux_bdy) then 
      avh(1)=0.
      avh(N)=0.
   else       ! Not needed any more for Dirichlet conditions, only 
! kept for "historical" reasons, see Burchard et al. [1998].
      if (tke_method .eq. 2) then  
         avh(1)=u_taub**4*2/(eps(0)+eps(1  ))
         avh(N)=u_taus**4*2/(eps(N)+eps(N-1))
      else   
         avh(1)= 0.5*(num(0)+numtke(1  ))
         avh(N)= 0.5*(num(N)+numtke(N-1))
      end if
   end if

   zz=0.
   do i=N-1,1,-1
      prod=P(i)
      buoyan=B(i)
      diss=eps(i)
      if (prod+buoyan.gt.0) then
         pplus(i)=prod+buoyan
         pminus(i)=diss
      else
         pplus(i)=prod
         pminus(i)=diss-buoyan
      end if
   end do
   i=N-1
 
   do i=1,N-1
      au(i)=-2.*dt*avh(i)/(h(i)+h(i+1))/h(i)
      cu(i)=-2.*dt*avh(i+1)/(h(i)+h(i+1))/h(i+1)
      bu(i)=1.-au(i)-cu(i)+pminus(i)*dt/tke(i)
      du(i)=(1+pplus(i)*dt/tke(i))*tke(i)
   end do

!  Surface flux of TKE due to surface wave breaking 
!  according to Craig and Banner 1994: 

   if (craig_banner) then 
      if (h(N) .gt. z0s .and.  warning) then
         write(0,*) 'WARNING: Surface roughness length smaller than'
         write(0,*) '         thickness of upper grid box. Calculations'
         write(0,*) '         might be inaccurate in this Craig and Banner'
         write(0,*) '         surface wave breaking parameterisation.'
         write(0,*) '         Computation is continued.' 
         warning=.false.   
      end if 
      du(N-1)=du(N-1)+cw*u_taus**3*dt/(0.5*(h(N)+h(N-1)))   &
                 *((0.5*h(N)+z0s)/z0s)**(-craig_m)                      
   end if 

   if (flux_bdy) then
!  +-------------------------------------------------------------+
!  | No-flux conditions for TKE                                  | 
!  +-------------------------------------------------------------+
      call tridiagonal(ObjGOTM,N,1,N-1,tke)
      tke(0) = u_taub*u_taub/sqrt(cm0*cde) 
      tke(N) = u_taus*u_taus/sqrt(cm0*cde) 
   else
!  +-------------------------------------------------------------+
!  | Dirichlet conditions for TKE                                | 
!  +-------------------------------------------------------------+
      if (craig_banner) then
          write(0,*) 'For the Craig and Banner wave breaking condition,' 
          write(0,*) 'flux boundary conditions should be used.'
          write(0,*) 'Please, change namelist gotmturb.inp !'
          write(0,*) 'Program aborted ...'
          stop
       end if   
      cu(1)=0.
      bu(1)=1.
      du(1)=u_taub*u_taub/sqrt(cm0*cde)
 
      bu(N-1)=1.
      au(N-1)=0.
      du(N-1)=u_taus*u_taus/sqrt(cm0*cde)

      call tridiagonal(ObjGOTM,N,1,N-1,tke)
   end if 

   where (tke .lt. k_min) tke = k_min

   !include 'GOTMVariables_out.F90'
   
   end subroutine tkeeq
!EOC




!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Canuto et al. [2000] version A non-equilibrium stability functions. 
! 
! !INTERFACE:
   subroutine cmue_ca(ObjGOTM,nlev)
!
! !DESCRIPTION:
!  This subroutine computes Canuto et al. [2000] version A non-equilibrium
!  stability functions.
!
! !USES:
!   use turbulence, only: an,as
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision, parameter  :: L1=0.1070
   double precision, parameter  :: L2=0.0032
   double precision, parameter  :: L3=0.0864
   double precision, parameter  :: L4=0.1200
   double precision, parameter  :: L5=11.9000
   double precision, parameter  :: L6=0.4000
   double precision, parameter  :: L7=0.0000
   double precision, parameter  :: L8=0.4800
   double precision, parameter  :: cm0_sf=0.5270
   double precision, parameter  :: tnmin=-12.27
   double precision, parameter  :: a2_cm03=2./cm0_sf**3
   double precision     :: s0,s1,s2,s4,s5,s6
   double precision     :: d,d0,d1,d2,d3,d4,d5
   double precision     :: tsmax
   double precision     :: tn,ts,sm,sh

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   s0 = 1.5*L1*L5*L5
   s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
   s2 = -0.375*L1*(L6*L6-L7*L7)
   s4 = 2.*L5
   s5 = 2.*L4
   s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7)

   d0 = 3.*L5*L5
   d1 = L5*(7.*L4+3.*L8)
   d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
   d3 = L4*(4.*L4+3.*L8)
   d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
   d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

   do i=1,nlev-1
      tn = 4./cm0_sf**6 * an(i)
      if (tn .lt. tnmin) tn = tnmin

      ts = 4./cm0_sf**6 * as(i)
      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)
      if (ts.gt.tsmax) ts = tsmax

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1(i) = a2_cm03 * sm
      cmue2(i) = a2_cm03 * sh
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_ca
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Canuto et al. [2000] version A quasi-equilibrium stability functions. 
!  
! !INTERFACE:
   subroutine cmue_caqe(ObjGOTM,nlev)
!
! !DESCRIPTION:
!  This subroutine computes Canuto et al. [2000] version A quasi-equilibrium
!  stability functions.
!
! !USES:
!   use turbulence, only: an
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team Q
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision, parameter  :: L1=0.1070
   double precision, parameter  :: L2=0.0032
   double precision, parameter  :: L3=0.0864
   double precision, parameter  :: L4=0.1200
   double precision, parameter  :: L5=11.9000
   double precision, parameter  :: L6=0.4000
   double precision, parameter  :: L7=0.0000
   double precision, parameter  :: L8=0.4800
   double precision, parameter  :: cm0_sf=0.5270
   double precision, parameter  :: tnmin=-12.27
   double precision, parameter  :: a2_cm03=2./cm0_sf**3
   double precision     :: s0,s1,s2,s4,s5,s6
   double precision     :: d,d0,d1,d2,d3,d4,d5
   double precision     :: tsmax
   double precision     :: tn,ts,sm,sh
   double precision     :: PP,QQ

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   s0 = 1.5*L1*L5*L5
   s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
   s2 = -0.375*L1*(L6*L6-L7*L7)
   s4 = 2.*L5
   s5 = 2.*L4
   s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7)

   d0 = 3.*L5*L5
   d1 = L5*(7.*L4+3.*L8)
   d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
   d3 = L4*(4.*L4+3.*L8)
   d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
   d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

   do i=1,nlev-1
      tn = 4./cm0_sf**6 * an(i)
      if (tn .lt. tnmin) tn = tnmin

      PP=(s0+(s1-s6)*tn-2.*(d2+d4*tn))/(s2-2.*d5)
      QQ=-(2.*(d0+d1*tn+d3*tn*tn)+(s4+s5*tn)*tn)/(s2-2.*d5)

      ts=-0.5*PP-sqrt(PP**2/4.-QQ)

      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1(i) = a2_cm03 * sm
      cmue2(i) = a2_cm03 * sh
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_caqe
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Canuto et al. [2000] version B non-equilibrium stability functions. 
! 
! !INTERFACE:
   subroutine cmue_cb(ObjGOTM,nlev)
!
! !DESCRIPTION:
!  This subroutine computes Canuto et al. [2000] version B non-equilibrium
!  stability functions.
!
! !USES:
!   use turbulence, only: an,as
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision, parameter  :: L1=0.1270
   double precision, parameter  :: L2=0.00336
   double precision, parameter  :: L3=0.0906
   double precision, parameter  :: L4=0.1010
   double precision, parameter  :: L5=11.2000
   double precision, parameter  :: L6=0.4000
   double precision, parameter  :: L7=0.0000
   double precision, parameter  :: L8=0.3180
   double precision, parameter  :: cm0_sf=0.5540
   double precision, parameter  :: tnmin=-12.27
   double precision, parameter  :: a2_cm03=2./cm0_sf**3
   double precision     :: s0,s1,s2,s4,s5,s6
   double precision     :: d,d0,d1,d2,d3,d4,d5
   double precision     :: tsmax
   double precision     :: tn,ts,sm,sh

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   s0 = 1.5*L1*L5*L5
   s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
   s2 = -0.375*L1*(L6*L6-L7*L7)
   s4 = 2.*L5
   s5 = 2.*L4
   s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7)

   d0 = 3.*L5*L5
   d1 = L5*(7.*L4+3.*L8)
   d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
   d3 = L4*(4.*L4+3.*L8)
   d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
   d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

   do i=1,nlev-1
      tn = 4./cm0_sf**6 * an(i)
      if (tn .lt. tnmin) tn = tnmin

      ts = 4./cm0_sf**6 * as(i)
      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)
      if (ts.gt.tsmax) ts = tsmax

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1(i) = a2_cm03 * sm
      cmue2(i) = a2_cm03 * sh
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_cb
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Canuto et al. [2000] version B quasi-equilibrium stability functions. 
! 
! !INTERFACE:
   subroutine cmue_cbqe(ObjGOTM,nlev)
!
! !DESCRIPTION:
!  This subroutine computes Canuto et al. [2000] version B quasi-equilibrium
!  stability functions.
!
! !USES:
!   use turbulence, only: an
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision, parameter  :: L1=0.1270
   double precision, parameter  :: L2=0.00336
   double precision, parameter  :: L3=0.0906
   double precision, parameter  :: L4=0.1010
   double precision, parameter  :: L5=11.2000
   double precision, parameter  :: L6=0.4000
   double precision, parameter  :: L7=0.0000
   double precision, parameter  :: L8=0.3180
   double precision, parameter  :: cm0_sf=0.5540
   double precision, parameter  :: tnmin=-12.27
   double precision, parameter  :: a2_cm03=2./cm0_sf**3
   double precision     :: s0,s1,s2,s4,s5,s6
   double precision     :: d,d0,d1,d2,d3,d4,d5
!kbk   double precision     :: tsmax
   double precision     :: tn,ts,sm,sh
   double precision     :: PP,QQ

   include 'GOTMVariables_in.F90'

!EOP
!-----------------------------------------------------------------------
!BOC
   s0 = 1.5*L1*L5*L5
   s1 = -L4*(L6+L7)+2.*L4*L5*(L1-1./3.*L2-L3)+1.5*L1*L5*L8
   s2 = -0.375*L1*(L6*L6-L7*L7)
   s4 = 2.*L5
   s5 = 2.*L4
   s6 = 2./3.*L5*(3.*L3*L3-L2*L2)-0.5*L5*L1*(3.*L3-L2)+0.75*L1*(L6-L7)

   d0 = 3.*L5*L5
   d1 = L5*(7.*L4+3.*L8)
   d2 = L5*L5*(3.*L3*L3-L2*L2)-0.75*(L6*L6-L7*L7)
   d3 = L4*(4.*L4+3.*L8)
   d4 = L4*(L2*L6-3.*L3*L7-L5*(L2*L2-L3*L3))+L5*L8*(3.*L3*L3-L2*L2)
   d5 = 0.25*(L2*L2-3*L3*L3)*(L6*L6-L7*L7)

   do i=1,nlev-1
      tn = 4./cm0_sf**6 * an(i)
      if (tn .lt. tnmin) tn = tnmin

      PP=(s0+(s1-s6)*tn-2.*(d2+d4*tn))/(s2-2.*d5)
      QQ=-(2.*(d0+d1*tn+d3*tn*tn)+(s4+s5*tn)*tn)/(s2-2.*d5)

      ts=-0.5*PP-sqrt(PP**2/4.-QQ)

!kbk      tsmax = (d0+d1*tn+d3*tn*tn)/(d2+d4*tn)

      d = d0 + d1*tn + d2*ts + d3*tn*tn + d4*tn*ts + d5*ts*ts

      sm = (s0 + s1*tn + s2*ts) / d
      sh = (s4 + s5*tn + s6*ts) / d

      cmue1(i) = a2_cm03 * sm
      cmue2(i) = a2_cm03 * sh
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_cbqe
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Kantha and Clayson [1994] non-equilibrium stability functions.
! 
! !INTERFACE:
   subroutine cmue_kc(ObjGOTM,nlev)
!
! !DESCRIPTION:
!  This subroutine computes Burchard and Baumert [1995] stability functions.
!
! !USES:
!   use turbulence, only: a1,a2,b2,c1,c2,c3
!   use turbulence, only: as,an
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision     :: gm,gh,sm,sh
   double precision     :: c11,c12,c13,c21,c22,c23

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      gm=0.5*as(i)             !Transformation to MY notation
      gh=-0.5*an(i)            !Transformation to MY notation
      if (gh.gt.0.029) gh=0.029

      if (gm.gt.0.825-25.0*gh) gm=0.825-25.0*gh

      c11=6*a1*a2*gm
      c12=1-3*a2*b2*(1-c3)*gh-12*a1*a2*gh
      c13=a2
      c21=1+6*a1*a1*gm-9*a1*a2*gh
      c22=-12*a1*a1*gh-9*a1*a2*(1-c2)*gh
      c23=a1*(1-3*c1)
      sm=(c12*c23-c22*c13)/(c12*c21-c22*c11)
      sh=(c21*c13-c11*c23)/(c12*c21-c22*c11)
      cmue1(i)=sqrt(2.)*sm     !Retransformation to GOTM notation
      cmue2(i)=sqrt(2.)*sh     !Retransformation to GOTM notation
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_kc
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Kantha and Clayson [1995] stability functions.
! 
! !INTERFACE:
   subroutine cmue_kcqe(ObjGOTM,nlev)
!
! !DESCRIPTION:
!  This subroutine computes Kantha and Clayson [19995] stability functions
!  including smoothing for convective conditions as discussed by Burchard
!  and Petersen [1999].
!
! !USES:
!   use turbulence, only: a1,a2,b1,b2,c2,c3,an
!   use turbulence, only: qesmooth,qeghmax,qeghmin,qeghcrit
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision         :: gh,sm,sh

   include 'GOTMVariables_in.F90'
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      gh=-0.5*an(i)           ! Transformation to MY notation 
      qeghmax=1/(a2*(b1+12*a1+3*b2*(1-c3)))
      if (qesmooth) then
         if (gh.gt.qeghcrit) gh=gh-(gh-qeghcrit)**2/(gh+qeghmax-2*qeghcrit)
      else 
         if (gh.gt.qeghmax) gh=qeghmax  
      end if     
      if (gh.lt.qeghmin)  gh = qeghmin  
      sh=a2*(1-6*a1/b1)/(1-3*a2*gh*(6*a1+b2*(1-c3)))
      sm=(b1**(-1./3.)+9*a1*(2*a1+a2*(1-c2))*sh*gh)/(1-9*a1*a2*gh)

      cmue1(i)=sqrt(2.)*sm    ! Retransformation to GOTM notation 
      cmue2(i)=sqrt(2.)*sh    ! Retransformation to GOTM notation
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_kcqe
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Munk-Anderson stability functions. 
! 
! !INTERFACE:
   subroutine cmue_ma(ObjGOTM,nlev,NN,SS)
!
! !DESCRIPTION: 
!  This subroutine computes stability functions 
!  according to Munk \& Anderson [1948]:  
!
!  \begin{equation}
!  \begin{array}{ll}
!  c_{\mu} = cm0 \\
!  c_{\mu}'= \frac{c_{\mu}}{P_r^0}
!  \frac{(1+10 R_i)^{1/2}}{(1+3.33 R_i)^{3/2}}, &  R_i \geq 0,\\
!  c_{\mu}'= c_{\mu}, &  R_i<0,
!  \end{array}
!  \end{equation}
!
! !USES:
!   use turbulence, only: cm0,Prandtl0
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision     :: Ri,Prandtl

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      Ri=NN(i)/(SS(i)+1e-8)   ! Gradient Richardson number 
      if (Ri.ge.1e-10) then 
         Prandtl=Prandtl0*(1.+3.33*Ri)**1.5/sqrt(1.+10.0*Ri)
      else
         Prandtl=Prandtl0
      end if
      cmue1(i)=cm0
      cmue2(i)=cm0/Prandtl
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Simple formula stability functions.
! 
! !INTERFACE:
   subroutine cmue_rf(ObjGOTM,nlev,NN,SS)
!
! !DESCRIPTION: 
!  This subroutine computes stability functions with some "simple" formulas
!
!  Stab=Isprastab
!  \begin{equation}
!  \begin{array}{ll}
!  c_{\mu} = cm0 \\
!  c_{\mu}'= \frac {1} {Prandtl0} (1-R_f)^{1/2}
!  \end{array}
!  \end{equation}
!  \begin{equation} 
!  1-R_f=(\sqrt{R_i^2+1}-R_i)^2 
!  \end{equation}
!  \begin{equation}
!  R_i=\frac {1} {2 Prandtl0} \frac {N^2} {S^2}
!  \end{equation}
!
!  $R_f$ is limited for supercritically stable stratification $1.8<R_f$.  
!
! !USES:
!   use turbulence, only: cm0,Prandtl0,xRF 
!   use turbulence, only: cmue1,cmue2 
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision     :: Ri,Prandtl_inv

   include 'GOTMVariables_in.F90'
! 
!EOP
!-----------------------------------------------------------------------
!BOC
!  Calculation of xRf=(1-Rf), where Rf is the flux Richardson number
   do i=1,nlev-1
      Ri=0.5/Prandtl0*NN(i)/(SS(i)+1e-8)
      xRf(i)=(sqrt(Ri*Ri+1)-Ri)**2 
      if (xRf(i) .gt. 2.) xRf(i)=2.
      Prandtl_inv=1/Prandtl0*sqrt(xRf(i))

      if (Prandtl_inv.lt.0.18) Prandtl_inv=0.18
      if (Prandtl_inv.gt.2.0)  Prandtl_inv=2.0
   
      cmue1(i)=cm0
      cmue2(i)=cm0*Prandtl_inv

   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_rf
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Schumann and Gerz [1995] parameterization of the Prandtl number. 
! 
! !INTERFACE:
   subroutine cmue_sg(ObjGOTM,nlev,NN,SS)
!
! !DESCRIPTION:
!  This subroutine computes Schumann and Gerz [1995] stability functions.
!
! \begin{equation}
! c_{\mu}=c_{\mu}^0,\qquad c'_{\mu}=\frac{c_{\mu}^0}{P_r}
! \end{equation}
!
! with constant $c_{\mu}^0$. We choose here for the Prandtl number $P_r$ a
! formulation suggested by {\it Schumann and Gerz} [1995]:
!
! \begin{equation}
! P_r=P_r^0\exp\left(-\frac{R_i}{P_r^0R_i^{\infty}}\right)
! -\frac{R_i}{R_i^{\infty}}
! \end{equation}
!
! with the neutral Prandtl number $P_r^0=0.74$ and $R_i^{\infty}=0.25$.
!
! !USES:
!   use turbulence, only: Prandtl0,cm0
!   use turbulence, only: cmue1,cmue2
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: nlev
   double precision, intent(in) :: NN(0:nlev),SS(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY: 
!  Original author(s): GOTM Team 
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   double precision     :: Ri,Prandtl

   include 'GOTMVariables_in.F90'
!EOP
!-----------------------------------------------------------------------
!BOC
   do i=1,nlev-1
      Ri=NN(i)/(SS(i)+1e-8)   ! Gradient Richardson number
      if (Ri.ge.1e-10) then
         Prandtl=Prandtl0*exp(-Ri/(Prandtl0*0.25))+Ri/0.25
      else
         Prandtl=Prandtl0
      end if

      cmue1(i)=cm0
      cmue2(i)=cm0/Prandtl

   end do

   !include 'GOTMVariables_out.F90'

   end subroutine cmue_sg
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: tridiagonal. For consistency, it should be include in ModuleTriangulation
!
! !INTERFACE:
   subroutine tridiagonal(ObjGOTM,N,fi,lt,value)
!
! !DESCRIPTION:
!
! A linear equation with tridiagonal matrix is solved here. The main
! diagonal is stored on bu, the upper diagonal on au, and the
! lower diagonal on cu, the right hand side is stored on du. The method
! used here is the simplified Gauss elimination, also called Thomas algorithm.  
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)  :: N,fi,lt
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
   double precision     :: value(0:N)
!
! !REVISION HISTORY:
!  Original author(s): GOTM Team
!
!  $Log$
!
! !LOCAL VARIABLES:
   integer      :: i
   
   include 'GOTMVariables_in.F90'

!
!EOP
!-----------------------------------------------------------------------
!BOC
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)
   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do

   !include 'GOTMVariables_out.F90'

   end subroutine tridiagonal

!EOC

end module ModuleGOTM 

!----------------------------------------------------------------------------------------------------------
!Copyright (C) 2000 - GOTM code.
!----------------------------------------------------------------------------------------------------------

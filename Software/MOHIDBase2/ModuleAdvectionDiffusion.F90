!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : ModuleAdvectionDiffusion 
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : Jun 2003
! REVISION      : Paulo Leitão - v4.0
! DESCRIPTION   : Module responsbile for computing advection/diffusion processes
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

Module ModuleAdvectionDiffusion
                                    
    use ModuleGlobalData
    use ModuleFunctions     , only : OrlanskiCelerity2D, THOMAS_3D, THOMASZ, ComputeAdvection1D_V2,     &
                                     SetMatrixValue, Chunk_J, Chunk_K, Chunk_I, T_THOMAS, T_VECGW, T_D_E_F
    use ModuleTime          , only : GetComputeCurrentTime, T_Time, KillComputeTime, &
                                     null_time, operator(+), operator(-),            &
                                     operator (==), operator (/=)
    use ModuleStopWatch     , only : StartWatch, StopWatch     
    use ModuleHorizontalGrid, only : KillHorizontalGrid, GetHorizontalGrid, UnGetHorizontalGrid
    use ModuleHorizontalMap , only : KillHorizontalMap, GetBoundaries, UnGetHorizontalMap
    use ModuleGeometry      , only : GetGeometrySize, GetGeometryKFloor, GetGeometryWaterColumn, &
                                     GetGeometryAreas, GetGeometryDistances, UnGetGeometry
    !$ use omp_lib

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartAdvectionDiffusion
    private ::      AllocateInstance
    private ::      AllocateVariables                                                


    !Selector
    public  :: GetAdvFlux
    public  :: GetDifFlux
    public  :: GetBoundaryConditionList

    public  :: UngetAdvectionDiffusion

    public  :: SetDischarges
    public  :: UnSetDischarges

    !Modifier
    public  :: AdvectionDiffusion                                                
    private ::      Set_Internal_State          
    private ::      Convert_Dif_Vertical

    private ::      AdvectionDiffusionIteration

    private ::          Convert_Visc_Dif_Horizontal                                      

    private ::          VolumeVariation                     !OK

    private ::          Discharges

    private ::          HorizontalDiffusion                 !OK                                         
    private ::              HorizontalDiffusionXX           !OK                             
    private ::                  CalcHorizontalDifFluxXX     !OK                                     
    private ::              HorizontalDiffusionYY           !OK
    private ::                  CalcHorizontalDifFluxYY     !OK                                     
                                      
    private ::          OpenBoundaryCondition                                                     
    private ::              FluxAtOpenBoundary                                               

    private ::          HorizontalAdvection                                             
    private ::              HorizontalAdvectionXX                                      
    private ::                  CalcHorizontalAdvFluxXX                                         
    private ::              HorizontalAdvectionYY                                      
    private ::                  CalcHorizontalAdvFluxYY                                         

    private ::          VerticalDiffusion                                                 
    private ::              CalcVerticalDifFlux                                              
    private ::          VerticalAdvection      
    private ::              CalcVerticalAdvFlux                                             

    private ::      FinishAdvectionDiffusionIt                            

    !Destructor
    public  ::  KillAdvectionDiffusion
    private ::      DeallocateInstance
                
                                        

    !Management
    private ::      Ready
    private ::          LocateObjAdvectionDiffusion

    
    !Interfaces----------------------------------------------------------------

    private :: UngetAdvectionDiffusion3Dreal4
    private :: UngetAdvectionDiffusion3Dreal8
    interface  UngetAdvectionDiffusion
        module procedure UngetAdvectionDiffusion3Dreal4
        module procedure UngetAdvectionDiffusion3Dreal8
    end interface  UngetAdvectionDiffusion

    !Parameter-----------------------------------------------------------------
    integer, parameter :: MassConservation_           = 1
    integer, parameter :: ImposedValue_               = 2
    integer, parameter :: NullGradient_               = 4
    integer, parameter :: SubModel_                   = 5
    integer, parameter :: Orlanski_                   = 6
    integer, parameter :: MassConservNullGrad_        = 7
    integer, parameter :: CyclicBoundary_             = 8

    !Maximum internal wave velocity allowed is 10 m/s.
    real,    parameter :: MaxInternalCelerity = 10. 

    !Types---------------------------------------------------------------------

    type       T_State
        logical :: VertAdv      = ON    !This state defines if the coeficients D E and F need to be recalculated
        logical :: HorAdv       = ON    !This state defines if the coeficients D E and F need to be recalculated
        logical :: VertDif      = ON    !This state defines if 
        logical :: HorDif       = ON
        logical :: CellFluxes   = OFF
        logical :: OpenBoundary = OFF
    end type T_State

    type       T_FluxCoef
        real   , pointer, dimension(: , : , :)  :: C_flux    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: D_flux    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: E_flux    !Coeficient to calculate AdvFlux and DifFlux
        real   , pointer, dimension(: , : , :)  :: F_flux    !Coeficient to calculate AdvFlux and DifFlux
    end type T_FluxCoef

    type       T_CellFluxes
        real(8), pointer, dimension(:,:,:)      :: AdvFluxX         !Former CFLUX 
        real(8), pointer, dimension(:,:,:)      :: AdvFluxY         !Former CFLUX 
        real(8), pointer, dimension(:,:,:)      :: AdvFluxZ         !Former CFLUX 

        real(8), pointer, dimension(:,:,:)      :: DifFluxX          !Former DFLUX
        real(8), pointer, dimension(:,:,:)      :: DifFluxY          !Former DFLUX
        real(8), pointer, dimension(:,:,:)      :: DifFluxZ          !Former DFLUX
        type(T_Time)                            :: LastFluxCalculation           
    end type T_CellFluxes
    

    type       T_External
        real,    pointer, dimension(:,:,:) :: PROP                          
        real,    pointer, dimension(:,:,:) :: ReferenceProp                 
        real,    pointer, dimension(:,:,:) :: PROPOld
    
        !Map    
        integer, pointer, dimension(:,:,:) :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:) :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:) :: ComputeFacesW3D
        integer, pointer, dimension(:,:,:) :: LandPoints3D
        integer, pointer, dimension(:,:,:) :: OpenPoints3D
        integer, pointer, dimension(:,:  ) :: BoundaryPoints2D

        integer, pointer, dimension(:,:  ) :: KFloorZ
        real(8), pointer, dimension(:,:,:) :: VolumeZ
        real(8), pointer, dimension(:,:,:) :: VolumeZOld
        real,    pointer, dimension(:,:  ) :: DUX
        real,    pointer, dimension(:,:  ) :: DVY 
        real,    pointer, dimension(:,:  ) :: DZX 
        real,    pointer, dimension(:,:  ) :: DZY
        real,    pointer, dimension(:,:,:) :: DWZ
        real,    pointer, dimension(:,:,:) :: DZZ
        real,    pointer, dimension(:,:,:) :: AreaU
        real,    pointer, dimension(:,:,:) :: AreaV

        real                               :: Schmidt_H = null_real
        real                               :: SchmidtCoef_V       = null_real
        real                               :: SchmidtBackground_V = null_real
        real                               :: DecayTime = null_real
        real                               :: DTProp    = null_real

        integer                            :: BoundaryCondition = null_int


        !Implicit-Explicit weight coeficients -> 1 = Implicit, 0 = Explicit
        real :: ImpExp_AdvXX = null_real             
        real :: ImpExp_AdvYY = null_real             
        real :: ImpExp_DifH  = null_real             !Presentlty horizontal is explicitly computed
        real :: ImpExp_AdvV  = null_real      
        real :: ImpExp_DifV  = null_real      



        !Hydrodynamic
        logical                            :: Nulldif=.false. !ppina    
        real(8), pointer, dimension(:,:,:) :: Wflux_X
        real(8), pointer, dimension(:,:,:) :: Wflux_Y
        real(8), pointer, dimension(:,:,:) :: Wflux_Z


        !Turbulence
        real,    pointer, dimension(:,:,:) :: Visc_H
        real,    pointer, dimension(:,:,:) :: Diff_V  
        logical, pointer, dimension(:,:  ) :: SmallDepths

        logical                            :: SmallDepthsPresent

        real                               :: VolumeRelMax
        integer                            :: AdvMethodH, TVDLimitationH
        integer                            :: AdvMethodV, TVDLimitationV
        logical                            :: Upwind2H, Upwind2V

        !Discharges
        real,    pointer, dimension(:)     :: DischFlow, DischConc
        integer, pointer, dimension(:)     :: DischI, DischJ, DischK, DischnCells
        integer, pointer, dimension(:)     :: DischVert
        integer                            :: DischNumber  = null_int
        logical                            :: DischON     
        logical, pointer, dimension(:)     :: IgnoreDisch


    end type T_External

    !For performance reasons the coeffiecient are order in the way the are allocated
    type      T_AdvectionDiffusion
        real , dimension(:,:,:), pointer        :: Diffusion_CoeficientX
        real , dimension(:,:,:), pointer        :: Diffusion_CoeficientY
        real , dimension(:,:,:), pointer        :: Diffusion_CoeficientZ
        type(T_CellFluxes)                      :: Fluxes
        type(T_D_E_F)                           :: COEF3                    !Former DCOEF3  ECOEF3  FCOEF3
        type(T_FluxCoef)                        :: COEF3_VertAdv            !Vertical    advection coeficients
        type(T_FluxCoef)                        :: COEF3_HorAdvXX           !Horinzontal advection coeficients
        type(T_FluxCoef)                        :: COEF3_HorAdvYY           !Horinzontal advection coeficients
        real, pointer, dimension(: , : , :)     :: TICOEF3       
        real(8), pointer, dimension(:,:,:)      :: WaterFluxOBoundary
        !griflet
        type(T_THOMAS), pointer                 :: THOMAS
        real(8), pointer, dimension(:)          :: VECG                     !Auxiliar thomas arrays 
        real(8), pointer, dimension(:)          :: VECW                     !Auxiliar thomas arrays     

        !griflet
        integer                                 :: MaxThreads

        type(T_External  )                      :: ExternalVar

        integer                                 :: InstanceID
        type(T_Size3D    )                      :: Size
        type(T_Size3D    )                      :: WorkSize
        type(T_State     )                      :: State

        type(T_Time      )                      :: Now
        type(T_Time      )                      :: LastCalc  

        logical                                 :: Vertical1D        = .false.
        logical                                 :: XZFlow            = .false.
    
   
        !Instance of ModuleHorizontalMap
        integer                                 :: ObjHorizontalMap   = 0
       
        !Instance of ModuleHorizontalGrid
        integer                                 :: ObjHorizontalGrid  = 0
       
        !Instance of ModuleGeometry
        integer                                 :: ObjGeometry        = 0
    
        !Instance of ModuleTime
        integer                                 :: ObjTime            = 0

        !Collection of instances
        type(T_AdvectionDiffusion),  pointer    :: Next               => null()

    end type T_AdvectionDiffusion


    !Global Module Variables
    type (T_AdvectionDiffusion), pointer            :: FirstAdvectionDiffusion
    type (T_AdvectionDiffusion), pointer            :: Me
    
    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartAdvectionDiffusion(AdvectionDiffusionID,                            &
                                       GeometryID,                                      &
                                       HorizontalMapID,                                 &
                                       HorizontalGridID,                                &
                                       TimeID,                                          &
                                       Vertical1D,                                      &
                                       XZFlow,                                          &
                                       STAT)

        !Arguments-------------------------------------------------------------

        integer                             :: AdvectionDiffusionID
        integer                             :: GeometryID     
        integer                             :: HorizontalMapID
        integer                             :: HorizontalGridID
        integer                             :: TimeID    
        logical, optional, intent(IN )      :: Vertical1D, XZFlow           
        integer, optional, intent(OUT)      :: STAT     

        !External--------------------------------------------------------------

        integer :: STAT_CALL
        integer :: ready_         

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
 
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mAdvectionDiffusion_)) then
            nullify (FirstAdvectionDiffusion)
            call RegisterModule (mAdvectionDiffusion_) 
        endif
        

        call Ready(AdvectionDiffusionID, ready_)

cd0 :   if (ready_ == OFF_ERR_) then
            
            call AllocateInstance 
            
            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_,  HorizontalMapID )
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_,       GeometryID      )

            !Gets the size from the Geometry
            call GetGeometrySize(Me%ObjGeometry,                                         &
                                 Size        = Me%Size,                                  &
                                 WorkSize    = Me%WorkSize,                              &
                                 STAT        = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'StartAdvectionDiffusion - ModuleAdvectionDiffusion - ERR01'


            !Allocates variables
            call AllocateVariables

            if (present(Vertical1D)) Me%Vertical1D = Vertical1D

            if (present(XZFlow    )) Me%XZFlow     = XZFlow

            call null_time(Me%Now)        

            STAT_ = SUCCESS_

           !Returns ID
            AdvectionDiffusionID    = Me%InstanceID

        else

            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'StartAdvectionDiffusion - ModuleAdvectionDiffusion - ERR02'

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartAdvectionDiffusion

    !--------------------------------------------------------------------------

    subroutine AllocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_AdvectionDiffusion), pointer           :: NewAdvectionDiffusion
        type (T_AdvectionDiffusion), pointer           :: PreviousAdvectionDiffusion


        !Allocates new instance
        allocate (NewAdvectionDiffusion)
        nullify  (NewAdvectionDiffusion%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstAdvectionDiffusion)) then
            FirstAdvectionDiffusion         => NewAdvectionDiffusion
            Me                              => NewAdvectionDiffusion
        else
            PreviousAdvectionDiffusion      => FirstAdvectionDiffusion
            Me                              => FirstAdvectionDiffusion%Next
            do while (associated(Me))
                PreviousAdvectionDiffusion  => Me
                Me                          => Me%Next
            enddo
            Me                              => NewAdvectionDiffusion
            PreviousAdvectionDiffusion%Next => NewAdvectionDiffusion
        endif

        Me%InstanceID = RegisterNewInstance (mADVECTIONDIFFUSION_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------
    ! This subroutine allocates the auxiliar variables need
    ! to compute the Advection & Diffusion of a generic propriety 
    subroutine AllocateVariables()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                                     :: IJKLB, IJKUB
        
        !griflet
        integer                                     :: m
        type(T_VECGW), pointer                      :: VECGW

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB

        IJKLB = min (ILB, JLB, KLB)
        IJKUB = max (IUB, JUB, KUB)


        allocate(Me%Diffusion_CoeficientX   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Diffusion_CoeficientY   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Diffusion_CoeficientZ   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Fluxes%AdvFluxX         (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Fluxes%AdvFluxY         (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Fluxes%AdvFluxZ         (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%Fluxes%DifFluxX         (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Fluxes%DifFluxY         (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%Fluxes%DifFluxZ         (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%COEF3%D                 (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3%E                 (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3%F                 (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%COEF3_VertAdv%C_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_VertAdv%D_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_VertAdv%E_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_VertAdv%F_Flux    (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%COEF3_HorAdvXX%C_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_HorAdvXX%D_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_HorAdvXX%E_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_HorAdvXX%F_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%COEF3_HorAdvYY%C_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_HorAdvYY%D_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_HorAdvYY%E_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))
        allocate(Me%COEF3_HorAdvYY%F_Flux   (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%TICOEF3                 (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%WaterFluxOBoundary      (ILB:IUB, JLB:JUB, KLB:KUB))

        allocate(Me%VECG                    (IJKLB:IJKUB))
        allocate(Me%VECW                    (IJKLB:IJKUB))

        !griflet: BEGIN this is the alternate version that allows parallel openmp

        Me%MaxThreads = 1
        !$ Me%MaxThreads = omp_get_max_threads()

        allocate(Me%THOMAS)
        allocate(Me%THOMAS%COEF3)
        allocate(Me%THOMAS%VEC(1:Me%MaxThreads))

        do m = 1, Me%MaxThreads

            VECGW => Me%THOMAS%VEC(m)

            allocate(VECGW%G(IJKLB:IJKUB))
            allocate(VECGW%W(IJKLB:IJKUB))

        enddo

        Me%THOMAS%COEF3 => Me%COEF3
        Me%THOMAS%TI => Me%TICOEF3

        !griflet: END

        Me%Diffusion_CoeficientX    = Null_real
        Me%Diffusion_CoeficientY    = Null_real
        Me%Diffusion_CoeficientZ    = Null_real
                            
        Me%Fluxes%AdvFluxX          = Null_real
        Me%Fluxes%AdvFluxY          = Null_real
        Me%Fluxes%AdvFluxZ          = Null_Real

        Me%Fluxes%DifFluxX          = Null_real
        Me%Fluxes%DifFluxY          = Null_real
        Me%Fluxes%DifFluxZ          = Null_Real

        Me%COEF3%D                  = Null_real
        Me%COEF3%E                  = Null_real
        Me%COEF3%F                  = Null_real

        Me%COEF3_VertAdv%C_Flux     = Null_real
        Me%COEF3_VertAdv%D_Flux     = Null_real
        Me%COEF3_VertAdv%E_Flux     = Null_real
        Me%COEF3_VertAdv%F_Flux     = Null_real

        Me%COEF3_HorAdvXX%C_Flux    = Null_real
        Me%COEF3_HorAdvXX%D_Flux    = Null_real
        Me%COEF3_HorAdvXX%E_Flux    = Null_real
        Me%COEF3_HorAdvXX%F_Flux    = Null_real

        Me%COEF3_HorAdvYY%C_Flux    = Null_real
        Me%COEF3_HorAdvYY%D_Flux    = Null_real
        Me%COEF3_HorAdvYY%E_Flux    = Null_real
        Me%COEF3_HorAdvYY%F_Flux    = Null_real

        Me%TICOEF3                  = Null_real
        Me%WaterFluxOBoundary       = Null_real

        Me%VECG                     = Null_real 
        Me%VECW                     = Null_real 


        !----------------------------------------------------------------------

    end subroutine AllocateVariables   


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine GetAdvFlux(AdvectionDiffusionID, AdvFluxX,                  &
                                                AdvFluxY,                  &
                                                AdvFluxZ, STAT)    !AdvFlux is the former CFLUX

        !Arguments-------------------------------------------------------------
        integer                                      :: AdvectionDiffusionID

        real(8), pointer, optional, dimension(:,:,:) :: AdvFluxX
        real(8), pointer, optional, dimension(:,:,:) :: AdvFluxY
        real(8), pointer, optional, dimension(:,:,:) :: AdvFluxZ

        integer,      optional, intent(OUT) :: STAT


        !External--------------------------------------------------------------

        integer :: ready_          
        integer :: STAT_CALL

        !Local-----------------------------------------------------------------

        integer :: STAT_   
           

        !----------------------------------------------------------------------


        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                            &
            (ready_ == READ_LOCK_ERR_)) then

            call GetComputeCurrentTime(Me%ObjTime, Me%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'GetAdvFlux - ModuleAdvectionDiffusion - ERR01'



cd5 :       if (Me%Now == Me%Fluxes%LastFluxCalculation) then

cd2 :           if (present(AdvFluxX)) then
                    call Read_Lock(mADVECTIONDIFFUSION_,  Me%InstanceID)   
                    AdvFluxX => Me%Fluxes%AdvFluxX
                end if cd2


cd3 :           if (present(AdvFluxY)) then
                    call Read_Lock(mADVECTIONDIFFUSION_,  Me%InstanceID)      
                    AdvFluxY => Me%Fluxes%AdvFluxY
                end if cd3


cd4 :           if (present(AdvFluxZ)) then
                    call Read_Lock(mADVECTIONDIFFUSION_,  Me%InstanceID)
                    AdvFluxZ => Me%Fluxes%AdvFluxZ
                end if cd4
                
                STAT_ = SUCCESS_
            else 
                STAT_ = TIME_ERR_
            end if cd5

            call null_time(Me%Now)
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetAdvFlux

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
 
    subroutine GetDifFlux(AdvectionDiffusionID, DifFluxX,                     &
                                                DifFluxY,                     &
                                                DifFluxZ, STAT)      !DifFlux is the former CFLUX

        !Arguments-------------------------------------------------------------
        integer                                      :: AdvectionDiffusionID

        real(8), pointer, optional, dimension(:,:,:) :: DifFluxX
        real(8), pointer, optional, dimension(:,:,:) :: DifFluxY
        real(8), pointer, optional, dimension(:,:,:) :: DifFluxZ

        integer,      optional, intent(OUT) :: STAT


        !External--------------------------------------------------------------

        integer :: ready_          
        integer :: STAT_CALL

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)
        
cd1 :   if ((ready_ == IDLE_ERR_     ) .OR.                                            &
            (ready_ == READ_LOCK_ERR_)) then

            call GetComputeCurrentTime(Me%ObjTime,                    &
                                       Me%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'GetDifFlux - ModuleAdvectionDiffusion - ERR01'

cd5 :       if (Me%Now ==                                           &
                Me%Fluxes%LastFluxCalculation) then
cd2 :           if (present(DifFluxX)) then
                    call Read_Lock(mADVECTIONDIFFUSION_,  Me%InstanceID)
                    DifFluxX => Me%Fluxes%DifFluxX
                end if cd2
           

cd3 :           if (present(DifFluxY)) then
                    call Read_Lock(mADVECTIONDIFFUSION_,  Me%InstanceID)
                    DifFluxY => Me%Fluxes%DifFluxY
                end if cd3
           

cd4 :           if (present(DifFluxZ)) then
                    call Read_Lock(mADVECTIONDIFFUSION_,  Me%InstanceID)   
                    DifFluxZ => Me%Fluxes%DifFluxZ
                end if cd4
                STAT_ = SUCCESS_
            else
                STAT_ = TIME_ERR_
            end if cd5
           
            call null_time(Me%Now)
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                               &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDifFlux

    !--------------------------------------------------------------------------

    subroutine GetBoundaryConditionList(MassConservation,                     &
                                        ImposedValue,                         &
                                        NullGradient, SubModel, Orlanski,     &
                                        MassConservNullGrad, CyclicBoundary)

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: MassConservation      
        integer, optional, intent(OUT) :: ImposedValue       
        integer, optional, intent(OUT) :: NullGradient        
        integer, optional, intent(OUT) :: SubModel        
        integer, optional, intent(OUT) :: Orlanski
        integer, optional, intent(OUT) :: MassConservNullGrad
        integer, optional, intent(OUT) :: CyclicBoundary

        !----------------------------------------------------------------------
     
        if (present(MassConservation          )) MassConservation           = MassConservation_
        if (present(ImposedValue              )) ImposedValue               = ImposedValue_
        if (present(NullGradient              )) NullGradient               = NullGradient_
        if (present(SubModel                  )) SubModel                   = SubModel_
        if (present(Orlanski                  )) Orlanski                   = Orlanski_                            
        if (present(MassConservNullGrad       )) MassConservNullGrad        = MassConservNullGrad_
        if (present(CyclicBoundary            )) CyclicBoundary             = CyclicBoundary_
        !----------------------------------------------------------------------

    end subroutine GetBoundaryConditionList

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine UngetAdvectionDiffusion3Dreal4(AdvectionDiffusionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: AdvectionDiffusionID

        real(4), pointer, dimension(:,:,:) :: Array

        integer, optional, intent (OUT)    :: STAT
   

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)

cd1 :   if (ready_ == READ_LOCK_ERR_) then
            nullify(Array)

            call Read_Unlock(mADVECTIONDIFFUSION_, Me%InstanceID,"UngetAdvectionDiffusion3Dreal4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAdvectionDiffusion3Dreal4

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine UngetAdvectionDiffusion3Dreal8(AdvectionDiffusionID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: AdvectionDiffusionID

        real(8), pointer, dimension(:,:,:) :: Array

        integer, optional, intent (OUT) :: STAT
   
        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)

cd1 :   if (ready_ == READ_LOCK_ERR_) then
            nullify(Array)

            call Read_Unlock(mADVECTIONDIFFUSION_, Me%InstanceID,"UngetAdvectionDiffusion3Dreal4")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetAdvectionDiffusion3Dreal8

        !----------------------------------------------------------------------

    subroutine SetDischarges (AdvectionDiffusionID, DischFlow, DischConc,               &
                              DischI, DischJ, DischK, DischVert, DischNumber,           &
                              IgnoreDisch, DischnCells, STAT)
    
    
        !Arguments--------------------------------------------------------------
        integer,           intent(IN )     :: AdvectionDiffusionID
        real,    pointer, dimension(:)     :: DischFlow, DischConc
        integer, pointer, dimension(:)     :: DischI, DischJ, DischK
        integer, pointer, dimension(:)     :: DischVert 
        integer,           intent(IN )     :: DischNumber
        logical, pointer, dimension(:)     :: IgnoreDisch
        integer, pointer, dimension(:)     :: DischnCells
        integer, optional, intent(OUT)     :: STAT

        !Local-----------------------------------------------------------------

        integer :: STAT_            
        integer :: ready_              

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)

cd1 :   if (ready_ == IDLE_ERR_) then

            Me%ExternalVar%DischON      = .true.

            Me%ExternalVar%DischNumber  =  DischNumber

            Me%ExternalVar%DischFlow    => DischFlow 
            Me%ExternalVar%DischConc    => DischConc 
            Me%ExternalVar%DischI       => DischI 
            Me%ExternalVar%DischJ       => DischJ 
            Me%ExternalVar%DischK       => DischK 
            Me%ExternalVar%DischVert    => DischVert
            Me%ExternalVar%IgnoreDisch  => IgnoreDisch
            Me%ExternalVar%DischnCells  => DischnCells


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                              &
            STAT = STAT_

    end subroutine SetDischarges 

        !----------------------------------------------------------------------

        !----------------------------------------------------------------------

    subroutine UnSetDischarges (AdvectionDiffusionID, STAT)
    
    
        !Arguments--------------------------------------------------------------
        integer,           intent(IN )     :: AdvectionDiffusionID
        integer, optional, intent(OUT)     :: STAT

        !Local-----------------------------------------------------------------

        integer :: STAT_            
        integer :: ready_              

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)

cd1 :   if (ready_ == IDLE_ERR_) then

            Me%ExternalVar%DischON      = .false.

            Me%ExternalVar%DischNumber  = null_int

            nullify(Me%ExternalVar%DischFlow)
            nullify(Me%ExternalVar%DischConc)
            nullify(Me%ExternalVar%DischI   )
            nullify(Me%ExternalVar%DischJ   )
            nullify(Me%ExternalVar%DischK   )
            nullify(Me%ExternalVar%DischVert)
            nullify(Me%ExternalVar%IgnoreDisch)
            nullify(Me%ExternalVar%DischnCells)


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                              &
            STAT = STAT_

    end subroutine UnSetDischarges 

        !----------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    !   AdvectionDiffusion Resolve a equacao de adveccao-difusao                      
    !                      na sua forma Euleriana com adveccao vertical e 
    !                      difusao vertical implicitas. Considera-se que a                      
    !                      propriedade encontra-se no centro das celulas e que os                     
    !                      fluxos de agua nas suas faces sao conhecidos.       


    subroutine AdvectionDiffusion(AdvectionDiffusionID,                         &
                                  PROP,                                         &    
                                  schmidt_H,                                    &
                                  SchmidtCoef_V,                                &
                                  SchmidtBackground_V,                          &
                                  AdvMethodH,                                   &
                                  TVDLimitationH,                               &
                                  AdvMethodV,                                   &
                                  TVDLimitationV,                               &
                                  Upwind2H,                                     &
                                  Upwind2V,                                     &
                                  VolumeRelMax,                                 &
                                  DTProp,                                       &
                                  ImpExp_AdvV, ImpExp_DifV,                     &
                                  ImpExp_AdvXX, ImpExp_AdvYY,                   &
                                  ImpExp_DifH, NullDif,                         &
                                  Wflux_X, Wflux_Y, Wflux_Z,                    &
                                  VolumeZOld, VolumeZ,                          &
                                  OpenPoints3D,                                 &
                                  LandPoints3D,                                 &
                                  ComputeFacesU3D, ComputeFacesV3D,             &
                                  ComputeFacesW3D,                              &
                                  Visc_H, Diff_V,                               &
                                  CellFluxes,                                   &
                                  ReferenceProp,                                &
                                  BoundaryCondition,                            &
                                  DecayTime,                                    &
                                  NumericStability,                             &   
                                  PROPOld,                                      &
                                  SmallDepths,                                  &
                                  STAT)

        !Arguments-------------------------------------------------------------
        integer                            :: AdvectionDiffusionID

        integer, optional, intent(OUT)     :: STAT

        logical, optional, intent(IN )     :: NumericStability

        real,                          pointer, dimension(:,:,:) :: PROP
        real,    optional,             pointer, dimension(:,:,:) :: ReferenceProp
        real,    optional,             pointer, dimension(:,:,:) :: PROPOld

        !Hydrodynamic
        real(8),                       pointer, dimension(:,:,:) :: Wflux_X, Wflux_Y, Wflux_Z
        real(8),                       pointer, dimension(:,:,:) :: VolumeZOld, VolumeZ
        real,                          pointer, dimension(:,:,:) :: Visc_H
        real,                          pointer, dimension(:,:,:) :: Diff_V
        integer, dimension(:, :, :), pointer                     :: OpenPoints3D
        integer, dimension(:, :, :), pointer                     :: LandPoints3D
        integer, dimension(:, :, :), pointer                     :: ComputeFacesU3D
        integer, dimension(:, :, :), pointer                     :: ComputeFacesV3D
        integer, dimension(:, :, :), pointer                     :: ComputeFacesW3D

        real,              intent(IN )     :: schmidt_H
        real,              intent(IN)      :: SchmidtCoef_V
        real,              intent(IN)      :: SchmidtBackground_V 
        integer,           intent(IN)      :: AdvMethodH, TVDLimitationH
        integer,           intent(IN)      :: AdvMethodV, TVDLimitationV
        logical,           intent(IN)      :: Upwind2H, Upwind2V
        real   ,           intent(IN)      :: VolumeRelMax 
        real,    optional, intent(IN )     :: DecayTime 
        real,              intent(IN )     :: DTProp
        real,              intent(IN )     :: ImpExp_DifV, ImpExp_AdvV 
        real,              intent(IN )     :: ImpExp_AdvXX, ImpExp_AdvYY, ImpExp_DifH
        Logical,           intent(IN )     :: NullDif
 
        logical, optional, intent(IN )     :: CellFluxes                    !State variable to calculate Difusive 
                                                                            ! and Advective fluxes among cells.

        integer, optional, intent(IN )     :: BoundaryCondition             !State variable to impose
                                                                            ! boundary conditions 

        logical, dimension(:, : ), pointer, optional :: SmallDepths

        !External--------------------------------------------------------------

        integer :: STAT_CALL
        integer :: ready_              

        logical :: NumericStability_
        logical :: CellFluxes_                

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable
        integer :: ILB, IUB, JLB, JUB, KLB, KUB

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)

cd1 :   if (ready_ == IDLE_ERR_) then

            if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusion")

            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime, Me%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_)                                        &
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR01'

            ILB = Me%Size%ILB
            IUB = Me%Size%IUB

            JLB = Me%Size%JLB 
            JUB = Me%Size%JUB

            KLB = Me%Size%KLB
            KUB = Me%Size%KUB

            if ((ImpExp_AdvXX == ImplicitScheme .or. ImpExp_AdvYY == ImplicitScheme) .and.  &
                (AdvMethodH   == UpwindOrder2   .or. AdvMethodH   == UpwindOrder3  )) then
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR100'
            endif

            if ((ImpExp_AdvV  == ImplicitScheme                                    ) .and.  &
                (AdvMethodV   == UpwindOrder2   .or. AdvMethodV   == UpwindOrder3  )) then
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR200'
            endif


cd12 :      if (present(ReferenceProp)) then
cd6  :          if (associated(ReferenceProp)) then
                    Me%ExternalVar%ReferenceProp => ReferenceProp
                else
                    nullify(Me%ExternalVar%ReferenceProp)
                end if cd6
            else
                nullify(Me%ExternalVar%ReferenceProp)
            end if cd12


cd62 :      if (present(PROPOld)) then
cd65  :          if (associated(PROPOld)) then
                    Me%ExternalVar%PROPOld => PROPOld
                else
                    nullify(Me%ExternalVar%PROPOld)
                end if cd65
            else 
                nullify(Me%ExternalVar%PROPOld)
            end if cd62



cd99 :      if (present(NumericStability)) then
                NumericStability_ = NumericStability
            else cd99
                NumericStability_ =.FALSE.
            end if cd99

            Me%ExternalVar%PROP               => PROP
            Me%ExternalVar%Wflux_X            => Wflux_X
            Me%ExternalVar%Wflux_Y            => Wflux_Y
            Me%ExternalVar%Wflux_Z            => Wflux_Z
            Me%ExternalVar%Visc_H             => Visc_H
            Me%ExternalVar%Diff_V             => Diff_V
            Me%ExternalVar%VolumeZOld         => VolumeZOld
            Me%ExternalVar%VolumeZ            => VolumeZ
            Me%ExternalVar%OpenPoints3D       => OpenPoints3D
            Me%ExternalVar%LandPoints3D       => LandPoints3D
            Me%ExternalVar%ComputeFacesU3D    => ComputeFacesU3D
            Me%ExternalVar%ComputeFacesV3D    => ComputeFacesV3D
            Me%ExternalVar%ComputeFacesW3D    => ComputeFacesW3D

            if (present(SmallDepths)) then
                Me%ExternalVar%SmallDepthsPresent = .true.
                Me%ExternalVar%SmallDepths        => SmallDepths
            else
                Me%ExternalVar%SmallDepthsPresent = .false.
            endif

            Me%ExternalVar%ImpExp_DifV      = ImpExp_DifV   
            Me%ExternalVar%ImpExp_AdvXX     = ImpExp_AdvXX 
            Me%ExternalVar%ImpExp_AdvYY     = ImpExp_AdvYY
            Me%ExternalVar%ImpExp_AdvV      = ImpExp_AdvV

            Me%ExternalVar%ImpExp_DifH      = ImpExp_DifH

            Me%ExternalVar%Upwind2H         = Upwind2H
            Me%ExternalVar%Upwind2V         = Upwind2V
            Me%ExternalVar%VolumeRelMax     = VolumeRelMax 


cd9 :       if (present(DecayTime)) then
                Me%ExternalVar%DecayTime = DecayTime
            else
                Me%ExternalVar%DecayTime = null_real
            end if cd9


cd10 :      if (present(BoundaryCondition)) then
                Me%ExternalVar%BoundaryCondition = BoundaryCondition
            else
                Me%ExternalVar%BoundaryCondition = null_int
            end if cd10



cd7 :       if (ImpExp_DifH  /= 0.0) then    !0 = Explicit
                write(*,*) 'Horizontal Diffusion must be explicit.'
                stop       'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR02'
            end if cd7


cd66 :      if (ImpExp_AdvXX == ImplicitScheme .and. ImpExp_AdvYY == ImplicitScheme) then    
                write(*,*) 'Horizontal Advection can not be implicit in both directions.'
                stop       'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR03'
            end if cd66


            !DUX, DVY, DZX, DZY
            call GetHorizontalGrid(Me%ObjHorizontalGrid,              &
                                   DUX = Me%ExternalVar%DUX,          &
                                   DVY = Me%ExternalVar%DVY,          &
                                   DZX = Me%ExternalVar%DZX,          &
                                   DZY = Me%ExternalVar%DZY,          &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR04'

            !BoundaryPoints2D
            call GetBoundaries(Me%ObjHorizontalMap,                   &
                               Me%ExternalVar%BoundaryPoints2D,       &  
                               STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR05'

            !KFloorZ
            call GetGeometryKFloor(Me%ObjGeometry,                    &
                                   Z = Me%ExternalVar%KFloorZ,        &
                                   STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR06'

            !AreaU, AreaV
            call GetGeometryAreas(Me%ObjGeometry,                     &
                                  AreaU = Me%ExternalVar%AreaU,       &
                                  AreaV = Me%ExternalVar%AreaV,       &
                                  STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR08'

            !DZZ, DWZ
            call GetGeometryDistances(Me%ObjGeometry,                 &
                                      DZZ = Me%ExternalVar%DZZ,       &
                                      DWZ = Me%ExternalVar%DWZ,       &
                                      STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                stop 'AdvectionDiffusion - ModuleAdvectionDiffusion - ERR09'

cd11 :      if (present(CellFluxes)) then
                CellFluxes_ = CellFluxes
            else
                CellFluxes_ = .FALSE.
            end if cd11

            call Set_Internal_State(DTProp, CellFluxes_, AdvMethodH, TVDLimitationH, AdvMethodV, TVDLimitationV, &
                                    SchmidtCoef_V, SchmidtBackground_V, Schmidt_H, NullDif)

            Me%ExternalVar%DTProp              = DTProp
            Me%ExternalVar%AdvMethodH          = AdvMethodH
            Me%ExternalVar%TVDLimitationH      = TVDLimitationH
            Me%ExternalVar%SchmidtCoef_V       = SchmidtCoef_V    
            Me%ExternalVar%SchmidtBackground_V = SchmidtBackground_V
            Me%ExternalVar%NullDif            =  NullDif
            Me%ExternalVar%Schmidt_H           = Schmidt_H
            Me%ExternalVar%DTProp              = DTProp
            
            if (Me%State%VertDif) then
                call Convert_Dif_Vertical()
            endif

            if (Me%State%HorDif) then
                call Convert_Visc_Dif_Horizontal()
            endif

cd8 :       if (Me%State%VertAdv) then

                Me%ExternalVar%AdvMethodV       = AdvMethodV
                Me%ExternalVar%TVDLimitationV   = TVDLimitationV

                call SetMatrixValue (Me%COEF3_VertAdv%C_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_VertAdv%D_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_VertAdv%E_flux, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3_VertAdv%F_flux, Me%Size, 0.0)
                    
            end if cd8

cd5 :       if (Me%State%CellFluxes) then
                
                Me%Fluxes%LastFluxCalculation = Me%Now  !Now is former TIMEI

                call SetMatrixValue (Me%Fluxes%AdvFluxX, Me%Size, dble(0.0))
                call SetMatrixValue (Me%Fluxes%AdvFluxY, Me%Size, dble(0.0))
                call SetMatrixValue (Me%Fluxes%AdvFluxZ, Me%Size, dble(0.0))

                call SetMatrixValue (Me%Fluxes%DifFluxX, Me%Size, dble(0.0))
                call SetMatrixValue (Me%Fluxes%DifFluxY, Me%Size, dble(0.0))
                call SetMatrixValue (Me%Fluxes%DifFluxZ, Me%Size, dble(0.0))

            end if cd5


            call AdvectionDiffusionIteration(ImpExp_AdvXX, ImpExp_AdvYY)                     !AdvectionDiffusion main cycle

            call FinishAdvectionDiffusionIt()                  !Nullify external variables

            Me%LastCalc = Me%Now

            call null_time   (Me%Now)

            STAT_ = SUCCESS_

            if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusion")

        else               
            STAT_ = ready_
        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine AdvectionDiffusion

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine AdvectionDiffusionIteration(ImpExp_AdvXX, ImpExp_AdvYY)

        !Arguments-------------------------------------------------------------
        real,    intent(IN)                 :: ImpExp_AdvXX, ImpExp_AdvYY

        !Local-----------------------------------------------------------------
        integer                             :: ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: i, j, k, di, dj
        integer                             :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration")

        ILBWS = Me%WorkSize%ILB
        IUBWS = Me%WorkSize%IUB

        JLBWS = Me%WorkSize%JLB
        JUBWS = Me%WorkSize%JUB

        KLBWS = Me%WorkSize%KLB
        KUBWS = Me%WorkSize%KUB

        ILB   = Me%Size%ILB
        IUB   = Me%Size%IUB

        JLB   = Me%Size%JLB 
        JUB   = Me%Size%JUB

        KLB   = Me%Size%KLB
        KUB   = Me%Size%KUB

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration_1")

        !Inicializacao dos coeficientes DCOEF3,ECOEF3,FCOEF3 e independent term
        call SetMatrixValue (Me%COEF3%D, Me%Size, 0.0)
        call SetMatrixValue (Me%COEF3%E, Me%Size, dble(1.0))
        call SetMatrixValue (Me%COEF3%F, Me%Size, 0.0)
        call SetMatrixValue (Me%TICOEF3, Me%Size, 0.0)
     
cd2 :   if (Me%State%HorAdv) then

            call SetMatrixValue (Me%COEF3_HorAdvXX%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvXX%F_flux, Me%Size, 0.0)

            call SetMatrixValue (Me%COEF3_HorAdvYY%C_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%D_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%E_flux, Me%Size, 0.0)
            call SetMatrixValue (Me%COEF3_HorAdvYY%F_flux, Me%Size, 0.0)

        end if cd2

        ! evolutive term + Volume Variation
        call VolumeVariation()

        if (Me%ExternalVar%DischON) call Discharges ()

        if (.not. Me%Vertical1D) then

            ! Calculo dos fluxos difusivos horizontais da propriedade e actualizacao                         
            ! imediata do termo independente
            call HorizontalDiffusion()

            ! Calculo dos fluxos advectivos da propriedade 
            call HorizontalAdvection(ImpExp_AdvXX, ImpExp_AdvYY)

        endif


        if (KUBWS > 1) then
            call VerticalDiffusion ()
            if (.not. Me%Vertical1D) call VerticalAdvection()
        endif
        

        ! This subroutine must be always the last to be called
        if (Me%State%OpenBoundary)                                                      &
            call OpenBoundaryCondition() 

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration_1")

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration_2")

        CHUNK = CHUNK_J(JLB,JUB)
        !$OMP PARALLEL PRIVATE(i,j,k)

        ! At this stage the variable TICOEF3 have null values in the land points
        ! We want that all proprieties in land points have the value of Null_Real.
dok7 :  do k = KLB, KUB
        !$OMP DO SCHEDULE (DYNAMIC,CHUNK)
doj7 :  do j = JLB, JUB
doi7 :  do i = ILB, IUB

            if (k == KUB .or. Me%ExternalVar%LandPoints3D(i,j,k) == 1) then
                Me%TICOEF3(i,j,KUB)= Null_Real
!            else
!                Me%TICOEF3(i,j,k) = Me%TICOEF3(i,j,k)                                       &
!                     * (1 - Me%ExternalVar%LandPoints3D(i,j,k))                             &
!                     +      Me%ExternalVar%LandPoints3D(i,j,k) * Null_Real
            endif

                 
        end do doi7
        end do doj7
        !$OMP END DO NOWAIT
        end do dok7

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration_2")

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration_TH")

cd3:    if (KUBWS == 1 .and. ImpExp_AdvXX == ImplicitScheme) then !ImplicitScheme = 0

            di = 0
            dj = 1

            !griflet: old call
            !call THOMAS_3D(ILBWS, IUBWS,                                                &
            !               JLBWS, JUBWS,                                                &
            !               KLBWS, KUBWS,                                                &
            !               di, dj,                                                      &
            !               Me%COEF3%D,                                                  &
            !               Me%COEF3%E,                                                  &
            !               Me%COEF3%F,                                                  &
            !               Me%TICOEF3,                                                  &
            !               Me%ExternalVar%PROP,                                         &
            !               Me%VECG,                                                     &
            !               Me%VECW)      
            !griflet: new  call
            call THOMAS_3D(ILBWS, IUBWS,                                                &
                           JLBWS, JUBWS,                                                &
                           KLBWS, KUBWS,                                                &
                           di, dj,                                                      &
                           Me%THOMAS,                                                   &
                           Me%ExternalVar%PROP)      
                            
        else if (KUBWS == 1 .and. ImpExp_AdvYY == ImplicitScheme) then cd3 !ImplicitScheme = 0

            di = 1
            dj = 0
            
            !griflet: old call
            !call THOMAS_3D(JLBWS, JUBWS,                                                &
            !               ILBWS, IUBWS,                                                &
            !               KLBWS, KUBWS,                                                &
            !               di, dj,                                                      &
            !               Me%COEF3%D,                                                  &
            !               Me%COEF3%E,                                                  &
            !               Me%COEF3%F,                                                  &
            !               Me%TICOEF3,                                                  &
            !               Me%ExternalVar%PROP,                                         &
            !               Me%VECG,                                                     &
            !               Me%VECW)      
            !griflet: new call                           
            call THOMAS_3D(JLBWS, JUBWS,                                                &
                           ILBWS, IUBWS,                                                &
                           KLBWS, KUBWS,                                                &
                           di, dj,                                                      &
                           Me%THOMAS,                                                   &
                           Me%ExternalVar%PROP)      
        else cd3
 
            ! If the model is 3D the vertical diffusion must be implicit so is necessary to 
            ! compute the vertical diffusion  implicitly
            
            !griflet: old call   
            !CALL THOMASZ(ILBWS, IUBWS,                                                  &
            !             JLBWS, JUBWS,                                                  &
            !             KLBWS, KUBWS,                                                  &
            !             Me%COEF3%D,                                                    &
            !             Me%COEF3%E,                                                    &
            !             Me%COEF3%F,                                                    &
            !             Me%TICOEF3,                                                    &
            !             Me%ExternalVar%PROP,                                           &
            !             Me%VECG,                                                       &
            !             Me%VECW)      
        
            !griflet: new call
            CALL THOMASZ(ILBWS, IUBWS,                                                  &
                         JLBWS, JUBWS,                                                  &
                         KLBWS, KUBWS,                                                  &
                         Me%THOMAS,                                                     &
                         Me%ExternalVar%PROP)

        endif cd3

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration_TH")


        if      (Me%ExternalVar%BoundaryCondition == NullGradient_) then

            call ImposeNullGradient()

        else if (Me%ExternalVar%BoundaryCondition == CyclicBoundary_) then

            call Prop_CyclicBoundary()

        endif


cd5 :   if (Me%State%CellFluxes) then
            if (Me%ExternalVar%ImpExp_AdvV > 0.0 .and. KUBWS > 1)                       &
                call CalcVerticalAdvFlux(Me%ExternalVar%ImpExp_AdvV)


            if (Me%ExternalVar%ImpExp_DifV > 0.0 .and. KUBWS > 1)                       &
                call CalcVerticalDifFlux (Me%ExternalVar%ImpExp_DifV )


            if (Me%ExternalVar%ImpExp_AdvXX == ImplicitScheme)                          &
                call CalcHorizontalAdvFluxXX(Me%ExternalVar%ImpExp_AdvXX)

            if (Me%ExternalVar%ImpExp_AdvYY == ImplicitScheme)                          &
                call CalcHorizontalAdvFluxYY(Me%ExternalVar%ImpExp_AdvYY)


        end if cd5

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "AdvectionDiffusionIteration")

        !----------------------------------------------------------------------

    end subroutine AdvectionDiffusionIteration

    !--------------------------------------------------------------------------

    subroutine ImposeNullGradient()


        !Arguments--------------------------------------------------------------

        !Local------------------------------------------------------------------
        real     :: BoundaryProp
        logical  :: Found
        integer  :: i, j, k, ILB, IUB, JLB, JUB, KLB, KUB
        integer  :: CHUNK

        !Begin------------------------------------------------------------------


        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KUB = Me%WorkSize%KUB

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "ImposeNullGradient")

        CHUNK = CHUNK_J(JLB,JUB)

        !$OMP PARALLEL PRIVATE(j,i,k,KLB,BoundaryProp,Found)

        !$OMP DO SCHEDULE(GUIDED,CHUNK)
do3 :   do  j = JLB, JUB
do2 :   do  i = ILB, IUB

cd1 :       if (Me%ExternalVar%BoundaryPoints2D(i,j) == 1) then 

                KLB = ABS(Me%ExternalVar%KFloorZ(i,j))
                

do1 :           do  k = KLB, KUB


                    call NullGradProp(BoundaryProp, i, j, k, Found)

                    if (Found) then

                        Me%ExternalVar%Prop(i, j, k) = BoundaryProp 

                    endif


                enddo do1


            endif cd1

        enddo do2
        enddo do3
        !$OMP END DO

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "ImposeNullGradient")

    end subroutine ImposeNullGradient


    subroutine NullGradProp(BoundaryProp, i, j, k, Found)

        !Arguments--------------------------------------------------------------
        integer, intent(IN)                 :: i, j, k  
        real,    intent(OUT)                :: BoundaryProp
        logical, intent(OUT)                :: Found

        !Local------------------------------------------------------------------
        integer  :: Aux

        !Begin------------------------------------------------------------------

        Aux = Me%ExternalVar%ComputeFacesV3D(i + 1, j, k)+&
              Me%ExternalVar%ComputeFacesV3D(i    , j, k)+&  
              Me%ExternalVar%ComputeFacesU3D(i, j + 1, k)+&  
              Me%ExternalVar%ComputeFacesU3D(i    , j, k)

cd2:    if (Aux > 0) then

            BoundaryProp =                                                               &
               (Me%ExternalVar%Prop           (i + 1,     j, k)    *  &
                Me%ExternalVar%ComputeFacesV3D(i + 1,     j, k)    +  &
                Me%ExternalVar%Prop           (i - 1,     j, k)    *  &
                Me%ExternalVar%ComputeFacesV3D(i    ,     j, k)    +  &
                Me%ExternalVar%Prop           (i    , j + 1, k)    *  &
                Me%ExternalVar%ComputeFacesU3D(i    , j + 1, k)    +  &
                Me%ExternalVar%Prop           (i    , j - 1, k)    *  &
                Me%ExternalVar%ComputeFacesU3D(i    ,     j, k))   /  &
                real(Aux)

            Found = .true.

        else 

            Found = .false.

        endif cd2


    end subroutine NullGradProp

    !--------------------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !
    ! This subroutine imposed Cyclic conditions in the open boundary                         !
    !                                                                                      !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Subroutine Prop_CyclicBoundary ()

        !Arguments------------------------------------------------------------

        !Local---------------------------------------------------------------------
        integer, dimension(:,:  ), pointer  :: BoundaryPoints2D, KFloorZ

        real,    dimension(:,:,:), pointer  :: Prop, PropRef

        real                                :: DT_Prop

        integer                             :: IUB, ILB, JUB, JLB, KUB, KLB
        integer                             :: i, j, k, kbottom

        !Begin----------------------------------------------------------------

        !Begin - Shorten variables name 

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        IUB = Me%WorkSize%IUB
        ILB = Me%WorkSize%ILB
        JUB = Me%WorkSize%JUB
        JLB = Me%WorkSize%JLB

        DT_Prop          =  Me%ExternalVar%DTProp

        BoundaryPoints2D => Me%ExternalVar%BoundaryPoints2D
        KFloorZ          => Me%ExternalVar%KFloorZ

        Prop             => Me%ExternalVar%Prop
        PropRef          => Me%ExternalVar%ReferenceProp


       
        !End   - Shorten variables name

do6:    do j = JLB, JUB     
do5:    do i = ILB, IUB     

            if (BoundaryPoints2D(i, j) == Boundary) Prop(i, j, KLB:KUB) = PropRef(i, j, KLB:KUB)

        enddo do5
        enddo do6
       
do1:    do  i = ILB + 1, IUB - 1

cd1:        if (BoundaryPoints2D(i, JLB) == Boundary .and.                               &
                BoundaryPoints2D(i, JUB) == Boundary) then

                kbottom = KFloorZ(i, JUB - 1)

do3:            do k = kbottom, KUB

                    Prop(i, JLB, k) = Prop(i, JUB-1, k)

                enddo do3

                kbottom = KFloorZ(i, JLB + 1)

do7:            do k = kbottom, KUB

                    Prop(i, JUB, k) = Prop(i, JLB+1, k)

                enddo do7

           
            endif cd1
                           
        enddo do1

        
do2:    do  j = JLB + 1, JUB - 1

cd2:        if (BoundaryPoints2D(ILB, j) == Boundary .and.                               &
                BoundaryPoints2D(IUB, j) == Boundary) then

                kbottom = KFloorZ(IUB-1, j)

do4:            do k=kbottom, KUB

                    Prop(ILB, j, k) = Prop(IUB-1, j, k)

                enddo do4

                kbottom = KFloorZ(ILB+1, j)

do8:            do k=kbottom, KUB

                    Prop(IUB, j, k) = Prop(ILB+1, j, k)

                enddo do8

            endif cd2
                           
        enddo do2

        !Nullify auxiliar variables
        nullify (Prop, PropRef, BoundaryPoints2D, KFloorZ)

        !----------------------------------------------------------------------

    End Subroutine Prop_CyclicBoundary

    !--------------------------------------------------------------------------


    subroutine FinishAdvectionDiffusionIt()

        !External--------------------------------------------------------------

        !----------------------------------------------------------------------
        integer                             :: STAT_CALL

        Me%ExternalVar%DecayTime = null_real

        !Nullifies pointer which entered per argument
        nullify(Me%ExternalVar%PROP              )
        nullify(Me%ExternalVar%Wflux_X           )
        nullify(Me%ExternalVar%Wflux_Y           )
        nullify(Me%ExternalVar%Wflux_Z           )
        nullify(Me%ExternalVar%Visc_H            )
        nullify(Me%ExternalVar%Diff_V            )
        nullify(Me%ExternalVar%VolumeZOld        )
        nullify(Me%ExternalVar%VolumeZ           )
        nullify(Me%ExternalVar%OpenPoints3D      )
        nullify(Me%ExternalVar%ComputeFacesU3D   )
        nullify(Me%ExternalVar%ComputeFacesV3D   )
        nullify(Me%ExternalVar%ComputeFacesW3D   )

        if (associated(Me%ExternalVar%ReferenceProp)) nullify(Me%ExternalVar%ReferenceProp)
        if (associated(Me%ExternalVar%PROPOld))       nullify(Me%ExternalVar%PROPOld)
            

        !DUX, DVY, DZY, DZX
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%DUX,                                    &
                                 STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR01'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%DVY,                                    &
                                 STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR02'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%DZY,                                    &
                                 STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR03'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%DZX,                                    &
                                 STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR04'

        !BoundaryPoints2D
        call UngetHorizontalMap(Me%ObjHorizontalMap,                                    &
                                Me%ExternalVar%BoundaryPoints2D,                        &
                                STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR05'

        !AreaU
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExternalVar%AreaU, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR08'

        !AreaV
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExternalVar%AreaV, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR09'

        !DZZ
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExternalVar%DZZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR10'

        !DWZ
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExternalVar%DWZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR11'

        !KFloorZ
        call UnGetGeometry(Me%ObjGeometry,                                              &
                           Me%ExternalVar%KFloorZ, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)                                                      &
            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR12'

        !LandPoints3D
!        call UngetMap(Me%AppMap,                                      &
!                      Me%ExternalVar%LandPoints3D, STAT = STAT_CALL) 
!        if (STAT_CALL /= SUCCESS_)                                                     &
!            stop 'FinishAdvectionDiffusionIt - ModuleAdvectionDiffusion - ERR17'


        !Implicit-Explicit weight coeficients -> 1 = Implicit, 0 = Explicit
        Me%ExternalVar%ImpExp_AdvXX  = null_real         
        Me%ExternalVar%ImpExp_AdvYY  = null_real         
        Me%ExternalVar%ImpExp_AdvV   = null_real      
        Me%ExternalVar%ImpExp_DifH   = null_real    
        Me%ExternalVar%ImpExp_DifV   = null_real      


        !----------------------------------------------------------------------

    end subroutine FinishAdvectionDiffusionIt

    !--------------------------------------------------------------------------



    
    !--------------------------------------------------------------------------
    !MAnolo. Modification to Convert_VIsc_Dif_Vertical
    ! Computes diffusivities for every property from turbulent diffusivity and 
    ! SchmidtCoef_V and SchmidtBackground_V (read in WaterProperties)
    ! Diffusivity= SchmidtCoef_V * TurbulentDiffusivity + SchmidtBackground_V
    ! Also does some preparation work for subroutine vertical diffusion
    !--------------------------------------------------------------------------

    subroutine Convert_Dif_Vertical()

        !Local-----------------------------------------------------------------
        integer :: ILB, IUB 
        integer :: JLB, JUB 
        integer :: KLB, KUB
        integer :: i, j, k
        integer :: CHUNK
        
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
         

        !Begin-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "Convert_Dif_Vertical")
        
        !CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)
        CHUNK = CHUNK_J(JLB, JUB)

        !$OMP PARALLEL PRIVATE(i,j,k)

do1:   do k = KLB + 1, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2:   do j = JLB, JUB
do3:   do i = ILB, IUB     

cd1 :       if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == 1) then   

                Me%Diffusion_CoeficientZ(i,j,k) =                                       &                    
                              (Me%ExternalVar%SchmidtCoef_V                             &
                               * Me%ExternalVar%Diff_V(i,j,k)                           &
                               + Me%ExternalVar%SchmidtBackground_V)

            endif  cd1

        end do do3
        end do do2
        !$OMP END DO
        end do do1
        !!$OMP END PARALLEL


        if (Me%ExternalVar%NullDif)then

            do k = KLB + 1, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB     

                if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == 1) then   

                    if (Me%ExternalVar%Wflux_Z(i,j,k)==0.) then

                        Me%Diffusion_CoeficientZ(i, j, k) = 0.

                    end if

                end if

            enddo    
            enddo
            !$OMP END DO NOWAIT
            enddo   
         
        endif          
 
        !$OMP END PARALLEL
 
        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "Convert_Dif_Vertical")
 
    end subroutine Convert_Dif_Vertical


    !--------------------------------------------------------------------------
    !
    ! Converte viscosidades turbulentas (definidas nos centros das células Z) 
    ! em difusões turbulentas (definidas nas faces das células Z)
    !
    !--------------------------------------------------------------------------

    subroutine Convert_Visc_Dif_Horizontal()

        !External--------------------------------------------------------------


        !Local-----------------------------------------------------------------
        
        integer :: ILB, IUB 
        integer :: JLB, JUB 
        integer :: KLB, KUB
        integer :: i, j, k
        integer :: CHUNK
        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB
        JLB = Me%WorkSize%JLB 
        JUB = Me%WorkSize%JUB
        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "Convert_Visc_Dif_Horizontal")

        CHUNK = CHUNK_J(JLB, JUB)

        !$OMP PARALLEL PRIVATE(i,j,k)

        do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            if (Me%ExternalVar%ComputeFacesU3D(i, j, k) == 1) then
 
                Me%Diffusion_CoeficientX(i, j, k) =                  &
                        Me%ExternalVar%Schmidt_H                     &
                        * (Me%ExternalVar%Visc_H(i,j,  k)            &
                        *  Me%ExternalVar%DUX(i,j-1)                 &
                        +  Me%ExternalVar%Visc_H(i,j-1,k)            &
                        *  Me%ExternalVar%DUX(i,j  ))                &
                        / (Me%ExternalVar%DUX(i,j    )               &
                        +  Me%ExternalVar%DUX(i,j-1))           

            endif

        enddo
        enddo
        !$OMP END DO NOWAIT
        enddo

        do k = KLB, KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
        do j = JLB, JUB
        do i = ILB, IUB

            If (Me%ExternalVar%ComputeFacesV3D(i, j, k) == 1) then


                Me%Diffusion_CoeficientY(i,j,k)  =                    &
                            Me%ExternalVar%Schmidt_H                  &
                         * (Me%ExternalVar%Visc_H(i,  j,k)            &
                         *  Me%ExternalVar%DVY(i-1,j)                 &
                         +  Me%ExternalVar%Visc_H(i-1,j,k)            &
                         *  Me%ExternalVar%DVY(i,  j))                &
                         / (Me%ExternalVar%DVY(i,j)                   &
                         +  Me%ExternalVar%DVY(i-1,j))

            Endif

        
        enddo
        enddo
        !$OMP END DO
        enddo


nulldif:If (Me%ExternalVar%NullDif) Then

            do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%ComputeFacesU3D(i, j, k) == 1) then

                    if (Me%ExternalVar%Wflux_X(i,j,k)==0.) then

                        Me%Diffusion_CoeficientX(i, j, k)=0.

                    endif

                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            enddo


            do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = JLB, JUB
            do i = ILB, IUB

                if (Me%ExternalVar%ComputeFacesV3D(i, j, k) == 1) then

                    if (Me%ExternalVar%Wflux_Y(i,j,k)==0.) then 

                        Me%Diffusion_CoeficientY(i,j,k)=0.

                    endif

                endif

            enddo
            enddo
            !$OMP END DO NOWAIT
            enddo


        end if nulldif

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "Convert_Visc_Dif_Horizontal")


        !----------------------------------------------------------------------

    end subroutine Convert_Visc_Dif_Horizontal

    !--------------------------------------------------------------------------


    logical function SmallDepthCell (i, j) 

        !Arguments-------------------------------------------------------------
        integer                             :: i, j

        !Begin-----------------------------------------------------------------

        if (Me%ExternalVar%SmallDepthsPresent) then

            SmallDepthCell = Me%ExternalVar%SmallDepths(i,j)

        else

            SmallDepthCell = .false.

        endif
        
        !--------------------------------------------------------------------------

    end function SmallDepthCell 

    !--------------------------------------------------------------------------
    !
    !This routine works either explicitly or implicitly, according to ImpExp_DifV
    !   1 = Implicit, 0 = Explicit
    !
    !--------------------------------------------------------------------------

    Subroutine VerticalDiffusion()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------               
        real(8) :: AuxK, Aux1, Aux2

        integer :: i, j, k

        integer :: CHUNK

        !----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "VerticalDiffusion")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
      
        !$OMP PARALLEL PRIVATE(i,j,k,AuxK,Aux1,Aux2)

do2 :   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%ComputeFacesW3D(i, j, k  ) == 1 .and.                    &
                .not. SmallDepthCell (i, j))  then
            
                ! [m^2/s * m * m / m]
                AuxK  =  Me%Diffusion_CoeficientZ(i,j,k  )                              &
                       * Me%ExternalVar%DUX      (i,j    )                              &
                       * Me%ExternalVar%DVY      (i,j    )                              &
                       / Me%ExternalVar%DZZ      (i,j,k-1)

                ![m^3/s * s / m^3]
                Aux1 = AuxK * dble(Me%ExternalVar%DTProp) / Me%ExternalVar%VolumeZ(i, j, k-1) 
                Aux2 = AuxK * dble(Me%ExternalVar%DTProp) / Me%ExternalVar%VolumeZ(i, j, k  ) 

                Me%COEF3%E(i,j,k-1) = Me%COEF3%E(i,j,k-1) + Aux1 * Me%ExternalVar%ImpExp_DifV

                Me%COEF3%F(i,j,k-1) = Me%COEF3%F(i,j,k-1) - Aux1 * Me%ExternalVar%ImpExp_DifV

                Me%TICOEF3(i,j,k-1) = Me%TICOEF3(i,j,k-1) + Aux1 *                                  &
                                     (Me%ExternalVar%PROP(i,j,k)-Me%ExternalVar%PROP(i,j,k-1)) *    &
                                     (1. - Me%ExternalVar%ImpExp_DifV)

                Me%COEF3%D(i,j,k  ) = Me%COEF3%D(i,j,k  ) - Aux2 * Me%ExternalVar%ImpExp_DifV

                Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k  ) + Aux2 * Me%ExternalVar%ImpExp_DifV

                Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k  ) - Aux2 *                                  &
                                     (Me%ExternalVar%PROP(i,j,k)-Me%ExternalVar%PROP(i,j,k-1)) *    &
                                     (1. - Me%ExternalVar%ImpExp_DifV)

            endif

        end do do1
        end do do3
        !$OMP END DO
        end do do2

        !$OMP END PARALLEL

        if (Me%State%CellFluxes .and. Me%ExternalVar%ImpExp_DifV < 1.)                              &
            call CalcVerticalDifFlux(1. - Me%ExternalVar%ImpExp_DifV)


        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "VerticalDiffusion")

        !----------------------------------------------------------------------

    End Subroutine VerticalDiffusion                                                          

    !--------------------------------------------------------------------------

    subroutine VerticalAdvection()

        !Local-----------------------------------------------------------------               
        real(8)             :: AdvFluxZ, DT1, DT2
        integer             :: i, j, k  
        integer             :: CHUNK                           

        !----------------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "VerticalAdvection")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxZ,DT1,DT2)

st:     if (Me%State%VertAdv) then

            !CHUNK = CHUNK_J(Me%Size%JLB, Me%Size%JUB)
        
            !!$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)

j1:         do j = Me%WorkSize%JLB, Me%WorkSize%JUB
i1:         do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%OpenPoints3D (i,j,Me%WorkSize%KUB) == 1) then
                    call ComputeAdvection1D_V2( Me%WorkSize%KLB+1,                          &
                                                Me%WorkSize%KUB+1,                          &
                                                Me%ExternalVar%DTProp,                      &
                                                Me%ExternalVar%DWZ          (i,j,:),        &
                                                Me%ExternalVar%PROP         (i,j,:),        &
                                                Me%ExternalVar%Wflux_Z      (i,j,:),        &
                                                Me%ExternalVar%VolumeZ      (i,j,:),        & 
                                                Me%ExternalVar%OpenPoints3D (i,j,:),        &
                                                Me%COEF3_VertAdv%C_flux     (i,j,:),        &
                                                Me%COEF3_VertAdv%D_flux     (i,j,:),        &
                                                Me%COEF3_VertAdv%E_flux     (i,j,:),        &
                                                Me%COEF3_VertAdv%F_flux     (i,j,:),        &
                                                Me%ExternalVar%AdvMethodV,                  &
                                                Me%ExternalVar%TVDLimitationV,              &
                                                Me%ExternalVar%VolumeRelMax,                &
                                                Me%ExternalVar%Upwind2V)
                endif

            end do i1
            end do j1
            
            !$OMP END DO
            !!$OMP END PARALLEL

        endif st

        

cd6:    if (Me%ExternalVar%ImpExp_AdvV == ExplicitScheme)  then !ExplicitScheme = 0

dok3 :      do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj3 :      do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi3 :      do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%ComputeFacesW3D(i, j  , k) == 1) then
                
                AdvFluxZ =    (Me%COEF3_VertAdv%C_flux(i,   j,   k)                     &
                            *  Me%ExternalVar%PROP    (i,   j, k-2)                     &
                            +  Me%COEF3_VertAdv%D_flux(i,   j,   k)                     &
                            *  Me%ExternalVar%PROP    (i,   j, k-1)                     &
                            +  Me%COEF3_VertAdv%E_flux(i,   j,   k)                     &
                            *  Me%ExternalVar%PROP    (i,   j,   k)                     &
                            +  Me%COEF3_VertAdv%F_flux(i,   j,   k)                     &
                            *  Me%ExternalVar%PROP    (i,   j, k+1))

                Me%TICOEF3(i,j,k-1) = Me%TICOEF3(i,j,k-1) - AdvFluxZ *                  &
                                      Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j,k-1)
                Me%TICOEF3(i,j,k  ) = Me%TICOEF3(i,j,k  ) + AdvFluxZ *                  &
                                      Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j,k  )


            endif

            end do doi3
            end do doj3
            !$OMP END DO
            end do dok3

        else if (Me%ExternalVar%ImpExp_AdvV == ImplicitScheme) then cd6 !ImplicitScheme = 1

dok4 :      do k = Me%WorkSize%KLB, Me%WorkSize%KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj4 :      do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi4 :      do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                if (Me%ExternalVar%ComputeFacesW3D(i, j  , k) == 1) then

                    DT1 = Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j,k-1)
                    DT2 = Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j,k  )

                    Me%COEF3%D(i,j,k  ) = Me%COEF3%D(i,j,k  ) - Me%COEF3_VertAdv%D_flux(i,   j, k) * DT2
                    Me%COEF3%E(i,j,k  ) = Me%COEF3%E(i,j,k  ) - Me%COEF3_VertAdv%E_flux(i,   j, k) * DT2

                    Me%COEF3%E(i,j,k-1) = Me%COEF3%E(i,j,k-1) + Me%COEF3_VertAdv%D_flux(i,   j, k) * DT1
                    Me%COEF3%F(i,j,k-1) = Me%COEF3%F(i,j,k-1) + Me%COEF3_VertAdv%E_flux(i,   j, k) * DT1

                endif

            end do doi4
            end do doj4
            !$OMP END DO           
            end do dok4

        else cd6

            stop 'sub. VerticalAdvection - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        !$OMP END PARALLEL

        !Fluxes among cells
        if (Me%State%CellFluxes .and. Me%ExternalVar%ImpExp_AdvV < 1.)    &
            call CalcVerticalAdvFlux(1. - Me%ExternalVar%ImpExp_AdvV)

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "VerticalAdvection")

        !----------------------------------------------------------------------

    end subroutine VerticalAdvection

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    subroutine CalcVerticalAdvFlux(Weigth)


        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "CalcVerticalAdvFlux")

        CHUNK = CHUNK_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(k,j,i)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

        if (Me%ExternalVar%ComputeFacesW3D(i, j, k) == 1) then

            Me%Fluxes%AdvFluxZ(i, j, k) =                                               &
                          Me%Fluxes%AdvFluxZ(i, j, k)                                   &
                        + Weigth                                                        &
                        *(Me%COEF3_VertAdv%C_flux(i,j,k)                                &
                        * Me%ExternalVar%PROP(i, j, k-2)                                &
                        + Me%COEF3_VertAdv%D_flux(i,j,k)                                &
                        * Me%ExternalVar%PROP(i, j, k-1)                                &
                        + Me%COEF3_VertAdv%E_flux(i,j,k)                                &
                        * Me%ExternalVar%PROP(i, j, k )                                 &
                        + Me%COEF3_VertAdv%F_flux(i,j,k)                                &
                        * Me%ExternalVar%PROP(i, j, k+1))

        endif

        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "CalcVerticalAdvFlux")

        !----------------------------------------------------------------------

    end subroutine CalcVerticalAdvFlux

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalAdvFluxXX(Weigth)

        !External--------------------------------------------------------------
    
        real, intent(IN) :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------

        integer :: i,     j,     k  
        integer :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "CalcHorizontalAdvFluxXX")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j,k)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExternalVar%ComputeFacesU3D(i, j, k) == 1) then

            Me%Fluxes%AdvFluxX(i, j, k) =                                               &
                          Me%Fluxes%AdvFluxX(i, j, k)                                   &
                        + Weigth                                                        &
                        * (Me%COEF3_HorAdvXX%C_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i, j-2, k)                          &
                        +  Me%COEF3_HorAdvXX%D_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i, j-1, k)                          &
                        +  Me%COEF3_HorAdvXX%E_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i,   j, k)                          &
                        +  Me%COEF3_HorAdvXX%F_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i, j+1, k))

        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "CalcHorizontalAdvFluxXX")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalAdvFluxXX

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine CalcHorizontalAdvFluxYY(Weigth)

        !External--------------------------------------------------------------
    
        real, intent(IN) :: Weigth !Refers to the wigth of Implicit-Explicit calculations


        !Local-----------------------------------------------------------------

        integer :: i,     j,     k  
        integer :: CHUNK
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "CalcHorizontalAdvFluxYY")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k)
        
dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExternalVar%ComputeFacesV3D(i  , j, k) == 1) then
            Me%Fluxes%AdvFluxY(i, j, k) =                                               &
                          Me%Fluxes%AdvFluxY(i, j, k)                                   &
                        + Weigth                                                        &
                        * (Me%COEF3_HorAdvYY%C_flux(i  , j, k)                          &
                        *  Me%ExternalVar%PROP     (i-2, j, k)                          &
                        +  Me%COEF3_HorAdvYY%D_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i-1, j, k)                          &
                        +  Me%COEF3_HorAdvYY%E_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i,   j, k)                          &
                        +  Me%COEF3_HorAdvYY%F_flux(i,   j, k)                          &
                        *  Me%ExternalVar%PROP     (i+1, j, k))
        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "CalcHorizontalAdvFluxYY")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalAdvFluxYY

    !--------------------------------------------------------------------------

    subroutine CalcVerticalDifFlux(Weigth)

        !Arguments-------------------------------------------------------------
        real, intent(IN)                            :: Weigth !Refers to the wigth of Implicit-Explicit calculations

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "CalcVerticalDifFlux")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i, j, k)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExternalVar%ComputeFacesW3D(i, j, k  ) == 1) then
            Me%Fluxes%DifFluxZ(i, j, k) =                                               &
                          Me%Fluxes%DifFluxZ(i, j, k)                                   &
                        - Weigth                                                        &
                        * Me%Diffusion_CoeficientZ(i,j,k)                               &
                        * Me%ExternalVar%DUX(i,j)                                       &
                        * Me%ExternalVar%DVY(i,j)                                       &
                        / Me%ExternalVar%DZZ(i,j,k-1)                                   &
                        *(Me%ExternalVar%PROP(i, j, k  )                                &
                        - Me%ExternalVar%PROP(i, j, k-1))
        endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "CalcVerticalDifFlux")

        !----------------------------------------------------------------------

    end subroutine CalcVerticalDifFlux

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalDifFluxXX()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK
        
        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "CalcHorizontalDifFluxXX")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)
        
        !$OMP PARALLEL PRIVATE(i,j,k)
        
dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
            if (Me%ExternalVar%ComputeFacesU3D(i, j  , k) == 1) then                           
                Me%Fluxes%DifFluxX(i, j, k) =                                           &
                              Me%Fluxes%DifFluxX      (i,  j, k)                        &
                            - Me%Diffusion_CoeficientX(i,  j, k)                        &
                            * Me%ExternalVar%AreaU    (i,  j, k)                        &
                            / Me%ExternalVar%DZX      (i,j-1   )                        &
                            *(Me%ExternalVar%PROP     (i,  j, k)                        &
                            - Me%ExternalVar%PROP     (i,j-1, k))

            endif
        end do doi1
        end do doj1
        !$OMP END DO NOWAIT
        end do dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "CalcHorizontalDifFluxXX")

        !----------------------------------------------------------------------

    end subroutine CalcHorizontalDifFluxXX

    !--------------------------------------------------------------------------

    subroutine CalcHorizontalDifFluxYY()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "CalcHorizontalDifFluxYY")

        CHUNK = CHUNK_J(Me%WorkSize%JLB, Me%WorkSize%JUB)

        !$OMP PARALLEL PRIVATE(i,j,k)

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB
        if (Me%ExternalVar%ComputeFacesV3D(i  , j, k) == 1) then
            Me%Fluxes%DifFluxY(i, j, k) =                                               &
                          Me%Fluxes%DifFluxY      (i  , j, k)                           &
                        - Me%Diffusion_CoeficientY(i  , j, k)                           &
                        * Me%ExternalVar%AreaV    (i  , j, k)                           &
                        / Me%ExternalVar%DZY      (i-1, j   )                           &
                        *(Me%ExternalVar%PROP     (i  , j, k)                           &
                        - Me%ExternalVar%PROP     (i-1, j, k))
        endif
        end do doi1
        end do doj1
        !$OMP END DO
        end do dok1
        
        !$OMP END PARALLEL
        
        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "CalcHorizontalDifFluxYY")
        
        !----------------------------------------------------------------------

    end subroutine CalcHorizontalDifFluxYY

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------
    !
    !This term appears due to (Ct+Dt.Vt+Dt - Ct.Vt)/Dt -> Ct.Vt/Vt+Dt
    !
    !--------------------------------------------------------------------------

    subroutine VolumeVariation()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8)                                     :: DT_V
        integer                                     :: i, j, k
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "VolumeVariation")

        CHUNK = Chunk_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j,k,DT_V) 

dok1:   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj1:   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
doi1:   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

cd1:       if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then

                Me%TICOEF3(i,j,k) =                                   &
                                   Me%TICOEF3(i,j,k)                  &
                                 + Me%ExternalVar%PROP(i,j,k)         &
                                 * (Me%ExternalVar%VolumeZOld(i,j,k)  &
                                 /  Me%ExternalVar%VolumeZ(i,j,k)   )

                

                ! if WFLUX_Z(i,j,KUB+1) different from zero corrects the mass fluxes 
                ! to maintain the global mass
                !In the surface layer the follow is true:
                ! Vreal = VolumeZ + dt*FluxZ(kub+1)
                ! d(VC)/dt = (C(t+dt)*Vreal(t+dt) - C(t)*V(t))/ dt = (C(t+dt) * (VolumeZ + dt *FluxZ(kub+1)) -C(t)*V(t))/ dt 
                ! C(t+dt)* (1 + dt * FluxZ(kub+1) / VolumeZ)  = C(t)*V(t) / VolumeZ + Transporte
                if (k == Me%WorkSize%KUB)  then
                    DT_V              = dble(Me%ExternalVar%DTProp) / Me%ExternalVar%VolumeZ(i,j,k)
                    Me%COEF3%E(i,j,k) = Me%COEF3%E(i,j,k) + DT_V * Me%ExternalVar%Wflux_Z(i,j,Me%WorkSize%KUB+1)
                endif
            
            else cd1 !This is important for the tidal flats points to maintain 
                     !their concentration when they are not covered 

                Me%TICOEF3(i,j,k) = Me%ExternalVar%PROP(i,j,k)

            endif cd1

         enddo doi1
         enddo doj1
         !$OMP END DO     
         enddo dok1

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "VolumeVariation")

        !----------------------------------------------------------------------
    
    end Subroutine VolumeVariation    

    !--------------------------------------------------------------------------

    subroutine Discharges ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real,   dimension(:,:), pointer             :: WaterColumnZ
        real(8)                                     :: DT_V
        real                                        :: Flow
        integer                                     :: i, j, k, kd, n, dis, nc, kmin, kmax
        integer                                     :: STAT_CALL

        !----------------------------------------------------------------------

        !WaterColumnZ
        call GetGeometryWaterColumn(Me%ObjGeometry, WaterColumn = WaterColumnZ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'Discharges - ModuleAdvectionDiffusion - ERR10'


        n = 0
        
ddis:   do dis = 1, Me%ExternalVar%DischNumber  

            if (Me%ExternalVar%IgnoreDisch(dis)) cycle 
            
dnc:        do   nc = 1, Me%ExternalVar%DischnCells(dis)     

                n = n + 1
            
                i  = Me%ExternalVar%DischI(n)
                j  = Me%ExternalVar%DischJ(n)
                kd = Me%ExternalVar%DischK(n)

                if (Me%ExternalVar%DischVert(dis) == DischUniform_) then

                    kmin = Me%ExternalVar%KFloorZ(i,j)
                    kmax = Me%WorkSize%KUB

                else
        
                    kmin = kd; kmax = kd

                endif

dk:             do k=kmin, kmax

                    DT_V  = dble(Me%ExternalVar%DTProp) / Me%ExternalVar%VolumeZ(i,j,k)

                    if (Me%ExternalVar%DischVert(dis)  == DischUniform_) then

                        Flow = Me%ExternalVar%DischFlow(n) *                            &
                               Me%ExternalVar%DWZ(i,j,k) / WaterColumnZ(i,j)

                    else

                        Flow = Me%ExternalVar%DischFlow(n)

                    endif

cd1:                if (Me%ExternalVar%OpenPoints3D(i, j, k) == OpenPoint) then
            
fl:                     if (Flow > 0.) then
                
                            Me%TICOEF3(i,j,k) = Me%TICOEF3(i,j,k)               +               &
                                                Flow * DT_V * Me%ExternalVar%DischConc(n)
                        else fl
                
                            Me%COEF3%E(i,j,k) = Me%COEF3%E(i,j,k)               -               &
                                                Flow * DT_V 
                        endif fl

                    else
                    !This for the case that tonly discharges changes the volume
                    ! In this case the cell is water point and is not open point
                    !In this case the subroutine VolumeVariation do not change the cell concentration (no dilution)
                        if (Flow > 0) then

                            Me%TICOEF3(i,j,k) = Me%ExternalVar%PROP(i, j, k)     *          &
                                                Me%ExternalVar%VolumeZOld(i,j,k) /          &
                                                Me%ExternalVar%VolumeZ   (i,j,k) +          &
                                                Flow * DT_V * Me%ExternalVar%DischConc(n)
                        endif
                    
                    endif cd1

                enddo dk

            enddo dnc
        
        enddo ddis


        !WaterColumnZ
        call UnGetGeometry(Me%ObjGeometry, WaterColumnZ, STAT = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)  stop 'Discharges - ModuleAdvectionDiffusion - ERR20'


    end subroutine Discharges 

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvection(ImpExp_AdvXX, ImpExp_AdvYY)

        !Arguments-------------------------------------------------------------
        real                                :: ImpExp_AdvXX, ImpExp_AdvYY

        !Local-----------------------------------------------------------------
        integer                             :: di,    dj    
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalAdvection")

        ILBWS = Me%WorkSize%ILB
        IUBWS = Me%WorkSize%IUB

        JLBWS = Me%WorkSize%JLB
        JUBWS = Me%WorkSize%JUB

        KLBWS = Me%WorkSize%KLB
        KUBWS = Me%WorkSize%KUB

        ILB   = Me%Size%ILB
        IUB   = Me%Size%IUB

        JLB   = Me%Size%JLB 
        JUB   = Me%Size%JUB

        KLB   = Me%Size%KLB
        KUB   = Me%Size%KUB

        
        call HorizontalAdvectionXX(ImpExp_AdvXX)

        if (.not. Me%XZFlow) call HorizontalAdvectionYY(ImpExp_AdvYY)

cd1:    if (ImpExp_AdvYY == ImplicitScheme .or. ImpExp_AdvXX == ImplicitScheme) then 


cd3D:       if (KUBWS > 1) then

                !if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorAdvection-THOMAS3D")


cd2:            if (ImpExp_AdvXX == ImplicitScheme) then 

                    di = 0
                    dj = 1

                    !griflet: old call
                    !call THOMAS_3D(ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS, di, dj,    &
                    !     Me%COEF3%D,                                                    &
                    !     Me%COEF3%E,                                                    &
                    !     Me%COEF3%F,                                                    &
                    !     Me%TICOEF3,                                                    &
                    !     Me%ExternalVar%PROP,                                           &
                    !     Me%VECG,                                                       &
                    !     Me%VECW)     
                    !griflet: new call 
                    call THOMAS_3D(ILBWS, IUBWS, JLBWS, JUBWS, KLBWS, KUBWS, di, dj,    &
                         Me%THOMAS,                                                     &
                         Me%ExternalVar%PROP)
                
                else if (ImpExp_AdvYY == ImplicitScheme) then cd2

                    di = 1
                    dj = 0

                    !griflet: old call
                    !call THOMAS_3D(JLBWS, JUBWS, ILBWS, IUBWS, KLBWS, KUBWS, di, dj,    &
                    !     Me%COEF3%D,                                                    &
                    !     Me%COEF3%E,                                                    &
                    !     Me%COEF3%F,                                                    &
                    !     Me%TICOEF3,                                                    &
                    !     Me%ExternalVar%PROP,                                           &
                    !     Me%VECG,                                                       &
                    !     Me%VECW)      
                    !griflet: new call
                    call THOMAS_3D(JLBWS, JUBWS, ILBWS, IUBWS, KLBWS, KUBWS, di, dj,    &
                         Me%THOMAS,                                                     &
                         Me%ExternalVar%PROP)      

                endif cd2

                
                call SetMatrixValue (Me%COEF3%D, Me%Size, 0.0)
                call SetMatrixValue (Me%COEF3%E, Me%Size, dble(1.0))
                call SetMatrixValue (Me%COEF3%F, Me%Size, 0.0)
                call SetMatrixValue (Me%TICOEF3, Me%Size, Me%ExternalVar%PROP)
!dok6 :          do k = KLB, KUB 
!doj6 :          do j = JLB, JUB
!doi6 :          do i = ILB, IUB
!                    Me%COEF3%D(i,j,k) = 0.0                                            
!                    Me%COEF3%E(i,j,k) = 1.0  
!                    Me%COEF3%F(i,j,k) = 0.0   
!                    Me%TICOEF3(i,j,k) = Me%ExternalVar%PROP(i,j,k)
!                end do doi6
!                end do doj6
!                end do dok6

                !if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorAdvection-THOMAS3D")

            endif cd3D
     

        endif cd1

        
        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalAdvection")

        !----------------------------------------------------------------------

    end subroutine HorizontalAdvection

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvectionXX(ImpExp_AdvXX)

        !External--------------------------------------------------------------

        real                                :: ImpExp_AdvXX

        !Local-----------------------------------------------------------------               

        real(8) :: AdvFluxX, DT1, DT2

        integer :: i,     j,     k                             
        integer :: ILB, IUB, JLB, JUB, KLB, KUB
        integer :: CHUNK
        !----------------------------------------------------------------------


        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalAdvectionXX")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB

        CHUNK = CHUNK_I(ILB, IUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxX,DT2,DT1)

st:     if (Me%State%HorAdv) then
            
            !CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)

            !!$OMP PARALLEL SHARED(CHUNK) PRIVATE(I,K)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)

k1:         do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
i1:         do i = ILB, IUB


                call ComputeAdvection1D_V2(JLB+1, JUB+1, Me%ExternalVar%DTProp,         &
                                        Me%ExternalVar%DUX          (i,:),              &
                                        Me%ExternalVar%PROP         (i,:,k),            &
                                        Me%ExternalVar%Wflux_X      (i,:,k),            &
                                        Me%ExternalVar%VolumeZ      (i,:,k),            & 
                                        Me%ExternalVar%OpenPoints3D (i,:,k),            &
                                        Me%COEF3_HorAdvXX%C_flux    (i,:,k),            &
                                        Me%COEF3_HorAdvXX%D_flux    (i,:,k),            &
                                        Me%COEF3_HorAdvXX%E_flux    (i,:,k),            &
                                        Me%COEF3_HorAdvXX%F_flux    (i,:,k),            &
                                        Me%ExternalVar%AdvMethodH,                      &
                                        Me%ExternalVar%TVDLimitationH,                  &
                                        Me%ExternalVar%VolumeRelMax,                    &
                                        Me%ExternalVar%Upwind2H)

            end do i1
            !$OMP END DO
            end do k1
        
        !!$OMP END PARALLEL

        endif st

        CHUNK = CHUNK_K(KLB, KUB)

cd6:    if (ImpExp_AdvXX == ExplicitScheme)  then !ExplicitScheme = 0

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
dok3 :      do k = KLB, KUB
doj3 :      do j = JLB, JUB
doi3 :      do i = ILB, IUB

            if (Me%ExternalVar%ComputeFacesU3D(i, j  , k) == 1) then

                AdvFluxX =    (Me%COEF3_HorAdvXX%C_flux(i,   j, k)                          &
                            *  Me%ExternalVar%PROP     (i, j-2, k)                          &
                            +  Me%COEF3_HorAdvXX%D_flux(i,   j, k)                          &
                            *  Me%ExternalVar%PROP     (i, j-1, k)                          &
                            +  Me%COEF3_HorAdvXX%E_flux(i,   j, k)                          &
                            *  Me%ExternalVar%PROP     (i,   j, k)                          &
                            +  Me%COEF3_HorAdvXX%F_flux(i,   j, k)                          &
                            *  Me%ExternalVar%PROP     (i, j+1, k))

                Me%TICOEF3(i,j  ,k) = Me%TICOEF3(i,j  ,k) + AdvFluxX * Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j  ,k)
                Me%TICOEF3(i,j-1,k) = Me%TICOEF3(i,j-1,k) - AdvFluxX * Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j-1,k)

            endif

            end do doi3
            end do doj3
            end do dok3
            !$OMP END DO NOWAIT 

        else if (ImpExp_AdvXX == ImplicitScheme) then cd6 !ImplicitScheme = 1

            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
dok4 :      do k = KLB, KUB
doj4 :      do j = JLB, JUB
doi4 :      do i = ILB, IUB


            if (Me%ExternalVar%ComputeFacesU3D(i, j  , k) == 1) then

                DT2 = Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j  ,k)
                DT1 = Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i,j-1,k)

                Me%COEF3%D(i,j  ,k) = Me%COEF3%D(i,j  ,k) - Me%COEF3_HorAdvXX%D_flux(i,   j, k) * DT2
                Me%COEF3%E(i,j  ,k) = Me%COEF3%E(i,j  ,k) - Me%COEF3_HorAdvXX%E_flux(i,   j, k) * DT2

                Me%COEF3%E(i,j-1,k) = Me%COEF3%E(i,j-1,k) + Me%COEF3_HorAdvXX%D_flux(i,   j, k) * DT1
                Me%COEF3%F(i,j-1,k) = Me%COEF3%F(i,j-1,k) + Me%COEF3_HorAdvXX%E_flux(i,   j, k) * DT1

!The module advection is only able to compute tridiagonal linear systems 
!
!                    Me%COEF3%C(i,j,k) =  not computed
!                    Me%COEF3%G(i,j,k) =  not computed

!                      The aux variables
!                      Me%COEF3_VertAdv%C_flux and are not used Me%COEF3_VertAdv%G_flux


            endif


            end do doi4
            end do doj4
            end do dok4
            !$OMP END DO NOWAIT 

        else cd6

            stop 'sub. HorizontalAdvectionXX - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalAdvectionXX")

        if (Me%State%CellFluxes .and. ImpExp_AdvXX == ExplicitScheme)                   &
            call CalcHorizontalAdvFluxXX(1. - ImpExp_AdvXX)

        !----------------------------------------------------------------------

    end subroutine HorizontalAdvectionXX

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine HorizontalAdvectionYY(ImpExp_AdvYY)

        !Arguments--------------------------------------------------------------
        real                                :: ImpExp_AdvYY

        !Local-----------------------------------------------------------------               
        real(8)                             :: AdvFluxY, DT1, DT2
        integer                             :: i,     j,     k                             
        integer                             :: ILB, IUB, JLB, JUB, KLB, KUB
        integer                             :: CHUNK

        !----------------------------------------------------------------------

        !if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalAdvectionYY")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KLB = Me%WorkSize%KLB
        KUB = Me%WorkSize%KUB
        
        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalAdvectionYY")
        
        !CHUNK = CHUNK_K(Me%Size%KLB, Me%Size%KUB)
        CHUNK = CHUNK_J(JLB, JUB)
      
        !$OMP PARALLEL PRIVATE(i,j,k,AdvFluxY,DT2,DT1)


st:     if (Me%State%HorAdv) then
        
            !!$OMP PARALLEL SHARED(CHUNK) PRIVATE(J,K)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)

k1:         do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
j1:         do j = JLB, JUB


                call ComputeAdvection1D_V2(ILB+1, IUB+1, Me%ExternalVar%DTProp,         &
                                        Me%ExternalVar%DVY          (:,j),              &
                                        Me%ExternalVar%PROP         (:,j,k),            &
                                        Me%ExternalVar%Wflux_Y      (:,j,k),            &
                                        Me%ExternalVar%VolumeZ      (:,j,k),            & 
                                        Me%ExternalVar%OpenPoints3D (:,j,k),            &
                                        Me%COEF3_HorAdvYY%C_flux    (:,j,k),            &
                                        Me%COEF3_HorAdvYY%D_flux    (:,j,k),            &
                                        Me%COEF3_HorAdvYY%E_flux    (:,j,k),            &
                                        Me%COEF3_HorAdvYY%F_flux    (:,j,k),            &
                                        Me%ExternalVar%AdvMethodH,                      &
                                        Me%ExternalVar%TVDLimitationH,                  &
                                        Me%ExternalVar%VolumeRelMax,                    &
                                        Me%ExternalVar%Upwind2H)

            end do j1
            !$OMP END DO
            end do k1
            
            !!$OMP END DO NOWAIT
            !!$OMP END PARALLEL

        endif st


cd6:    if (ImpExp_AdvYY == ExplicitScheme)  then !ExplicitScheme = 0

dok3 :      do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj3 :      do j = JLB, JUB
doi3 :      do i = ILB, IUB

            if (Me%ExternalVar%ComputeFacesV3D(i, j  , k) == 1) then

                AdvFluxY =    (Me%COEF3_HorAdvYY%C_flux(  i, j, k)                          &
                            *  Me%ExternalVar%PROP     (i-2, j, k)                          &
                            +  Me%COEF3_HorAdvYY%D_flux(  i, j, k)                          &
                            *  Me%ExternalVar%PROP     (i-1, j, k)                          &
                            +  Me%COEF3_HorAdvYY%E_flux(  i, j, k)                          &
                            *  Me%ExternalVar%PROP     (  i, j, k)                          &
                            +  Me%COEF3_HorAdvYY%F_flux(  i, j, k)                          &
                            *  Me%ExternalVar%PROP     (i+1, j, k))

                Me%TICOEF3(i  ,j,k) = Me%TICOEF3(i  ,j,k) + AdvFluxY * Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i  ,j,k)
                Me%TICOEF3(i-1,j,k) = Me%TICOEF3(i-1,j,k) - AdvFluxY * Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i-1,j,k)

            endif

            end do doi3
            end do doj3
            !$OMP END DO NOWAIT
            end do dok3

        else if (ImpExp_AdvYY == ImplicitScheme) then cd6 !ImplicitScheme = 1

dok4 :      do k = KLB, KUB
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj4 :      do j = JLB, JUB
doi4 :      do i = ILB, IUB


            if (Me%ExternalVar%ComputeFacesV3D(i, j  , k) == 1) then

                DT2 = Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i  ,j  ,k)
                DT1 = Me%ExternalVar%DTProp / Me%ExternalVar%VolumeZ(i-1,j  ,k)

                Me%COEF3%D(i,j  ,k) = Me%COEF3%D(i,j  ,k) - Me%COEF3_HorAdvYY%D_flux(i,   j, k) * DT2
                Me%COEF3%E(i,j  ,k) = Me%COEF3%E(i,j  ,k) - Me%COEF3_HorAdvYY%E_flux(i,   j, k) * DT2

                Me%COEF3%E(i-1,j,k) = Me%COEF3%E(i-1,j,k) + Me%COEF3_HorAdvYY%D_flux(i,   j, k) * DT1
                Me%COEF3%F(i-1,j,k) = Me%COEF3%F(i-1,j,k) + Me%COEF3_HorAdvYY%E_flux(i,   j, k) * DT1

!The module advection is only able to compute tridiagonal linear systems 
!
!                    Me%COEF3%C(i,j,k) =  not computed
!                    Me%COEF3%G(i,j,k) =  not computed

!                      The aux variables
!                      Me%COEF3_VertAdv%C_flux and are not used Me%COEF3_VertAdv%G_flux


            endif


            end do doi4
            end do doj4
            !$OMP END DO NOWAIT
            end do dok4

        else cd6

            stop 'sub. HorizontalAdvectionYY - ModuleAdvectionDiffusion - ERR01'
        
        endif cd6

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalAdvectionYY")

        if (Me%State%CellFluxes .and. ImpExp_AdvYY == ExplicitScheme)                   &
            call CalcHorizontalAdvFluxYY(1. - ImpExp_AdvYY)
        
        !if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalAdvectionYY")

        !----------------------------------------------------------------------

    end subroutine HorizontalAdvectionYY

    
    !--------------------------------------------------------------------------


    subroutine HorizontalDiffusion()

        !External--------------------------------------------------------------

        !----------------------------------------------------------------------
        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalDiffusion")

        call HorizontalDiffusionXX()
        
        if (.not. Me%XZFlow) call HorizontalDiffusionYY()
                
        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalDiffusion")

        !----------------------------------------------------------------------

    end subroutine HorizontalDiffusion

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine HorizontalDiffusionXX()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8)                                     :: DTPropDouble 
        real(8)                                     :: AuxJ
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        DTPropDouble = dble(Me%ExternalVar%DTProp)

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalDiffusionXX")

        CHUNK = Chunk_K(Me%WorkSize%KLB,Me%WorkSize%KUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AuxJ) 
 
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
do2 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%ComputeFacesU3D(i, j, k) == 1) then

                !DTV1 = DTPropDouble / Me%ExternalVar%VolumeZ(i, j-1, k)
                !DTV2 = DTPropDouble / Me%ExternalVar%VolumeZ(i, j  , k)
                        
                AuxJ = Me%Diffusion_CoeficientX (i,j  ,k)                               &
                       * Me%ExternalVar%AreaU   (i,j  ,k)                               &
                       / Me%ExternalVar%DZX     (i,j-1  )                    

                Me%TICOEF3(i,j-1,k) = Me%TICOEF3(i,j-1,k) + AuxJ * DTPropDouble /       &
                                      Me%ExternalVar%VolumeZ(i, j-1, k) *               &
                                     (Me%ExternalVar%PROP(i,j,k) - Me%ExternalVar%PROP(i,j-1,k))


                Me%TICOEF3(i,j  ,k) = Me%TICOEF3(i,j  ,k) - AuxJ * DTPropDouble /       &
                                      Me%ExternalVar%VolumeZ(i, j  , k) *               &
                                     (Me%ExternalVar%PROP(i,j,k) - Me%ExternalVar%PROP(i,j-1,k))


            endif

        end do do1
        end do do2
        end do do3
        !$OMP END DO
            
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalDiffusionXX")
    
        if (Me%State%CellFluxes) call CalcHorizontalDifFluxXX()

        !----------------------------------------------------------------------

    end subroutine HorizontalDiffusionXX

    !--------------------------------------------------------------------------

    subroutine HorizontalDiffusionYY()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        real(8)                                     :: DTPropDouble 
        real(8)                                     :: AuxI
        integer                                     :: i, j, k  
        integer                                     :: CHUNK

        !----------------------------------------------------------------------

        DTPropDouble = dble(Me%ExternalVar%DTProp) 

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "HorizontalDiffusionYY")

        CHUNK = Chunk_J(Me%WorkSize%JLB,Me%WorkSize%JUB)
        !$OMP PARALLEL PRIVATE(i,j,k,AuxI) 

do3 :   do k = Me%WorkSize%KLB, Me%WorkSize%KUB
        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do2 :   do j = Me%WorkSize%JLB, Me%WorkSize%JUB
do1 :   do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            if (Me%ExternalVar%ComputeFacesV3D(i  , j, k) == 1) then

                !DTV1 = DTPropDouble / Me%ExternalVar%VolumeZ(i-1, j, k)
                !DTV2 = DTPropDouble / Me%ExternalVar%VolumeZ(i  , j, k)
                        
                AuxI = Me%Diffusion_CoeficientY (i  ,j,k)                               &
                       * Me%ExternalVar%AreaV   (i  ,j,k)                               &
                       / Me%ExternalVar%DZY     (i-1,j  )                    

                Me%TICOEF3(i-1,j,k) = Me%TICOEF3(i-1,j,k) + AuxI * DTPropDouble /       &
                                      Me%ExternalVar%VolumeZ(i-1, j, k) *               &
                                     (Me%ExternalVar%PROP(i,j,k) - Me%ExternalVar%PROP(i-1,j,k))


                Me%TICOEF3(i,j  ,k) = Me%TICOEF3(i,j  ,k) - AuxI * DTPropDouble /       &
                                      Me%ExternalVar%VolumeZ(i  , j, k) *               &
                                     (Me%ExternalVar%PROP(i,j,k) - Me%ExternalVar%PROP(i-1,j,k))
            endif
        end do do1
        end do do2
        !$OMP END DO
        end do do3

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "HorizontalDiffusionYY")

        if (Me%State%CellFluxes) call CalcHorizontalDifFluxYY()

        !----------------------------------------------------------------------

    end subroutine HorizontalDiffusionYY

    !--------------------------------------------------------------------------

    subroutine OpenBoundaryCondition()

        !Arguments-------------------------------------------------------------

        !local-----------------------------------------------------------------

        integer :: i, j, k  
        integer :: ILB, IUB
        integer :: JLB, JUB
        integer :: KLB, KUB
        integer :: ILBSize, IUBSize, JLBSize, JUBSize

        integer :: BoundaryCondition

        real    :: TdecAux
        real    :: ExteriorProp, ImpExp_DifV
        real(8) :: DT_V, DTPropDouble
        real    :: Atotal, A1, A2, A3, A4
        real    :: P1, P2, P3, P4, InteriorProp, VelBound
        logical :: EastNorthBoundary
        integer :: di, dj, iext, jext
        real    :: LimitMax, BoundaryProp
        logical :: Found
        
        !integer :: CHUNK

        !----------------------------------------------------------------------

        !if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "OpenBoundaryCondition")

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        ILBSize = Me%WorkSize%ILB
        IUBSize = Me%WorkSize%IUB

        JLBSize = Me%WorkSize%JLB
        JUBSize = Me%WorkSize%JUB


        KUB = Me%WorkSize%KUB

        BoundaryCondition = Me%ExternalVar%BoundaryCondition

        ImpExp_DifV       = Me%ExternalVar%ImpExp_DifV

        DTPropDouble      = dble(Me%ExternalVar%DTProp)

        TdecAux           =  1.0 / (1.0 + Me%ExternalVar%DecayTime &
                             / Me%ExternalVar%DTProp)    

        if (BoundaryCondition == MassConservation_ .or. BoundaryCondition == Orlanski_ .or. &
            BoundaryCondition == MassConservNullGrad_ )                                     &
            call FluxAtOpenBoundary()


        if (BoundaryCondition == Orlanski_ )                                             &
            Me%ExternalVar%PROPOld(:,:,:) = Me%ExternalVar%PROP(:,:,:)

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "OpenBoundaryCondition")

        !CHUNK = CHUNK_J(JLB,JUB)

        !!$OMP PARALLEL PRIVATE(i,j,k,KLB,DT_V,A1,A2,A3,A4,P1,P2,P3,P4,Atotal,InteriorProp) &
        !!$OMP PRIVATE(ExteriorProp,VelBound,di,dj,iext,jext,LimitMax,Found,Me%ExternalVar%PROP)

        !!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do3 :   do  j = JLB, JUB
do1 :   do  i = ILB, IUB

cd9 :       if (Me%ExternalVar%BoundaryPoints2D(i,j) == 1) then 

                KLB = ABS(Me%ExternalVar%KFloorZ(i,j))


do2 :           do  k = KLB, KUB
                                         
cd3:                if (Me%ExternalVar%OpenPoints3D(i, j, k) == 1) then
                
                        DT_V = DTPropDouble / Me%ExternalVar%VolumeZ(i,j,k)    


                        if (BoundaryCondition == ImposedValue_ .or.                      &
                            BoundaryCondition == SubModel_) then

                            !Valor da propriedade a entrar no dominio em funcao do                      
                            !  valor imposto no exterior e do tempo de decaimento
                            !A expressao do tempo de decaimento e' resolvida para 
                            !  trás no tempo por razoes de estabilidade (Tese de doutoramento 
                            !  do Aires pág. 129

                            !Neighbours points of the "boundary point" that are interior.                         
                            A1 = Me%ExternalVar%OpenPoints3D(i + 1, j    , k) * &
                                 (1 - Me%ExternalVar%BoundaryPoints2D(i + 1, j    ))

                            A2 = Me%ExternalVar%OpenPoints3D(i - 1, j    , k) * &
                                 (1 - Me%ExternalVar%BoundaryPoints2D(i - 1, j    ))

                            A3 = Me%ExternalVar%OpenPoints3D(i    , j + 1, k) * &
                                 (1 - Me%ExternalVar%BoundaryPoints2D(i    , j + 1))

                            A4 = Me%ExternalVar%OpenPoints3D(i    , j - 1, k) * &
                                 (1 - Me%ExternalVar%BoundaryPoints2D(i    , j - 1))

                            !Property values of Neighbours points
                            P1 = Me%ExternalVar%PROP        (i + 1, j    , k)
                            P2 = Me%ExternalVar%PROP        (i - 1, j    , k)
                            P3 = Me%ExternalVar%PROP        (i    , j + 1, k)
                            P4 = Me%ExternalVar%PROP        (i    , j - 1, k)

                            Atotal = A1 + A2 + A3 + A4

                            if (Atotal > 0) then

                                InteriorProp = (P1 * A1 + P2 * A2 + P3 * A3 + P4 * A4) / Atotal

                                !Property imposed in the exterior using a decay time 
                                ExteriorProp = InteriorProp * (1.0 - TdecAux)            &
                                           +   Me%ExternalVar%ReferenceProp(i,j,k) *  TdecAux
                            else 

                                ExteriorProp = Me%ExternalVar%ReferenceProp(i,j,k)

                            endif

                        endif

                                               
cd1:                    if (BoundaryCondition == MassConservation_ .or.                  &
                            BoundaryCondition == Orlanski_         .or.                  &
                            BoundaryCondition == MassConservNullGrad_ ) then 

                            if (BoundaryCondition == Orlanski_) then


                                ![]      = [m^3/s] * [s/m^3]
                                VelBound = Me%WaterFluxOBoundary(i,j,k) * DT_V 

                                if      (i == ILB .or. i == IUB) then

                                    di = 1 
                                    dj = 0

                                else if (j == JLB .or. j == JUB) then

                                    di = 0 
                                    dj = 1

                                else

                                    Stop 'Orlanski Advection 2'

                                endif

                                if      (i == ILB .or. j == JLB) then

                                    EastNorthBoundary = .false.

                                    iext = i - di
                                    jext = j - dj

                                else if (i == IUB .or. j == JUB) then

                                    EastNorthBoundary = .true.

                                    iext = i + di
                                    jext = j + dj

                                else

                                    Stop 'Orlanski Advection 1'

                                endif

                                LimitMax = MaxInternalCelerity * Me%ExternalVar%DTProp / &
                                            (Me%ExternalVar%DUX(i, j) + &
                                             Me%ExternalVar%DVY(i, j)) * 2


                                call OrlanskiCelerity2D(NewField          = Me%ExternalVar%PROP,          &
                                                        OldField          = Me%ExternalVar%PROPOld,       &
                                                        ComputePoints     = Me%ExternalVar%OpenPoints3D,  &
                                                        ReferenceField    = Me%ExternalVar%ReferenceProp, &
                                                        Imin              = ILBSize,                                         &
                                                        Imax              = IUBSize,                                         &
                                                        Jmin              = JLBSize,                                         &
                                                        Jmax              = JUBSize,                                         &
                                                        di                = di,                                              &
                                                        dj                = dj,                                              &
                                                        i                 = iext,                                            &
                                                        j                 = jext,                                            &
                                                        k                 =  k,                                              &
                                                        LimitMax          = LimitMax,                                        &
                                                        EastNorthBoundary = EastNorthBoundary,                               &
                                                        DT                = Me%ExternalVar%DTProp,        &
                                                        FlowVelX          = VelBound,                                        &
                                                        NewValue          = ExteriorProp)

!                                call ExteriorValuesRadiation(Me, VelBound, i, j, k, BP2D)

                            else if (BoundaryCondition == MassConservation_ ) then

                                !In this option the interior point is also the boundary point
                                InteriorProp = Me%ExternalVar%PROP(i, j, k)


                                ! dP/dt = (Preference -P ) / Tdec
                                ! Tdec - infinity => TdecAux = 0 => ExteriorProp = P
                                ! Tdec - zero     => TdecAux = 1 => ExteriorProp = Preference

                                !Property imposed in the exterior using a decay time 
                                ExteriorProp = InteriorProp * (1.0 - TdecAux)            &
                                           +   Me%ExternalVar%ReferenceProp(i,j,k) *  TdecAux

                            endif


                            ! This method don't conserve mass if the Diffusion fluxe between the 
                            ! frontier and the interior is of the same magnitude of the advective fluxe
                        
cd2 :                       if (Me%WaterFluxOBoundary(i,j,k) < 0.0) Then ! water is flowing in    


                                if (.not. BoundaryCondition == MassConservNullGrad_) then
                                !Mass is being Added to the domain by default the flow 
                                ! is negative when is going in and this process is explicit computed
                                    Me%TICOEF3(i, j, k) = Me%TICOEF3(i, j, k)  &               
                                          -  Me%WaterFluxOBoundary(i,j,k)                         &
                                          *  ExteriorProp * DT_V

                                else if (BoundaryCondition == MassConservNullGrad_) then

                                    call NullGradProp(BoundaryProp, i, j, k, Found)

                                    if (Found) then

                                        Me%TICOEF3(i, j, k) = BoundaryProp

                                    else 
                                        
                                        Me%TICOEF3(i, j, k) = Me%ExternalVar%Prop(i, j, k)

                                    endif

                                    Me%COEF3%E(i,j,k) = 1.
                                    Me%COEF3%D(i,j,k) = 0.
                                    Me%COEF3%F(i,j,k) = 0.


                                endif



                            else   cd2
                                !Mass is being taken from the domain by default the 
                                ! flow is positive when is  going out and this 
                                ! process is implicit computed 
                                Me%COEF3%E(i,j,k) =  Me%COEF3%E(i,j,k)   &
                                      + Me%WaterFluxOBoundary(i,j,k) * DT_V 
                            end if cd2                                                      
                            
                         else if (BoundaryCondition == ImposedValue_ .or. BoundaryCondition == SubModel_) then cd1


                            !Property in the exterior constant in time and
                            ! boundary property equal to the exterior

                            Me%TICOEF3(i, j, k) = ExteriorProp

                            Me%COEF3%D(i,j,k) = 0.0
                            Me%COEF3%E(i,j,k) = 1.0
                            Me%COEF3%F(i,j,k) = 0.0

                         else if (BoundaryCondition == NullGradient_ .or. BoundaryCondition == CyclicBoundary_) then

                            Me%TICOEF3(i, j, k) = Me%ExternalVar%PROP(i,j,k)

                            Me%COEF3%D(i,j,k) = 0.0
                            Me%COEF3%E(i,j,k) = 1.0
                            Me%COEF3%F(i,j,k) = 0.0

                         else if (BoundaryCondition /= NullGradient_ .and. BoundaryCondition /= CyclicBoundary_) then
                            
                            write(*,*) ' This is not a valid boundary Condition'
                            stop 'sub. OpenBoundaryCondition - ModuleAdvectionDiffusion - ERR02'
                   
                        end if cd1
                    end if cd3
                end do do2
            end if cd9
        end do do1
        end do do3
        !!$OMP END DO
               
        !!$OMP END PARALLEL   
                
        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "OpenBoundaryCondition")
                                                                   
        !----------------------------------------------------------------------

    end subroutine OpenBoundaryCondition

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine FluxAtOpenBoundary()

        !Arguments-------------------------------------------------------------

        !local-----------------------------------------------------------------

        integer                            :: i, j, k  
        integer                            :: ILB, IUB
        integer                            :: JLB, JUB
        integer                            :: KLB, KUB
        integer, pointer, dimension(:,:,:) :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:) :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:) :: ComputeFacesW3D
        integer                            :: CHUNK

        !----------------------------------------------------------------------

        ILB = Me%WorkSize%ILB
        IUB = Me%WorkSize%IUB

        JLB = Me%WorkSize%JLB
        JUB = Me%WorkSize%JUB

        KUB = Me%WorkSize%KUB

        ComputeFacesU3D => Me%ExternalVar%ComputeFacesU3D
        ComputeFacesV3D => Me%ExternalVar%ComputeFacesV3D
        ComputeFacesW3D => Me%ExternalVar%ComputeFacesW3D

        if (MonitorPerformance) call StartWatch ("ModuleAdvectionDiffusion", "FluxAtOpenBoundary")

        CHUNK = CHUNK_J(JLB,JUB)
        
        !$OMP PARALLEL PRIVATE(i,j,k,KLB)

        !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
do3 :   do j = JLB, JUB
do1 :   do i = ILB, IUB
cd1 :       if (Me%ExternalVar%BoundaryPoints2D(i,j) == 1) then 
                KLB = ABS(Me%ExternalVar%KFloorZ(i,j))

do2 :           do k = KLB, KUB 

                    Me%WaterFluxOBoundary(i, j, k) =                                           &            
                        Me%ExternalVar%Wflux_X(i  , j  , k)   * ComputeFacesU3D(i  , j  , k  ) &
                      - Me%ExternalVar%Wflux_X(i  , j+1, k)   * ComputeFacesU3D(i  , j+1, k  ) &            
                      + Me%ExternalVar%Wflux_Y(i  , j  , k)   * ComputeFacesV3D(i  , j  , k  ) &
                      - Me%ExternalVar%Wflux_Y(i+1, j  , k)   * ComputeFacesV3D(i+1, j  , k  ) &
                      + Me%ExternalVar%Wflux_Z(i  , j  , k)   * ComputeFacesW3D(i  , j  , k  ) &            
                      - Me%ExternalVar%Wflux_Z(i  , j  , k+1) * ComputeFacesW3D(i  , j  , k+1) &            
                      -(Me%ExternalVar%VolumeZ(i  , j  , k)                                    &
                       -Me%ExternalVar%VolumeZOld(i, j, k))                                    &
                      / Me%ExternalVar%DTProp
                end do do2

            end if cd1
        end do do1
        end do do3
        !$OMP END DO

        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleAdvectionDiffusion", "FluxAtOpenBoundary")

        nullify(ComputeFacesU3D)
        nullify(ComputeFacesV3D)
        nullify(ComputeFacesW3D)


        !----------------------------------------------------------------------

    end subroutine FluxAtOpenBoundary

    !--------------------------------------------------------------------------

    subroutine Set_Internal_State(DTProp, CellFluxes, AdvMethodH, TVDLimitationH,        &
                                  AdvMethodV, TVDLimitationV, SchmidtCoef_V,             &
                                  SchmidtBackground_V, Schmidt_H, NullDif)

        !Arguments-------------------------------------------------------------
        real,       intent(IN)                      :: DTProp
        integer,    intent(IN)                      :: AdvMethodH, TVDLimitationH
        integer,    intent(IN)                      :: AdvMethodV, TVDLimitationV
        logical,    intent(IN)                      :: CellFluxes    
        real,       intent(IN)                      :: SchmidtCoef_V, SchmidtBackground_V
        real,       intent(IN)                      :: Schmidt_H
        logical,    intent(IN)                      :: NullDif

        !----------------------------------------------------------------------

        !If DTProp is different from previous DT Prop or 
        !the time the last time advection diffusion was called is different
        !from the current time, State%VertAdv & HorAdv must be set true 

cd2 :   if ((Me%ExternalVar%DTProp    /= DTProp) .OR.                                           &
            (Me%LastCalc              /= Me%Now)) then

            Me%State%VertAdv = ON
            Me%State%HorAdv  = ON

            Me%State%VertDif = ON
            Me%State%HorDif  = ON
        else

            !If parameters for AdvMethodH are diferent 
            if ((Me%ExternalVar%AdvMethodH      /= AdvMethodH    )              .OR.            &
                (Me%ExternalVar%TVDLimitationH  /= TVDLimitationH)              .OR.            &
                (Me%ExternalVar%AdvMethodH      == P2_TVD)) then 
                Me%State%HorAdv  = ON
            else
                Me%State%HorAdv  = OFF
            endif

            if ((Me%ExternalVar%AdvMethodV      /= AdvMethodV    )              .OR.            &
                (Me%ExternalVar%TVDLimitationV  /= TVDLimitationV)              .OR.            &
                (Me%ExternalVar%AdvMethodV      == P2_TVD)) then 
                Me%State%VertAdv  = ON
            else
                Me%State%VertAdv  = OFF
            endif

            if ((Me%ExternalVar%SchmidtCoef_V       /= SchmidtCoef_V)           .OR.            &
                (Me%ExternalVar%SchmidtBackground_V /= SchmidtBackground_V)     .OR.            &
                 NullDif) then
                Me%State%VertDif = ON
            else
                Me%State%VertDif = OFF
            endif

            if ((Me%ExternalVar%Schmidt_H           /= Schmidt_H)               .OR.            &
                NullDif) then
                Me%State%HorDif = ON
            else
                Me%State%HorDif = OFF
            endif

        end if cd2

cd3 :   if (CellFluxes) then
            Me%State%CellFluxes = ON
        else
            Me%State%CellFluxes = OFF
        end if cd3


cd6 :   if (associated(Me%ExternalVar%ReferenceProp)) then
            Me%State%OpenBoundary = ON

            if ((Me%ExternalVar%BoundaryCondition /= MassConservation_          ) .AND. &
                (Me%ExternalVar%BoundaryCondition /= ImposedValue_              ) .AND. &
                (Me%ExternalVar%BoundaryCondition /= SubModel_                  ) .AND. &
                (Me%ExternalVar%BoundaryCondition /= Orlanski_                  ) .AND. &
                (Me%ExternalVar%BoundaryCondition /= NullGradient_              ) .AND. &
                (Me%ExternalVar%BoundaryCondition /= CyclicBoundary_            ) .AND. &
                (Me%ExternalVar%BoundaryCondition /= MassConservNullGrad_       ))      &

                stop 'Set_Internal_State - ModuleAdvectionDiffusion - ERR01'
        else
            Me%State%OpenBoundary = OFF
        end if cd6


        !----------------------------------------------------------------------

    end subroutine Set_Internal_State

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR  

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillAdvectionDiffusion(AdvectionDiffusionID, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: AdvectionDiffusionID
        integer, optional, intent(OUT) :: STAT


        !External--------------------------------------------------------------

        integer :: ready_             
        integer :: STAT_CALL

        !Local-----------------------------------------------------------------

        integer :: STAT_, nUsers
        !griflet
        integer                 :: p
        type(T_VECGW), pointer  :: VECGW

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(AdvectionDiffusionID, ready_)

cd1 :   if (ready_ /= OFF_ERR_) then
            
            nUsers = DeassociateInstance(mADVECTIONDIFFUSION_,  Me%InstanceID)

            if (nUsers == 0) then

                nUsers = DeassociateInstance (mTIME_,           Me%ObjTime)
                if (nUsers == 0) stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR01'

                nUsers = DeassociateInstance (mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR02'
                
                nUsers = DeassociateInstance (mHORIZONTALMAP_,  Me%ObjHorizontalMap)
                if (nUsers == 0) stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR03'
                
                nUsers = DeassociateInstance (mGEOMETRY_,       Me%ObjGeometry)
                if (nUsers == 0) stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR04'



                deallocate(Me%Diffusion_CoeficientX, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR06'
                nullify   (Me%Diffusion_CoeficientX)

                deallocate(Me%Diffusion_CoeficientY, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR06a'
                nullify   (Me%Diffusion_CoeficientY)

                deallocate(Me%Diffusion_CoeficientZ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR06b'
                nullify   (Me%Diffusion_CoeficientZ)

                deallocate(Me%Fluxes%AdvFluxX, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR07'
                nullify   (Me%Fluxes%AdvFluxX   ) 


                deallocate(Me%Fluxes%AdvFluxY, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR08'
                nullify   (Me%Fluxes%AdvFluxY   ) 


                deallocate(Me%Fluxes%AdvFluxZ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR09'
                nullify   (Me%Fluxes%AdvFluxZ) 


                deallocate(Me%Fluxes%DifFluxX, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR10'
                nullify   (Me%Fluxes%DifFluxX) 

 
                deallocate(Me%Fluxes%DifFluxY, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR11'
                nullify   (Me%Fluxes%DifFluxY) 

 
                deallocate(Me%Fluxes%DifFluxZ, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR12'
                nullify   (Me%Fluxes%DifFluxZ) 

 
                deallocate(Me%WaterFluxOBoundary, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR13'
                nullify   (Me%WaterFluxOBoundary)


                deallocate(Me%COEF3%D, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR14a'
                nullify   (Me%COEF3%D) 

                deallocate(Me%COEF3%E, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR14b'
                nullify   (Me%COEF3%E) 

                deallocate(Me%COEF3%F, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR14c'
                nullify   (Me%COEF3%F) 

                deallocate(Me%COEF3_VertAdv%D_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR16f'
                nullify   (Me%COEF3_VertAdv%D_Flux) 

                deallocate(Me%COEF3_VertAdv%E_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR16g'
                nullify   (Me%COEF3_VertAdv%E_Flux)

                deallocate(Me%COEF3_VertAdv%C_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR16h'
                nullify   (Me%COEF3_VertAdv%C_Flux) 

                deallocate(Me%COEF3_VertAdv%F_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR16i'
                nullify   (Me%COEF3_VertAdv%F_Flux)


                deallocate(Me%COEF3_HorAdvXX%C_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR18h'
                nullify   (Me%COEF3_HorAdvXX%C_Flux) 

                deallocate(Me%COEF3_HorAdvXX%D_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR18f'
                nullify   (Me%COEF3_HorAdvXX%D_Flux) 

                deallocate(Me%COEF3_HorAdvXX%E_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR18g'
                nullify   (Me%COEF3_HorAdvXX%E_Flux) 

                deallocate(Me%COEF3_HorAdvXX%F_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR18i'
                nullify   (Me%COEF3_HorAdvXX%F_Flux) 

                deallocate(Me%COEF3_HorAdvYY%D_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR20f'
                nullify   (Me%COEF3_HorAdvYY%D_Flux) 

                deallocate(Me%COEF3_HorAdvYY%E_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR20g'
                nullify   (Me%COEF3_HorAdvYY%E_Flux) 

                deallocate(Me%COEF3_HorAdvYY%C_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR20h'
                nullify   (Me%COEF3_HorAdvYY%C_Flux) 

                deallocate(Me%COEF3_HorAdvYY%F_Flux, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR20i'
                nullify   (Me%COEF3_HorAdvYY%F_Flux) 



                deallocate(Me%TICOEF3, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR21'
                nullify   (Me%TICOEF3) 

                !griflet
                do p = 1, Me%MaxThreads
                    VECGW => Me%THOMAS%VEC(p)
                    deallocate(VECGW%G)
                    deallocate(VECGW%W)
                enddo 
                deallocate(Me%THOMAS%VEC)
                deallocate(Me%THOMAS%COEF3)
                deallocate(Me%THOMAS)

                deallocate(Me%VECW, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR22'
                nullify   (Me%VECW)


                deallocate(Me%VECG, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)                                             &
                    stop 'KillAdvectionDiffusion - ModuleAdvectionDiffusion - ERR23'
                nullify   (Me%VECG)

                call DeallocateInstance()

                AdvectionDiffusionID = 0
                STAT_                = SUCCESS_


            endif

        else cd1

            STAT_ = UNKNOWN_

        end if cd1


        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillAdvectionDiffusion

        !--------------------------------------------------------------------------

    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_AdvectionDiffusion), pointer           :: AuxAdvectionDiffusion
        type (T_AdvectionDiffusion), pointer           :: PreviousAdvectionDiffusion

        !Updates pointers
        if (Me%InstanceID == FirstAdvectionDiffusion%InstanceID) then
            FirstAdvectionDiffusion => FirstAdvectionDiffusion%Next
        else
            PreviousAdvectionDiffusion => FirstAdvectionDiffusion
            AuxAdvectionDiffusion      => FirstAdvectionDiffusion%Next
            do while (AuxAdvectionDiffusion%InstanceID /= Me%InstanceID)
                PreviousAdvectionDiffusion => AuxAdvectionDiffusion
                AuxAdvectionDiffusion      => AuxAdvectionDiffusion%Next
            enddo

            !Now update linked list
            PreviousAdvectionDiffusion%Next => AuxAdvectionDiffusion%Next


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

    subroutine Ready (AdvectionDiffusionID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: AdvectionDiffusionID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (AdvectionDiffusionID > 0) then
            call LocateObjAdvectionDiffusion (AdvectionDiffusionID)
            ready_ = VerifyReadLock (mADVECTIONDIFFUSION_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1


        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjAdvectionDiffusion (AdvectionDiffusionID)

        !Arguments-------------------------------------------------------------
        integer                                     :: AdvectionDiffusionID

        !Local-----------------------------------------------------------------

        Me => FirstAdvectionDiffusion
        do while (associated (Me))
            if (Me%InstanceID == AdvectionDiffusionID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                                        &
            stop 'ModuleAdvectionDiffusion - LocateObjAdvectionDiffusion - ERR01'

    end subroutine LocateObjAdvectionDiffusion



end module ModuleAdvectionDiffusion

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------


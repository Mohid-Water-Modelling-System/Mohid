!-------------------------------------------------------------------------
!        IST/MARETEC, Marine Modelling Group, Mohid2000 modelling system
!-------------------------------------------------------------------------
!BOI
! !TITLE: Mohid2000 hydrodynamic model 
! !AUTHORS: Manuel Ruiz Villarreal
! !AFFILIATION: IST/MARETEC, Marine Modelling Group
! !DATE: 2001
! !INTRODUCTION: An interface to subroutines for computing turbulence coefficients 
!                using a one or two equation turbulence closure. Subroutines are taken 
!                from GOTM (General Ocean Turbulence Model, http://www.gotm.net)
!
!EOI
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

Module ModuleTurbGOTM

    use ModuleGlobalData
    use ModuleTime  
    use ModuleEnterData,        only : ReadFileName
    use ModuleFunctions,        only : TimeToString, ChangeSuffix, CHUNK_J, CHUNK_I, CHUNK_K
    use ModuleHorizontalMap,    only : GetBoundaries, GetWaterPoints2D, UnGetHorizontalMap
    use ModuleGeometry
    use ModuleMap
    use ModuleGOTM
    use ModuleStopWatch,        only : StartWatch, StopWatch         

    Implicit none

    private 


    !Subroutines & Functions---------------------------------------------------

    !Constructor
    public  :: StartTurbGOTM        
    private ::      AllocateInstance        
    private ::      AllocateVariables


    !Selector
    public  :: GetViscosityTurbGOTM
    public  :: GetDiffusivityTurbGOTM
    public  :: GetMixingLenghTurbGOTM
    public  :: GetTurbGOTM_TurbEq

    public  :: UngetTurbGOTM
    public  :: UnGetTurbGOTM_TurbEq

    public  :: SetTurbGOTMBottomRugosity
    public  :: SetTurbGOTMSurfaceRugosity
    public  :: SetTurbGOTMBottomShearVelocity
    public  :: SetTurbGOTMWindShearVelocity


    !Modifier
    public  :: TurbGOTM  
    private :: Read_Final_Turbulence_File      
    public  :: Write_Final_Turbulence_File                                     
    
    !Destructor
    public  ::  KillTurbGOTM   
    private ::      DeAllocateInstance

    
    !Management
    private ::      Ready

    !Interfaces----------------------------------------------------------------

    private :: UngetTurb1
    interface  UngetTurbGOTM
        module procedure UngetTurb1
    end interface UngetTurbGOTM

    
    !Types---------------------------------------------------------------------

    type       T_ID
        integer                       :: IDnumber
        character(LEN = StringLength) :: name
    end type T_ID

    type       T_Files
         character(len=PathLength)    :: ConstructData
         character(len=PathLength)    :: FinalTurbulence
         character(len=PathLength)    :: InitialTurbulence
    end type T_Files

    type       T_External
        real,    pointer, dimension(:,:,:) :: DWZ,dzz
        real,    pointer, dimension(:,:,:) :: NN       !SQUARED FREQ. BRUNT-VAISALLA N=SQRT(G/ROo*dB/dZ)
        real,    pointer, dimension(:,:,:) :: SS       !SQUARED FREQ. PRANDTL        M=dU/dZ
        integer, pointer, dimension(:,:,:) :: WaterPoints3D
        integer, pointer, dimension(:,:  ) :: WaterPoints2D
        integer, pointer, dimension(:,:,:) :: ComputeFacesU3D
        integer, pointer, dimension(:,:,:) :: ComputeFacesV3D
        integer, pointer, dimension(:,:,:) :: LandPoints
        integer, pointer, dimension(:,:,:) :: OpenPoints3D       
        real,    pointer, dimension(:,:  ) :: HT
        real,    pointer, dimension(:,:  ) :: u_taub  !shear bottom velocity
        real,    pointer, dimension(:,:  ) :: u_taus  !shear surface velocity
        integer, pointer, dimension(:,:  ) :: KFloorZ
        real,    pointer, dimension(:,:  ) :: SurfaceRugosity
        real,    pointer, dimension(:,:  ) :: BottomRugosity
        type(T_Time)                       :: Now, BeginTime, EndTime
    end type T_External


    type      T_TurbGOTM
        integer                                     :: InstanceID
        type(T_Size3D  )                            :: Size            !Allocated size -> former IIM,  JJM,  KKM
        type(T_Size3D  )                            :: WorkSize        !Work      size -> former IMAX, JMAX, KMAX
        real(8), dimension(:, :, :),  pointer       :: Matrix
        type(T_TurbGOTM), pointer                   :: Next
        real,    pointer, dimension(:,:,:)          :: TKE
        real,    pointer, dimension(:,:,:)          :: L
        real,    pointer, dimension(:,:,:)          :: EPS
        real,    pointer, dimension(:,:,:)          :: NUM  
        real,    pointer, dimension(:,:,:)          :: NUH
        !real,    pointer, dimension(:,:,:)          :: Lupward
        !real,    pointer, dimension(:,:,:)          :: Ldownward
        real,    pointer, dimension(:,:,:)          :: P
        real,    pointer, dimension(:,:,:)          :: B
        real(8), pointer, dimension(:    )          :: P_1D,B_1D,NN_1D,SS_1D,h
!        double precision, pointer, dimension(:)     :: tke_1D,L_1D,eps_1D,num_1D,nuh_1D
        real                           :: DT            = null_real
        !type(T_ID      ) :: ID
        type(T_External) :: ExternalVar
        type(T_Files   ) :: Files
        type (T_Gotm), pointer :: ObjGotm

        !Information associated with the Gotm model that does not change 
        !in time and in space 
        type(T_Gotmparameters   ), pointer      :: ObjGOTMparameters
        !Instance of ModuleTime
        integer                                 :: ObjTime            = 0

        !Instance of ModuleBathymetry
        integer                                 :: ObjGridData        = 0

     
        !Instance of ModuleHorizontalMap
        integer                                 :: ObjHorizontalMap   = 0
       
        !Instance of ModuleHorizontalGrid
        integer                                 :: ObjHorizontalGrid  = 0
       
        !Instance of ModuleMap
        integer                                 :: ObjMap             = 0
 
        !Instance of ModuleGeometry
        integer                                 :: ObjGeometry        = 0

        
    end type T_TurbGOTM
     
    !Global Module Variables
    type (T_TurbGOTM), pointer                      :: FirstObjTurbGOTM
    type (T_TurbGOTM), pointer                      :: Me

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine StartTurbGOTM(TurbGOTMID,                                        &
                             TimeID,                                            &
                             GridDataID,                                        &
                             MapID,                                             &
                             HorizontalMapID,                                   &
                             HorizontalGridID,                                  &
                             GeometryID,                                        &
                             Continuous_Compute,  STAT)

   
        !Arguments---------------------------------------------------------------
        integer                                         :: TurbGOTMID 
        integer                                         :: TimeID
        integer                                         :: GridDataID
        integer                                         :: MapID
        integer                                         :: HorizontalMapID
        integer                                         :: HorizontalGridID
        integer                                         :: GeometryID
        logical                                         :: Continuous_Compute
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_        
        integer                                         :: STAT_CALL
        type(T_Time)                                    :: BeginTime, EndTime
        
        !Local-------------------------------------------------------------------
        integer                                         :: STAT_        
        integer                                         :: GOTM_UNIT

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mTurbGOTM_)) then
            nullify (FirstObjTurbGOTM)
            call RegisterModule (mTurbGOTM_) 
        endif

        call Ready(TurbGOTMID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_          , TimeID          )
            Me%ObjGridData       = AssociateInstance (mGRIDDATA_      , GridDataID      )
            Me%ObjMap            = AssociateInstance (mMAP_           , MapID           )
            Me%ObjHorizontalMap  = AssociateInstance (mHORIZONTALMAP_ , HorizontalMapID )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            Me%ObjGeometry       = AssociateInstance (mGEOMETRY_      , GeometryID      )
            
            call GetComputeTimeLimits(Me%ObjTime,                               &
                                      EndTime   = EndTime,                      &
                                      BeginTime = BeginTime,                    &
                                      STAT      = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR01'     

            !Stores the begin and the end of the water properties compute
            Me%ExternalVar%BeginTime = BeginTime
            Me%ExternalVar%EndTime   = EndTime

            !Actualizes the time
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR02'     

            call GetComputeTimeStep(Me%ObjTime, Me%DT, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR03'     

                                             
            call GetGeometrySize(Me%ObjGeometry, Size = Me%Size,                &
                                 WorkSize = Me%WorkSize, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR03'     

                 
            !Opens the data file with GOTM turbulence information
            call ReadFileName('TURB_GOTM', Me%Files%ConstructData,              &
                               Message = "GOTM Data File", STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR04'     

            if (Continuous_Compute) then
                
                !Data From previous run
                call ReadFileName('TURB_INI', Me%Files%InitialTurbulence,         &
                                  Message = "GOTM Initial Conditions",            &
                                  STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR05'     

            end if


            call ReadFileName('TURB_FIN', Me%Files%FinalTurbulence,             &
                              Message = "GOTM Final Conditions",                &
                              STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR06'     
        
            call UnitsManager (GOTM_UNIT, OPEN_FILE) 

            !Initializes the value of the parameters in the turbulence module GOTM
            call init_turbulence_parameters(Me%ObjGotmParameters,                        &
                                            GOTM_UNIT, Me%Files%ConstructData,           &
                                            STAT=STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'StartTurbGOTM - ModuleTurbGOTM - ERR05'     
           
            call AllocateVariables 
            
            if (Continuous_Compute) then
                call Read_Final_Turbulence_File
            end if

            !Returns ID
            TurbGOTMID          = Me%InstanceID
            STAT_               = SUCCESS_

        else cd0
            
            stop 'ModuleTurbGOTM - StartTurbGOTM - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartTurbGOTM

    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_TurbGOTM), pointer                      :: NewObjTurbGOTM
        type (T_TurbGOTM), pointer                      :: PreviousObjTurbGOTM


        !Allocates new instance
        allocate (NewObjTurbGOTM)
        nullify  (NewObjTurbGOTM%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjTurbGOTM)) then
            FirstObjTurbGOTM         => NewObjTurbGOTM
            Me                       => NewObjTurbGOTM
        else
            PreviousObjTurbGOTM      => FirstObjTurbGOTM
            Me                       => FirstObjTurbGOTM%Next
            do while (associated(Me))
                PreviousObjTurbGOTM  => Me
                Me                   => Me%Next
            enddo
            Me                       => NewObjTurbGOTM
            PreviousObjTurbGOTM%Next => NewObjTurbGOTM
        endif

        Me%InstanceID = RegisterNewInstance (mTurbGOTM_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    subroutine AllocateVariables

        !Arguments-------------------------------------------------------------


        !External--------------------------------------------------------------        
           
        integer :: STAT_CALL

        !Local-----------------------------------------------------------------

        integer :: ILB, IUB 
        integer :: JLB, JUB
        integer :: KLB, KUB
        !----------------------------------------------------------------------

        ILB = Me%Size%ILB
        IUB = Me%Size%IUB

        JLB = Me%Size%JLB
        JUB = Me%Size%JUB

        KLB = Me%Size%KLB
        KUB = Me%Size%KUB


        allocate(Me%TKE         (ILB:IUB, JLB:JUB, KLB:KUB))         
        allocate(Me%L           (ILB:IUB, JLB:JUB, KLB:KUB))         
        allocate(Me%EPS         (ILB:IUB, JLB:JUB, KLB:KUB))         
        allocate(Me%NUM         (ILB:IUB, JLB:JUB, KLB:KUB))         
        allocate(Me%NUH         (ILB:IUB, JLB:JUB, KLB:KUB))         
        allocate(Me%P           (ILB:IUB, JLB:JUB, KLB:KUB))         
        allocate(Me%B           (ILB:IUB, JLB:JUB, KLB:KUB))         
        !allocate(Me%Lupward     (ILB:IUB, JLB:JUB, KLB:KUB))         
        !allocate(Me%Ldownward   (ILB:IUB, JLB:JUB, KLB:KUB))         

        Me%TKE                  = Me%ObjGOTMParameters%K_MIN 
        Me%L                    = Me%ObjGOTMParameters%L_MIN 
        Me%EPS                  = Me%ObjGOTMParameters%EPS_MIN     
        Me%NUM                  = 0.
        Me%NUH                  = 0.
        Me%P                    = null_real
        Me%B                    = null_real
        !Me%Lupward              = null_real
        !Me%Ldownward            = null_real
                         
        ! Allocates matrices for GOTM 1D column. 
        ! The local matrices in GOTM are allocated with the maximum numbers of layers, 
        ! although the computation limits are 0:nlev and nlev is different for every point.

        !Variables transferred to ModuleGOTM as arguments
        allocate(Me%NN_1D   (KLB:KUB))         
        allocate(Me%SS_1D   (KLB:KUB))         
        allocate(Me%P_1D    (KLB:KUB))         
        allocate(Me%B_1D    (KLB:KUB))         
        allocate(Me%h       (KLB:KUB))         
        allocate(Me%ObjGotm         )

        Me%NN_1D                    = null_real
        Me%SS_1D                    = null_real
        Me%P_1D                    = null_real
        Me%B_1D                    = null_real
        Me%h                       = null_real
                
        !Here 1D arrays of turbulent magnitudes(tke,eps,num,nuh,L) are allocated 
        call init_turbulence(Me%ObjGotm, Me%ObjGOTMParameters, KLB, KUB, STAT=STAT_CALL) 
        if (STAT_CALL .NE. SUCCESS_)  stop 'AllocateVariables - ModuleTurbGOTM - ERR01'    

        Me%ObjGOTM%Parameters => Me%ObjGOTMParameters

        nullify(Me%ExternalVar%u_taus)
        allocate(Me%ExternalVar%u_taus (ILB:IUB, JLB:JUB), STAT = STAT_CALL)         
        if (STAT_CALL .NE. SUCCESS_) stop 'AllocateVariables - ModuleTurbGOTM - ERR02'
        Me%ExternalVar%u_taus = 0.

    end subroutine AllocateVariables 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MO 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine TurbGOTM(TurbGOTMID, NN, SS, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: TurbGOTMID
        real,    pointer, dimension(:,:,:)  :: NN       !SQUARE FREQ. BRUNT-VAISALLA N=SQRT(G/ROo*dB/dZ)
        real,    pointer, dimension(:,:,:)  :: SS       !SQUARE FREQ. PRANDTL        M=dU/dZ
        integer, optional,     intent(OUT)  :: STAT

        !Local-----------------------------------------------------------------
        integer,    pointer, dimension(:,:) :: BoundaryPoints2D
        real(8),    pointer, dimension(:)   :: h,NN_1D,SS_1D,P_1D,B_1D
        real(8),    pointer, dimension(:)   :: tke_1D,eps_1D,L_1D,num_1D,nuh_1D   
        real(8)                             :: depth,u_taub,u_taus,dt,z0s,z0b
        integer                             :: KBottom,nlev
        integer                             :: i, j, k, ILB, IUB, JLB, JUB, KLB, KUB                            
        integer                             :: STAT_                
        integer                             :: STAT_CALL
        integer                             :: ready_            
        integer                             :: CHUNK            

        !----------------------------------------------------------------------                         

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)    

        if (MonitorPerformance) &
            call StartWatch ("ModuleTurbGOTM", "TurbGOTM")

cd1 :   if (ready_ .EQ. IDLE_ERR_) then

            !Time Properties - Actualises CurrentTime
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR01'
                
            Me%ExternalVar%NN   => NN
            Me%ExternalVar%SS   => SS

            NN_1D               => Me%NN_1D
            SS_1D               => Me%SS_1D
            P_1D                => Me%P_1D
            B_1D                => Me%B_1D
            h                   => Me%h            

            ILB = Me%WorkSize%ILB   
            IUB = Me%WorkSize%IUB
            JLB = Me%WorkSize%JLB 
            JUB = Me%WorkSize%JUB
            KLB = Me%WorkSize%KLB
            KUB = Me%WorkSize%KUB

            !WaterPoints3D
            call GetWaterPoints3D(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR02'

            !ComputeFacesU3D
            call GetComputeFaces3D(Me%ObjMap,                                           &
                                   ComputeFacesU3D = Me%ExternalVar%ComputeFacesU3D,    &  
                                   STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR03'
 
            !ComputeFacesV3D
            call GetComputeFaces3D(Me%ObjMap,                                           &
                                   ComputeFacesV3D = Me%ExternalVar%ComputeFacesV3D,    &
                                   STAT = STAT_CALL)   
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR04'

            !LandPoints
            call GetLandPoints3D(Me%ObjMap, Me%ExternalVar%LandPoints, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR05'
            
            !OpenPoints3D
            call GetOpenPoints3D(Me%ObjMap, Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR06'

            call GetWaterPoints2D(Me%ObjHorizontalMap,                                  &
                                  Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR07'

            call GetBoundaries(Me%ObjHorizontalMap, BoundaryPoints2D, STAT= STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR08'
                       
            call GetGeometryKFloor(Me%ObjGeometry, Z = Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR09'

            call GetGeometryDistances(Me%ObjGeometry, DWZ = Me%ExternalVar%DWZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR10'


            call GetGeometryWaterColumn(Me%ObjGeometry,                                 &
                                        WaterColumn = Me%ExternalVar%HT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR11'
                
  
            !In principle time step for turbulence is the same as in the model
            call GetComputeTimeStep(Me%ObjTime, Me%DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR12'

            dt = Me%DT

            ! Production terms for turbulence equations are computed here.
            ! We consider them as a hydrodynamic determined property.
            call Production
            
            !Call to module GOTM, 1D column model for turbulence   
            

            ! We only allocate ONCE space for these 1D variables. This is done in init_turbulence in ModuleGOTM
            ! Here we associate these variables to those allocated there.
            call Associate_to_ExportGOTM (Me%ObjGOTM, Tke_1D, L_1D, eps_1D, num_1D, nuh_1D)
            
            !Doesn't work CHUNK = CHUNK_J(JLB,JUB)
            !Parallelize with private variables that are arrays doesn't work.
            !Doesn't work !$OMP PARALLEL PRIVATE(I,J,k,Kbottom,Depth,u_taus,u_taub,z0b,z0s,nlev,NN_1D,SS_1D,P_1D,B_1D,Tke_1D,L_1D,eps_1D,num_1D,nuh_1D,h) 
            ! We don't compute turbulence coefficients at the limits of the domain. 
            !Doesn't work !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
doj :       do J = JLB+1, JUB-1
doi :       do I = ILB+1, IUB-1 
      
ifwp :          if (Me%ExternalVar%OpenPoints3D(i,j,KUB) .EQ. Openpoint) then       
                        
                    Kbottom = Me%ExternalVar%KFloorZ        (i,j) 
                    Depth   = Me%ExternalVar%HT             (i,j)
                    u_taus  = Me%ExternalVar%u_taus         (i,j)
                    u_taub  = Me%ExternalVar%u_taub         (i,j)
                    z0b     = Me%ExternalVar%BottomRugosity (i,j)
                    z0s     = Me%ExternalVar%SurfaceRugosity(i,j) 

                    nlev = KUB - KBottom + 1


                    ! We must be careful with the limits of the matrix, as GOTM computes from 0 to nlev
                    ! Kbottom face correspond to 0 in GOTM notation
                    ! Level 0 is needed for boundary confitions of turbulent magnitudes...
                    
dok :               do k=KBottom,KUB+1
                        NN_1D (k-KBottom) = dble(Me%ExternalVar%NN(i,j,k))
                        SS_1D (k-KBottom) = dble(Me%ExternalVar%SS(i,j,k))
                        P_1D  (k-KBottom) = dble(Me%P     (i,j,k))
                        B_1D  (k-KBottom) = dble(Me%B     (i,j,k))
                        Tke_1D(k-KBottom) = dble(Me%Tke   (i,j,k))
                        L_1D  (k-KBottom) = dble(Me%L     (i,j,k))
                        eps_1D(k-KBottom) = dble(Me%eps   (i,j,k))
                        num_1D(k-KBottom) = dble(Me%num   (i,j,k))
                        nuh_1D(k-KBottom) = dble(Me%nuh   (i,j,k))
                    end do dok

                    ! h is the distance between vertical faces (i.e. where turbulent quantities are defined)
                    ! See GOTM report page 34. In MOHID notation h = DWZ
dok2:               do k=KBottom,KUB
                     h(k+1-KBottom) = Me%ExternalVar%DWZ(i,j,k)
                    end do dok2
            
                    !The information read in from 3d variables (like Me%tke) is transferred to ModuleGOTM
                    call AssociateExportGOTM (Me%ObjGOTM, Tke_1D, L_1D, eps_1D, num_1D, nuh_1D)

                    !The resolution of turbulence equations is done here (ModuleGOTM.f90)
                    call do_turbulence(Me%ObjGOTM, nlev,dt,depth,u_taus,u_taub,z0s,z0b, &
                                                   h,NN_1D,SS_1D,P_1D,B_1D)
                                    
            
dok4 :              do k=KBottom,KUB+1
                        Me%Tke   (i,j,k) = real(Tke_1D(k-KBottom))  
                        Me%L     (i,j,k) = real(L_1D(k-KBottom  ))  
                        Me%eps   (i,j,k) = real(eps_1D(k-KBottom))  
                        Me%num   (i,j,k) = real(num_1D(k-KBottom))  
                        Me%nuh   (i,j,k) = real(nuh_1D(k-KBottom)) 
                    end do dok4

                end if ifwp

            end do doi
            end do doj
            !Doesn't work !$OMP END DO

            !Doesn't work !$OMP END PARALLEL          
            
            !Values at open boundary points at the limits of the domain 
            ! are set to the values of the nearest interior point. Null_gradient
            i=IUB
                       
            CHUNK = CHUNK_J(JLB,JUB)
            !$OMP PARALLEL PRIVATE(J,k)
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j=JLB,JUB
                if(BoundaryPoints2D(i,j) == 1) then
                    Kbottom = Me%ExternalVar%KFloorZ(I,J)
                    do k=Kbottom,KUB+1
                        Me%num   (i,j,k) = Me%num(i-1,j,k)  
                        Me%nuh   (i,j,k) = Me%nuh(i-1,j,k)
                    end do 
                end if 
            end do
            !$OMP END DO

            !$OMP END PARALLEL          

            i=ILB
            CHUNK = CHUNK_J(JLB,JUB)
            !$OMP PARALLEL PRIVATE(j,k)
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j=JLB,JUB
                if(BoundaryPoints2D(i,j) == 1) then
                    Kbottom = Me%ExternalVar%KFloorZ(I,J)
                    do k=Kbottom,KUB+1
                        Me%num   (i,j,k) = Me%num(i+1,j,k)  
                        Me%nuh   (i,j,k) = Me%nuh(i+1,j,k)
                    end do 
                end if 
            end do
            !$OMP END DO

            !$OMP END PARALLEL          

            j=JUB
            CHUNK = CHUNK_I(ILB,IUB)
            !$OMP PARALLEL PRIVATE(i,k)
            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i=ILB,IUB
                if(BoundaryPoints2D(i,j) == 1) then
                    Kbottom = Me%ExternalVar%KFloorZ(I,J)
                    do k=Kbottom,KUB+1
                        Me%num   (i,j,k) = Me%num(i,j-1,k)  
                        Me%nuh   (i,j,k) = Me%nuh(i,j-1,k)
                    end do 
                end if 
            end do
            !$OMP END DO

            !$OMP END PARALLEL          

            j=JLB
            CHUNK = CHUNK_I(ILB,IUB)
            !$OMP PARALLEL PRIVATE(i,k)            
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i=ILB,IUB
                if(BoundaryPoints2D(i,j) == 1) then
                    Kbottom = Me%ExternalVar%KFloorZ(I,J)
                    do k=Kbottom,KUB+1
                        Me%num   (i,j,k) = Me%num(i,j+1,k)  
                        Me%nuh   (i,j,k) = Me%nuh(i,j+1,k)
                    end do 
                end if 
            end do
            !$OMP END DO
            !$OMP END PARALLEL          

            !WaterPoints3D
            call UngetMap(Me%ObjMap, Me%ExternalVar%WaterPoints3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR13'

            !ComputeFacesU3D
            call UngetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesU3D, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR14'
            
            !ComputeFacesV3D
            call UngetMap(Me%ObjMap, Me%ExternalVar%ComputeFacesV3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR15'

            !LandPoints
            call UngetMap(Me%ObjMap,Me%ExternalVar%LandPoints, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR16'
            
            !OpenPoints3D
            call UngetMap(Me%ObjMap,Me%ExternalVar%OpenPoints3D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR17'

            !WaterPoints2D
            call UngetHorizontalMap(Me%ObjHorizontalMap,                                &
                                    Me%ExternalVar%WaterPoints2D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR18'

            !BoundaryPoints2D
            call UngetHorizontalMap(Me%ObjHorizontalMap, BoundaryPoints2D, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR19'

            !KFloorZ
            call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%KFloorZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR20'

            !DWZ
            call UnGetGeometry(Me%ObjGeometry, Me%ExternalVar%DWZ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR21'

            !HT
            call UnGetGeometry(Me%ObjGeometry,Me%ExternalVar%HT, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'TurbGOTM - ModuleTurbGOTM - ERR22'
           
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
            
        if (MonitorPerformance) &
            call StopWatch ("ModuleTurbGOTM", "TurbGOTM")

    end subroutine TurbGOTM

    !--------------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Various useful variables. 
!
! !INTERFACE:
!   subroutine production(nlev,alpha,num,nuh,P,B)
    subroutine Production
!
! !DESCRIPTION:
!  This subroutine calculates different parameters needed for the
!  turbulence equations such as:
!
!  \begin{itemize}
!  \item P : shear production of turbulent kinetic energy
!  \item B : buoyancy production of turbulent kinetic energy
!  \end{itemize}
!
!  xP is an extra production term which might come from e.g. seagrass friction.
!
! !USES:
!   use meanflow, only: NN,SS,xP,P,B
!   IMPLICIT NONE
!
! !INPUT PARAMETERS:
!   integer, intent(in)  :: nlev
!   double precision, intent(in) :: alpha
!   double precision, intent(in) :: num(0:nlev),nuh(0:nlev)
!
! !INPUT/OUTPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): GOTM code
!
!  $Log$
!
!EOP
!-----------------------------------------------------------------------
!BOC
! In GOTM code
!   P=num*(SS+alpha*NN) + xP
!   B=-nuh*NN 

       !Arguments-------------------------------------------------------------

       !Local-----------------------------------------------------------------
          
 
!        integer :: JLB, JUB 
!        integer :: KLB, KUB
!        integer :: I, J, K 
!        integer :: Kbottom
        real :: alpha

        real,    pointer, dimension(:,:,:) :: NUM, NUH,NN,SS

!        integer, pointer, dimension(:,:  ) :: WaterPoints2D
!        integer, pointer, dimension(:,:  ) :: KFloorZ

        !----------------------------------------------------------------------


        alpha  = Me%ObjGOTMParameters%alpha
           
        NN              => Me%ExternalVar%NN
        SS              => Me%ExternalVar%SS
        NUM             => Me%NUM
        NUH             => Me%NUH

        Me%P = (SS +alpha*NN)*NUM
        Me%B = -NN*NUH

        nullify(NUM            )
        nullify(NUH            )
        nullify(SS             )
        nullify(NN             )

    end subroutine Production
    !----------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !If the user wants to use the values of a previous   
    ! run the read the property values form the final      
    ! results file of a previous run. By default this      
    ! file is in HDF format                                

    subroutine Read_Final_Turbulence_File


        !Local-----------------------------------------------------------------
        real                    :: Year_File, Month_File, Day_File
        real                    :: Hour_File, Minute_File, Second_File
        real                    :: DT_error
        type (T_Time)           :: BeginTime, EndTimeFile, EndTime
        integer                 :: IUB, JUB, KUB, ILB, JLB, KLB
        integer                 :: InitialFile, i, j, k
        logical                 :: exists
        integer                 :: STAT_CALL

        !----------------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        KLB = Me%Size%KLB 
        KUB = Me%Size%KUB 

        inquire (FILE=Me%Files%InitialTurbulence, EXIST = exists) 
        if (.not. exists) &
            stop 'Read_Final_Turbulence_File - ModuleTurbGOTM - ERR01'
               
        call UnitsManager(InitialFile, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)  &
            stop 'Read_Final_Turbulence_File - ModuleTurbGOTM - ERR02' 

        open(Unit = InitialFile, File = Me%Files%InitialTurbulence, &
             Form = 'UNFORMATTED', status = 'OLD', IOSTAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'Read_Final_Turbulence_File - ModuleTurbGOTM - ERR03' 

        !Start reading the end file of the previous run 
        read(InitialFile) Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File

        call SetDate(EndTimeFile, Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File)

        call GetComputeTimeLimits(Me%ObjTime, BeginTime = BeginTime, &
                                  EndTime = EndTime, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  &
            stop 'Read_Final_Turbulence_File - ModuleTurbGOTM - ERR04' 
        
        DT_error = EndTimeFile - BeginTime

        if (abs(DT_error) > 1.e-5) then

            write(*,*) 'The end time of the previous run is different from the start time of this run'
            write(*,*) 'File Time ', Year_File, Month_File, Day_File, Hour_File, Minute_File, Second_File
            write(*,*) 'DT_error ', DT_error
            stop 'Read_Final_Turbulence_File - ModuleTurbGOTM - ERR05'   

        endif

        !Eddy viscosities and diffusivities
        read(InitialFile) (((Me%NUM(I,J,K),     i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
        read(InitialFile) (((Me%NUH(I,J,K),     i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
                   
        !TKE
        read(InitialFile) (((Me%TKE(i, j, k),   i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)

        !L  
        read(InitialFile) (((Me%L(i, j, k),     i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)

        !eps
        read(InitialFile) (((Me%EPS(i, j, k),   i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
                               

        call UnitsManager(InitialFile, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)  &
            stop 'Read_Final_Turbulence_File - ModuleTurbGOTM - ERR06' 

    end subroutine Read_Final_Turbulence_File

    !--------------------------------------------------------------------------

    ! The model writes the final turbulence properties in a binary file
    subroutine Write_Final_Turbulence_File(TurbGOTMID, Overwrite, STAT)

        !Arguments-------------------------------------------------------------
        integer, optional, intent(IN)               :: TurbGOTMID
        logical, optional, intent(IN)               :: Overwrite
        integer, optional, intent(OUT )             :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: Year, Month, Day, Hour, Minute, Second
        integer                                     :: IUB, JUB, KUB, ILB, JLB, KLB
        integer                                     :: FinalFile, i, j, k, STAT_, ready_
        logical                                     :: WriteOK, Overwrite_
        character (Len = Pathlength)                :: filename

        !Begin-----------------------------------------------------------------

        ILB = Me%Size%ILB 
        IUB = Me%Size%IUB
        JLB = Me%Size%JLB 
        JUB = Me%Size%JUB 
        KLB = Me%Size%KLB 
        KUB = Me%Size%KUB 

        WriteOK = .false. 

        if(present(TurbGOTMID))then

            STAT_ = UNKNOWN_

            call Ready(TurbGOTMID, ready_)

            if (ready_ .EQ. IDLE_ERR_) then

                WriteOK = .true.
                STAT_ = SUCCESS_

            else
                
                WriteOK = .false.
                STAT_ = ready_

            end if 

            if (present(STAT)) STAT = STAT_

        else

            WriteOK = .true.
            
        endif

        if(WriteOK)then

            if(present(Overwrite))then
                Overwrite_  = Overwrite
            else    
                Overwrite_  = .true.
            endif

            !Checks if it's at the end of the run 
            !or !if it's supposed to overwrite the final HDF file
            if (Overwrite_) then

                filename = trim(Me%Files%FinalTurbulence)

            else

                filename =  ChangeSuffix(Me%Files%FinalTurbulence,                 &
                                "_"//trim(TimeToString(Me%ExternalVar%Now))//".fin")

            endif

            call UnitsManager(FinalFile, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_Final_Turbulence_File - ModuleTurbGOTM - ERR00'

            open(Unit = FinalFile, File = trim(filename), &
                 Form = 'UNFORMATTED', status = 'UNKNOWN', IOSTAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_Final_Turbulence_File - ModuleTurbGOTM - ERR10'

            !Time Properties - Actualizes CurrentTime
            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Write_Final_Turbulence_File - ModuleTurbGOTM - ERR20'

            !Start writting the final turbulence conditions file
            call ExtractDate(Me%ExternalVar%Now, Year, Month, Day, Hour, Minute, Second)

            write(FinalFile) Year, Month, Day, Hour, Minute, Second

            !Eddy viscosities and diffusivities
            write(FinalFile) (((Me%NUM(I,J,K),      i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
            write(FinalFile) (((Me%NUH(I,J,K),      i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
                   
            !TKE
            write(FinalFile) (((Me%TKE(i, j, k),    i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
                               
            !L  
            write(FinalFile) (((Me%L(i, j, k),      i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)

            !eps
            write(FinalFile) (((Me%EPS(i, j, k),    i = ILB, IUB), j = JLB, JUB), k = KLB, KUB)
                               
            call UnitsManager(FinalFile, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Write_Final_Turbulence_File - ModuleTurbGOTM - ERR30'

        end if

    end subroutine Write_Final_Turbulence_File


    !--------------------------------------------------------------------------
!EOP

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetMixingLenghTurbGOTM(TurbGOTMID, Lupward, Ldownward, STAT)    

        !Arguments-------------------------------------------------------------

        integer                                      :: TurbGOTMID
        real, optional, pointer, dimension(: , :, :) :: Lupward   
        real, optional, pointer, dimension(: , :, :) :: Ldownward   
        integer, optional, intent(OUT) :: STAT

   
        !External--------------------------------------------------------------

        integer :: ready_  

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


cd3 :       if (present(Lupward)) then
                call Read_Lock(mTurbGOTM_, TurbGOTMID)

                Lupward => Me%L
            end if cd3


cd2 :       if (present(Ldownward)) then
                call Read_Lock(mTurbGOTM_, TurbGOTMID)

                Ldownward => Me%L
            end if cd2


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetMixingLenghTurbGOTM

    !--------------------------------------------------------------------------



    !--------------------------------------------------------------------------

    subroutine GetViscosityTurbGOTM(TurbGOTMID, Viscosity, STAT)    

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: STAT

        real, pointer, dimension(:,:,:) :: Viscosity   

        integer                         :: TurbGOTMID

        !External--------------------------------------------------------------

        integer :: ready_  

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
                call Read_Lock(mTurbGOTM_, TurbGOTMID)
                Viscosity => Me%NUM


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetViscosityTurbGOTM

  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------

    subroutine GetDiffusivityTurbGOTM(TurbGOTMID, Diffusivity, STAT)    

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: STAT

        real, pointer, dimension(:,:,:) :: Diffusivity   

        integer                         :: TurbGOTMID

        !External--------------------------------------------------------------
 
        integer :: ready_  

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
                call Read_Lock(mTurbGOTM_, TurbGOTMID)
                Diffusivity => Me%NUH


            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

    !----------------------------------------------------------------------

    end subroutine GetDiffusivityTurbGOTM

    !--------------------------------------------------------------------------

    subroutine GetTurbGOTM_TurbEq(TurbGOTMID, TKE, eps, L, P, B, STAT)    

        !Arguments-------------------------------------------------------------

        integer, optional, intent(OUT) :: STAT

        real, pointer, dimension(:,:,:) :: TKE, eps, L, P, B  

        integer                         :: TurbGOTMID

        !External--------------------------------------------------------------
 
        integer :: ready_  

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
                call Read_Lock(mTurbGOTM_, TurbGOTMID)
                TKE => Me%TKE
                eps => Me%eps
                L   => Me%L
                P   => Me%P
                B   => Me%B

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

    !----------------------------------------------------------------------

    end subroutine GetTurbGOTM_TurbEq

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine UngetTurbGOTM_TurbEq(TurbGOTMID, tke, eps, L, P, B, STAT)

        !Arguments-------------------------------------------------------------

        integer, optional, intent (OUT) :: STAT

        integer                         :: TurbGOTMID

        real, pointer, dimension(:,:,:) :: TKE, eps, L, P, B

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(tke)
            nullify(eps)
            nullify(L)
            nullify(P)
            nullify(B)

            call Read_Unlock(mTurbGOTM_,TurbGOTMID,'UngetTurbGOTM_TurbEq')

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetTurbGOTM_TurbEq

    !--------------------------------------------------------------------------

    subroutine UngetTurb1(TurbGOTMID, Array, STAT)

        !Arguments-------------------------------------------------------------
   
        real, pointer, dimension(:,:,:) :: Array

        integer, optional, intent (OUT) :: STAT

        integer                         :: TurbGOTMID  

        !External--------------------------------------------------------------

        integer :: ready_   

        !Local-----------------------------------------------------------------

        integer :: STAT_              !Auxiliar local variable

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array)

            call Read_Unlock(mTurbGOTM_, TurbGOTMID, 'UngetTurb1')

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UngetTurb1


    !----------------------------------------------------------------------

    subroutine SetTurbGOTMBottomRugosity(TurbGOTMID, BottomRugosity, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: TurbGOTMID
        real, dimension(:,:), pointer   :: BottomRugosity
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then


            Me%ExternalVar%BottomRugosity => BottomRugosity            
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetTurbGOTMBottomRugosity
    
    !----------------------------------------------------------------------

    subroutine SetTurbGOTMSurfaceRugosity(TurbGOTMID, SurfaceRugosity, STAT)

        !Arguments-------------------------------------------------------------
        integer                         :: TurbGOTMID
        real, dimension(:,:), pointer   :: SurfaceRugosity
        integer, optional, intent(OUT)  :: STAT

        !External--------------------------------------------------------------
        integer                         :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                         :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then


            Me%ExternalVar%SurfaceRugosity => SurfaceRugosity
            
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetTurbGOTMSurfaceRugosity


    !----------------------------------------------------------------------

    subroutine SetTurbGOTMWindShearVelocity(TurbGOTMID, WindShearVelocity, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: TurbGOTMID
        real,    pointer, dimension(:,:)    :: WindShearVelocity
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                             :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then


            Me%ExternalVar%u_taus => WindShearVelocity
            
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetTurbGOTMWindShearVelocity

    !----------------------------------------------------------------------
    
    subroutine SetTurbGOTMBottomShearVelocity(TurbGOTMID, BottomShearVelocity, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: TurbGOTMID
        real,    pointer, dimension(:,:)    :: BottomShearVelocity
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              
        
        !Local-----------------------------------------------------------------
        integer                             :: STAT_            

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)  
        
cd1 :   if (ready_ == IDLE_ERR_)then


            Me%ExternalVar%u_taub => BottomShearVelocity
            
   
            STAT_ = SUCCESS_

        else cd1

            STAT_ = ready_

        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine SetTurbGOTMBottomShearVelocity

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillTurbGOTM(TurbGOTMID, STAT)

        !Arguments-------------------------------------------------------------
        integer, optional, intent(OUT)  :: STAT
        integer                         :: TurbGOTMID  

        !External--------------------------------------------------------------
        integer                         :: ready_             

        !Local-----------------------------------------------------------------
        integer                         :: STAT_, nUsers

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(TurbGOTMID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mTurbGOTM_,  Me%InstanceID)
           
            if (nUsers == 0) then

                !Deassociates External Instances
                nUsers = DeassociateInstance (mTIME_          , Me%ObjTime          )
                nUsers = DeassociateInstance (mGRIDDATA_      , Me%ObjGridData      )
                nUsers = DeassociateInstance (mMAP_           , Me%ObjMap           )
                nUsers = DeassociateInstance (mHORIZONTALMAP_ , Me%ObjHorizontalMap )
                nUsers = DeassociateInstance (mHORIZONTALGRID_, Me%ObjHorizontalGrid)
                nUsers = DeassociateInstance (mGEOMETRY_      , Me%ObjGeometry      )

                call Write_Final_Turbulence_File
                call KillGotmMemory

                deallocate(Me%ExternalVar%u_taus)

                !Deallocates Instance
                call DeallocateInstance

                TurbGOTMID = 0

                STAT_      = SUCCESS_

            else 
                STAT_ = NBUSERS_ERR_
            end if 
       
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine KillTurbGOTM

    !------------------------------------------------------------------------
    
    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_TurbGOTM), pointer          :: AuxObjTurbGOTM
        type (T_TurbGOTM), pointer          :: PreviousObjTurbGOTM

        !Updates pointers
        if (Me%InstanceID == FirstObjTurbGOTM%InstanceID) then
            FirstObjTurbGOTM => FirstObjTurbGOTM%Next
        else
            PreviousObjTurbGOTM => FirstObjTurbGOTM
            AuxObjTurbGOTM      => FirstObjTurbGOTM%Next
            do while (AuxObjTurbGOTM%InstanceID /= Me%InstanceID)
                PreviousObjTurbGOTM => AuxObjTurbGOTM
                AuxObjTurbGOTM      => AuxObjTurbGOTM%Next
            enddo

            !Now update linked list
            PreviousObjTurbGOTM%Next => AuxObjTurbGOTM%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine KillGotmMemory

        !Arguments-------------------------------------------------------------

        !Begin-----------------------------------------------------------------
     
        deallocate(Me%TKE       )
        deallocate(Me%L         )    
        deallocate(Me%EPS       )
        deallocate(Me%NUM       )
        deallocate(Me%NUH       )
        deallocate(Me%P         )
        deallocate(Me%B         )
        !deallocate(Me%Lupward   )         
        !deallocate(Me%Ldownward )         
        deallocate(Me%P_1D      )
        deallocate(Me%B_1D      )
        deallocate(Me%NN_1D     )
        deallocate(Me%SS_1D     )
        deallocate(Me%h         )
        deallocate(Me%ObjGotm   )


    end subroutine KillGotmMemory


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready(ObjTurbGOTM_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTurbGOTM_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjTurbGOTM_ID > 0) then
            call LocateObjTurbGOTM (ObjTurbGOTM_ID)
            ready_ = VerifyReadLock (mTurbGOTM_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

        !--------------------------------------------------------------------------

    subroutine LocateObjTurbGOTM (ObjTurbGOTMID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjTurbGOTMID

        !Local-----------------------------------------------------------------

        Me => FirstObjTurbGOTM
        do while (associated (Me))
            if (Me%InstanceID == ObjTurbGOTMID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleTurbGOTM - LocateObjTurbGOTM - ERR01'

    end subroutine LocateObjTurbGOTM

    !--------------------------------------------------------------------------


end module ModuleTurbGOTM

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------









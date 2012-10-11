!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 2
! MODULE        : BoxDif
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to integrate values in 2D and 3D monitoring boxes
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

Module ModuleBoxDif

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions,        only: RodaXY, IsOdd
    use ModuleGridData,         only: WriteGridData
    use ModuleHorizontalGrid,   only: GetZCoordinates,  GetGridBorderPolygon,           &
                                      GetHorizontalGridSize, UnGetHorizontalGrid,       &
                                      GetCornersCoordinates  
    use ModuleTimeSerie,        only: StartTimeSerie, WriteTimeSerieLine, KillTimeSerie
    use ModuleDrawing

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartBoxDif
    private ::      AllocateInstance
    private ::      ConstructBoxes2D
    private ::      ConstructBoxes3D
    private ::          ReadGlobalOptions
    private ::          AllocateVariables2D
    private ::          AllocateVariables3D
    private ::          InitializeEnvironmentBox2D
    private ::          InitializeEnvironmentBox3D
    private ::          ReadBoxes2D
    private ::          ReadBoxes3D
    private ::              AddBox
    private ::              CountNumberOfBoxes2D
    private ::              CountNumberOfLayers
    private ::              ConstructBox2D
    private ::                  SetLimitsPolygon
    private ::              ConstructBox3D
    private ::          ConvertBoxesToMap2D
    private ::          ConvertBoxesToMap3D
    private ::          FindAdjacentBoxesBoundaries2D
    private ::          FindAdjacentBoxesBoundaries3D
    private ::          WriteBoxes
    private ::      ConstructOutputFluxes2D
    private ::      ConstructOutputFluxes3D
    private ::          FluxesTimeSerieHeader2D
    private ::          FluxesTimeSerieHeader3D
    private ::          AddFluxesTimeSerie
    private ::      ConstructOutputScalar2D
    private ::      ConstructOutputScalar3D
    private ::          ScalarTimeSerieHeader2D
    private ::          ScalarTimeSerieHeader3D
    private ::          AddScalarTimeSerie
    public  :: UpdateBoxDif !soffs


    !Selector
    public  :: GetBoxes
    public  :: GetDTBoxes
    public  :: GetNumberOfBoxes
    public  :: CheckIfInsideBox
    public  :: GetIfBoxInsideDomain
    public  :: UnGetBoxDif
    
    private ::          ReadLockGridInformation
    private ::          ReadUnLockGridInformation
    
    !Modifier
    public  :: BoxDif
    private ::      OutputFluxesTimeSerie2D
    private ::      OutputFluxesTimeSerie3D
    private ::      OutputScalarTimeSerie

    !Destructor
    public  :: KillBoxDif                                                     
    

    !Management
    private ::      Ready
    private ::          LocateObjBoxDif 
    
    !Interfaces----------------------------------------------------------------
    private :: StartBoxDif2D
    private :: StartBoxDif3D
    interface  StartBoxDif
        module procedure StartBoxDif2D
        module procedure StartBoxDif3D
    end interface  StartBoxDif

    private :: ConvertBoxesToMap2D_2D
    private :: ConvertBoxesToMap2D_3D
    interface  ConvertBoxesToMap2D
        module procedure ConvertBoxesToMap2D_2D
        module procedure ConvertBoxesToMap2D_3D
    end interface  ConvertBoxesToMap2D
    
    private :: BoxDifFluxes3D
    private :: BoxDifFluxes2D
    private :: BoxDifScalar3D_R4
    private :: BoxDifScalar3D_R8
    private :: BoxDifScalar2D_R4
    private :: BoxDifScalar2D_R8
    interface  BoxDif
        module procedure BoxDifFluxes3D
        module procedure BoxDifFluxes2D
        module procedure BoxDifScalar3D_R4
        module procedure BoxDifScalar3D_R8
        module procedure BoxDifScalar2D_R4
        module procedure BoxDifScalar2D_R8
    end interface  BoxDif

    private :: GetBoxes2D
    private :: GetBoxes3D
    interface  GetBoxes
        module procedure GetBoxes2D
        module procedure GetBoxes3D
    end interface  GetBoxes

    private :: UnGetBoxDif3D_I
    private :: UnGetBoxDif2D_I
    private :: UnGetBoxDif3D_R8
    private :: UnGetBoxDif2D_R8
    interface  UnGetBoxDif
        module procedure UnGetBoxDif3D_I
        module procedure UnGetBoxDif2D_I
        module procedure UnGetBoxDif3D_R8
        module procedure UnGetBoxDif2D_R8
    end interface  UnGetBoxDif
    
    !Parameter-----------------------------------------------------------------
    integer, parameter                                      :: Grid_Coordinates     = 1
    integer, parameter                                      :: Geo_Coordinates      = 2


    !Types---------------------------------------------------------------------
    
    private :: T_ID
    type       T_ID
        integer                                             :: Number               = 0
        character(LEN = StringLength)                       :: Name
    end type   T_ID                                                 
    

    private :: T_Box
    type       T_Box
        integer                                             :: MainID
        type(T_ID),      dimension(:  ), pointer            :: ID
        integer                                             :: NumberOfLayers
        integer,         dimension(:,:), pointer            :: Layers
        type(T_Polygon),                 pointer            :: Polygon
        type(T_Box),                     pointer            :: Next
        logical                                             :: InsideDomain
    end type   T_Box                                        
                                                            
                                                            
    private :: T_External                                   
    type       T_External                                   
        type(T_Time)                                        :: Now
        real,    dimension(:,:  ), pointer                  :: ZCoordX, ZCoordY
        real,    dimension(:,:  ), pointer                  :: CornersX, CornersY
        real(8), dimension(:,:  ), pointer                  :: FluxX2D
        real(8), dimension(:,:  ), pointer                  :: FluxY2D
        real(8), dimension(:,:,:), pointer                  :: FluxX3D
        real(8), dimension(:,:,:), pointer                  :: FluxY3D
        real(8), dimension(:,:,:), pointer                  :: FluxZ3D
        real(8), dimension(:,:  ), pointer                  :: Scalar2D_R8
        real(4), dimension(:,:  ), pointer                  :: Scalar2D_R4
        real(8), dimension(:,:,:), pointer                  :: Scalar3D_R8
        real(4), dimension(:,:,:), pointer                  :: Scalar3D_R4
        integer, dimension(:,:  ), pointer                  :: WaterPoints2D
        integer, dimension(:,:,:), pointer                  :: WaterPoints3D
        character(LEN = StringLength)                       :: CurrentTimeSerieName
    end type T_External


    private :: T_BoxTimeSerie
    type       T_BoxTimeSerie
        type(T_ID)                                          :: ID
        integer                                             :: ObjTimeSerie         = 0
        type(T_BoxTimeSerie), pointer                       :: Next
    end type   T_BoxTimeSerie

    private :: T_BoxDif
    type       T_BoxDif
        integer                                             :: InstanceID
        character(LEN = PathLength)                         :: BoxesFilePath
        type(T_Box),                  pointer               :: FirstBox
        integer                                             :: NumberOfBoxes2D      = 0
        integer                                             :: NumberOfBoxes3D      = 0
        real                                                :: DT
        integer, dimension(:, :   ),  pointer               :: Boxes2D
        integer, dimension(:, :, :),  pointer               :: Boxes3D
        real,    dimension(:, :   ),  pointer               :: Fluxes2D
        real,    dimension(:, :   ),  pointer               :: Fluxes3D
        real,    dimension(:      ),  pointer               :: FluxesTimeSerieLine
        real,    dimension(:      ),  pointer               :: ScalarTimeSerieLine
        integer, dimension(:,:    ),  pointer               :: AdjacentBoxesBoundaries2D
        integer, dimension(:,:    ),  pointer               :: AdjacentBoxesBoundaries3D
        integer                                             :: nAdjacentBoxesBoundaries2D = 0
        integer                                             :: nAdjacentBoxesBoundaries3D = 0
        logical, dimension(:, :, :),  pointer               :: BoundaryFace3DX
        logical, dimension(:, :, :),  pointer               :: BoundaryFace3DY
        logical, dimension(:, :, :),  pointer               :: BoundaryFace3DZ
        logical, dimension(:, :   ),  pointer               :: BoundaryFace2DX
        logical, dimension(:, :   ),  pointer               :: BoundaryFace2DY
        type(T_External)                                    :: ExternalVar
        type(T_Size2D)                                      :: Size2D
        type(T_Size2D)                                      :: WorkSize2D
        type(T_Size3D)                                      :: Size3D
        type(T_Size3D)                                      :: WorkSize3D
        character(len=StringLength), dimension(:), pointer  :: FluxesTimeSerieHeader
        character(len=StringLength), dimension(:), pointer  :: ScalarTimeSerieHeader
        type(T_BoxTimeSerie), pointer                       :: FirstFluxesTimeSerie
        type(T_BoxTimeSerie), pointer                       :: FirstScalarTimeSerie
        integer                                             :: CoordinateType, GridType
        logical                                             :: WriteBoxes
        character(LEN = PathLength)                         :: BoxesOutputFile
        integer                                             :: NextBoxMainID        = 1
        integer                                             :: NextBoxID            = 1
        integer                                             :: ObjTime              = 0
        integer                                             :: ObjHorizontalGrid    = 0
        integer                                             :: ObjEnterData         = 0
        type(T_BoxDif), pointer                             :: Next
    end type  T_BoxDif

    !Global Module Variables
    type (T_BoxDif), pointer                                :: FirstObjBoxDif
    type (T_BoxDif), pointer                                :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartBoxDif2D(BoxDifID, TimeID, HorizontalGridID,                &
                             BoxesFilePath, FluxesOutputList, ScalarOutputList, &
                             WaterPoints2D, STAT)
                           
        !Arguments---------------------------------------------------------------
        integer                                                 :: BoxDifID
        integer                                                 :: TimeID
        integer                                                 :: HorizontalGridID
        character(LEN = *),           intent(IN )               :: BoxesFilePath
        character(LEN = *), optional, pointer, dimension(:  )   :: FluxesOutputList
        character(LEN = *), optional, pointer, dimension(:  )   :: ScalarOutputList
        integer,                      pointer, dimension(:,:)   :: WaterPoints2D
        integer,            optional, intent(OUT)               :: STAT     

        !External----------------------------------------------------------------
        integer                                             :: ready_, STAT_CALL     

        !Local-------------------------------------------------------------------
        integer                                             :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBoxDif_)) then
            nullify (FirstObjBoxDif)
            call RegisterModule (mBoxDif_) 
        endif

        call Ready(BoxDifID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            
            Me%BoxesFilePath     = trim(BoxesFilePath)
            Me%ExternalVar%WaterPoints2D => WaterPoints2D

            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                            &
                                       Size        = Me%Size2D,                         &
                                       WorkSize    = Me%WorkSize2D,                     &
                                       STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartBoxDif2D - ModuleBoxDif - ERR01'


            call ConstructBoxes2D

            if(present(FluxesOutputList)) call ConstructOutputFluxes2D  (FluxesOutputList)
                
            if(present(ScalarOutputList)) call ConstructOutputScalar2D   (ScalarOutputList)

            !Returns ID
            BoxDifID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleBoxDif - StartBoxDif2D - ERR99' 

        end if cd0



        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartBoxDif2D
    
    !--------------------------------------------------------------------------
    
    subroutine StartBoxDif3D(BoxDifID, TimeID, HorizontalGridID,                &
                             BoxesFilePath, FluxesOutputList, ScalarOutputList, &
                             WaterPoints3D, Size3D, WorkSize3D, STAT)
                           
        !Arguments---------------------------------------------------------------
        integer                                                 :: BoxDifID
        integer                                                 :: TimeID
        integer                                                 :: HorizontalGridID
        character(LEN = *),           intent(IN )               :: BoxesFilePath
        character(LEN = *), optional, pointer, dimension(:)     :: FluxesOutputList
        character(LEN = *), optional, pointer, dimension(:)     :: ScalarOutputList
        integer,                      pointer, dimension(:,:,:) :: WaterPoints3D
        type(T_Size3D)                                          :: Size3D, WorkSize3D
        integer,            optional, intent(OUT)               :: STAT     

        !External----------------------------------------------------------------
        integer                                             :: ready_, STAT_CALL     

        !Local-------------------------------------------------------------------
        integer                                             :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBoxDif_)) then
            nullify (FirstObjBoxDif)
            call RegisterModule (mBoxDif_) 
        endif

        call Ready(BoxDifID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           TimeID          )
            Me%ObjHorizontalGrid = AssociateInstance (mHORIZONTALGRID_, HorizontalGridID)
            
            Me%ExternalVar%WaterPoints3D => WaterPoints3D
            Me%BoxesFilePath             = trim(BoxesFilePath)
            Me%Size3D                    = Size3D  
            Me%WorkSize3D                = WorkSize3D

            call GetHorizontalGridSize(Me%ObjHorizontalGrid,                            &
                                       Size        = Me%Size2D,                         &
                                       WorkSize    = Me%WorkSize2D,                     &
                                       STAT        = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartBoxDif3D - ModuleBoxDif - ERR02'


            call ConstructBoxes3D

            if(present(FluxesOutputList)) call ConstructOutputFluxes3D  (FluxesOutputList)
                
            if(present(ScalarOutputList)) call ConstructOutputScalar3D  (ScalarOutputList)

            !Returns ID
            BoxDifID = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleBoxDif - StartBoxDif3D - ERR99' 

        end if cd0



        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine StartBoxDif3D

    !--------------------------------------------------------------------------

    subroutine UpdateBoxDif (BoxDifID, NewFluxesOutputList, NewScalarOutputList, nDimensions,STAT)
                           
        !Arguments---------------------------------------------------------------
        integer                                                 :: BoxDifID
        character(LEN = *), optional, pointer, dimension(:  )   :: NewFluxesOutputList
        character(LEN = *), optional, pointer, dimension(:  )   :: NewScalarOutputList
        integer, intent(IN)                                     :: nDimensions
        integer,            optional, intent(OUT)               :: STAT     

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                                 :: STAT_, ready_

        !------------------------------------------------------------------------
        
        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            STAT_ = UNKNOWN_
            
            if (nDimensions .eq. 2) then
                call ConstructOutputScalar2D (NewScalarOutputList)
                call ConstructOutputFluxes2D (NewFluxesOutputList)
            else 
                call ConstructOutputScalar3D (NewScalarOutputList)
                call ConstructOutputFluxes3D (NewFluxesOutputList)
            end if
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UpdateBoxDif
    
    !----------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_BoxDif), pointer                         :: NewObjBoxDif
        type (T_BoxDif), pointer                         :: PreviousObjBoxDif


        !Allocates new instance
        allocate (NewObjBoxDif)
        nullify  (NewObjBoxDif%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjBoxDif)) then
            FirstObjBoxDif          => NewObjBoxDif
            Me                      => NewObjBoxDif
        else
            PreviousObjBoxDif       => FirstObjBoxDif
            Me                      => FirstObjBoxDif%Next
            do while (associated(Me))
                PreviousObjBoxDif   => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjBoxDif
            PreviousObjBoxDif%Next  => NewObjBoxDif
        endif

        Me%InstanceID       = RegisterNewInstance (mBOXDIF_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ConstructBoxes2D 
        
        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL      
     
        !Local-----------------------------------------------------------------

        call ConstructEnterData (Me%ObjEnterData, Me%BoxesFilePath, STAT = STAT_CALL)

        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBoxes - ModuleBoxDif - ERR02'

        call ReadGlobalOptions

        call ReadLockGridInformation
        
        call ReadBoxes2D

        call AllocateVariables2D

        call InitializeEnvironmentBox2D

        call ConvertBoxesToMap2D(Me%ExternalVar%WaterPoints2D)

        call FindAdjacentBoxesBoundaries2D

        if(Me%WriteBoxes) call WriteBoxes

        call ReadUnLockGridInformation

        call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBoxes - ModuleBoxDif - ERR03'


    end subroutine ConstructBoxes2D
    
    !--------------------------------------------------------------------------

    subroutine ConstructBoxes3D 
        
        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL      
     
        !Local-----------------------------------------------------------------


        call ConstructEnterData (Me%ObjEnterData, Me%BoxesFilePath, STAT = STAT_CALL)

        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBoxes - ModuleBoxDif - ERR02'

        call ReadGlobalOptions

        call ReadLockGridInformation
        
        call ReadBoxes3D

        call AllocateVariables2D

        call AllocateVariables3D

        call InitializeEnvironmentBox3D

        call ConvertBoxesToMap2D(Me%ExternalVar%WaterPoints3D)

        call ConvertBoxesToMap3D

        call FindAdjacentBoxesBoundaries3D

        if(Me%WriteBoxes) call WriteBoxes

        call ReadUnLockGridInformation

        call KillEnterData      (Me%ObjEnterData, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBoxes - ModuleBoxDif - ERR03'


    end subroutine ConstructBoxes3D
    
    !--------------------------------------------------------------------------


    subroutine ReadGlobalOptions 
        
        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                     :: iflag, STAT_CALL       
                      
        !Local-----------------------------------------------------------------

        call GetData(Me%CoordinateType,                                         &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      ='TYPE',                                      &
                     Default      = Grid_Coordinates,                           &
                     ClientModule = MohidModules(mBOXDIF_)%Name,                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadGlobalOptions - ModuleBoxDif - ERR10'
        
        call GetData(Me%WriteBoxes,                                             &
                     Me%ObjEnterData, iflag,                                    &
                     SearchType   = FromFile,                                   &
                     keyword      ='WRITE_BOXES',                               &
                     Default      = OFF,                                        &
                     ClientModule = MohidModules(mBOXDIF_)%Name,                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadGlobalOptions - ModuleBoxDif - ERR20'

        call GetData(Me%DT,                                 &
                     Me%ObjEnterData, iflag,                &
                     SearchType   = FromFile,               &
                     keyword      = 'DT_BOXES',             &
                     Default      = 300.,                  &
                     ClientModule = 'ModuleBoxDif',        &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadGlobalOptions - ModuleBoxDif - ERR40'

        if(Me%WriteBoxes)then

            call GetData(Me%BoxesOutputFile,                                    &
                         Me%ObjEnterData, iflag,                                &
                         SearchType   = FromFile,                               &
                         keyword      ='OUTPUT_FILE',                           &
                         ClientModule = MohidModules(mBOXDIF_)%Name,            &
                         STAT         = STAT_CALL)        
            if(STAT_CALL .ne. SUCCESS_)stop 'ReadGlobalOptions - ModuleBoxDif - ERR30'

            if(iflag==0)then
                write(*,*)'Option WRITE_BOXES is activated in boxes file'
                write(*,*)trim(Me%BoxesFilePath)
                write(*,*)'Please define keyword OUTPUT_FILE'
                stop 'ReadGlobalOptions - ModuleBoxDif - ERR40'
            endif

        end if

    end subroutine ReadGlobalOptions
    
    !--------------------------------------------------------------------------

    subroutine CountNumberOfBoxes2D

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL       
        integer                                         :: ClientNumber
        logical                                         :: BlockFound

        !Local-----------------------------------------------------------------
        integer                                         :: NumberOfBoxes2D
        
        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                                        ClientNumber,                                &
                                        '<beginpolygon>',                            &
                                        '<endpolygon>',                              &
                                        BlockFound,                                  &
                                        STAT = STAT_CALL)

cd1 :           if      (STAT_CALL .EQ. SUCCESS_ ) then    

cd2 :               if (BlockFound) then

                        !BlockInBlock
                        call RewindBlock(Me%ObjEnterData, ClientNumber = ClientNumber, STAT = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)stop 'CountNumberOfBoxes2D - ModuleBoxDif - ERR01'

                        NumberOfBoxes2D = NumberOfBoxes2D + 1
                    else

                        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                        if(STAT_CALL .ne. SUCCESS_)stop 'CountNumberOfBoxes2D - ModuleBoxDif - ERR04'

                        exit do1    !No more blocks

                    end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1

                if(STAT_CALL .ne. SUCCESS_)stop 'CountNumberOfBoxes2D - ModuleBoxDif - ERR05'

            end if cd1
        end do do1


    end subroutine CountNumberOfBoxes2D

    !--------------------------------------------------------------------------

    subroutine CountNumberOfLayers (ClientNumber, NumberOfLayers)

        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber 
        integer,         intent(out)                    :: NumberOfLayers

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL       
        logical                                         :: BlockInBlockFound

        !Local-----------------------------------------------------------------
        integer                                         :: Line, FirstLine, LastLine
        
        !Begin-----------------------------------------------------------------
        
        NumberOfLayers = 0

        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,        &
                                   '<<beginverticallayer>>',             &
                                   '<<endverticallayer>>',               &
                                   BlockInBlockFound,                    &
                                   FirstLine = FirstLine,                &
                                   LastLine  = LastLine,                 &
                                   STAT      = STAT_CALL)
        if      (STAT_CALL .EQ. SUCCESS_) then    
            
            if (BlockInBlockFound) then
            
                if     (((LastLine + 1) - (FirstLine - 1)) .GE. 1) then
                    
                    do Line = FirstLine + 1, LastLine - 1    
                        NumberOfLayers = NumberOfLayers + 1
                    end do

                    call RewindBlock(Me%ObjEnterData, ClientNumber = ClientNumber, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)stop 'CountNumberOfLayers - ModuleBoxDif - ERR01'

                elseif (((LastLine + 1) - (FirstLine - 1)) .EQ. 0) then 
                    
                    NumberOfLayers = NumberOfLayers + 1

                else
                    write(*,*)  
                    write(*,*) 'Error counting boxes. '
                    if(STAT_CALL .ne. SUCCESS_)stop 'CountNumberOfLayers - ModuleBoxDif - ERR02'
                end if 

            else
                !It assumes that at least one 3D box exixts
                if (NumberOfLayers .EQ. 0) &
                    NumberOfLayers = 1
            end if 

        else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            if(STAT_CALL .ne. SUCCESS_)stop 'CountNumberOfLayers - ModuleBoxDif - ERR03'
        end if 



    end subroutine CountNumberOfLayers
    
    !--------------------------------------------------------------------------
    
    subroutine ReadBoxes2D

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL       
        integer                                         :: ClientNumber
        logical                                         :: BlockFound, BlockInBlockFound

        !Local-----------------------------------------------------------------
        integer                                         :: FirstLine, LastLine
        type (T_Box),    pointer                        :: NewBox

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                                        ClientNumber,                                &
                                        '<beginpolygon>',                            &
                                        '<endpolygon>',                              &
                                        BlockFound,                                  &
                                        STAT = STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_ )then    

cd2 :           if (BlockFound) then

                    !BlockInBlock
                    call RewindBlock(Me%ObjEnterData, ClientNumber = ClientNumber, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes2D - ModuleBoxDif - ERR01'
                        
                    call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,        &
                                               '<<beginvertix>>',                    &
                                               '<<endvertix>>',                      &
                                               BlockInBlockFound,                    &
                                               FirstLine = FirstLine,                &
                                               LastLine  = LastLine,                 &
                                               STAT      = STAT_CALL)
                    if(STAT_CALL .EQ. SUCCESS_)then    
                    
                        if (BlockInBlockFound) then

                            call AddBox             (NewBox, 1)

                            call ConstructBox2D     (NewBox, FirstLine, LastLine)

                        end if

                    else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                        write(*,*)  
                        write(*,*) 'Error calling ExtractBlockFromBuffer. '
                        if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes2D - ModuleBoxDif - ERR03'
                    end if

                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes2D - ModuleBoxDif - ERR04'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1

                if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes2D - ModuleBoxDif - ERR05'

            end if cd1
        end do do1


    end subroutine ReadBoxes2D

    !--------------------------------------------------------------------------
    
    subroutine ReadBoxes3D

        !Arguments-------------------------------------------------------------

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL       
        integer                                         :: ClientNumber
        logical                                         :: BlockFound, BlockInBlockFound

        !Local-----------------------------------------------------------------
        integer                                         :: FirstLine, LastLine
        integer                                         :: NumberOfLayers
        type (T_Box),    pointer                        :: NewBox

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                             &
                                        ClientNumber,                                &
                                        '<beginpolygon>',                            &
                                        '<endpolygon>',                              &
                                        BlockFound,                                  &
                                        STAT = STAT_CALL)

cd1 :       if(STAT_CALL .EQ. SUCCESS_ )then    

cd2 :           if (BlockFound) then

                    !BlockInBlock
                    call RewindBlock(Me%ObjEnterData, ClientNumber = ClientNumber, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes3D - ModuleBoxDif - ERR01'
                        
                    call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,        &
                                               '<<beginvertix>>',                    &
                                               '<<endvertix>>',                      &
                                               BlockInBlockFound,                    &
                                               FirstLine = FirstLine,                &
                                               LastLine  = LastLine,                 &
                                               STAT      = STAT_CALL)
                    if(STAT_CALL .EQ. SUCCESS_)then    
                    
                        if (BlockInBlockFound) then

                            call CountNumberOfLayers(ClientNumber, NumberOfLayers)

                            call AddBox             (NewBox, NumberOfLayers)

                            call ConstructBox2D     (NewBox, FirstLine, LastLine)

                            call ConstructBox3D     (NewBox, ClientNumber)

                        end if

                    else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                        write(*,*)  
                        write(*,*) 'Error calling ExtractBlockFromBuffer. '
                        if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes3D - ModuleBoxDif - ERR03'
                    end if

                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes3D - ModuleBoxDif - ERR04'

                    exit do1    !No more blocks

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1

                if(STAT_CALL .ne. SUCCESS_)stop 'ReadBoxes3D - ModuleBoxDif - ERR05'

            end if cd1
        end do do1


    end subroutine ReadBoxes3D

    !--------------------------------------------------------------------------

    subroutine ConstructBox2D (NewBox, FirstLine, LastLine)

        !Arguments-------------------------------------------------------------
        type (T_Box),               pointer         :: NewBox
        integer                                     :: FirstLine, LastLine

        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag, FromBlock       
        
        !Local-----------------------------------------------------------------
        type (T_Polygon), pointer                   :: ModelDomainPolygon
        integer                                     :: i, Line, Ic, Jc
        real,    dimension(:),   pointer            :: VertixLocation
        real                                        :: X, Y


        !Begin-----------------------------------------------------------------


        !Always has one vertix more(last equals the first)             
        NewBox%Polygon%Count = LastLine - FirstLine

        if(NewBox%Polygon%Count .le. 2)then
            write(*,*)
            write(*,*) 'No valid vertix block found for box'
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBox2D - ModuleBoxDif - ERR10'
        end if

        allocate(NewBox%Polygon%Vertices (1:NewBox%Polygon%Count))
        allocate(NewBox%Polygon%VerticesF(1:NewBox%Polygon%Count))

        allocate(VertixLocation (1:2))

        i = 1

        do Line = FirstLine + 1, LastLine - 1
            
            call GetExtractType(FromBlock = FromBlock)

            call GetData(VertixLocation,                                 &
                         EnterDataID = Me%ObjEnterData,                  &
                         flag        = iflag,                            &
                         SearchType  = FromBlock,                        &
                         Buffer_Line = Line,                             &
                         STAT        = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBox2D - ModuleBoxDif - ERR20'

            select case (Me%CoordinateType)

                case (Grid_Coordinates)
                    
                    NewBox%Polygon%Vertices(i)%J = int(VertixLocation(1))
                    NewBox%Polygon%Vertices(i)%I = int(VertixLocation(2))

                       
                    Ic = NewBox%Polygon%Vertices(i)%I
                    Jc = NewBox%Polygon%Vertices(i)%J

                    X = Me%ExternalVar%ZCoordX(Ic, Jc)
                    Y = Me%ExternalVar%ZCoordY(Ic, Jc)
                        
                        
                    if (X<FillValueReal/4. .or.  Y<FillValueReal/4.) stop 'ConstructBox2D - ModuleBoxDif - ERR30'
                 

                    NewBox%Polygon%VerticesF(i)%X = X
                    NewBox%Polygon%VerticesF(i)%Y = Y

                case (Geo_Coordinates)

                    NewBox%Polygon%VerticesF(i)%X = VertixLocation(1)
                    NewBox%Polygon%VerticesF(i)%Y = VertixLocation(2)

            end select

            i = i + 1

        end do
        
        !Box has always one vertix more(last equals the first)
        NewBox%Polygon%Vertices (NewBox%Polygon%Count)%I = NewBox%Polygon%Vertices (1)%I
        NewBox%Polygon%Vertices (NewBox%Polygon%Count)%J = NewBox%Polygon%Vertices (1)%J
        NewBox%Polygon%VerticesF(NewBox%Polygon%Count)%X = NewBox%Polygon%VerticesF(1)%X
        NewBox%Polygon%VerticesF(NewBox%Polygon%Count)%Y = NewBox%Polygon%VerticesF(1)%Y

        
        call SetLimitsPolygon(NewBox%Polygon)

        call GetGridBorderPolygon(Me%ObjHorizontalGrid, ModelDomainPolygon, STAT = STAT_CALL)
        
        if(STAT_CALL /= SUCCESS_)stop 'ConstructBox2D - ModuleBoxDif - ERR30'

        NewBox%InsideDomain = VertPolygonInsidePolygon(NewBox%Polygon, ModelDomainPolygon)

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid, ModelDomainPolygon, STAT = STAT_CALL)
        
        if(STAT_CALL /= SUCCESS_)stop 'ConstructBox2D - ModuleBoxDif - ERR50'


        deallocate(VertixLocation)


    end subroutine ConstructBox2D


    !--------------------------------------------------------------------------

    subroutine ConstructBox3D (NewBox, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Box),    pointer                        :: NewBox
        integer, intent(in)                             :: ClientNumber 

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL       
        logical                                         :: BlockInBlockFound
        integer                                         :: iflag, FromBlock

        !Local-----------------------------------------------------------------
        integer                                         :: i, Line, FirstLine, LastLine
        integer, dimension(:), pointer                  :: LayersLimits
        !Begin-----------------------------------------------------------------
        
        allocate(LayersLimits (1:2)                         )
        allocate(NewBox%Layers(1:NewBox%NumberOfLayers, 1:2))

        call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,               &
                                   '<<beginverticallayer>>',                    &
                                   '<<endverticallayer>>',                      &
                                   BlockInBlockFound,                           &
                                   FirstLine = FirstLine,                       &
                                   LastLine  = LastLine,                        &
                                   STAT      = STAT_CALL)
        if      (STAT_CALL .EQ. SUCCESS_) then    
            
            if (BlockInBlockFound) then
            
                if     (((LastLine + 1) - (FirstLine - 1)) .GE. 1) then
                
                    i = 1

                    do Line = FirstLine + 1, LastLine - 1
            
                        call GetExtractType(FromBlock = FromBlock)

                        call GetData(LayersLimits,                                   &
                                     EnterDataID = Me%ObjEnterData,                  &
                                     flag        = iflag,                            &
                                     SearchType  = FromBlock,                        &
                                     Buffer_Line = Line,                             &
                                     STAT        = STAT_CALL)
                        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBox3D - ModuleBoxDif - ERR01'
                        
                        NewBox%Layers(i,1) = LayersLimits(1)
                        NewBox%Layers(i,2) = LayersLimits(2)

                        i = i + 1

                    end do

                    elseif (((LastLine + 1) - (FirstLine - 1)) .EQ. 0) then 
                        
                        NewBox%Layers(1,1) = Me%WorkSize3D%KLB
                        NewBox%Layers(1,2) = Me%WorkSize3D%KUB
                     
                    else
                        write(*,*)  
                        write(*,*) 'Error counting boxes. '
                        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBox3D - ModuleBoxDif - ERR02'
                    end if 

                else

                    NewBox%Layers(1,1) = Me%WorkSize3D%KLB
                    NewBox%Layers(1,2) = Me%WorkSize3D%KUB

                end if 

        else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
            write(*,*)  
            write(*,*) 'Error calling ExtractBlockFromBuffer. '
            if(STAT_CALL .ne. SUCCESS_)stop 'ConstructBox3D - ModuleBoxDif - ERR03'
        end if 



    end subroutine ConstructBox3D
    
    
    !--------------------------------------------------------------------------
    
    
    subroutine ReadLockGridInformation
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL       

        !Begin-----------------------------------------------------------------

        call GetZCoordinates(Me%ObjHorizontalGrid, Me%ExternalVar%ZCoordX, Me%ExternalVar%ZCoordY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockGridInformation - ModuleBoxDif - ERR10'

        call GetCornersCoordinates(Me%ObjHorizontalGrid, Me%ExternalVar%CornersX, Me%ExternalVar%CornersY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLockGridInformation - ModuleBoxDif - ERR10'

    end subroutine ReadLockGridInformation
    
    
    !--------------------------------------------------------------------------

    subroutine ReadUnLockGridInformation 

        
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL       
                      
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%ZCoordX,                                &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockGridInformation - ModuleBoxDif - ERR10'
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%ZCoordY,                                &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockGridInformation - ModuleBoxDif - ERR20'

        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%CornersX,                                &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockGridInformation - ModuleBoxDif - ERR30'
        
        call UnGetHorizontalGrid(Me%ObjHorizontalGrid,                                  &
                                 Me%ExternalVar%CornersY,                                &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadUnLockGridInformation - ModuleBoxDif - ERR40'


    end subroutine ReadUnLockGridInformation
    
    
    !--------------------------------------------------------------------------
    
    
    subroutine ConvertBoxesToMap2D_2D(WaterPoints2D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:), pointer                :: WaterPoints2D


        !Local-----------------------------------------------------------------
        type (T_Box),    pointer                        :: CurrentBox
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB
        integer                                         :: i, j, ii, jj
        type(T_PointF), pointer                         :: Point
        logical                                         :: NoPointInsidePolygon =.true.

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize2D%ILB
        WIUB = Me%WorkSize2D%IUB
        WJLB = Me%WorkSize2D%JLB
        WJUB = Me%WorkSize2D%JUB


        allocate(Point)

        CurrentBox => Me%FirstBox

        do while(associated(CurrentBox))
                
            NoPointInsidePolygon  =.true.                
                    
            do i = WILB,  WIUB
            do j = WJLB , WJUB

                if(WaterPoints2D(i,j).eq. WaterPoint)then

                    Point%X = Me%ExternalVar%ZCoordX(I, J  )
                    Point%Y = Me%ExternalVar%ZCoordY(I, J  )
                    
                    if(IsPointInsidePolygon(Point, CurrentBox%Polygon))then

                        Me%Boxes2D(i,j) = CurrentBox%MainID
                        
                         NoPointInsidePolygon  =.false.

                    end if
                end if
            enddo
            enddo
            
            if (NoPointInsidePolygon) then
            
                do i = WILB,  WIUB
                do j = WJLB , WJUB

                    if(WaterPoints2D(i,j).eq. WaterPoint)then
                        do ii=i,i+1
                        do jj=j,j+1
    
                            Point%X = Me%ExternalVar%CornersX(ii,jj)
                            Point%Y = Me%ExternalVar%CornersY(ii,jj)
                        
                            if(IsPointInsidePolygon(Point, CurrentBox%Polygon))then
                                Me%Boxes2D(i,j) = CurrentBox%MainID
                            end if
                            
                        enddo
                        enddo
                    end if
                enddo
                enddo            
            
            endif        

            CurrentBox => CurrentBox%Next
            
        end do

        deallocate(Point)


    end subroutine ConvertBoxesToMap2D_2D


   !--------------------------------------------------------------------------
    
    
    subroutine ConvertBoxesToMap2D_3D(WaterPoints3D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:,:,:) , pointer             :: WaterPoints3D

        !Local-----------------------------------------------------------------
        type (T_Box),    pointer                        :: CurrentBox
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB, WKUB
        integer                                         :: i, j, ii, jj
        logical                                         :: NoPointInsidePolygon
        type(T_PointF), pointer                         :: Point

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize3D%ILB
        WIUB = Me%WorkSize3D%IUB
        WJLB = Me%WorkSize3D%JLB
        WJUB = Me%WorkSize3D%JUB
        WKUB = Me%WorkSize3D%KUB


        allocate(Point)

        CurrentBox => Me%FirstBox

        do while(associated(CurrentBox))
                
            NoPointInsidePolygon  =.true.                
                    
            do i = WILB,  WIUB
            do j = WJLB , WJUB

                if(WaterPoints3D(i,j, WKUB).eq. WaterPoint)then

                    Point%X = Me%ExternalVar%ZCoordX(I, J  )
                    Point%Y = Me%ExternalVar%ZCoordY(I, J  )
                    
                    if(IsPointInsidePolygon(Point, CurrentBox%Polygon))then

                        Me%Boxes2D(i,j) = CurrentBox%MainID
                        
                         NoPointInsidePolygon  =.false.

                    end if
                end if
            enddo
            enddo
            
            if (NoPointInsidePolygon) then
            
                do i = WILB,  WIUB
                do j = WJLB , WJUB

                    if(WaterPoints3D(i,j, WKUB).eq. WaterPoint)then
                        do ii=i,i+1
                        do jj=j,j+1
    
                            Point%X = Me%ExternalVar%CornersX(ii,jj)
                            Point%Y = Me%ExternalVar%CornersY(ii,jj)
                        
                            if(IsPointInsidePolygon(Point, CurrentBox%Polygon))then
                                Me%Boxes2D(i,j) = CurrentBox%MainID
                            end if
                            
                        enddo
                        enddo
                    end if
                enddo
                enddo            
            
            endif        

            CurrentBox => CurrentBox%Next
            
        end do

        deallocate(Point)


    end subroutine ConvertBoxesToMap2D_3D

    !--------------------------------------------------------------------------


    subroutine ConvertBoxesToMap3D

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Box),    pointer                        :: CurrentBox
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB
        integer                                         :: WKLB, WKUB
        integer                                         :: i, j, k, Layer
        integer                                         :: bottomlayer, toplayer

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize3D%ILB
        WIUB = Me%WorkSize3D%IUB
        WJLB = Me%WorkSize3D%JLB
        WJUB = Me%WorkSize3D%JUB
        WKLB = Me%WorkSize3D%KLB
        WKUB = Me%WorkSize3D%KUB


        do k = WKLB, WKUB
        do i = WILB, WIUB
        do j = WJLB, WJUB

            CurrentBox => Me%FirstBox

            do while(associated(CurrentBox))

                if(Me%Boxes2D(i,j) .eq. CurrentBox%MainID)then

                    do Layer = 1, CurrentBox%NumberOfLayers
                    
                        if (CurrentBox%Layers(Layer,1) < CurrentBox%Layers(Layer,2) ) then
                            bottomlayer = CurrentBox%Layers(Layer,1)
                            toplayer = CurrentBox%Layers(Layer,2)
                        else
                            bottomlayer = CurrentBox%Layers(Layer,2)
                            toplayer = CurrentBox%Layers(Layer,1)                                                    
                        end if

                        if(k >= bottomlayer .and. k <= toplayer)then
                            Me%Boxes3D(i,j,k) = CurrentBox%ID(Layer)%Number
                        end if

                    end do

                end if

                CurrentBox => CurrentBox%Next

            end do

        enddo
        enddo
        enddo



    end subroutine ConvertBoxesToMap3D

    !--------------------------------------------------------------------------

    subroutine InitializeEnvironmentBox2D

        !Local-----------------------------------------------------------------
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB
        integer                                         :: i, j

        !Begin-----------------------------------------------------------------


        WILB = Me%WorkSize2D%ILB
        WIUB = Me%WorkSize2D%IUB
        WJLB = Me%WorkSize2D%JLB
        WJUB = Me%WorkSize2D%JUB

        do i = WILB, WIUB
        do j = WJLB, WJUB

            if((Me%ExternalVar%WaterPoints2D(i,j) .eq. WaterPoint) .and. &
               (Me%Boxes2D(i,j) < 0))then

                Me%Boxes2D(i,j) = 0
            else
                Me%Boxes2D(i,j) = null_int
            end if


        enddo
        enddo


    end subroutine InitializeEnvironmentBox2D
    
    !--------------------------------------------------------------------------

    subroutine InitializeEnvironmentBox3D

        !Local-----------------------------------------------------------------
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB
        integer                                         :: WKLB, WKUB
        integer                                         :: i, j, k

        !Begin-----------------------------------------------------------------


        WILB = Me%WorkSize3D%ILB
        WIUB = Me%WorkSize3D%IUB
        WJLB = Me%WorkSize3D%JLB
        WJUB = Me%WorkSize3D%JUB
        WKLB = Me%WorkSize3D%KLB
        WKUB = Me%WorkSize3D%KUB

        do k = WKLB, WKUB
        do i = WILB, WIUB
        do j = WJLB, WJUB

            if((Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint) .and. &
               (Me%Boxes3D(i,j,k) < 0))then

                Me%Boxes2D(i,j  ) = 0
                Me%Boxes3D(i,j,k) = 0
            else
                Me%Boxes2D(i,j  ) = null_int
                Me%Boxes3D(i,j,k) = null_int
            end if


        enddo
        enddo
        enddo


    end subroutine InitializeEnvironmentBox3D
    
    !--------------------------------------------------------------------------
    
    subroutine AllocateVariables2D

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB
        integer                                         :: JLB, JUB

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size2D%ILB
        IUB = Me%Size2D%IUB
        JLB = Me%Size2D%JLB
        JUB = Me%Size2D%JUB

        nullify (Me%Boxes2D)
        allocate(Me%Boxes2D(ILB:IUB,JLB:JUB))
        Me%Boxes2D(:,:  )  = null_int
        
        nullify (Me%Fluxes2D)
        allocate(Me%Fluxes2D(0:Me%NumberOfBoxes2D, 0:Me%NumberOfBoxes2D))
        Me%Fluxes2D = 0.

        nullify (Me%AdjacentBoxesBoundaries2D)
        allocate(Me%AdjacentBoxesBoundaries2D(0:Me%NumberOfBoxes2D, 0:Me%NumberOfBoxes2D))
        Me%AdjacentBoxesBoundaries2D(:,:) = 0

        nullify (Me%BoundaryFace2DX)
        allocate(Me%BoundaryFace2DX(ILB:IUB,JLB:JUB         ))
        Me%BoundaryFace2DX(:,:  )  = OFF

        nullify (Me%BoundaryFace2DY)
        allocate(Me%BoundaryFace2DY(ILB:IUB,JLB:JUB         ))
        Me%BoundaryFace2DY(:,:  )  = OFF

    end subroutine AllocateVariables2D
    
    !--------------------------------------------------------------------------

    subroutine AllocateVariables3D

        !Local-----------------------------------------------------------------
        integer                                         :: ILB, IUB
        integer                                         :: JLB, JUB
        integer                                         :: KLB, KUB

        !Begin-----------------------------------------------------------------
        
        ILB = Me%Size3D%ILB
        IUB = Me%Size3D%IUB
        JLB = Me%Size3D%JLB
        JUB = Me%Size3D%JUB
        KLB = Me%Size3D%KLB
        KUB = Me%Size3D%KUB
        
        nullify (Me%Boxes3D)
        allocate(Me%Boxes3D(ILB:IUB,JLB:JUB, KLB:KUB))
        Me%Boxes3D(:,:,:)  = null_int
        
        nullify (Me%Fluxes3D)
        allocate(Me%Fluxes3D(0:Me%NumberOfBoxes3D, 0:Me%NumberOfBoxes3D))
        Me%Fluxes3D = 0.

        nullify (Me%AdjacentBoxesBoundaries3D)
        allocate(Me%AdjacentBoxesBoundaries3D(0:Me%NumberOfBoxes3D, 0:Me%NumberOfBoxes3D))
        Me%AdjacentBoxesBoundaries3D(:,:) = 0

        nullify (Me%BoundaryFace3DX)
        allocate(Me%BoundaryFace3DX(ILB:IUB,JLB:JUB, KLB:KUB))
        Me%BoundaryFace3DX(:,:,:)  = OFF

        nullify (Me%BoundaryFace3DY)
        allocate(Me%BoundaryFace3DY(ILB:IUB,JLB:JUB, KLB:KUB))
        Me%BoundaryFace3DY(:,:,:)  = OFF

        nullify (Me%BoundaryFace3DZ)
        allocate(Me%BoundaryFace3DZ(ILB:IUB,JLB:JUB, KLB:KUB))
        Me%BoundaryFace3DZ(:,:,:)  = OFF

    end subroutine AllocateVariables3D
    
    !--------------------------------------------------------------------------

    subroutine FindAdjacentBoxesBoundaries2D
        
        !Local-----------------------------------------------------------------
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB
        integer                                         :: i, j

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize2D%ILB
        WIUB = Me%WorkSize2D%IUB
        WJLB = Me%WorkSize2D%JLB
        WJUB = Me%WorkSize2D%JUB
        
        do j = WJLB , WJUB
        do i = WILB , WIUB

            if(Me%Boxes2D(i,j) > -55 .and.                                &
               Me%ExternalVar%WaterPoints2D(i,j) .eq. WaterPoint)then

                !XX axis
                if ((Me%Boxes2D(i,j) .ne. Me%Boxes2D(i,j+1)) .and.    &
                    (Me%Boxes2D(i,j+1) .gt. -55)) then

                    Me%AdjacentBoxesBoundaries2D(Me%Boxes2D(i,j), Me%Boxes2D(i,j+1)) = Boundary
                    Me%BoundaryFace2DX(i,j)                                        = ON           

                end if

                !YY axis
                if ((Me%Boxes2D(i,  j) .ne. Me%Boxes2D(i+1,j)) .and.    &
                    (Me%Boxes2D(i+1,j) .gt. -55)) then

                    Me%AdjacentBoxesBoundaries2D(Me%Boxes2D(i,j), Me%Boxes2D(i+1,j)) = Boundary
                    Me%BoundaryFace2DY(i,j)                                        = ON           

                end if

            end if

        enddo
        enddo

        Me%nAdjacentBoxesBoundaries2D = sum(Me%AdjacentBoxesBoundaries2D)


    end subroutine FindAdjacentBoxesBoundaries2D

    !--------------------------------------------------------------------------

    subroutine FindAdjacentBoxesBoundaries3D
        
        !Local-----------------------------------------------------------------
        integer                                         :: WILB, WIUB
        integer                                         :: WJLB, WJUB
        integer                                         :: WKLB, WKUB
        integer                                         :: i, j, k

        !Begin-----------------------------------------------------------------
        
        WILB = Me%WorkSize3D%ILB
        WIUB = Me%WorkSize3D%IUB
        WJLB = Me%WorkSize3D%JLB
        WJUB = Me%WorkSize3D%JUB
        WKLB = Me%WorkSize3D%KLB
        WKUB = Me%WorkSize3D%KUB

        
        do k = WKLB,  WKUB
        do j = WJLB , WJUB
        do i = WILB,  WIUB

            if(Me%Boxes3D(i,j,k) > -55 .and.                                &
               Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then

                !XX axis
                if ((Me%Boxes3D(i,j,  k) .ne. Me%Boxes3D(i,j+1,k)) .and.    &
                    (Me%Boxes3D(i,j+1,k) .gt. -55)) then

                    Me%AdjacentBoxesBoundaries3D(Me%Boxes3D(i,j,k), Me%Boxes3D(i,j+1,k)) = Boundary
                    Me%BoundaryFace2DX(i,j  )                                          = ON           
                    Me%BoundaryFace3DX(i,j,k)                                          = ON

                end if

                !YY axis
                if ((Me%Boxes3D(i,  j,k) .ne. Me%Boxes3D(i+1,j,k)) .and.    &
                    (Me%Boxes3D(i+1,j,k) .gt. -55)) then

                    Me%AdjacentBoxesBoundaries3D(Me%Boxes3D(i,j,k), Me%Boxes3D(i+1,j,k)) = Boundary
                    Me%BoundaryFace2DY(i,j  )                                          = ON           
                    Me%BoundaryFace3DY(i,j,k)                                          = ON

                end if

                !ZZ axis
                if ((Me%Boxes3D(i,j,k  ) .ne. Me%Boxes3D(i,j,k+1)) .and.    &
                    (Me%Boxes3D(i,j,k+1) .gt. -55)) then

                    Me%AdjacentBoxesBoundaries3D(Me%Boxes3D(i,j,k), Me%Boxes3D(i,j,k+1)) = Boundary
                    Me%BoundaryFace3DZ(i,j,k)                                          = ON

                end if

            end if

        enddo
        enddo
        enddo

        Me%nAdjacentBoxesBoundaries3D = sum(Me%AdjacentBoxesBoundaries3D)


    end subroutine FindAdjacentBoxesBoundaries3D
    
    !--------------------------------------------------------------------------
    
    subroutine WriteBoxes

        !External-------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !Begin-----------------------------------------------------------------
        


        call WriteGridData (FileName            = Me%BoxesOutputFile,                   &
                            COMENT1             = "Defined 2D Boxes from file",         &
                            COMENT2             = Me%BoxesFilePath,                     &
                            HorizontalGridID    = Me%ObjHorizontalGrid,                 &
                            FillValue           = -99.,                                 &
                            Overwrite           = ON,                                   &
                            GridData2D_Int      = Me%Boxes2D,                           &
                            STAT                = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteBoxes - ModuleBoxDif - ERR10'


    end subroutine WriteBoxes


    !--------------------------------------------------------------------------


    subroutine AddBox (Box, NumberOfLayers)

        !Arguments-------------------------------------------------------------
        type (T_Box),    pointer                    :: Box
        integer                                     :: NumberOfLayers

        !Local-----------------------------------------------------------------
        type (T_Box),    pointer                    :: NewBox
        type (T_Box),    pointer                    :: PreviousBox
        integer                                     :: i
        character(len = StringLength)               :: Number

        
        !Allocates new instance
        allocate (NewBox)
        allocate (NewBox%Polygon)
        nullify  (NewBox%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%FirstBox)) then
            Me%FirstBox         => NewBox
            Box                 => NewBox
        else
            PreviousBox         => Me%FirstBox
            Box                 => Me%FirstBox%Next
            do while (associated(Box))
                PreviousBox     => Box
                Box             => Box%Next
            enddo
            Box                 => NewBox
            PreviousBox%Next    => NewBox
        endif
        
        Me%NumberOfBoxes2D      = Me%NumberOfBoxes2D + 1
        Me%NumberOfBoxes3D      = Me%NumberOfBoxes3D + NumberOfLayers
        Box%NumberOfLayers      = NumberOfLayers

        Box%MainID              = Me%NextBoxMainID
        Me%NextBoxMainID        = Me%NextBoxMainID + 1

        allocate(Box%ID(1:Box%NumberOfLayers))

        do i = 1, Box%NumberOfLayers

            Box%ID(i)%Number    = Me%NextBoxID
            write(Number,*)Box%ID(i)%Number
            Number              = trim(adjustl(Number))
            Box%ID(i)%Name      = trim("Box_"//Number)
            Me%NextBoxID        = Me%NextBoxID + 1

        end do


    end subroutine AddBox
        
    !--------------------------------------------------------------------------
    
    subroutine ConstructOutputFluxes2D (FluxesOutputList)

        !Arguments-------------------------------------------------------------
        character(len = *), optional, pointer, dimension(:) :: FluxesOutputList
        
        !External----------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                             :: NumberOfFluxes, i
        type (T_BoxTimeSerie), pointer                      :: NewFluxesTimeSerie

        !Begin-----------------------------------------------------------------
        
        nullify(Me%FirstFluxesTimeSerie)

        !Number of properties which fluxes are to be integrated
        NumberOfFluxes = size(FluxesOutputList)

        call FluxesTimeSerieHeader2D

        do i = 1, NumberOfFluxes

            call AddFluxesTimeSerie(NewFluxesTimeSerie, FluxesOutputList(i))

            call StartTimeSerie(NewFluxesTimeSerie%ObjTimeSerie, Me%ObjTime,        &
                                Me%BoxesFilePath,                                   &
                                PropertyList    = Me%FluxesTimeSerieHeader,         &
                                Extension       = 'BXF',                            &
                                ResultFileName  = NewFluxesTimeSerie%ID%Name,       &
                                STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutputFluxes2D - ModuleBoxDif - ERR02'

        end do


    end subroutine ConstructOutputFluxes2D

    !--------------------------------------------------------------------------

    subroutine ConstructOutputFluxes3D (FluxesOutputList)

        !Arguments-------------------------------------------------------------
        character(len = *), optional, pointer, dimension(:) :: FluxesOutputList
        
        !External----------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                             :: NumberOfFluxes, i
        type (T_BoxTimeSerie), pointer                      :: NewFluxesTimeSerie

        !Begin-----------------------------------------------------------------
        
        nullify(Me%FirstFluxesTimeSerie)


        !Number of properties which fluxes are to be integrated
        NumberOfFluxes = size(FluxesOutputList)

        call FluxesTimeSerieHeader3D

        do i = 1, NumberOfFluxes

            call AddFluxesTimeSerie(NewFluxesTimeSerie, FluxesOutputList(i))

            call StartTimeSerie(NewFluxesTimeSerie%ObjTimeSerie, Me%ObjTime,        &
                                Me%BoxesFilePath,                                   &
                                PropertyList    = Me%FluxesTimeSerieHeader,         &
                                WaterPoints3D   = Me%ExternalVar%WaterPoints3D,     &
                                Extension       = 'BXF',                            &
                                ResultFileName  = NewFluxesTimeSerie%ID%Name,       &
                                STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutputFluxes3D - ModuleBoxDif - ERR02'

        end do


    end subroutine ConstructOutputFluxes3D

    !--------------------------------------------------------------------------
    
    subroutine AddFluxesTimeSerie(FluxesTimeSerie, FluxesName)
        
        !Arguments-------------------------------------------------------------
        type (T_BoxTimeSerie), pointer                      :: FluxesTimeSerie
        character(len = *)                                  :: FluxesName

        !Local-----------------------------------------------------------------
        type (T_BoxTimeSerie), pointer                      :: NewFluxesTimeSerie
        type (T_BoxTimeSerie), pointer                      :: PreviousFluxesTimeSerie
        integer, save                                       :: NextID = 1
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewFluxesTimeSerie)
        nullify  (NewFluxesTimeSerie%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%FirstFluxesTimeSerie)) then
            Me%FirstFluxesTimeSerie       => NewFluxesTimeSerie
            FluxesTimeSerie               => NewFluxesTimeSerie
        else
            PreviousFluxesTimeSerie       => Me%FirstFluxesTimeSerie
            FluxesTimeSerie               => Me%FirstFluxesTimeSerie%Next
            do while (associated(FluxesTimeSerie))
                PreviousFluxesTimeSerie   => FluxesTimeSerie
                FluxesTimeSerie           => FluxesTimeSerie%Next
            enddo
            FluxesTimeSerie               => NewFluxesTimeSerie
            PreviousFluxesTimeSerie%Next  => NewFluxesTimeSerie
        endif

        FluxesTimeSerie%ID%Number = NextID
        NextID                    = NextID + 1

        FluxesTimeSerie%ID%Name = trim(FluxesName)
    
    
    end subroutine AddFluxesTimeSerie

    !--------------------------------------------------------------------------

    subroutine FluxesTimeSerieHeader3D

        !Arguments-------------------------------------------------------------
        
        !External----------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        integer                                             :: i, j, Column
        character(len = 3)                                  :: Char_I, Char_J
        
        !Begin-----------------------------------------------------------------

        nullify (Me%FluxesTimeSerieLine  )
        allocate(Me%FluxesTimeSerieLine  (1:Me%nAdjacentBoxesBoundaries3D))
        Me%FluxesTimeSerieLine = 0.

        nullify (Me%FluxesTimeSerieHeader)
        allocate(Me%FluxesTimeSerieHeader(1:Me%nAdjacentBoxesBoundaries3D))
        Me%FluxesTimeSerieHeader = null_str

        Column = 0

        do i = 0, Me%NumberOfBoxes3D
        do j = 0, Me%NumberOfBoxes3D

            if(Me%AdjacentBoxesBoundaries3D(i,j) .eq. Boundary)then

                Column = Column + 1

                write(Char_I, '(i3)') i
                write(Char_J, '(i3)') j

                Me%FluxesTimeSerieHeader(Column) = 'Flux_'//trim(adjustl(Char_I))//&
                                                       '_'//trim(adjustl(Char_J))

            end if

        enddo
        enddo

    end subroutine FluxesTimeSerieHeader3D

    !--------------------------------------------------------------------------

    subroutine FluxesTimeSerieHeader2D

        !Local-----------------------------------------------------------------
        integer                                             :: i, j, Column
        character(len = 3)                                  :: Char_I, Char_J
        
        !Begin-----------------------------------------------------------------

        nullify (Me%FluxesTimeSerieLine  )
        allocate(Me%FluxesTimeSerieLine  (1:Me%nAdjacentBoxesBoundaries2D))
        Me%FluxesTimeSerieLine = 0.

        nullify (Me%FluxesTimeSerieHeader)
        allocate(Me%FluxesTimeSerieHeader(1:Me%nAdjacentBoxesBoundaries2D))
        Me%FluxesTimeSerieHeader = null_str

        Column = 0

        do i = 0, Me%NumberOfBoxes2D
        do j = 0, Me%NumberOfBoxes2D

            if(Me%AdjacentBoxesBoundaries2D(i,j) .eq. Boundary)then

                Column = Column + 1

                write(Char_I, '(i3)') i
                write(Char_J, '(i3)') j

                Me%FluxesTimeSerieHeader(Column) = 'Flux_'//trim(adjustl(Char_I))//&
                                                       '_'//trim(adjustl(Char_J))

            end if

        enddo
        enddo

    end subroutine FluxesTimeSerieHeader2D

    !--------------------------------------------------------------------------

    subroutine AddScalarTimeSerie(ScalarTimeSerie, ScalarName)
        
        !Arguments-------------------------------------------------------------
        type (T_BoxTimeSerie), pointer                      :: ScalarTimeSerie
        character(len = *)                                  :: ScalarName

        !Local-----------------------------------------------------------------
        type (T_BoxTimeSerie), pointer                      :: NewScalarTimeSerie
        type (T_BoxTimeSerie), pointer                      :: PreviousScalarTimeSerie
        integer, save                                       :: NextID = 1
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewScalarTimeSerie)
        nullify  (NewScalarTimeSerie%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%FirstScalarTimeSerie)) then
            Me%FirstScalarTimeSerie       => NewScalarTimeSerie
            ScalarTimeSerie               => NewScalarTimeSerie
        else
            PreviousScalarTimeSerie       => Me%FirstScalarTimeSerie
            ScalarTimeSerie               => Me%FirstScalarTimeSerie%Next
            do while (associated(ScalarTimeSerie))
                PreviousScalarTimeSerie   => ScalarTimeSerie
                ScalarTimeSerie           => ScalarTimeSerie%Next
            enddo
            ScalarTimeSerie               => NewScalarTimeSerie
            PreviousScalarTimeSerie%Next  => NewScalarTimeSerie
        endif

        ScalarTimeSerie%ID%Number = NextID
        NextID                    = NextID + 1

        ScalarTimeSerie%ID%Name = trim(ScalarName)
    
    
    end subroutine AddScalarTimeSerie
    
    !--------------------------------------------------------------------------
    
    subroutine ConstructOutputScalar2D (ScalarOutputList)

        !Arguments-------------------------------------------------------------
        character(len = *), optional, pointer, dimension(:) :: ScalarOutputList
        
        !External----------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                             :: NumberOfScalares, i
        type (T_BoxTimeSerie), pointer                      :: NewScalarTimeSerie

        !Begin-----------------------------------------------------------------

        
        call ScalarTimeSerieHeader2D

        !Number of properties which values are to be integrated
        NumberOfScalares = size(ScalarOutputList)

        do i = 1, NumberOfScalares

            call AddScalarTimeSerie(NewScalarTimeSerie, ScalarOutputList(i))

            call StartTimeSerie(NewScalarTimeSerie%ObjTimeSerie, Me%ObjTime,        &
                                Me%BoxesFilePath,                                   &
                                PropertyList    = Me%ScalarTimeSerieHeader,         &
                                Extension       = 'BXM',                            &
                                ResultFileName  = NewScalarTimeSerie%ID%Name,       &
                                STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutputScalar2D - ModuleBoxDif - ERR02'
        end do


    end subroutine ConstructOutputScalar2D

    !--------------------------------------------------------------------------
    
    subroutine ConstructOutputScalar3D (ScalarOutputList)

        !Arguments-------------------------------------------------------------
        character(len = *), optional, pointer, dimension(:) :: ScalarOutputList
        
        !External----------------------------------------------------------------
        integer                                             :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                             :: NumberOfScalares, i
        type (T_BoxTimeSerie), pointer                      :: NewScalarTimeSerie

        !Begin-----------------------------------------------------------------

        
        call ScalarTimeSerieHeader3D


        !Number of properties which values are to be integrated
        NumberOfScalares = size(ScalarOutputList)

        do i = 1, NumberOfScalares

            call AddScalarTimeSerie(NewScalarTimeSerie, ScalarOutputList(i))

            call StartTimeSerie(NewScalarTimeSerie%ObjTimeSerie, Me%ObjTime,        &
                                Me%BoxesFilePath,                                   &
                                PropertyList    = Me%ScalarTimeSerieHeader,         &
                                WaterPoints3D   = Me%ExternalVar%WaterPoints3D,     &
                                Extension       = 'BXM',                            &
                                ResultFileName  = NewScalarTimeSerie%ID%Name,       &
                                STAT            = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructOutputScalar3D - ModuleBoxDif - ERR02'
        end do


    end subroutine ConstructOutputScalar3D

    !--------------------------------------------------------------------------

    subroutine ScalarTimeSerieHeader2D

        !Local-----------------------------------------------------------------
        integer                                             :: i, Column
        character(len = 3)                                  :: Char_I
        
        !Begin-----------------------------------------------------------------
        
        nullify (Me%ScalarTimeSerieLine)
        allocate(Me%ScalarTimeSerieLine  (1:Me%NumberOfBoxes2D + 1))
        Me%ScalarTimeSerieLine = 0.

        nullify (Me%ScalarTimeSerieHeader)
        allocate(Me%ScalarTimeSerieHeader(1:Me%NumberOfBoxes2D + 1)) !count with environment box
        Me%ScalarTimeSerieHeader = null_str

        Column = 0

        do i = 0, Me%NumberOfBoxes2D

            Column = Column + 1

            write(Char_I,'(i3)') i

            Me%ScalarTimeSerieHeader(Column) = 'Box_'//trim(adjustl(Char_I))

        enddo


    end subroutine ScalarTimeSerieHeader2D

    !--------------------------------------------------------------------------

    subroutine ScalarTimeSerieHeader3D

        !Local-----------------------------------------------------------------
        integer                                             :: i, Column
        character(len = 3)                                  :: Char_I
        
        !Begin-----------------------------------------------------------------
        
        nullify (Me%ScalarTimeSerieLine)
        allocate(Me%ScalarTimeSerieLine  (1:Me%NumberOfBoxes3D + 1))
        Me%ScalarTimeSerieLine = 0.

        nullify (Me%ScalarTimeSerieHeader)
        allocate(Me%ScalarTimeSerieHeader(1:Me%NumberOfBoxes3D + 1)) !count with environment box
        Me%ScalarTimeSerieHeader = null_str

        Column = 0

        do i = 0, Me%NumberOfBoxes3D

            Column = Column + 1

            write(Char_I,'(i3)') i

            Me%ScalarTimeSerieHeader(Column) = 'Box_'//trim(adjustl(Char_I))

        enddo


    end subroutine ScalarTimeSerieHeader3D

    !--------------------------------------------------------------------------
    
    subroutine SetLimitsPolygon(Polygon)

        !Arguments-------------------------------------------------------------
        type (T_Polygon)                       :: Polygon
        
        !Begin-----------------------------------------------------------------

        Polygon%Limits%Left   = minval(Polygon%VerticesF%X)
        Polygon%Limits%Right  = maxval(Polygon%VerticesF%X)
        Polygon%Limits%Bottom = minval(Polygon%VerticesF%Y)
        Polygon%Limits%Top    = maxval(Polygon%VerticesF%Y)

    end subroutine SetLimitsPolygon


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------
    subroutine GetBoxes2D (BoxDifID, Boxes2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        integer, dimension(:, :),  pointer              :: Boxes2D
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBOXDIF_, Me%InstanceID)

            Boxes2D => Me%Boxes2D

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoxes2D
    
    
    !--------------------------------------------------------------------------
    
    subroutine GetBoxes3D (BoxDifID, Boxes3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        integer, dimension(:, :, :),  pointer           :: Boxes3D
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBOXDIF_, Me%InstanceID)

            Boxes3D => Me%Boxes3D

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetBoxes3D
    
    !--------------------------------------------------------------------------

    
    subroutine GetDTBoxes(Boxes_ID, DT, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: Boxes_ID
        real,  intent(OUT)                  :: DT
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Boxes_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            DT = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetDTBoxes

    !--------------------------------------------------------------------------
    
    subroutine GetNumberOfBoxes (BoxDifID, NumberOfBoxes2D, NumberOfBoxes3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        integer, intent(OUT), optional                  :: NumberOfBoxes2D
        integer, intent(OUT), optional                  :: NumberOfBoxes3D
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if(present(NumberOfBoxes2D))NumberOfBoxes2D = Me%NumberOfBoxes2D
            if(present(NumberOfBoxes3D))NumberOfBoxes3D = Me%NumberOfBoxes3D

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetNumberOfBoxes

    !--------------------------------------------------------------------------

    
    !--------------------------------------------------------------------------
    
    logical function CheckIfInsideBox (BoxDifID, BoxNumber, Xcoord, Ycoord, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: BoxDifID, BoxNumber
        real,    intent(IN)                             :: Xcoord, Ycoord
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        type (T_Box),    pointer                        :: CurrentBox
        type(T_PointF), pointer                         :: Point

        integer                                         :: STAT_, ready_, Layer
        logical                                         :: NotFound

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            allocate(Point)

            Point%X = Xcoord 

            Point%Y = Ycoord

            CheckIfInsideBox = .false. 
            NotFound         = .true. 
        
            CurrentBox       => Me%FirstBox

            do while(associated(CurrentBox) .and. NotFound)

                do Layer = 1, CurrentBox%NumberOfLayers

                    if(BoxNumber == CurrentBox%ID(Layer)%Number)then

                        if(IsPointInsidePolygon(Point, CurrentBox%Polygon))then

                            CheckIfInsideBox = .true. 

                        endif

                        NotFound = .false.

                        exit

                    end if

                enddo

                 CurrentBox => CurrentBox%Next
            end do

            deallocate(Point)

            nullify(CurrentBox, Point)

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end function CheckIfInsideBox

    !--------------------------------------------------------------------------
    
    logical function GetIfBoxInsideDomain (BoxDifID, BoxNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                             :: BoxDifID
        integer, intent(IN)                             :: BoxNumber
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        type (T_Box),    pointer                        :: CurrentBox

        integer                                         :: STAT_, ready_, Layer

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

      
            CurrentBox => Me%FirstBox

            do while(associated(CurrentBox))

                do Layer = 1, CurrentBox%NumberOfLayers

                    if(BoxNumber == CurrentBox%ID(Layer)%Number)then

                        GetIfBoxInsideDomain = CurrentBox%InsideDomain
                        exit

                    end if

                enddo

                CurrentBox => CurrentBox%Next
            end do

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end function GetIfBoxInsideDomain

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    subroutine UnGetBoxDif3D_I(BoxDifID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        integer, dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mBOXDIF_, Me%InstanceID, "UnGetBoxDif3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBoxDif3D_I

    !--------------------------------------------------------------------------

    subroutine UnGetBoxDif2D_I(BoxDifID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        integer, dimension(:,:), pointer                :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mBOXDIF_, Me%InstanceID, "UnGetBoxDif2D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBoxDif2D_I
    
    !--------------------------------------------------------------------------
    
    subroutine UnGetBoxDif3D_R8(BoxDifID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        real(8), dimension(:, :, :), pointer            :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mBOXDIF_, Me%InstanceID, "UnGetBoxDif3D_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBoxDif3D_R8

    !--------------------------------------------------------------------------

    subroutine UnGetBoxDif2D_R8(BoxDifID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BoxDifID
        real(8), dimension(:,:), pointer                :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_UnLock(mBOXDIF_, Me%InstanceID, "UnGetBoxDif2D_R8")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBoxDif2D_R8

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine BoxDifFluxes3D(BoxDifID, FluxX3D, FluxY3D, FluxZ3D, FluxName, &
                              WaterPoints3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID
        real(8), dimension(:,:,:), pointer          :: FluxX3D
        real(8), dimension(:,:,:), pointer          :: FluxY3D
        real(8), dimension(:,:,:), pointer          :: FluxZ3D
        character(len = *)                          :: FluxName
        integer, dimension(:,:,:), pointer          :: WaterPoints3D
        integer, intent(OUT), optional              :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: WILB, WIUB
        integer                                     :: WJLB, WJUB
        integer                                     :: WKLB, WKUB
        integer                                     :: i, j, k
        integer                                     :: IN, OUT

        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_
        
        
        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            WILB = Me%WorkSize3D%ILB
            WIUB = Me%WorkSize3D%IUB
            WJLB = Me%WorkSize3D%JLB
            WJUB = Me%WorkSize3D%JUB
            WKLB = Me%WorkSize3D%KLB
            WKUB = Me%WorkSize3D%KUB

            Me%ExternalVar%WaterPoints3D => WaterPoints3D


            Me%ExternalVar%FluxX3D              => FluxX3D
            Me%ExternalVar%FluxY3D              => FluxY3D
            Me%ExternalVar%FluxZ3D              => FluxZ3D
            Me%ExternalVar%CurrentTimeSerieName =  FluxName


            Me%Fluxes3D(:,:) = 0.

            do k = WKLB,  WKUB
            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then
                    
                    !XX axis
                    if (Me%BoundaryFace3DX(i,j,k)) then
                   
                       OUT = Me%Boxes3D(i,j,k  )
                       IN  = Me%Boxes3D(i,j+1,k)

                       Me%Fluxes3D(OUT,IN ) = Me%Fluxes3D(OUT,IN) + FluxX3D(i,j+1,k)
                       Me%Fluxes3D(IN, OUT) = -Me%Fluxes3D(OUT,IN)

                    end if
                end if
            enddo
            enddo
            enddo

            do k = WKLB,  WKUB
            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then
                    !YY axis
                    if (Me%BoundaryFace3DY(i,j,k))then
                    
                        OUT = Me%Boxes3D(i,j,k  )
                        IN  = Me%Boxes3D(i+1,j,k)

                        Me%Fluxes3D(OUT,IN ) =  Me%Fluxes3D(OUT,IN) + FluxY3D(i+1,j,k)
                        Me%Fluxes3D(IN, OUT) = -Me%Fluxes3D(OUT,IN)

                    end if
                end if
            enddo
            enddo
            enddo

            do k = WKLB,  WKUB
            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then
                    !ZZ axis
                    if (Me%BoundaryFace3DZ(i,j,k))then
                    
                        OUT = Me%Boxes3D(i,j,k  )
                        IN  = Me%Boxes3D(i,j,k+1)

                        Me%Fluxes3D(OUT,IN)  =  Me%Fluxes3D(OUT,IN) + FluxZ3D(i,j,k+1)
                        Me%Fluxes3D(IN, OUT) = -Me%Fluxes3D(OUT,IN)

                    end if
                end if
            enddo
            enddo
            enddo
            
            call OutputFluxesTimeSerie3D

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine BoxDifFluxes3D

    !--------------------------------------------------------------------------

    subroutine BoxDifFluxes2D(BoxDifID, FluxX2D, FluxY2D, FluxName, WaterPoints2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID
        real(8), dimension(:,:), pointer            :: FluxX2D
        real(8), dimension(:,:), pointer            :: FluxY2D
        character(len = *)                          :: FluxName
        integer, dimension(:,:), pointer            :: WaterPoints2D
        integer, intent(OUT), optional              :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: WILB, WIUB
        integer                                     :: WJLB, WJUB
        integer                                     :: i, j
        integer                                     :: IN, OUT

        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            WILB = Me%WorkSize2D%ILB
            WIUB = Me%WorkSize2D%IUB
            WJLB = Me%WorkSize2D%JLB
            WJUB = Me%WorkSize2D%JUB

            Me%ExternalVar%WaterPoints2D => WaterPoints2D

            Me%ExternalVar%FluxX2D              => FluxX2D
            Me%ExternalVar%FluxY2D              => FluxY2D
            Me%ExternalVar%CurrentTimeSerieName =  FluxName


            Me%Fluxes2D(:,:) = 0.

            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%ExternalVar%WaterPoints2D(i,j) .eq. WaterPoint)then

                    !XX axis
                    if (Me%BoundaryFace2DX(i,j)) then
                    
                        OUT = Me%Boxes2D(i,j  )
                        IN  = Me%Boxes2D(i,j+1)

                        Me%Fluxes2D(OUT,IN ) =  Me%Fluxes2D(OUT,IN) + FluxX2D(i,j+1)
                        Me%Fluxes2D(IN, OUT) = -Me%Fluxes2D(OUT,IN)

                    end if

                    !YY axis
                    if (Me%BoundaryFace2DY(i,j)) then
                    
                        OUT = Me%Boxes2D(i,j  )
                        IN  = Me%Boxes2D(i+1,j)

                        Me%Fluxes2D(OUT,IN ) =  Me%Fluxes2D(OUT,IN) + FluxY2D(i+1,j)
                        Me%Fluxes2D(IN, OUT) = -Me%Fluxes2D(OUT,IN)

                    end if
                end if

            enddo
            enddo

            call OutputFluxesTimeSerie2D

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine BoxDifFluxes2D

    !--------------------------------------------------------------------------
    
    subroutine BoxDifScalar3D_R4(BoxDifID, Scalar3D_R4, ScalarName,WaterPoints3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID
        real(4), dimension(:,:,:), pointer          :: Scalar3D_R4
        character(len = *)                          :: ScalarName
        integer, dimension(:,:,:), pointer          :: WaterPoints3D
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: WILB, WIUB
        integer                                     :: WJLB, WJUB
        integer                                     :: WKLB, WKUB
        integer                                     :: i, j, k
        integer                                     :: BoxNumber

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            WILB = Me%WorkSize3D%ILB
            WIUB = Me%WorkSize3D%IUB
            WJLB = Me%WorkSize3D%JLB
            WJUB = Me%WorkSize3D%JUB
            WKLB = Me%WorkSize3D%KLB
            WKUB = Me%WorkSize3D%KUB

            Me%ExternalVar%WaterPoints3D => WaterPoints3D

            Me%ExternalVar%Scalar3D_R4          => Scalar3D_R4
            Me%ExternalVar%CurrentTimeSerieName =  ScalarName

            Me%ScalarTimeSerieLine(:) = 0.

            do k = WKLB,  WKUB
            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%Boxes3D(i,j,k) > -55 .and.                        &
                   Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then
                    
                    BoxNumber                             = Me%Boxes3D(i,j,k)
                    Me%ScalarTimeSerieLine(BoxNumber + 1) = Me%ScalarTimeSerieLine(BoxNumber + 1) + &
                                                            Scalar3D_R4(i,j,k)

                end if

            enddo
            enddo
            enddo

            call OutputScalarTimeSerie


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine BoxDifScalar3D_R4

    !--------------------------------------------------------------------------
    
    subroutine BoxDifScalar3D_R8(BoxDifID, Scalar3D_R8, ScalarName, WaterPoints3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID
        real(8), dimension(:,:,:), pointer          :: Scalar3D_R8
        character(len = *)                          :: ScalarName
        integer, dimension(:,:,:), pointer          :: WaterPoints3D
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: WILB, WIUB
        integer                                     :: WJLB, WJUB
        integer                                     :: WKLB, WKUB
        integer                                     :: i, j, k
        integer                                     :: BoxNumber

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        

        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            WILB = Me%WorkSize3D%ILB
            WIUB = Me%WorkSize3D%IUB
            WJLB = Me%WorkSize3D%JLB
            WJUB = Me%WorkSize3D%JUB
            WKLB = Me%WorkSize3D%KLB
            WKUB = Me%WorkSize3D%KUB

            Me%ExternalVar%WaterPoints3D => WaterPoints3D


            Me%ExternalVar%Scalar3D_R8          => Scalar3D_R8
            Me%ExternalVar%CurrentTimeSerieName =  ScalarName

            Me%ScalarTimeSerieLine(:) = 0.

            do k = WKLB,  WKUB
            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%Boxes3D(i,j,k) > -55 .and.                        &
                   Me%ExternalVar%WaterPoints3D(i,j,k) .eq. WaterPoint)then
                    
                    BoxNumber                             = Me%Boxes3D(i,j,k)
                    Me%ScalarTimeSerieLine(BoxNumber + 1) = Me%ScalarTimeSerieLine(BoxNumber + 1) + &
                                                            Scalar3D_R8(i,j,k)

                end if

            enddo
            enddo
            enddo

            call OutputScalarTimeSerie


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine BoxDifScalar3D_R8

    !--------------------------------------------------------------------------

    subroutine BoxDifScalar2D_R8(BoxDifID, Scalar2D_R8, ScalarName, WaterPoints2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID
        real(8), dimension(:,:), pointer            :: Scalar2D_R8
        character(len = *)                          :: ScalarName
        integer, dimension(:,:), pointer            :: WaterPoints2D
        integer, intent(OUT), optional              :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: WILB, WIUB
        integer                                     :: WJLB, WJUB
        integer                                     :: i, j
        integer                                     :: BoxNumber

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            WILB = Me%WorkSize2D%ILB
            WIUB = Me%WorkSize2D%IUB
            WJLB = Me%WorkSize2D%JLB
            WJUB = Me%WorkSize2D%JUB

            Me%ExternalVar%WaterPoints2D => WaterPoints2D


            Me%ExternalVar%Scalar2D_R8          => Scalar2D_R8
            Me%ExternalVar%CurrentTimeSerieName =  ScalarName

            Me%ScalarTimeSerieLine(:) = 0.

            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%Boxes2D(i,j) > -55 .and.                        &
                   Me%ExternalVar%WaterPoints2D(i,j) .eq. WaterPoint)then

                    BoxNumber                             = Me%Boxes2D(i,j)
                    Me%ScalarTimeSerieLine(BoxNumber + 1) = Me%ScalarTimeSerieLine(BoxNumber + 1) + &
                                                            Scalar2D_R8(i,j)

                end if

            enddo
            enddo
            
            call OutputScalarTimeSerie

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine BoxDifScalar2D_R8


    !--------------------------------------------------------------------------

    subroutine BoxDifScalar2D_R4(BoxDifID, Scalar2D_R4, ScalarName, WaterPoints2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID
        real(4), dimension(:,:), pointer            :: Scalar2D_R4
        character(len = *)                          :: ScalarName
        integer, dimension(:,:), pointer            :: WaterPoints2D
        integer, intent(OUT), optional              :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: WILB, WIUB
        integer                                     :: WJLB, WJUB
        integer                                     :: i, j
        integer                                     :: BoxNumber

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_
        
        call Ready(BoxDifID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            WILB = Me%WorkSize2D%ILB
            WIUB = Me%WorkSize2D%IUB
            WJLB = Me%WorkSize2D%JLB
            WJUB = Me%WorkSize2D%JUB

            Me%ExternalVar%WaterPoints2D => WaterPoints2D


            Me%ExternalVar%Scalar2D_R4          => Scalar2D_R4
            Me%ExternalVar%CurrentTimeSerieName =  ScalarName

            Me%ScalarTimeSerieLine(:) = 0.

            do j = WJLB , WJUB
            do i = WILB,  WIUB

                if(Me%Boxes2D(i,j) > -55 .and.                        &
                   Me%ExternalVar%WaterPoints2D(i,j) .eq. WaterPoint)then

                    BoxNumber                             = Me%Boxes2D(i,j)
                    Me%ScalarTimeSerieLine(BoxNumber + 1) = Me%ScalarTimeSerieLine(BoxNumber + 1) + &
                                                            Scalar2D_R4(i,j)

                end if

            enddo
            enddo
            
            call OutputScalarTimeSerie

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine BoxDifScalar2D_R4
    
    !--------------------------------------------------------------------------

    subroutine OutputFluxesTimeSerie2D
        
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type(T_BoxTimeSerie), pointer                   :: FluxTimeSerie
        character(LEN = StringLength)                   :: CurrentTimeSerieName
        integer                                         :: Column
        integer                                         :: IN, OUT, i, j
        !Begin-----------------------------------------------------------------


        CurrentTimeSerieName = Me%ExternalVar%CurrentTimeSerieName

        Column = 0

        do i  = 0, Me%NumberOfBoxes2D
        do j  = 0, Me%NumberOfBoxes2D

            if(Me%AdjacentBoxesBoundaries2D(i,j) .eq. Boundary)then

                Column = Column + 1
                OUT    = i
                IN     = j

                Me%FluxesTimeSerieLine (Column) = Me%Fluxes2D(OUT, IN)

            end if

        enddo
        enddo

        FluxTimeSerie => Me%FirstFluxesTimeSerie

        do while(associated(FluxTimeSerie))

            if(trim(FluxTimeSerie%ID%Name) .eq. trim(CurrentTimeSerieName))then

                call WriteTimeSerieLine(FluxTimeSerie%ObjTimeSerie,             &
                                        DataLine = Me%FluxesTimeSerieLine,      &
                                        STAT      = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)stop 'OutputFluxesTimeSerie - ModuleBoxDif -ERR01'

                Me%FluxesTimeSerieLine = null_real

            end if

            FluxTimeSerie => FluxTimeSerie%Next
        end do




    end subroutine OutputFluxesTimeSerie2D

    !--------------------------------------------------------------------------

    subroutine OutputFluxesTimeSerie3D
        
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL
        
        !Local-----------------------------------------------------------------
        type(T_BoxTimeSerie), pointer                   :: FluxTimeSerie
        character(LEN = StringLength)                   :: CurrentTimeSerieName
        integer                                         :: Column
        integer                                         :: IN, OUT, i, j
        !Begin-----------------------------------------------------------------


        CurrentTimeSerieName = Me%ExternalVar%CurrentTimeSerieName

        Column = 0

        do i  = 0, Me%NumberOfBoxes3D
        do j  = 0, Me%NumberOfBoxes3D

            if(Me%AdjacentBoxesBoundaries3D(i,j) .eq. Boundary)then

                Column = Column + 1
                OUT    = i
                IN     = j

                Me%FluxesTimeSerieLine (Column) = Me%Fluxes3D(OUT, IN)

            end if

        enddo
        enddo

        FluxTimeSerie => Me%FirstFluxesTimeSerie

        do while(associated(FluxTimeSerie))

            if(trim(FluxTimeSerie%ID%Name) .eq. trim(CurrentTimeSerieName))then

                call WriteTimeSerieLine(FluxTimeSerie%ObjTimeSerie,             &
                                        DataLine = Me%FluxesTimeSerieLine,      &
                                        STAT      = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)stop 'OutputFluxesTimeSerie - ModuleBoxDif -ERR01'

                Me%FluxesTimeSerieLine = null_real

            end if

            FluxTimeSerie => FluxTimeSerie%Next
        end do

    end subroutine OutputFluxesTimeSerie3D

    !--------------------------------------------------------------------------
   
    subroutine OutputScalarTimeSerie
        
        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL

        !Local-----------------------------------------------------------------
        type(T_BoxTimeSerie), pointer                   :: ScalarTimeSerie
        character(LEN = StringLength)                   :: CurrentTimeSerieName

        !Begin-----------------------------------------------------------------

        CurrentTimeSerieName = Me%ExternalVar%CurrentTimeSerieName

        ScalarTimeSerie => Me%FirstScalarTimeSerie

        do while(associated(ScalarTimeSerie))

            if(trim(ScalarTimeSerie%ID%Name) .eq. trim(CurrentTimeSerieName))then
                
                call WriteTimeSerieLine(ScalarTimeSerie%ObjTimeSerie,           &
                                        DataLine    = Me%ScalarTimeSerieLine,   &
                                        STAT        = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_)stop 'OutputFluxesTimeSerie - ModuleBoxDif -ERR01'

                Me%ScalarTimeSerieLine = null_real

            end if

            ScalarTimeSerie => ScalarTimeSerie%Next
        end do

    end subroutine OutputScalarTimeSerie

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillBoxDif(BoxDifID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: BoxDifID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_, STAT_CALL          

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers  
        type(T_BoxTimeSerie), pointer       :: FluxTimeSerie
        type(T_BoxTimeSerie), pointer       :: ScalarTimeSerie
        type(T_Box         ), pointer       :: CurrentBox

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BoxDifID, ready_)    


cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBOXDIF_,  Me%InstanceID)

            if (nUsers == 0) then

                nUsers = DeassociateInstance(mTIME_,            Me%ObjTime)
                if (nUsers == 0) stop 'KillBoxDif - ModuleBoxDif -ERR01'

                nUsers = DeassociateInstance(mHORIZONTALGRID_,  Me%ObjHorizontalGrid)
                if (nUsers == 0) stop 'KillBoxDif - ModuleBoxDif -ERR03'

                !Kill fluxes time series
                FluxTimeSerie => Me%FirstFluxesTimeSerie

                do while(associated(FluxTimeSerie))

                    call KillTimeSerie(FluxTimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)stop 'KillBoxDif - ModuleBoxDif -ERR04'

                    FluxTimeSerie => FluxTimeSerie%Next
                end do

                nullify(Me%FirstFluxesTimeSerie)

                !Kill Scalar time series
                ScalarTimeSerie => Me%FirstScalarTimeSerie

                do while(associated(ScalarTimeSerie))

                    call KillTimeSerie(ScalarTimeSerie%ObjTimeSerie, STAT = STAT_CALL)
                    if(STAT_CALL .ne. SUCCESS_)stop 'KillBoxDif - ModuleBoxDif -ERR05'

                    ScalarTimeSerie => ScalarTimeSerie%Next
                end do

                nullify(Me%FirstScalarTimeSerie)

                !Kill boxes
                CurrentBox => Me%FirstBox

                do while(associated(CurrentBox))

                    deallocate(CurrentBox%Polygon%Vertices )
                    deallocate(CurrentBox%Polygon%VerticesF)

                    CurrentBox => CurrentBox%Next
                end do

                nullify(Me%FirstBox)

                if(associated(Me%Boxes2D))                   deallocate  (Me%Boxes2D)
                nullify      (Me%Boxes2D)
                
                if(associated(Me%Boxes3D))                   deallocate  (Me%Boxes3D)
                nullify      (Me%Boxes3D)
                
                if(associated(Me%Fluxes2D))                  deallocate  (Me%Fluxes2D)
                nullify      (Me%Fluxes2D)

                if(associated(Me%Fluxes3D))                  deallocate  (Me%Fluxes3D)
                nullify      (Me%Fluxes3D)


                if(associated(Me%FluxesTimeSerieLine))       deallocate  (Me%FluxesTimeSerieLine)
                nullify      (Me%FluxesTimeSerieLine)

                if(associated(Me%ScalarTimeSerieLine))       deallocate  (Me%ScalarTimeSerieLine)
                nullify      (Me%ScalarTimeSerieLine)

                if(associated(Me%FluxesTimeSerieHeader))     deallocate  (Me%FluxesTimeSerieHeader)
                nullify      (Me%FluxesTimeSerieHeader)
                
                if(associated(Me%ScalarTimeSerieHeader))     deallocate  (Me%ScalarTimeSerieHeader)
                nullify      (Me%ScalarTimeSerieHeader)
                
                if(associated(Me%AdjacentBoxesBoundaries2D)) deallocate  (Me%AdjacentBoxesBoundaries2D)
                nullify      (Me%AdjacentBoxesBoundaries2D)
                
                if(associated(Me%AdjacentBoxesBoundaries3D)) deallocate  (Me%AdjacentBoxesBoundaries3D)
                nullify      (Me%AdjacentBoxesBoundaries3D)
                
                if(associated(Me%BoundaryFace3DX))           deallocate  (Me%BoundaryFace3DX)
                nullify      (Me%BoundaryFace3DX)

                if(associated(Me%BoundaryFace3DY))           deallocate  (Me%BoundaryFace3DY)
                nullify      (Me%BoundaryFace3DY)

                if(associated(Me%BoundaryFace3DZ))           deallocate  (Me%BoundaryFace3DZ)
                nullify      (Me%BoundaryFace3DZ)

                if(associated(Me%BoundaryFace2DX))           deallocate  (Me%BoundaryFace2DX)
                nullify      (Me%BoundaryFace2DX)

                if(associated(Me%BoundaryFace2DY))           deallocate  (Me%BoundaryFace2DY)
                nullify      (Me%BoundaryFace2DY)

                !Deallocates Instance
                call DeallocateInstance
                BoxDifID = 0

                STAT_ = SUCCESS_

            end if
        
        else 
            STAT_ = ready_
        end if cd1


        if (present(STAT)) STAT = STAT_

        !------------------------------------------------------------------------

    end subroutine KillBoxDif
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_BoxDif), pointer          :: AuxObjBoxDif
        type (T_BoxDif), pointer          :: PreviousObjBoxDif

        !Updates pointers
        if (Me%InstanceID == FirstObjBoxDif%InstanceID) then
            FirstObjBoxDif      => FirstObjBoxDif%Next
        else
            PreviousObjBoxDif   => FirstObjBoxDif
            AuxObjBoxDif        => FirstObjBoxDif%Next
            do while (AuxObjBoxDif%InstanceID /= Me%InstanceID)
                PreviousObjBoxDif => AuxObjBoxDif
                AuxObjBoxDif      => AuxObjBoxDif%Next
            enddo

            !Now update linked list
            PreviousObjBoxDif%Next => AuxObjBoxDif%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjBoxDif_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBoxDif_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBoxDif_ID > 0) then
            call LocateObjBoxDif (ObjBoxDif_ID)
            ready_ = VerifyReadLock (mBOXDIF_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjBoxDif (BoxDifID)

        !Arguments-------------------------------------------------------------
        integer                                     :: BoxDifID

        !Local-----------------------------------------------------------------

        Me => FirstObjBoxDif
        do while (associated (Me))
            if (Me%InstanceID == BoxDifID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me))                          &
            stop 'ModuleBoxDif - LocateObjBoxDif - ERR01'

    end subroutine LocateObjBoxDif

    !--------------------------------------------------------------------------

end module ModuleBoxDif

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------










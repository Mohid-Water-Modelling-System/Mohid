!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : MovingProfile
! PROGRAM       : MainMovingProfile
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : January 2005
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to extract moving time series from MOHID HDF5 files
!
!------------------------------------------------------------------------------

!DataFile: needs a nomfich.dat to redirect to input data file
!   IN_MODEL                : char                  [-]         !Path to input data file
!   ROOT_SRT                : char                  [-]         !Path to folder where output files are written
!
!DataFile: Keywords in the input data file
!   MOVING_TIMESERIE        : 0/1                   [1]         !Write a time serie given by a track line 
!                                                               !(e.g probe results on a moving boat)
!   MOVING_TIMESERIE_FILE   : char                  [-]         !File with the track line
!   X_COLUMN                : int                   [2]         !Column with coordinates in XX axis (meters or degrees)
!   Y_COLUMN                : int                   [3]         !Column with coordinates in YY axis (meters or degrees)
!   USE_LATLON              : 0/1                   [1]         !Use geographic coordinates or metric coordinates
!                                                               !to convert XYZ to IJK. Must be in agreement with
!                                                               !the coordinate type given by the track line.
!   CORRECT_TRAJECTORY      : 0/1                   [0]         !Correct trajectory due to GPS errors, by interpolating 
!                                                               !points not trustable
!   CORRECTION_METHOD       : int                   [1]         !Method to correct trajectory. 1 - Deviation 2-Forecast trajectory
!   TOLERANCE               : real              [5e-4 or 50]    !Tolerance for deviation from trajectory in moving time serie
!                                                               !Default value varies on type of coordinate used
!   TOLERANCE_INCREASE      : real                 [0.1]        !Factor to increase tolerance if next point 
!   MAXIMUM_TOLERANCE       : real             [10xTolerance]   !Maximum tolerance to be used
!   MAXIMUM_FORECAST        : int                  [10]         !Maximum number of forecasts. If exceeded, resumes next point
!                                                               !is not inside tolerance
!   <BeginHDF5File>                                                                 
!   NAME                    : char                  [-]         !Name of HDF5 file 
!   <EndHDF5File>                                               !(specify one block for each HDF5 file to use)
!
!   <BeginParameter>                                            
!   HDF_GROUP               : char                  [-]         !Path of the HDF5 group in HDF5 file for the property
!   PROPERTY                : char                  [-]         !Property name (should be equal to the one specified 
!   <EndParameter>                                              !in ModuleGlobalData)
!                                                               !(specify one block for each property)                                                                                 
!   <BeginTimeSerie>
!   NAME                    : char                  [-]         !Name for output time serie file
!   LOCALIZATION_I          : int                   [-]         !Redundant information
!   LOCALIZATION_J          : int                   [-]         !Redundant information
!   LOCALIZATION_K          : int                   [-]         !Redundant information
!   <EndTimeSerie>                                              !(specify one block for each output time serie file)

program MohidMovingProfile

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleHorizontalGrid
    use ModuleProfile
    use ModuleHDF5
    use ModuleDrawing
    use ModuleTimeSerie

    implicit none

    !CPU Time variables
    type (T_Time)                                           :: InitialSystemTime, FinalSystemTime
    real                                                    :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                                   :: F95Time
    
    !Parameters
    integer, parameter                                      :: Deviation            = 1
    integer, parameter                                      :: ForecastTrajectory   = 2

    type       T_Parameter
        character(len=StringLength)                         :: Name
        character(len=PathLength)                           :: Group
        integer                                             :: Rank
        character(len=StringLength)                         :: Units
        type(T_Parameter), pointer                          :: Next
    end type  T_Parameter

    type     T_Grid
        real,    dimension(:,:  ), pointer                  :: XX, YY
        real,    dimension(:,:  ), pointer                  :: Bathymetry
        real,    dimension(:,:,:), pointer                  :: VerticalZ
        integer, dimension(:,:  ), pointer                  :: WaterPoints2D
        integer, dimension(:,:,:), pointer                  :: WaterPoints3D
    end type T_Grid

    type       T_HDF5File
        integer                                             :: ObjHDF5              = 0
        character(len=StringLength)                         :: Name
        integer                                             :: FirstSignificantInstant
        integer                                             :: LastSignificantInstant
        integer                                             :: nSignificantInstants = 0
        integer                                             :: CurrentInstant       = 1
        type(T_Time), dimension(:), pointer                 :: InstantsArray 
        logical                                             :: IsSignificant        = .true.

        type(T_HDF5File), pointer                           :: Next
    end type  T_HDF5File


    type T_MovingProfile
        character(PathLength)                               :: MovingProfileFile
        character(PathLength)                               :: DataFile
        integer                                             :: ObjTimeSerie         = 0
        integer                                             :: ObjProfile           = 0
        integer                                             :: X_Column
        integer                                             :: Y_Column
        logical                                             :: Use_LatLon           = .true.
        logical                                             :: CorrectTrajectory    = .false.
        integer                                             :: Correction_Method
        real                                                :: Tolerance
        real                                                :: MaxTolerance
        real                                                :: ToleranceIncrease
        real                                                :: MaxNumberOfForecasts
        integer,        dimension(:    ), pointer           :: I, J
        real,           dimension(:    ), pointer           :: X, Y
        type(T_Time),   dimension(:    ), pointer           :: Times
        real,           dimension(:    ), pointer           :: Values
        type(T_Size3D)                                      :: Size3D
        type(T_Time)                                        :: InitialData
        character(len=StringLength)                         :: TimeUnits
        integer                                             :: nPositions
        integer                                             :: nParameters
        integer                                             :: nColumns
        real,           dimension(:,:  ), pointer           :: Values2D
        real,           dimension(:,:,:), pointer           :: Values3D
        real,           dimension(:,:,:), pointer           :: PreviousField3D, PreviousSZZ
        real,           dimension(:,:,:), pointer           :: NextField3D,     NextSZZ
        real,           dimension(:,:,:), pointer           :: SZZ
        real                                                :: DT
        type(T_Time)                                        :: BeginTime
        type(T_Time)                                        :: EndTime
        type(T_Time)                                        :: CurrentTime
        type(T_Grid)                                        :: Grid
        type(T_HDF5File),  pointer                          :: FirstHDF5File
        type(T_Parameter), pointer                          :: FirstParameter
        integer                                             :: ObjTime           = 0
        integer                                             :: ObjHorizontalGrid = 0
        integer                                             :: ObjEnterData      = 0
        character(len=StringLength), dimension(:,:), pointer:: PropertyList
    end type T_MovingProfile

    type(T_MovingProfile), pointer                        :: Me


    call ConstructMohidMovingProfile
    call ModifyMohidMovingProfile
    call KillMohidMovingProfile

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidMovingProfile

        !Local-----------------------------------------------------------------
        integer             :: STAT_CALL

        !Begin-----------------------------------------------------------------

        allocate(Me)
        
        call StartUpMohid("MohidMovingProfile")

        call StartCPUTime

        call ReadGlobalData

        call OrganizeHDF5Files

        call ComputeLocations

        call AllocateFields

        call StartComputeTime(Me%ObjTime, Me%BeginTime, Me%EndTime, DT = Me%DT, &
                              VariableDT = .false., STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructMohidMovingProfile MohidMovingProfile - ERR01'

        call ConstructPropertyList

        call StartProfile(ProfileID         = Me%ObjProfile,                    &
                          ObjTime           = Me%ObjTime,                       &
                          ProfileDataFile   = trim(Me%DataFile),                &
                          WaterPoints2D     = Me%Grid%WaterPoints2D,            &
                          nProperties       = Me%nParameters,                   &
                          KUB               = Me%Size3D%KUB,                    &
                          PropertyList      = Me%PropertyList,                  &
                          ClientName        = "MovingProfile",                  &
                          VerifyLocation    = .false.,                          &
                          OutTime           = Me%Times,                         &
                          STAT              = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                             &
            stop 'ConstructMohidMovingProfile - MohidMovingProfile - ERR02'


    end subroutine ConstructMohidMovingProfile
    
    !--------------------------------------------------------------------------


    subroutine AllocateFields


        allocate(Me%Values3D       (Me%Size3D%ILB:Me%Size3D%IUB, &
                                    Me%Size3D%JLB:Me%Size3D%JUB, &
                                    Me%Size3D%KLB:Me%Size3D%KUB))
        
        allocate(Me%PreviousField3D(Me%Size3D%ILB:Me%Size3D%IUB, &
                                    Me%Size3D%JLB:Me%Size3D%JUB, &
                                    Me%Size3D%KLB:Me%Size3D%KUB))

        allocate(Me%NextField3D    (Me%Size3D%ILB:Me%Size3D%IUB, &
                                    Me%Size3D%JLB:Me%Size3D%JUB, &
                                    Me%Size3D%KLB:Me%Size3D%KUB))


        allocate(Me%PreviousSZZ    (Me%Size3D%ILB  :Me%Size3D%IUB, &
                                    Me%Size3D%JLB  :Me%Size3D%JUB, &
                                    Me%Size3D%KLB-1:Me%Size3D%KUB))

        allocate(Me%NextSZZ        (Me%Size3D%ILB  :Me%Size3D%IUB, &
                                    Me%Size3D%JLB  :Me%Size3D%JUB, &
                                    Me%Size3D%KLB-1:Me%Size3D%KUB))


        allocate(Me%SZZ            (Me%Size3D%ILB  :Me%Size3D%IUB, &
                                    Me%Size3D%JLB  :Me%Size3D%JUB, &
                                    Me%Size3D%KLB-1:Me%Size3D%KUB))


    end subroutine AllocateFields
    
    !--------------------------------------------------------------------------

    subroutine ConstructPropertyList

        !Local-----------------------------------------------------------------
        type(T_Parameter), pointer                      :: ObjParameter
        integer                                         :: n

        !Begin-----------------------------------------------------------------

        n = 1

        allocate(Me%PropertyList(1:Me%nParameters, 1:2))

        ObjParameter => Me%FirstParameter

        do while(associated(ObjParameter))

            Me%PropertyList(n, 1) = ObjParameter%Name
            Me%PropertyList(n, 2) = ObjParameter%Units

            n = n + 1

            ObjParameter => ObjParameter%Next

        end do

    end subroutine ConstructPropertyList
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidMovingProfile
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        integer                                     :: STAT_CALL
        integer                                     :: iPosition
        type(T_HDF5File),  pointer                  :: HDF5File
        type(T_Time)                                :: PreviousTime
        type(T_Time)                                :: NextTime
        type(T_Parameter), pointer                  :: ObjParameter
        integer                                     :: OutI, OutJ
        integer                                     :: UpdateCounter
        integer                                     :: FivePercent
        real                                        :: CPUTime, LastCPUTime
        !Begin-----------------------------------------------------------------
        
        Running         = .true.
        Me%CurrentTime  = Me%BeginTime
        
        write(*,*)
        write(*,*)'Getting results based on file:'
        write(*,*)
        write(*,*)'------->   ', trim(Me%MovingProfileFile)
        write(*,*)
        write(*,*)'Total number of positions :', Me%nPositions
        write(*,*)

        FivePercent     = int(Me%nPositions/20)
        UpdateCounter   = 0

        write(*,*)'Started conversion...'
        write(*,*)

        do iPosition = 1, Me%nPositions

            OutI = Me%I(iPosition)
            OutJ = Me%J(iPosition)


            HDF5File => Me%FirstHDF5File

            do while(associated(HDF5File))

                PreviousTime = HDF5File%InstantsArray(HDF5File%CurrentInstant)

                do while(PreviousTime .lt. Me%CurrentTime)

                    HDF5File%CurrentInstant = HDF5File%CurrentInstant + 1
                    
                    PreviousTime = HDF5File%InstantsArray(HDF5File%CurrentInstant)

                enddo

                PreviousTime = HDF5File%InstantsArray(HDF5File%CurrentInstant-1)
                NextTime     = HDF5File%InstantsArray(HDF5File%CurrentInstant  )
                
                ObjParameter => Me%FirstParameter

                do while(associated(ObjParameter))

                    call ReadSZZ(HDF5File, HDF5File%CurrentInstant, PreviousTime, NextTime)

                    call ReadPreviousAndNextField(HDF5File, ObjParameter, HDF5File%CurrentInstant)

                    call InterpolateMatrix3DInTime(Me%CurrentTime,                                  &
                                                   Me%Size3D, PreviousTime, Me%PreviousField3D,     &
                                                   NextTime, Me%NextField3D,                        &
                                                   Me%Values3D)

                    call WriteProfile(Me%ObjProfile,                                                &
                                      Me%Values3D,                                                  &
                                      SZZ    = Me%SZZ,                                              &
                                      LocI   = OutI,                                                &
                                      LocJ   = OutJ,                                                &
                                      STAT   = STAT_CALL)


                    ObjParameter => ObjParameter%Next

                end do


                HDF5File => HDF5File%Next

            enddo

            if(UpdateCounter .ge. FivePercent)then

                UpdateCounter = 0
                write(*,*)'Just exported position ', iPosition, ' of ', Me%nPositions

            else

                UpdateCounter = UpdateCounter + 1

            end if


            if(iPosition .lt. Me%nPositions)then
                Me%DT = Me%Times(iPosition+1) - Me%Times(iPosition)
            end if
            
            Me%CurrentTime = Me%CurrentTime + Me%DT

            call ActualizeCurrentTime(Me%ObjTime, Me%DT, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMohidMovingProfile - ERR01'

            if (abs(Me%CurrentTime - Me%EndTime) > Me%DT / 10.) then
                Running = .true.
            else
                Running = .false.
            endif


            call CPU_TIME(CPUTime)
            if (CPUTime - LastCPUTime > 60.) then
                LastCPUTime = CPUTime
                call PrintProgress(Me%ObjTime, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidLand - MohidLand - ERR05'
            endif

        enddo
    
    
    
    end subroutine ModifyMohidMovingProfile
    
    
    subroutine ReadGlobalData

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        real                                        :: DefaultTolerance
        
        !Begin-----------------------------------------------------------------

        call ReadFileName('IN_MODEL', Me%DataFile, "MohidMovingProfile", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidMovingProfile - ERR01'

        call ConstructEnterData (Me%ObjEnterData, Me%DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidMovingProfile - ERR02'

        !Moving Time Serie File
        call GetData(Me%MovingProfileFile, Me%ObjEnterData, iflag,          &
                     keyword      = 'MOVING_PROFILE_FILE',                  &
                     SearchType   = FromFile,                               &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR03'
            
        call GetData(Me%X_Column, Me%ObjEnterData, iflag,                   &
                     keyword      = 'X_COLUMN',                             &
                     SearchType   = FromFile,                               &
                     Default      = 2,                                      &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR04'


        call GetData(Me%Y_Column, Me%ObjEnterData, iflag,                   &
                     keyword      = 'Y_COLUMN',                             &
                     SearchType   = FromFile,                               &
                     Default      = 3,                                      &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0)                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR05'

        call GetData(Me%Use_LatLon, Me%ObjEnterData, iflag,                 &
                     keyword      = 'USE_LATLON',                           &
                     SearchType   = FromFile,                               &
                     Default      = .true.,                                 &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR07'

        call GetData(Me%CorrectTrajectory, Me%ObjEnterData, iflag,          &
                     keyword      = 'CORRECT_TRAJECTORY',                   &
                     SearchType   = FromFile,                               &
                     Default      = .false.,                                &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR09'

        call GetData(Me%Correction_Method, Me%ObjEnterData, iflag,          &
                     keyword      = 'CORRECTION_METHOD',                    &
                     SearchType   = FromFile,                               &
                     Default      = Deviation,                              &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR10'

            

        if(Me%Use_LatLon)then

            DefaultTolerance = 5e-4

        else
            
            DefaultTolerance = 50.

        end if

        call GetData(Me%Tolerance, Me%ObjEnterData, iflag,                  &
                     keyword      = 'TOLERANCE',                            &
                     SearchType   = FromFile,                               &
                     Default      = DefaultTolerance,                       &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR11'


        call GetData(Me%ToleranceIncrease, Me%ObjEnterData, iflag,          &
                     keyword      = 'TOLERANCE_INCREASE',                   &
                     SearchType   = FromFile,                               &
                     Default      = 0.1,                                    &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR12'

        call GetData(Me%MaxTolerance, Me%ObjEnterData, iflag,               &
                     keyword      = 'MAXIMUM_TOLERANCE',                    &
                     SearchType   = FromFile,                               &
                     Default      = DefaultTolerance * 10,                  &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR13'


        call GetData(Me%MaxNumberOfForecasts, Me%ObjEnterData, iflag,       &
                     keyword      = 'MAXIMUM_FORECAST',                     &
                     SearchType   = FromFile,                               &
                     Default      = 0.1,                                    &
                     ClientModule = 'MohidMovingProfile',                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                          &
            stop 'ReadGlobalData - MohidMovingProfile - ERR14'

        call ReadMovingProfileFile

        call ReadParameters

        call ConstructHDF5FilesList

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidMovingProfile - ERR03'

    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------

    subroutine ReadMovingProfileFile
        
        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL
        integer                                 :: ObjEnterData = 0
        integer                                 :: n
        real                                    :: MinDT, CurrentDT
        real, dimension(:  ), pointer           :: CumulativeTime
        real, dimension(:,:), pointer           :: FieldData

        !Begin-----------------------------------------------------------------

        call ConstructEnterData (ObjEnterData, Me%MovingProfileFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR01'
            
        call StartTimeSerieInput(Me%ObjTimeSerie,                               &
                                 Me%MovingProfileFile,                          &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR02'

        call GetTimeSerieDataMatrix(Me%ObjTimeSerie, FieldData, STAT = STAT_CALL)
                                    
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR03'

        call GetTimeSerieDataValues (Me%ObjTimeSerie, Me%nPositions, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR03'
        
        call GetTimeSerieDataColumns (Me%ObjTimeSerie, Me%nColumns, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR03'

        call GetTimeSerieInitialData (Me%ObjTimeSerie, Me%InitialData, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR03'
        
        call GetTimeSerieTimeUnits (Me%ObjTimeSerie, Me%TimeUnits, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR03'

        allocate(Me%X           (1:Me%nPositions))
        allocate(Me%Y           (1:Me%nPositions))
        allocate(CumulativeTime (1:Me%nPositions))

        Me%X            (1:Me%nPositions) = FieldData(:, Me%X_Column)
        Me%Y            (1:Me%nPositions) = FieldData(:, Me%Y_Column)
        CumulativeTime  (1:Me%nPositions) = FieldData(:, 1)
        
        allocate(Me%Times (1:Me%nPositions))
        allocate(Me%I     (1:Me%nPositions)) 
        allocate(Me%J     (1:Me%nPositions))

        do n = 1, Me%nPositions

            Me%Times(n) = Me%InitialData + CumulativeTime(n) 

        end do

        call WritePositions(CumulativeTime)

        if(Me%CorrectTrajectory) call CorrectPositions(CumulativeTime, FieldData)

        Me%BeginTime = Me%Times(1)
        Me%EndTime   = Me%Times(Me%nPositions)

        if(Me%nPositions == 1)then
            MinDT = 10.
        else
            MinDT = 1e32
        end if

        do n = 2, Me%nPositions

            CurrentDT = CumulativeTime(n) - CumulativeTime(n-1)

            if(CurrentDT<0)then
                write(*,*)'Inconsistency found in time serie input file'
                write(*,*)'Time1 greater than Time2'
                stop 'ReadMovingProfileFile - MohidMovingProfile - ERR04'
            end if

            if(CurrentDT < MinDT)then
                MinDT = CurrentDT
            endif

        end do

        deallocate(CumulativeTime)

        Me%DT = MinDT/2.

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidMovingProfile - ERR05'


    end subroutine ReadMovingProfileFile
    
    !--------------------------------------------------------------------------
    
    
    subroutine CorrectPositions(CumulativeTime, FieldData)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:  ), pointer   :: CumulativeTime
        real, dimension(:,:), pointer   :: FieldData

        !Local-----------------------------------------------------------------
        character(len=StringLength)     :: XYZFileName, LinesFileName, CorrectedFileName
        integer                         :: XYZUnit, LinesUnit, CorrectedFileUnit, STAT_CALL
        integer                         :: NameLength, n
        logical                         :: NeedsCorrection
        real,    dimension(:), pointer  :: NewX, NewY
        logical, dimension(:), pointer  :: Valid
        type(T_PointF), pointer         :: NextPoint, PreviousPoint, ForecastPoint
        real                            :: DT, PreviousDT
        real                            :: DistanceX, DistanceY, CurrentTolerance
        integer                         :: nForecasts
        character(len=2)                :: CharFormat               = '  '

        !Begin-----------------------------------------------------------------
        
        allocate(Valid(1:Me%nPositions))
        allocate(NewX(1:Me%nPositions))
        allocate(NewY(1:Me%nPositions))
        allocate(NextPoint, PreviousPoint, ForecastPoint)
       
        NewX(:) = Me%X(:)
        NewY(:) = Me%Y(:)


        PreviousDT          = 1e32 
        NeedsCorrection     = .false.
       
        select case(Me%Correction_Method)

            case(Deviation)

                do n = 1, Me%nPositions - 1
            
                    PreviousPoint%X = NewX(n)
                    PreviousPoint%Y = NewY(n)
                    NextPoint%X     = NewX(n+1)
                    NextPoint%Y     = NewY(n+1)

                    DT = CumulativeTime(n+1) - CumulativeTime(n)

                    !if time difference between 2 points is higher that last DT 
                    !then it is considerec that there is a gap in the data. 
                    !Correction is resumed for the next point
                    if(DT .le. PreviousDT)then

                        PreviousDT = DT

                        if(IsPointInsideCircle(NextPoint, PreviousPoint, Me%Tolerance))then
                            Valid(n+1) = .true.
                        else
                            Valid(n+1) = .false.
                        end if
                    
                    else
                        
                        Valid(n+1) = .true.

                    end if

                end do


            case(ForecastTrajectory)

                Valid(:) = .true.

                nForecasts = 0

                CurrentTolerance = Me%Tolerance

                do n = 1, Me%nPositions - 1
            
                    PreviousPoint%X = NewX(n)
                    PreviousPoint%Y = NewY(n)
                    NextPoint%X     = NewX(n+1)
                    NextPoint%Y     = NewY(n+1)

                    DT = CumulativeTime(n+1) - CumulativeTime(n)

                    !if time difference between 2 points is higher that last DT 
                    !then it is considerec that there is a gap in the data. 
                    !Correction is resumed for the next point
                    if(DT .le. PreviousDT)then

                        PreviousDT = DT

                        if(IsPointInsideCircle(NextPoint, PreviousPoint, CurrentTolerance))then

                            DistanceX        = (NextPoint%X - PreviousPoint%X)
                            DistanceY        = (NextPoint%Y - PreviousPoint%Y)

                            ForecastPoint%X  = NextPoint%X + DistanceX
                            ForecastPoint%Y  = NextPoint%Y + DistanceY

                            CurrentTolerance = Me%Tolerance

                        else

                            if(nForecasts .gt. Me%MaxNumberOfForecasts)then

                                nForecasts = 0
                        
                                write(*,*)'MaxNumberOfForecasts EXCEEDED! n = ', n

                                NeedsCorrection = .false.

                            else

                                NeedsCorrection = .true.


                                nForecasts = nForecasts + 1

                                NewX(n+1)       = ForecastPoint%X
                                NewY(n+1)       = ForecastPoint%Y

                                ForecastPoint%X = ForecastPoint%X + DistanceX
                                ForecastPoint%Y = ForecastPoint%Y + DistanceY

                                CurrentTolerance = CurrentTolerance + CurrentTolerance * Me%ToleranceIncrease

                                if(CurrentTolerance .gt. Me%MaxTolerance)then

                                    CurrentTolerance = Me%Tolerance

                                    write(*,*)'Exceeded maximum tolerance! n = ', n

                                end if

                            endif



                        end if


                    end if

                end do

        end select


        NeedsCorrection = .true.

        deallocate(NextPoint, PreviousPoint, ForecastPoint)

        if(NeedsCorrection)then

            NameLength          = len_trim(Me%MovingProfileFile)
            XYZFileName         = Me%MovingProfileFile(:NameLength-4)//'corrected.xyz'
            LinesFileName       = Me%MovingProfileFile(:NameLength-4)//'corrected.lin'
            CorrectedFileName   = Me%MovingProfileFile(:NameLength-4)//'corrected.srm'

            call UnitsManager(XYZUnit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR01'

            open(XYZUnit, File = trim(XYZFileName), Form='FORMATTED', &
                 status = 'UNKNOWN',IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR02'

            call UnitsManager(LinesUnit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR03'

            open(LinesUnit, File = trim(LinesFileName), Form='FORMATTED', &
                 status = 'UNKNOWN',IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR04'

            call UnitsManager(CorrectedFileUnit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR03'

            open(CorrectedFileUnit, File = trim(CorrectedFileName), Form='FORMATTED', &
                 status = 'UNKNOWN',IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR04'



            write(XYZUnit,   *)'<begin_xyz>'
            write(LinesUnit, *)'<begin_line>'
            call WriteDataLine(CorrectedFileUnit,  'SERIE_INITIAL_DATA', Me%InitialData)
            call WriteDataLine(CorrectedFileUnit,  'TIME_UNITS', Me%TimeUnits)
            write(CorrectedFileUnit, *)
            write(CorrectedFileUnit, *)
            write(CorrectedFileUnit, *)'<BeginTimeSerie>'

            write(CharFormat, '(i2)')Me%nColumns-1

            do n = 1, Me%nPositions

                if(Valid(n))then
            
                    write(XYZUnit,  *)NewX(n), NewY(n), CumulativeTime(n)
                    write(LinesUnit,*)NewX(n), NewY(n)
                    write(CorrectedFileUnit, '(i10,'//CharFormat//'f14.6)')int(FieldData(n, 1)), FieldData(n, 2:)

                end if

            end do

            write(CorrectedFileUnit, *)'<EndTimeSerie>'
            write(XYZUnit,   *)'<end_xyz>'
            write(LinesUnit, *)'<end_line>'

            call UnitsManager(XYZUnit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR05'
        
            call UnitsManager(LinesUnit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR06'
            
            
            call UnitsManager(CorrectedFileUnit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'CorrectPositions - MohidMovingProfile - ERR05'
 
            write(*,*)
            write(*,*)'GPS positions were corrected!'
            write(*,*)'Please compare corrected XYZ and line files with the original ones'
            write(*,*)
            stop 'CorrectPositions - MohidMovingProfile - ERR06'
        end if

        deallocate(NewX, NewY, Valid)

    end subroutine CorrectPositions
    
    !--------------------------------------------------------------------------

    subroutine WritePositions(CumulativeTime)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:  ), pointer   :: CumulativeTime

        !Local-----------------------------------------------------------------
        character(len=StringLength)     :: XYZFileName, LinesFileName
        integer                         :: XYZUnit, LinesUnit, STAT_CALL
        integer                         :: NameLength, n

        !Begin-----------------------------------------------------------------
        
        NameLength    = len_trim(Me%MovingProfileFile)
        XYZFileName   = Me%MovingProfileFile(:NameLength-4)//'.xyz'
        LinesFileName = Me%MovingProfileFile(:NameLength-4)//'.lin'

        call UnitsManager(XYZUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WritePositions - MohidMovingProfile - ERR01'

        open(XYZUnit, File = trim(XYZFileName), Form='FORMATTED', &
             status = 'UNKNOWN',IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WritePositions - MohidMovingProfile - ERR02'

        call UnitsManager(LinesUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WritePositions - MohidMovingProfile - ERR03'

        open(LinesUnit, File = trim(LinesFileName), Form='FORMATTED', &
             status = 'UNKNOWN',IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WritePositions - MohidMovingProfile - ERR04'


        write(XYZUnit,   *)'<begin_xyz>'
        write(LinesUnit, *)'<begin_line>'

        do n = 1, Me%nPositions
            
            write(XYZUnit,  *)Me%X(n), Me%Y(n), CumulativeTime(n)
            write(LinesUnit,*)Me%X(n), Me%Y(n)

        end do

        write(XYZUnit,   *)'<end_xyz>'
        write(LinesUnit, *)'<end_line>'

        call UnitsManager(XYZUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WritePositions - MohidMovingProfile - ERR05'
        
        call UnitsManager(LinesUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WritePositions - MohidMovingProfile - ERR06'


    end subroutine WritePositions

    !--------------------------------------------------------------------------

    subroutine ConstructHDF5FilesList

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, ClientNumber
        type (T_HDF5File),       pointer        :: NewHDF5File
        logical                                 :: BlockFound
        logical                                 :: AtLeastOneBlock = .false.
        integer                                 :: iflag

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginHDF5File>',    &
                                        block_end       = '<EndHDF5File>',      &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddHDF5File       (NewHDF5File)

                    
                    call GetData(NewHDF5File%Name,                                  &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'NAME',                             &
                                 ClientModule = 'MohidMovingProfile',                &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'ConstructHDF5File - MohidMovingProfile - ERR00'

                    nullify(NewHDF5File)

                else cd2
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                                      

                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'ConstructHDF5FilesList - MohidMovingProfile - ERR01'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructHDF5FilesList - MohidMovingProfile - ERR02'
            else cd1
                stop 'ConstructHDF5FilesList - MohidMovingProfile - ERR03'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No HDF5 file block is indicated in input file. '
            stop 'ConstructHDF5FilesList - MohidMovingProfile - ERR04'
        end if

    end subroutine ConstructHDF5FilesList

    !--------------------------------------------------------------------------

    subroutine AddHDF5File(ObjHDF5File)

        !Arguments-------------------------------------------------------------
        type (T_HDF5File),     pointer           :: ObjHDF5File

        !Local-----------------------------------------------------------------
        type (T_HDF5File),     pointer           :: PreviousHDF5File
        type (T_HDF5File),     pointer           :: NewHDF5File

        !Begin-----------------------------------------------------------------

        !Allocates new HDF5File
        allocate (NewHDF5File)
        nullify  (NewHDF5File%Next)

        !Insert new Parameter into list and makes current ?? point to it
        if (.not. associated(Me%FirstHDF5File)) then
            Me%FirstHDF5File         => NewHDF5File
            ObjHDF5File              => NewHDF5File
        else
            PreviousHDF5File         => Me%FirstHDF5File
            ObjHDF5File              => Me%FirstHDF5File%Next
            do while (associated(ObjHDF5File))
                PreviousHDF5File     => ObjHDF5File
                ObjHDF5File          => ObjHDF5File%Next
            enddo
            ObjHDF5File              => NewHDF5File
            PreviousHDF5File%Next    => NewHDF5File
        end if

    end subroutine AddHDF5File

    !--------------------------------------------------------------------------

    subroutine ConstructGridFromHDF5(ObjHDF5File)
    
        !Arguments-------------------------------------------------------------
        type (T_HDF5File),     pointer              :: ObjHDF5File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer, dimension(7)                       :: Dimensions
        integer                                     :: Rank, nItems, n
        character(len=StringLength)                 :: GroupName
        type(T_Size3D)                              :: Size3D
        integer                                     :: i, j, k

        !Begin-----------------------------------------------------------------


        call ReadParameterRanks(ObjHDF5File)


        call GetHDF5GroupNumberOfItems(ObjHDF5File%ObjHDF5, "/Grid", nItems, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR01'
            
        do n = 1, nItems

            call GetHDF5GroupID(ObjHDF5File%ObjHDF5, FatherGroupName = "/Grid",     &
                                GroupPosition = n, GroupName =  GroupName,          &
                                Rank = Rank, Dimensions = Dimensions,               &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR02'

            select case(Rank)

                case(2)

                    Size3D%ILB = 1
                    Size3D%IUB = Dimensions(1)
                    Size3D%JLB = 1
                    Size3D%JUB = Dimensions(2)

                    call HDF5SetLimits (ObjHDF5File%ObjHDF5, Size3D%ILB, Size3D%IUB,    &
                                        Size3D%JLB,Size3D%JUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR03'

                case(3)

                    Size3D%ILB = 1
                    Size3D%IUB = Dimensions(1)
                    Size3D%JLB = 1
                    Size3D%JUB = Dimensions(2)
                    Size3D%KLB = 1
                    Size3D%KUB = Dimensions(3)

                    call HDF5SetLimits (ObjHDF5File%ObjHDF5, Size3D%ILB,                &
                                        Size3D%IUB, Size3D%JLB, Size3D%JUB,             &   
                                        Size3D%KLB, Size3D%KUB, STAT = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR04'


            end select


            select case(GroupName)

                case ('Bathymetry')

                    Me%Size3D%ILB = Size3D%ILB
                    Me%Size3D%IUB = Size3D%IUB
                    Me%Size3D%JLB = Size3D%JLB
                    Me%Size3D%JUB = Size3D%JUB

                    allocate(Me%Grid%Bathymetry(Size3D%ILB:Size3D%IUB, Size3D%JLB:Size3D%JUB))

                    call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid",                 &
                                      "Bathymetry",                                 &
                                      Array2D   = Me%Grid%Bathymetry,               &
                                      STAT      = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR05'


                case('ConnectionX', 'Longitude')

                    allocate(Me%Grid%XX(Size3D%ILB:Size3D%IUB, Size3D%JLB:Size3D%JUB))

                    if(Me%Use_LatLon)then


                        call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid",                     &
                                          "Longitude",                                      &
                                          Array2D   = Me%Grid%XX,                           &
                                          STAT      = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR06'

                    else


                        call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid",                     &
                                          "ConnectionX",                                    &
                                          Array2D   = Me%Grid%XX,                           &
                                          STAT      = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR07'

                    end if

                case('ConnectionY', 'Latitude')

                    allocate(Me%Grid%YY(Size3D%ILB:Size3D%IUB,                              &
                                        Size3D%JLB:Size3D%JUB))

                    if(Me%Use_LatLon)then

                        call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid",                     &
                                          "Latitude",                                       &
                                          Array2D   = Me%Grid%YY,                           &
                                          STAT      = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                        &
                            stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR08'

                    else

                        call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid",                     &
                                          "ConnectionY",                                    &
                                          Array2D   = Me%Grid%YY,                           &
                                          STAT      = STAT_CALL)
                        if (STAT_CALL .NE. SUCCESS_)                                        &
                            stop 'ConstructGridFromHDF5 - MohidMovingProfile - ERR09'

                    end if


                case('WaterPoints3D')

                    Me%Size3D%KLB = Size3D%KLB
                    Me%Size3D%KUB = Size3D%KUB

                    allocate(Me%Grid%WaterPoints3D(Size3D%ILB:Size3D%IUB,                   &
                                                   Size3D%JLB:Size3D%JUB,                   &
                                                   Size3D%KLB:Size3D%KUB))

                    allocate(Me%Grid%WaterPoints2D(Size3D%ILB:Size3D%IUB,                   &
                                                   Size3D%JLB:Size3D%JUB))
                    
                    do j = Size3D%JLB,Size3D%JUB
                    do i = Size3D%ILB,Size3D%IUB
                    do k = Size3D%JLB,Size3D%KUB

                        Me%Grid%WaterPoints2D(i,j) =  Me%Grid%WaterPoints3D(i,j,k)

                    enddo
                    enddo
                    enddo



                case default

            end select

        end do

        
    end subroutine ConstructGridFromHDF5
    
    !--------------------------------------------------------------------------
    
    subroutine ReadParameterRanks(ObjHDF5File)
        
        !Arguments-------------------------------------------------------------
        type(T_HDF5File), pointer                   :: ObjHDF5File
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_Parameter), pointer                  :: ObjParameter
        character(len=StringLength)                 :: ParameterName

        !Begin-----------------------------------------------------------------
        
        ObjParameter => Me%FirstParameter

        do while(associated(ObjParameter))

            ParameterName = null_str
        
            call GetHDF5GroupID(ObjHDF5File%ObjHDF5, ObjParameter%Group,  1,    &
                                ParameterName,                                  &
                                ObjParameter%Units, ObjParameter%Rank,          &
                                STAT = STAT_CALL)                                
            if (STAT_CALL .NE. SUCCESS_) stop 'ReadParameterRanks - MohidMovingProfile - ERR01'

            ObjParameter => ObjParameter%Next

        enddo 
        
        nullify(ObjParameter)

    end subroutine ReadParameterRanks
    
    
    
    !--------------------------------------------------------------------------

    subroutine ComputeLocations

        !Local-----------------------------------------------------------------
        integer                         :: n
        !Begin-----------------------------------------------------------------


        do n = 1, Me%nPositions

            call GetIJLocation(Me%I(n), Me%J(n), n)

        end do


    end subroutine ComputeLocations

    !--------------------------------------------------------------------------

    subroutine  GetIJLocation(ICell, JCell, CurrentLocation)
       
        !Arguments-------------------------------------------------------------
        integer, intent(out)                        :: ICell, JCell
        integer                                     :: CurrentLocation

        !Local-----------------------------------------------------------------
        real,   dimension(:)  , pointer             :: XX, YY

        !Begin-----------------------------------------------------------------

        ICell = 0
        JCell = 0

        XX => Me%Grid%XX(1,:)

        YY => Me%Grid%YY(:,1)

        call LocateCell (XX, YY,                    &
                         Me%X(CurrentLocation),     &
                         Me%Y(CurrentLocation),     &
                         Me%Size3D%ILB,             &
                         Me%Size3D%IUB+1,           &
                         Me%Size3D%JLB,             &
                         Me%Size3D%JUB+1,           &
                         ICell, JCell)

    end subroutine  GetIJLocation

    !--------------------------------------------------------------------------

    subroutine ReadParameters

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, ClientNumber
        type (T_Parameter),       pointer       :: NewParameter
        logical                                 :: BlockFound
        logical                                 :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

        Me%nParameters = 0

        ! Obtain Parameters for the Time Serie
do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<BeginParameter>',   &
                                        block_end       = '<EndParameter>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then

                    AtLeastOneBlock = .true.
                    
                    call AddParameter                    (NewParameter)

                    call ConstructTimeSerieParameters    (NewParameter)

                    nullify(NewParameter)

                else cd2
                    call Block_Unlock(Me%ObjEnterData,                          & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                    stop 'ReadParameters - MohidMovingProfile - ERR01'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadParameters - MohidMovingProfile - ERR02'
            else cd1
                stop 'ReadParameters - MohidMovingProfile - ERR03'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No property block is indicated in input file. '
            stop 'ReadParameters - MohidMovingProfile - ERR04'
        end if

    end subroutine ReadParameters


    subroutine AddParameter (ObjParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),     pointer           :: ObjParameter

        !Local-----------------------------------------------------------------
        type (T_Parameter),     pointer           :: PreviousParameter
        type (T_Parameter),     pointer           :: NewParameter

        !Begin-----------------------------------------------------------------

        !Allocates new Parameter
        allocate (NewParameter)
        nullify  (NewParameter%Next)

        if (.not. associated(Me%FirstParameter)) then
            Me%FirstParameter         => NewParameter
            ObjParameter              => NewParameter
        else
            PreviousParameter         => Me%FirstParameter
            ObjParameter              => Me%FirstParameter%Next
            do while (associated(ObjParameter))
                PreviousParameter     => ObjParameter
                ObjParameter          => ObjParameter%Next
            enddo
            ObjParameter              => NewParameter
            PreviousParameter%Next    => NewParameter
        end if

        Me%nParameters = Me%nParameters + 1

    end subroutine AddParameter

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSerieParameters (NewParameter)

        !Arguments-------------------------------------------------------------
        type (T_Parameter),      pointer          :: NewParameter

        !Local-----------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------
        
        ! Obtain parameter name
        call GetData(NewParameter%Name,                         &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'PROPERTY',                 &
                     ClientModule = 'MohidMovingProfile',       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                            &
            stop 'ConstructTimeSerieParameters - MohidMovingProfile - ERR01'

        if (.not.CheckPropertyName(NewParameter%Name)) then
            write(*,*)
            write(*,*) 'The property name is not recognised by the model.'
            write(*,*) 'ConstructTimeSerieParameters - MohidMovingProfile - WARN02' 
        end if

        ! Obtain parameter group
        call GetData(NewParameter%Group,                        &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'HDF_GROUP',                &
                     ClientModule = 'MohidMovingProfile',       &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                            &
            stop 'ConstructTimeSerieParameters - MohidMovingProfile - ERR02'

    end subroutine ConstructTimeSerieParameters 

    !--------------------------------------------------------------------------

    subroutine OrganizeHDF5Files
       
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i
        type (T_HDF5File), pointer                  :: HDF5File
        logical                                     :: exist
        integer                                     :: HDF5_READ
        type(T_Time), dimension(:), pointer         :: AuxInstantsArray 
        logical                                     :: FoundFirstInstant    = .false.
        logical                                     :: FoundLastInstant     = .false.
        type(T_Time)                                :: FirstTime
        type(T_Time)                                :: LastTime
        integer                                     :: NumberOfInstants     = 0
        logical                                     :: HasReadParameters    = .false.

        !Begin-----------------------------------------------------------------
        
        HasReadParameters    = .false.


        HDF5File => Me%FirstHDF5File
        
        do while (associated(HDF5File))

            NumberOfInstants = 0

            inquire(FILE = HDF5File%Name, EXIST = exist)
            if (.not. exist) then
                write(*,*)'HDF5 file does not exist:'//trim(HDF5File%Name)
                stop 'OrganizeHDF5Files - MohidMovingProfile - ERR01'
            endif

            call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

            call ConstructHDF5 (HDF5File%ObjHDF5, trim(HDF5File%Name),      &
                                HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                    &
                stop 'OrganizeHDF5Files - MohidMovingProfile - ERR02'

            !Obtain start and end times of HDF5 file
            !(obtain number of instants) 
            call GetHDF5GroupNumberOfItems(HDF5File%ObjHDF5, "/Time",       &
                                           NumberOfInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                    & 
                stop 'OrganizeHDF5Files - MohidMovingProfile - ERR03'


            !(obtain HDF5 start time)
            FirstTime = HDF5TimeInstant(1, HDF5File)

            !(obtain HDF5 end time)
            LastTime  = HDF5TimeInstant(NumberOfInstants, HDF5File)

            if(FirstTime > Me%EndTime .or. LastTime < Me%BeginTime)then

                HDF5File%IsSignificant = .false.

            else

                HDF5File%IsSignificant = .true.

                allocate(AuxInstantsArray(1:NumberOfInstants))

                HDF5File%nSignificantInstants = 0
                FoundFirstInstant    = .false.
                FoundLastInstant     = .false.

                do i = 1, NumberOfInstants

                    AuxInstantsArray(i) = HDF5TimeInstant(i, HDF5File)

                    if(.not. FoundFirstInstant .and. AuxInstantsArray(i) .ge. Me%BeginTime)then

                        FoundFirstInstant                = .true.

                        HDF5File%FirstSignificantInstant = i - 1   !get the instant before begin time

                    end if

                    if(.not. FoundLastInstant  .and. AuxInstantsArray(i) .gt. Me%EndTime)then

                        FoundLastInstant                = .true.

                        HDF5File%LastSignificantInstant = i        !get the instant after begin time

                    end if

                end do

                HDF5File%nSignificantInstants = HDF5File%LastSignificantInstant - &
                                                HDF5File%FirstSignificantInstant + 1

                allocate(HDF5File%InstantsArray(1:HDF5File%nSignificantInstants))

                do i = 1, HDF5File%nSignificantInstants

                    HDF5File%InstantsArray(i) = AuxInstantsArray(i + HDF5File%FirstSignificantInstant - 1) 

                enddo


                deallocate(AuxInstantsArray)


                if(HDF5File%nSignificantInstants > 0 .and. .not. HasReadParameters) then

                    HasReadParameters = .true.

                    call ConstructGridFromHDF5(HDF5File)

                end if

            end if

            HDF5File => HDF5File%Next

        end do

    end subroutine OrganizeHDF5Files

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant, ObjHDF5File)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        type(T_HDF5File), pointer               :: ObjHDF5File
        
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL

        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5File%ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjHDF5File%ObjHDF5,      &
                             GroupName      = "/Time",                  &
                             Name           = "Time",                   &
                             Array1D        = TimeVector,               &
                             OutputNumber   = Instant,                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                      &
            stop 'HDF5TimeInstant - MohidMovingProfile - ERR01'

        call SetDate(HDF5TimeInstant, Year  = TimeVector(1),            &
                     Month  = TimeVector(2), Day      = TimeVector(3),  &
                     Hour   = TimeVector(4), Minute   = TimeVector(5),  &
                     Second = TimeVector(6))

        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------

    subroutine ReadPreviousAndNextField(ObjHDF5File, ObjParameter, CurrentInstant)
        
        !Arguments-------------------------------------------------------------
        type(T_HDF5File ), pointer              :: ObjHDF5File
        type(T_Parameter), pointer              :: ObjParameter
        integer, intent(in)                     :: CurrentInstant

        !Local-----------------------------------------------------------------
        integer                                 :: PreviousInstant, NextInstant
        integer, dimension(3)                   :: lower_bound, upper_bound

        !Begin-----------------------------------------------------------------
        
        PreviousInstant = CurrentInstant  - 1
        
        NextInstant     = CurrentInstant  + ObjHDF5File%FirstSignificantInstant - 1
        PreviousInstant = PreviousInstant + ObjHDF5File%FirstSignificantInstant - 1

        lower_bound(1)  = Me%Size3D%ILB
        lower_bound(2)  = Me%Size3D%JLB
        lower_bound(3)  = Me%Size3D%KLB

        upper_bound(1)  = Me%Size3D%IUB
        upper_bound(2)  = Me%Size3D%JUB
        upper_bound(3)  = Me%Size3D%KUB

        
        select case(ObjParameter%Rank)

            case(2)

                stop 'Cannot output profiles from 2D variables'

            case(3)

                call HDF5ReadHyperSlab (ObjHDF5File%ObjHDF5, trim(ObjParameter%Group),  &
                                        trim(ObjParameter%Name),                        &
                                        lower_bound, upper_bound,                       &
                                        Array3D      = Me%PreviousField3D,              &
                                        OutputNumber = PreviousInstant)


                call HDF5ReadHyperSlab (ObjHDF5File%ObjHDF5, trim(ObjParameter%Group),  &
                                        trim(ObjParameter%Name),                        &
                                        lower_bound, upper_bound,                       &
                                        Array3D      = Me%NextField3D,                  &
                                        OutputNumber = NextInstant)

        end select


    end subroutine ReadPreviousAndNextField
    
    !--------------------------------------------------------------------------

    subroutine ReadSZZ(ObjHDF5File, CurrentInstant, PreviousTime, NextTime)
        
        !Arguments-------------------------------------------------------------
        type(T_HDF5File ), pointer              :: ObjHDF5File
        integer, intent(in)                     :: CurrentInstant
        type(T_Time)                            :: PreviousTime, NextTime

        !Local-----------------------------------------------------------------
        integer                                 :: PreviousInstant, NextInstant
        integer, dimension(3)                   :: lower_bound, upper_bound
        type(T_Size3D)                          :: Size3D
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        PreviousInstant = CurrentInstant  - 1
        
        NextInstant     = CurrentInstant  + ObjHDF5File%FirstSignificantInstant - 1
        PreviousInstant = PreviousInstant + ObjHDF5File%FirstSignificantInstant - 1


        Size3D          = Me%Size3D
        Size3D%KLB      = Size3D%KLB - 1

        lower_bound(1)  = Size3D%ILB
        lower_bound(2)  = Size3D%JLB
        lower_bound(3)  = Size3D%KLB 

        upper_bound(1)  = Size3D%IUB
        upper_bound(2)  = Size3D%JUB
        upper_bound(3)  = Size3D%KUB 


        call HDF5SetLimits   (ObjHDF5File%ObjHDF5, Size3D%ILB, Size3D%IUB, Size3D%JLB,      &
                              Size3D%JUB, Size3D%KLB, Size3D%KUB, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR02'
        
        call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid/VerticalZ",               &
                          "Vertical", Array3D = Me%PreviousSZZ,                 &
                          OutputNumber = PreviousInstant, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR03'
        
        
        call HDF5ReadData(ObjHDF5File%ObjHDF5, "/Grid/VerticalZ",               &
                          "Vertical", Array3D = Me%NextSZZ,                     &
                          OutputNumber = NextInstant, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ReadHDF5 - ERR03'


        call InterpolateMatrix3DInTime(Me%CurrentTime,                          &
                                       Size3D, PreviousTime,                    &
                                       Me%PreviousSZZ,                          &
                                       NextTime, Me%NextSZZ,                    &
                                       Me%SZZ)

    end subroutine ReadSZZ

    !--------------------------------------------------------------------------

    subroutine KillMohidMovingProfile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_HDF5File), pointer                   :: HDF5File
        !Begin-----------------------------------------------------------------


        call DeallocateFields

        HDF5File => Me%FirstHDF5File
        
        do while (associated(HDF5File))

            call KillHDF5(HDF5File%ObjHDF5, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'KillMohidMovingProfile - ERR01'

            HDF5File => HDF5File%Next

        end do


        call KillProfile(Me%ObjProfile, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'KillMohidMovingProfile - ERR02'

        call StopCPUTime

        call ShutdownMohid ("MohidMovingProfile", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidMovingProfile
    
    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        
        call cpu_time(TotalCPUTime)

        ElapsedSeconds = FinalSystemTime - InitialSystemTime

    end subroutine StopCPUTime
    
    !--------------------------------------------------------------------------

    subroutine DeallocateFields

        deallocate(Me%Values3D       )
        deallocate(Me%PreviousField3D)
        deallocate(Me%NextField3D    )
        deallocate(Me%SZZ            )
        deallocate(Me%PreviousSZZ    )
        deallocate(Me%NextSZZ        )

    end subroutine DeallocateFields
   
    !--------------------------------------------------------------------------


end program MohidMovingProfile

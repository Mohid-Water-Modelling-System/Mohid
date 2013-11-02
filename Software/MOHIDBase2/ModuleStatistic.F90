!-------------------------------------------------------------------------
!        IST/MARETEC, Marine Modelling Group, Mohid2000 modelling system
!-------------------------------------------------------------------------
!BOI
! !TITLE: Mohid2000 Statistics module 
! !AUTHORS: Paulo Chambel, Frank Braunschweig
! !AFFILIATION: IST/MARETEC, Marine Modelling Group
! !DATE: May2002
! !INTRODUCTION: Module which makes basic statistic operations over matrixes
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

Module ModuleStatistic

    use ModuleGlobalData
    use ModuleTime                  
    use ModuleStopWatch, only: StartWatch, StopWatch
    use ModuleEnterData
    use ModuleFunctions
    use ModuleHDF5             

    implicit none 

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  ::  ConstructStatistic
    private ::      AllocateInstance
    private ::      ReadDataFile
    private ::      AllocateStatisticMatrixes
    
    !Modifier
    public  ::  ModifyStatistic
    private ::      StatisticAnalysis3D
    private ::          ModifyGlobalStatistic
    private ::          ModifyDailyStatistic
    private ::          ModifyMonthlyStatistic
    private ::          ModifySpecificHourStatistic
    private ::          ModifyClassification
    private ::      StatisticAnalysis2D
    private ::          ModifyGlobalStatistic2D
    private ::          ModifyDailyStatistic2D
    private ::          ModifyMonthlyStatistic2D
    private ::          ModifySpecificHourStatistic2D
    private ::          ModifyClassification2D
    private ::      WriteValuesToFileHDF5

    public  ::  AddStatisticLayers
    private ::      AverageValueBetweenDepths
    private ::      AverageValueBetweenLayers

                
    !Selector
    public  ::  GetStatisticMethod
    public  ::  GetStatisticLayerDef
    public  ::  GetStatisticParameters
    public  ::  GetStatisticLayersNumber
    public  ::  GetStatisticClasses
    public  ::  GetStatisticClassesNumber
    public  ::  UnGetStatistic

    !Destructor
    public  :: KillStatistic
    private ::      DeallocateInstance
    private ::      DeallocateStatisticMatrixes

    !Interfaces----------------------------------------------------------------

    private :: UngetStatistic2D
    interface  UngetStatistic
        module procedure UngetStatistic2D
    end interface UngetStatistic

    !Parameter-----------------------------------------------------------------
    integer, parameter :: Value3DStat3D_ = 1, Value3DStatLayers_ = 2, Value2DStat2D_ = 3
    integer, parameter :: Depth_ = 1, Layer_ = 2

    type T_Classification
        logical                                     :: On            = .false. !inicialization: Carina
        real                                        :: Percentil     = null_real !inicialization: Carina
        integer                                     :: nClasses      = null_int !inicialization: Carina
        real, dimension(:, :),       pointer        :: Classes       => null() !inicialization: Carina
        real, dimension(:, :, :, :), pointer        :: Frequency     => null() !inicialization: Carina
        real, dimension(:, :   , :), pointer        :: Frequency2D   => null() !inicialization: Carina
        real                                        :: RunPeriod     = null_real !inicialization: Carina
        type (T_Time)                               :: LastCalculation
    end type T_Classification

    type T_SimpleStatistic
        logical                                     :: On           = .false.
        
        real, dimension(:, :, :), pointer           :: Minimum      => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: Maximum      => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: Average      => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: SquareAverage  => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: StandardDeviation => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: GeomAverage       => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: SquareGeomAverage => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: GeomStandardDeviation  => null() !inicialization: Carina
        real, dimension(:, :, :), pointer           :: Accumulated   => null() !inicialization: Carina
        !guillaume juan
        real, dimension(:, :, :), pointer           :: PercentBelowCriticalValue  => null() !inicialization: Carina

        real, dimension(:, :   ), pointer           :: Minimum2D            => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: Maximum2D            => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: Average2D            => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: SquareAverage2D      => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: StandardDeviation2D  => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: GeomAverage2D        => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: SquareGeomAverage2D  => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: GeomStandardDeviation2D => null() !inicialization: Carina
        real, dimension(:, :   ), pointer           :: Accumulated2D       => null() !inicialization: Carina
        !guillaume juan
        real, dimension(:, :   ), pointer           :: PercentBelowCriticalValue2D => null() !inicialization: Carina

        real                                        :: RunPeriod    = null_real !inicialization: Carina
        integer                                     :: OutputNumber = null_int !inicialization: Carina
        type (T_Time)                               :: LastCalculation
        type (T_Time)                               :: NextOutputTime
    end type T_SimpleStatistic


    type T_Layers

        integer                                     :: Number = null_int !inicialization: Carina

        integer, dimension(:), pointer              :: Definition  => null() !inicialization: Carina !1 - Matrix of layer numbers
                                                                                                     !2 - Matrix of depths

        ! These two matrixes are depths and are allocated only when Methodology = 2 (Values 3D and Statistic 2D)
        real,    dimension(:,:,: ), pointer         :: UpperDepth  => null() !inicialization: Carina   
        real,    dimension(:,:,: ), pointer         :: LowerDepth  => null() !inicialization: Carina   
                                                                                                      
        !These two matrixes are layers number and are allocated only when Methodology = 2 (Values 3D and Statistic 2D)
        real,    dimension(:,:,: ), pointer         :: UpperLayer  => null() !inicialization: Carina 
        real,    dimension(:,:,: ), pointer         :: LowerLayer  => null() !inicialization: Carina

        !The average value in the layer
        real,    dimension(:,:,: ), pointer         :: Value  => null() !inicialization: Carina

        !if exist one water point between the layer limits than WaterPoints2D (i, j) = 1
        integer, dimension(:,:,: ), pointer         :: WaterPoints  => null() !inicialization: Carina

        integer, dimension(:), pointer              :: MinLayer  => null() !inicialization: Carina
        integer, dimension(:), pointer              :: MaxLayer  => null() !inicialization: Carina
        real,    dimension(:), pointer              :: MinDepth  => null() !inicialization: Carina
        real,    dimension(:), pointer              :: MaxDepth  => null() !inicialization: Carina

    end type T_Layers

    type T_ExternalVar
        type (T_Time    )                           :: Now
        type (T_Size3D  )                           :: Size
        type (T_Size3D  )                           :: WorkSize
        real,    dimension(:, :, :), pointer        :: Value3D  => null() !inicialization: Carina
        real,    dimension(:, :   ), pointer        :: Value2D  => null() !inicialization: Carina
    end type T_ExternalVar


    type T_Statistic
        integer                                     :: InstanceID = null_int !inicialization: Carina
        type (T_ExternalVar    )                    :: ExternalVar
        integer                                     :: Methodology = null_int !inicialization: Carina 
                                                                              !1 - Values 3D and Statistic 3D
                                                                              !2 - Values 3D and Statistic 2D (2D layer)
                                                                              !3 - Values 2D and Statistic 2D
        type (T_Layers )                            :: Layers 
        type (T_SimpleStatistic)                    :: Global
        type (T_SimpleStatistic)                    :: Daily
        type (T_SimpleStatistic)                    :: Monthly
        type (T_SimpleStatistic)                    :: SpecificHour
        integer                                     :: SpecificHourValue = null_int !inicialization: Carina
        type (T_Classification )                    :: Classification
        integer                                     :: ObjTime          = 0
        integer                                     :: ObjHDF5          = 0
        integer                                     :: ObjEnterData     = 0
        character(len=StringLength)                 :: Name             = null_str !inicialization: Carina
        character(len=StringLength)                 :: GroupName        = null_str !inicialization: Carina
        logical                                     :: Accumulated      = .false.  !Accumulated Values
        logical                                     :: GeomMean         = .false. !Geometric statistics
        !guillaume juan
        logical                                     :: Critical = .false.
        real                                        :: CriticalValue  = null_real !inicialization: Carina
        
        logical                                     :: NormalizeFreq = .false.  !inicialization: Carina
        
        type (T_Statistic      ), pointer           :: Next  => null() !inicialization: Carina
    end type T_Statistic


    type (T_Statistic), pointer                     :: FirstStatistic   => null()
    type (T_Statistic), pointer                     :: Me               => null()

    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine ConstructStatistic (StatisticID, ObjTime, ObjHDF5,                        &
                                   Size, WorkSize, DataFile, Name, GroupName, Rank, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: StatisticID
        integer                                     :: ObjTime
        integer                                     :: ObjHDF5
        type (T_Size3D)                             :: Size, WorkSize
        character(len=*)                            :: DataFile
        character(len=*)                            :: Name
        character(len=*), optional                  :: GroupName
        integer, optional                           :: Rank
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: STAT_, ready_
        type (T_Time)                               :: AuxTime
        real                                        :: Year,Month,Day

        !----------------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mStatistic_)) then
            nullify (FirstStatistic)
            call RegisterModule (mStatistic_) 
        endif

        write(*,*) 'constructing...'
        call Ready (StatisticID, ready_)
        write(*,*) 'constructed...'
        if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            !Name
            Me%Name = Name

            if (present(GroupName)) then
                Me%GroupName = "/Statistics/"//GroupName
            else
                Me%GroupName = "/Statistics/"
            endif

            !Associates External Instances
            Me%ObjTime           = AssociateInstance (mTIME_,           ObjTime          )
            Me%ObjHDF5           = AssociateInstance (mHDF5_,           ObjHDF5          )


            ! Actualized the time
            call GetComputeCurrentTime(Me%ObjTime,                              &
                                       Me%ExternalVar%Now, STAT = STAT_CALL)                    
            if (STAT_CALL /= SUCCESS_) stop 'ConstructStatistic - ModuleStatistic - ERR03'

            !Sets Size
            Me%ExternalVar%Size     = Size
            Me%ExternalVar%WorkSize = WorkSize

            !Reads Data File            
            call ReadDataFile (DataFile)
            
            if (present(Rank)) then
                if     (Rank == 2) then
                
                    Me%Methodology = Value2DStat2D_
                
                elseif (Rank == 3) then
                
                    Me%Methodology = Value3DStat3D_
                
                endif
            
            endif
            

            if (Me%Methodology == Value3DStatLayers_) then
                call AllocateLayerMatrixes
            endif

            !Allocates Variables 3D
            if (Me%Global%On)                                                  &
                call AllocateStatisticMatrixes (Me%Global)

            if (Me%Daily%On) then
                call AllocateStatisticMatrixes (Me%Daily)

                call ExtractDate (Time1=Me%ExternalVar%Now, Year=Year, Month=Month, Day=Day)
                
                call SetDate(Time1=AuxTime, Year=Year, Month=Month, Day=Day,          &
                             Hour=0.0, Minute=0.0, Second=0.0)
                
                Me%Daily%NextOutputTime = AuxTime + 24.*3600.
            endif
    
            if (Me%Monthly%On) then
                call AllocateStatisticMatrixes (Me%Monthly)
            endif

            if (Me%SpecificHour%On)                                            &
                call AllocateStatisticMatrixes (Me%SpecificHour)

            if (Me%Classification%On)                                          &
                Call AllocateFrequencyMatrixes

            !Returns Statistic InstanceID
            StatisticID = Me%InstanceID

        else
            write(*,*) Name

            stop 'ModuleStatistic - ConstructStatistic - ERR99' 

        endif

        STAT_ = SUCCESS_

        if (present(STAT)) STAT = STAT_


    end subroutine ConstructStatistic

    !--------------------------------------------------------------------------

    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Statistic), pointer                 :: NewStatistic
        type (T_Statistic), pointer                 :: PreviousStatistic


        !Allocates new instance
        allocate (NewStatistic)
        nullify  (NewStatistic%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstStatistic)) then
            FirstStatistic          => NewStatistic
            Me                      => NewStatistic
        else
            PreviousStatistic       => FirstStatistic
            Me                      => FirstStatistic%Next
            do while (associated(Me))
                PreviousStatistic   => Me
                Me                  => Me%Next
            enddo
            Me                      => NewStatistic
            PreviousStatistic%Next  => NewStatistic
        endif

        Me%InstanceID = RegisterNewInstance (mSTATISTIC_)

    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine ReadDataFile (DataFile)

        !Arguments-------------------------------------------------------------
        character (len=*)                           :: DataFile

        !Local-----------------------------------------------------------------
        integer                                     :: FromFile, iflag
        integer                                     :: STAT_CALL
        integer                                     :: ClientNumber    
        logical                                     :: BlockFound
        integer                                     :: FirstLine, LastLine, iLine
        integer                                     :: nClasses
        integer                                     :: ILB, IUB    
        integer                                     :: JLB, JUB    
        integer                                     :: KLB, KUB    

        !Shorten
        ILB = Me%ExternalVar%Size%ILB
        IUB = Me%ExternalVar%Size%IUB
        JLB = Me%ExternalVar%Size%JLB
        JUB = Me%ExternalVar%Size%JUB
        KLB = Me%ExternalVar%Size%KLB
        KUB = Me%ExternalVar%Size%KUB

        call GetExtractType    (FromFile = FromFile)


        !Constructs the data file
        call ConstructEnterData (Me%ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR10')

        !1 - Values 3D and Statistic 3D
        !2 - Values 3D and Statistic 2D (2D layer)
        !3 - Values 2D and Statistic 2D
        call GetData(Me%Methodology,                                                    &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'METHOD_STATISTIC',                                 &
                     default      = Value3DStat3D_,                                     &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR20')



        !Global Statistic
        call GetData(Me%Global%On,                                                      &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GLOBAL_STATISTIC',                                 &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR30')


        !Daily Statistic
        call GetData(Me%Daily%On,                                                       &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DAILY_STATISTIC',                                  &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR40')


        !Monthly Statistic
        call GetData(Me%Monthly%On,                                                     &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MONTHLY_STATISTIC',                                &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR50')

        !Statistics for a specific hour of day
        call GetData(Me%SpecificHour%On,                                                &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SPECIFIC_HOUR_STATISTIC',                          &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR60')

        if (Me%SpecificHour%On) then
            !Statistics for a specific hour of day
            call GetData(Me%SpecificHourValue,                                          &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'SPECIFIC_HOUR',                                &
                         default      = 12,                                             &
                         ClientModule = 'ModuleStatistic',                              &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                  &
                call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR70')
        
        endif

        !Accumulated values
        call GetData(Me%Accumulated,                                                    &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'ACCUMULATED',                                      &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR80')


        !Geometric Average and Standard Deviation
        call GetData(Me%GeomMean,                                                       &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'GEOMETRIC_MEAN',                                   &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR90')

        !guillaume juan
        call GetData(Me%Critical,                                                       &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'CRITICAL',                                         &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR100')

        !guillaume juan
        if (Me%Critical) then
            !Statistics for a specific hour of day
            call GetData(Me%CriticalValue,                                              &
                         Me%ObjEnterData,                                               &
                         iflag,                                                         &
                         SearchType   = FromFile,                                       &
                         keyword      = 'CRITICAL_VALUE',                               &
                         default      = 0.02,                                           &
                         ClientModule = 'ModuleStatistic',                              &
                         STAT         = STAT_CALL)

            if (STAT_CALL /= SUCCESS_)                                                  &
                call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR110')
        
        endif


        call GetData(Me%NormalizeFreq,                                                  &
                     Me%ObjEnterData,                                                   &
                     iflag,                                                             &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NORMALIZE_FREQ',                                   &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleStatistic',                                  &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                      &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR120')        

        
        Me%Classification%On = .false.
        nullify (Me%Classification%Classes)
        nullify (Me%Classification%Frequency  )
        nullify (Me%Classification%Frequency2D)



        !Classification
        !The infinit loop is necessary to scan the file again. When
        !the block ("<BeginClass>", "<EndClass>") is not found (BlockFound = .false.)
        !all the variables are initialize and the Object EnterData is read to be scan again
        !(ex: scan the block "<beginlayer>, <endlayer>).
do1 :   do
        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                   &
                                    "<BeginClass>", "<EndClass>",                    &
                                    BlockFound,                                      &
                                    FirstLine = FirstLine,                           &
                                    LastLine  = LastLine,                            &
                                    STAT      = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                   &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR07')


        !Reads Classification Data
cd1:    if (BlockFound) then
            Me%Classification%On       = .true.
            nClasses                                 = (LastLine-1) - (FirstLine+1) + 1
            Me%Classification%nClasses = nClasses
            allocate (Me%Classification%Classes(nClasses, 2))
            nClasses = 0
            do iLine = FirstLine+1, LastLine-1
                nClasses = nClasses + 1
                call GetData(Me%Classification%Classes(nClasses, :),                 &
                             Me%ObjEnterData, iflag,                                 &
                             Buffer_Line  = iLine,                                   & 
                             STAT         = STAT_CALL)
                if (iflag     /=        2) stop 'ModuleStatistic - ReadDataFile - ERR08'
                if (STAT_CALL /= SUCCESS_) stop 'ModuleStatistic - ReadDataFile - ERR09'
            enddo

        else  cd1
            
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)                                               &
                call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR10')

            exit do1    !No more blocks
            
        endif cd1

        enddo do1

        if (Me%Classification%On) then

            !Percentil
            call GetData(Me%Classification%Percentil,                                    &
                         Me%ObjEnterData,                                                &
                         iflag,                                                          &
                         SearchType   = FromFile,                                        &
                         keyword      = 'PERCENTILE',                                    &
                         default      = 90.,                                             &
                         ClientModule = 'ModuleStatistic',                               &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                                   &
                call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR11')

        endif


        if (Me%Methodology == Value3DStatLayers_)                                        &
            call Construct_Layers

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                       &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - ReadDataFile - ERR12')


    
    end subroutine ReadDataFile

    !--------------------------------------------------------------------------

    subroutine AllocateLayerMatrixes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i    
        integer                                     :: JLB, JUB    

        !Shorten
        ILB = Me%ExternalVar%Size%ILB
        IUB = Me%ExternalVar%Size%IUB
        JLB = Me%ExternalVar%Size%JLB
        JUB = Me%ExternalVar%Size%JUB

 
        allocate (Me%Layers%Value      (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
        allocate (Me%Layers%WaterPoints(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
        

        Me%Layers%Value        (:,:,:) = FillValueReal
        Me%Layers%WaterPoints  (:,:,:) = FillValueInt

        allocate (Me%Layers%UpperDepth(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
        allocate (Me%Layers%LowerDepth(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))


        allocate (Me%Layers%UpperLayer(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
        allocate (Me%Layers%LowerLayer(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))


        do i=1, Me%Layers%Number

            Me%Layers%UpperDepth(:,:,i) = Me%Layers%MinDepth(i)
            Me%Layers%LowerDepth(:,:,i) = Me%Layers%MaxDepth(i) 

            Me%Layers%UpperLayer(:,:,i) = Me%Layers%MaxLayer(i)
            Me%Layers%LowerLayer(:,:,i) = Me%Layers%MinLayer(i)

        enddo

    end subroutine AllocateLayerMatrixes

    !--------------------------------------------------------------------------

    subroutine AllocateStatisticMatrixes (Statistic)

        !Arguments-------------------------------------------------------------
        type (T_SimpleStatistic)                    :: Statistic

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB    
        integer                                     :: JLB, JUB    
        integer                                     :: KLB, KUB    

        !Shorten
        ILB = Me%ExternalVar%Size%ILB
        IUB = Me%ExternalVar%Size%IUB
        JLB = Me%ExternalVar%Size%JLB
        JUB = Me%ExternalVar%Size%JUB
        KLB = Me%ExternalVar%Size%KLB
        KUB = Me%ExternalVar%Size%KUB

        Statistic%RunPeriod       = 0.
        Statistic%OutputNumber      = 1.
        Statistic%LastCalculation   = Me%ExternalVar%Now


        if   (Me%Methodology == Value3DStat3D_) then

            allocate (Statistic%Minimum          (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate (Statistic%Maximum          (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate (Statistic%Average          (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate (Statistic%SquareAverage    (ILB:IUB, JLB:JUB, KLB:KUB))
            allocate (Statistic%StandardDeviation(ILB:IUB, JLB:JUB, KLB:KUB))

            Statistic%Minimum               = - FillValueReal
            Statistic%Maximum               = + FillValueReal
            Statistic%Average               = + FillValueReal
            Statistic%SquareAverage         = + FillValueReal
            Statistic%StandardDeviation     = + FillValueReal

            if (Me%Accumulated) then           
                allocate (Statistic%Accumulated (ILB:IUB, JLB:JUB, KLB:KUB))  
                Statistic%Accumulated = 0.0
            endif
            
            if (Me%GeomMean) then !Geometric Average is to be calculated
                allocate (Statistic%GeomAverage          (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Statistic%SquareGeomAverage    (ILB:IUB, JLB:JUB, KLB:KUB))
                allocate (Statistic%GeomStandardDeviation(ILB:IUB, JLB:JUB, KLB:KUB))

                Statistic%GeomAverage           = - FillValueReal !must be positive because log
                Statistic%SquareGeomAverage     = + FillValueReal
                Statistic%GeomStandardDeviation = + FillValueReal
            endif

            !guillaume juan
            if (Me%Critical) then           
                allocate (Statistic%PercentBelowCriticalValue (ILB:IUB, JLB:JUB, KLB:KUB))  
                Statistic%PercentBelowCriticalValue = + FillValueReal
            endif

        else if   (Me%Methodology == Value3DStatLayers_) then


            allocate (Statistic%Minimum          (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
            allocate (Statistic%Maximum          (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
            allocate (Statistic%Average          (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
            allocate (Statistic%SquareAverage    (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
            allocate (Statistic%StandardDeviation(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))

            Statistic%Minimum           = - FillValueReal
            Statistic%Maximum           = + FillValueReal
            Statistic%Average           = + FillValueReal
            Statistic%SquareAverage     = + FillValueReal
            Statistic%StandardDeviation = + FillValueReal

            if (Me%Accumulated) then           
                allocate (Statistic%Accumulated (ILB:IUB, JLB:JUB, 1:Me%Layers%Number)) 
                Statistic%Accumulated = 0.0
            endif

            if (Me%GeomMean) then !Geometric Average is to be calculated
                allocate (Statistic%GeomAverage          (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
                allocate (Statistic%SquareGeomAverage    (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))
                allocate (Statistic%GeomStandardDeviation(ILB:IUB, JLB:JUB, 1:Me%Layers%Number))

                Statistic%GeomAverage           = - FillValueReal !must be positive because log
                Statistic%SquareGeomAverage     = + FillValueReal
                Statistic%GeomStandardDeviation = + FillValueReal
            endif

            !guillaume juan
            if (Me%Critical) then           
                allocate (Statistic%PercentBelowCriticalValue (ILB:IUB, JLB:JUB, 1:Me%Layers%Number))  
                Statistic%PercentBelowCriticalValue = + FillValueReal
            endif

        else if   (Me%Methodology == Value2DStat2D_) then

            allocate (Statistic%Minimum2D          (ILB:IUB, JLB:JUB))
            allocate (Statistic%Maximum2D          (ILB:IUB, JLB:JUB))
            allocate (Statistic%Average2D          (ILB:IUB, JLB:JUB))
            allocate (Statistic%SquareAverage2D    (ILB:IUB, JLB:JUB))
            allocate (Statistic%StandardDeviation2D(ILB:IUB, JLB:JUB))

            Statistic%Minimum2D           = - FillValueReal
            Statistic%Maximum2D           = + FillValueReal
            Statistic%Average2D           = + FillValueReal
            Statistic%SquareAverage2D     = + FillValueReal
            Statistic%StandardDeviation2D = + FillValueReal

            if (Me%Accumulated) then           
                allocate (Statistic%Accumulated2D (ILB:IUB, JLB:JUB))
                Statistic%Accumulated2D = 0.0
            endif
 
            if (Me%GeomMean) then !Geometric Average is to be calculated
                allocate (Statistic%GeomAverage2D          (ILB:IUB, JLB:JUB))
                allocate (Statistic%SquareGeomAverage2D    (ILB:IUB, JLB:JUB))
                allocate (Statistic%GeomStandardDeviation2D(ILB:IUB, JLB:JUB))

                Statistic%GeomAverage2D           = - FillValueReal !must be positive because log
                Statistic%SquareGeomAverage2D     = + FillValueReal
                Statistic%GeomStandardDeviation2D = + FillValueReal
            endif

            !guillaume juan
            if (Me%Critical) then           
                allocate (Statistic%PercentBelowCriticalValue2D (ILB:IUB, JLB:JUB))  
                Statistic%PercentBelowCriticalValue2D = + FillValueReal
            endif

        endif

    end subroutine AllocateStatisticMatrixes

    !--------------------------------------------------------------------------

    subroutine AllocateFrequencyMatrixes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB    
        integer                                     :: JLB, JUB    
        integer                                     :: KLB, KUB, nClasses    

        !Shorten
        ILB = Me%ExternalVar%Size%ILB
        IUB = Me%ExternalVar%Size%IUB
        JLB = Me%ExternalVar%Size%JLB
        JUB = Me%ExternalVar%Size%JUB
        KLB = Me%ExternalVar%Size%KLB
        KUB = Me%ExternalVar%Size%KUB

        nClasses = Me%Classification%nClasses


        Me%Classification%LastCalculation   = Me%ExternalVar%Now
        Me%Classification%RunPeriod       = 0.


        if   (Me%Methodology == Value3DStat3D_) then

            allocate (Me%Classification%Frequency  (ILB:IUB, JLB:JUB, KLB:KUB, 1:nClasses))
                                                             

            Me%Classification%Frequency(:,:,:,:) = 0.


        else if   (Me%Methodology == Value3DStatLayers_) then


            allocate (Me%Classification%Frequency  (ILB:IUB, JLB:JUB, 1:Me%Layers%Number, 1:nClasses))
                                     

            Me%Classification%Frequency(:,:,:,:) = 0.

        else if   (Me%Methodology == Value2DStat2D_) then

                allocate (Me%Classification%Frequency2D(ILB:IUB, JLB:JUB, 1:nClasses))
                                                                 

                Me%Classification%Frequency2D(:,:,:) = 0.
        endif

    end subroutine AllocateFrequencyMatrixes

    !--------------------------------------------------------------------------

    subroutine Construct_Layers

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------

        character(LEN = StringLength), parameter :: block_begin = '<beginlayer>'
        character(LEN = StringLength), parameter :: block_end   = '<endlayer>'

        real, dimension(:), pointer              :: Aux1, Aux2, Aux3, Aux4, Aux5 

        integer                                  :: STAT_CALL
        integer                                  :: ClientNumber    
        logical                                  :: BlockFound
        integer                                  :: LayerNumber

        !Begin-------------------------------------------------------------------


        LayerNumber = 0

   do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = block_begin,          &
                                        block_end       = block_end,            &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)

cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :            if (BlockFound) then              
                    LayerNumber = LayerNumber + 1

                    if (associated (Me%Layers%MinDepth))                       &
                        deallocate (Me%Layers%MinDepth)

                    if (associated (Me%Layers%MaxDepth))                       &
                        deallocate (Me%Layers%MaxDepth) 

                    if (associated (Me%Layers%MaxLayer))                       &
                        deallocate (Me%Layers%MaxLayer)

                    if (associated (Me%Layers%MinLayer))                       &
                        deallocate(Me%Layers%MinLayer)

                    if (associated (Me%Layers%Definition))                     &
                        deallocate(Me%Layers%Definition)

                    allocate(Me%Layers%MinDepth  (LayerNumber))
                    allocate(Me%Layers%MaxDepth  (LayerNumber)) 
                    allocate(Me%Layers%MaxLayer  (LayerNumber))
                    allocate(Me%Layers%MinLayer  (LayerNumber))
                    allocate(Me%Layers%Definition(LayerNumber))

                    if (LayerNumber>1) then

                        Me%Layers%MinDepth  (1:LayerNumber-1) = Aux1(1:LayerNumber-1)
                        Me%Layers%MaxDepth  (1:LayerNumber-1) = Aux2(1:LayerNumber-1) 
                        Me%Layers%MinLayer  (1:LayerNumber-1) = Aux3(1:LayerNumber-1)
                        Me%Layers%MaxLayer  (1:LayerNumber-1) = Aux4(1:LayerNumber-1) 
                        Me%Layers%Definition(1:LayerNumber-1) = Aux5(1:LayerNumber-1) 

                        deallocate(Aux1)
                        deallocate(Aux2) 
                        deallocate(Aux3)
                        deallocate(Aux4)
                        deallocate(Aux5)

                    endif

                    ! Read the characteristics of a New Layer 
                    Call Read_Layer(LayerNumber)

                    allocate(Aux1(LayerNumber))
                    allocate(Aux2(LayerNumber)) 
                    allocate(Aux3(LayerNumber))
                    allocate(Aux4(LayerNumber))
                    allocate(Aux5(LayerNumber))

                    Aux1(:) = Me%Layers%MinDepth  (:) 
                    Aux2(:) = Me%Layers%MaxDepth  (:) 
                    Aux3(:) = Me%Layers%MinLayer  (:) 
                    Aux4(:) = Me%Layers%MaxLayer  (:) 
                    Aux5(:) = Me%Layers%Definition(:) 


                else cd2

                    if (associated (Aux1)) deallocate(Aux1)
                    if (associated (Aux2)) deallocate(Aux2) 
                    if (associated (Aux3)) deallocate(Aux3)
                    if (associated (Aux4)) deallocate(Aux4)
                    if (associated (Aux5)) deallocate(Aux5)


                    Me%Layers%Number = LayerNumber

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Subroutine Construct_Layers; ModuleStatistic. ERR01.'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop       'Construct_Layers; ModuleStatistic. ERR02'
            else cd1
                stop       'Construct_Layers; ModuleStatistic. ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_Layers

    !--------------------------------------------------------------------------

    subroutine Read_Layer(LayerNumber)

        !Arguments---------------------------------------------------------------

        integer                                  :: LayerNumber

        !Local----------------------------------------------------------------
        integer                                  :: FromBlock, iflag
        integer                                  :: STAT_CALL


        !Begin-------------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)


        !2 - Matrix of layer numbers 
        !1 - depth 
        call GetData(Me%Layers%Definition(LayerNumber),                              &                    
                     Me%ObjEnterData,                                                &
                     iflag,                                                          &
                     SearchType   = FromBlock,                                       &
                     keyword      = 'LAYER_DEFINITION',                              &
                     default      = Depth_,                                          &
                     ClientModule = 'ModuleStatistic',                               &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                   &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - Read_Layer - ERR01')


        call GetData(Me%Layers%MinLayer(LayerNumber),                                &                    
                     Me%ObjEnterData,                                                &
                     iflag,                                                          &
                     SearchType   = FromBlock,                                       &
                     keyword      = 'MIN_LAYER',                                     &
                     default      = Me%ExternalVar%WorkSize%KLB,                     &
                     ClientModule = 'ModuleStatistic',                               &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                   &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - Read_Layer - ERR02')

        call GetData(Me%Layers%MaxLayer(LayerNumber),                                &
                     Me%ObjEnterData,                                                &
                     iflag,                                                          &
                     SearchType   = FromBlock,                                       &
                     keyword      = 'MAX_LAYER',                                     &
                     default      = Me%ExternalVar%WorkSize%KUB,                     &
                     ClientModule = 'ModuleStatistic',                               &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                   &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - Read_Layer - ERR03')


        call GetData(Me%Layers%MinDepth(LayerNumber),                                &
                     Me%ObjEnterData,                                                &
                     iflag,                                                          &
                     SearchType   = FromBlock,                                       &
                     keyword      = 'MIN_DEPTH',                                     &
                     default      = 0.,                                              &
                     ClientModule = 'ModuleStatistic',                               &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                   &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - Read_Layer - ERR04')

        call GetData(Me%Layers%MaxDepth(LayerNumber),                                &
                     Me%ObjEnterData,                                                &
                     iflag,                                                          &
                     SearchType   = FromBlock,                                       &
                     keyword      = 'MAX_DEPTH',                                     &
                     default      = 10000.,                                          &
                     ClientModule = 'ModuleStatistic',                               &
                     STAT         = STAT_CALL)

        if (STAT_CALL /= SUCCESS_)                                                   &
            call SetError (FATAL_, KEYWORD_, 'ModuleStatistic - Read_Layer - ERR05')

    end subroutine Read_Layer


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyStatistic (StatisticID, Value2D, Value3D, WaterPoints2D,            &
                                WaterPoints3D, STAT)

        !Arguments---------------------------------------------------------------
        integer                                        :: StatisticID
        real,    dimension(:, :   ), pointer, optional :: Value2D
        real,    dimension(:, :, :), pointer, optional :: Value3D
        integer, dimension(:, :   ), pointer, optional :: WaterPoints2D
        integer, dimension(:, :, :), pointer, optional :: WaterPoints3D
        integer,    optional                           :: STAT

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, ready_   
        integer                                         :: STAT_CALL

        !Monitores Performance
        if (MonitorPerformance) call StartWatch ("ModuleStatistic", "ModifyStatistic")

        STAT_ = UNKNOWN_

        call Ready (StatisticID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)              
            if (STAT_CALL /= SUCCESS_) stop 'ModifyStatistic - ModuleStatistic - ERR01'

            if (Me%Methodology == Value2DStat2D_) then
            
                if (present(Value2D)) then
                
                    call StatisticAnalysis2D (Value2D, WaterPoints2D)
                    
                    STAT_ = SUCCESS_
                else

                    STAT_ = UNKNOWN_

                endif

            endif

            if (Me%Methodology == Value3DStat3D_) then
            
                if (present(Value3D)) then

                    Me%ExternalVar%Value3D => Value3D
                
                    call StatisticAnalysis3D (WaterPoints3D)

                    nullify(Me%ExternalVar%Value3D)

                    STAT_ = SUCCESS_
                else

                    STAT_ = UNKNOWN_

                endif

            endif


            if (Me%Methodology == Value3DStatLayers_) then
            
                if (present(Value3D)) then
                
                    call StatisticAnalysis3D (WaterPoints3D)

                    STAT_ = SUCCESS_
                else

                    STAT_ = UNKNOWN_

                endif

                
            endif

        else

            STAT_ = ready_

        endif

        if (MonitorPerformance) call StopWatch ("ModuleStatistic", "ModifyStatistic")
        if (present(STAT)) STAT = STAT_


    end subroutine ModifyStatistic

    !--------------------------------------------------------------------------
            
    subroutine AddStatisticLayers (StatisticID, Value3D, WaterPoints3D, DZ3D,            &
                                   LayerNumber, UpperDepth, LowerDepth, UpperLayer,      &
                                   LowerLayer, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: StatisticID
        real,    dimension(:, :, :), pointer            :: Value3D
        integer, dimension(:, :, :), pointer            :: WaterPoints3D
        real,    dimension(:, :, :), pointer            :: DZ3D
        integer                                         :: LayerNumber
        real,    dimension(:, :   ), pointer, optional  :: UpperDepth, LowerDepth
        integer, dimension(:, :   ), pointer, optional  :: UpperLayer, LowerLayer

        integer,    optional                            :: STAT

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, ready_    
        integer                                         :: STAT_CALL

        !Monitores Performance
        if (MonitorPerformance) call StartWatch ("ModuleStatistic", "ModifyStatistic")

        STAT_ = UNKNOWN_

        call Ready (StatisticID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            call GetComputeCurrentTime(Me%ObjTime, Me%ExternalVar%Now, STAT = STAT_CALL)              
            if (STAT_CALL /= SUCCESS_) stop 'AddStatisticLayers - ModuleStatistic - ERR01'

            if (Me%Methodology == Value3DStatLayers_) then

                    if      (Me%Layers%Definition(LayerNumber) == Depth_) then
                        
                        if (present(UpperDepth))                                         &
                            Me%Layers%UpperDepth(:,:, LayerNumber) = UpperDepth(:,:)
                        if (present(LowerDepth))                                         &
                            Me%Layers%LowerDepth(:,:, LayerNumber) = LowerDepth(:,:) 

                        call AverageValueBetweenDepths(Value3D, WaterPoints3D, DZ3D, LayerNumber)

                    else if (Me%Layers%Definition(LayerNumber) == Layer_) then

                        if (present(UpperLayer))                                         &
                            Me%Layers%UpperLayer(:,:, LayerNumber) = UpperLayer(:,:)
                        if (present(LowerLayer))                                         &
                            Me%Layers%LowerLayer(:,:, LayerNumber) = LowerLayer(:,:) 

                        call AverageValueBetweenLayers(Value3D, WaterPoints3D, DZ3D, LayerNumber)

                    endif

                STAT_ = SUCCESS_

            else

                STAT_ = UNKNOWN_

                
            endif
            
        else

            STAT_ = ready_

        endif

        if (MonitorPerformance) call StopWatch ("ModuleStatistic", "ModifyStatistic")
        if (present(STAT)) STAT = STAT_


    end subroutine AddStatisticLayers

    !--------------------------------------------------------------------------
            
    subroutine AverageValueBetweenDepths(Value3D, WaterPoints3D, DZ3D, LayerNumber)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer           :: Value3D
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real, dimension(:, :, :), pointer           :: DZ3D
        integer                                     :: LayerNumber

        !Local-----------------------------------------------------------------
        real                                        :: DepthTopLayer, DepthBottomLayer,  &
                                                       DZ_Total, DZ, DepthMin, DepthMax
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k

        !Begin-----------------------------------------------------------------

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB
        KLB = Me%ExternalVar%WorkSize%KLB
        KUB = Me%ExternalVar%WorkSize%KUB


        do j = JLB, JUB
        do i = ILB, IUB
        
            DepthTopLayer    = 0.
            DepthBottomLayer = 0.
            DZ_Total         = 0.
            Me%Layers%WaterPoints (i, j, LayerNumber) = 0
            Me%Layers%Value       (i, j, LayerNumber) = 0


            do k = KUB, KLB, -1

            DepthMin = Me%Layers%UpperDepth(i, j, LayerNumber)
            DepthMax = Me%Layers%LowerDepth(i, j, LayerNumber)


            if (WaterPoints3D(i,j,k) == WaterPoint) then

                DepthTopLayer    = DepthBottomLayer
                DepthBottomLayer = DepthBottomLayer + DZ3D(i, j, k)

                !If the layer intercept the vertical interval predefined than an average value can be compute
                if (DepthTopLayer <= DepthMax .and. DepthBottomLayer >= DepthMin) then

                    Me%Layers%WaterPoints(i, j, LayerNumber) = 1

                    if (DepthTopLayer >= DepthMin .and. DepthBottomLayer <= DepthMax) DZ = DZ3D(i, j, k)

                    if (DepthTopLayer <  DepthMin .and. DepthBottomLayer <= DepthMax) DZ = DepthBottomLayer - DepthMin

                    if (DepthTopLayer >= DepthMin .and. DepthBottomLayer >  DepthMax) DZ = DepthMax - DepthTopLayer

                    if ((DZ_Total + DZ)>0 ) then

                         Me%Layers%Value(i, j, LayerNumber) =                   &
                        (Me%Layers%Value(i, j, LayerNumber) * DZ_Total +        &
                         Value3D(i, j, k) * DZ) / (DZ_Total + DZ)

                    else

                        Me%Layers%Value(i, j, LayerNumber) = Value3D(i, j, k)

                    endif

                endif

                
                DZ_Total = DZ_Total + DZ

            endif

            enddo

        enddo
        enddo


    end subroutine AverageValueBetweenDepths

    !--------------------------------------------------------------------------

    subroutine AverageValueBetweenLayers(Value3D, WaterPoints3D, DZ3D, LayerNumber)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :, :), pointer           :: Value3D
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        real, dimension(:, :, :), pointer           :: DZ3D
        integer                                     :: LayerNumber

        !Local-----------------------------------------------------------------
        real                                        :: DepthTopLayer, DepthBottomLayer,  &
                                                       DZ_Total
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KUB, k, kup, klo

        !Begin-----------------------------------------------------------------

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB
        KUB = Me%ExternalVar%WorkSize%KUB


        do j = JLB, JUB
        do i = ILB, IUB

            if (WaterPoints3D(i,j,KUB) == WaterPoint) then
        
                DepthTopLayer    = 0.
                DepthBottomLayer = 0.
                DZ_Total         = 0.
                Me%Layers%WaterPoints (i, j, LayerNumber) = 0
                Me%Layers%Value       (i, j, LayerNumber) = 0

                kup = Me%Layers%UpperLayer(i, j, LayerNumber)
                klo = Me%Layers%LowerLayer(i, j, LayerNumber)

                do k = kup, klo, -1

                if (WaterPoints3D(i,j,k) == WaterPoint) then
                
                    Me%Layers%WaterPoints(i, j, LayerNumber) = 1

                    Me%Layers%Value(i, j, LayerNumber) =                           &
                        (Me%Layers%Value(i, j, LayerNumber) * DZ_Total +           &
                         Value3D(i, j, k) * DZ3D(i, j, k)) / (DZ_Total + DZ3D(i, j, k))

            
                    DZ_Total = DZ_Total + DZ3D(i, j, k)

                endif

                enddo

            endif

        enddo
        enddo

    end subroutine AverageValueBetweenLayers

    !--------------------------------------------------------------------------

    subroutine StatisticAnalysis3D (WaterPoints3D)

        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :), pointer        :: WaterPoints3D

        !Local-----------------------------------------------------------------
        real,    dimension(:, :, :), pointer        :: Value3D
        integer                                     :: KLB, KUB
        integer, dimension(:, :, :), pointer        :: lWaterPoints3D

        !Begin-----------------------------------------------------------------

            if (Me%Methodology == Value3DStat3D_) then

                Value3D        => Me%ExternalVar%Value3D
                lWaterPoints3D => WaterPoints3D
                KLB            =  Me%ExternalVar%WorkSize%KLB
                KUB            =  Me%ExternalVar%WorkSize%KUB

            else if (Me%Methodology == Value3DStatLayers_) then

                Value3D        => Me%Layers%Value
                lWaterPoints3D => Me%Layers%WaterPoints
                KLB            =  1
                KUB            =  Me%Layers%Number

            endif

            !Global
            if (Me%Global%On) then
                call ModifyGlobalStatistic  (Value3D, lWaterPoints3D, KLB, KUB)
            endif

            !Daily
            if (Me%Daily%On) then
                call ModifyDailyStatistic   (Value3D, lWaterPoints3D, KLB, KUB)
            endif
            
            !Monthly
            if (Me%Monthly%On) then
                call ModifyMonthlyStatistic (Value3D, lWaterPoints3D, KLB, KUB)
            endif

           !SpecificHour
            if (Me%SpecificHour%On) then
                call ModifySpecificHourStatistic (Value3D, lWaterPoints3D, KLB, KUB)
            endif

            !Classification
            if (Me%Classification%On) then
                call ModifyClassification   (Value3D, lWaterPoints3D, KLB, KUB)
            endif

            nullify(Value3D)
            nullify(lWaterPoints3D)

    end subroutine StatisticAnalysis3D

    !--------------------------------------------------------------------------

    subroutine StatisticAnalysis2D (Value2D, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :), pointer              :: Value2D
        integer, dimension(:, :), pointer           :: WaterPoints2D

        !Local-----------------------------------------------------------------

            !Global
            if (Me%Global%On) then
                call ModifyGlobalStatistic2D  (Value2D, WaterPoints2D)
            endif

            !Daily
            if (Me%Daily%On) then
                call ModifyDailyStatistic2D   (Value2D, WaterPoints2D)
            endif
            
            !Monthly
            if (Me%Monthly%On) then
                call ModifyMonthlyStatistic2D (Value2D, WaterPoints2D)
            endif

            !SpecificHour
            if (Me%SpecificHour%On) then
                call ModifySpecificHourStatistic2D (Value2D, WaterPoints2D)
            endif
            !Classification
            if (Me%Classification%On) then
                call ModifyClassification2D   (Value2D, WaterPoints2D)
            endif

    end subroutine StatisticAnalysis2D

    !--------------------------------------------------------------------------

    subroutine ModifyGlobalStatistic (Value, WaterPoints3D, KLB, KUB)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer        :: Value
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        integer                                     :: KLB, KUB

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: k
        real                                        :: DT, DX, AuxValue

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB


        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Global%LastCalculation

cd1:    if (DT>0) then

        !Loops
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints3D(i, j, k) == WaterPoint) then

                !Minimum Value
                if (Value (i, j, k) < Me%Global%Minimum (i, j, k))          &
                    Me%Global%Minimum (i, j, k) = Value (i, j, k)

                !Maximum Value
                if (Value (i, j, k) > Me%Global%Maximum (i, j, k))          &
                    Me%Global%Maximum (i, j, k) = Value (i, j, k)

                !Average
                Me%Global%Average (i, j, k) =                               &
                    (Me%Global%Average (i, j, k) *                          &
                     Me%Global%RunPeriod      +                             &
                     Value (i, j, k) * DT) / (Me%Global%RunPeriod + DT)

                !Square Average
                Me%Global%SquareAverage (i, j, k) =                         &
                    (Me%Global%SquareAverage (i, j, k) *                    &
                     Me%Global%RunPeriod      +                             &
                     Value (i, j, k)**2 * DT) / (Me%Global%RunPeriod + DT)

                !Standard deviation
                DX = Me%Global%SquareAverage (i, j, k) -                    &
                     Me%Global%Average       (i, j, k) ** 2

                DX = abs(DX) 

                Me%Global%StandardDeviation(i, j, k) = sqrt(DX)

                !Accumulated Values
                if (Me%Accumulated)                                                     &
                    Me%Global%Accumulated (i,j,k) = Me%Global%Accumulated (i,j,k)       &
                                                  + Value (i, j, k) 

                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value (i, j, k)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifyGlobalStatistic - ModuleStatistic - ERR01'
                    endif
                 
                    Me%Global%GeomAverage (i, j, k) = 10**                  &
                        ((LOG10(Me%Global%GeomAverage (i, j, k)) *          &
                        Me%Global%RunPeriod      +                          &
                        LOG10(AuxValue) * DT) / (Me%Global%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%Global%SquareGeomAverage (i, j, k) =                 &
                        (Me%Global%SquareGeomAverage (i, j, k) *            &
                        Me%Global%RunPeriod      +                          &
                        AuxValue**2 * DT) / (Me%Global%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%Global%SquareGeomAverage (i, j, k) -            &
                         Me%Global%GeomAverage       (i, j, k) ** 2

                    DX = abs(DX) 

                    Me%Global%GeomStandardDeviation(i, j, k) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical .and. (Value (i, j, k) > FillValueReal / 2.)) then
                    if (Value (i, j, k) < Me%CriticalValue) then
                        Me%Global%PercentBelowCriticalValue (i, j, k) =                               &
                            (Me%Global%PercentBelowCriticalValue (i, j, k) *                          &
                             Me%Global%RunPeriod + DT) / (Me%Global%RunPeriod + DT)
                    else
                        Me%Global%PercentBelowCriticalValue (i, j, k) =                               &
                            (Me%Global%PercentBelowCriticalValue (i, j, k) *                          &
                             Me%Global%RunPeriod + 0.) / (Me%Global%RunPeriod + DT)
                    endif
                endif

            endif
        enddo
        enddo
        enddo

        !Updates Time
        Me%Global%RunPeriod     = Me%Global%RunPeriod + DT
        Me%Global%LastCalculation = Me%ExternalVar%Now

        endif cd1

    end subroutine ModifyGlobalStatistic

    !--------------------------------------------------------------------------

    subroutine ModifyDailyStatistic (Value, WaterPoints3D, KLB, KUB)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer        :: Value
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        integer                                     :: KLB, KUB

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: k
        real                                        :: DT, DX, AuxValue
        real                                        :: OldDay, PresentDay    

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB


        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Daily%LastCalculation

cd1:    if (DT>0) then

        !Loops
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints3D(i, j, k) == WaterPoint) then

                !Minimum Value
                if (Value (i, j, k) < Me%Daily%Minimum (i, j, k))           &
                    Me%Daily%Minimum (i, j, k) = Value (i, j, k)

                !Maximum Value
                if (Value (i, j, k) > Me%Daily%Maximum (i, j, k))           &
                    Me%Daily%Maximum (i, j, k) = Value (i, j, k)

                !Average
                Me%Daily%Average (i, j, k) =                                &
                    (Me%Daily%Average (i, j, k) *                           &
                     Me%Daily%RunPeriod      +                              &
                     Value (i, j, k) * DT) / (Me%Daily%RunPeriod + DT)

                !Square Average
                Me%Daily%SquareAverage (i, j, k) =                          &
                    (Me%Daily%SquareAverage (i, j, k) *                     &
                     Me%Daily%RunPeriod      +                              &
                     Value (i, j, k)**2 * DT) / (Me%Daily%RunPeriod + DT)

                !Standard deviation
                DX = Me%Daily%SquareAverage (i, j, k) -                     &
                     Me%Daily%Average       (i, j, k) ** 2

                DX = abs(DX) 

                Me%Daily%StandardDeviation(i, j, k) = sqrt(DX)

                !Accumulated Values
                if (Me%Accumulated)                                                     &
                    Me%Daily%Accumulated (i,j,k) = Me%Daily%Accumulated (i,j,k)         &
                                                  + Value (i, j, k) 

                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value (i, j, k)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifyDailyStatistic - ModuleStatistic - ERR01'
                    endif
                  
                    if (Me%Daily%GeomAverage (i, j, k) == 0.0) then
                        Me%Daily%GeomAverage (i, j, k) = 1.0
                    endif

                    Me%Daily%GeomAverage (i, j, k) = 10**                   &
                        ((LOG10(Me%Daily%GeomAverage (i, j, k)) *           &
                        Me%Daily%RunPeriod      +                           &
                        LOG10(AuxValue) * DT) / (Me%Daily%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%Daily%SquareGeomAverage (i, j, k) =                  &
                        (Me%Daily%SquareGeomAverage (i, j, k) *             &
                        Me%Daily%RunPeriod      +                           &
                        AuxValue**2 * DT) / (Me%Daily%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%Daily%SquareGeomAverage (i, j, k) -             &
                         Me%Daily%GeomAverage       (i, j, k) ** 2

                    DX = abs(DX) 

                    Me%Daily%GeomStandardDeviation(i, j, k) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical .and. (Value (i, j, k) > FillValueReal / 2.)) then
                    if (Value (i, j, k) < Me%CriticalValue) then
                        Me%Daily%PercentBelowCriticalValue (i, j, k) =                               &
                            (Me%Daily%PercentBelowCriticalValue (i, j, k) *                          &
                             Me%Daily%RunPeriod + DT) / (Me%Daily%RunPeriod + DT)
                    else
                        Me%Daily%PercentBelowCriticalValue (i, j, k) =                               &
                            (Me%Daily%PercentBelowCriticalValue (i, j, k) *                          &
                             Me%Daily%RunPeriod + 0.) / (Me%Daily%RunPeriod + DT)
                    endif
                endif

            endif
        enddo
        enddo
        enddo

        !Verifies if the present time is a new output
        call ExtractDate (Time1=Me%ExternalVar%Now,       Day = PresentDay)
        call ExtractDate (Time1=Me%Daily%LastCalculation, Day = OldDay)
        if (int(PresentDay) /= int(OldDay)) then
            call WriteValuesToFileHDF5 (.false., .true., .false., .false., .false.)
            Me%Daily%Minimum           = Value
            Me%Daily%Maximum           = Value
            Me%Daily%Average           = Value
            Me%Daily%SquareAverage     = Value
            Me%Daily%StandardDeviation = Value
  
            if (Me%Accumulated) Me%Daily%Accumulated = 0.0
            
            if (Me%GeomMean) then !Geometric Average to calculate
                Me%Daily%GeomAverage           = Value
                Me%Daily%SquareGeomAverage     = Value
                Me%Daily%GeomStandardDeviation = Value
            endif

            !guillaume juan
            if (Me%Critical) then
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if ( (WaterPoints3D(i, j, k) == WaterPoint) .and. (Value (i, j, k) > FillValueReal / 2.)) then
                        if(Value(i, j, k) < Me%CriticalValue) then
                            Me%Daily%PercentBelowCriticalValue (i, j, k) = 1.
                        else
                            Me%Daily%PercentBelowCriticalValue (i, j, k) = 0.
                        endif
                    endif
                enddo
                enddo
                enddo
            endif

            Me%Daily%RunPeriod  = 0.
        else
            Me%Daily%RunPeriod  = Me%Daily%RunPeriod + DT
        endif

        !Updates Time
        Me%Daily%LastCalculation  = Me%ExternalVar%Now


        endif cd1


    end subroutine ModifyDailyStatistic

    !--------------------------------------------------------------------------

    subroutine ModifyMonthlyStatistic (Value, WaterPoints3D, KLB, KUB)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer        :: Value
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        integer                                     :: KLB, KUB

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: k
        real                                        :: DT, DX, AuxValue
        real                                        :: OldMonth, PresentMonth

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB


        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Monthly%LastCalculation

cd1:    if (DT>0) then

        !Loops
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints3D(i, j, k) == WaterPoint) then

                !Minimum Value
                if (Value (i, j, k) < Me%Monthly%Minimum (i, j, k))         &
                    Me%Monthly%Minimum (i, j, k) = Value (i, j, k)

                !Maximum Value
                if (Value (i, j, k) > Me%Monthly%Maximum (i, j, k))         &
                    Me%Monthly%Maximum (i, j, k) = Value (i, j, k)

                !Average
                Me%Monthly%Average (i, j, k) =                              &
                    (Me%Monthly%Average (i, j, k) *                         &
                     Me%Monthly%RunPeriod      +                            &
                     Value (i, j, k) * DT) / (Me%Monthly%RunPeriod + DT)

                !Square Average
                Me%Monthly%SquareAverage (i, j, k) =                        &
                    (Me%Monthly%SquareAverage (i, j, k) *                   &
                     Me%Monthly%RunPeriod      +                            &
                     Value (i, j, k)**2 * DT) / (Me%Monthly%RunPeriod + DT)

                !Standard deviation
                DX = Me%Monthly%SquareAverage (i, j, k) -                   &
                     Me%Monthly%Average       (i, j, k) ** 2

                DX = abs(DX) 


                Me%Monthly%StandardDeviation(i, j, k) = sqrt(DX)
   
   
                !Accumulated Values
                if (Me%Accumulated)                                                     &
                    Me%Monthly%Accumulated (i,j,k) = Me%Monthly%Accumulated (i,j,k)     &
                                                   + Value (i, j, k) 


                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value (i, j, k)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifyMonthlyStatistic - ModuleStatistic - ERR01'
                    endif

                    if (Me%Monthly%GeomAverage (i, j, k) == 0.0) then
                        Me%Monthly%GeomAverage (i, j, k) = 1.0
                    endif
                  
                    Me%Monthly%GeomAverage (i, j, k) = 10**                 &
                        ((LOG10(Me%Monthly%GeomAverage (i, j, k)) *         &
                        Me%Monthly%RunPeriod      +                         &
                        LOG10(AuxValue) * DT) / (Me%Monthly%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%Monthly%SquareGeomAverage (i, j, k) =                &
                        (Me%Monthly%SquareGeomAverage (i, j, k) *           &
                        Me%Monthly%RunPeriod      +                         &
                        AuxValue**2 * DT) / (Me%Monthly%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%Monthly%SquareGeomAverage (i, j, k) -           &
                         Me%Monthly%GeomAverage       (i, j, k) ** 2

                    DX = abs(DX) 

                    Me%Monthly%GeomStandardDeviation(i, j, k) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical .and. (Value (i, j, k) > FillValueReal / 2.)) then
                    if (Value (i, j, k) < Me%CriticalValue) then
                        Me%Monthly%PercentBelowCriticalValue (i, j, k) =                               &
                            (Me%Monthly%PercentBelowCriticalValue (i, j, k) *                          &
                             Me%Monthly%RunPeriod + DT) / (Me%Monthly%RunPeriod + DT)
                    else
                        Me%Monthly%PercentBelowCriticalValue (i, j, k) =                               &
                            (Me%Monthly%PercentBelowCriticalValue (i, j, k) *                          &
                             Me%Monthly%RunPeriod + 0.) / (Me%Monthly%RunPeriod + DT)
                    endif
                endif

            endif
        enddo
        enddo
        enddo

        !Verifies if the present time is a new output
        call ExtractDate (Time1=Me%ExternalVar%Now,         Month = PresentMonth)
        call ExtractDate (Time1=Me%Monthly%LastCalculation, Month = OldMonth)
        if (int(PresentMonth) /= int(OldMonth)) then
            call WriteValuesToFileHDF5 (.false., .false., .true., .false., .false.)
            Me%Monthly%Minimum           = Value
            Me%Monthly%Maximum           = Value
            Me%Monthly%Average           = Value
            Me%Monthly%SquareAverage     = Value
            Me%Monthly%StandardDeviation = Value

            if (Me%Accumulated) Me%Monthly%Accumulated = 0.0

            if (Me%GeomMean) then !Geometric Average to calculate
                Me%Monthly%GeomAverage           = Value
                Me%Monthly%SquareGeomAverage     = Value
                Me%Monthly%GeomStandardDeviation = Value
            endif

            !guillaume juan
            if (Me%Critical) then
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                    if ((WaterPoints3D(i, j, k) == WaterPoint) .and. (Value (i, j, k) > FillValueReal / 2.)) then
                        if(Value(i, j, k) < Me%CriticalValue) then
                            Me%Monthly%PercentBelowCriticalValue (i, j, k) = 1.
                        else
                            Me%Monthly%PercentBelowCriticalValue (i, j, k) = 0.
                        endif
                    endif
                enddo
                enddo
                enddo
            endif

            Me%Monthly%RunPeriod  = 0.
        else
            Me%Monthly%RunPeriod  = Me%Monthly%RunPeriod + DT
        endif

        !Updates Time
        Me%Monthly%LastCalculation  = Me%ExternalVar%Now

        endif cd1


    end subroutine ModifyMonthlyStatistic

    !--------------------------------------------------------------------------

    subroutine ModifySpecificHourStatistic (Value, WaterPoints3D, KLB, KUB)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer        :: Value
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        integer                                     :: KLB, KUB

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: k
        real                                        :: DT, DX, AuxValue
        real                                        :: PresentHour

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

        call ExtractDate (Time1=Me%ExternalVar%Now, Hour = PresentHour)
        
if1:    if (int(PresentHour) == int(Me%SpecificHourValue)) then

            !Time since last calculation 
            DT  = Me%ExternalVar%Now - Me%SpecificHour%LastCalculation

cd1:        if (DT>0) then

            !Loops
            do k = KLB, KUB
            do j = JLB, JUB
            do i = ILB, IUB
                if (WaterPoints3D(i, j, k) == WaterPoint) then

                    !Minimum Value
                    if (Value (i, j, k) < Me%SpecificHour%Minimum (i, j, k))         &
                        Me%SpecificHour%Minimum (i, j, k) = Value (i, j, k)

                    !Maximum Value
                    if (Value (i, j, k) > Me%SpecificHour%Maximum (i, j, k))         &
                        Me%SpecificHour%Maximum (i, j, k) = Value (i, j, k)

                    !Average
                    Me%SpecificHour%Average (i, j, k) =                              &
                        (Me%SpecificHour%Average (i, j, k) *                         &
                         Me%SpecificHour%RunPeriod      +                            &
                         Value (i, j, k) * DT) / (Me%SpecificHour%RunPeriod + DT)

                    !Square Average
                    Me%SpecificHour%SquareAverage (i, j, k) =                        &
                        (Me%SpecificHour%SquareAverage (i, j, k) *                   &
                         Me%SpecificHour%RunPeriod      +                            &
                         Value (i, j, k)**2 * DT) / (Me%SpecificHour%RunPeriod + DT)

                    !Standard deviation
                    DX = Me%SpecificHour%SquareAverage (i, j, k) -                   &
                         Me%SpecificHour%Average       (i, j, k) ** 2

                    DX = abs(DX) 


                    Me%SpecificHour%StandardDeviation(i, j, k) = sqrt(DX)

                    !Accumulated Values
                    if (Me%Accumulated)                                                             &
                        Me%SpecificHour%Accumulated (i,j,k) = Me%SpecificHour%Accumulated (i,j,k)   &
                                                            + Value (i, j, k) 

                    if (Me%GeomMean) then !Geometric Average to calculate
                        !Geometric Average
                        AuxValue = Value (i, j, k)
                        if (AuxValue == 0.0) then
                            AuxValue = 1.0                   
                        elseif (AuxValue .lt. 0.0) then                   
                            write(*,*) 'Negative valued property.'
                            write(*,*) 'Geometric Average cannot be calculated.'
                            stop 'ModifySpecificHourStatistic - ModuleStatistic - ERR01'
                        endif

                        if (Me%SpecificHour%GeomAverage (i, j, k) == 0.0) then
                            Me%SpecificHour%GeomAverage (i, j, k) = 1.0
                        endif
                  
                        Me%SpecificHour%GeomAverage (i, j, k) = 10**                 &
                            ((LOG10(Me%SpecificHour%GeomAverage (i, j, k)) *         &
                            Me%SpecificHour%RunPeriod      +                         &
                            LOG10(AuxValue) * DT) / (Me%SpecificHour%RunPeriod + DT))

                        !Squared Geometric Average
                        Me%SpecificHour%SquareGeomAverage (i, j, k) =                &
                            (Me%SpecificHour%SquareGeomAverage (i, j, k) *           &
                            Me%SpecificHour%RunPeriod      +                         &
                            AuxValue**2 * DT) / (Me%SpecificHour%RunPeriod + DT)

                        !Geometric Standard Deviation
                        DX = Me%SpecificHour%SquareGeomAverage (i, j, k) -           &
                             Me%SpecificHour%GeomAverage       (i, j, k) ** 2

                        DX = abs(DX) 

                        Me%SpecificHour%GeomStandardDeviation(i, j, k) = sqrt(DX)

                    endif

                    !guillaume juan
                    if (Me%Critical  .and. (Value (i, j, k) > FillValueReal / 2.)) then
                        if (Value (i, j, k) < Me%CriticalValue) then
                            Me%SpecificHour%PercentBelowCriticalValue (i, j, k) =                               &
                                (Me%SpecificHour%PercentBelowCriticalValue (i, j, k) *                          &
                                 Me%SpecificHour%RunPeriod + DT) / (Me%SpecificHour%RunPeriod + DT)
                        else
                            Me%SpecificHour%PercentBelowCriticalValue (i, j, k) =                               &
                                (Me%SpecificHour%PercentBelowCriticalValue (i, j, k) *                          &
                                 Me%SpecificHour%RunPeriod + 0.) / (Me%SpecificHour%RunPeriod + DT)
                        endif
                    endif

                endif
            enddo
            enddo
            enddo

            
            !Updates Time
            Me%SpecificHour%RunPeriod       = Me%SpecificHour%RunPeriod + DT          
            Me%SpecificHour%LastCalculation = Me%ExternalVar%Now

        endif cd1

    endif if1

    end subroutine ModifySpecificHourStatistic

    !--------------------------------------------------------------------------

    subroutine ModifyClassification (Value, WaterPoints3D, KLB, KUB)

        !Arguments-------------------------------------------------------------
        real,    dimension(:, :, :), pointer        :: Value
        integer, dimension(:, :, :), pointer        :: WaterPoints3D
        integer                                     :: KLB, KUB


        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: k
        integer                                     :: iClass   
        real                                        :: DT, Aux

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB


        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Classification%LastCalculation

cd1:    if (DT>0) then
        
        !Loops
        do k = KLB, KUB
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints3D(i, j, k) == WaterPoint) then
doClass:        do iClass = 1, Me%Classification%nClasses
                    if (Value(i, j, k) >= Me%Classification%Classes(iClass, 1) .and.    &
                        Value(i, j, k)  < Me%Classification%Classes(iClass, 2)) then
                        Aux = DT
                    else
                        Aux = 0
                    endif
                    
                    Me%Classification%Frequency     (i, j, k, iClass) =                 &
                        (Me%Classification%Frequency(i, j, k, iClass) *                 &
                         Me%Classification%RunPeriod + Aux          ) /                 &
                        (Me%Classification%RunPeriod + DT)
                enddo doClass
            endif
        enddo
        enddo
        enddo


        Me%Classification%RunPeriod       = Me%Classification%RunPeriod + DT

        Me%Classification%LastCalculation = Me%ExternalVar%Now

        endif cd1

    end subroutine ModifyClassification

    !--------------------------------------------------------------------------

    subroutine ModifyGlobalStatistic2D (Value2D, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:),   pointer          :: Value2D
        integer, dimension(:,:),   pointer          :: WaterPoints2D

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        real                                        :: DT, DX, AuxValue

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

       
        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Global%LastCalculation

cd1:    if (DT>0) then

        !Loops
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint) then


                !Minimum Value2D
                if (Value2D (i, j) < Me%Global%Minimum2D (i, j))            &
                    Me%Global%Minimum2D (i, j) = Value2D (i, j)

                !Maximum Value2D
                if (Value2D (i, j) > Me%Global%Maximum2D (i, j))            &
                    Me%Global%Maximum2D (i, j) = Value2D (i, j)

                !Average
                Me%Global%Average2D (i, j) =                                &
                    (Me%Global%Average2D (i, j) *                           &
                     Me%Global%RunPeriod      +                             &
                     Value2D (i, j) * DT) / (Me%Global%RunPeriod + DT)

                !Square Average
                Me%Global%SquareAverage2D (i, j) =                          &
                    (Me%Global%SquareAverage2D (i, j) *                     &
                     Me%Global%RunPeriod      +                             &
                     Value2D (i, j)**2 * DT) / (Me%Global%RunPeriod + DT)

                !Standard deviation
                DX = Me%Global%SquareAverage2D (i, j) -                     &
                     Me%Global%Average2D       (i, j) ** 2

                DX = abs(DX) 

                Me%Global%StandardDeviation2D(i, j) = sqrt(DX)


                !Accumulated Values
                if (Me%Accumulated)                                                     &
                    Me%Global%Accumulated2D (i,j) = Me%Global%Accumulated2D (i,j)       &
                                                  + Value2D (i, j) 

                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value2D (i, j)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifyGlobalStatistic2D - ModuleStatistic - ERR01'
                    endif
                
                    Me%Global%GeomAverage2D (i, j) = 10**                   &
                        ((LOG10(Me%Global%GeomAverage2D(i, j)) *            &
                        Me%Global%RunPeriod      +                          &
                        LOG10(AuxValue) * DT) / (Me%Global%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%Global%SquareGeomAverage2D (i, j) =                  &
                        (Me%Global%SquareGeomAverage2D (i, j) *             &
                        Me%Global%RunPeriod      +                          &
                        AuxValue**2 * DT) / (Me%Global%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%Global%SquareGeomAverage2D (i, j) -             &
                         Me%Global%GeomAverage2D       (i, j) ** 2

                    DX = abs(DX) 

                    Me%Global%GeomStandardDeviation2D(i, j) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical  .and. (Value2D (i, j) > FillValueReal / 2.)) then
                    if (Value2D (i, j) < Me%CriticalValue) then
                        Me%Global%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%Global%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%Global%RunPeriod + DT) / (Me%Global%RunPeriod + DT)                        
                    else
                        Me%Global%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%Global%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%Global%RunPeriod + 0.) / (Me%Global%RunPeriod + DT)                        
                    endif
                endif

            endif
        enddo
        enddo


        !Updates Time
        Me%Global%RunPeriod     = Me%Global%RunPeriod + DT
        Me%Global%LastCalculation = Me%ExternalVar%Now

      
        endif cd1

    end subroutine ModifyGlobalStatistic2D

    !--------------------------------------------------------------------------

    subroutine ModifyDailyStatistic2D (Value2D, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:),   pointer          :: Value2D
        integer, dimension(:,:),   pointer          :: WaterPoints2D

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        real                                        :: DT, DX, AuxValue
!        real                                        :: OldDay, PresentDay    

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB
        
        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Daily%LastCalculation

cd1:    if (DT>0) then


        !Loops
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint) then

                !Minimum Value2D
                if (Value2D (i, j) < Me%Daily%Minimum2D (i, j))             &
                    Me%Daily%Minimum2D (i, j) = Value2D (i, j)

                !Maximum Value2D
                if (Value2D (i, j) > Me%Daily%Maximum2D (i, j))             &
                    Me%Daily%Maximum2D (i, j) = Value2D (i, j)

                !Average
                Me%Daily%Average2D (i, j) =                                 &
                    (Me%Daily%Average2D (i, j) *                            &
                     Me%Daily%RunPeriod      +                              &
                     Value2D (i, j) * DT) / (Me%Daily%RunPeriod + DT)

                !Square Average
                Me%Daily%SquareAverage2D (i, j) =                           &
                    (Me%Daily%SquareAverage2D (i, j) *                      &
                     Me%Daily%RunPeriod      +                              &
                     Value2D (i, j)**2 * DT) / (Me%Daily%RunPeriod + DT)

                !Standard deviation
                DX = Me%Daily%SquareAverage2D (i, j) -                      &
                     Me%Daily%Average2D       (i, j) ** 2

                DX = abs(DX) 

                Me%Daily%StandardDeviation2D(i, j) = sqrt(DX)


                !Accumulated Values
                if (Me%Accumulated)                                                     &
                    Me%Daily%Accumulated2D (i,j) = Me%Daily%Accumulated2D (i,j)         &
                                                  + Value2D (i, j) 

                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value2D (i, j)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifyDailyStatistic2D - ModuleStatistic - ERR01'
                    endif

                    if (Me%Daily%GeomAverage2D (i, j) == 0.0) then
                        Me%Daily%GeomAverage2D (i, j) = 1.0
                    endif
                 
                    Me%Daily%GeomAverage2D (i, j) = 10**                    &
                        ((LOG10(Me%Daily%GeomAverage2D (i, j)) *            &
                        Me%Daily%RunPeriod      +                           &
                        LOG10(AuxValue) * DT) / (Me%Daily%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%Daily%SquareGeomAverage2D (i, j) =                   &
                        (Me%Daily%SquareGeomAverage2D (i, j) *              &
                        Me%Daily%RunPeriod      +                           &
                        AuxValue**2 * DT) / (Me%Daily%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%Daily%SquareGeomAverage2D (i, j) -              &
                         Me%Daily%GeomAverage2D       (i, j) ** 2

                    DX = abs(DX) 

                    Me%Daily%GeomStandardDeviation2D(i, j) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical .and. (Value2D (i, j) > FillValueReal / 2.)) then
                    if (Value2D (i, j) < Me%CriticalValue) then
                        Me%Daily%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%Daily%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%Daily%RunPeriod + DT) / (Me%Daily%RunPeriod + DT)                        
                    else
                        Me%Daily%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%Daily%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%Daily%RunPeriod + 0.) / (Me%Daily%RunPeriod + DT)                        
                    endif
                endif

            endif
        enddo
        enddo

        !Verifies if the present time is a new output

!        call ExtractDate (Me%ExternalVar%Now,       Day = PresentDay)
!        call ExtractDate (Me%Daily%LastCalculation, Day = OldDay)
!        if (int(PresentDay) /= int(OldDay)) then

        if (Me%ExternalVar%Now .GT. Me%Daily%NextOutputTime) then

            call WriteValuesToFileHDF5 (.false., .true., .false., .false., .false.)
            Me%Daily%Minimum2D           = Value2D
            Me%Daily%Maximum2D           = Value2D
            Me%Daily%Average2D           = Value2D
            Me%Daily%SquareAverage2D     = Value2D
            Me%Daily%StandardDeviation2D = Value2D

            if (Me%Accumulated) Me%Daily%Accumulated2D = 0.0

            if (Me%GeomMean) then !Geometric Average to calculate
                Me%Daily%GeomAverage2D           = Value2D
                Me%Daily%SquareGeomAverage2D     = Value2D
                Me%Daily%GeomStandardDeviation2D = Value2D
            endif

            !guillaume juan
            if (Me%Critical) then
                do j = JLB, JUB
                do i = ILB, IUB
                    if ( (WaterPoints2D(i, j) == WaterPoint) .and. (Value2D (i, j) > FillValueReal / 2.)) then
                        if(Value2D(i, j) < Me%CriticalValue) then
                            Me%Daily%PercentBelowCriticalValue2D (i, j) = 1.
                        else
                            Me%Daily%PercentBelowCriticalValue2D (i, j) = 0.
                        endif
                    endif
                enddo
                enddo
            endif

            Me%Daily%RunPeriod  = 0.
            Me%Daily%NextOutputTime = Me%Daily%NextOutputTime + 24*3600.
        else
            Me%Daily%RunPeriod  = Me%Daily%RunPeriod + DT
        endif

        !Updates Time
        Me%Daily%LastCalculation  = Me%ExternalVar%Now


        endif cd1


    end subroutine ModifyDailyStatistic2D

    !--------------------------------------------------------------------------

    subroutine ModifyMonthlyStatistic2D (Value2D, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:),   pointer          :: Value2D
        integer, dimension(:,:),   pointer          :: WaterPoints2D

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        real                                        :: DT, DX, AuxValue
        real                                        :: OldMonth, PresentMonth

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB



        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Monthly%LastCalculation

cd1:    if (DT>0) then

        !Loops
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint) then

                !Minimum Value2D
                if (Value2D (i, j) < Me%Monthly%Minimum2D (i, j))           &
                    Me%Monthly%Minimum2D (i, j) = Value2D (i, j)

                !Maximum Value2D
                if (Value2D (i, j) > Me%Monthly%Maximum2D (i, j))           &
                    Me%Monthly%Maximum2D (i, j) = Value2D (i, j)

                !Average
                Me%Monthly%Average2D (i, j) =                               &
                    (Me%Monthly%Average2D (i, j) *                          &
                     Me%Monthly%RunPeriod      +                            &
                     Value2D (i, j) * DT) / (Me%Monthly%RunPeriod + DT)

                !Square Average
                Me%Monthly%SquareAverage2D (i, j) =                         &
                    (Me%Monthly%SquareAverage2D (i, j) *                    &
                     Me%Monthly%RunPeriod      +                            &
                     Value2D (i, j)**2 * DT) / (Me%Monthly%RunPeriod + DT)

                !Standard deviation
                DX = Me%Monthly%SquareAverage2D (i, j) -                    &
                     Me%Monthly%Average2D       (i, j) ** 2

                DX = abs(DX) 

                Me%Monthly%StandardDeviation2D(i, j) = sqrt(DX)

                !Accumulated Values
                if (Me%Accumulated)                                                     &
                    Me%Monthly%Accumulated2D (i,j) = Me%Monthly%Accumulated2D (i,j)     &
                                                   + Value2D (i, j) 

                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value2D (i, j)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifyMonthlyStatistic2D - ModuleStatistic - ERR01'
                    endif

                    if (Me%Monthly%GeomAverage2D (i, j) == 0.0) then
                        Me%Monthly%GeomAverage2D (i, j) = 1.0
                    endif
                                   
                    Me%Monthly%GeomAverage2D (i, j) = 10**                  &
                        ((LOG10(Me%Monthly%GeomAverage2D (i, j)) *          &
                        Me%Monthly%RunPeriod      +                         &
                        LOG10(AuxValue) * DT) / (Me%Monthly%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%Monthly%SquareGeomAverage2D (i, j) =                 &
                        (Me%Monthly%SquareGeomAverage2D (i, j) *            &
                        Me%Monthly%RunPeriod      +                         &
                        AuxValue**2 * DT) / (Me%Monthly%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%Monthly%SquareGeomAverage2D (i, j) -            &
                         Me%Monthly%GeomAverage2D       (i, j) ** 2

                    DX = abs(DX) 

                    Me%Monthly%GeomStandardDeviation2D(i, j) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical .and. (Value2D (i, j) > FillValueReal / 2.)) then
                    if (Value2D (i, j) < Me%CriticalValue) then
                        Me%Monthly%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%Monthly%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%Monthly%RunPeriod + DT) / (Me%Monthly%RunPeriod + DT)                        
                    else
                        Me%Monthly%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%Monthly%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%Monthly%RunPeriod + 0.) / (Me%Monthly%RunPeriod + DT)                        
                    endif
                endif

            endif
        enddo
        enddo


        !Verifies if the present time is a new output
        call ExtractDate (Time1=Me%ExternalVar%Now,         Month = PresentMonth)
        call ExtractDate (Time1=Me%Monthly%LastCalculation, Month = OldMonth)
        if (int(PresentMonth) /= int(OldMonth)) then
            call WriteValuesToFileHDF5 (.false., .false., .true., .false., .false.)
            Me%Monthly%Minimum2D           = Value2D
            Me%Monthly%Maximum2D           = Value2D
            Me%Monthly%Average2D           = Value2D
            Me%Monthly%SquareAverage2D     = Value2D
            Me%Monthly%StandardDeviation2D = Value2D

            if (Me%Accumulated) Me%Monthly%Accumulated2D = 0.0
            
            if (Me%GeomMean) then !Geometric Average to calculate
                Me%Monthly%GeomAverage2D           = Value2D
                Me%Monthly%SquareGeomAverage2D     = Value2D
                Me%Monthly%GeomStandardDeviation2D = Value2D
            endif

            !guillaume juan
            if (Me%Critical) then
                do j = JLB, JUB
                do i = ILB, IUB
                    if ( (WaterPoints2D(i, j) == WaterPoint) .and. (Value2D (i, j) > FillValueReal / 2.)) then
                        if(Value2D(i, j) < Me%CriticalValue) then
                            Me%Monthly%PercentBelowCriticalValue2D (i, j) = 1.
                        else
                            Me%Monthly%PercentBelowCriticalValue2D (i, j) = 0.
                        endif
                    endif
                enddo
                enddo
            endif

            Me%Monthly%RunPeriod  = 0.
        else
            Me%Monthly%RunPeriod  = Me%Monthly%RunPeriod + DT
        endif

        !Updates Time
        Me%Monthly%LastCalculation  = Me%ExternalVar%Now


        endif cd1


    end subroutine ModifyMonthlyStatistic2D

    !--------------------------------------------------------------------------

    subroutine ModifySpecificHourStatistic2D (Value2D, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:),   pointer          :: Value2D
        integer, dimension(:,:),   pointer          :: WaterPoints2D

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        real                                        :: DT, DX, AuxValue
        real                                        :: PresentHour

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

        call ExtractDate (Time1=Me%ExternalVar%Now,         Hour = PresentHour)
        
if1:    if (int(PresentHour) == int(Me%SpecificHourValue)) then
       
        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%SpecificHour%LastCalculation

cd1:    if (DT>0) then

        !Loops
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint) then

                !Minimum Value2D
                if (Value2D (i, j) < Me%SpecificHour%Minimum2D (i, j))            &
                    Me%SpecificHour%Minimum2D (i, j) = Value2D (i, j)

                !Maximum Value2D
                if (Value2D (i, j) > Me%SpecificHour%Maximum2D (i, j))            &
                    Me%SpecificHour%Maximum2D (i, j) = Value2D (i, j)

                !Average
                Me%SpecificHour%Average2D (i, j) =                                &
                    (Me%SpecificHour%Average2D (i, j) *                           &
                     Me%SpecificHour%RunPeriod      +                             &
                     Value2D (i, j) * DT) / (Me%SpecificHour%RunPeriod + DT)

                !Square Average
                Me%SpecificHour%SquareAverage2D (i, j) =                          &
                    (Me%SpecificHour%SquareAverage2D (i, j) *                     &
                     Me%SpecificHour%RunPeriod      +                             &
                     Value2D (i, j)**2 * DT) / (Me%SpecificHour%RunPeriod + DT)

                !Standard deviation
                DX = Me%SpecificHour%SquareAverage2D (i, j) -                     &
                     Me%SpecificHour%Average2D       (i, j) ** 2

                DX = abs(DX) 

                Me%SpecificHour%StandardDeviation2D(i, j) = sqrt(DX)

                !Accumulated Values
                if (Me%Accumulated)                                                             &
                    Me%SpecificHour%Accumulated2D (i,j) = Me%SpecificHour%Accumulated2D (i,j)   &
                                                        + Value2D (i, j) 

                if (Me%GeomMean) then !Geometric Average to calculate
                    !Geometric Average
                    AuxValue = Value2D (i, j)
                    if (AuxValue == 0.0) then
                        AuxValue = 1.0                   
                    elseif (AuxValue .lt. 0.0) then                   
                        write(*,*) 'Negative valued property.'
                        write(*,*) 'Geometric Average cannot be calculated.'
                        stop 'ModifySpecificHourStatistic2D - ModuleStatistic - ERR01'
                    endif
                
                    Me%SpecificHour%GeomAverage2D (i, j) = 10**                   &
                        ((LOG10(Me%SpecificHour%GeomAverage2D(i, j)) *            &
                        Me%SpecificHour%RunPeriod      +                          &
                        LOG10(AuxValue) * DT) / (Me%SpecificHour%RunPeriod + DT))

                    !Squared Geometric Average
                    Me%SpecificHour%SquareGeomAverage2D (i, j) =                  &
                        (Me%SpecificHour%SquareGeomAverage2D (i, j) *             &
                        Me%SpecificHour%RunPeriod      +                          &
                        AuxValue**2 * DT) / (Me%SpecificHour%RunPeriod + DT)

                    !Geometric Standard Deviation
                    DX = Me%SpecificHour%SquareGeomAverage2D (i, j) -             &
                         Me%SpecificHour%GeomAverage2D       (i, j) ** 2

                    DX = abs(DX) 

                    Me%SpecificHour%GeomStandardDeviation2D(i, j) = sqrt(DX)

                endif

                !guillaume juan
                if (Me%Critical .and. (Value2D (i, j) > FillValueReal / 2.)) then
                    if (Value2D (i, j) < Me%CriticalValue) then
                        Me%SpecificHour%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%SpecificHour%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%SpecificHour%RunPeriod + DT) / (Me%SpecificHour%RunPeriod + DT)                        
                    else
                        Me%SpecificHour%PercentBelowCriticalValue2D (i, j) =                               &
                            (Me%SpecificHour%PercentBelowCriticalValue2D (i, j) *                          &
                             Me%SpecificHour%RunPeriod + 0.) / (Me%SpecificHour%RunPeriod + DT)                        
                    endif
                endif

            endif
        enddo
        enddo


        !Updates Time
        Me%SpecificHour%RunPeriod     = Me%SpecificHour%RunPeriod + DT
        Me%SpecificHour%LastCalculation = Me%ExternalVar%Now

      
        endif cd1
        endif if1

    end subroutine ModifySpecificHourStatistic2D
    !--------------------------------------------------------------------------

    subroutine ModifyClassification2D (Value2D, WaterPoints2D)

        !Arguments-------------------------------------------------------------
        real,    dimension(:,:),   pointer          :: Value2D
        integer, dimension(:,:),   pointer          :: WaterPoints2D

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: iClass    
        real                                        :: DT, Aux

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB


        !Time since last calculation 
        DT  = Me%ExternalVar%Now - Me%Classification%LastCalculation

cd1:    if (DT>0) then


        !Loops
        do j = JLB, JUB
        do i = ILB, IUB
            if (WaterPoints2D(i, j) == WaterPoint) then
doClass:        do iClass = 1, Me%Classification%nClasses
                    if (Value2D(i, j) >= Me%Classification%Classes(iClass, 1) .and. &
                        Value2D(i, j)  < Me%Classification%Classes(iClass, 2)) then
                        Aux = DT
                    else
                        Aux = 0
                    endif

                    Me%Classification%Frequency2D     (i, j, iClass) =                   &
                        (Me%Classification%Frequency2D(i, j, iClass) *                   &
                         Me%Classification%RunPeriod + Aux)          /                   & 
                        (Me%Classification%RunPeriod + DT )

                enddo doClass
            endif
        enddo
        enddo

        Me%Classification%RunPeriod     = Me%Classification%RunPeriod + DT

        Me%Classification%LastCalculation = Me%ExternalVar%Now

        endif cd1


    end subroutine ModifyClassification2D

    !--------------------------------------------------------------------------
    
    subroutine WriteValuesToFileHDF5 (WriteGlobal, WriteDaily,             &
                                      WriteMonthly, WriteClassification,   &
                                      WriteSpecificHour)

        !Arguments-------------------------------------------------------------
        logical                                     :: WriteGlobal
        logical                                     :: WriteDaily
        logical                                     :: WriteMonthly
        logical                                     :: WriteClassification
        logical                                     :: WriteSpecificHour

        !Local-----------------------------------------------------------------
        integer                                     :: ILB, IUB, i
        integer                                     :: JLB, JUB, j
        integer                                     :: KLB, KUB, k
        integer                                     :: iClass 
        integer                                     :: STAT_CALL   
        character(len=StringLength)                 :: AuxChar, AuxChar1, AuxChar2
        real, dimension(:, :, :), pointer           :: AuxMatrix3D
        real, dimension(:, :   ), pointer           :: AuxMatrix2D
        real                                        :: Aux, d1, d2
        real(8), dimension(:), allocatable          :: P,C
        real(8)                                     :: Px,Cx, sumFreq
        integer                                     :: nc, n

        !Shorten
        ILB = Me%ExternalVar%WorkSize%ILB
        IUB = Me%ExternalVar%WorkSize%IUB
        JLB = Me%ExternalVar%WorkSize%JLB
        JUB = Me%ExternalVar%WorkSize%JUB

        if (Me%Methodology==Value3DStat3D_) then

            KLB = Me%ExternalVar%WorkSize%KLB
            KUB = Me%ExternalVar%WorkSize%KUB

        else if (Me%Methodology==Value3DStatLayers_) then

            KLB = 1
            KUB = Me%Layers%Number

        endif


        call HDF5SetLimits (Me%ObjHDF5, ILB, IUB, JLB, JUB, KLB, KUB,                                       &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR01'


        if (Me%Global%On .and. WriteGlobal) then

            !Minimum
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Minimum", &
                                      "Minimum","-", Array3D = Me%Global%Minimum,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR02'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Minimum", &
                                      "Minimum","-", Array2D = Me%Global%Minimum2D,                         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR03'

            endif

            !Maximum
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Maximum", &
                                      "Maximum","-", Array3D = Me%Global%Maximum,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR04'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Maximum", &
                                      "Maximum","-", Array2D = Me%Global%Maximum2D,                         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR05'

            endif

            !Average
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Average", &
                                      "Average", "-", Array3D = Me%Global%Average,                          &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR06'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Average", &
                                      "Average","-", Array2D = Me%Global%Average2D,                         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR07'

            endif

            !Standard Deviation
            if (Me%Methodology==Value3DStat3D_ .or.  &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/StandDev",&
                                      "StandDev", "-", Array3D = Me%Global%StandardDeviation,               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR08'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/StandDev",&
                                      "StandDev","-", Array2D = Me%Global%StandardDeviation2D,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR09'

            endif

           !Accumulated
           if (Me%Accumulated) then
           
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
            
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Accumulated", &
                                            "Accumulated", "-", Array3D = Me%Global%Accumulated,                    &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR10'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/Accumulated", &
                                            "Accumulated","-", Array2D = Me%Global%Accumulated2D,                   &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR11'

                endif           
           
            endif

            if (Me%GeomMean) then

                !Geometric Average
                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/GeomAverage", &
                                          "GeomAverage",                                                    & 
                                          "-", Array3D = Me%Global%GeomAverage,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR10'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/GeomAverage", & 
                                          "GeomAverage",                                                    &
                                          "-", Array2D = Me%Global%GeomAverage2D,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR11'

                endif

                !Geometric Standard Deviation
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/GeomStandDev",& 
                                          "GeomStandDev",                                                   &
                                          "-", Array3D = Me%Global%GeomStandardDeviation,                   &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR12'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/GeomStandDev",& 
                                          "GeomStandDev",                                                   &
                                          "-", Array2D = Me%Global%GeomStandardDeviation2D,                 &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR13'

                endif

            endif

            !guillaume e juan
            if (Me%Critical) then

                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/CriticalValue", &
                                          "CriticalValue",                                                    & 
                                          "-", Array3D = Me%Global%PercentBelowCriticalValue,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR10a'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/CriticalValue", & 
                                          "CriticalValue",                                                    &
                                          "-", Array2D = Me%Global%PercentBelowCriticalValue2D,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR11a'

                endif

            endif

            if (Me%Methodology==Value3DStatLayers_) then
    
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Global"//"/OpenPoints",& 
                                      "OpenPoints",                                                    &
                                      "-", Array3D = Me%Layers%WaterPoints,                            &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR13a'

            endif

            Me%Global%OutputNumber = Me%Global%OutputNumber + 1

        endif


        if (Me%Daily%On .and. WriteDaily) then

            !Minimum
            if (Me%Methodology==Value3DStat3D_ .or.                                                         &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Minimum", "Minimum",   &
                                      "-", Array3D = Me%Daily%Minimum,                                      &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR14'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Minimum", "Minimum",   &
                                      "-", Array2D = Me%Daily%Minimum2D,                                    &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR15'

            endif

            !Maximum
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Maximum", "Maximum",   &
                                      "-", Array3D = Me%Daily%Maximum,                                      &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR16'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Maximum", "Maximum",   &
                                      "-", Array2D = Me%Daily%Maximum2D,                                    &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR17'

            endif



            !Average
            if (Me%Methodology==Value3DStat3D_ .or.                                                         &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Average", "Average",   &
                                      "-", Array3D = Me%Daily%Average,                                      &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR18'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Average", "Average",   &
                                      "-", Array2D = Me%Daily%Average2D,                                    &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR19'

            endif

            !Standard Deviation
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/StandDev", "StandDev",  &
                                      "-", Array3D = Me%Daily%StandardDeviation,                            &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR20'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/StandDev", "StandDev",  &
                                      "-", Array2D = Me%Daily%StandardDeviation2D,                          &
                                      OutputNumber = Me%Daily%OutputNumber,                                 &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR21'

            endif

            !Accumulated
            if (Me%Accumulated) then
                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
            
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Accumulated",  &
                                          "Accumulated",  "-", Array3D = Me%Daily%Accumulated,                      &
                                          OutputNumber = Me%Daily%OutputNumber,                                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR20'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/Accumulated",  &
                                          "Accumulated", "-", Array2D = Me%Daily%Accumulated2D,                     &
                                          OutputNumber = Me%Daily%OutputNumber,                                     &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR21'

                endif
            endif

            if (Me%GeomMean) then

                !Geometric Average
                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/GeomAverage",          &
                                          "GeomAverage",                                                    & 
                                          "-", Array3D = Me%Daily%GeomAverage,                              &
                                          OutputNumber = Me%Daily%OutputNumber,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR22'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/GeomAverage",          & 
                                          "GeomAverage",                                                    &
                                          "-", Array2D = Me%Daily%GeomAverage2D,                            &
                                          OutputNumber = Me%Daily%OutputNumber,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR23'

                endif

                !Geometric Standard Deviation
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/GeomStandDev",          & 
                                          "GeomStandDev",                                                   &
                                          "-", Array3D = Me%Daily%GeomStandardDeviation,                    &
                                          OutputNumber = Me%Daily%OutputNumber,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR24'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/GeomStandDev",          & 
                                          "GeomStandDev",                                                   &
                                          "-", Array2D = Me%Daily%GeomStandardDeviation2D,                  &
                                          OutputNumber = Me%Daily%OutputNumber,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR25'

                endif

            endif


            !guillaume e juan
            if (Me%Critical) then

                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/CriticalValue", &
                                          "CriticalValue",                                                    & 
                                          "-", Array3D = Me%Daily%PercentBelowCriticalValue,                             &
                                          OutputNumber = Me%Daily%OutputNumber,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR10a'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/CriticalValue", & 
                                          "CriticalValue",                                                    &
                                          "-", Array2D = Me%Daily%PercentBelowCriticalValue2D,                           &
                                          OutputNumber = Me%Daily%OutputNumber,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR11a'

                endif

            endif

            if (Me%Methodology==Value3DStatLayers_) then
    
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Daily"//"/OpenPoints",& 
                                      "OpenPoints",                                                    &
                                      "-", Array3D = Me%Layers%WaterPoints,                            &
                                      OutputNumber = Me%Daily%OutputNumber,                            &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR25a'

            endif

            Me%Daily%OutputNumber = Me%Daily%OutputNumber + 1

        endif

        if (Me%Monthly%On .and. WriteMonthly) then

            !Minimum
            if (Me%Methodology==Value3DStat3D_ .or.  &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Minimum", "Minimum", &
                                      "-", Array3D = Me%Monthly%Minimum,                                    &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR26'

            else if (Me%Methodology==Value2DStat2D_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Minimum", "Minimum", &
                                      "-", Array2D = Me%Monthly%Minimum2D,                                  &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR27'

            endif

            !Maximum
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Maximum", "Maximum", &
                                      "-", Array3D = Me%Monthly%Maximum,                                    &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR28'

            else if (Me%Methodology==Value2DStat2D_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Maximum", "Maximum", &
                                      "-", Array2D = Me%Monthly%Maximum2D,                                  &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR29'

            endif

            !Average
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Average", "Average", &
                                      "-", Array3D = Me%Monthly%Average,                                    &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR30'

            else if (Me%Methodology==Value2DStat2D_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Average", "Average", &
                                      "-", Array2D = Me%Monthly%Average2D,                                  &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR31'

            endif


            !Standard Deviation
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/StandDev", "StandDev", &
                                      "-", Array3D = Me%Monthly%StandardDeviation,                          &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR32'

            else if (Me%Methodology==Value2DStat2D_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/StandDev", "StandDev", &
                                      "-", Array2D = Me%Monthly%StandardDeviation2D,                        &
                                      OutputNumber = Me%Monthly%OutputNumber,                               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR33'

            endif

           !Accumulated
           if (Me%Accumulated) then
           
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
            
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Accumulated",    &
                                            "Accumulated", "-", Array3D = Me%Monthly%Accumulated,                       &
                                            OutputNumber = Me%Monthly%OutputNumber,                                     &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR34'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/Accumulated",    &
                                            "Accumulated","-", Array2D = Me%Monthly%Accumulated2D,                      &
                                            OutputNumber = Me%Monthly%OutputNumber,                                     &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR35'

                endif           
           
           endif

            if (Me%GeomMean) then

                !Geometric Average
                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/GeomAverage",        &
                                          "GeomAverage",                                                    & 
                                          "-", Array3D = Me%Monthly%GeomAverage,                            &
                                          OutputNumber = Me%Monthly%OutputNumber,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR34'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/GeomAverage",        & 
                                          "GeomAverage",                                                    &
                                          "-", Array2D = Me%Monthly%GeomAverage2D,                          &
                                          OutputNumber = Me%Monthly%OutputNumber,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR35'

                endif

                !Geometric Standard Deviation
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/GeomStandDev",        & 
                                          "GeomStandDev",                                                   &
                                          "-", Array3D = Me%Monthly%GeomStandardDeviation,                  &
                                          OutputNumber = Me%Monthly%OutputNumber,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR36'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/GeomStandDev",        & 
                                          "GeomStandDev",                                                   &
                                          "-", Array2D = Me%Monthly%GeomStandardDeviation2D,                &
                                          OutputNumber = Me%Monthly%OutputNumber,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR37'

                endif


            endif

            !guillaume e juan
            if (Me%Critical) then

                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/CriticalValue", &
                                          "CriticalValue",                                                    & 
                                          "-", Array3D = Me%Monthly%PercentBelowCriticalValue,                             &
                                          OutputNumber = Me%Monthly%OutputNumber,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR10a'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/CriticalValue", & 
                                          "CriticalValue",                                                    &
                                          "-", Array2D = Me%Monthly%PercentBelowCriticalValue2D,                           &
                                          OutputNumber = Me%Monthly%OutputNumber,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR11a'

                endif

            endif

            if (Me%Methodology==Value3DStatLayers_) then
    
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Monthly"//"/OpenPoints",& 
                                      "OpenPoints",                                                    &
                                      "-", Array3D = Me%Layers%WaterPoints,                            &
                                      OutputNumber = Me%Daily%OutputNumber,                            &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR37a'

            endif


            Me%Monthly%OutputNumber = Me%Monthly%OutputNumber + 1

        endif

        if (Me%Classification%On .and. WriteClassification) then

            do iClass = 1, Me%Classification%nClasses


                write(AuxChar1, fmt=*)Me%Classification%Classes(iClass, 1)
                write(AuxChar2, fmt=*)Me%Classification%Classes(iClass, 2)

                AuxChar = trim(adjustl(AuxChar1))//"_"//trim(adjustl(AuxChar2))

            
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then

                    !Allocates auxiliar matrix 3D
                    allocate (AuxMatrix3D(                               &
                        Me%ExternalVar%Size%ILB:Me%ExternalVar%Size%IUB, &
                        Me%ExternalVar%Size%JLB:Me%ExternalVar%Size%JUB, &
                        Me%ExternalVar%Size%KLB:Me%ExternalVar%Size%KUB))


                            
                    do k = KLB, KUB
                    do j = JLB, JUB
                    do i = ILB, IUB
                        AuxMatrix3D(i, j, k) = 100. * Me%Classification%Frequency(i, j, k, iClass) 
                    enddo
                    enddo
                    enddo

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Classes",        &
                                          trim(adjustl(AuxChar)),                                           &
                                          "-", Array3D = AuxMatrix3D,                                       &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR38'

                else if (Me%Methodology==Value2DStat2D_) then

                    !Allocates auxiliar matrix 2D
                    allocate (AuxMatrix2D(Me%ExternalVar%Size%ILB:Me%ExternalVar%Size%IUB,                  &
                                          Me%ExternalVar%Size%JLB:Me%ExternalVar%Size%JUB))

                        

                    do j = JLB, JUB
                    do i = ILB, IUB
                        AuxMatrix2D(i, j) = 100. * Me%Classification%Frequency2D(i, j, iClass) 
                    enddo
                    enddo

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Classes",        &
                                          trim(adjustl(AuxChar)),                                           &
                                          "-", Array2D = AuxMatrix2D,                                       &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR39'

                endif

            enddo

            write(AuxChar1, fmt=*)Me%Classification%Percentil
            AuxChar = "Percentil_"//trim(adjustl(AuxChar1))//"_Class"
            
            if (Me%Methodology==Value3DStat3D_ .or.                                                         &
                Me%Methodology==Value3DStatLayers_) then
            
                !By default all domain belongs to the first class
                AuxMatrix3D(:,:,:) = 0.
                
                nc = Me%Classification%nClasses
                
                allocate(C(nc+1),P(nc+1))
                
                C(1) = Me%Classification%Classes(1, 1)
                
                P(1) = 0. 
                
                do iClass = 1, nc 
                
                    C(iClass+1) = Me%Classification%Classes(iClass, 2)
                    
                enddo
                
                Px = Me%Classification%Percentil                
                
                do k = KLB, KUB
                do j = JLB, JUB
                do i = ILB, IUB
                
                if (Me%NormalizeFreq) then
                
                    sumFreq = 0.

                    do iClass = 1, nc

                        sumFreq = sumFreq + Me%Classification%Frequency(i, j, k, iClass)

                    enddo
                    
                    if (sumFreq > 0.) then

                        do iClass = 1, nc

                            P(iClass + 1) = P(iClass) + Me%Classification%Frequency(i, j, k, iClass) / sumFreq * 100.

                        enddo

                    else                

                        do iClass = 1, nc

                            P(iClass + 1) = P(iClass) + Me%Classification%Frequency(i, j, k, iClass) * 100.

                        enddo
                        
                    endif
                    
                else
                
doClass1:           do iClass = 1, nc

                        P(iClass + 1) = P(iClass) + Me%Classification%Frequency(i, j, k, iClass) * 100.

                    enddo doClass1                
                    
                endif                    
                
                
                !Land point 
                if (P(nc + 1) == 0.) then
                    Cx = FillValueReal                
                    
                else if (Px > P(nc + 1)) then
                    Cx = C(nc + 1)
                    write(*,*) 'WriteValuesToFileHDF5 - ModuleStatistic - WRN35'
                    write(*,*) 'Percentil out of range please add more classes' 
                    write(*,*) 'Cell i=',i, ' j', j, ' k', k

                else                
                    ! polynomail interpolation 
                    !call polint(P(1:nc+1),C(1:nc+1),nc+1,Px,Cx,Error, STAT = status)
                    
                    !if large incertanty or a error is return - linear interpolation 
                    !if (Error > abs(C(nc+1) - C(1))/100. .or. status /= SUCCESS_) then
                        do n = 1,nc+1
                            if (Px>=P(n) .and. Px<=P(n+1)) then
                                d1 = Px - P(n)
                                d2 = P(n+1) - Px
                                Cx = (C(n) * d2 + C(n+1) * d1) / (d2 + d1)
                                exit
                            endif 
                        enddo
                    !endif
                    
                endif
                
                AuxMatrix3D(i, j, k) = Cx 
                
                enddo
                enddo
                enddo

                deallocate(C,P)

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Classes",            &
                                      trim(adjustl(AuxChar)),                                               &
                                      "-", Array3D = AuxMatrix3D,                                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR40'

                if (Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Classes"//"/OpenPoints",& 
                                          "OpenPoints",                                                    &
                                          "-", Array3D = Me%Layers%WaterPoints,                            &
                                          OutputNumber = Me%Daily%OutputNumber,                            &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR407a'

                endif

                !Deallocates auxiliar matrix
                deallocate (AuxMatrix3D)

            else if (Me%Methodology==Value2DStat2D_) then


                !By default all domain belongs to the first class
                AuxMatrix2D(:,:) = 0.
                
                nc = Me%Classification%nClasses
                
                allocate(C(nc+1),P(nc+1))
                
                C(1) = Me%Classification%Classes(1, 1)
                
                P(1) = 0. 
                
                do iClass = 1, nc 
                
                    C(iClass+1) = Me%Classification%Classes(iClass, 2)
                    
                enddo
                
                Px = Me%Classification%Percentil  

                do j = JLB, JUB
                do i = ILB, IUB
                Aux = 0.
                
!doClass2:       do iClass = 1, Me%Classification%nClasses

!                    Aux = Aux + Me%Classification%Frequency2D(i, j, iClass) / TotalFreq * 100.

!                    if (Aux <= Me%Classification%Percentil) then
!                        AuxMatrix2D(i, j) = float(iClass)
!                    else
!                        exit doClass2
!                    endif

                if (Me%NormalizeFreq) then
                
                    sumFreq = 0.

                    do iClass = 1, nc

                        sumFreq = sumFreq + Me%Classification%Frequency2D(i, j, iClass)

                    enddo
                    
                    if (sumFreq > 0.) then

                        do iClass = 1, nc

                            P(iClass + 1) = P(iClass) + Me%Classification%Frequency2D(i, j, iClass) / sumFreq * 100.

                        enddo

                    else

                        do iClass = 1, nc

                            P(iClass + 1) = P(iClass) + Me%Classification%Frequency2D(i, j, iClass) * 100.

                        enddo
                        
                    endif                
                else
                    
doClass2:           do iClass = 1, nc

                        P(iClass + 1) = P(iClass) + Me%Classification%Frequency2D(i, j, iClass) * 100.

                    enddo doClass2
                
                endif
                
                
                !Land point 
                if (P(nc + 1) == 0.) then
                    Cx = FillValueReal                
                    
                else if (Px > P(nc + 1)) then
                    Cx = C(nc + 1)
                    write(*,*) 'WriteValuesToFileHDF5 - ModuleStatistic - WRN35'
                    write(*,*) 'Percentil out of range please add more classes' 
                    write(*,*) 'Cell i=',i, ' j', j

                else                
                    ! polynomail interpolation 
                    !call polint(P(1:nc+1),C(1:nc+1),nc+1,Px,Cx,Error, STAT = status)
                    
                    !if large incertanty or a error is return - linear interpolation 
                    !if (Error > abs(C(nc+1) - C(1))/100. .or. status /= SUCCESS_) then
                        do n = 1,nc+1
                            if (Px>=P(n) .and. Px<=P(n+1)) then
                                d1 = Px - P(n)
                                d2 = P(n+1) - Px
                                Cx = (C(n) * d2 + C(n+1) * d1) / (d2 + d1)
                                exit
                            endif 
                        enddo
                    !endif
                    
                endif


                AuxMatrix2D(i, j) = Cx                     
                
                enddo
                enddo

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/Classes",            &
                                      trim(adjustl(AuxChar)),                                               &
                                      "-", Array2D = AuxMatrix2D,                                           &
                                      STAT = STAT_CALL)

                !Deallocates auxiliar matrix
                deallocate (AuxMatrix2D)
                deallocate (C,P)

            endif

        endif

        if (Me%SpecificHour%On .and. WriteSpecificHour) then

            !Minimum
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Minimum", &
                                      "Minimum","-", Array3D = Me%SpecificHour%Minimum,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR51'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Minimum", &
                                      "Minimum","-", Array2D = Me%SpecificHour%Minimum2D,                         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR52'

            endif

            !Maximum
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Maximum", &
                                      "Maximum","-", Array3D = Me%SpecificHour%Maximum,                           &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR53'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Maximum", &
                                      "Maximum","-", Array2D = Me%SpecificHour%Maximum2D,                         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR54'

            endif

            !Average
            if (Me%Methodology==Value3DStat3D_ .or. &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Average", &
                                      "Average", "-", Array3D = Me%SpecificHour%Average,                          &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR55'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Average", &
                                      "Average","-", Array2D = Me%SpecificHour%Average2D,                         &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR56'

            endif

            !Standard Deviation
            if (Me%Methodology==Value3DStat3D_ .or.  &
                Me%Methodology==Value3DStatLayers_) then
            
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/StandDev",&
                                      "StandDev", "-", Array3D = Me%SpecificHour%StandardDeviation,               &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR57'

            else if (Me%Methodology==Value2DStat2D_) then

                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/StandDev",&
                                      "StandDev","-", Array2D = Me%SpecificHour%StandardDeviation2D,              &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR58'

            endif

           !Accumulated
           if (Me%Accumulated) then
           
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
            
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Accumulated",   &
                                            "Accumulated", "-", Array3D = Me%SpecificHour%Accumulated,                      &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR34'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/Accumulated",   &
                                            "Accumulated","-", Array2D = Me%SpecificHour%Accumulated2D,                     &
                                            STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR35'

                endif           
           
           endif


            if (Me%GeomMean) then

                !Geometric Average
                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/GeomAverage", &
                                          "GeomAverage",                                                    & 
                                          "-", Array3D = Me%SpecificHour%GeomAverage,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR60'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/GeomAverage", & 
                                          "GeomAverage",                                                    &
                                          "-", Array2D = Me%SpecificHour%GeomAverage2D,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR61'

                endif

                !Geometric Standard Deviation
                if (Me%Methodology==Value3DStat3D_ .or.  &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/GeomStandDev",& 
                                          "GeomStandDev",                                                   &
                                          "-", Array3D = Me%SpecificHour%GeomStandardDeviation,                   &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR62'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/GeomStandDev",& 
                                          "GeomStandDev",                                                   &
                                          "-", Array2D = Me%SpecificHour%GeomStandardDeviation2D,                 &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR63'

                endif

            endif

            !guillaume e juan
            if (Me%Critical) then

                if (Me%Methodology==Value3DStat3D_ .or. &
                    Me%Methodology==Value3DStatLayers_) then
        
                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/CriticalValue", &
                                          "CriticalValue",                                                    & 
                                          "-", Array3D = Me%SpecificHour%PercentBelowCriticalValue,                             &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR10a'

                else if (Me%Methodology==Value2DStat2D_) then

                    call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/CriticalValue", & 
                                          "CriticalValue",                                                    &
                                          "-", Array2D = Me%SpecificHour%PercentBelowCriticalValue2D,                           &
                                          STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR11a'

                endif

            endif

            if (Me%Methodology==Value3DStatLayers_) then
    
                call HDF5WriteData   (Me%ObjHDF5, trim(Me%GroupName)//trim(Me%Name)//"/SpecificHour"//"/OpenPoints",& 
                                      "OpenPoints",                                                    &
                                      "-", Array3D = Me%Layers%WaterPoints,                            &
                                      STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFileHDF5 - ModuleStatistic - ERR63a'

            endif

            Me%SpecificHour%OutputNumber = Me%SpecificHour%OutputNumber + 1

        endif


        if (STAT_CALL /= SUCCESS_) stop 'WriteValuesToFile - ModuleStatistic - ERR99'

    end subroutine WriteValuesToFileHDF5


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine GetStatisticParameters(StatisticID, Value3DStat3D, Value3DStatLayers,     &
                                      Value2DStat2D, Depth, Layer, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: StatisticID
        integer, optional, intent(OUT) :: Value3DStat3D, Value3DStatLayers,              &
                                          Value2DStat2D, Depth, Layer
        integer, optional, intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: ready_              
        integer                         :: STAT_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(Value3DStat3D    )) Value3DStat3D    = Value3DStat3D_
            if (present(Value3DStatLayers)) Value3DStatLayers= Value3DStatLayers_
            if (present(Value2DStat2D    )) Value2DStat2D    = Value2DStat2D_
            if (present(Depth            )) Depth            = Depth_
            if (present(Layer            )) Layer            = Layer_

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetStatisticParameters

    !--------------------------------------------------------------------------

    subroutine GetStatisticLayerDef(StatisticID, LayerNumber, LayerDefinition, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: StatisticID
        integer,           intent(IN)  :: LayerNumber
        integer,           intent(OUT) :: LayerDefinition
        integer, optional, intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: ready_              
        integer                         :: STAT_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            LayerDefinition  = Me%Layers%Definition(LayerNumber)
            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetStatisticLayerDef

    !--------------------------------------------------------------------------

    subroutine GetStatisticLayersNumber(StatisticID, LayersNumber, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: StatisticID
        integer,           intent(OUT) :: LayersNumber
        integer, optional, intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: ready_              
        integer                         :: STAT_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then
            LayersNumber  = Me%Layers%Number

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetStatisticLayersNumber

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine GetStatisticClasses(StatisticID, Classes, STAT)

        !Arguments-------------------------------------------------------------
        integer                                  :: StatisticID
        real, dimension(:,:),          pointer   :: Classes
        integer, optional,           intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                  :: ready_              
        integer                                  :: STAT_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mSTATISTIC_, Me%InstanceID)

            Classes  => Me%Classification%Classes

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetStatisticClasses

    !--------------------------------------------------------------------------

    subroutine GetStatisticClassesNumber(StatisticID, nClasses, STAT)

        !Arguments-------------------------------------------------------------
        integer                                  :: StatisticID
        integer,                     intent(OUT) :: nClasses
        integer, optional,           intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                                  :: ready_              
        integer                                  :: STAT_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            nClasses  = Me%Classification%nClasses

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetStatisticClassesNumber

    !--------------------------------------------------------------------------



    subroutine GetStatisticMethod(StatisticID, MethodStatistic, STAT)

        !Arguments-------------------------------------------------------------
        integer                        :: StatisticID
        integer,           intent(OUT) :: MethodStatistic
        integer, optional, intent(OUT) :: STAT

        !Local-----------------------------------------------------------------
        integer                         :: ready_              
        integer                         :: STAT_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_) 
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            MethodStatistic  = Me%Methodology

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine GetStatisticMethod

    !--------------------------------------------------------------------------

    !----------------------------------------------------------------------
    
    subroutine UnGetStatistic2D(StatisticID, Array2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: StatisticID
        real, pointer, dimension(:,:)               :: Array2D
        integer, optional, intent(OUT)              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: ready_          
        integer                                     :: STAT_ 

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(StatisticID, ready_)    

cd1 :   if (ready_ .EQ. READ_LOCK_ERR_) then
            nullify(Array2D)

            call Read_UnLock(mSTATISTIC_, Me%InstanceID, "UngetStatistic2D")

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))                                                    &
            STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine UnGetStatistic2D

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCT 

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine KillStatistic (StatisticID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: StatisticID
        integer, optional                           :: STAT

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_, ready_    
        integer                                     :: nUsers

        STAT_ = UNKNOWN_

        call Ready (StatisticID, ready_)
           
cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mSTATISTIC_,  Me%InstanceID)

            if (nUsers == 0) then
                
                !Writes the final values to the HDF file
                call WriteValuesToFileHDF5 (.true., .true., .true., .true., .true.)

                !Associates External Instances
                nUsers = DeassociateInstance (mTIME_,          Me%ObjTime         )
                if (nUsers == 0) stop 'KillStatistic - ModuleStatistic - ERR01'
                nUsers = DeassociateInstance (mHDF5_,          Me%ObjHDF5         )
                if (nUsers == 0) stop 'KillStatistic - ModuleStatistic - ERR02'

                if (Me%Global%On)                                                   &
                    call DeallocateStatisticMatrixes (Me%Global)

                if (Me%SpecificHour%On)                                             &
                    call DeallocateStatisticMatrixes (Me%SpecificHour)
                
                if (Me%Daily%On)                                                    & 
                    call DeallocateStatisticMatrixes (Me%Daily)
        
                if (Me%Monthly%On)                                                  &
                    call DeallocateStatisticMatrixes (Me%Monthly)

                if (Me%Classification%On)  then

                    deallocate (Me%Classification%Classes)

                    if (Me%Methodology==Value3DStat3D_ .or.                         &
                        Me%Methodology==Value3DStatLayers_) then

                        deallocate (Me%Classification%Frequency)

                    else if (Me%Methodology==Value2DStat2D_) then

                        deallocate (Me%Classification%Frequency2D)

                    endif


                    if (Me%Methodology == Value3DStatLayers_) then

                        call DeallocateLayerMatrixes

                    endif

                endif

                !Deallocates Instance of Statistic
                call DeallocateInstance 

                StatisticID = 0
                STAT_       = SUCCESS_

            endif

        else cd1

            STAT_ = ready_

        endif cd1

        if (present(STAT)) STAT = STAT_


    end subroutine KillStatistic

    !--------------------------------------------------------------------------

    subroutine DeallocateInstance

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_Statistic), pointer                 :: AuxStatistic
        type (T_Statistic), pointer                 :: PreviousStatistic

        !Updates pointers
        if (Me%InstanceID == FirstStatistic%InstanceID) then
            FirstStatistic => FirstStatistic%Next
        else
            PreviousStatistic => FirstStatistic
            AuxStatistic      => FirstStatistic%Next
            do while (AuxStatistic%InstanceID /= Me%InstanceID)
                PreviousStatistic => AuxStatistic
                AuxStatistic      => AuxStatistic%Next
            enddo

            !Now update linked list
            PreviousStatistic%Next => AuxStatistic%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------

    subroutine DeallocateStatisticMatrixes (Statistic)

        !Arguments-------------------------------------------------------------
        type (T_SimpleStatistic)                    :: Statistic

        !Local-----------------------------------------------------------------

        if (Me%Methodology==Value3DStat3D_ .or.                                &
            Me%Methodology==Value3DStatLayers_) then

            deallocate (Statistic%Minimum)
            deallocate (Statistic%Maximum)    
            deallocate (Statistic%Average)    
            deallocate (Statistic%SquareAverage) 
            deallocate (Statistic%StandardDeviation) 

            if (Me%Accumulated) then           
                deallocate (Statistic%Accumulated)  
            endif

            if (Me%GeomMean) then !Geometric Average is to be calculated
                deallocate (Statistic%GeomAverage)
                deallocate (Statistic%SquareGeomAverage)
                deallocate (Statistic%GeomStandardDeviation)
            endif

            !guillaume juan
            if (Me%Critical) then
                deallocate (Statistic%PercentBelowCriticalValue) 
            endif

        else if (Me%Methodology==Value2DStat2D_) then

            deallocate (Statistic%Minimum2D)
            deallocate (Statistic%Maximum2D)    
            deallocate (Statistic%Average2D)    
            deallocate (Statistic%SquareAverage2D) 
            deallocate (Statistic%StandardDeviation2D) 

            if (Me%Accumulated) then           
                deallocate (Statistic%Accumulated2D)  
            endif

            if (Me%GeomMean) then !Geometric Average is to be calculated
                deallocate (Statistic%GeomAverage2D)
                deallocate (Statistic%SquareGeomAverage2D)
                deallocate (Statistic%GeomStandardDeviation2D)
            endif

            !guillaume juan
            if (Me%Critical) then
                deallocate (Statistic%PercentBelowCriticalValue2D) 
            endif

        endif

    end subroutine DeallocateStatisticMatrixes


    !--------------------------------------------------------------------------

    subroutine DeallocateLayerMatrixes

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------

 
        deallocate (Me%Layers%Value)
        deallocate (Me%Layers%WaterPoints)


        deallocate (Me%Layers%UpperDepth)
        deallocate (Me%Layers%LowerDepth)

        deallocate (Me%Layers%UpperLayer)
        deallocate (Me%Layers%LowerLayer)

    end subroutine DeallocateLayerMatrixes



    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !--------------------------------------------------------------------------

    subroutine Ready (StatisticID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: StatisticID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (StatisticID > 0) then
            call LocateObjStatistic (StatisticID)
            ready_ = VerifyReadLock (mSTATISTIC_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjStatistic (StatisticID)

        !Arguments-------------------------------------------------------------
        integer                                     :: StatisticID

        !Local-----------------------------------------------------------------

        Me => FirstStatistic
        do while (associated (Me))
            if (Me%InstanceID == StatisticID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) then
            stop 'ModuleStatistic - LocateObjStatistic - ERR01'
        endif

    end subroutine LocateObjStatistic

    !--------------------------------------------------------------------------


End Module ModuleStatistic

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------

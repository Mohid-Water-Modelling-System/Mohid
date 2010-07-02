!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : JoinTimeSeries
! PROGRAM       : MainJoinTimeSeries
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : January 2005
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to extract moving time series from MOHID HDF5 files
!
!------------------------------------------------------------------------------

program MohidJoinTimeSeries

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleTimeSerie
    use ModuleFunctions

    implicit none

    !CPU Time variables
    type (T_Time)                                           :: InitialSystemTime, FinalSystemTime
    real                                                    :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                                   :: F95Time
    
    type       T_TimeSeriesFile
        integer                                             :: ObjTimeSerie         = 0
        character(len=PathLength)                           :: FileName
        character(len=line_length)                          :: Residual = null_str
        type(T_Time)                                        :: BeginTime
        character(len=256)                                  :: TimeUnits
        real, dimension(:,:), pointer                       :: DataMatrix
        integer                                             :: nValues, nColumns
        type(T_TimeSeriesFile), pointer                     :: Next
    end type  T_TimeSeriesFile


    type T_JoinTimeSeries
        character(PathLength)                               :: DataFile
        type(T_Size3D)                                      :: Size3D
        type(T_Time)                                        :: BeginTime
        type(T_Time)                                        :: EndTime
        type(T_Time)                                        :: CurrentTime
        type(T_TimeSeriesFile),  pointer                    :: FirstTimeSeriesFile
        integer                                             :: ObjEnterData         = 0
        integer                                             :: nTimeSeries          = 0
        character(len=PathLength)                           :: OutputFileName
        type(T_Time), dimension(:), pointer                 :: InitialTimes, OrganizedTimes
        real                                                :: DT, AddModelDT
        integer                                             :: LocationI, LocationJ, LocationK
        logical                                             :: WriteResiduals
        logical                                             :: FirstEqualsLast
    end type T_JoinTimeSeries

    type(T_JoinTimeSeries), pointer                         :: Me


    call ConstructMohidJoinTimeSeries
    call ModifyMohidJoinTimeSeries
    call KillMohidJoinTimeSeries

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidJoinTimeSeries

        !Begin-----------------------------------------------------------------

        allocate(Me)
        
        call StartUpMohid("MohidJoinTimeSeries")

        call StartCPUTime

        call ReadGlobalData

        call OrganizeTimeSeriesFiles

        call ReadLocation

    end subroutine ConstructMohidJoinTimeSeries
    
    !--------------------------------------------------------------------------


    subroutine ModifyMohidJoinTimeSeries
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j, OutputUnit
        type(T_TimeSeriesFile),  pointer            :: TimeSeriesFile
        real                                        :: SecondsSinceFirst
        character(len=4)                            :: CharFormat = '   '
        integer                                     :: n_columns
        character(len=line_length)                  :: Header

        !Begin-----------------------------------------------------------------

        n_columns = Me%FirstTimeSeriesFile%nColumns

        if    (n_columns .lt. 10                             )then
            write(CharFormat, '(i1)')n_columns
        elseif(n_columns .ge. 10   .and. n_columns .lt. 100  )then
            write(CharFormat, '(i2)')n_columns
        elseif(n_columns .ge. 100  .and. n_columns .lt. 1000 )then
            write(CharFormat, '(i3)')n_columns
        elseif(n_columns .ge. 1000 .and. n_columns .lt. 10000)then
            write(CharFormat, '(i4)')n_columns
        else
            stop 'Number of volumes limited to 9999.'
        endif

        call GetTimeSerieHeader(Me%FirstTimeSeriesFile%ObjTimeSerie, Header, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidJoinTimeSeries - MohidJoinTimeSeries - ERR00'
       
        call UnitsManager(OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidJoinTimeSeries - MohidJoinTimeSeries - ERR01'

        open(unit = OutputUnit, file = trim(Me%OutputFileName ), status = 'replace')


        call WriteDataLine(OutputUnit, "SERIE_INITIAL_DATA", Me%OrganizedTimes(1))
        call WriteDataLine(OutputUnit, "TIME_UNITS", Me%FirstTimeSeriesFile%TimeUnits)

        call WriteDataLine(OutputUnit, "LOCALIZATION_I", Me%LocationI)
        call WriteDataLine(OutputUnit, "LOCALIZATION_J", Me%LocationJ)
        call WriteDataLine(OutputUnit, "LOCALIZATION_K", Me%LocationK)


        call WriteDataLine(OutputUnit, Header)
        call WriteDataLine(OutputUnit, "<BeginTimeSerie>")

        SecondsSinceFirst = Me%AddModelDT 

        do i = 1, Me%nTimeSeries

            TimeSeriesFile => Me%FirstTimeSeriesFile

            do while(associated(TimeSeriesFile))

                if(TimeSeriesFile%BeginTime .eq. Me%OrganizedTimes(i))then

                    write(*,*)TimeToString(TimeSeriesFile%BeginTime)

                    do j = 1, TimeSeriesFile%nValues

                        if(j == 1 .and. i > 1 .and. Me%FirstEqualsLast)then

                        else

                            write(OutputUnit, '('//CharFormat//'f25.5, 1x)') SecondsSinceFirst, TimeSeriesFile%DataMatrix(j, 2:)

                            SecondsSinceFirst = SecondsSinceFirst + Me%DT

                        endif



                    enddo

                end if

                TimeSeriesFile => TimeSeriesFile%Next

            enddo

        enddo

        call WriteDataLine(OutputUnit, "<EndTimeSerie>")

        if(Me%WriteResiduals)then
            
            do i = 1, Me%nTimeSeries

                TimeSeriesFile => Me%FirstTimeSeriesFile

                do while(associated(TimeSeriesFile))

                    if(TimeSeriesFile%BeginTime .eq. Me%OrganizedTimes(i))then
                    
                        call WriteDataLine(OutputUnit, "")
                        call WriteDataLine(OutputUnit, "<BeginResidual>")
                        write(OutputUnit, '(A)') TimeSeriesFile%Residual
                        call WriteDataLine(OutputUnit, "<EndResidual>")
                    
                    end if

                    TimeSeriesFile => TimeSeriesFile%Next

                enddo

            enddo

        end if
        
    
    end subroutine ModifyMohidJoinTimeSeries

    !--------------------------------------------------------------------------
    

    subroutine ReadGlobalData

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        
        !Begin-----------------------------------------------------------------

        call ConstructEnterData (Me%ObjEnterData, 'JoinTimeSeries.dat', STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR01'
    
        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'OUT_FILE',                         &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR02'

            
        call GetData(Me%DT,                                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'DT',                               &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR03'

        call GetData(Me%AddModelDT,                                     &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'ADD_MODEL_DT',                     &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR04'
        
        call GetData(Me%WriteResiduals,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'WRITE_RESIDUALS',                  &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR05'

        call GetData(Me%FirstEqualsLast,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'FIRST_EQUALS_LAST',                &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR06'

        

        call ConstructTimeSeriesFilesList

        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidJoinTimeSeries - ERR07'

    end subroutine ReadGlobalData

    !--------------------------------------------------------------------------

    subroutine ConstructTimeSeriesFilesList

        !Arguments-------------------------------------------------------------
          
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, ClientNumber
        type (T_TimeSeriesFile), pointer        :: NewTimeSeriesFile
        logical                                 :: BlockFound
        logical                                 :: AtLeastOneBlock = .false.
        integer                                 :: iflag

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                            &
                                        ClientNumber    = ClientNumber,             &
                                        block_begin     = '<BeginTimeSeriesFile>',  &
                                        block_end       = '<EndTimeSeriesFile>',    &
                                        BlockFound      = BlockFound,               &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddTimeSeriesFile (NewTimeSeriesFile)

                    call GetData(NewTimeSeriesFile%FileName,                        &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'NAME',                             &
                                 ClientModule = 'MohidJoinTimeSeries',              &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'ConstructTimeSeriesFile - MohidJoinTimeSeries - ERR00'

                    if(Me%WriteResiduals)then
                        call ReadResiduals(NewTimeSeriesFile)
                    end if

                    nullify(NewTimeSeriesFile)

                else cd2
                    
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                    &
                        stop 'ConstructTimeSeriesFilesList - MohidJoinTimeSeries - ERR01'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ConstructTimeSeriesFilesList - MohidJoinTimeSeries - ERR02'
            else cd1
                stop 'ConstructTimeSeriesFilesList - MohidJoinTimeSeries - ERR03'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No TimeSeriesFile block is indicated in input file. '
            stop 'ConstructTimeSeriesFilesList - MohidJoinTimeSeries - ERR04'
        end if

    end subroutine ConstructTimeSeriesFilesList

    !--------------------------------------------------------------------------

    subroutine AddTimeSeriesFile(ObjTimeSeriesFile)

        !Arguments-------------------------------------------------------------
        type (T_TimeSeriesFile),     pointer           :: ObjTimeSeriesFile

        !Local-----------------------------------------------------------------
        type (T_TimeSeriesFile),     pointer           :: PreviousTimeSeriesFile
        type (T_TimeSeriesFile),     pointer           :: NewTimeSeriesFile

        !Begin-----------------------------------------------------------------

        !Allocates new TimeSeriesFile
        allocate (NewTimeSeriesFile)
        nullify  (NewTimeSeriesFile%Next)

        !Insert new Parameter into list and makes current ?? point to it
        if (.not. associated(Me%FirstTimeSeriesFile)) then
            Me%FirstTimeSeriesFile         => NewTimeSeriesFile
            ObjTimeSeriesFile              => NewTimeSeriesFile
        else
            PreviousTimeSeriesFile         => Me%FirstTimeSeriesFile
            ObjTimeSeriesFile              => Me%FirstTimeSeriesFile%Next
            do while (associated(ObjTimeSeriesFile))
                PreviousTimeSeriesFile     => ObjTimeSeriesFile
                ObjTimeSeriesFile          => ObjTimeSeriesFile%Next
            enddo
            ObjTimeSeriesFile              => NewTimeSeriesFile
            PreviousTimeSeriesFile%Next    => NewTimeSeriesFile
        end if

        Me%nTimeSeries = Me%nTimeSeries + 1

    end subroutine AddTimeSeriesFile
    
    !--------------------------------------------------------------------------

    subroutine OrganizeTimeSeriesFiles
       
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i, j, k
        type (T_TimeSeriesFile),    pointer         :: TimeSeriesFile
        logical                                     :: exist
        type(T_Time)                                :: OldestTime
        logical                                     :: Done 

        !Begin-----------------------------------------------------------------
        
        allocate(Me%InitialTimes  (1:Me%nTimeSeries))
        allocate(Me%OrganizedTimes(1:Me%nTimeSeries))

        TimeSeriesFile => Me%FirstTimeSeriesFile

        i = 1
        
        do while (associated(TimeSeriesFile))

            inquire(FILE = TimeSeriesFile%FileName, EXIST = exist)
            if (.not. exist) then
                write(*,*)'Time series file does not exist:'//trim(TimeSeriesFile%FileName)
                stop 'OrganizeTimeSeriesFiles - MohidJoinTimeSeries - ERR01'
            endif

            call StartTimeSerieInput(TimeSeriesFile%ObjTimeSerie,                       &
                                     TimeSeriesFile%FileName,                           &
                                     STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OrganizeTimeSeriesFiles - MohidJoinTimeSeries - ERR02'

            call GetTimeSerieInitialData (TimeSeriesFile%ObjTimeSerie, TimeSeriesFile%BeginTime, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_)stop 'OrganizeTimeSeriesFiles - MohidJoinTimeSeries - ERR03'

            call GetTimeSerieDataMatrix(TimeSeriesFile%ObjTimeSerie, TimeSeriesFile%DataMatrix, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidJoinTimeSeries - ERR03'


            call GetTimeSerieTimeUnits(TimeSeriesFile%ObjTimeSerie, TimeSeriesFile%TimeUnits, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidJoinTimeSeries - ERR03'
            
            call GetTimeSerieDataValues(TimeSeriesFile%ObjTimeSerie, TimeSeriesFile%nValues, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidJoinTimeSeries - ERR03'

            call GetTimeSerieDataColumns(TimeSeriesFile%ObjTimeSerie, TimeSeriesFile%nColumns, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ReadMovingProfileFile - MohidJoinTimeSeries - ERR03'

            Me%InitialTimes(i) = TimeSeriesFile%BeginTime

            i = i + 1 

            TimeSeriesFile => TimeSeriesFile%Next

        end do
        
        i = 1
        j = 1
        Done = .false.

        do while (.not. Done)

            call SetDate(OldestTime, 3000,1,1,0,0,0)
           
            do i = 1, Me%nTimeSeries
                if(Me%InitialTimes(i) .lt. OldestTime)then
                    OldestTime              = Me%InitialTimes(i)
                    Me%OrganizedTimes(j)    = Me%InitialTimes(i)
                    k                       = i
                endif
            end do

            call SetDate(Me%InitialTimes(k), 3000,1,1,0,0,0)

            if(j .eq. Me%nTimeSeries)then
                Done = .true.
            end if

            j = j + 1

        end do

        write(*,*)'organized files'

    end subroutine OrganizeTimeSeriesFiles

    !--------------------------------------------------------------------------

    subroutine ReadLocation
        
        !Local-----------------------------------------------------------------
        integer                                     :: ObjEnterData = 0
        integer                                     :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------


        call ConstructEnterData (ObjEnterData, Me%FirstTimeSeriesFile%FileName, &
                                 STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLocation - MohidJoinTimeSeries - ERR01'

        call GetData(Me%LocationI,                                      &
                     ObjEnterData, iflag,                               &
                     SearchType   = FromFile,                           &
                     keyword      = 'LOCALIZATION_I',                   &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     Default      = -9999,                              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLocation - MohidJoinTimeSeries - ERR02'


        call GetData(Me%LocationJ,                                      &
                     ObjEnterData, iflag,                               &
                     SearchType   = FromFile,                           &
                     keyword      = 'LOCALIZATION_J',                   &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     Default      = -9999,                              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLocation - MohidJoinTimeSeries - ERR03'

        call GetData(Me%LocationJ,                                      &
                     ObjEnterData, iflag,                               &
                     SearchType   = FromFile,                           &
                     keyword      = 'LOCALIZATION_K',                   &
                     ClientModule = 'MohidJoinTimeSeries',              &
                     Default      = -9999,                              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLocation - MohidJoinTimeSeries - ERR04'


        call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadLocation - MohidJoinTimeSeries - ERR05'


    end subroutine ReadLocation

    !--------------------------------------------------------------------------

    subroutine ReadResiduals(TimeSeriesFile)

        !Arguments-------------------------------------------------------------
        type (T_TimeSeriesFile), pointer            :: TimeSeriesFile

        !Local-----------------------------------------------------------------
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber, STAT_CALL
        integer                                     :: StartLine, EndLine
        logical                                     :: BlockFound
        !Begin-----------------------------------------------------------------


        call ConstructEnterData (ObjEnterData, TimeSeriesFile%FileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadResiduals - MohidJoinTimeSeries - ERR01'

        call ExtractBlockFromBuffer(ObjEnterData,                               &
                                    ClientNumber    = ClientNumber,             &
                                    block_begin     = '<BeginResidual>',        &
                                    block_end       = '<EndResidual>',          &
                                    BlockFound      = BlockFound,               &
                                    STAT            = STAT_CALL)
        if(STAT_CALL .EQ. SUCCESS_)then
            
            if (BlockFound) then
                
                !Gets the number of values in the block
                call GetBlockSize(ObjEnterData, ClientNumber, StartLine, EndLine)

                call GetFullBufferLine(ObjEnterData, StartLine+1, TimeSeriesFile%Residual, STAT = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadResiduals - MohidJoinTimeSeries - ERR02'
            
            else 
                    
                call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                if (STAT_CALL .NE. SUCCESS_) stop 'ReadResiduals - MohidJoinTimeSeries - ERR03'
            
            end if

        elseif (STAT_CALL .EQ. BLOCK_END_ERR_) then
                
            write(*,*)  
            write(*,*) 'Error calling ReadResiduals. '
            stop 'ReadResiduals - MohidJoinTimeSeries - ERR04'
        
        else
            
            stop 'ReadResiduals - MohidJoinTimeSeries - ERR05'

        end if


        call KillEnterData(ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadResiduals - MohidJoinTimeSeries - ERR06'


    end subroutine ReadResiduals

    !--------------------------------------------------------------------------


    subroutine KillMohidJoinTimeSeries

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type(T_TimeSeriesFile), pointer             :: TimeSeriesFile
        !Begin-----------------------------------------------------------------


        TimeSeriesFile => Me%FirstTimeSeriesFile
        
        do while (associated(TimeSeriesFile))

            call KillTimeSerie(TimeSeriesFile%ObjTimeSerie, STAT = STAT_CALL)
            if(STAT_CALL .ne. SUCCESS_)stop 'KillMohidJoinTimeSeries - ERR01'

            TimeSeriesFile => TimeSeriesFile%Next

        end do

        call StopCPUTime

        call ShutdownMohid ("MohidJoinTimeSeries", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidJoinTimeSeries
    
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


end program MohidJoinTimeSeries

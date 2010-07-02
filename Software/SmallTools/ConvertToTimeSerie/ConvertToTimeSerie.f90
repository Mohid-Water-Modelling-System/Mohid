!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToTimeSerie
! PROGRAM       : ConvertToTimeSerie
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : April 2005
! REVISION      : Angela Canas - v4.0
! DESCRIPTION   : ConvertToTimeSerie to convert files in IH and Toga tide 
!                 gauge format to TimeSerie format
!
!------------------------------------------------------------------------------

!Nomfich.dat:
!
!   IN_MODEL                : char                  [-]     !Path to the input file with
!                                                           !user's instructions (DataFile)                                                                
!   ROOT_SRT                : char                  [-]     !Path of folder where the
!                                                           !files to convert are
!                                                           !and where output files
!                                                           !will appear                   
!DataFile:
!
!   INPUT_FILE_TYPE         : char                  [IH]    !Type of file to convert
!                                                           !(options: IH (Inst. Hidrográfico)
!                                                           !          TOGA (Toga Center))
!   INPUT_FILENAME          : char                  [-]     !Path to the file to convert if IH
!
!   <BeginTogaFile>
!   INPUT_FILENAME          : char                  [-]     !Path to the file to convert if TOGA 
!   <EndTogaFile>                                                                                          
!
!   START                   : YYYY MM DD HH MM SS   [-]     !Start date of time series
!   END                     : YYYY MM DD HH MM SS   [-]     !End date of time series
!                                                           !(only day is relevant for end)
!   DT                      : integer               [-]     !Time interval (seconds) 
!
!   COMPUTE_RESIDUAL        : 0/1                   [1]     !Residual appearance in time series: 
!                                                           !0 - Appears not ; 1 - Appears   
!                                                           !(residual is meanless in time series produced: 
!                                                           !should always be 0)
!                                                                                
!   <BeginTimeSerie>
!   NAME                    : char                  [-]     !Name for output time serie file
!   LOCALIZATION_I          : int                   [-]     !Put always 1 (for only one spot)
!   LOCALIZATION_J          : int                   [-]     !Put always 1
!   LOCALIZATION_K          : int                   [-]     !Put always 1
!   <EndTimeSerie>                                              

program ConvertToTimeSerie

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleTimeSerie

    implicit none

    type(T_Time)                                        :: BeginTime, EndTime
    type(T_Time)                                        :: CurrentTime, DataTime
    real                                                :: DT
    logical                                             :: VariableDT
    type (T_Time)                                       :: InitialSystemTime, FinalSystemTime
    real                                                :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                               :: F95Time
    character(PathLength)                               :: FileName
    integer, dimension(24)                              :: DataArray
    character(PathLength)                               :: DataFile
    character(len=StringLength)                         :: FileType
    real, dimension(:,:  ),     pointer                 :: Values2D
    real                                                :: Year, Month = 1, Day = 1
    real                                                :: Hour = 0, Minute = 0, Second = 0
    integer                                             :: FileYear
    integer                                             :: UnitConvert
    real                                                :: ValuesSum
    integer                                             :: ValuesCount
    real                                                :: ValuesAverage
    integer                                             :: ObjTime = 0
    integer                                             :: STAT_CALL
    integer                                             :: ObjTimeSerie = 0
    character(len=StringLength), dimension(:), pointer  :: ParameterList

    ! Definition of type T_TogaFile
    type       T_TogaFile
        integer                                         :: TFileID = 0
        character(len=StringLength)                     :: Name
        type(T_Time)                                    :: StartTime
        type(T_Time)                                    :: EndTime
        type(T_TogaFile), pointer                       :: Next
    end type  T_TogaFile

    type(T_TogaFile), pointer                           :: FirstTogaFile
    type(T_TogaFile), pointer                           :: FirstTSTogaFile, LastTSTogaFile


    call ConstructConvertToTimeSerie
    call ModifyConvertToTimeSerie
    call KillConvertToTimeSerie

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructConvertToTimeSerie

        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        call StartUpMohid("ConvertToTimeSerie")

        call StartCPUTime

        call ReadKeywords

        call StartComputeTime(ObjTime, BeginTime, EndTime, DT, VariableDT, STAT = STAT_CALL)

    end subroutine ConstructConvertToTimeSerie
    
    !--------------------------------------------------------------------------

    subroutine ModifyConvertToTimeSerie
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        integer                                     :: i
        type(T_TogaFile), pointer                   :: ObjectTogaFile

        !Begin-----------------------------------------------------------------

        !Allocates Parameter list
        allocate(ParameterList(1))
        ParameterList(1) = trim(adjustl('water level'))

        !Get date from files, check if date is according with specified
        if (FileType == 'IH') then

            !Opens the data file   
            open(unit = 99, file = FileName, status = 'unknown')

            read(unit = 99, fmt = 15) FileYear

            15 FORMAT(3/, 72X, I2)

            close(unit = 99)
            
            FileYear = 1900 + FileYear
            
            Year = real(FileYear)

            call SetDate(DataTime, Year, Month, Day, Hour, Minute, Second)

            if ((DataTime .lt. BeginTime) .OR. (DataTime .gt. EndTime )) then
                write(*,*) 'Data file does not include data for specified times'
                stop 'ModifyConvertToTimeSerie - ConvertToTimeSerie - ERR01'
            endif  

        else !It is 'TOGA' data type

            !Cycle all files, get date and check is date is correct, order files by date 
            call OpenAndDateFiles 

        endif
        
        !Constructs Time Serie Header
        call StartTimeSerie(ObjTimeSerie, ObjTime,                                  &
                            trim(DataFile),                                         &
                            ParameterList, "cts",                                   &
                            STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                  &
        stop 'ModifyConvertToTimeSerie - ConvertToTimeSerie - ERR02'

        write(*,*)'Reading data from file...'
        
        write(*,*)'Reading '//trim(ParameterList(1))//' fields'

        Running      = .true.
        CurrentTime  = BeginTime

        !Read values for first day
        if (FileType == 'IH') then

            open(unit = 99, file = FileName, status = 'unknown')

            read(unit = 99, fmt = 10) (DataArray(i), i = 1, 24)

            10 FORMAT (3/, 24I3)
            
            UnitConvert = 100 !(Water level must be in meters in TimeSerie)        

        else !It is 'TOGA' data type

            ObjectTogaFile => FirstTSTogaFile

            !Opens the data file   
            open(unit = 99, file = ObjectTogaFile%Name, status = 'unknown')

            read(unit = 99, fmt = 12) (DataArray(i), i = 1, 24)

            12 FORMAT ((/, 21X, 12(I4,1X)),(/, 21X, 12(I4,1X)))
            
            UnitConvert = 1000 !(Water level must be in meters in TimeSerie)                   

        endif 

        !allocate data field
        nullify (Values2D)
        allocate(Values2D(1,1))

        !initializate sum and count of valid values
        ValuesSum = 0.0
        ValuesCount = 0

        do while (Running)

            do i = 1, 24

                !Write data in file only if it is valid
                if (DataArray(i) /= 9999) then

                    !Transfer data to Values2D
                    Values2D(1,1) = real(DataArray(i))/UnitConvert
                    !(Water level must be in meters in TimeSerie)

                    !Actualize sum and count of valid values
                    ValuesSum = ValuesSum + Values2D(1,1)
                    ValuesCount = ValuesCount + 1

                    !Water level is a 2D parameter 
                    call WriteTimeSerie(ObjTimeSerie,                               &
                                Data2D = Values2D,                                  &
                                STAT = STAT_CALL)
                endif

                CurrentTime = CurrentTime + DT

                call ActualizeCurrentTime(ObjTime, DT, STAT = STAT_CALL)
                if(STAT_CALL .ne. SUCCESS_) stop 'ModifyConvertToTimeSerie - ERR03'

            end do

            if ((abs(CurrentTime - EndTime) > DT / 10.) .AND.                       &
                (CurrentTime <= EndTime)) then

                Running = .true.

               !Read next line of data
               if (FileType == 'IH') then

                    read(unit = 99, fmt = 11) (DataArray(i), i = 1, 24)

                    11 FORMAT ((24I3))

                else !It is 'TOGA' data type

                    if (CurrentTime > ObjectTogaFile%EndTime) then
                        !Data is in the next file
                        if (associated(ObjectTogaFile%Next)) then
                            ObjectTogaFile => ObjectTogaFile%Next
                        else
                            exit
                        endif

                        close(unit = 99)
                        open(unit = 99, file = ObjectTogaFile%Name, status = 'unknown')
                        read(unit = 99, fmt = 16) (DataArray(i), i = 1, 24)

                        16 FORMAT ((/, 21X, 12(I4,1X)),(/, 21X, 12(I4,1X)))

                    else !Remain in the same Toga file

                        read(unit = 99, fmt = 13) (DataArray(i), i = 1, 24)

                        13 FORMAT ((21X, 12(I4,1X)),(/, 21X, 12(I4,1X)))

                    endif

                endif

            else
                Running = .false.
            endif

        enddo

        close(unit = 99)
        
        call KillTimeSerie(ObjTimeSerie, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)                                                 &
        stop 'ModifyConvertToTimeSerie - ConvertToTimeSerie - ERR04'
        
        !Calculation of average
        !(average is printed in TimeSerie file)
        ValuesAverage = ValuesSum/(real(ValuesCount))
        write(*, *) 'Average value:'
        write(*, fmt = 20) ValuesAverage
        write(*, *) 'Values count:'
        write(*, fmt = 21) ValuesCount
        20 format (1x, f6.4)
        21 format (1x, I10)
   
    end subroutine ModifyConvertToTimeSerie
    
    !--------------------------------------------------------------------------

    subroutine KillConvertToTimeSerie

        call StopCPUTime

        call ShutdownMohid ("ConvertToTimeSerie", ElapsedSeconds, TotalCPUTime)

    end subroutine KillConvertToTimeSerie
    
    !--------------------------------------------------------------------------

    subroutine StartCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)), &
                                              float(F95Time(3)), float(F95Time(5)), &
                                              float(F95Time(6)), float(F95Time(7))+ &
                                              F95Time(8)/1000.)

    end subroutine StartCPUTime
    
    !--------------------------------------------------------------------------

    subroutine StopCPUTime

        call date_and_time(Values = F95Time)
        
        call SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)), &
                                              float(F95Time(3)), float(F95Time(5)), &
                                              float(F95Time(6)), float(F95Time(7))+ &
                                              F95Time(8)/1000.)
        
        call cpu_time(TotalCPUTime)

        ElapsedSeconds = FinalSystemTime - InitialSystemTime

    end subroutine StopCPUTime
    
    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, iflag
        integer                                     :: ObjEnterData = 0
        integer                                     :: FromFile

        !Begin-----------------------------------------------------------------

        call ReadFileName('IN_MODEL', DataFile, "ConvertToTimeSerie", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ConvertToTimeSerie - ERR01'

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ConvertToTimeSerie - ERR02'

        call GetExtractType     (FromFile = FromFile)

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,    &
                                 VariableDT, "ConvertToTimeSerie")

        ! Obtain input file type
        call GetData(FileType,                                                      &
                     ObjEnterData, iflag,                                           &
                     SearchType   = FromFile,                                       &
                     keyword      ='INPUT_FILE_TYPE',                               &
                     ClientModule ='ConvertToTimeSerie',                            &
                     default      = 'IH',                                           &
                     STAT         = STAT_CALL)             
        if (STAT_CALL /= SUCCESS_)  stop 'ReadKeywords - ConvertToTimeSerie - ERR03'

        FileType = trim(adjustl(FileType))
        if ((FileType /= 'IH') .AND. (FileType /= 'TOGA')) then
           write(*,*) 'Type '//trim(FileType)//' not a valid data file type'
           stop 'ReadKeywords - ConvertToTimeSerie - ERR04'
        endif

        if (FileType == 'IH') then

            ! Obtain data file name
            call GetData(FileName,                                                  &
                         ObjEnterData, iflag,                                       &
                         SearchType   = FromFile,                                   &
                         keyword      = 'INPUT_FILENAME',                           &
                         ClientModule = 'ConvertToTimeSerie',                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_)                                            &
            stop 'ReadKeywords - ConvertToTimeSerie - ERR05'

        else !It is 'TOGA' data type

            ! Obtain data files names: several files are admited
            call ReadTogaFileName(ObjEnterData)

        endif

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - ConvertToTimeSerie - ERR06'

    end subroutine ReadKeywords

    !--------------------------------------------------------------------------
    
    subroutine ReadTogaFileName(ObjEnterData)

        !Arguments-------------------------------------------------------------
         
        integer                                     :: ObjEnterData
          
        !Local-----------------------------------------------------------------
        integer                                 :: STAT_CALL, ClientNumber
        type (T_TogaFile),       pointer        :: NewTogaFile
        logical                                 :: BlockFound
        logical                                 :: AtLeastOneBlock = .false.

        !Begin-----------------------------------------------------------------

do1 :   do
            call ExtractBlockFromBuffer(ObjEnterData,                               &
                                        ClientNumber    = ClientNumber,             &
                                        block_begin     = '<BeginTogaFile>',        &
                                        block_end       = '<EndTogaFile>',          &
                                        BlockFound      = BlockFound,               &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    AtLeastOneBlock = .true.

                    call AddTogaFile                     (NewTogaFile)

                    call ConstructTogaFile               (NewTogaFile, ObjEnterData)

                    nullify(NewTogaFile)

                else cd2
                    call Block_Unlock(ObjEnterData,                                 & 
                                      ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                    &
                    stop 'ReadTogaFileName - ConvertToTimeSerie - ERR01'

                    exit do1

                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'ReadTogaFileName - ConvertToTimeSerie - ERR02'
            else cd1
                stop 'ReadTogaFileName - ConvertToTimeSerie - ERR03'
            end if cd1

        end do do1

        if (.not. AtLeastOneBlock) then                                            
            write(*,*) 'No data file block is indicated in input file. '
            stop 'ReadTogaFileName - ConvertToTimeSerie - ERR04'
        end if

    end subroutine ReadTogaFileName

    !--------------------------------------------------------------------------

    subroutine AddTogaFile(ObjTogaFile)

        !Arguments-------------------------------------------------------------
        type (T_TogaFile),     pointer           :: ObjTogaFile

        !Local-----------------------------------------------------------------
        type (T_TogaFile),     pointer           :: PreviousTogaFile
        type (T_TogaFile),     pointer           :: NewTogaFile

        !Begin-----------------------------------------------------------------

        !Allocates new TogaFile
        allocate (NewTogaFile)
        nullify  (NewTogaFile%Next)

        !Insert new file into list and makes current ?? point to it
        if (.not. associated(FirstTogaFile)) then
            FirstTogaFile            => NewTogaFile
            ObjTogaFile              => NewTogaFile
        else
            PreviousTogaFile         => FirstTogaFile
            ObjTogaFile              => FirstTogaFile%Next
            do while (associated(ObjTogaFile))
                PreviousTogaFile     => ObjTogaFile
                ObjTogaFile          => ObjTogaFile%Next
            enddo
            ObjTogaFile              => NewTogaFile
            PreviousTogaFile%Next    => NewTogaFile
        end if

    end subroutine AddTogaFile

    !--------------------------------------------------------------------------

    subroutine ConstructTogaFile(NewTogaFile,ObjEnterData)

        !Arguments-------------------------------------------------------------
        type (T_TogaFile),      pointer           :: NewTogaFile
        integer                                     :: ObjEnterData

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        
        ! Obtain data file name
        call GetData(NewTogaFile%Name,                                              &
                     ObjEnterData, iflag,                                           &
                     SearchType   = FromBlock,                                      &
                     keyword      = 'INPUT_FILENAME',                               &
                     ClientModule = 'ConvertToTimeSerie',                           &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_)                                                &
        stop 'ConstructTogaFile - ConvertToTimeSerie - ERR01'

    end subroutine ConstructTogaFile

    !--------------------------------------------------------------------------

    subroutine OpenAndDateFiles

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_TogaFile), pointer                  :: TogaFileX
        logical                                     :: exist, FirstTime
        integer                                     :: fileunit
        real                                        :: MonthEnd = 12, DayEnd = 31 
        real                                        :: HourEnd = 23
      
        !Begin-----------------------------------------------------------------

        FirstTime = .true.

        TogaFileX => FirstTogaFile
        
        fileunit = 99

        !In a DO cycle open all HDF5 files provided by the user
        do while (associated(TogaFileX))

            !Verifies if file exists
            inquire(FILE = TogaFileX%Name, EXIST = exist)
            if (.not. exist) then
                write(*,*)'Toga file does not exist:'//trim(TogaFileX%Name)
                stop 'OpenAndDateFiles - ConvertToTimeSerie - ERR01'
            endif

            TogaFileX%TFileID = fileunit

            !Opens the data file   
            open(unit = TogaFileX%TFileID, file = TogaFileX%Name, status = 'unknown')

            read(unit = TogaFileX%TFileID, fmt = 14) FileYear

            14 FORMAT(44X,I4)

            close(unit = TogaFileX%TFileID)

            Year = real(FileYear)

            call SetDate(TogaFileX%StartTime, Year, Month, Day, Hour, Minute, Second)

            if ((TogaFileX%StartTime .lt. BeginTime) .OR.                           & 
                (TogaFileX%StartTime .gt. EndTime )) then
                write(*,*) 'Data file does not include data for specified times'
                stop 'OpenAndDateFiles - ConvertToTimeSerie - ERR02'
            endif
            
            call SetDate(TogaFileX%EndTime, Year, MonthEnd, DayEnd, HourEnd, Minute, Second)
            
            !Add file to list of relevant files
            call AddTSTogaFile(TogaFileX)

            TogaFileX => TogaFileX%Next

        end do

    end subroutine OpenAndDateFiles

    !--------------------------------------------------------------------------

    subroutine AddTSTogaFile(TogaFileX)

        !Arguments-------------------------------------------------------------
        type(T_TogaFile), pointer                   :: TogaFileX

        !Local-----------------------------------------------------------------
        type(T_TogaFile), pointer                   :: TogaFileAux
        type(T_TogaFile), pointer                   :: PreviousTogaFile, LastTogaFile
      
        !Begin-----------------------------------------------------------------

        if (.not. associated(FirstTSTogaFile)) then

            call CreateTSTogaFile(FirstTSTogaFile, TogaFileX)
            call CreateTSTogaFile(LastTSTogaFile, TogaFileX)
            deallocate(FirstTSTogaFile%Next) 
            nullify(FirstTSTogaFile%Next)

        else

            if (TogaFileX%StartTime < FirstTSTogaFile%StartTime) then
                !current file should be the first file in list

                !save the previous list 
                allocate(TogaFileAux)
                call CreateTSTogaFile(TogaFileAux, FirstTSTogaFile)
                TogaFileAux%Next => FirstTSTogaFile%Next              

                !make the first element in the list of relevant files equal to current file
                call CreateTSTogaFile(FirstTSTogaFile, TogaFileX)
                FirstTSTogaFile%Next => TogaFileAux

            else
                !check next files in list

                !locate previous file in the first file
                allocate(PreviousTogaFile)
                PreviousTogaFile => FirstTSTogaFile                   

                do while(associated(PreviousTogaFile))

                    if (.not. associated(PreviousTogaFile%Next)) then
        
                        !current file is the last file in the list of relevant files
                        call CreateTSTogaFile(PreviousTogaFile%Next, TogaFileX)
                        allocate(LastTogaFile)
                        LastTogaFile => PreviousTogaFile%Next
                        deallocate(LastTogaFile%Next)
                        nullify(LastTogaFile%Next)
                        call CreateTSTogaFile(LastTSTogaFile, TogaFileX)

                        !current file was added to list
                        exit

                    else

                        !check if current file should be located before the next file
                        if (TogaFileX%StartTime < PreviousTogaFile%Next%StartTime) then
                            !current file should be located before next file

                            !save the previous list begining in PreviousHDF5File%Next
                            allocate(LastTogaFile)
                            LastTogaFile => PreviousTogaFile%Next
                            allocate(TogaFileAux)

                            call CreateTSTogaFile(TogaFileAux, LastTogaFile)
                            TogaFileAux%Next => LastTogaFile%Next

                            call CreateTSTogaFile(LastTogaFile, TogaFileX)

                            PreviousTogaFile%Next => LastTogaFile
                            LastTogaFile%Next => TogaFileAux

                            !current file was added to list
                            exit

                        end if

                    end if

                    !check next file in list
                    PreviousTogaFile => PreviousTogaFile%Next     

                end do

            end if

        end if

    end subroutine AddTSTogaFile

    !--------------------------------------------------------------------------

    subroutine CreateTSTogaFile(TogaFileNew, TogaFileX)

        !Arguments-------------------------------------------------------------
        type(T_TogaFile), pointer                   :: TogaFileX
        type(T_TogaFile), pointer                   :: TogaFileNew

        !Local-----------------------------------------------------------------
      
        !Begin-----------------------------------------------------------------

        !This subroutine atributes the values of the fields of a HDF5File to another HDF5File
        allocate(TogaFileNew)

        TogaFileNew%Name             =  TogaFileX%Name
        TogaFileNew%StartTime        =  TogaFileX%StartTime
        TogaFileNew%EndTime          =  TogaFileX%EndTime

    end subroutine CreateTSTogaFile

    !--------------------------------------------------------------------------

end program ConvertToTimeSerie

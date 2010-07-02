!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Swat Link
! PROJECT       : SwatLink
! PROGRAM       : MainSwatLink
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2005
! REVISION      : Frank Braunschweig
! DESCRIPTION   : Interpolates modified SWAT output files in time and creates 
!                 discharges locations
!
!------------------------------------------------------------------------------
! Keywords read in the Data File
!
! Keyword                   : Data Type             Default     !Comment
!
! START                     : YYYY MM DD hh mm ss   -           !Start Date of Conversion
! END                       : YYYY MM DD hh mm ss   -           !End Date of Conversion
! DT                        : sec                   -           !Dummy for DT (not used)
! SWAT_FILES_PATH           : char                  -           !Path to folder where SWAT results are located
                                                                !must end with "\"
! MOHID_FILES_PATH          : char                  -           !Path to folder where MOHID input files are to be placed
                                                                !must end with "\"
! NETWORK_FILE              : char                  -           !Path to file which defines the Drainage Network
! DISCHARGES_FILE           : char                  -           !Discharges File to be written for Posterior Usage with
                                                                !MOHID River Network
! XYZ_FILE                  : char                  -           !XYZ File which defines the location of the discharges
! RAIN_DATA_FILE            : char                  -           !Path to Time Series file with hourly Rainfall Data
! RAIN_COLUMN               : int                   -           !Column in which rain values are stored
! RUN_WATER_QUALITY         : logical               0           !It to run Water Quality




program MohidSwatLink

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleTimeSerie
    use ModuleDrawing

    implicit none

    integer, parameter                              :: SwatColOverLand_     = 19
    integer, parameter                              :: SwatColGroundWater_  =  9

    integer, parameter                              :: MohidColOverland_    =  2
    integer, parameter                              :: MohidColGroundWater_ =  3

    integer, parameter                              :: FirstPropertyInSwat  = 12
    !integer, parameter                              :: SwatColGroundW_      = 1
    !integer, parameter                              :: SwatColNitrate_      = 
    !integer, parameter                              :: SwatColNitrate_      = 

    type T_Node
        integer                                     :: ID                       = null_int
        real                                        :: X                        = null_real
        real                                        :: Y                        = null_real
        type (T_Node), pointer                      :: Next
    end type T_Node

    type T_SwatValue
        real                                        :: OverlandFlow
        real                                        :: GroundWaterFlow
        real, dimension(6)                          :: PropertyValues
        type (T_TIME)                               :: Date
    end type T_SwatValue                            
                                                    
    type T_DischargeProperty                        
        type (T_PropertyID)                         :: ID
        logical                                     :: VariableOverLandConc
        integer                                     :: SwatSeriesColumn
        real                                        :: SwatToMohidFactor
        real                                        :: OLConcentration 
        real                                        :: GWConcentration
        type (T_DischargeProperty), pointer         :: Next => null()
    end type T_DischargeProperty                    
          
                                                    
    type(T_Time)                                    :: BeginTime, EndTime, CurrentTime
    real                                            :: DT, MaxDT
    logical                                         :: VariableDT
    character(PathLength)                           :: ModelName            = 'SwatLink'
                                                    
    type (T_Time)                                   :: InitialSystemTime, FinalSystemTime
    real                                            :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)                           :: F95Time
    integer                                         :: ObjTime = 0
                                                    
    character(PathLength), parameter                :: ProjectFile  = "MohidSwatLink.dat"
    character(PathLength)                           :: DNetFile, XYZLocFile, DischargesFile
    character(PathLength)                           :: SwatFilesPath, MohidFilesPath

    character(PathLength)                           :: RainDataFile
    integer                                         :: RainDataColumn

    logical                                         :: RunWaterQuality
                                                    
    integer                                         :: NumberOfDays
    type (T_DischargeProperty), pointer             :: FirstProperty                => null()
    type (T_XYZPoints), pointer                     :: DischargesLocations          => null()
    type (T_SwatValue), dimension(:), allocatable   :: SwatValues        
    real, dimension(:), allocatable                 :: DailyRainValues


    type (T_Node), pointer                          :: FirstNode                    => null()

    call ConstructMohidSwatLink
    call ModifyMohidSwatLink
    call KillMohidSwatLink

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidSwatLink
        

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        call StartUpMohid("MohidSwatLink")

        call StartCPUTime

        call ReadKeywords

        call StartComputeTime(ObjTime, BeginTime, EndTime, DT,     &
                              VariableDT, MaxDT, STAT = STAT_CALL)   
        if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - Mohid - ERR0'


        !Reads Discharges Location File
        call New(DischargesLocations, XYZLocFile)

        !Calculates Number Of Days
        NumberOfDays = int((EndTime - BeginTime) / 86400.0) + 1

    end subroutine ConstructMohidSwatLink
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidSwatLink
        
        !Local-----------------------------------------------------------------
        integer                                     :: iP, iProp
        integer                                     :: STAT_CALL
        integer, dimension(:), allocatable          :: ObjTimeSeriesIn
        integer, dimension(:), allocatable          :: UnitsTimeSeriesOut
        character(5)                                :: AuxString
        type (T_TIME)                               :: RainInitialData
        type (T_TIME)                               :: Time1, Time2
        real                                        :: Value1, Value2, AuxValue, Year, Month, Day
        logical                                     :: TimeCycle
        integer                                     :: iDay, nRainValues, iValue
        real                                        :: CurrentDay, LastDay, CurrentFlowValue
        real, dimension(:,:), pointer               :: RainValuesMatrix
        real, dimension(:), allocatable             :: AuxPropConc
        character(20), dimension(:), allocatable    :: PropertyHeader
        integer                                     :: nVariableProps, iAux, iCh
        type (T_DischargeProperty), pointer         :: CurrProperty
        real                                        :: OverlandFlow
        integer                                     :: i

        write(*, *)"-------------------------- MOHID -------------------------"
        write(*, *)
        write(*, *)"Running MOHID SwatLink, please wait..."
        write(*, *)                    


        call ReadRainTimeSerie (nRainValues, RainValuesMatrix, RainInitialData)

        call WriteDischargeFile


        allocate (SWATValues     (NumberOfDays))
        SWATValues(:)%OverlandFlow    = 0.0
        SWATValues(:)%GroundWaterFlow = 0.0

        allocate (ObjTimeSeriesIn     (DischargesLocations%Count))    
        allocate (UnitsTimeSeriesOut  (DischargesLocations%Count))    
        ObjTimeSeriesIn  = 0
        UnitsTimeSeriesOut = 0
        nVariableProps = 0

        if (.not. RunWaterQuality) then
            allocate(PropertyHeader(3))
            PropertyHeader(1) = "Seconds"
            PropertyHeader(2) = "Overland_Flow"
            PropertyHeader(3) = "GroundWater_Flow"
        else

            !See how many properties have variable overland conc
            CurrProperty => FirstProperty
            do while (associated(CurrProperty))
                if (CurrProperty%VariableOverLandConc) then
                    nVariableProps = nVariableProps + 1
                endif
                CurrProperty => CurrProperty%Next
            enddo

            allocate(PropertyHeader(3 + nVariableProps))
            PropertyHeader(1) = "Seconds"
            PropertyHeader(2) = "Overland_Flow"
            PropertyHeader(3) = "GroundWater_Flow"

            nVariableProps = 0
            CurrProperty => FirstProperty
            do while (associated(CurrProperty))
                if (CurrProperty%VariableOverLandConc) then
                    nVariableProps = nVariableProps + 1
                    PropertyHeader(3 + nVariableProps) = CurrProperty%ID%Name
                endif
                CurrProperty => CurrProperty%Next
            enddo
        endif

        !Changes white spaces to underscores in the property list (Turns Input to Excel more easy)
        do iP = 1, nVariableProps + 3
            do iCh = 1, len_trim(PropertyHeader(iP))
                if (PropertyHeader(iP)(iCh:iCh) == ' ') then
                    PropertyHeader(iP)(iCh:iCh) =  '_'
                endif
            enddo
        enddo

        allocate (AuxPropConc(nVariableProps))

        do iP = 1, DischargesLocations%Count

            write(AuxString, fmt='(i5)') iP
            write(*,*)'Writing ',trim(adjustl(AuxString))//".dis",'...'

            !Opens Discharges File produced by SWAT
            call StartTimeSerieInput (ObjTimeSeriesIn(iP), trim(adjustl(SwatFilesPath))//trim(adjustl(AuxString))//".dis", ObjTime, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR02'

            !Reads SWAT Flow Values
            CurrentTime = BeginTime
            iDay        = 1
            do while (CurrentTime <= EndTime)

                !Gets Overland Flow values
                call GetTimeSerieValue (ObjTimeSeriesIn(iP), CurrentTime, SwatColOverLand_,     &
                                        Time1, Value1, Time2, Value2, TimeCycle, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR03'

                call InterpolateValueInTime(CurrentTime, Time1, Value1, Time2, Value2, SWATValues(iDay)%OverlandFlow)

                !Gets GroundWater Flow values
                call GetTimeSerieValue (ObjTimeSeriesIn(iP), CurrentTime, SwatColGroundWater_,  &
                                        Time1, Value1, Time2, Value2, TimeCycle, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR03'

                call InterpolateValueInTime(CurrentTime, Time1, Value1, Time2, Value2, SWATValues(iDay)%GroundWaterFlow)

                !Gets Property Concentrations
                if (RunWaterQuality) then
                    do iProp = 1, 6
                        call GetTimeSerieValue (ObjTimeSeriesIn(iP), CurrentTime, iProp + FirstPropertyInSwat,  &
                                                Time1, Value1, Time2, Value2, TimeCycle, STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR04'

                        call InterpolateValueInTime(CurrentTime, Time1, Value1, Time2, Value2, SWATValues(iDay)%PropertyValues(iProp))
                    enddo
                endif

                SWATValues(iDay)%Date = CurrentTime

                CurrentTime = CurrentTime + 86400
                iDay        = iDay        + 1

                !call ExtractDate(CurrentTime, Year, Month, Day)                
                !open (unit=111, file="a.dat")
                !write (111,*) iDay, Year, Month, Day


            enddo



            !Opens Final (Interpolated) Discharge File
            call UnitsManager (UnitsTimeSeriesOut(iP), OPEN_FILE)
            
            write(AuxString, fmt='(i5)') iP

            open (unit = UnitsTimeSeriesOut(iP), file = trim(adjustl(MohidFilesPath))//trim(adjustl(AuxString))//".ssl")

            call WriteDataLine(UnitsTimeSeriesOut(iP), "Interpolated Time Series")

            call WriteDataLine(UnitsTimeSeriesOut(iP), "LOCALIZATION_I", 1)
            call WriteDataLine(UnitsTimeSeriesOut(iP), "LOCALIZATION_J", 1)
            call WriteDataLine(UnitsTimeSeriesOut(iP), "SERIE_INITIAL_DATA", BeginTime)
            call WriteDataLine(UnitsTimeSeriesOut(iP), "TIME_UNITS", "SECONDS")
            call WriteDataLine(UnitsTimeSeriesOut(iP), " ")

            if (.not. RunWaterQuality) then
                write(UnitsTimeSeriesOut(iP), fmt=999)(PropertyHeader(iProp), iProp = 1, 3)
            else
                write(UnitsTimeSeriesOut(iP), fmt=999)(PropertyHeader(iProp), iProp = 1, 3 + nVariableProps)
            endif
            
            call WriteDataLine(UnitsTimeSeriesOut(iP), "<BeginTimeSerie>")


            LastDay = -99
            iDay    = 0

            do iValue = 1, nRainValues

                CurrentTime = RainInitialData + RainValuesMatrix(iValue, 1)

                if (CurrentTime >= BeginTime .and. CurrentTime <= EndTime) then

                    call ExtractDate(CurrentTime, Day = CurrentDay)

                    if (CurrentDay == LastDay) then
                        !Do Nothing
                    else
                        LastDay               = CurrentDay
                        iDay                  = iDay + 1
                    endif
                    
                    OverlandFlow = SWATValues(iDay)%OverlandFlow
                    
                    if (SWATValues(iDay)%Date /= CurrentTime) then
                        do i=1,NumberOfDays
                            if (SWATValues(i)%Date == CurrentTime) then
                                OverlandFlow = SWATValues(i)%OverlandFlow
                            endif
                        enddo
                    endif

                    !call ExtractDate(CurrentTime, Year, Month, Day)                
                    !open (unit=112, file="b.dat")
                    !write (112,*) iDay, Year, Month, Day

                    !Calculates Hourly Flow
                    if (DailyRainValues(iDay) /= 0.0) then
                        CurrentFlowValue = 86400 / 3600 * OverlandFlow * RainValuesMatrix(iValue, RainDataColumn) / DailyRainValues(iDay)
                    else
                        CurrentFlowValue = 0.0
                    endif

                    if (.not. RunWaterQuality) then
                        !Write DT since beginning, OverlandFlow & GroundWaterFlow
                        write(UnitsTimeSeriesOut(iP), fmt=1000)CurrentTime - BeginTime, CurrentFlowValue,   &
                                                               SWATValues(iDay)%GroundWaterFlow
                    else

                        !Converts SWAT properties to MOHID properties
                        nVariableProps = 0
                        CurrProperty => FirstProperty
                        do while (associated(CurrProperty))
                            if (CurrProperty%VariableOverLandConc) then
                                nVariableProps = nVariableProps + 1
                                iAux           = CurrProperty%SwatSeriesColumn - FirstPropertyInSwat
                                AuxPropConc(nVariableProps) = SWATValues(iDay)%PropertyValues(iAux) *       &
                                                              CurrProperty%SwatToMohidFactor
                            endif
                            CurrProperty => CurrProperty%Next
                        enddo

                        write(UnitsTimeSeriesOut(iP), fmt=1000)CurrentTime - BeginTime, CurrentFlowValue,   &
                                                               SWATValues(iDay)%GroundWaterFlow,            &
                                                               (AuxPropConc(iProp), iProp = 1, nVariableProps)
                    endif

                endif
            enddo    

            call WriteDataLine(UnitsTimeSeriesOut(iP), "<EndTimeSerie>")
            call UnitsManager (UnitsTimeSeriesOut(iP), CLOSE_FILE)

        enddo

    999 format(2x,   A13, 1x, 1000(1x, A20))
             
   1000 format(1x, f13.2, 1x, 1000(1x, e20.12e3))

   
    end subroutine ModifyMohidSwatLink

    !--------------------------------------------------------------------------

    subroutine ReadRainTimeSerie (nRainValues, RainValuesMatrix, RainInitialData)

        !Arguments-------------------------------------------------------------
        integer                                     :: nRainValues
        real, dimension(:,:), pointer               :: RainValuesMatrix
        type (T_TIME)                               :: RainInitialData

        !Local---------------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ObjRainTimeSerie = 0
        real                                        :: LastDay, CurrentDay
        integer                                     :: iDay, iValue

        allocate (DailyRainValues(NumberOfDays))

        DailyRainValues             = 0.0

        !Reads rain values
        call StartTimeSerieInput(ObjRainTimeSerie, RainDataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR01'

        call GetTimeSerieDataMatrix (ObjRainTimeSerie, RainValuesMatrix, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR03'

        !Gets Number of rows
        call GetTimeSerieDataValues (ObjRainTimeSerie, nRainValues, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR04'

        !Gets initial date
        call GetTimeSerieInitialData(ObjRainTimeSerie, RainInitialData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'MohidSwatLink - ModifyMohidSwatLink - ERR05'

        LastDay = -99
        iDay    = 0
        do iValue = 1, nRainValues

            CurrentTime = RainInitialData + RainValuesMatrix(iValue, 1)

            if (CurrentTime >= BeginTime .and. CurrentTime <= EndTime) then
                call ExtractDate(CurrentTime, Day = CurrentDay)
                if (CurrentDay == LastDay) then
                    DailyRainValues(iDay) = DailyRainValues(iDay) + RainValuesMatrix(iValue, RainDataColumn)
                else
                    LastDay               = CurrentDay
                    iDay                  = iDay + 1
                    DailyRainValues(iDay) = RainValuesMatrix(iValue, RainDataColumn)
                endif
            endif
        enddo


    end subroutine ReadRainTimeSerie
    
    !--------------------------------------------------------------------------

    subroutine WriteDischargeFile

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: ObjNetworkFile = 0
        logical                                     :: BlockFound  
        integer                                     :: ClientNumber  
        integer                                     :: STAT_CALL, flag, UnitID
        type (T_Node), pointer                      :: NewNode, CurrNode
        real, dimension(2)                          :: AuxCoord
        integer                                     :: iP    
        character(5)                                :: AuxString
        integer                                     :: NearNodeID
        real                                        :: Dist, MinDist
        integer                                     :: FirstLine, LastLine
        type (T_DischargeProperty), pointer         :: CurrProperty
        integer                                     :: iAux


        !Reads Drainage Network File
        call ConstructEnterData (ObjNetworkFile, DNetFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'WriteDischargeFile - MohidSwatLink - ERR01'

        do 

            call ExtractBlockFromBuffer(ObjNetworkFile, ClientNumber,     &
                                        '<BeginNode>', '<EndNode>', BlockFound,                  &
                                        FirstLine, LastLine, STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'WriteDischargeFile - MohidSwatLink - ERR02'

            if (BlockFound) then                 
                    
                allocate    (NewNode)
                nullify     (NewNode%Next)

                !Gets ID
                call GetData(NewNode%ID,                                                &
                             ObjNetworkFile, flag,                                      & 
                             keyword      = 'ID',                                       &
                             ClientModule = 'MohidSwatLink',                            &
                             SearchType   = FromBlock,                                  &
                             STAT         = STAT_CALL)                                  
                if (STAT_CALL .NE. SUCCESS_) stop 'WriteDischargeFile - MohidSwatLink - ERR03'

                !Gets X / Y
                call GetData(AuxCoord,                                                  &
                             ObjNetworkFile, flag,                                      &
                             keyword      = 'COORDINATES',                              &
                             ClientModule = 'MohidSwatLink',                            &
                             SearchType   = FromBlock,                                  &
                             STAT         = STAT_CALL)                                  
                if (STAT_CALL .NE. SUCCESS_) stop 'WriteDischargeFile - ConstructNode - ERR04'

                NewNode%X = AuxCoord(1)                                                                                
                NewNode%Y = AuxCoord(2)                                                                                
                call AddNodeToList (NewNode)

            else

                exit

            endif

        end do 


        !Opens Discharge Data File
        call UnitsManager (UnitID, OPEN_FILE)
        open (file = DischargesFile, unit = UnitID)
    
        call WriteDataLine(UnitID, "Automaticly generated discharge file")

        do iP = 1, DischargesLocations%Count

            MinDist = -null_real
            CurrNode => FirstNode
            do while (associated(CurrNode))
                
                Dist = sqrt((DischargesLocations%X(iP) - CurrNode%X) ** 2.0 + (DischargesLocations%Y(iP) - CurrNode%Y) ** 2.0)
                
                if (Dist < MinDist) then
                    MinDist    = Dist
                    NearNodeID = CurrNode%ID
                endif

                CurrNode => CurrNode%Next
            enddo



            !Writes Overland Flow discharge
            call WriteDataLine(UnitID, "<begindischarge>")
            write(AuxString, fmt='(i5)') NearNodeID
            call WriteDataLine(UnitID, "NAME",          "NODE_"//trim(adjustl(AuxString)))
            call WriteDataLine(UnitID, "DESCRIPTION",   "Overland & Lateral Discharge")
            call WriteDataLine(UnitID, "NODE_ID",        NearNodeID)
            write(AuxString, fmt='(i5)')int(DischargesLocations%Z(iP))
            call WriteDataLine(UnitID, "DATA_BASE_FILE", trim(adjustl(MohidFilesPath))//trim(adjustl(AuxString))//".ssl")
            call WriteDataLine(UnitID, "FLOW_COLUMN",    MohidColOverland_)
            
            if (RunWaterQuality) then
                
                iAux = 3
                CurrProperty => FirstProperty
                do while (associated(CurrProperty))
                    call WriteDataLine(UnitID, "")
                    call WriteDataLine(UnitID, "<<beginproperty>>")
                    call WriteDataLine(UnitID, "NAME",          CurrProperty%ID%Name)
                    call WriteDataLine(UnitID, "UNITS",         "SI Units")
                    call WriteDataLine(UnitID, "DESCRIPTION",   "No description given")
                    call WriteDataLine(UnitID, "DEFAULTVALUE",  CurrProperty%OLConcentration)
                    if (CurrProperty%VariableOverLandConc) then
                        iAux = iAux + 1
                        call WriteDataLine(UnitID, "TIME_SERIE_COLUMN", iAux)
                    endif
                    call WriteDataLine(UnitID, "<<endproperty>>")
                    CurrProperty => CurrProperty%Next
                enddo

            endif
            call WriteDataLine(UnitID, "<enddischarge>")
            call WriteDataLine(UnitID, "")


            !Writes GW Flow discharge
            call WriteDataLine(UnitID, "<begindischarge>")
            write(AuxString, fmt='(i5)') NearNodeID
            call WriteDataLine(UnitID, "NAME",          "NODE_"//trim(adjustl(AuxString)))
            call WriteDataLine(UnitID, "DESCRIPTION",   "GW Discharge")
            call WriteDataLine(UnitID, "NODE_ID",        NearNodeID)
            write(AuxString, fmt='(i5)')int(DischargesLocations%Z(iP))
            call WriteDataLine(UnitID, "DATA_BASE_FILE", trim(adjustl(MohidFilesPath))//trim(adjustl(AuxString))//".ssl")
            call WriteDataLine(UnitID, "FLOW_COLUMN",    MohidColGroundWater_)
            
            if (RunWaterQuality) then

                CurrProperty => FirstProperty
                do while (associated(CurrProperty))
                    call WriteDataLine(UnitID, "")
                    call WriteDataLine(UnitID, "<<beginproperty>>")
                    call WriteDataLine(UnitID, "NAME",          CurrProperty%ID%Name)
                    call WriteDataLine(UnitID, "UNITS",         "SI Units")
                    call WriteDataLine(UnitID, "DESCRIPTION",   "No description given")
                    call WriteDataLine(UnitID, "DEFAULTVALUE",  CurrProperty%GWConcentration)
                    call WriteDataLine(UnitID, "<<endproperty>>")
                    CurrProperty => CurrProperty%Next
                enddo

            endif
            call WriteDataLine(UnitID, "<enddischarge>")
            call WriteDataLine(UnitID, "")

        enddo

        call UnitsManager (UnitID, CLOSE_FILE)
    
    end subroutine WriteDischargeFile

    !--------------------------------------------------------------------------

    subroutine AddNodeToList (NewNode)

        !Arguments-------------------------------------------------------------
        type (T_Node), pointer                      :: NewNode

        !Local-----------------------------------------------------------------
        type (T_Node), pointer                      :: PreviousNode
        type (T_Node), pointer                      :: CurrNode

        !Begin-----------------------------------------------------------------

        if (.not. associated(FirstNode)) then
            FirstNode            => NewNode
        else
            PreviousNode         => FirstNode
            CurrNode             => FirstNode%Next
            do while (associated(CurrNode))
                PreviousNode     => CurrNode
                CurrNode         => CurrNode%Next
            enddo
            PreviousNode%Next    => NewNode
        end if

    end subroutine AddNodeToList
    
    !--------------------------------------------------------------------------

    subroutine KillMohidSwatLink

        call StopCPUTime

        call ShutdownMohid ("MohidSwatLink", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidSwatLink
    
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
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: ClientNumber, flag
        logical                                     :: BlockFound
        type (T_DischargeProperty), pointer         :: NewProperty

        !Opens Data File
        call ConstructEnterData (ObjEnterData, ProjectFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR01'

        !Reads Time Keywords
        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, ModelName, MaxDT)

        !Gets Path to Drainage Network File
        call GetData    (DNetFile, ObjEnterData, flag,                                      &
                         keyword      = 'NETWORK_FILE',                                     &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR01'

        !Gets Path to XYZ file which contains information about subbasin outlets
        call GetData    (XYZLocFile, ObjEnterData, flag,                                    &
                         keyword      = 'XYZ_FILE',                                         &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR02'

        !Gets Path to Discharge file to produce
        call GetData    (DischargesFile, ObjEnterData, flag,                                &
                         keyword      = 'DISCHARGES_FILE',                                  &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR03'

        !Gets Path to time series file with rain information
        call GetData    (RainDataFile, ObjEnterData, flag,                                  &
                         keyword      = 'RAIN_DATA_FILE',                                   &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)                                          
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR05'              
                                                                                            
        !Gets column where rain data values are stored                                      
        call GetData    (RainDataColumn, ObjEnterData, flag,                                &
                         keyword      = 'RAIN_COLUMN',                                      &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR06'

        !Gets Path to place where Swat result files are located
        call GetData    (SwatFilesPath, ObjEnterData, flag,                                 &
                         keyword      = 'SWAT_FILES_PATH',                                  &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR06a'


        !Gets Path to place where Swat result files are located
        call GetData    (MohidFilesPath, ObjEnterData, flag,                                &
                         keyword      = 'MOHID_FILES_PATH',                                 &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR06a'


        !Verifies if the user wants to run water quality
        call GetData    (RunWaterQuality, ObjEnterData, flag,                               &
                         keyword      = 'RUN_WATER_QUALITY',                                &
                         ClientModule = 'MohidSwatLink',                                    &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR07'

        if (RunWaterQuality) then
        
            !Searches for properties
            do
                call ExtractBlockFromBuffer(ObjEnterData,                                   &
                                            ClientNumber    = ClientNumber,                 &
                                            block_begin     = '<beginproperty>',            &
                                            block_end       = '<endproperty>',              &
                                            BlockFound      = BlockFound,                   &
                                            STAT            = STAT_CALL)
                if (STAT_CALL == SUCCESS_) then
                    
                    if (BlockFound) then
                        call AddDischargeProperty       (NewProperty)
                        call ConstructDischargeProperty (NewProperty, ObjEnterData)
                    else
                        call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                        if (STAT_CALL .NE. SUCCESS_) stop 'ReadKeywords - MohidSwatLink - ERR08'
                        exit
                    end if

                else if (STAT_CALL .EQ. BLOCK_END_ERR_) then

                    write(*,*)  
                    write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop 'ReadKeywords - MohidSwatLink - ERR09'
            
                else
                
                    stop 'ReadKeywords - MohidSwatLink - ERR10'

                end if

            end do

        endif




    end subroutine ReadKeywords

    !--------------------------------------------------------------------------

    subroutine AddDischargeProperty (ObjProperty)

        !Arguments-------------------------------------------------------------
        type (T_DischargeProperty), pointer         :: ObjProperty

        !Local-----------------------------------------------------------------
        type (T_DischargeProperty), pointer         :: PreviousProperty
        type (T_DischargeProperty), pointer         :: NewProperty

        !Begin-----------------------------------------------------------------

        !Allocates new Parameter
        allocate (NewProperty)
        nullify  (NewProperty%Next)

        if (.not. associated(FirstProperty)) then
            FirstProperty            => NewProperty
            ObjProperty              => NewProperty
        else
            PreviousProperty         => FirstProperty
            ObjProperty              => FirstProperty%Next
            do while (associated(ObjProperty))
                PreviousProperty     => ObjProperty
                ObjProperty          => ObjProperty%Next
            enddo
            ObjProperty              => NewProperty
            PreviousProperty%Next    => NewProperty
        end if

    end subroutine AddDischargeProperty

    !--------------------------------------------------------------------------

    subroutine ConstructDischargeProperty (ObjProperty, ObjEnterData)

        !Arguments-------------------------------------------------------------
        type (T_DischargeProperty), pointer         :: ObjProperty
        integer                                     :: ObjEnterData
        integer                                     :: flag, STAT_CALL

        !Local-----------------------------------------------------------------

        call ConstructPropertyID (ObjProperty%ID, ObjEnterData, FromBlock)

        call GetData    (ObjProperty%VariableOverLandConc, ObjEnterData, flag,              &
                         keyword      = 'VARIABLE_OVERAND_CONC',                            &
                         ClientModule = 'MohidSwatLink',                                    &
                         SearchType   = FromBlock,                                          &
                         default      = .false.,                                            &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDischargeProperty - MohidSwatLink - ERR01'

        if (ObjProperty%VariableOverLandConc) then

            call GetData    (ObjProperty%SwatSeriesColumn, ObjEnterData, flag,              &
                             keyword      = 'SWAT_SERIES_COLUMN',                           &
                             ClientModule = 'MohidSwatLink',                                &
                             SearchType   = FromBlock,                                      &
                             STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructDischargeProperty - MohidSwatLink - ERR02'

            call GetData    (ObjProperty%SwatToMohidFactor, ObjEnterData, flag,             &
                             keyword      = 'SWAT_TO_MOHID_FACTOR',                         &
                             ClientModule = 'MohidSwatLink',                                &
                             SearchType   = FromBlock,                                      &
                             default      = 1.0,                                            &
                             STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructDischargeProperty - MohidSwatLink - ERR02a'

        else
            
            call GetData    (ObjProperty%OLConcentration, ObjEnterData, flag,               &
                             keyword      = 'OVERLAND_CONC',                                &
                             ClientModule = 'MohidSwatLink',                                &
                             SearchType   = FromBlock,                                      &
                             STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructDischargeProperty - MohidSwatLink - ERR03'

        endif

        call GetData    (ObjProperty%GWConcentration, ObjEnterData, flag,                   &
                         keyword      = 'GROUNDWATER_CONC',                                 &
                         ClientModule = 'MohidSwatLink',                                    &
                         SearchType   = FromBlock,                                          &
                         STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDischargeProperty - MohidSwatLink - ERR04'

    end subroutine ConstructDischargeProperty

end program MohidSwatLink

!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2006. Instituto Superior Técnico, Technical University of Lisbon. 

!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : LagHDF5_2_TimeSeries
! PROGRAM       : MainLagHDF5_2_TimeSeries
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Paulo Leitao
! DESCRIPTION   : LagHDF5_2_TimeSeries to create main program to use MOHID modules
!
!------------------------------------------------------------------------------

program MohidLagHDF5_2_TimeSeries

    use ModuleGlobalData
    use ModuleTime      
    use ModuleStopWatch 
    use ModuleEnterData    
    use ModuleFunctions    
    use ModuleHDF5    
    use ModuleDrawing     


    implicit none
    
    character(LEN = StringLength), parameter    :: Poly_TS_begin              = '<BeginPoly_TS>'
    character(LEN = StringLength), parameter    :: Poly_TS_end                = '<EndPoly_TS>'    
    character(LEN = StringLength), parameter    :: individual_Poly_TS_begin    = '<<BeginIndividualPoly_TS>>'
    character(LEN = StringLength), parameter    :: individual_Poly_TS_end      = '<<EndIndividualPoly_TS>>'    
    
    !Input / Output
    integer, parameter :: FileOpen = 1, FileClose = 0    
    
    type T_IndividualPoly_TS
        type (T_Polygon),        pointer        :: Polygon              => null()
        real                                    :: AgeMin               = null_real
        real                                    :: AgeMax               = null_real
        real                                    :: Factor               = null_real
        character(Len=StringLength)             :: Name                 = null_str
        character(Len=StringLength)             :: Description          = null_str
        character(Len=PathLength)               :: TS_Out_Filename      = null_str
        integer                                 :: TS_Out_Obj           = null_int
    end type T_IndividualPoly_TS
    
    type T_Poly_TS
        integer                                 :: Number               = null_int
        type (T_IndividualPoly_TS), dimension(:), pointer :: Individual    => null()
        integer                                 :: ObjEnterData         = 0
        character(len=PathLength)               :: FileNameHDF5In   
        integer                                 :: ObjHDF5              = 0
        integer                                 :: iTotal               = 0
        character(PathLength)                   :: GroupX, GroupY, GroupP
        character(PathLength)                   :: FieldX, FieldY, FieldP    
        real, dimension(:), pointer             :: X, Y, P
        integer                                 :: Size1D
    end type T_Poly_TS    

    type(T_Time)                :: BeginTime, SecondTime, EndTime
    real                        :: DT
    logical                     :: VariableDT
    type (T_Time)               :: InitialSystemTime, FinalSystemTime
    real                        :: TotalCPUTime, ElapsedSeconds
    integer, dimension(8)       :: F95Time

    integer                     :: ObjTime = 0
    integer                     :: STAT_CALL
    
    type(T_Poly_TS), pointer    :: Me


    allocate(Me)
    
    call ConstructMohidLagHDF5_2_TimeSeries
    call ModifyMohidLagHDF5_2_TimeSeries
    call KillMohidLagHDF5_2_TimeSeries
    
    deallocate(Me)

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidLagHDF5_2_TimeSeries
    
        !Variables------------------------------------------------
        character(PathLength)                       :: DataFile    
        
        call StartUpMohid("MohidLagHDF5_2_TimeSeries")

        call StartCPUTime
        
        call ReadFileName('IN_MODEL', DataFile, "MohidLagHDF5_2_TimeSeries", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR01'        
        
        call ConstructEnterData (Me%ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidLagHDF5_2_TimeSeries - MohidLagHDF5_2_TimeSeries - ERR10'

        call ReadKeywords
        
        call ReadTime

        call StartComputeTime(ObjTime, BeginTime, BeginTime, EndTime, DT = DT, VariableDT = VariableDT, STAT = STAT_CALL)
                
        call ReadPoly_TS
        
        call KillEnterData (Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidLagHDF5_2_TimeSeries - MohidLagHDF5_2_TimeSeries - ERR20'
        

    end subroutine ConstructMohidLagHDF5_2_TimeSeries
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidLagHDF5_2_TimeSeries
        
        !Local-----------------------------------------------------------------
        integer                                     :: i
        type(T_Time)                                :: CurrentTime


        do i=1, Me%iTotal
            
            CurrentTime = HDF5TimeInstant(i, Me%ObjHDF5)
            
            call Read_HDF5_XYP(i)
            
            call WritePolyTimeSeries(CurrentTime)
            


            !call PrintProgress(ObjTime, STAT = STAT_CALL)
            !if(STAT_CALL .ne. SUCCESS_)stop 'ModifyMohidLagHDF5_2_TimeSeries - ERR02'

        enddo
    
    
    
    end subroutine ModifyMohidLagHDF5_2_TimeSeries
    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine Read_HDF5_XYP(i)
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, imax, STAT_CALL, NDim
        
        
        
        call GetHDF5ArrayDimensions(HDF5ID      = Me%ObjHDF5,                           & 
                                    GroupName   = Me%GroupX,                            &
                                    ItemName    = Me%FieldX,                            &
                                    OutputNumber= i,                                    &
                                    NDim        = NDim,                                 &
                                    Imax        = imax,                                 &
                                    STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_XYP - MohidLagHDF5_2_TimeSeries - ERR10'

        call HDF5SetLimits  (Me%ObjHDF5, 1, imax, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_XYP - MohidLagHDF5_2_TimeSeries - ERR20'        
        
        allocate(Me%X(1:imax))
        allocate(Me%Y(1:imax))        
        allocate(Me%P(1:imax))
        
        Me%Size1D = imax

        call HDF5ReadWindow (HDF5ID         = Me%ObjHDF5,                               &
                             GroupName      = Me%GroupX,                                &
                             Name           = Me%FieldX,                                &
                             Array1D        = Me%X,                                     &
                             OutputNumber   = i,                                        &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_XYP - MohidLagHDF5_2_TimeSeries - ERR30'
    
        call HDF5ReadWindow (HDF5ID         = Me%ObjHDF5,                               &
                             GroupName      = Me%GroupY,                                &
                             Name           = Me%FieldY,                                &
                             Array1D        = Me%Y,                                     &
                             OutputNumber   = i,                                        &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_XYP - MohidLagHDF5_2_TimeSeries - ERR40'
    
        call HDF5ReadWindow (HDF5ID         = Me%ObjHDF5,                               &
                             GroupName      = Me%GroupP,                                &
                             Name           = Me%FieldP,                                &
                             Array1D        = Me%P,                                     &
                             OutputNumber   = i,                                        &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Read_HDF5_XYP - MohidLagHDF5_2_TimeSeries - ERR50'
        
    
    end subroutine Read_HDF5_XYP
    
    !--------------------------------------------------------------------------    
    !--------------------------------------------------------------------------

    subroutine WritePolyTimeSeries(CurrentTime)
        
        !Local-----------------------------------------------------------------
        type(T_Time)                                :: CurrentTime
        integer                                     :: n, nPoly_TS, np, STAT_CALL, Unit
        type (T_PointF)   , pointer                 :: Point
        type (T_Polygon)  , pointer                 :: Polygons
        real                                        :: Days, Index, Sum, Average
        character(len=256)                          :: String256
        
        
        allocate (Point) 
        
DONB:   do nPoly_TS = 1, Me%Number
    
            np = 0
            Sum = 0
    
            do n=1, Me%Size1D
                
                Point%X = Me%X(n)
                Point%Y = Me%Y(n)
                
                Polygons => Me%Individual(nPoly_TS)%Polygon
                
                if (IsVisible(Polygons, Point)) then
                    
                    if (Me%P(n) <  Me%Individual(nPoly_TS)%AgeMax .and.                 &
                        Me%P(n) >= Me%Individual(nPoly_TS)%AgeMin) then
                        
                        np = np + 1
                        
                        Sum = Sum + Me%P(n)
                        
                    endif
                    
                endif

            
            enddo
            
            Unit = Me%Individual(nPoly_TS)%TS_Out_Obj
             
            
            if (CurrentTime == BeginTime) then

                call WriteDataLine(Unit, 'NAME'              , 'MOHID Parcels' )
                call WriteDataLine(Unit, 'LOCALIZATION_I'    , '-999999')
                call WriteDataLine(Unit, 'LOCALIZATION_J'    , '-999999')
                call WriteDataLine(Unit, 'LOCALIZATION_K'    , '-999999')
                call WriteDataLine(Unit, 'SERIE_INITIAL_DATA', CurrentTime)
                call WriteDataLine(Unit, 'TIME_UNITS        ', 'DAYS')

                write (String256,'(A)')  'Days Number Index Average'
                write (Unit,'(A)') trim(adjustl(String256))
                write (Unit,*)  '<BeginTimeSerie>'                            
            
            endif
            
            Days = (CurrentTime - BeginTime) / 86400.
            
            Index = real(np) * Me%Individual(nPoly_TS)%Factor
            
            if (np > 0) then
                Average = Sum /  real(np)
            else
                Average = 0.
            endif
            
            
            write(Unit, fmt=100) Days, np, Index, Average

100         format(1x, f16.2, 1x, i16, 1x, f16.2, 1x, f16.2)   
            if (CurrentTime == EndTime) then            
                write (Unit,*)  '<EndTimeSerie>' 
                
                call UnitsManager(Unit, FileClose, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'WritePolyTimeSeries - MohidLagHDF5_2_TimeSeries - ERR10'                   
            endif
                                
        enddo DONB            
        
        
        deallocate(Point) 
        
        deallocate(Me%X)
        deallocate(Me%Y)        
        deallocate(Me%P)

    
    end subroutine WritePolyTimeSeries
    
    !--------------------------------------------------------------------------    

    subroutine KillMohidLagHDF5_2_TimeSeries
    
        !Local--------------------------------------
        integer                 :: nPoly_TS, STAT_CALL
    
        !Kill HDF5 File
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'KillMohidLagHDF5_2_TimeSeries - MohidLagHDF5_2_TimeSeries - ERR10'    
        endif
        
DONB:   do nPoly_TS = 1, Me%Number
                
           deallocate(Me%Individual(nPoly_TS)%Polygon)
            
        enddo DONB        

        deallocate(Me%Individual)

        call StopCPUTime

        call ShutdownMohid ("MohidLagHDF5_2_TimeSeries", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidLagHDF5_2_TimeSeries
    
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
        integer                                     :: FromFile, flag
        integer                                     :: HDF5_READ


        call GetExtractType     (FromFile = FromFile)

        !call ReadTimeKeyWords   (Me%ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
        !                         VariableDT, "MohidLagHDF5_2_TimeSeries")

        call GetData(Me%FileNameHDF5In,                                                     &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='FILENAME_LAG_HDF',                                      &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR03'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR04'    
        
        call GetData(Me%GroupX,                                                             &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='GROUP_HDF_X',                                           &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR05'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR06'  
        
        
        call GetData(Me%GroupY,                                                             &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='GROUP_HDF_Y',                                           &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR07'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR08'  
        
        
        call GetData(Me%GroupP,                                                             &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='GROUP_HDF_P',                                           &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR09'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR10'  

        call GetData(Me%FieldX,                                                             &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='FIELD_HDF_X',                                           &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR05'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR06'  
        
        
        call GetData(Me%FieldY,                                                             &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='FIELD_HDF_Y',                                           &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR07'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR08'  
        
        
        call GetData(Me%FieldP,                                                             &
                     Me%ObjEnterData,                                                       &
                     flag,                                                                  &
                     SearchType   = FromFile,                                               &
                     keyword      ='FIELD_HDF_P',                                           &
                     ClientModule ='MohidLagHDF5_2_TimeSeries',                             &
                     STAT         = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR09'
        if (flag==0)               stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR10'         
        
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%FileNameHDF5In, HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidLagHDF5_2_TimeSeries - ERR11'    
        



    end subroutine ReadKeywords
    
    !--------------------------------------------------------------------------

    subroutine ReadTime

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iTotal


        call GetHDF5GroupNumberOfItems(Me%ObjHDF5, "/Time", iTotal, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'ReadTime - MohidLagHDF5_2_TimeSeries - ERR60'
        endif
        
        BeginTime = HDF5TimeInstant(1     , Me%ObjHDF5)        
        DT        = HDF5TimeInstant(2     , Me%ObjHDF5) - BeginTime
        EndTime   = HDF5TimeInstant(iTotal, Me%ObjHDF5)        
        
        Me%iTotal = iTotal


    end subroutine ReadTime
    
    !--------------------------------------------------------------------------

    subroutine ReadPoly_TS

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: Poly_TSFound
        integer                                     :: STAT_CALL, ClientNumber

        !----------------------------------------------------------------------


        call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,                      &   
                                   Poly_TS_begin, Poly_TS_end,                  &
                                   Poly_TSFound, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadPoly_TS - MohidLagHDF5_2_TimeSeries - ERR10'

if1:    if (Poly_TSFound) then

            call CountIndividualPoly_TS(ClientNumber)

            call ReadIndividualPoly_TS (ClientNumber)
            

        endif if1

        call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadPoly_TS - MohidLagHDF5_2_TimeSeries - ERR20'
        

        !Prepares file for a new block search throughout the entire file
        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadPoly_TS - MohidLagHDF5_2_TimeSeries - ERR30'
        

    end subroutine ReadPoly_TS

    !--------------------------------------------------------------------------

    subroutine CountIndividualPoly_TS(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_CALL, NPoly_TS
        logical                                         :: Poly_TSFound

        !Begin-----------------------------------------------------------------

        NPoly_TS = 0

DOPROP: do 
            
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                   &
                                       individual_Poly_TS_begin, individual_Poly_TS_end,&
                                       Poly_TSFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CountIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR10'
            
i1:         if (Poly_TSFound) then

                NPoly_TS = NPoly_TS + 1
 
            else i1
            
                call RewindBlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)  
                if (STAT_CALL /= SUCCESS_) stop 'CountIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR20'
                exit
            endif i1

        enddo DOPROP
        
        Me%Number = NPoly_TS
        
        if (Me%Number<1) stop 'CountIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR30'                

        allocate(Me%Individual(Me%Number))

    end subroutine CountIndividualPoly_TS
    !--------------------------------------------------------------------------

    subroutine ReadIndividualPoly_TS(ClientNumber)

        !Arguments-------------------------------------------------------------
        integer                                         :: ClientNumber

        !Local-----------------------------------------------------------------
        character (Len = PathLength)                    :: PolyFileName
        logical                                         :: Poly_TSFound
        integer                                         :: STAT_CALL, nPoly_TS, flag
        !Begin-----------------------------------------------------------------
        
        
        
DONB:   do nPoly_TS = 1, Me%Number
 
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                   &
                                       individual_Poly_TS_begin, individual_Poly_TS_end,&
                                       Poly_TSFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR10'
            
i1:         if (Poly_TSFound) then

                call GetData(Me%Individual(nPoly_TS)%Name,                      &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='NAME',                                      &
                             default      = 'Poly_TS X',                                &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR20'
 
                call GetData(Me%Individual(nPoly_TS)%Description,               &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='DESCRIPTION',                               &
                             default      = 'Poly_TS to contain a oil spill',           &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR30'
 
                call GetData(Me%Individual(nPoly_TS)%AgeMax,                            &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                            !Age in days
                             keyword      ='AGE_MAX',                                   &
                             default      = 1e6,                                        &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR40'

                call GetData(Me%Individual(nPoly_TS)%AgeMin,                            &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                            !Age in days
                             keyword      ='AGE_MIN',                                   &
                             default      = 0.,                                         &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR50'
                
                call GetData(Me%Individual(nPoly_TS)%Factor,                            &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                            !factor
                             keyword      ='FACTOR',                                    &
                             default      = 1.,                                         &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR50'                
                
                call GetData(Me%Individual(nPoly_TS)%TS_Out_Filename,                   &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='TS_OUT',                                    &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR55'   
                if (flag == 0            ) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR57'  
                
                call UnitsManager(Me%Individual(nPoly_TS)%TS_Out_Obj, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR59'   

                open(unit   = Me%Individual(nPoly_TS)%TS_Out_Obj,                       &
                     file   = Me%Individual(nPoly_TS)%TS_Out_Filename,                  &
                     form   = 'FORMATTED',                                              &
                     status = 'UNKNOWN')
                                
                call GetData(PolyFileName,                                               &
                             Me%ObjEnterData,                                           &
                             flag,                                                      &
                             SearchType   = FromBlockInBlock,                           &
                             keyword      ='FILENAME',                                  &
                             ClientModule ='MohidLagHDF5_2_TimeSeries',                 &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR60'
                
                if (flag == 0) then
                    write(*,*) 'Poly_TS named ', trim(Me%Individual(nPoly_TS)%Name),' needs a ".xy" file'
                    stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR70'                
                endif
                
                nullify(Me%Individual(nPoly_TS)%Polygon)
                
                call New(Me%Individual(nPoly_TS)%Polygon, PolyFileName)                                
                
              
                
            else i1
            
                stop 'ReadIndividualPoly_TS - MohidLagHDF5_2_TimeSeries - ERR90'

            endif i1
            
        enddo DONB
        


    end subroutine ReadIndividualPoly_TS
                
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    type(T_Time) function HDF5TimeInstant(Instant,ObjHDF5)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        integer                                 :: ObjHDF5
        !External--------------------------------------------------------------
        integer                                 :: STAT_CALL
        !Local-----------------------------------------------------------------
        real, dimension(:), pointer             :: TimeVector
        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID     = ObjHDF5,               &
                             GroupName  = "/Time", Name = "Time",                   &
                             Array1D    = TimeVector, OutputNumber   = Instant, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleUpscaleHDF5 - ERR10'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))
        deallocate(TimeVector)
        nullify   (TimeVector)

    end function HDF5TimeInstant

    !--------------------------------------------------------------------------    

end program MohidLagHDF5_2_TimeSeries

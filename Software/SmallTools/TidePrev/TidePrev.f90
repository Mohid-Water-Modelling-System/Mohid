!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Tide Preview
! PROJECT       : TidePreview
! PROGRAM       : MainTidePreview
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : August 2004
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Program to output tidal elevations preview based on gauge files
!
!------------------------------------------------------------------------------

!Data file - default name 'TidePrevInput.dat' (must be placed in working directory)

!   START                       : YYYY MM DD HH MM SS -         !Start time to compute water level
!   END                         : YYYY MM DD HH MM SS -         !End time to compute water level
!   DT                          : real              -           !Time step to compute water level
!   EXPORT_TO_XYZ               : 0/1               0           !Create a XYZ file with gauge locations
!   XYZ_FILE                    : char              -           !Name of XYZ file to be created

!<begintideprev>
!   IN_TIDES                    : char              -           !Path to gauge file
!   OUT_FILE                    : char              -           !Path to output water level time serie file 
!<endtideprev>

program MohidTidePreview

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleGauge
    !use lsq
	USE nrtype; USE nrutil, ONLY : arth,assert,poly
	USE nr, ONLY : lubksb,ludcmp

    implicit none

    !Time variables
    type(T_Time)                        :: BeginTime, EndTime, CurrentTime
    type(T_Time)                        :: InitialSystemTime, FinalSystemTime
    logical                             :: VariableDT
    real                                :: DT, TotalCPUTime, CPUTime, LastCPUTime, ElapsedSeconds
    integer, dimension(8)               :: F95Time

    !Types
    type T_TidePrev
        integer                         :: ObjGauge             = 0
        character(len=PathLength)       :: GaugeFile
        character(len=PathLength)       :: OutPutFileName
        character(len=PathLength)       :: OutFileNameHighLow
        character(len=StringLength)     :: Name
        integer                         :: OutPutFileUnit
        integer                         :: OutPutHighLowUnit        
        real                            :: Longitude, Latitude
        real, dimension(:), pointer     :: ReferenceLevel
        type(T_TidePrev), pointer       :: Next
        real                            :: Aux1,Aux2,Aux3,Aux4,Aux5
        integer                         :: counter
        real                            :: PreviousTime
        logical                         :: SmoothTimeSerie = .true. 
    end type T_TidePrev
    

    integer                             :: ObjTime              = 0
    type(T_TidePrev), pointer           :: FirstTidePrev
    character(len=PathLength)           :: DataFile             = 'TidePrevInput.dat'
    logical                             :: ExportToXYZ          = OFF
    character(len=PathLength)           :: XYZFile


    call ConstructMohidTidePreview
    call ModifyMohidTidePreview
    call KillMohidTidePreview

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidTidePreview

        !Local-----------------------------------------------------------------
        integer                                     :: NGauges, TotalTidePrevs
        integer                                     :: ObjEnterData         = 0
        integer                                     :: iflag, STAT_CALL, ClientNumber
        logical                                     :: BlockFound
        type(T_TidePrev),   pointer                 :: NewTidePrev
        real, dimension(:), pointer                 :: XLocation, YLocation

        !Begin-----------------------------------------------------------------

        nullify(NewTidePrev, FirstTidePrev)
        
        call StartUpMohid("MohidTidePreview")

        call StartCPUTime
        
        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR10'


        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,        &
                                 VariableDT, "MohidTidePreview")
        
        call StartComputeTime(ObjTime, InitialSystemTime, BeginTime, EndTime, DT = DT,  &
                              VariableDT = VariableDT, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR20'

        call GetData(ExportToXYZ,                                                       &
                     ObjEnterData,  iflag,                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'EXPORT_TO_XYZ',                                    &
                     Default      = OFF,                                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR30'


        call GetData(XYZFile,                                                           &
                     ObjEnterData,  iflag,                                              &
                     SearchType   = FromFile,                                           &
                     keyword      = 'XYZ_FILE',                                         &
                     Default      = 'GaugeLocation.xyz',                                &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR40'


        TotalTidePrevs = 0

        write(*,*)
        write(*,*)'Reading gauge(s) information'
        
        do 
            call ExtractBlockFromBuffer(ObjEnterData, ClientNumber,                         &
                                        '<begintideprev>', '<endtideprev>', BlockFound,     &
                                        STAT = STAT_CALL)
            if (STAT_CALL .EQ. SUCCESS_) then
                
                if (BlockFound) then

                    call AddTidePrev (FirstTidePrev, NewTidePrev)
        

                    call GetData(NewTidePrev%Name,                                          &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'NAME',                                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR50'

                    call GetData(NewTidePrev%OutPutFileName,                                &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'OUT_FILE',                                 &
                                 default      = 'TidePrev.srh',                             &                                 
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR60'

                    call GetData(NewTidePrev%OutFileNameHighLow,                            &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'OUT_FILE_HIGH_LOW',                        &
                                 default      = 'TidePrevHighLowTide.srh',                  &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR70'
                    

                    call GetData(NewTidePrev%SmoothTimeSerie,                               &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'SMOOTH_TIME_SERIE',                        &
                                 default      = .true.,                                     &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR75'

                    call GetData(NewTidePrev%GaugeFile,                                     &
                                 ObjEnterData,  iflag,                                      &
                                 SearchType   = FromBlock,                                  &
                                 keyword      = 'IN_TIDES',                                 &
                                 STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR80'
                    if (iflag     == 0       ) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR90'


                    call ConstructGauges(GaugeID    = NewTidePrev%ObjGauge,                 &
                                         TimeID     = ObjTime,                              &
                                         GaugeFile  = NewTidePrev%GaugeFile,                &
                                         STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR100'


                    !Get the number of gauges in use
                    call GetNGauges(NewTidePrev%ObjGauge, NGauges, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR110'
        
                    if(NGauges .ne. 1)then
                        write(*,*)'Can only compute one gauge per tide preview at the time.'
                        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR120'
                    end if

                    allocate(NewTidePrev%ReferenceLevel(NGauges))

                    !Get the current elevation at the gauges
                    call GetReferenceLevel(NewTidePrev%ObjGauge,                            &
                                           NewTidePrev%ReferenceLevel,                      &
                                           STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR130'

                    allocate(XLocation(NGauges), YLocation(NGauges))

                    call GetGaugeLocation(NewTidePrev%ObjGauge, XLocation, YLocation, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR140'

                    NewTidePrev%Longitude = XLocation(1)
                    NewTidePrev%Latitude  = YLocation(1)
                    
                    TotalTidePrevs = TotalTidePrevs + 1
                else
                    exit     !No more blocks
                end if 
            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then 
                stop 'ConstructMohidTidePreview - MohidTidePreview - ERR150'
            end if
        end do

        if(ExportToXYZ) call ExportGaugeLocations

        call Block_Unlock(ObjEnterData, ClientNumber, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR160'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructMohidTidePreview - MohidTidePreview - ERR170'

        write(*,*)
        write(*,*)'Total number of gauges to preview :', TotalTidePrevs

    end subroutine ConstructMohidTidePreview
    
    !--------------------------------------------------------------------------

    subroutine ExportGaugeLocations
        
        !Local-----------------------------------------------------------------
        type(T_TidePrev), pointer                   :: TidePrev
        integer                                     :: STAT_CALL, XYZUnit
        integer                                     :: ID

        !Begin-----------------------------------------------------------------

        ID = 1

        call UnitsManager (XYZUnit, OPEN_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExportGaugeLocations - MohidTidePreview - ERR01'

        open(XYZUnit, file = trim(XYZFile), Status = 'unknown')

        write(XYZUnit, *)"<begin_xyz>"

        TidePrev => FirstTidePrev
        
        do while(associated(TidePrev))

            write(XYZUnit, *) TidePrev%Longitude, TidePrev%Latitude, ID, trim(TidePrev%Name)

            ID = ID + 1

            TidePrev => TidePrev%Next
        end do

        write(XYZUnit, *)"<end_xyz>"

        call UnitsManager (XYZUnit, CLOSE_FILE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ExportGaugeLocations - MohidTidePreview - ERR02'


    end subroutine ExportGaugeLocations
   
    !--------------------------------------------------------------------------

    subroutine ModifyMohidTidePreview
        
        !Local-----------------------------------------------------------------
        logical                                     :: Running
        integer                                     :: STAT_CALL
        real                                        :: TotalTime, PreviousTime, AuxDT
        real,       dimension(:), pointer           :: OpenPoints, WaterLevel
        real                                        :: WL
        type(T_TidePrev), pointer                   :: TidePrev
        real                                        :: DTAux, WL_out

        !Begin-----------------------------------------------------------------


        Running      = .true.
        CurrentTime  = BeginTime

        allocate(OpenPoints(1), WaterLevel(1))
        OpenPoints = 1

        
        TidePrev => FirstTidePrev
        
        do while(associated(TidePrev))

            call UnitsManager (TidePrev%OutPutFileUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR10'

            open(TidePrev%OutPutFileUnit, file = TidePrev%OutPutFileName, Status = 'unknown')

            call WriteDataLine(TidePrev%OutPutFileUnit, "SERIE_INITIAL_DATA", CurrentTime)
            write(TidePrev%OutPutFileUnit, *) "TIME_UNITS              : SECONDS"
            write(TidePrev%OutPutFileUnit, *) "LONGITUDE               : ", TidePrev%Longitude
            write(TidePrev%OutPutFileUnit, *) "LATITUDE                : ", TidePrev%Latitude

            write(TidePrev%OutPutFileUnit, *)
            write(TidePrev%OutPutFileUnit, *)
            write(TidePrev%OutPutFileUnit, *) 'Seconds   Elevation'
            write(TidePrev%OutPutFileUnit, *) '<BeginTimeSerie>'
            
            call UnitsManager (TidePrev%OutPutHighLowUnit, OPEN_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR20'

            open(TidePrev%OutPutHighLowUnit, file = TidePrev%OutFileNameHighLow, Status = 'unknown')

            call WriteDataLine(TidePrev%OutPutHighLowUnit, "SERIE_INITIAL_DATA", CurrentTime)
            write(TidePrev%OutPutHighLowUnit, *) "TIME_UNITS              : SECONDS"
            write(TidePrev%OutPutHighLowUnit, *) "LONGITUDE               : ", TidePrev%Longitude
            write(TidePrev%OutPutHighLowUnit, *) "LATITUDE                : ", TidePrev%Latitude

            write(TidePrev%OutPutHighLowUnit, *)
            write(TidePrev%OutPutHighLowUnit, *)
            write(TidePrev%OutPutHighLowUnit, *) 'Seconds   Elevation'
            write(TidePrev%OutPutHighLowUnit, *) '<BeginTimeSerie>'

            TidePrev%counter = 0
            
            TidePrev%PreviousTime  = 0.

            TidePrev => TidePrev%Next
        end do


        TotalTime      = 0.
        PreviousTime   = 0.

        do while (Running)

            TidePrev => FirstTidePrev
            
            do while(associated(TidePrev))

                call GaugeLevel(TidePrev%ObjGauge,                          &
                                WaterLevel,                                 &
                                OpenPoints,                                 &
                                CurrentTime,                                &
                                ReferenceLevel = TidePrev%ReferenceLevel,   &
                                STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR30'
                
                WL = WaterLevel(1) + TidePrev%ReferenceLevel(1)


                if (TidePrev%SmoothTimeSerie) then

                    call FindLowAndHighTideSmoothTS(TidePrev, WL, PreviousTime)
                
                else
                
                    call FindLowAndHighTideNoisyTS(TidePrev, WL, PreviousTime, TotalTime, WL_Out)  
                    
                    WL = WL_Out
                
                endif
                
                write(TidePrev%OutPutFileUnit, *) TotalTime, WL                
                
                TidePrev => TidePrev%Next

            end do

            call cpu_time(CPUTime)
            !if (CPUTime - LastCPUTime > 10.) then
            !    LastCPUTime = CPUTime
            !    call PrintProgress(ObjTime, STAT = STAT_CALL)
            !    if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR40'
            !endif

            
            CurrentTime     = CurrentTime + DT
            PreviousTime    = TotalTime
            TotalTime       = TotalTime   + DT

            call ActualizeCurrentTime(ObjTime, DT, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR50'


            if (abs(CurrentTime - EndTime) > DT / 10.) then
                Running = .true.
            else
                Running = .false.
            endif

        enddo
    



        TidePrev => FirstTidePrev

        do while(associated(TidePrev))

            write(TidePrev%OutPutFileUnit, *) '<EndTimeSerie>'

            call UnitsManager (TidePrev%OutPutFileUnit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR60'
            
            write(TidePrev%OutPutHighLowUnit, *) '<EndTimeSerie>'

            call UnitsManager (TidePrev%OutPutHighLowUnit, CLOSE_FILE, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyMohidTidePreview - MohidTidePreview - ERR70'
            

            TidePrev => TidePrev%Next

        end do

     
    end subroutine ModifyMohidTidePreview
    
    !--------------------------------------------------------------------------
    
    subroutine FindLowAndHighTideSmoothTS(TidePrev, WaterLevel, PreviousTime)

        !Arguments-------------------------------------------------------------
        real                                        :: PreviousTime
        real                                        :: WaterLevel
        type(T_TidePrev), pointer                   :: TidePrev
        
        !Local-----------------------------------------------------------------
        real                                        :: DTAux

        !Begin-----------------------------------------------------------------

        TidePrev%counter    = TidePrev%counter + 1
        TidePrev%Aux3       = TidePrev%Aux2
        TidePrev%Aux2       = TidePrev%Aux1
        TidePrev%Aux1       = WaterLevel
        
        DTAux = PreviousTime - TidePrev%PreviousTime
        
        if (TidePrev%counter > 3 .and. DTAux > 10800.) then
            !low tide 
            if ((TidePrev%Aux2 <= TidePrev%Aux1 .and. TidePrev%Aux2 <= TidePrev%Aux3) .or. &
            !high tide 
                (TidePrev%Aux2 >= TidePrev%Aux1 .and. TidePrev%Aux2 >= TidePrev%Aux3)) then
                write(TidePrev%OutPutHighLowUnit, *) PreviousTime, TidePrev%Aux2
                TidePrev%PreviousTime = PreviousTime
            endif
        endif
    
    end subroutine FindLowAndHighTideSmoothTS                    
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    subroutine FindLowAndHighTideNoisyTS(TidePrev, WaterLevel, PreviousTime, TotalTime, WL_Out)

        !Arguments-------------------------------------------------------------
        real                                        :: PreviousTime, TotalTime
        real                                        :: WaterLevel
        type(T_TidePrev), pointer                   :: TidePrev
        real                                        :: WL_Out
        
        !Local-----------------------------------------------------------------
        integer,  parameter                         :: NS         = 1000
        integer,  parameter                         :: Nx         = 20
        real,     parameter                         :: DTanalysis = 10800
        real,     dimension(1:NS)                   :: X    = FillValueReal
        real,     dimension(1:NS)                   :: Y    = FillValueReal         
        real,     dimension(1:NS)                   :: X1   = FillValueReal
        real,     dimension(1:NS)                   :: Y1   = FillValueReal
        REAL(SP), dimension(1:NS)                   :: C    = FillValueReal        
	    INTEGER(I4B)                                :: nl,nrr,ld,m
        real                                        :: DTAux
        integer                                     :: i, iS, NxX, ii
        
        !Begin-----------------------------------------------------------------

        X(2:NS) = X(1:NS-1) 
        Y(2:NS) = Y(1:NS-1) 
        X(1)     = TotalTime
        Y(1)     = WaterLevel
        
        do i=2,NS
            if (X(1) - X(i) > DTanalysis) then
                iS = i-1
                exit
            else
                iS = NS
            endif
        enddo
        
        TidePrev%counter    = TidePrev%counter + 1
        TidePrev%Aux5       = TidePrev%Aux4
        TidePrev%Aux4       = TidePrev%Aux3
        TidePrev%Aux3       = TidePrev%Aux2
        TidePrev%Aux2       = TidePrev%Aux1
        
        if (TidePrev%counter > 3) then
            !call fit_poly(3,X(1:iS),Y(1:iS),iS,X(1),TidePrev%Aux1) 
             TidePrev%Aux1 = 0
             nl=0;nrr=iS-1;ld = 0;m = 2
             c(1:iS) = savgol(nl=nl,nrr=nrr,ld = ld,m = m)
             do i=1, iS
                TidePrev%Aux1 = TidePrev%Aux1 + C(i) * Y(i)
             enddo
        else
        
            TidePrev%Aux1 = WaterLevel
        
        endif    
        
        WL_Out = TidePrev%Aux1
        
        DTAux = PreviousTime - TidePrev%PreviousTime

        if (TidePrev%counter > 3 .and. DTAux > DTanalysis) then
            !low tide 
            if ((TidePrev%Aux2 <= TidePrev%Aux1 .and. TidePrev%Aux2 <= TidePrev%Aux3 .and. &
                 TidePrev%Aux2 <= TidePrev%Aux4 .and. TidePrev%Aux2 <= TidePrev%Aux5) .or. &
            !high tide 
                (TidePrev%Aux2 >= TidePrev%Aux1 .and. TidePrev%Aux2 >= TidePrev%Aux3 .and. &
                 TidePrev%Aux2 >= TidePrev%Aux4 .and. TidePrev%Aux2 >= TidePrev%Aux5)) then
                write(TidePrev%OutPutHighLowUnit, *) PreviousTime, TidePrev%Aux2
                TidePrev%PreviousTime = PreviousTime
            endif
        endif
    
    end subroutine FindLowAndHighTideNoisyTS                    
    !--------------------------------------------------------------------------    
    
	FUNCTION savgol(nl,nrr,ld,m)
	INTEGER(I4B), INTENT(IN) :: nl,nrr,ld,m
	REAL(SP), DIMENSION(nl+nrr+1) :: savgol
	INTEGER(I4B) :: imj,ipj,mm,np
	INTEGER(I4B), DIMENSION(m+1) :: indx
	REAL(SP) :: d,sm
	REAL(SP), DIMENSION(m+1) :: b
	REAL(SP), DIMENSION(m+1,m+1) :: a
	INTEGER(I4B) :: irng(nl+nrr+1)
	call assert(nl >= 0, nrr >= 0, ld <= m, nl+nrr >= m, 'savgol args')
	do ipj=0,2*m
		sm=sum(arth(1.0_sp,1.0_sp,nrr)**ipj)+&
			sum(arth(-1.0_sp,-1.0_sp,nl)**ipj)
		if (ipj == 0) sm=sm+1.0_sp
		mm=min(ipj,2*m-ipj)
		do imj=-mm,mm,2
			a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sm
		end do
	end do
	call ludcmp(a(:,:),indx(:),d)
	b(:)=0.0
	b(ld+1)=1.0
	call lubksb(a(:,:),indx(:),b(:))
	savgol(:)=0.0
	irng(:)=arth(-nl,1,nrr+nl+1)
	np=nl+nrr+1
	savgol(mod(np-irng(:),np)+1)=poly(real(irng(:),sp),b(:))
	END FUNCTION savgol
    
    !-----------------------------------------------------------------------------------------    
    !subroutine fit_poly(m, ArrayX, ArrayY, n, X, Y)
    ! A simple program to fit a polynomial in one variable.
    ! Data must be store in the form of pairs, either (x,y) or (y,x)
    ! Polynomial to be fitted:

    ! Y = a(0) + a(1).X + a(2).X^2 + ... + a(m).X^m
    
        !Arguments------------------------------------------------------------------
     !   integer,               intent(IN) :: m, n
     !   real,   dimension(:),  intent(IN) :: ArrayX, ArrayY
     !   real,                  intent(IN) :: X         
     !   real,                  intent(OUT):: Y

        !Local----------------------------------------------------------------------
     !   CHARACTER (LEN=50)                :: fname
     !   CHARACTER (LEN= 1)                :: ans
     !   REAL (dp), dimension(0:20)        :: xrow, beta, sterr
     !   REAL (dp), dimension(231)         :: covmat
     !   REAL (dp)                         :: wt = 1.0_dp, var, totalSS, Yaux
     !   INTEGER                           :: i, ier, iostatus, j
     !   LOGICAL,   dimension(0:20)        :: lindep
     !   LOGICAL                           :: fit_const = .TRUE., xfirst

        !Begin----------------------------------------------------------------------

        ! Least-squares calculations
        !CALL startup(m, fit_const)
        !DO i = 1, n
        !  xrow(0) = 1.0_dp
        !  DO j = 1, m
        !    xrow(j) = ArrayX(i) * xrow(j-1)
        !  END DO
        !  Yaux = ArrayY(i)
        !  CALL includ(wt, xrow, Yaux)
        !END DO

       ! CALL sing(lindep, ier)
       ! IF (ier /= 0) THEN
       !   DO i = 0, m
       !     IF (lindep(i)) WRITE(*, '(a, i3)') ' Singularity detected for power: ', i
       !     IF (lindep(i)) WRITE(9, '(a, i3)') ' Singularity detected for power: ', i
       !   END DO
       ! END IF

        ! Calculate progressive residual sums of squares
      !  CALL ss()
      !  var = rss(m+1) / (n - m - 1)

        ! Calculate least-squares regn. coeffs.
      !  CALL regcf(beta, m+1, ier)

        ! Calculate covariance matrix, and hence std. errors of coeffs.
      !  CALL cov(m+1, var, covmat, 231, sterr, ier)

        !WRITE(*, *) 'Least-squares coefficients & std. errors'
        !WRITE(9, *) 'Least-squares coefficients & std. errors'
        !WRITE(*, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'
        !WRITE(9, *) 'Power  Coefficient          Std.error      Resid.sum of sq.'
        !DO i = 0, m
        !  WRITE(*, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &
        !        i, beta(i), sterr(i), rss(i+1)
        !  WRITE(9, '(i4, g20.12, "   ", g14.6, "   ", g14.6)')  &
        !        i, beta(i), sterr(i), rss(i+1)
        !END DO

        !WRITE(*, *)
        !WRITE(9, *)
        !WRITE(*, '(a, g20.12)') ' Residual standard deviation = ', SQRT(var)
        !WRITE(9, '(a, g20.12)') ' Residual standard deviation = ', SQRT(var)
        !totalSS = rss(1)
        !WRITE(*, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS
        !WRITE(9, '(a, g20.12)') ' R^2 = ', (totalSS - rss(m+1))/totalSS

     !   Y = 0.
     !   do i = 0, m
     !       Y = Y + beta(i) * X ** real(i)
     !   enddo
        !Y = a(0) + a(1).X + a(2).X^2 + ... + a(m).X^m

    !end subroutine fit_poly
    
    !--------------------------------------------------------------------------    


    subroutine KillMohidTidePreview

        call StopCPUTime

        call ShutdownMohid ("MohidTidePreview", ElapsedSeconds, TotalCPUTime)

    end subroutine KillMohidTidePreview
    
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

    subroutine AddTidePrev (FirstTidePrev, ObjTidePrev)

        !Arguments-------------------------------------------------------------
        type (T_TidePrev), pointer                   :: FirstTidePrev
        type (T_TidePrev), pointer                   :: ObjTidePrev

        !Local-----------------------------------------------------------------
        type (T_TidePrev), pointer                   :: NewTidePrev
        type (T_TidePrev), pointer                   :: PreviousTidePrev
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewTidePrev)
        nullify  (NewTidePrev%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstTidePrev)) then
            FirstTidePrev         => NewTidePrev
            ObjTidePrev           => NewTidePrev
        else
            PreviousTidePrev      => FirstTidePrev
            ObjTidePrev           => FirstTidePrev%Next
            do while (associated(ObjTidePrev))
                PreviousTidePrev  => ObjTidePrev
                ObjTidePrev       => ObjTidePrev%Next
            enddo
            ObjTidePrev           => NewTidePrev
            PreviousTidePrev%Next => NewTidePrev
        endif


    end subroutine AddTidePrev

end program MohidTidePreview

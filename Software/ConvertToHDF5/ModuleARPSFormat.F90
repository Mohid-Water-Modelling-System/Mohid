!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid ConvertToHDF5
! MODULE        : ARPSFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Pablo Carracedo Garcia
! DESCRIPTION   : Module to convert ARPS format files into HDF5 format.
!
!------------------------------------------------------------------------------


Module ModuleARPSFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertARPSFormat
    private ::      ReadOptions
    private ::      OpenAndReadARPSFile
    private ::          AddField
    private ::          AddDate
    private ::      OutputFields
    private ::          Open_HDF5_OutPut_File
    private ::      KillARPSFormat
    

    !Parameters----------------------------------------------------------------
    integer, parameter                          :: DotGrid      = 1
    integer, parameter                          :: CrossGrid    = 2

    !Types---------------------------------------------------------------------
    
    private :: T_Date
    type       T_Date
        type(T_Time)                            :: Date
        type(T_Date), pointer                   :: Next
    end type  T_Date


    private :: T_Field
    type       T_Field
        character(len=StringLength)             :: Name
        character(len=StringLength)             :: Units
        integer                                 :: GridLocation
        type(T_Time)                            :: Date
        real, dimension(:,:),       pointer     :: Values2D
        integer                                 :: OutputNumber         = 1
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: Next
    end type  T_Field

    
    private :: T_MM5Format
    type       T_MM5Format
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: Unit
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:  ),       pointer     :: XX
        real, dimension(:  ),       pointer     :: YY
        real                                    :: Xorig
        real                                    :: Yorig
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: FirstField         
        type(T_Date),               pointer     :: FirstDate     
    end type  T_MM5Format

    type(T_MM5Format), pointer                  :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertARPSFormat(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions

        call OpenAndReadARPSFile

        call ConstructGrid

        call OutputFields

        call KillARPSFormat

        STAT = SUCCESS_


    end subroutine ConvertARPSFormat

    !------------------------------------------------------------------------

    subroutine ReadOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------
       
        call GetData(Me%FileName,                                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleMM5Format - ERR02'
        end if

        call GetData(Me%GridFileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR03'


        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleMM5Format - ERR04'

    end subroutine ReadOptions


    !--------------------------------------------------------------------------

    subroutine OpenAndReadARPSFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist

        !Others----------------------------------------------------------------
        real, pointer,     dimension(:,:)           :: u, v, rlon, rlat
        character(len=8)                            :: fecha
        integer                                     :: Year, Month, Day, Hour, Minute, Second
        integer                                     :: imax, jmax, i, j, gInst, n
        integer                                     :: imedio, jmedio
        logical                                     :: First   = .TRUE.
        type(T_Field), pointer                      :: NewField
        type(T_Date ), pointer                      :: ActualTime
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: WILB, WIUB, WJLB, WJUB

        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading ARPS output file...'
        
        nullify(NewField)
        nullify(Me%FirstField)
        nullify(Me%FirstDate )
        nullify(Me%Bathymetry)
        nullify(Me%XX        )
        nullify(Me%YY        )

        !Verifies if file exists
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'ARPS File does not exist'
            stop 'OpenAndReadMM5File - ModuleMM5Format - ERR01'
        endif


        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR02'

        open(Unit   = Me%Unit,          &
             File   = Me%FileName,      &
             STATUS = 'OLD',            &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR03'

    !Reads user Data
    write (*,fmt=8)
  8 format('Grid cells in X-Dimension :')
    read  (*,*)jmax
    write (*,fmt=9)
  9 format('Grid cells in Y-Dimension :')
    read  (*,*)imax

    allocate(u(1:imax , 1:jmax))
    allocate(v(1:imax , 1:jmax))
    allocate(rlon(1:imax , 1:jmax))
    allocate(rlat(1:imax , 1:jmax))

    gInst = 0
    !Reads file Input file
    do 

        read(Me%Unit,900, end=10) fecha, n

        do i=1,imax
        do j=1,jmax
            read(Me%Unit,901) rlat(i,j),rlon(i,j),u(i,j),v(i,j)
        enddo
        enddo

        read(fecha(1:4), *)Year
        read(fecha(5:6), *)Month
        read(fecha(7:8), *)Day
        Hour   = n
        Minute = 0.
        Second = 0.
        
        call AddDate(Me%FirstDate, ActualTime)
        call SetDate(ActualTime%Date, Year, Month, Day, Hour, Minute, Second)
        Actualtime%Date = Actualtime%Date + 0.0

        gInst = gInst + 1

!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        call AddField(Me%FirstField, NewField)

        NewField%WorkSize%ILB   = 1
        NewField%WorkSize%IUB   = imax
        NewField%WorkSize%JLB   = 1
        NewField%WorkSize%JUB   = jmax

        NewField%Size%ILB       = NewField%WorkSize%ILB - 1
        NewField%Size%IUB       = NewField%WorkSize%IUB + 1
        NewField%Size%JLB       = NewField%WorkSize%JLB - 1
        NewField%Size%JUB       = NewField%WorkSize%IUB + 1

        WILB                    = NewField%WorkSize%ILB 
        WIUB                    = NewField%WorkSize%IUB 
        WJLB                    = NewField%WorkSize%JLB 
        WJUB                    = NewField%WorkSize%JUB 
                                            
        ILB                     = NewField%Size%ILB
        IUB                     = NewField%Size%IUB
        JLB                     = NewField%Size%JLB
        JUB                     = NewField%Size%JUB

        NewField%Name   = GetPropertyName(WindVelocityX_)
        NewField%Date   = ActualTime%Date
        NewField%Units  = 'm/s'
                    
        allocate(NewField%Values2D(ILB:IUB, JLB:JUB))
        NewField%Values2D(WILB:WIUB,WJLB:WJUB) = u(WILB:WIUB,WJLB:WJUB)

!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        call AddField(Me%FirstField, NewField)

        NewField%WorkSize%ILB   = 1
        NewField%WorkSize%IUB   = imax
        NewField%WorkSize%JLB   = 1
        NewField%WorkSize%JUB   = jmax

        NewField%Size%ILB       = NewField%WorkSize%ILB - 1
        NewField%Size%IUB       = NewField%WorkSize%IUB + 1
        NewField%Size%JLB       = NewField%WorkSize%JLB - 1
        NewField%Size%JUB       = NewField%WorkSize%IUB + 1

        WILB                    = NewField%WorkSize%ILB 
        WIUB                    = NewField%WorkSize%IUB 
        WJLB                    = NewField%WorkSize%JLB 
        WJUB                    = NewField%WorkSize%JUB 
                                            
        ILB                     = NewField%Size%ILB
        IUB                     = NewField%Size%IUB
        JLB                     = NewField%Size%JLB
        JUB                     = NewField%Size%JUB

        NewField%Name   = GetPropertyName(WindVelocityY_)
        NewField%Date   = ActualTime%Date
        NewField%Units  = 'm/s'
                    
        allocate(NewField%Values2D(ILB:IUB, JLB:JUB))
        NewField%Values2D(WILB:WIUB,WJLB:WJUB) = v(WILB:WIUB,WJLB:WJUB)
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

        !Calcula a malha:
        if (First) then

            Me%WorkSize%ILB = 1
            Me%WorkSize%IUB = imax
            Me%WorkSize%JLB = 1
            Me%WorkSize%JUB = jmax

            Me%Size%ILB     = Me%WorkSize%ILB - 1
            Me%Size%IUB     = Me%WorkSize%IUB + 1
            Me%Size%JLB     = Me%WorkSize%JLB - 1
            Me%Size%JUB     = Me%WorkSize%IUB + 1

            WILB            = Me%WorkSize%ILB 
            WIUB            = Me%WorkSize%IUB 
            WJLB            = Me%WorkSize%JLB 
            WJUB            = Me%WorkSize%JUB 

            ILB             = Me%Size%ILB
            IUB             = Me%Size%IUB
            JLB             = Me%Size%JLB
            JUB             = Me%Size%JUB

            allocate(Me%Bathymetry(ILB:IUB, JLB:JUB))
            Me%Bathymetry(ILB:IUB, JLB:JUB) = 0.

            !Calculates XX
            allocate(Me%XX(JLB:JUB))

            imedio = imax / 2
            Me%Xorig = rlon(imedio, 1) - (rlon(imedio, 2) - rlon(imedio, 1))
            Me%XX(1) = 0.0
            do j = 1, jmax - 1
                Me%XX(j+1) = (rlon(imedio, j+1)+rlon(imedio, j))/2.0 - Me%Xorig
            enddo
            Me%XX(jmax+1) = rlon(imedio, jmax) + (rlon(imedio, jmax)-rlon(imedio, jmax-1))/2.0 - Me%Xorig

            !Calculates YY
            allocate(Me%YY(ILB:IUB))

            jmedio = jmax / 2
            Me%Yorig = rlat(1, jmedio) - (rlat(2, jmedio) - rlat(1, jmedio))
            Me%YY(1) = 0.0
            do i = 1, imax-1
                Me%YY(i+1) = (rlat(i+1, jmedio)+rlat(i, jmedio))/2.0 - Me%Yorig
            enddo
            Me%YY(imax+1) = rlat(imax, jmedio) + (rlat(imax, jmedio)-rlat(imax-1, jmedio))/2.0 - Me%Yorig
                                   
            First = .false.

        endif


    enddo

10 continue


900         format(a8,2x,i2)
901         format(4(f9.5,1x))

800     format(i2)
801     format(1600(4(f9.5,1x)))


100     write(*,*)'Finished reading ARPS output file.'

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleMM5Format - ERR04'



    end subroutine OpenAndReadARPSFile
    
    
    !------------------------------------------------------------------------

    
    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        real                        :: Latitude, Longitude
        integer                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'


        call WriteGridData (FileName        = Me%GridFileName,              &
                            XX              = Me%XX,                        &
                            YY              = Me%YY,                        &
                            COMENT1         = " MM5 Grid based on file :",  &
                            COMENT2         = Me%FileName,                  &
                            WorkSize        = Me%WorkSize,                  & 
                            CoordType       = 4,                            &
                            Xorig           = Me%Xorig,                     &
                            Yorig           = Me%Yorig,                     &
                            Zone            = 0,                            &
                            Grid_Angle      = 0.,                           &
                            Latitude        = Latitude,                     &
                            Longitude       = Longitude,                    &
                            GridData2D_Real = Me%Bathymetry,                &
                            Overwrite       = ON,                           &
                            FillValue       = -99.,                         &
                            STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMM5Format - ERR01'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleMM5Format - ERR02'


    end subroutine ConstructGrid


    
    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                   :: FirstField
        type (T_Field), pointer                   :: ObjField

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        type (T_Field), pointer                   :: PreviousField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
            FirstField         => NewField
            ObjField           => NewField
        else
            PreviousField      => FirstField
            ObjField           => FirstField%Next
            do while (associated(ObjField))
                PreviousField  => ObjField
                ObjField       => ObjField%Next
            enddo
            ObjField           => NewField
            PreviousField%Next => NewField
        endif


    end subroutine AddField
    
    
    !------------------------------------------------------------------------


    subroutine AddDate (FirstDate, ObjDate)

        !Arguments-------------------------------------------------------------
        type (T_Date), pointer                   :: FirstDate
        type (T_Date), pointer                   :: ObjDate

        !Local-----------------------------------------------------------------
        type (T_Date), pointer                   :: NewDate
        type (T_Date), pointer                   :: PreviousDate
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewDate)
        nullify  (NewDate%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstDate)) then
            FirstDate         => NewDate
            ObjDate           => NewDate
        else
            PreviousDate      => FirstDate
            ObjDate           => FirstDate%Next
            do while (associated(ObjDate))
                PreviousDate  => ObjDate
                ObjDate       => ObjDate%Next
            enddo
            ObjDate           => NewDate
            PreviousDate%Next => NewDate
        endif


    end subroutine AddDate
    
    
    !------------------------------------------------------------------------
    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, OutputNumber
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate

        !Begin-----------------------------------------------------------------
        
        write(*,*)
        write(*,*)'Writing HDF5 file...'

        call Open_HDF5_OutPut_File

        OutputNumber = 1
        CurrentDate  => Me%FirstDate

        do while(associated(CurrentDate))

            call ExtractDate   (CurrentDate%Date,                           &
                                AuxTime(1), AuxTime(2), AuxTime(3),         &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR01'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR02'


            Field => Me%FirstField

            do while(associated(Field))
            
                
                if(Field%Date == CurrentDate%Date)then
                
                    Field%OutputNumber = OutputNumber
            
                    call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                       Field%WorkSize%JLB, Field%WorkSize%JUB, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR05'

                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//Field%Name,                         &
                                       Field%Name,                                      &
                                       Field%Units,                                     &
                                       Array2D      = Field%Values2D,                   &
                                       OutputNumber = Field%OutputNumber,               &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR06'

                end if

                Field => Field%Next

            end do

            OutputNumber = OutputNumber + 1

            CurrentDate => CurrentDate%Next

        end do

        write(*,*)
        write(*,*)'Closing HDF5 file...'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR07'



    end subroutine OutputFields


    !----------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        integer,    dimension(:,:,:), pointer       :: WaterPoints3D

        !----------------------------------------------------------------------

        
        
        allocate(WaterPoints3D(Me%Size%ILB:Me%Size%IUB,&
                               Me%Size%JLB:Me%Size%JUB, 1))

        WaterPoints3D = 1

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
        
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR02'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR04'            
   

        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, 1,1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints", "-",    &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR08'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR09'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine Open_HDF5_OutPut_File


    !--------------------------------------------------------------------------
    
    subroutine KillARPSFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMM5Format - ModuleMM5Format - ERR01'
        
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMM5Format - ModuleMM5Format - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillMM5Format - ModuleMM5Format - ERR03'


        deallocate(Me%Bathymetry)
        deallocate(Me%XX)
        deallocate(Me%YY)
        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillARPSFormat

    !--------------------------------------------------------------------------

end module ModuleARPSFormat

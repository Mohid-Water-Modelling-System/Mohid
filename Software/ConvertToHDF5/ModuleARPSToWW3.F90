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


Module ModuleARPSToWW3

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
    public  :: ConvertARPSWW3Format
    private ::      ReadOptions
    private ::      OpenAndReadARPSFile
    private ::          AddField
    private ::          AddDate
    private ::      OutputFields
    private ::          Open_HDF5_OutPut_File
    private ::      KillARPSFormat
    

    !Parameters----------------------------------------------------------------

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
        character(len=PathLength)               :: ARPSGridFileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:,:),       pointer     :: ConnectionX
        real, dimension(:,:),       pointer     :: ConnectionY
        real                                    :: Xorig
        real                                    :: Yorig
        logical                                 :: InputGridInMohidFormat = .false.
        type(T_Size2D)                          :: Size, WorkSize
        type(T_Time  )                          :: StartDate, EndDate
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

    subroutine ConvertARPSWW3Format(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions

        call ConstructGrid

        call OpenAndReadARPSFile

        call OutputFields

        call KillARPSFormat

        STAT = SUCCESS_


    end subroutine ConvertARPSWW3Format

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
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleARPSToWW3 - ERR02'
        end if


        call GetData(Me%InputGridInMohidFormat,                         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'INPUTGRIDINMOHIDFORMAT',           &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR01'


        
        call GetData(Me%ARPSGridFileName,                               &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'ARPS_GRID_FILENAME',               &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleARPSToWW3 - ERR02'
        end if


        call GetData(Me%GridFileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR03'


        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR04'

        call GetData(Me%StartDate,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR05'

        call GetData(Me%EndDate,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleARPSToWW3 - ERR06'



    end subroutine ReadOptions


    !--------------------------------------------------------------------------

    subroutine OpenAndReadARPSFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist

        !Others----------------------------------------------------------------
        real, pointer,     dimension(:,:)           :: U, V
        character(len=8)                            :: StrDate
        character(len=6)                            :: StrTime
        character(len=2048)                         :: string
        integer                                     :: Year, Month, Day, Hour, Minute, Second
        integer                                     :: i, j, Instant, istr
        type(T_Field), pointer                      :: NewField
        type(T_Date ), pointer                      :: ActualTime
        logical                                     :: NewTime

        !Begin-----------------------------------------------------------------
        
        write(*,*)'---------------------------'
        write(*,*)
        write(*,*)'Reading ARPS output file...', trim(Me%FileName)
        
        nullify(NewField)
        nullify(Me%FirstField)
        nullify(Me%FirstDate )

        !Verifies if file exists
        inquire(file = trim(Me%FileName), exist = exist)
        if (.not. exist) then
            write(*,*)'ARPS File does not exist'
            stop 'OpenAndReadMM5File - ModuleARPSToWW3 - ERR01'
        endif


        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleARPSToWW3 - ERR02'

        open(Unit   = Me%Unit,          &
             File   = Me%FileName,      &
             STATUS = 'OLD',            &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleARPSToWW3 - ERR03'

        allocate(U(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        allocate(V(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        500 format(a8, 2x, a6)

        rewind(Me%Unit)

        NewTime = .false.

        Instant = 0
        !Reads file Input file
        do

            read(Me%Unit, fmt=500, end=900)StrDate, StrTime

            read(StrDate(1:4),*)Year
            read(StrDate(5:6),*)Month
            read(StrDate(7:8),*)Day

            read(StrTime(1:2),*)Hour
            read(StrTime(3:4),*)Minute
            read(StrTime(5:6),*)Second

            call AddDate(Me%FirstDate, ActualTime)
            call SetDate(ActualTime%Date, Year, Month, Day, Hour, Minute, Second)
            ActualTime%Date = ActualTime%Date + 0.

            call AddField(Me%FirstField, NewField)
            allocate(NewField%Values2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            NewField%Name   = GetPropertyName(WindVelocityX_)
            NewField%Date   = ActualTime%Date
            NewField%Units  = 'm/s'

            do i=Me%WorkSize%ILB, Me%WorkSize%IUB+1
            
                read(Me%Unit, '(A)')string

                istr = 1

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB+1

                    read(string(istr:istr+7), *)NewField%Values2D(i,j)

                    istr = istr+7
                enddo

            enddo
                
            call AddField(Me%FirstField, NewField)
            allocate(NewField%Values2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            NewField%Name   = GetPropertyName(WindVelocityY_)
            NewField%Date   = ActualTime%Date
            NewField%Units  = 'm/s'

            do i=Me%WorkSize%ILB, Me%WorkSize%IUB+1

                read(Me%Unit, '(A)')string

                istr = 1

                do j = Me%WorkSize%JLB, Me%WorkSize%JUB+1

                    read(string(istr:istr+7), *)NewField%Values2D(i,j)

                    istr = istr+7
                enddo

            enddo

        enddo

        900 write(*,*)'Finished reading ARPS output file.'

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadMM5File - ModuleARPSToWW3 - ERR04'



    end subroutine OpenAndReadARPSFile
    
    
    !------------------------------------------------------------------------

    
    subroutine ConstructGrid
        
        !Local-----------------------------------------------------------------
        real                        :: Latitude, Longitude
        integer                     :: STAT_CALL
        integer                     :: BufferSize, Unit, stringsize
        character(len=2048)         :: string
        integer                     :: i, j, istr, ix, iy

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Constructing grid...'


        nullify(Me%ConnectionY)
        nullify(Me%ConnectionY)

        if(Me%InputGridInMohidFormat) then

            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%ARPSGridFileName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR04'

            call GetHorizontalGridSize(Me%ObjHorizontalGrid, Size = Me%Size, WorkSize = Me%WorkSize, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR04'
           
            allocate(Me%Bathymetry (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            Me%Bathymetry = 0.


        else

            call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR01'

            open(Unit   = Unit,                     &
                 File   = Me%ARPSGridFileName,      &
                 STATUS = 'OLD',                    &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR02'

            BufferSize = 0

            do

                BufferSize = BufferSize + 1   

                read(Unit, '(A)', end = 100)string
                
            end do


            100 write(*,*) "Found grid limits..."
            write(*,*)

            stringsize = len_trim(string) / 7
    
            Me%Size%ILB = 0;                    Me%Size%IUB = BufferSize/2
            Me%Size%JLB = 0;                    Me%Size%JUB = stringsize
        
            Me%WorkSize%ILB = Me%Size%ILB + 1;  Me%WorkSize%IUB = Me%Size%IUB - 1
            Me%WorkSize%JLB = Me%Size%JLB + 1;  Me%WorkSize%JUB = Me%Size%JUB - 1
        
            allocate(Me%ConnectionX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Me%ConnectionY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            allocate(Me%Bathymetry (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            Me%Bathymetry = 0.; Me%ConnectionX = null_real; Me%ConnectionY = null_real

            rewind(Unit)

            i = 0; ix = 0; iy = 0

            do
                i = i + 1   
                
                read(Unit, '(A)', end = 200)string

                if(i .le. BufferSize/2)then
                
                    istr = 1

                    iy = iy + 1

                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB + 1

                        read(string(istr:istr+7), *)Me%ConnectionY(iy,j)

                        istr = istr+7

                    end do

                else

                    istr = 1

                    ix = ix + 1


                    do j = Me%WorkSize%JLB, Me%WorkSize%JUB + 1

                        read(string(istr:istr+7), *)Me%ConnectionX(ix,j)

                        istr = istr+7

                    end do
            
                end if
             
            end do

            200 write(*,*)"Writing grid in Mohid format"

            call WriteGridData (FileName        = Me%GridFileName,              &
                                ConnectionX     = Me%ConnectionX,               &
                                ConnectionY     = Me%ConnectionY,               &
                                COMENT1         = " ARPS Grid based on file :", &
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
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR03'

            call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR04'

            call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleARPSToWW3 - ERR05'


        end if




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

            if(CurrentDate%Date .ge. Me%StartDate .and. CurrentDate%Date .le. Me%EndDate)then

                call ExtractDate   (CurrentDate%Date,                           &
                                    AuxTime(1), AuxTime(2), AuxTime(3),         &
                                    AuxTime(4), AuxTime(5), AuxTime(6))

                TimePtr => AuxTime

                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleARPSToWW3 - ERR01'


                call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",             &
                                     Array1D = TimePtr,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleARPSToWW3 - ERR02'


                Field => Me%FirstField

                do while(associated(Field))
            
                
                    if(Field%Date == CurrentDate%Date)then
                
                        Field%OutputNumber = OutputNumber
            
                        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                                           Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleARPSToWW3 - ERR05'

                        call HDF5WriteData(Me%ObjHDF5,                                      &
                                           "/Results/"//Field%Name,                         &
                                           Field%Name,                                      &
                                           Field%Units,                                     &
                                           Array2D      = Field%Values2D,                   &
                                           OutputNumber = OutputNumber,                     &
                                           STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleARPSToWW3 - ERR06'

                    end if

                    Field => Field%Next

                end do

                OutputNumber = OutputNumber + 1

            end if

            CurrentDate => CurrentDate%Next

        end do

        write(*,*)
        write(*,*)'Closing HDF5 file...'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleARPSToWW3 - ERR07'



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
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR02'

        
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR03'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR04'            
   

        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, 1,1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR07'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints3D", "-",    &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR08'
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleARPSToWW3 - ERR09'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine Open_HDF5_OutPut_File


    !--------------------------------------------------------------------------
    
    subroutine KillARPSFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMM5Format - ModuleARPSToWW3 - ERR01'
        
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillMM5Format - ModuleARPSToWW3 - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillMM5Format - ModuleARPSToWW3 - ERR03'


        deallocate(Me%Bathymetry)
        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillARPSFormat

    !--------------------------------------------------------------------------

end module ModuleARPSToWW3

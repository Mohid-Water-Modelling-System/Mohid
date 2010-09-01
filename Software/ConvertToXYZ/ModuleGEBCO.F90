!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : GEBCO
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : June 2004
! REVISION      : Luis Fernandes
! DESCRIPTION   : Module to convert GEBCO file into XYZ format within a defined window.
!                 Source on the Internet: http://www.ngdc.noaa.gov/mgg/gebco/gebco.html
!
!------------------------------------------------------------------------------


!   INPUT_FILENAME              : char              -           !Path to input file to convert
!   OUTPUT_FILENAME             : char              -           !Path to XYZ file generated
!   WRITE_AS_BATHYMETRY         : 0/1               -           !Write XYZ values multiplied by -1.
!   CONVERT_TO_HDF5             : 0/1               -           !Convert to HDF5

Module ModuleGEBCO

    use ModuleGlobalData
    use ModuleDrawing
    use ModuleEnterData
    use ModuleHDF5
    use netcdf90

    implicit none

    private 
    
    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertGEBCO
    private ::   ReadFileAndWriteOutput
    private ::   KillGEBCO

    !Parameters----------------------------------------------------------------
    integer, parameter                  :: iLongitude   = 1201
    integer, parameter                  :: iLatitude    = 1201   

    !Types---------------------------------------------------------------------
    type      T_GEBCO
        integer                         :: Unit
        integer                         :: OutputUnit
        type(T_Limits)                  :: Window
        type(T_Limits)                  :: OriginalWindow
        type(T_Size2D)                  :: Size
        real(8)                         :: SpacingX, SpacingY
        logical                         :: WriteAsBathymetry, ConvertToHDF5
        character(len=StringLength)     :: FileName, OutputFileName
        integer, dimension(:), pointer  :: Topography
        real,    dimension(:), pointer  :: Latitude
        real,    dimension(:), pointer  :: Longitude
        integer                         :: nPoints
        integer                         :: ObjEnterData = 0
    end type  T_GEBCO

    type(T_GEBCO), pointer              :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertGEBCO(Window, EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        type(T_Limits),    intent(IN )                  :: Window
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_

        nullify (Me)
        allocate(Me)

        Me%Window               = Window
        Me%ObjEnterData         = AssociateInstance(mEnterData_, EnterDataID)

        call ReadOptions

        call ReadFileAndWriteOutput

        call KillGEBCO

        STAT = SUCCESS_

    end subroutine ConvertGEBCO
    
    !--------------------------------------------------------------------------

    subroutine ReadOptions
        
        !Local-----------------------------------------------------------------
        integer                 :: iflag, STAT_CALL

        !Begin---------------------------------------------------------

        call GetData(Me%FileName,                                   &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'INPUT_FILENAME',               &
                     ClientModule = 'ModuleGEBCO',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGEBCO - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
             stop 'ReadOptions - ModuleGEBCO - ERR02'
        end if



        call GetData(Me%OutputFileName,                             &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'OUTPUT_FILENAME',              &
                     ClientModule = 'ModuleGEBCO',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGEBCO - ERR03'

        call GetData(Me%WriteAsBathymetry,                          &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'WRITE_AS_BATHYMETRY',          &
                     Default      = .false.,                        &
                     ClientModule = 'ModuleGEBCO',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGEBCO - ERR04'

        call GetData(Me%ConvertToHDF5,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'CONVERT_TO_HDF5',              &
                     Default      = .false.,                        &
                     ClientModule = 'ModuleGEBCO',                  &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGEBCO - ERR05'


    end subroutine ReadOptions

    !--------------------------------------------------------------------------

    
    subroutine ReadFileAndWriteOutput

        !Local-----------------------------------------------------------------
        integer                                     :: stat, iPoint, i, j, STAT_CALL
        logical                                     :: exist
        integer                                     :: FileID, nDimensions, nVariables, nAttributes
        integer                                     :: VariableID, VariableType
        integer, dimension(4)                       :: DimensionsID
        integer, dimension(:), pointer              :: DimensionsArray
        real(8), dimension(:), pointer              :: SpacingArray, XRange, YRange
        integer                                     :: DimensionSize
        character(len=80)                           :: VariableName, DimensionName

        !Begin---------------------------------------------------------

        nullify(XRange, YRange, SpacingArray, DimensionsArray)
        nullify(Me%Longitude, Me%Latitude, Me%Topography)

        call UnitsManager(Me%OutputUnit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR02'

        open(Unit   = Me%OutputUnit,                &
             File   = Me%OutputFileName,            &
             STATUS = 'UNKNOWN',                    &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR03'


        write(*,*)
        write(*,*)"Opening GEBCO file..."

        !Verifies if file exists
        inquire(file = Me%FileName, exist = exist)
        if (.not. exist) then
            write(*,*)'GEBCO file does not exist'
            stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR01'
        endif

        stat = nf90_open(Me%FileName, NF90_NOWRITE, FileID)
        if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR02'

        stat = nf90_inquire(FileID, nDimensions, nVariables, nAttributes)
        if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR03'


        do VariableID = 1, nVariables
        
            stat = nf90_inquire_variable(FileID, VariableID, VariableName, &
                                         VariableType, nDimensions, DimensionsID)
            if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR04'

            stat = nf90_inquire_dimension(FileID, DimensionsID(1), DimensionName, DimensionSize)

            select case(trim(VariableName))


                case('x_range')

                    allocate(XRange(DimensionSize))

                    stat = nf90_get_var(FileID, VariableID, XRange)
                    if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR05'

                    Me%OriginalWindow%Left      = XRange(1)
                    Me%OriginalWindow%Right     = XRange(2)


                case('y_range')

                    allocate(YRange(DimensionSize))

                    stat = nf90_get_var(FileID, VariableID, YRange)
                    if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR05'

                    Me%OriginalWindow%Bottom    = YRange(1)
                    Me%OriginalWindow%Top       = YRange(2)

                case('spacing')

                    allocate(SpacingArray(DimensionSize))

                    stat = nf90_get_var(FileID, VariableID, SpacingArray)
                    if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR06'

                    Me%SpacingX = SpacingArray(1)
                    Me%SpacingY = SpacingArray(2)
                    

                case('dimension')

                    allocate(DimensionsArray(DimensionSize))

                    stat = nf90_get_var(FileID, VariableID, DimensionsArray)
                    if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR07'

                    Me%Size%ILB = 1
                    Me%Size%JLB = 1
                    Me%Size%IUB = DimensionsArray(1)
                    Me%Size%JUB = DimensionsArray(2)

                    allocate(Me%Longitude(1:Me%Size%JUB))
                    allocate(Me%Latitude (1:Me%Size%IUB))


                case('z')

                    Me%nPoints = DimensionSize

                    allocate(Me%Topography(Me%nPoints))

                    stat = nf90_get_var(FileID, VariableID, Me%Topography)
                    if (stat /= nf90_noerr) stop 'ReadFileAndWriteOutput - ModuleGEBCO - ERR08'

            end select

        end do

        write(*,*)
        write(*,*)"Opened GEBCO file..."

        write(*,*)
        write(*,*)"Converting GEBCO file. Please wait..."

        Me%Longitude(1) = Me%OriginalWindow%Left
        Me%Latitude(1)  = Me%OriginalWindow%Top

        do iPoint = Me%Size%JLB + 1, Me%Size%JUB

            Me%Longitude(iPoint) = Me%Longitude(iPoint - 1) + Me%SpacingX

        end do

        
        do iPoint = Me%Size%ILB + 1, Me%Size%IUB

            Me%Latitude (iPoint) = Me%Latitude (iPoint - 1) - Me%SpacingY

        end do


        write(Me%OutputUnit,*)"<begin_xyz>"

        iPoint = 1

        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB

            if(Me%Topography(iPoint) .ge. Me%Window%MinimumValue .and. &
               Me%Topography(iPoint) .le. Me%Window%MaximumValue)then

                if(Me%Longitude(i) .ge. Me%Window%Left  .and. &
                   Me%Longitude(i) .le. Me%Window%Right .and. &
                   Me%Latitude (j) .le. Me%Window%Top   .and. &
                   Me%Latitude (j) .ge. Me%Window%Bottom)then
            
                    if(Me%WriteAsBathymetry)then
                        write(Me%OutputUnit,99) Me%Longitude(i), Me%Latitude(j),-Me%Topography(iPoint)
                    else
                        write(Me%OutputUnit,99) Me%Longitude(i), Me%Latitude(j), Me%Topography(iPoint)
                    endif

                 endif

            end if

            iPoint = iPoint + 1

        enddo
        enddo


        if(Me%ConvertToHDF5)then

            call ConvertToHDF5

        end if


        write(Me%OutputUnit,*)'<end_xyz>'

        99 format(1x, f10.5, 1x, f10.5,1x, i6)

        deallocate(XRange, YRange, SpacingArray, DimensionsArray)
        deallocate(Me%Longitude, Me%Latitude, Me%Topography)


        write(*,*)'Finished converting'


    end subroutine ReadFileAndWriteOutput
    
    !------------------------------------------------------------------

    subroutine ConvertToHDF5
        
        !Local-----------------------------------------------------------------
        integer,    dimension(:,:), pointer         :: WaterPoints2D
        real,       dimension(:,:), pointer         :: Bathymetry
        real,       dimension(:,:), pointer         :: ConnectionX, ConnectionY
        integer                                     :: iPoint, i, j, HDF5_CREATE
        integer                                     :: ObjHDF5, STAT_CALL

        !Begin-----------------------------------------------------------------


        write(*,*)'Converting to HDF5'

        nullify(WaterPoints2D, Bathymetry, ConnectionX, ConnectionY)

        allocate(WaterPoints2D(Me%Size%ILB:Me%Size%IUB,   Me%Size%JLB:Me%Size%JUB  ))
        allocate(Bathymetry   (Me%Size%ILB:Me%Size%IUB,   Me%Size%JLB:Me%Size%JUB  ))
        allocate(ConnectionX  (Me%Size%ILB:Me%Size%IUB+1, Me%Size%JLB:Me%Size%JUB+1))
        allocate(ConnectionY  (Me%Size%ILB:Me%Size%IUB+1, Me%Size%JLB:Me%Size%JUB+1))

        iPoint = 1

        Bathymetry = 0.; ConnectionX = 0.; ConnectionY = 0.
        

        !Compute ConnectionX
        ConnectionX(:, Me%Size%JLB) = Me%Longitude(Me%Size%JLB) - Me%SpacingX/2.
        
        do j = Me%Size%JLB, Me%Size%JUB
        do i = Me%Size%ILB, Me%Size%IUB+1

            ConnectionX(i, j+1) = Me%Longitude(j) + Me%SpacingX/2.
        
        end do
        end do

        ConnectionX(:, Me%Size%JUB+1) = Me%Longitude(Me%Size%JUB) + Me%SpacingX/2.



        !Compute ConnectionX
        ConnectionY(Me%Size%ILB,: )   = Me%Latitude (Me%Size%IUB) - Me%SpacingY/2.

        do i = Me%Size%ILB, Me%Size%IUB
        do j = Me%Size%JLB, Me%Size%JUB+1

            ConnectionY(i+1, j)   = Me%Latitude(Me%Size%IUB - i + 1) + Me%SpacingY/2.

        end do
        end do
        
        ConnectionY(Me%Size%IUB+1, :) = Me%Latitude (Me%Size%JLB)  + Me%SpacingY/2.


        !Fill Bathymetry
        
        do i = Me%Size%IUB, Me%Size%ILB, -1
        do j = Me%Size%JLB, Me%Size%JUB

            if(Me%Topography(iPoint) .ge. Me%Window%MinimumValue .and. &
               Me%Topography(iPoint) .le. Me%Window%MaximumValue)then

                WaterPoints2D(i,j) = 1
                    
                if(Me%WriteAsBathymetry)then
                    Bathymetry(i,j) = -Me%Topography(iPoint)
                else
                    Bathymetry(i,j) =  Me%Topography(iPoint)
                endif

            else

                WaterPoints2D(i,j) = 0

            end if

            iPoint = iPoint + 1

        enddo
        enddo



        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(ObjHDF5, trim(Me%OutputFileName)//".hdf5", HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR01'
        
        call HDF5SetLimits  (ObjHDF5, Me%Size%ILB, Me%Size%IUB, &
                             Me%Size%JLB,Me%Size%JUB,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR02'

        call HDF5WriteData   (ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR03'            
        
        call HDF5WriteData   (ObjHDF5, "/Grid", "WaterPoints", "-",       &
                              Array2D =  WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR04'            

        call HDF5SetLimits  (ObjHDF5, Me%Size%ILB, Me%Size%IUB+1, &
                             Me%Size%JLB,Me%Size%JUB+1,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR05'


        call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionX", "-",       &
                              Array2D =  ConnectionX, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR06'            

        call HDF5WriteData   (ObjHDF5, "/Grid", "ConnectionY", "-",       &
                              Array2D =  ConnectionY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR07'            

        call HDF5FlushMemory (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR08'      
       
        call KillHDF5 (ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConvertToHDF5 - ModuleGEBCO - ERR09'      

        write(*,*)

        deallocate(WaterPoints2D, Bathymetry, ConnectionX, ConnectionY)
       
        write(*,*)'Finished converting to HDF5'

    end subroutine ConvertToHDF5
        
    !------------------------------------------------------------------

    subroutine KillGEBCO
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnitsManager(Me%OutputUnit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'KillGEBCO - ModuleGEBCO - ERR02'


        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0)           stop 'KillGEBCO - ModuleGEBCO - ERR03'


        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillGEBCO

    !--------------------------------------------------------------------------

 
end module ModuleGEBCO










!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : TestBoxDif
! PROGRAM       : ConvertProfiles
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2005
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Converts TZ profile HDF outputs to a Matlab speed load format
!
!------------------------------------------------------------------------------

program ConvertProfiles

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData
    use ModuleFunctions
    use ModuleHDF5

    implicit none

    integer                             :: ObjHDF5              = 0
    integer                             :: ObjOutHDF5           = 0
    integer                             :: ObjEnterData         = 0
    character(PathLength)               :: DataFile = 'Profile.dat'
    character(PathLength)               :: HDF5File, OutHDF5File
    real, dimension(:  ),  pointer      :: Z, Profile
    integer                             :: STAT_CALL, iflag, HDF5_READ, HDF5_CREATE


    call StartConvertProfiles

    contains
    
    !--------------------------------------------------------------------------

    subroutine StartConvertProfiles

        !Local-----------------------------------------------------------------
        character(len=StringLength)                 :: PropertyName, nPropertyName
        character(len=StringLength)                 :: GroupName, Units 
        integer                                     :: i, n, nInstants, nProperties
        integer, dimension(7)                       :: Dimensions
        type(T_Time)                                :: CurrentTime
        integer                                     :: nLayers
        real(8), dimension(:,:),  pointer           :: Field2D
        real(8), dimension(:,:),  pointer           :: Z2D
        real(8), dimension(:,:),  pointer           :: Time2D
        integer                                     :: iJulDay
        real(8)                                     :: JDay
        real                                        :: Hours, Minutes, Seconds

        !Begin-----------------------------------------------------------------

        call StartUpMohid("ConvertProfiles")

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertProfiles - ERR10'


        call GetData(HDF5File,                                                          &
                     ObjEnterData,iflag,                                                &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF5_PROFILE_FILE',                                &
                     ClientModule = 'ConvertProfiles',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertProfiles - ERR20'

        call GetData(OutHDF5File,                                                       &
                     ObjEnterData,iflag,                                                &
                     SearchType   = FromFile,                                           &
                     keyword      = 'HDF5_TZ_OUTFILE',                                  &
                     ClientModule = 'ConvertProfiles',                                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertProfiles - ERR30'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConvertProfiles - ERR40'


        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ, HDF5_CREATE = HDF5_CREATE)

        call ConstructHDF5 (ObjHDF5, trim(HDF5File), HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR50'


        call ConstructHDF5 (ObjOutHDF5, trim(OutHDF5File), HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR60'


        call GetHDF5GroupNumberOfItems(ObjHDF5, "/Time", nInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR70'

        call GetHDF5GroupNumberOfItems(ObjHDF5, "/Results", nProperties, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR80'
 
        if(nProperties < 1)then
            stop 'No properties - ConvertProfiles - ERR90'
        end if
        
        do n = 1, nProperties

            !get property name
            call GetHDF5GroupID(ObjHDF5, "/Results", n, PropertyName, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR100'

            write(*,*)'Converting property: ', trim(PropertyName)

            GroupName     = "/Results/"//trim(adjustl(PropertyName))

            nPropertyName = PropertyName

            call GetHDF5GroupID(ObjHDF5, GroupName, 1, nPropertyName,                   &
                                Units = Units, Dimensions = Dimensions,                 &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR110'

            nLayers = Dimensions(1)

            allocate(Profile(1:nLayers            ))
            allocate(Field2D(1:nLayers,1:nInstants))
            
            call HDF5SetLimits  (ObjHDF5, 1, nLayers, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR120'


            do i = 1, nInstants

                call HDF5ReadData   (HDF5ID         = ObjHDF5,                          &
                                     GroupName      = trim(GroupName),                  &
                                     Name           = trim(PropertyName),               &
                                     Array1D        = Profile,                          &
                                     OutputNumber   = i,                                &
                                     STAT           = STAT_CALL)
                if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR130'
                
                Field2D(:,i) = Profile(:)

            enddo

            call HDF5SetLimits (ObjOutHDF5, 1, nLayers, 1, nInstants, STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR140'

            call HDF5WriteData (ObjOutHDF5, trim(GroupName), trim(PropertyName),        &
                                trim(Units), Array2D = Field2D,                         &
                                STAT = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR150'

            deallocate(Profile, Field2D)

            
        end do

        allocate(Time2D(1:nLayers,1:nInstants))
        
        write(*,*)'...'
        write(*,*)'Converting time'

        do i = 1, nInstants

            CurrentTime = HDF5TimeInstant(i)

            call JulianDay(CurrentTime, iJulDay)

            JDay = dble(iJulDay) - 1.

            call ExtractDate(CurrentTime, Hour = Hours, Minute = Minutes, Second = Seconds)

            JDay = JDay + Hours/24. + Minutes/(1440.) + Seconds/86400.

            Time2D(:, i) = JDay

        enddo

        
        
        call HDF5SetLimits (ObjOutHDF5, 1, nLayers, 1, nInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR160'

        call HDF5WriteData (ObjOutHDF5, '/Time', 'Time',                            &
                            'seconds', Array2D = Time2D,                            &
                            STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR170'

        deallocate(Time2D)


        write(*,*)'...'
        write(*,*)'Converting Z'

        allocate(Z  (1:nLayers+1            ))
        allocate(Z2D(1:nLayers  ,1:nInstants))

        call HDF5SetLimits (ObjHDF5, 1, nLayers+1, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR180'

        do i = 1, nInstants

            call HDF5ReadData   (HDF5ID         = ObjHDF5,                              &
                                 GroupName      = '/Grid/VerticalZ',                    &
                                 Name           = "Vertical",                           &
                                 Array1D        = Z,                                    &
                                 OutputNumber   = i,                                    &
                                 STAT           = STAT_CALL)
            if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR190'

            do n = 1, nLayers

                Z2D(n,i) =  0.5*(Z(n)+Z(n+1))

            enddo            
             
        enddo

        call HDF5SetLimits (ObjOutHDF5, 1, nLayers, 1, nInstants, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR200'

        call HDF5WriteData (ObjOutHDF5, '/Grid/VerticalZ', 'VerticalZ',                 &
                            'm', Array2D = Z2D,STAT = STAT_CALL)
                            
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR210'

        deallocate(Z2D)

        call KillHDF5(ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConvertProfiles - ERR999'


    end subroutine StartConvertProfiles



    type(T_Time) function HDF5TimeInstant(Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        

        !Local-----------------------------------------------------------------
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ConvertProfiles - ERR01'

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjHDF5,                                  &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ConvertProfiles - ERR02'


        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

    end function HDF5TimeInstant


end program ConvertProfiles

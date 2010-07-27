!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : ConvertToXYZ
! MODULE        : ConvertToXYZ
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2004
! REVISION      : Luis Fernandes
! DESCRIPTION   : Module to convert different base file format to MohidGIS format
! SPECIFICATIONS: Uses hdf5.lib, netcdf.lib, netcf_.lib to compile and link
!                 Mandatory that netcdf.dll is placed alongside with the executable for this to run 
!------------------------------------------------------------------------------

!DataFile
!<begin_window>
!   FILE_TYPE                   : int               -           !Type of file to convert
!   LEFT                        : real              -           !Window left coordinate
!   RIGHT                       : real              -           !Window right coordinate
!   TOP                         : real              -           !Window top coordinate
!   BOTTOM                      : real              -           !Window bottom coordinate
!   MINIMUM                     : real              -           !Minimum value to be converted
!   MAXIMUM                     : real              -           !Maximum value to be converted
!   ---specific keywords from each module---
!<begin_window>

program ConvertToXYZ

    use ModuleGlobalData
    use ModuleEnterData
    use ModuleDrawing
    use ModuleEtopo2
    use ModuleEtopo5
    use ModuleNASA
    use ModuleGEBCO
    use ModuleNOAA_ShoreLine
    use ModuleASCII
    use ModuleSRTM30

    implicit none

    !Types--------------------------------------------------------------
    type T_File
        integer                                 :: Unit
        character(len=StringLength)             :: FileName
    end type T_File

    type     T_ConvertToXYZ
        integer                                 :: FileType
        type(T_Limits)                          :: Window
        type(T_File  )                          :: Input
        type(T_File  )                          :: XYZOutput
        integer                                 :: ObjEnterData = 0
        logical                                 :: WriteAsBathymetry    = .false.
    end type T_ConvertToXYZ

    type(T_ConvertToXYZ)                        :: Me
    
    !Parameters---------------------------------------------------------
    integer, parameter                          :: NASA                 = 1
    integer, parameter                          :: Etopo2               = 2
    integer, parameter                          :: Etopo5               = 3
    integer, parameter                          :: NOAA_ShoreLine       = 4
    integer, parameter                          :: GEBCO                = 5
    integer, parameter                          :: ASCII                = 6
    integer, parameter                          :: SRTM30               = 7

    character(len=StringLength), parameter      :: DataFile             = 'ConvertToXYZ.dat'

    !Begin--------------------------------------------------------------
    
    call StartUpMohid("ConvertToXYZ")

    call StartConverting

    contains
    
    !------------------------------------------------------------------
    
    subroutine StartConverting

        !Local---------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer                                     :: ClientNumber
        logical                                     :: BlockFound

        !Begin---------------------------------------------------------

        call ConstructEnterData(Me%ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR010'


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData, ClientNumber,              &
                                        '<begin_window>', '<end_window>',           &
                                        BlockFound, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR020'

if1 :       if(STAT_CALL .EQ. SUCCESS_) then 
   
if2 :           if (BlockFound) then

                    call GetData(Me%FileType,                                       &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'FILE_TYPE',                        &
                                 default      = -99999,                             &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR030'

                    if ( Me%FileType == -99999 ) then

                        call GetData(Me%FileType,                                       &
                                     Me%ObjEnterData, iflag,                            &
                                     SearchType   = FromBlock,                          &
                                     keyword      = 'FYLE_TYPE',                        &
                                     ClientModule = 'ConvertToXYZ',                     &
                                     STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR040'
                        write(*,*)
                        write(*,*) 'Warning! The correct keyword syntax is FILE_TYPE, not FYLE_TYPE!'
                        write(*,*) 'Please correct ConvertToXYZ.dat file. '

                    endif

                    call GetData(Me%Window%Left,                                    &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'LEFT',                             &
                                 Default      = -99999.,                            &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR050'

                    call GetData(Me%Window%Right,                                   &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'RIGHT',                            &
                                 Default      =  99999.,                            &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR060'

                    call GetData(Me%Window%Top,                                     &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'TOP',                              &
                                 Default      =  99999.,                            &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR070'

                    call GetData(Me%Window%Bottom,                                  &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'BOTTOM',                           &
                                 Default      =  -99999.,                           &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR080'

                    call GetData(Me%Window%MinimumValue,                            &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'MINIMUM',                          &
                                 Default      = -9999999.,                          &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR090'

                    call GetData(Me%Window%MaximumValue,                            &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'MAXIMUM',                          &
                                 Default      = 9999999.,                           &
                                 ClientModule = 'ConvertToXYZ',                     &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR100'


                    call ConvertFile

                    call ClearMe

                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if(STAT_CALL .ne. SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR110'
                        
                    exit do1    !No more blocks

                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'StartConverting - ConvertToXYZ - ERR120'
                    
            end if if1
        end do do1


        call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'StartConverting - ConvertToXYZ - ERR130'


    end subroutine StartConverting

    !------------------------------------------------------------------
        
    subroutine ConvertFile

        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL

        !Begin---------------------------------------------------------

        select case(Me%FileType)

            case(NASA)
                
                call ConvertNASA  (Me%Window, Me%ObjEnterData, STAT = STAT_CALL)

            case(Etopo2)

                call ConvertEtopo2(Me%Window, Me%ObjEnterData, STAT = STAT_CALL)
            
            case(Etopo5)

                call ConvertEtopo5(Me%Window, Me%ObjEnterData, STAT = STAT_CALL)


            case(GEBCO)

                call ConvertGEBCO (Me%Window, Me%ObjEnterData, STAT = STAT_CALL)


            case(NOAA_ShoreLine)

                call ConvertNOAA_ShoreLine(Me%ObjEnterData, STAT = STAT_CALL)

            case(ASCII)

                call ConvertASCII (Me%Window, Me%ObjEnterData, STAT = STAT_CALL)

            case(SRTM30)

                call ConvertSRTM30 (Me%Window, Me%ObjEnterData, STAT = STAT_CALL)

            case default
                
                stop 'ConvertFile - ConvertToXYZ - ERR01'

        end select


    end subroutine ConvertFile

    !------------------------------------------------------------------

    subroutine ClearMe

        Me%FileType             = null_int
        Me%Window%Left          = null_real
        Me%Window%Right         = null_real
        Me%Window%Top           = null_real
        Me%Window%Left          = null_real
        Me%Input%FileName       = null_str
        Me%XYZOutput%FileName   = null_str

    end subroutine ClearMe


end program ConvertToXYZ






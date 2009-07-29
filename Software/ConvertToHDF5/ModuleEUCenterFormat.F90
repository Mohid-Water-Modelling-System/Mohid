!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : EUCenterFormat
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to convert EUCenterFormat files into HDF5 format
!
!------------------------------------------------------------------------------


Module ModuleEUCenterFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleTime

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertEUCenterFormat
    private ::      AllocateVariables
    private ::      ReadOptions
    private ::      OpenAndReadEUCenterFile
    private ::          EUCenterFieldName
    private ::          EUCenterFieldUnits
    private ::      ConstructEUCenterGrid
    private ::      OutputFields
    private ::          Open_HDF5_OutPut_File

    !Parameters----------------------------------------------------------------
    
    !Properties ID-------------------------------------------------------------
    integer, parameter                          :: SensibleHeat_ID         = 146
    integer, parameter                          :: LatentHeat_ID           = 147
    integer, parameter                          :: WindStressX_ID          = 180
    integer, parameter                          :: WindStressY_ID          = 181
    integer, parameter                          :: SolarRadiation_ID       = 176
    integer, parameter                          :: NetLongWaveRadiation_ID = 177
    integer, parameter                          :: Precipitation_ID        = 228
    
    !Grid---------------------------------------------------------------------
    integer, parameter                          :: EUCenterGrid_ILB      = 1
    integer, parameter                          :: EUCenterGrid_IUB      = 57
    integer, parameter                          :: EUCenterGrid_JLB      = 1
    integer, parameter                          :: EUCenterGrid_JUB      = 63
    real,    parameter                          :: EUCenterGrid_OriginX  = -26.1565
    real,    parameter                          :: EUCenterGrid_OriginY  = 28.9405
    real,    parameter                          :: EUCenterGrid_DX       = 0.563
    real,    parameter                          :: EUCenterGrid_DY       = 0.563
    
    !Types---------------------------------------------------------------------
    
    private :: T_EUCenterField
    type       T_EUCenterField
        character(len=StringLength)             :: Name
        character(len=StringLength)             :: Units
        integer                                 :: IDNumber
        type(T_Time)                            :: Date
        real, dimension(:,:),       pointer     :: Scalar
        type(T_EUCenterField),      pointer     :: Next
    end type  T_EUCenterField
    

    private :: T_EUCenterFormat
    type       T_EUCenterFormat
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjEUCenterGrid      = 0
        integer                                 :: ObjNewGrid           = 0
        integer                                 :: ObjHDF5              = 0
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: EUCenterGridFileName
        character(len=PathLength)               :: OutputFileName
        integer                                 :: Unit
        integer                                 :: NumberOfFields
        integer, dimension(:,:,:),  pointer     :: WaterPoints3D
        integer, dimension(:,:  ),  pointer     :: Bathymetry
        type(T_Size2D)                          :: Size
        type(T_EUCenterField),      pointer     :: FirstEUCenterField
    end type  T_EUCenterFormat

    type(T_EUCenterFormat),         pointer     :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertEUCenterFormat(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !External----------------------------------------------------------------
        
        
        !Local-------------------------------------------------------------------
        integer                                         :: nUsers

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions

        call ConstructEUCenterGrid

        call AllocateVariables

        call OpenAndReadEUCenterFile

        call OutputFields

        call KillEUCenterFormat

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'ConvertEUCenterFormat - ModuleEUCenterFormat - ERR01' 

        STAT = SUCCESS_

        deallocate(Me)
        nullify   (Me)



    end subroutine ConvertEUCenterFormat

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
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleEUCenterFormat - ERR02'
        end if

        call GetData(Me%EUCenterGridFileName,                           &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'EUCENTER_GRID_FILENAME',           &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR03'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleEUCenterFormat - ERR04'
        end if

        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR05'

        call GetData(Me%NumberOfFields,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NUMBER_OF_FIELDS',                 &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR06'


    end subroutine ReadOptions

    !------------------------------------------------------------------------

    subroutine AllocateVariables


        allocate(Me%WaterPoints3D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB, 0:2))
        Me%WaterPoints3D = 1

        allocate(Me%Bathymetry(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        Me%Bathymetry = 0.


    end subroutine AllocateVariables

    !------------------------------------------------------------------------

    subroutine OpenAndReadEUCenterFile

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist
        real, dimension(6)                          :: TimeVector
        real                                        :: Dummy
        type(T_EUCenterField), pointer              :: NewEUCenterField
        integer                                     :: Field, i, j

        !Begin-----------------------------------------------------------------
       
        nullify(Me%FirstEUCenterField)
        nullify(NewEUCenterField)

        !Verifies if file exists
        inquire(FILE = Me%FileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'EU Center File does not exist'
            stop 'OpenAndReadEUCenterFile - ModuleEUCenterFormat - ERR01'
        endif
        
        write(*,*)
        write(*,*)'Reading file: '//trim(Me%FileName)

        !opens file
        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadEUCenterFile - ModuleEUCenterFormat - ERR02'

        open(Me%Unit, File = Me%FileName, Form='FORMATTED', &
             status = 'OLD', Action = 'READ', IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadEUCenterFile - ModuleEUCenterFormat - ERR03'

        write(*,*) 'Reading fields. Please wait...'

        do Field = 1, Me%NumberOfFields

            call AddEUCenterField(NewEUCenterField)

            nullify (NewEUCenterField%Scalar)
            allocate(NewEUCenterField%Scalar(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

            read(Me%Unit, *)NewEUCenterField%IDNumber
            read(Me%Unit, *)Dummy, TimeVector(1),TimeVector(2),TimeVector(3),TimeVector(4)
        

            TimeVector(1) = TimeVector(1) + 1900.
            TimeVector(5) = 0.
            TimeVector(6) = 0.

            NewEUCenterField%Name   = EUCenterFieldName  (NewEUCenterField%IDNumber)
            NewEUCenterField%Units  = EUCenterFieldUnits (NewEUCenterField%IDNumber)


            call SetDate(NewEUCenterField%Date, TimeVector(1),TimeVector(2),TimeVector(3),&
                                                TimeVector(4),TimeVector(5),TimeVector(6))

            do i = Me%Size%IUB, Me%Size%ILB, -1
            do j = Me%Size%JLB, Me%Size%JUB
                
                read(Me%Unit,*) Dummy, Dummy, NewEUCenterField%Scalar(i,j)

            enddo
            enddo

            NewEUCenterField%Scalar = NewEUCenterField%Scalar/86400.

        end do


        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadEUCenterFile - ModuleEUCenterFormat - ERR04'

        close(Me%Unit, IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'OpenAndReadEUCenterFile - ModuleEUCenterFormat - ERR05'



    end subroutine OpenAndReadEUCenterFile

    
    !------------------------------------------------------------------------


    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, OutputNumber
        type (T_EUCenterField), pointer                 :: EUCenterField
        type(T_Time)                                    :: CurrentDate

        !Begin-----------------------------------------------------------------
        
        call Open_HDF5_OutPut_File

        OutputNumber = 1
        CurrentDate  = Me%FirstEUCenterField%Date
        
        write(*,*)'Writing field number: ', OutputNumber

        call ExtractDate   (CurrentDate,                                &
                            AuxTime(1), AuxTime(2), AuxTime(3),         &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR01'


        call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                             "Time", "YYYY/MM/DD HH:MM:SS",             &
                             Array1D = TimePtr,                         &
                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR02'

        EUCenterField => Me%FirstEUCenterField

        do while(associated(EUCenterField))

            if(EUCenterField%Date .gt. CurrentDate)then

                CurrentDate = EUCenterField%Date

                OutputNumber = OutputNumber + 1

                write(*,*)'Writing field number: ', OutputNumber

                call ExtractDate   (CurrentDate,                                &
                                    AuxTime(1), AuxTime(2), AuxTime(3),         &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime
            
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR03'
                

                call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",             &
                                     Array1D = TimePtr,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR04'

            end if


            if(EUCenterField%Date .eq. CurrentDate)then

                call HDF5SetLimits  (Me%ObjHDF5, Me%Size%ILB, Me%Size%IUB,&
                                     Me%Size%JLB, Me%Size%JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR05'

                call HDF5WriteData(Me%ObjHDF5,                                      &
                                   "/Results/"//EUCenterField%Name,                 &
                                   EUCenterField%Name,                              &
                                   EUCenterField%Units,                             &
                                   Array2D      = EUCenterField%Scalar,             &
                                   OutputNumber = OutPutNumber,                     &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR06'

            end if

            EUCenterField => EUCenterField%Next

        end do

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR07'

    end subroutine OutputFields


    !------------------------------------------------------------------------


    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleEUCenterFormat - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%Size%ILB, Me%Size%IUB,&
                             Me%Size%JLB, Me%Size%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleEUCenterFormat - ERR02'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleEUCenterFormat - ERR03'            

        call WriteHorizontalGrid (Me%ObjEUCenterGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleEUCenterFormat - ERR04'


        call HDF5SetLimits  (Me%ObjHDF5, Me%Size%ILB, Me%Size%IUB,&
                             Me%Size%JLB, Me%Size%JUB, 1, 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleEUCenterFormat - ERR05'


        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints", "-",    &
                              Array3D = Me%WaterPoints3D,                 &
                              STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleEUCenterFormat - ERR06'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleEUCenterFormat - ERR07'

    end subroutine Open_HDF5_OutPut_File

    
    !------------------------------------------------------------------------
    
    
    subroutine ConstructEUCenterGrid
        
        !External--------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Local-----------------------------------------------------------------
        logical                                     :: exist

        !Begin-----------------------------------------------------------------

        
        !Verifies if file exists
        inquire(FILE = Me%EUCenterGridFileName, EXIST = exist)
        if (.not. exist) then
            write(*,*)'Grid file does not exist'
            stop 'ConstructEUCenterGrid -  ModuleEUCenterFormat - ERR01'
        endif

        call ConstructHorizontalGrid(Me%ObjEUCenterGrid, Me%EUCenterGridFileName, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructEUCenterGrid -  ModuleEUCenterFormat - ERR02'

        call GetHorizontalGridSize(Me%ObjEUCenterGrid, WorkSize = Me%Size, STAT = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'ConstructEUCenterGrid -  ModuleEUCenterFormat - ERR03'


    end subroutine ConstructEUCenterGrid
    
    


    !------------------------------------------------------------------------

    character(len=StringLength) function EUCenterFieldName(FieldID)

        !Arguments---------------------------------------------------------------
        integer,                    intent(IN )         :: FieldID

        !------------------------------------------------------------------------

        select case(FieldID)

            case(SensibleHeat_ID)
                
                EUCenterFieldName = GetPropertyName (SensibleHeat_)

            case(LatentHeat_ID)

                EUCenterFieldName = GetPropertyName (LatentHeat_)
                
            case(WindStressX_ID)
            
                EUCenterFieldName = GetPropertyName (WindStressX_)
            
            case(WindStressY_ID)

                EUCenterFieldName = GetPropertyName (WindStressY_)

            case(SolarRadiation_ID)
            
                EUCenterFieldName = GetPropertyName (SolarRadiation_)

            case(NetLongWaveRadiation_ID)

                EUCenterFieldName = GetPropertyName (NetLongWaveRadiation_)

            case(Precipitation_ID)
                
                EUCenterFieldName = GetPropertyName (Precipitation_)

            case default

                write(*,*)'Unknown EUCenter field name'
                stop 'EUCenterFieldName -  ModuleEUCenterFormat - ERR01'


        end select


    end function EUCenterFieldName


    character(len=StringLength) function EUCenterFieldUnits(FieldID)

        !Arguments---------------------------------------------------------------
        integer,                    intent(IN )         :: FieldID

        !------------------------------------------------------------------------

        select case(FieldID)

            case(SensibleHeat_ID)
                
                EUCenterFieldUnits = 'W/m2'

            case(LatentHeat_ID)

                EUCenterFieldUnits = 'W/m2'
                
            case(WindStressX_ID)
            
                EUCenterFieldUnits = 'N/m2'
            
            case(WindStressY_ID)

                EUCenterFieldUnits = 'N/m2'

            case(SolarRadiation_ID)
            
                EUCenterFieldUnits = 'W/m2'

            case(NetLongWaveRadiation_ID)

                EUCenterFieldUnits = 'W/m2'

            case(Precipitation_ID)
                
                EUCenterFieldUnits = 'm/s'

            case default

                write(*,*)'Unknown EUCenter field units'
                stop 'EUCenterFieldUnits -  ModuleEUCenterFormat - ERR01'


        end select


    end function EUCenterFieldUnits


    !------------------------------------------------------------------------


    subroutine AddEUCenterField (ObjEUCenterField)

        !Arguments-------------------------------------------------------------
        type (T_EUCenterField), pointer                   :: ObjEUCenterField

        !Local-----------------------------------------------------------------
        type (T_EUCenterField), pointer                   :: NewEUCenterField
        type (T_EUCenterField), pointer                   :: PreviousEUCenterField
        
        !Begin-----------------------------------------------------------------
        
        !Allocates new instance
        allocate (NewEUCenterField)
        nullify  (NewEUCenterField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(Me%FirstEUCenterField)) then
            Me%FirstEUCenterField      => NewEUCenterField
            ObjEUCenterField           => NewEUCenterField
        else
            PreviousEUCenterField      => Me%FirstEUCenterField
            ObjEUCenterField           => Me%FirstEUCenterField%Next
            do while (associated(ObjEUCenterField))
                PreviousEUCenterField  => ObjEUCenterField
                ObjEUCenterField       => ObjEUCenterField%Next
            enddo
            ObjEUCenterField           => NewEUCenterField
            PreviousEUCenterField%Next => NewEUCenterField
        endif


    end subroutine AddEUCenterField

    !------------------------------------------------------------------------

    subroutine KillEUCenterFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjEUCenterGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillEUCenterFormat - ModuleEUCenterFormat - ERR01'
        
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillEUCenterFormat - ModuleEUCenterFormat - ERR02'

    
    end subroutine KillEUCenterFormat


    !------------------------------------------------------------------------
 
end module ModuleEUCenterFormat










!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : NETCDF
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2006
! REVISION      : Luis Fernandes - v4.0
! DESCRIPTION   : Module to Write/Read NETCDF Files
!                 See more on http://www.mohid.com/wiki/index.php?title=ModuleNETCDF
!
!------------------------------------------------------------------------------

Module ModuleNETCDF

    use ModuleGlobalData
    ! Manages NetCDF files
#ifdef _USE_NIX
    use netcdf
#else
    use netcdf90
#endif

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructNETCDF
    private ::      AllocateInstance

    !Selector
    public   :: GetNCDFFileAccess
                     
    
    !Modifier
    public  :: NETCDFWriteData
    public  :: NETCDFWriteHeader
    public  :: NETCDFWriteTime
    public  :: NETCDFSetDimensions
    public  :: NETCDFWriteLatLon
    public  :: NETCDFWriteVert
    !public  :: NETCDFReadData
    private :: NETCDFWriteAttributes

    !Destructor
    public  :: KillNETCDF                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjNETCDF 
    
    !Interfaces----------------------------------------------------------------
    interface NETCDFWriteData
        module procedure NETCDFWriteDataR4_2D
        module procedure NETCDFWriteDataR4_3D
        module procedure NETCDFWriteDataR8_2D
        module procedure NETCDFWriteDataR8_3D
        module procedure NETCDFWriteDataI4_2D
        module procedure NETCDFWriteDataI4_3D
    end interface
    
    interface NETCDFWriteLatLon    
        module procedure NETCDFWriteLatLon_1D
        module procedure NETCDFWriteLatLon_2D            
    end interface

    !interface NETCDFReadData
        !module procedure NETCDFReadDataR4_1D
        !module procedure NETCDFReadDataR4_2D
        !module procedure NETCDFReadDataR4_3D
        !module procedure NETCDFReadDataR8_1D
        !module procedure NETCDFReadDataR8_2D
        !module procedure NETCDFReadDataR8_3D
        !module procedure NETCDFReadDataI4_1D
        !module procedure NETCDFReadDataI4_2D
        !module procedure NETCDFReadDataI4_3D
    !end interface
    
    !Parameters----------------------------------------------------------------
    integer, parameter                              :: NCDF_CREATE_     = 1
    integer, parameter                              :: NCDF_READ_       = 2
    integer, parameter                              :: NCDF_READWRITE_  = 3

    !nf90_byte   = 1; nf90_char   = 2; nf90_short  = 3
    !nf90_int    = 4; nf90_float  = 5; nf90_double = 6

    !Types---------------------------------------------------------------------
    type     T_ID
        character(len=StringLength)                 :: Name
        character(len=StringLength)                 :: LongName
        integer                                     :: Number
        character(len=StringLength)                 :: Units
    end type T_ID
   
    type     T_Attribute
        type(T_ID)                                  :: ID
    end type T_Attribute


    type     T_Dimension
        type(T_ID)                                  :: ID
        integer                                     :: LB
        integer                                     :: UB
        integer                                     :: VarID
    end type T_Dimension

    type T_NETCDF
        character(len=PathLength)                   :: FileName
        character(len=StringLength)                 :: Title            
        character(len=StringLength)                 :: Convention       
        character(len=StringLength)                 :: Version          
        character(len=StringLength)                 :: History          
        character(len=StringLength)                 :: Source           
        character(len=StringLength)                 :: Institution      
        character(len=StringLength)                 :: References       

        type(T_Dimension), dimension(6)             :: Dims

        integer                                     :: Idate
        integer                                     :: nAttributes      = null_int     
        integer                                     :: nVariables       = null_int
        integer                                     :: nDimensions      = null_int  
        integer                                     :: ncid
        integer                                     :: InstanceID
        type(T_NETCDF), pointer                     :: Next
    end type T_NETCDF

    !Global Module Variables
    type (T_NETCDF), pointer                        :: FirstObjNETCDF
    type (T_NETCDF), pointer                        :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructNETCDF(ObjNCDFID, FileName, Access, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjNCDFID 
        character(len=*)                                :: FileName
        integer                                         :: Access
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_
        integer                                         :: STAT_CALL
        logical                                         :: exist

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mNETCDF_)) then
            nullify (FirstObjNETCDF)
            call RegisterModule (mNETCDF_) 
        endif

        call Ready(ObjNCDFID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%FileName = trim(FileName)

            if      (Access == NCDF_CREATE_) then

                STAT_CALL = nf90_create(path=trim(FileName), cmode = NF90_CLOBBER , ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'ConstructNETCDF - ModuleNETCDF - ERR01'

                STAT_CALL = nf90_close(Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'ConstructNETCDF - ModuleNETCDF - ERR02'

                STAT_CALL = nf90_open(path = trim(Me%FileName), mode = nf90_write, ncid = Me%ncid)  
                if(STAT_CALL /= nf90_noerr) stop 'ConstructNETCDF - ModuleNETCDF - ERR03'

            elseif  (Access == NCDF_READ_) then

                inquire(file = Me%FileName, exist = exist)
                if (.not. exist) then
                    write(*,*)
                    write(*,*)'NETCDF file does not exist'
                    write(*,*)trim(Me%FileName)
                    stop 'ConstructNETCDF - ModuleNETCDF - ERR04'
                endif

                stat = nf90_open(Me%FileName, nf90_nowrite, Me%ncid)
                if (stat /= nf90_noerr) stop 'ConstructNETCDF - ModuleNETCDF - ERR05'

                stat = nf90_inquire(Me%ncid, Me%nDimensions, Me%nVariables, Me%nAttributes)
                if (stat /= nf90_noerr) stop 'ConstructNETCDF - ModuleNETCDF - ERR06'

            elseif  (Access == NCDF_READWRITE_) then


            else

                stop 'ModuleNETCDF - ConstructNETCDF - ERR07' 

            endif


            !Returns ID
            ObjNCDFID  = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleNETCDF - ConstructNETCDF - ERR99' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructNETCDF
 
    !--------------------------------------------------------------------------
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_NETCDF), pointer                         :: NewObjNETCDF
        type (T_NETCDF), pointer                         :: PreviousObjNETCDF


        !Allocates new instance
        allocate (NewObjNETCDF)
        nullify  (NewObjNETCDF%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjNETCDF)) then
            FirstObjNETCDF          => NewObjNETCDF
            Me                      => NewObjNETCDF
        else
            PreviousObjNETCDF       => FirstObjNETCDF
            Me                      => FirstObjNETCDF%Next
            do while (associated(Me))
                PreviousObjNETCDF   => Me
                Me                  => Me%Next
            enddo
            Me                      => NewObjNETCDF
            PreviousObjNETCDF%Next  => NewObjNETCDF
        endif

        Me%InstanceID = RegisterNewInstance (mNETCDF_)


    end subroutine AllocateInstance


    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine GetNCDFFileAccess (NCDF_CREATE, NCDF_READ, NCDF_READWRITE)

        !Arguments-------------------------------------------------------------
        integer, optional                               :: NCDF_CREATE
        integer, optional                               :: NCDF_READ
        integer, optional                               :: NCDF_READWRITE

        !Local-----------------------------------------------------------------

        if (present(NCDF_CREATE     )) NCDF_CREATE      = NCDF_CREATE_
        if (present(NCDF_READ       )) NCDF_READ        = NCDF_READ_
        if (present(NCDF_READWRITE  )) NCDF_READWRITE   = NCDF_READWRITE_

    end subroutine GetNCDFFileAccess
    
    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    !--------------------------------------------------------------------------

    subroutine NETCDFWriteHeader(NCDFID, Title, Convention, Version, History, &
                                iDate, Source, Institution, References, STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID
        character(len=*)                            :: Title, Convention
        character(len=*)                            :: Version, History
        integer                                     :: iDate
        character(len=*)                            :: Source, Institution
        character(len=*)                            :: References
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%Title        = trim(Title)
            Me%Convention   = trim(Convention )
            Me%Version      = trim(Version    )
            Me%History      = trim(History    )
            Me%Source       = trim(Source     )
            Me%Institution  = trim(Institution)
            Me%References   = trim(References )
            Me%iDate        = iDate
            
            !Enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR00' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'Title',             Me%Title)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR01' 
            
            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'Conventions',       Me%Convention)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR02' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'netcdf_version_id', Me%Version)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR03' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'history',           Me%History)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR04' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'date',              Me%iDate)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR05' 
            
            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'source',            Me%Source)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR06' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'institution',       Me%Institution)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR07' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'references',        Me%References) 
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR08' 

            !Exit definition mode
            STAT_CALL = nf90_enddef(Me%ncid) 
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR09' 

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteHeader

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteAttributes (VarID, LongName, StandardName, Units, Positive, &
                                             Calendar, ScaleFactor, FillValue,        &
                                             MissingValue, ValidMin, ValidMax,        &
                                             Maximum, Minimum, Add_Offset, Step, iFillValue)

        !Arguments-------------------------------------------------------------
        integer                                     :: VarID
        character(len=*), optional                  :: LongName
        character(len=*), optional                  :: StandardName
        character(len=*), optional                  :: Units
        character(len=*), optional                  :: Positive
        character(len=*), optional                  :: Calendar
        real            , optional                  :: ScaleFactor
        real            , optional                  :: FillValue
        integer         , optional                  :: iFillValue
        real            , optional                  :: MissingValue
        real            , optional                  :: ValidMin
        real            , optional                  :: ValidMax
        real            , optional                  :: Minimum
        real            , optional                  :: Maximum
        real            , optional                  :: Add_Offset
        real            , optional                  :: Step

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        real                                        :: LastMax, LastMin
       
        !Begin-----------------------------------------------------------------

        if(present(LongName))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'long_name',     trim(LongName))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR01' 
        endif

        if(present(StandardName))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'standard_name', trim(StandardName))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR02' 
        endif

        if(present(Units))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'units',         trim(Units))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR03' 
        endif

        if(present(ScaleFactor))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'scale_factor', ScaleFactor) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR04' 
        end if

        if(present(FillValue))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, '_FillValue', FillValue) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR05' 
        end if
        
        if(present(MissingValue))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'missing_value', MissingValue) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR06' 
        end if

        if(present(ValidMin))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'valid_min', ValidMin) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR07' 
        end if
        
        if(present(ValidMax))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'valid_max', ValidMax) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR08' 
        end if
        
         if(present(Add_Offset))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'add_offset',     Add_Offset)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR09' 
        endif

        if(present(Step))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'step', Step) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR10' 
        end if

        if(present(Maximum))then
  
            STAT_CALL = nf90_inquire_attribute(Me%ncid, VarID, 'maximum')

            if    (STAT_CALL == nf90_noerr)then
                STAT_CALL = nf90_get_att(Me%ncid, VarID, 'maximum', LastMax)

                if(Maximum > LastMax)then
                    STAT_CALL = nf90_put_att(Me%ncid, VarID, 'maximum', Maximum) 
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR11' 
                end if

            elseif(STAT_CALL == nf90_enotatt)then
                
                STAT_CALL = nf90_put_att(Me%ncid, VarID, 'maximum', Maximum) 
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR12' 

            elseif(STAT_CALL /= nf90_noerr)then 

                stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR13'

            end if

        end if

        if(present(Minimum))then


            STAT_CALL = nf90_inquire_attribute(Me%ncid, VarID, 'minimum')

            if    (STAT_CALL == nf90_noerr)then
                STAT_CALL = nf90_get_att(Me%ncid, VarID, 'minimum', LastMin)

                if(Minimum < LastMin)then
                    STAT_CALL = nf90_put_att(Me%ncid, VarID, 'minimum', Minimum) 
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR14' 
                end if

            elseif(STAT_CALL == nf90_enotatt)then
                
                STAT_CALL = nf90_put_att(Me%ncid, VarID, 'minimum', Minimum) 
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR15' 

            elseif(STAT_CALL /= nf90_noerr)then 

                stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR16'

            end if

        end if

        if(present(Positive))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'positive',     trim(Positive))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR17' 
        endif

        if(present(Calendar))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'calendar',     trim(Calendar))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR18' 
        endif

        if(present(iFillValue))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, '_FillValue', iFillValue) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR19' 
        end if


    end subroutine NETCDFWriteAttributes

    !--------------------------------------------------------------------------

    subroutine NETCDFSetDimensions (NCDFID, IUB, JUB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID
        integer                                     :: IUB
        integer                                     :: JUB
        integer                                     :: KUB
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            Me%Dims(1)%ID%Name = "lon"
            Me%Dims(1)%LB      = 1
            Me%Dims(1)%UB      = JUB

            Me%Dims(2)%ID%Name = "lat"
            Me%Dims(2)%LB      = 1
            Me%Dims(2)%UB      = IUB

            if(KUB > 0)then
                Me%Dims(3)%ID%Name = "depth"
                Me%Dims(3)%LB      = 1
                Me%Dims(3)%UB      = KUB
            endif

            Me%Dims(5)%ID%Name = "lon_staggered"
            Me%Dims(5)%LB      = 1
            Me%Dims(5)%UB      = JUB + 1

            Me%Dims(6)%ID%Name = "lat_staggered"
            Me%Dims(6)%LB      = 1
            Me%Dims(6)%UB      = IUB + 1

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR00'

            !define longitude as dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(1)%ID%Name), Me%Dims(1)%UB, Me%Dims(1)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR01'

            !define latitude as dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(2)%ID%Name), Me%Dims(2)%UB, Me%Dims(2)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR02'
            
            if(KUB > 0)then
                !define depth as dimension
                STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(3)%ID%Name), Me%Dims(3)%UB, Me%Dims(3)%ID%Number)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR03'
            end if

            !define lon_staggered as dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(5)%ID%Name), Me%Dims(5)%UB, Me%Dims(5)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR04'

            !define lat_staggered as dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(6)%ID%Name), Me%Dims(6)%UB, Me%Dims(6)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR05'

            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR06'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFSetDimensions

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR4_2D (NCDFID, Name, LongName, StandardName, Units,           &
                                     ValidMin, ValidMax, MinValue, MaxValue, MissingValue,  &
                                     Positive, OutputNumber, Array2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName, Positive
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        integer,                      optional          :: OutputNumber
        real(4), dimension(:, :) , pointer              :: Array2D
        integer,                      optional          :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(2)                           :: Dims2ID
        integer, dimension(3)                           :: Dims3ID

        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            if     (STAT_CALL == nf90_enotvar)then

                Dims2ID(1) = Me%Dims(1)%ID%Number  !x
                Dims2ID(2) = Me%Dims(2)%ID%Number  !y
            
                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR00'

                if(present(OutputNumber))then
                    Dims3ID(1:2) = Dims2ID(1:2)
                    Dims3ID(3)   = Me%Dims(4)%ID%Number  !time

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_float, Dims3ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR01'
                else
                    

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_float, Dims2ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR02'

                endif

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           LongName       = trim(LongName),             &
                                           StandardName   = trim(StandardName),         &
                                           Units          = trim(Units),                &
                                           FillValue      = FillValueReal,              &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue,               &
                                           Positive       = Positive)
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR03'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR04'

            end if
            
            if(present(OutputNumber))then
                
                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue)

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array2D,                       &
                                         start = (/ 1, 1, OutputNumber /),              &
                                         count = (/ Me%Dims(1)%UB, Me%Dims(2)%UB, 1 /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR05' 

            else

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array2D)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR06' 

            end if


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataR4_2D

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR4_3D (NCDFID, Name, LongName, StandardName, Units,            &
                                     ValidMin, ValidMax,  MinValue, MaxValue, MissingValue,  &
                                     Positive, OutputNumber, Array3D, STAT)
                                   
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName, Positive
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        integer,                      optional          :: OutputNumber
        real(4), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(3)                           :: Dims3ID
        integer, dimension(4)                           :: Dims4ID

        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            if     (STAT_CALL == nf90_enotvar)then
            
                Dims3ID(1) = Me%Dims(1)%ID%Number  !x
                Dims3ID(2) = Me%Dims(2)%ID%Number  !y
                Dims3ID(3) = Me%Dims(3)%ID%Number  !z                            
                
                if(present(OutputNumber))then
                
                    if (OutputNumber == 1) then

                        !enter definition mode
                        STAT_CALL = nf90_redef(ncid = Me%ncid)
                        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR10'
                    
                        Dims4ID(1:3) = Dims3ID(1:3)
                        Dims4ID(4)   = Me%Dims(4)%ID%Number  !time

                        STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_float, Dims4ID, VarID)
                        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR20'
                        
                        call NETCDFWriteAttributes(VarID          = VarID,                      &
                                                   LongName       = trim(LongName),             &
                                                   StandardName   = trim(StandardName),         &
                                                   Units          = trim(Units),                &
                                                   FillValue      = FillValueReal,              &
                                                   ValidMin       = ValidMin,                   &
                                                   ValidMax       = ValidMax,                   &
                                                   Minimum        = MinValue,                   &
                                                   Maximum        = MaxValue,                   &
                                                   MissingValue   = MissingValue,               &
                                                   Positive       = Positive)
                        !exit definition mode
                        STAT_CALL = nf90_enddef(ncid = Me%ncid)
                        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR30'

                        
                        
                   endif
                else

                    !enter definition mode
                    STAT_CALL = nf90_redef(ncid = Me%ncid)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR40'
                    

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_float, Dims3ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR50'

                    call NETCDFWriteAttributes(VarID          = VarID,                      &
                                               LongName       = trim(LongName),             &
                                               StandardName   = trim(StandardName),         &
                                               Units          = trim(Units),                &
                                               FillValue      = FillValueReal,              &
                                               ValidMin       = ValidMin,                   &
                                               ValidMax       = ValidMax,                   &
                                               Minimum        = MinValue,                   &
                                               Maximum        = MaxValue,                   &
                                               MissingValue   = MissingValue,               &
                                               Positive       = Positive)
                    !exit definition mode
                    STAT_CALL = nf90_enddef(ncid = Me%ncid)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR60'


                endif

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR70'

            end if
            
            if(present(OutputNumber))then

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue)

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array3D,                       &
                                         start = (/ 1, 1, 1, OutputNumber /),           &
                                         count = (/ Me%Dims(1)%UB, Me%Dims(2)%UB, Me%Dims(3)%UB, 1 /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR80' 

            else

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array3D)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_3D - ModuleNETCDF - ERR90' 

            end if


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataR4_3D

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR8_2D (NCDFID, Name, LongName, StandardName, Units,           &
                                     ValidMin, ValidMax, MinValue, MaxValue, MissingValue,  &
                                     OutputNumber, Array2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        integer,                      optional          :: OutputNumber
        real(8), dimension(:, :) ,    pointer           :: Array2D
        integer,                      optional          :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(2)                           :: Dims2ID
        integer, dimension(3)                           :: Dims3ID

        
        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            if     (STAT_CALL == nf90_enotvar)then
            
                Dims2ID(1) = Me%Dims(1)%ID%Number  !x
                Dims2ID(2) = Me%Dims(2)%ID%Number  !y

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR00'

                if(present(OutputNumber))then
                    Dims3ID(1:2) = Dims2ID(1:2)                    
                    Dims3ID(3  ) = Me%Dims(4)%ID%Number  !time

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_double, Dims3ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR01'
                else
                    
                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_double, Dims2ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR02'

                endif

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           LongName       = trim(LongName),             &
                                           StandardName   = trim(StandardName),         &
                                           Units          = trim(Units),                &
                                           FillValue      = FillValueReal,              &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue)

                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR03'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR04'

            end if
            
            if(present(OutputNumber))then

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue)

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array2D,                       &
                                         start = (/ 1, 1, OutputNumber /),              &
                                         count = (/ Me%Dims(1)%UB, Me%Dims(2)%UB, 1 /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR05' 

            else

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array2D)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_2D - ModuleNETCDF - ERR06' 

            end if


            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataR8_2D
    
    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR8_3D (NCDFID, Name, LongName, StandardName, Units,           &
                                     ValidMin, ValidMax,  MinValue, MaxValue, MissingValue, &
                                     OutputNumber, Array3D, STAT)
                                   
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        integer,                      optional          :: OutputNumber
        real(8), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(3)                           :: Dims3ID
        integer, dimension(4)                           :: Dims4ID

        
        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            Dims3ID(1) = Me%Dims(1)%ID%Number  !z
            Dims3ID(2) = Me%Dims(2)%ID%Number  !y
            Dims3ID(3) = Me%Dims(3)%ID%Number  !z
            
            if     (STAT_CALL == nf90_enotvar)then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR00'

                if(present(OutputNumber))then
                    Dims4ID(1:3) = Dims3ID(1:3)
                    Dims4ID(4)   = Me%Dims(4)%ID%Number  !time

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_double, Dims4ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR01'
                else
                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_double, Dims3ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR02'

                endif

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           LongName       = trim(LongName),             &
                                           StandardName   = trim(StandardName),         &
                                           Units          = trim(Units),                &
                                           FillValue      = FillValueReal,              &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue)
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR03'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR04'

            end if
        
            if(present(OutputNumber))then

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue)

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array3D,                       &
                                         start = (/ 1, 1, 1, OutputNumber /),           &
                                         count = (/ Me%Dims(1)%UB, Me%Dims(2)%UB, Me%Dims(3)%UB, 1 /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR05' 

            else

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array3D)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_3D - ModuleNETCDF - ERR06' 

            end if


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif
            
        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataR8_3D
    
    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataI4_2D (NCDFID, Name, LongName, StandardName, Units, &
                                     ValidMin, ValidMax,  MinValue, MaxValue,     &
                                     MissingValue, OutputNumber, Array2D, STAT)
                                   
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax
        real,             intent(IN), optional          :: MinValue, MaxValue, MissingValue
        integer,                      optional          :: OutputNumber
        integer, dimension(:, :), pointer               :: Array2D
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(2)                           :: Dims2ID
        integer, dimension(3)                           :: Dims3ID

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)

            Dims2ID(1) = Me%Dims(1)%ID%Number  !x
            Dims2ID(2) = Me%Dims(2)%ID%Number  !y
        
            if     (STAT_CALL == nf90_enotvar)then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR00'

                if(present(OutputNumber))then
                    
                    Dims3ID(1:2) = Dims2ID(1:2)                    
                    Dims3ID(3  ) = Me%Dims(4)%ID%Number  !time

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_int, Dims3ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR01'
                else
                    
                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_int, Dims2ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR02'

                endif

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           LongName       = trim(LongName),             &
                                           StandardName   = trim(StandardName),         &
                                           Units          = trim(Units),                &
                                           iFillValue     = FillValueInt,               &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue)
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR03'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR04'

            end if
            
            if(present(OutputNumber))then

            
                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue)

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array2D,                       &
                                         start = (/ 1, 1, OutputNumber /),              &
                                         count = (/ Me%Dims(1)%UB, Me%Dims(2)%UB, Me%Dims(3)%UB, 1 /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR05' 
            else
                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array2D)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_2D - ModuleNETCDF - ERR06' 
            end if


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataI4_2D
    
    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataI4_3D (NCDFID, Name, LongName, StandardName, Units, &
                                     ValidMin, ValidMax,  MinValue, MaxValue,     &
                                     MissingValue, OutputNumber, Array3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        integer,                      optional          :: OutputNumber
        integer, dimension(:, :, :),  pointer           :: Array3D
        integer,                      optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(3)                           :: Dims3ID
        integer, dimension(4)                           :: Dims4ID

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            Dims3ID(1) = Me%Dims(1)%ID%Number  !x
            Dims3ID(2) = Me%Dims(2)%ID%Number  !y
            Dims3ID(3) = Me%Dims(3)%ID%Number  !z
        
            if     (STAT_CALL == nf90_enotvar)then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR00'

                if(present(OutputNumber))then
                    Dims4ID(1:3) = Dims3ID(1:3)                    
                    Dims4ID(4  ) = Me%Dims(4)%ID%Number  !time

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_int, Dims4ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR01'
                else

                    STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_int, Dims3ID, VarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR02'

                endif

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           LongName       = trim(LongName),             &
                                           StandardName   = trim(StandardName),         &
                                           Units          = trim(Units),                &
                                           iFillValue     = FillValueInt,               &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue)
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR03'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR04'

            end if
            
            if(present(OutputNumber))then
            
                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue)

                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array3D,                       &
                                         start = (/ 1, 1, 1, OutputNumber /),           &
                                         count = (/ Me%Dims(1)%UB, Me%Dims(2)%UB, Me%Dims(3)%UB, 1 /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR05' 
            else
                STAT_CALL = nf90_put_var(Me%ncid, VarID, Array3D)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataI4_3D - ModuleNETCDF - ERR06' 
            end if

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataI4_3D
    
    !--------------------------------------------------------------------------

    subroutine NETCDFWriteTime(NCDFID, InitialDate, nInstants, Times, StartInstant, DefDimTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=19), intent(IN)                   :: InitialDate
        integer                                         :: nInstants
        real(8), dimension(:), pointer                  :: Times
        integer, optional                               :: StartInstant
        logical, optional                               :: DefDimTime
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer, dimension(1)                           :: TimeDimID

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%Dims(4)%ID%Name = "time"
            Me%Dims(4)%LB      = 1
            Me%Dims(4)%UB      = nInstants
            
            DefDimTime = .true.
            
            if (present(StartInstant)) then
                if (StartInstant > 1) DefDimTime = .false.
            endif
            
            if (DefDimTime) then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR00'

                !define time as dimension
                STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(4)%ID%Name), NF90_UNLIMITED, Me%Dims(4)%ID%Number)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR01'

                TimeDimID(1) = Me%Dims(4)%ID%Number

                !define time as variable 
                STAT_CALL = nf90_def_var(Me%ncid, trim(Me%Dims(4)%ID%Name), NF90_DOUBLE, TimeDimID, Me%Dims(4)%VarID)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR02'

                !write variable attributes
                call NETCDFWriteAttributes(Me%Dims(4)%VarID, Units          = "seconds since "//InitialDate, &
                                                             LongName       = "time",                        &
                                                             StandardName   = "time") 
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR03'
                
            endif                
            
            !write time
            if (present(StartInstant)) then
                STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(4)%VarID, Times,                  &
                                             start = (/ StartInstant /),                    &
                                             count = (/ nInstants    /))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR04'
            else
                STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(4)%VarID, Times)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR05'
            endif

            STAT_ = SUCCESS_
            
        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteTime

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteLatLon_1D(NCDFID, Lat, Lon, Lat_Stag, Lon_Stag, GeoCoordinates, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:  ), pointer                   :: Lat, Lon
        real, dimension(:  ), pointer                   :: Lat_Stag, Lon_Stag
        logical                                         :: GeoCoordinates
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR00'

            !define latitude as variable
            STAT_CALL = nf90_def_var(Me%ncid, "lat", NF90_FLOAT, Me%Dims(2)%ID%Number, Me%Dims(2)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR01'

            !define longitude as variable 
            STAT_CALL = nf90_def_var(Me%ncid, "lon", NF90_FLOAT, Me%Dims(1)%ID%Number, Me%Dims(1)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR02'

            !define latitude staggered as variable
            STAT_CALL = nf90_def_var(Me%ncid, "lat_staggered", NF90_FLOAT, Me%Dims(6)%ID%Number, Me%Dims(6)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR03'

            !define longitude staggered as variable 
            STAT_CALL = nf90_def_var(Me%ncid, "lon_staggered", NF90_FLOAT, Me%Dims(5)%ID%Number, Me%Dims(5)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR04'

            if(GeoCoordinates)then

                call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "latitude",             &
                                                             StandardName = "latitude",             &
                                                             Units        = "degrees_north",        &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -90.,                   &
                                                             ValidMax     = 90.,                    &
                                                             MissingValue = FillValueReal)
               
                call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "longitude",            &
                                                             StandardName = "longitude",            &
                                                             Units        = "degrees_east",         &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -180.,                  &
                                                             ValidMax     = 180.,                   &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(5)%VarID, LongName     = "latitude staggered",   &
                                                             StandardName = "latitude",             &
                                                             Units        = "degrees_north",        &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -90.,                   &
                                                             ValidMax     = 90.,                    &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(6)%VarID, LongName     = "longitude staggered",  &
                                                             StandardName = "longitude",            &
                                                             Units        = "degrees_east",         &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -180.,                  &
                                                             ValidMax     = 180.,                   &
                                                             MissingValue = FillValueReal)
            else

                call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "metric X coordinate",  &
                                                             StandardName = "longitude",            &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "metric Y coordinate",  &
                                                             StandardName = "latitude",             &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(5)%VarID, LongName     = "metric X coordinate staggered",  &
                                                             StandardName = "longitude",            &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(6)%VarID, LongName     = "metric Y coordinate staggered", &
                                                             StandardName = "latitude",             &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)
            end if

            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR05'

            !write lat
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(2)%VarID, Lat)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR06'

            !write lon
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(1)%VarID, Lon)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR07'

            !write lat staggered
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(6)%VarID, Lat_Stag)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR08'

            !write lon staggered
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(5)%VarID, Lon_Stag)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_1D - ModuleNETCDF - ERR09'


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteLatLon_1D

    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteLatLon_2D(NCDFID, Lat, Lon, Lat_Stag, Lon_Stag, SphericX, SphericY, GeoCoordinates, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:,:), pointer                   :: Lat, Lon
        real, dimension(:,:), pointer                   :: Lat_Stag, Lon_Stag
        real(8), dimension(:,:), pointer                :: SphericX, SphericY       
        logical                                         :: GeoCoordinates
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        integer, dimension(2)                           :: Dims2ID, Dims2IDStag
        integer                                         :: SphericXVarID, SphericYVarID
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR10'
            
            Dims2ID(1) = Me%Dims(2)%ID%Number  !y
            Dims2ID(2) = Me%Dims(1)%ID%Number  !x
            
            Dims2IDStag(1) = Me%Dims(6)%ID%Number  !y
            Dims2IDStag(2) = Me%Dims(5)%ID%Number  !x            
    

            !define latitude as variable
            STAT_CALL = nf90_def_var(Me%ncid, "lat", NF90_FLOAT, Dims2ID, Me%Dims(2)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR20'

            !define longitude as variable 
            STAT_CALL = nf90_def_var(Me%ncid, "lon", NF90_FLOAT, Dims2ID, Me%Dims(1)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR30'

            !define latitude staggered as variable
            STAT_CALL = nf90_def_var(Me%ncid, "lat_staggered", NF90_FLOAT, Dims2IDStag, Me%Dims(6)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR40'

            !define longitude staggered as variable 
            STAT_CALL = nf90_def_var(Me%ncid, "lon_staggered", NF90_FLOAT, Dims2IDStag, Me%Dims(5)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR50'

            !define spherical mercator coordinates adopted by Google and Bing X staggered as variable 
            if (associated(SphericX)) then
                STAT_CALL = nf90_def_var(Me%ncid, "googlemaps_x", NF90_FLOAT, Dims2IDStag, SphericXVarID)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR60'
            endif

            !define spherical mercator Y staggered as variable 
            if (associated(SphericY)) then
                STAT_CALL = nf90_def_var(Me%ncid, "googlemaps_y", NF90_FLOAT, Dims2IDStag, SphericYVarID)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR70'
            endif

            if(GeoCoordinates)then

                call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "latitude",             &
                                                             StandardName = "latitude",             &
                                                             Units        = "degrees_north",        &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -90.,                   &
                                                             ValidMax     = 90.,                    &
                                                             MissingValue = FillValueReal)
               
                call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "longitude",            &
                                                             StandardName = "longitude",            &
                                                             Units        = "degrees_east",         &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -180.,                  &
                                                             ValidMax     = 180.,                   &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(5)%VarID, LongName     = "latitude staggered",   &
                                                             StandardName = "latitude",             &
                                                             Units        = "degrees_north",        &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -90.,                   &
                                                             ValidMax     = 90.,                    &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(6)%VarID, LongName     = "longitude staggered",  &
                                                             StandardName = "longitude",            &
                                                             Units        = "degrees_east",         &
                                                             FillValue    = FillValueReal,          &
                                                             ValidMin     = -180.,                  &
                                                             ValidMax     = 180.,                   &
                                                             MissingValue = FillValueReal)


                if (associated(SphericX)) then
                                                             
                    call NETCDFWriteAttributes(SphericXVarID,    LongName     = "spherical mercator - google maps - x staggered",  &
                                                                 StandardName = "spherical mercator - google maps - x",            &
                                                                 Units        = "paper meters",         &
                                                                 FillValue    = FillValueReal,          &
                                                                 ValidMin     = -20037508.34,           &
                                                                 ValidMax     =  20037508.34,           &
                                                                 MissingValue = FillValueReal)
                endif
                
                if (associated(SphericY)) then
                                                             
                    call NETCDFWriteAttributes(SphericYVarID,    LongName     = "spherical mercator - google maps - y staggered",  &
                                                                 StandardName = "spherical mercator - google maps - y",            &
                                                                 Units        = "paper meters",         &
                                                                 FillValue    = FillValueReal,          &
                                                                 ValidMin     = -20037508.34,           &
                                                                 ValidMax     =  20037508.34,           &
                                                                 MissingValue = FillValueReal)
                endif                
                                                                                 
            else

                call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "metric X coordinate",  &
                                                             StandardName = "longitude",            &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "metric Y coordinate",  &
                                                             StandardName = "latitude",             &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(5)%VarID, LongName     = "metric X coordinate staggered",  &
                                                             StandardName = "longitude",            &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)

                call NETCDFWriteAttributes(Me%Dims(6)%VarID, LongName     = "metric Y coordinate staggered", &
                                                             StandardName = "latitude",             &
                                                             Units        = "m",                    &
                                                             FillValue    = FillValueReal,          &
                                                             MissingValue = FillValueReal)
            end if

            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR80'

            !write lat
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(2)%VarID, Lat)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR90'

            !write lon
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(1)%VarID, Lon)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR100'

            !write lat staggered
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(6)%VarID, Lat_Stag)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR110'

            !write lon staggered
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(5)%VarID, Lon_Stag)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR120'

            !write spheric X staggered
            if (associated(SphericX)) then
                STAT_CALL = nf90_put_var(Me%ncid, SphericXVarID, SphericX)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR130'
            endif
            
            !write spheric Y staggered
            if (associated(SphericY)) then
                STAT_CALL = nf90_put_var(Me%ncid, SphericYVarID, SphericY)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon_2D - ModuleNETCDF - ERR140'
            endif            

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteLatLon_2D

    !--------------------------------------------------------------------------    

    subroutine NETCDFWriteVert(NCDFID, Vert, VertCoordinate, OffSet, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:  ), pointer                   :: Vert
        logical                                         :: VertCoordinate
        real, intent(in), optional                      :: OffSet
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVert - ModuleNETCDF - ERR00'
            
            !define depth as variable
            STAT_CALL = nf90_def_var(Me%ncid, "depth", NF90_FLOAT, Me%Dims(3)%ID%Number, Me%Dims(3)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVert - ModuleNETCDF - ERR01'

            !Is it a sigma grid?
            if(VertCoordinate)then

                call NETCDFWriteAttributes(Me%Dims(3)%VarID, LongName     = "depth",             &
                                                             StandardName = "depth",             &
                                                             Units        = "",                  &
                                                             Positive     = "down",              &
                                                             FillValue    = FillValueReal,       &
                                                             ValidMin     = 0.,                  &
                                                             ValidMax     = 1.,                  & 
                                                             Add_Offset   = OffSet,              &
                                                             MissingValue = FillValueReal)
               
            !No? Then it must be a z-coordinate grid
            else

                call NETCDFWriteAttributes(Me%Dims(3)%VarID, LongName     = "depth",             &
                                                             StandardName = "depth",             &
                                                             Units        = "m",                 &
                                                             Positive     = "down",              &
                                                             FillValue    = FillValueReal,       &
                                                             ValidMin     = -10.,                &
                                                             ValidMax     = 10000.,              &
                                                             Add_Offset   = OffSet,              &
                                                             MissingValue = FillValueReal)

            end if

            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVert - ModuleNETCDF - ERR05'

            !write Vertical coordinate
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(3)%VarID, Vert)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVert - ModuleNETCDF - ERR06'

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteVert

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillNETCDF(ObjNCDFID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjNCDFID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers
        integer                             :: STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjNCDFID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mNETCDF_,  Me%InstanceID)

            if (nUsers == 0) then

                STAT_CALL =nf90_close(Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - KillNETCDF - ERR01' 

                !Deallocates Instance
                call DeallocateInstance ()

                ObjNCDFID = 0
                STAT_       = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillNETCDF
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_NETCDF), pointer          :: AuxObjNETCDF
        type (T_NETCDF), pointer          :: PreviousObjNETCDF

        !Updates pointers
        if (Me%InstanceID == FirstObjNETCDF%InstanceID) then
            FirstObjNETCDF          => FirstObjNETCDF%Next
        else
            PreviousObjNETCDF       => FirstObjNETCDF
            AuxObjNETCDF            => FirstObjNETCDF%Next
            do while (AuxObjNETCDF%InstanceID /= Me%InstanceID)
                PreviousObjNETCDF   => AuxObjNETCDF
                AuxObjNETCDF        => AuxObjNETCDF%Next
            enddo

            !Now update linked list
            PreviousObjNETCDF%Next  => AuxObjNETCDF%Next

        endif

        !Deallocates instance
        deallocate (Me)
        nullify    (Me) 

            
    end subroutine DeallocateInstance

    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEMENT MANAGEME

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !--------------------------------------------------------------------------

    subroutine Ready (ObjNETCDF_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjNETCDF_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjNETCDF_ID > 0) then
            call LocateObjNETCDF (ObjNETCDF_ID)
            ready_ = VerifyReadLock (mNETCDF_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjNETCDF (ObjNCDFID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjNCDFID

        !Local-----------------------------------------------------------------

        Me => FirstObjNETCDF
        do while (associated (Me))
            if (Me%InstanceID == ObjNCDFID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleNETCDF - LocateObjNETCDF - ERR01'

    end subroutine LocateObjNETCDF

    !--------------------------------------------------------------------------

end module ModuleNETCDF

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------









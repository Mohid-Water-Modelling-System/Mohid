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
    public  :: NETCDFReadTime
    public  :: NETCDFReadGrid2D    
    public  :: NETCDFSetDimensions
    public  :: NETCDFSet1D_Dimension
    public  :: NETCDFGetDimensions    
    public  :: NETCDFWriteLatLon
    public  :: NETCDFWriteLatLon1D    
    public  :: NETCDFWriteVert
    public  :: NETCDFWriteVertStag    
    public  :: NETCDFReadVert    
    public  :: NETCDFWithVert
    public  :: NETCDFReadData
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
        module procedure NETCDFWriteDataR4_1D                
        module procedure NETCDFWriteDataR8_1D        
        !module procedure NETCDFWriteDataI4_1D                
    end interface
    


    interface NETCDFReadData
        module procedure NETCDFReadDataR4_2D
        module procedure NETCDFReadDataR4_3D
        module procedure NETCDFReadDataR8_2D
        module procedure NETCDFReadDataR8_3D
        module procedure NETCDFReadDataI4_2D
        module procedure NETCDFReadDataI4_3D
    end interface
    
    !Parameters----------------------------------------------------------------
    integer, parameter                              :: NCDF_CREATE_     = 1
    integer, parameter                              :: NCDF_READ_       = 2
    integer, parameter                              :: NCDF_READWRITE_  = 3

    !nf90_byte   = 1; nf90_char   = 2; nf90_short  = 3
    !nf90_int    = 4; nf90_float  = 5; nf90_double = 6
    
    !Grid
    character(len=StringLength), parameter          :: Lat_Name        = "lat"
    character(len=StringLength), parameter          :: Lon_Name        = "lon"    
    character(len=StringLength), parameter          :: Lat_Stag_Name   = "lat_staggered"
    character(len=StringLength), parameter          :: Lon_Stag_Name   = "lon_staggered"    
    character(len=StringLength), parameter          :: gmaps_x_Name    = "googlemaps_x"
    character(len=StringLength), parameter          :: gmaps_y_Name    = "googlemaps_y"

    character(len=StringLength), parameter          :: depth_Name      = "depth"
    character(len=StringLength), parameter          :: depth_Stag_Name = "depth_staggered"

    
    !Time
    character(len=StringLength), parameter          :: Time_Name                = "time"    
    !Dimensions
    character(len=StringLength), parameter          :: Column_Name              = "column_c"
    character(len=StringLength), parameter          :: Line_Name                = "line_c"
    character(len=StringLength), parameter          :: Layer_Name               = "depth"
    !horizontal vertices dimension
    character(len=StringLength), parameter          :: HorizCell_Vertices_Names = "nvh"
    !layers vertices dimension
    character(len=StringLength), parameter          :: Layer_Vertices_Names     = "nvv"    
    
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

        type(T_Dimension), dimension(7)             :: Dims

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

    subroutine ConstructNETCDF(NCDFID, FileName, Access, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: NCDFID 
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

        call Ready(NCDFID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            Me%FileName = trim(FileName)

            if      (Access == NCDF_CREATE_) then

                STAT_CALL = nf90_create(path  = trim(FileName),                         &
!                                        cmode = NF90_CLOBBER,                           &
                                         cmode = NF90_HDF5,                           &
                                        ncid  = Me%ncid)
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
            NCDFID  = Me%InstanceID

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

    subroutine NETCDFWriteHeader(NCDFID, Title, Convention, Version, History,           &
                                iDate, Source, Institution, References,                 & 
                                geospatial_lat_min, geospatial_lat_max,                 &
                                geospatial_lon_min, geospatial_lon_max,                 & 
                                CoordSysBuilder, contact, field_type, bulletin_date,    &
                                bulletin_type, comment, MetadataAtt, MetadataLink,      & 
                                STAT)


        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID
        character(len=*)                            :: Title, Convention
        character(len=*)                            :: Version, History
        integer                                     :: iDate
        character(len=*)                            :: Source, Institution
        character(len=*)                            :: References
        real            , optional                  :: geospatial_lat_min
        real            , optional                  :: geospatial_lat_max
        real            , optional                  :: geospatial_lon_min
        real            , optional                  :: geospatial_lon_max
        character(len=*), optional                  :: CoordSysBuilder
        character(len=*), optional                  :: contact
        character(len=*), optional                  :: field_type
        character(len=*), optional                  :: bulletin_date        
        character(len=*), optional                  :: bulletin_type        
        character(len=*), optional                  :: comment        
        character(len=*), optional                  :: MetadataAtt
        character(len=*), optional                  :: MetadataLink
        integer         , optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%Title                = trim(Title)
            Me%Convention           = trim(Convention )
            Me%Version              = trim(Version    )
            Me%History              = trim(History    )
            Me%Source               = trim(Source     )
            Me%Institution          = trim(Institution)
            Me%References           = trim(References )
            Me%iDate                = iDate
        
            
            !Enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR10' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'Title',             Me%Title)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR20' 
            
            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'Conventions',       Me%Convention)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR30' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'netcdf_version_id', Me%Version)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR40' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'history',           Me%History)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR50' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'date',              Me%iDate)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR60' 
            
            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'source',            Me%Source)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR70' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'institution',       Me%Institution)  
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR80' 

            STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'references',        Me%References) 
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR90'
            
            if (present(geospatial_lat_min)) then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'geospatial_lat_min', geospatial_lat_min) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR100' 
            endif                
            
            
            if (present(geospatial_lat_max)) then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'geospatial_lat_max', geospatial_lat_max) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR110' 
            endif

            if (present(geospatial_lon_min)) then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'geospatial_lon_min', geospatial_lon_min) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR120'             
            endif
            
            if (present(geospatial_lon_max)) then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'geospatial_lon_max', geospatial_lon_max) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR130'             
            endif
            
            if (present(CoordSysBuilder))   then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, '_CoordSysBuilder',   CoordSysBuilder) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR140'         
            endif
            
            if (present(contact))           then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'contact',            contact) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR150'
            endif
            
            if (present(field_type))        then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'field_type',         field_type) 
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR160'            
            endif

            if (present(bulletin_date))     then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'bulletin_date',      bulletin_date)
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR170'                        
            endif
            
            if (present(bulletin_type))     then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'bulletin_type',      bulletin_type)
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR180'               
            endif
            
            if (present(comment))           then
                STAT_CALL = nf90_put_att(Me%ncid, nf90_global, 'comment',            comment)
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR190'                           
            endif
            
            if (present(MetadataAtt) .and. present(MetadataLink))   then
                if (trim(MetadataAtt) /= trim(null_str)) then
                    STAT_CALL = nf90_put_att(Me%ncid, nf90_global, trim(MetadataAtt), trim(MetadataLink))
                    if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR200'                           
                endif                    
            endif            

            !Exit definition mode
            STAT_CALL = nf90_enddef(Me%ncid) 
            if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - NETCDFWriteHeader - ERR210' 

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteHeader

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteAttributes (VarID, LongName, StandardName, Units, Positive,   &
                                             Calendar, ScaleFactor, FillValue,          &
                                             MissingValue, ValidMin, ValidMax,          &
                                             Maximum, Minimum, Add_Offset, Step,        &
                                             iFillValue, iMissingValue,                 &
                                             coordinates, bounds, axis, reference,      &
                                             CoordinateAxisType,  CoordinateZisPositive)
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
        integer         , optional                  :: iMissingValue        
        real            , optional                  :: ValidMin
        real            , optional                  :: ValidMax
        real            , optional                  :: Minimum
        real            , optional                  :: Maximum
        real            , optional                  :: Add_Offset
        real            , optional                  :: Step
        character(len=*), optional                  :: coordinates
        character(len=*), optional                  :: bounds
        character(len=*), optional                  :: axis
        character(len=*), optional                  :: reference 
        character(len=*), optional                  :: CoordinateAxisType
        character(len=*), optional                  :: CoordinateZisPositive               

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, xtype
        real                                        :: LastMax, LastMin
        real(4)                                     :: AuxR4
        real(8)                                     :: AuxR8        

       
        !Begin-----------------------------------------------------------------

        if(present(LongName))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'long_name',     trim(LongName))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR01' 
        endif

        if(present(StandardName))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'standard_name', trim(StandardName))
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR02' 
        endif

        if(present(Units)) then
            if (units/=null_str) then
                STAT_CALL = nf90_put_att(Me%ncid, VarID, 'units',         trim(Units))
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR03' 
            endif                
        endif

        if(present(ScaleFactor))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'scale_factor', ScaleFactor) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR04' 
        end if

        if(present(FillValue))then
            if(present(MissingValue))then
                FillValue = MissingValue
            end if            
            STAT_CALL = nf90_inquire_variable(ncid = Me%ncid, varid = VarID, xtype = xtype)        
            if (xtype == nf90_double) then
                !STAT_CALL = nf90_put_att(ncid = Me%ncid, varid = VarID, name = '_FillValue', values = NF90_FILL_DOUBLE) 
                AuxR8 = FillValue            
                STAT_CALL = nf90_put_att(ncid = Me%ncid, varid = VarID, name = '_FillValue', values = AuxR8)             
                if(STAT_CALL /= nf90_noerr) then
                    write(*,*) trim(NF90_STRERROR(STAT_CALL))
                    stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR05' 
                endif    
            elseif (xtype == nf90_float) then                
                !STAT_CALL = nf90_put_att(ncid = Me%ncid, varid = VarID, name = '_FillValue', values = NF90_FILL_FLOAT) 
                AuxR4 = FillValue
                STAT_CALL = nf90_put_att(ncid = Me%ncid, varid = VarID, name = '_FillValue', values = AuxR4)             
                if(STAT_CALL /= nf90_noerr) then
                    write(*,*) trim(NF90_STRERROR(STAT_CALL))    
                    stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR05a'                     
                endif                    
                
            endif
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
            if(present(iMissingValue))then
                iFillValue = iMissingValue
            end if              
            STAT_CALL = nf90_put_att(Me%ncid, VarID, '_FillValue', iFillValue) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR19' 
        end if
        
        if(present(iMissingValue))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'missing_Value', iMissingValue) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR20' 
        end if        


        if(present(bounds))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'bounds', bounds) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR30' 
        end if

        if(present(coordinates))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'coordinates', coordinates) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR40' 
        endif
        


        if(present(axis))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'axis', axis) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR50' 
        end if


        if(present(reference))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, 'reference', reference) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR60' 
        end if
        

        if(present(CoordinateAxisType))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, '_CoordinateAxisType', CoordinateAxisType) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR70' 
        end if
                
        if(present(CoordinateZisPositive))then
            STAT_CALL = nf90_put_att(Me%ncid, VarID, '_CoordinateZisPositive', CoordinateZisPositive) 
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteAttributes - ModuleNETCDF - ERR80' 
        end if
                
        

    end subroutine NETCDFWriteAttributes

    !--------------------------------------------------------------------------

    subroutine NETCDFSetDimensions (NCDFID, IUB, JUB, KUB, SimpleGrid, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID
        integer                                     :: IUB
        integer                                     :: JUB
        integer                                     :: KUB
        logical, optional                           :: SimpleGrid
        integer,          optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL
        logical                                     :: SimpleGrid_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(SimpleGrid)) then
                SimpleGrid_ = SimpleGrid
            else                
                SimpleGrid_ = .false. 
            endif


            if (SimpleGrid_) then            
                Me%Dims(1)%ID%Name = trim(Lon_Name)
            else
                Me%Dims(1)%ID%Name = trim(Column_Name)
            endif                    
            Me%Dims(1)%LB      = 1
            Me%Dims(1)%UB      = JUB

            if (SimpleGrid_) then            
                Me%Dims(2)%ID%Name = trim(Lat_Name)
            else
                Me%Dims(2)%ID%Name = trim(Line_Name)
            endif                    

            Me%Dims(2)%LB      = 1
            Me%Dims(2)%UB      = IUB

            if(KUB > 0)then
                Me%Dims(3)%ID%Name = trim(Layer_Name)
                Me%Dims(3)%LB      = 1
                Me%Dims(3)%UB      = KUB
                
                Me%Dims(6)%ID%Name = trim(Layer_Vertices_Names)
                Me%Dims(6)%LB      = 1
                Me%Dims(6)%UB      = 2             
            endif

            Me%Dims(5)%ID%Name = trim(HorizCell_Vertices_Names)
            Me%Dims(5)%LB      = 1
            Me%Dims(5)%UB      = 4

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR00'

            !define matrixes columns as a dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(1)%ID%Name), Me%Dims(1)%UB, Me%Dims(1)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR01'

            !define matrixes lines as a dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(2)%ID%Name), Me%Dims(2)%UB, Me%Dims(2)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR02'
            
            if (.not. SimpleGrid_) then  
            
            !define the four horizontal cell vertices as a dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(5)%ID%Name), Me%Dims(5)%UB, Me%Dims(5)%ID%Number)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR04'
            
            endif
            
            if(KUB > 0)then
                !define layers as dimension
                STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(3)%ID%Name), Me%Dims(3)%UB, Me%Dims(3)%ID%Number)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR03'
                
                if (.not. SimpleGrid_) then  
                !define the two layers vertices as a dimension                
                STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(6)%ID%Name), Me%Dims(6)%UB, Me%Dims(6)%ID%Number)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFSetDimensions - ModuleNETCDF - ERR03'
                endif
                
            end if

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

    !--------------------------------------------------------------------------

    subroutine NETCDFSet1D_Dimension (NCDFID, IUB, DimName, Dim1DID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID
        integer                                     :: IUB
        character(Len=*)                            :: DimName
        integer                                     :: Dim1DID
        integer,          optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSet1D_Dimension - ModuleNETCDF - ERR10'

            !define matrixes columns as a dimension
            STAT_CALL = nf90_def_dim(Me%ncid, trim(DimName), IUB, Dim1DID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSet1D_Dimension - ModuleNETCDF - ERR20'

            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFSet1D_Dimension - ModuleNETCDF - ERR30'
            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFSet1D_Dimension

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR4_2D (NCDFID, Name, LongName, StandardName, Units,           &
                                     ValidMin, ValidMax, MinValue, MaxValue, MissingValue,  &
                                     Positive, coordinates, OutputNumber, Array2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName, Positive
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        character(len=*), intent(IN), optional          :: coordinates        
        integer,                      optional          :: OutputNumber
        real(4), dimension(:, :) , pointer              :: Array2D
        integer,                      optional          :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: FillValue, MissingValue_
        character(len=StringLength)                     :: coordinates_
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(2)                           :: Dims2ID
        integer, dimension(3)                           :: Dims3ID

        !Begin-----------------------------------------------------------------
        
        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(coordinates)) then
                coordinates_=coordinates
            else
                coordinates_="lon lat"
            endif        
            
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
                
                FillValue       = FillValueReal
                MissingValue_   = FillValue

                if(present(MissingValue))then
                    FillValue       = MissingValue
                    MissingValue_   = MissingValue
                end if                
                

                call NETCDFWriteAttributes(VarID          = VarID,                      &
                                           LongName       = trim(LongName),             &
                                           StandardName   = trim(StandardName),         &
                                           Units          = trim(Units),                &
                                           FillValue      = FillValue,                  &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue_,              &
                                           Positive       = Positive,                   &
                                           coordinates    = "lon lat")
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
                if(STAT_CALL /= nf90_noerr) then
                    print *, trim(nf90_strerror(STAT_CALL))
                    stop 'NETCDFWriteDataR4_2D - ModuleNETCDF - ERR06' 
                endif
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
                                     Positive, coordinates, OutputNumber, Array3D, STAT)
                                   
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName, Positive
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        character(len=*), intent(IN), optional          :: coordinates
        integer,                      optional          :: OutputNumber
        real(4), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: FillValue, MissingValue_
        character(len=StringLength)                     :: coordinates_
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

            if (present(coordinates)) then
                coordinates_=coordinates
            else
                coordinates_="lon lat"
            endif        
            
            FillValue       = FillValueReal
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = MissingValue
                MissingValue_   = MissingValue
            end if                
                
        
            
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
                                                   FillValue      = FillValue,                  &
                                                   ValidMin       = ValidMin,                   &
                                                   ValidMax       = ValidMax,                   &
                                                   Minimum        = MinValue,                   &
                                                   Maximum        = MaxValue,                   &
                                                   MissingValue   = MissingValue_,              &
                                                   Positive       = Positive,                   &
                                                   coordinates    = "lon lat")
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
                                               FillValue      = FillValue,                  &
                                               ValidMin       = ValidMin,                   &
                                               ValidMax       = ValidMax,                   &
                                               Minimum        = MinValue,                   &
                                               Maximum        = MaxValue,                   &
                                               MissingValue   = MissingValue_,              &
                                               Positive       = Positive,                   &
                                               coordinates    = "lon lat")
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
                                     coordinates, OutputNumber, Array2D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        character(len=*), intent(IN), optional          :: coordinates        
        integer,                      optional          :: OutputNumber
        real(8), dimension(:, :) ,    pointer           :: Array2D
        integer,                      optional          :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: FillValue, MissingValue_
        character(len=StringLength)                     :: coordinates_
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
            
            if (present(coordinates)) then
                coordinates_=coordinates
            else
                coordinates_="lon lat"
            endif       
            
            FillValue       = FillValueReal
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = MissingValue
                MissingValue_   = MissingValue
            end if        
            
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
                                           FillValue      = FillValue,                  &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue_,              &
                                           coordinates    = "lon lat")

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
                                     coordinates, OutputNumber, Array3D, STAT)
                                   
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        character(len=*), intent(IN), optional          :: coordinates        
        integer,                      optional          :: OutputNumber
        real(8), dimension(:, :, :), pointer            :: Array3D
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        real                                            :: FillValue, MissingValue_
        character(len=StringLength)                     :: coordinates_        
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
            
            FillValue       = FillValueReal
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = MissingValue
                MissingValue_   = MissingValue
            end if        
            
            if (present(coordinates)) then
                coordinates_=coordinates
            else
                coordinates_="lon lat"
            endif        
            
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
                                           FillValue      = FillValue,                  &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           MissingValue   = MissingValue_,              &
                                           coordinates    = "lon lat")
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

    subroutine NETCDFWriteDataI4_2D (NCDFID, Name, LongName, StandardName, Units,       &
                                     ValidMin, ValidMax,  MinValue, MaxValue,           &
                                     MissingValue, coordinates, OutputNumber,           &
                                     Array2D, STAT)
    
    

                                   
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax
        real,             intent(IN), optional          :: MinValue, MaxValue, MissingValue
        character(len=*), intent(IN), optional          :: coordinates        
        integer,                      optional          :: OutputNumber
        integer, dimension(:, :), pointer               :: Array2D
        integer, optional                               :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: FillValue, MissingValue_
        character(len=StringLength)                     :: coordinates_        
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
            
            if (present(coordinates)) then
                coordinates_=coordinates
            else
                coordinates_="lon lat"
            endif            
            
            FillValue       = FillValueInt
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = int(MissingValue)
                MissingValue_   = int(MissingValue)
            end if        
        
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
                                           iFillValue     = FillValue,                  &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           iMissingValue  = MissingValue_,              &
                                           coordinates    = "lon lat")
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

    subroutine NETCDFWriteDataI4_3D (NCDFID, Name, LongName, StandardName, Units,       &
                                     ValidMin, ValidMax,  MinValue, MaxValue,           &
                                     MissingValue, coordinates, OutputNumber,           &
                                     Array3D, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        character(len=*), intent(IN), optional          :: coordinates        
        integer,                      optional          :: OutputNumber
        integer, dimension(:, :, :),  pointer           :: Array3D
        integer,                      optional          :: STAT
        
        !Local-----------------------------------------------------------------
        character(len=StringLength)                     :: coordinates_        
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(3)                           :: Dims3ID
        integer, dimension(4)                           :: Dims4ID
        integer                                         :: FillValue, MissingValue_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            Dims3ID(1) = Me%Dims(1)%ID%Number  !x
            Dims3ID(2) = Me%Dims(2)%ID%Number  !y
            Dims3ID(3) = Me%Dims(3)%ID%Number  !z

            if (present(coordinates)) then
                coordinates_=coordinates
            else
                coordinates_="lon lat"
            endif        
            
            FillValue       = FillValueInt
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = int(MissingValue)
                MissingValue_   = int(MissingValue)
            end if                    
            
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
                                           iFillValue     = FillValue,                  &
                                           ValidMin       = ValidMin,                   &
                                           ValidMax       = ValidMax,                   &
                                           Minimum        = MinValue,                   &
                                           Maximum        = MaxValue,                   &
                                           iMissingValue  = MissingValue_,              &
                                           coordinates    = coordinates_)
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

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR4_1D (NCDFID, Name, LongName, StandardName, Units,       &
                                     ValidMin, ValidMax,  MinValue, MaxValue,           &
                                     MissingValue, Array1D, DimID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        real(4), dimension(:),        pointer           :: Array1D
        integer, optional                               :: DimID
        integer,                      optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(1)                           :: Dims1ID
        integer                                         :: FillValue, MissingValue_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            if (present(DimID)) then
                Dims1ID(1) = DimID
            else
                Dims1ID(1) = Me%Dims(1)%ID%Number  ! ID
            endif

            FillValue       = FillValueInt
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = int(MissingValue)
                MissingValue_   = int(MissingValue)
            end if                    
            
            if     (STAT_CALL == nf90_enotvar)then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_1D - ModuleNETCDF - ERR10'

                STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_double, Dims1ID, VarID)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_1D - ModuleNETCDF - ERR20'
                
                !call NETCDFWriteAttributes(VarID          = VarID,                      &
                !                           LongName       = trim(LongName),             &
                !                           StandardName   = trim(StandardName),         &
                !                           Units          = trim(Units),                &
                !                           iFillValue     = FillValue,                  &
                !                           ValidMin       = ValidMin,                   &
                !                           ValidMax       = ValidMax,                   &
                !                           Minimum        = MinValue,                   &
                !                           Maximum        = MaxValue,                   &
                !                           iMissingValue  = MissingValue_)
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_1D - ModuleNETCDF - ERR30'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataR4_1D - ModuleNETCDF - ERR40'

            end if
            
            STAT_CALL = nf90_put_var(Me%ncid, VarID, Array1D)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR4_1D - ModuleNETCDF - ERR50' 

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataR4_1D
    
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteDataR8_1D (NCDFID, Name, LongName, StandardName, Units,       &
                                     ValidMin, ValidMax,  MinValue, MaxValue,           &
                                     MissingValue, Array1D, DimID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        character(len=*), intent(IN)                    :: Name
        character(len=*), intent(IN)                    :: Units
        character(len=*), intent(IN), optional          :: LongName, StandardName
        real,             intent(IN), optional          :: ValidMin, ValidMax, MissingValue
        real,             intent(IN), optional          :: MinValue, MaxValue
        real(8), dimension(:),        pointer           :: Array1D
        integer, optional                               :: DimID
        integer,                      optional          :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VarID
        integer, dimension(1)                           :: Dims1ID
        integer                                         :: FillValue, MissingValue_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            STAT_CALL= nf90_inq_varid (Me%ncid, trim(Name), VarID)
            
            if (present(DimID)) then
                Dims1ID(1) = DimID
            else
                Dims1ID(1) = Me%Dims(1)%ID%Number  ! ID
            endif

            FillValue       = FillValueInt
            MissingValue_   = FillValue

            if(present(MissingValue))then
                FillValue       = int(MissingValue)
                MissingValue_   = int(MissingValue)
            end if                    
            
            if     (STAT_CALL == nf90_enotvar)then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_1D - ModuleNETCDF - ERR10'

                STAT_CALL= nf90_def_var(Me%ncid, trim(Name), nf90_double, Dims1ID, VarID)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_1D - ModuleNETCDF - ERR20'
                
                !call NETCDFWriteAttributes(VarID          = VarID,                      &
                !                           LongName       = trim(LongName),             &
                !                           StandardName   = trim(StandardName),         &
                !                           Units          = trim(Units),                &
                !                           iFillValue     = FillValue,                  &
                !                           ValidMin       = ValidMin,                   &
                !                           ValidMax       = ValidMax,                   &
                !                           Minimum        = MinValue,                   &
                !                           Maximum        = MaxValue,                   &
                !                           iMissingValue  = MissingValue_)
                !exit definition mode
                STAT_CALL = nf90_enddef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_1D - ModuleNETCDF - ERR30'

            elseif(STAT_CALL /= nf90_noerr)then

                stop 'NETCDFWriteDataR8_1D - ModuleNETCDF - ERR40'

            end if
            
            STAT_CALL = nf90_put_var(Me%ncid, VarID, Array1D)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteDataR8_1D - ModuleNETCDF - ERR50' 

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFWriteDataR8_1D
    
    !--------------------------------------------------------------------------
    !
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
        logical                                         :: DefDimTime_
        integer, dimension(1)                           :: TimeDimID

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%Dims(4)%ID%Name = "time"
            Me%Dims(4)%LB      = 1
            Me%Dims(4)%UB      = nInstants
            
            if (present(DefDimTime))then
                DefDimTime_ = DefDimTime
            else
                DefDimTime_ =.true.
            endif
            
            if (present(StartInstant)) then
                if (StartInstant > 1) DefDimTime_ = .false.
            endif
            
            if (DefDimTime_) then

                !enter definition mode
                STAT_CALL = nf90_redef(ncid = Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR00'

                !define time as dimension
                STAT_CALL = nf90_def_dim(Me%ncid, trim(Me%Dims(4)%ID%Name), NF90_UNLIMITED, Me%Dims(4)%ID%Number)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR01'

                TimeDimID(1) = Me%Dims(4)%ID%Number

                !define time as variable 
                STAT_CALL = nf90_def_var(Me%ncid, trim(Time_Name), NF90_DOUBLE, TimeDimID, Me%Dims(4)%VarID)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteTime - ModuleNETCDF - ERR02'

                !write variable attributes
                call NETCDFWriteAttributes(Me%Dims(4)%VarID, Units          = "seconds since "//InitialDate, &
                                                             LongName       = trim(Time_Name),               &
                                                             StandardName   = trim(Time_Name),               &
                                                             Calendar       = "standard") 
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

    subroutine NETCDFReadTime(NCDFID, InitialDate, nInstants, Instants, TimeName, TimeDimName, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real,    dimension(6)                           :: InitialDate
        integer                                         :: nInstants
        real(8), dimension(:), pointer                  :: Instants
        character(*), intent(IN), optional              :: TimeName, TimeDimName        
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real                                            :: UnitsFactor
        character(len=StringLength)                     :: TimeName_, TimeDimName_
        integer                                         :: STAT_, ready_
        character(Len=StringLength)                     :: ref_date
        integer                                         :: n, status, dimid, i, tmax
        logical                                         :: ReadTime

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
            
            if (present(TimeName)) then
                TimeName_ = TimeName
            else
                TimeName_ = Time_Name
            endif
            
            if (present(TimeDimName)) then
                TimeDimName_ = TimeDimName
            else
                TimeDimName_ = Time_Name
            endif
                   
            status=NF90_INQ_DIMID(Me%ncid, trim(TimeDimName_),dimid)
            if (status /= nf90_noerr) stop 'NETCDFReadTime - ModuleNETCDF - ERR10'

            status=NF90_INQUIRE_DIMENSION(Me%ncid, dimid, len = nInstants)
            if (status /= nf90_noerr) stop 'NETCDFReadTime - ModuleNETCDF - ERR20'
            
            allocate(Instants(nInstants))

            status = nf90_inq_varid(Me%ncid, trim(TimeName_), n)
            if (status /= nf90_noerr) stop 'ReadTimeNetCDF - ModuleNetCDFCF_2_HDF5MOHID - ERR30'
            
            status = NF90_GET_VAR(Me%ncid,n,Instants)
            if (status /= nf90_noerr) stop 'NETCDFReadTime - ModuleNETCDF - ERR40'


            !CF convention 
            status=NF90_GET_ATT(Me%ncid,n,"units", ref_date)
            if (status /= nf90_noerr) stop 'NETCDFReadTime - ModuleNETCDF - ERR50'
            
            tmax = len_trim(ref_date)

            ReadTime =.false.
            
            UnitsFactor = 3600.

            do i=1,tmax-5
                if (ref_date(i:i+5)== "second") then
                    ReadTime =.true.
                    UnitsFactor = 1.
                    exit
                endif
            enddo
            
            do i=1,tmax-2
                if (ref_date(i:i+2)== "day") then
                    ReadTime =.true.
                    UnitsFactor = 86400.
                    exit
                endif
            enddo            

            do i=1,tmax-5
                if (ref_date(i:i+5)== "minute") then
                    ReadTime =.true.
                    UnitsFactor = 60.
                    exit
                endif
            enddo            


            do i=1,tmax-4            
                if (ref_date(i:i+4)== "since") then
                    ref_date = ref_date(i+5:tmax)
                    exit
                endif
                
            enddo

            ReadTime = .false.
            do i=1,len_trim(ref_date)

                if (ref_date(i:i) ==':') then
                    ReadTime = .true.
                endif

            enddo  

            do i=1,len_trim(ref_date)

                if (ref_date(i:i) =='_'.or.ref_date(i:i) ==':'.or. ref_date(i:i) =='-'&
                    .or. ref_date(i:i) =='Z'.or. ref_date(i:i) =='T') then
                    ref_date(i:i) = ' '
                endif

            enddo  
            
            ref_date(1:19) = trim(adjustl(ref_date))
            
            InitialDate(:) = 0.
            
            read(ref_date(1:4),*) InitialDate(1)    
            read(ref_date(6:7),*) InitialDate(2)                
            read(ref_date(9:10),*) InitialDate(3)
            if (ReadTime) then                            
                read(ref_date(12:13),*) InitialDate(4) 
                read(ref_date(15:16),*) InitialDate(5)                                                                   
                read(ref_date(18:19),*) InitialDate(6)  
            endif
            
            Me%Dims(4)%ID%Name = trim(TimeName_)
            Me%Dims(4)%LB      = 1
            Me%Dims(4)%UB      = nInstants            
            
            STAT_ = SUCCESS_
            
        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadTime


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteLineColumn()

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        real                                            :: FillValue, MissingValue
        real, dimension(:  ), pointer                   :: Line, Column
        integer                                         :: STAT_CALL, j, i
 
        !Begin-----------------------------------------------------------------

        !enter definition mode
        STAT_CALL = nf90_redef(ncid = Me%ncid)
        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLineColumn - ModuleNETCDF - ERR10'

        !define matrix lines as variable
        STAT_CALL = nf90_def_var(Me%ncid, trim(Line_Name), NF90_FLOAT, Me%Dims(2)%ID%Number, Me%Dims(2)%VarID)
        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLineColumn - ModuleNETCDF - ERR20'

        !define longitude as variable 
        STAT_CALL = nf90_def_var(Me%ncid, trim(Column_Name), NF90_FLOAT, Me%Dims(1)%ID%Number, Me%Dims(1)%VarID)
        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLineColumn - ModuleNETCDF - ERR30'

        FillValue    = FillValueReal
        MissingValue = FillValueReal
        
        call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "X",                &
                                                     !StandardName = "x direction in the canonical space",   &
                                                     Units        = " ",       &
!                                                     FillValue    = FillValue,          &
!                                                     MissingValue = MissingValue,       &
                                                     axis           ="X")

        call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "Y",                &
                                                     !StandardName = "y direction in the canonical space",   &
                                                     Units        = " ",       &
!                                                     FillValue    = FillValue,          &
!                                                     MissingValue = MissingValue,       &
                                                     axis           ="Y")
                                                     
        allocate(Column(Me%Dims(1)%LB:Me%Dims(1)%UB))
        allocate(Line  (Me%Dims(2)%LB:Me%Dims(2)%UB))  
        
        do j=Me%Dims(1)%LB,Me%Dims(1)%UB
            Column  (j)=j
        enddo        
        
        do i=Me%Dims(2)%LB,Me%Dims(2)%UB
            Line    (i)=i
        enddo            

        !exit definition mode
        STAT_CALL = nf90_enddef(ncid = Me%ncid)
        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLineColumn - ModuleNETCDF - ERR40'

        !write column indexes
        STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(1)%VarID, Column)
        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLineColumn - ModuleNETCDF - ERR50'

        !write line indexes
        STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(2)%VarID, Line)
        if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLineColumn - ModuleNETCDF - ERR60'

        deallocate(Column)
        deallocate(Line  ) 

    end subroutine NETCDFWriteLineColumn

    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteLatLon(NCDFID, Lat, Lon, Lat_Stag, Lon_Stag,                  &
                                 SphericX, SphericY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:,:), pointer                   :: Lat, Lon
        real, dimension(:,:), pointer                   :: Lat_Stag, Lon_Stag
        real(8), dimension(:,:), pointer, optional      :: SphericX, SphericY     
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:,:), pointer              :: AuxStag
        integer, dimension(2)                           :: Dims2ID
        integer, dimension(3)                           :: Dims3IDStag
        real                                            :: FillValue, MissingValue        
        integer                                         :: LatStagID, LonStagID
        integer                                         :: SphericXVarID, SphericYVarID
        integer                                         :: i, j, nv, di, dj
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR10'
            
            Dims2ID(1) = Me%Dims(1)%ID%Number  !x
            Dims2ID(2) = Me%Dims(2)%ID%Number  !y
            
            Dims3IDStag(1)  = Me%Dims(5)%ID%Number  ! nvh
            Dims3IDStag(2:3)= Dims2ID(1:2)  

            !define latitude as variable
            STAT_CALL = nf90_def_var(Me%ncid, trim(Lat_Name), NF90_FLOAT, Dims2ID, Me%Dims(2)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR20'

            !define longitude as variable 
            STAT_CALL = nf90_def_var(Me%ncid, trim(Lon_Name), NF90_FLOAT, Dims2ID, Me%Dims(1)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR30'
            
            !define latitude staggered as variable
            STAT_CALL = nf90_def_var(Me%ncid, trim(Lat_Stag_Name), NF90_FLOAT, Dims3IDStag, LatStagID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR40'

            !define longitude staggered as variable 
            STAT_CALL = nf90_def_var(Me%ncid, trim(Lon_Stag_Name), NF90_FLOAT, Dims3IDStag, LonStagID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR50'

            !define spherical mercator coordinates adopted by Google and Bing X staggered as variable 
            if (present(SphericX)) then
                if (associated(SphericX)) then
                    STAT_CALL = nf90_def_var(Me%ncid, trim(gmaps_x_Name), NF90_FLOAT, Dims3IDStag, SphericXVarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR60'
                endif                    
            endif

            !define spherical mercator Y staggered as variable 
            if (present(SphericY)) then
                if (associated(SphericY)) then
                    STAT_CALL = nf90_def_var(Me%ncid, trim(gmaps_y_Name), NF90_FLOAT, Dims3IDStag, SphericYVarID)
                    if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR70'
                endif                    
            endif
            
            FillValue       = FillValueReal
            MissingValue    = FillValueReal


            call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "latitude",             &
                                                         StandardName = "latitude",             &
                                                         Units        = "degrees_north",        &
!                                                         FillValue    = FillValue,              &
                                                         ValidMin     = -90.,                   &
                                                         ValidMax     = 90.,                    &
                                                         bounds       = trim(Lat_Stag_Name)) !,    &
!                                                         MissingValue = MissingValue)
           
            call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "longitude",            &
                                                         StandardName = "longitude",            &
                                                         Units        = "degrees_east",         &
!                                                         FillValue    = FillValue,              &
                                                         ValidMin     = -180.,                  &
                                                         ValidMax     = 180.,                   &
                                                         bounds       = trim(Lon_Stag_Name)) !,    &
!                                                         MissingValue = MissingValue)

            if (present(SphericX)) then
                if (associated(SphericX)) then
                                                         
                    call NETCDFWriteAttributes(SphericXVarID,    LongName     = "spherical mercator - google maps - x staggered",  &
                                                                 !StandardName = "spherical mercator - google maps - x",            &
                                                                 Units        = "meters",         &
!                                                                 FillValue    = FillValue,              &
                                                                 ValidMin     = -20037508.34,           &
                                                                 ValidMax     =  20037508.34)
!                                                                 MissingValue = MissingValue)
                endif                    
            endif
            
            if (present(SphericY)) then
                if (associated(SphericY)) then

                    call NETCDFWriteAttributes(SphericYVarID,    LongName     = "spherical mercator - google maps - y staggered",  &
                                                                 !StandardName = "spherical mercator - google maps - y",            &
                                                                 Units        = "meters",         &
!                                                                 FillValue    = FillValue,              & 
                                                                 ValidMin     = -20037508.34,           &
                                                                 ValidMax     =  20037508.34)
!                                                                 MissingValue = MissingValue)
                endif
            endif                
            
            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR80'

            !write lat
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(2)%VarID, Lat)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR90'

            !write lon
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(1)%VarID, Lon)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR100'
            
            allocate(AuxStag(Me%Dims(5)%LB:Me%Dims(5)%UB,Me%Dims(1)%LB:Me%Dims(1)%UB,Me%Dims(2)%LB:Me%Dims(2)%UB))
            
            do j = Me%Dims(1)%LB,Me%Dims(1)%UB
            do i = Me%Dims(2)%LB,Me%Dims(2)%UB
                do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                    if (nv==1) then
                        di=0;dj=0
                    endif
                    if (nv==2) then
                        di=0;dj=1
                    endif
                    if (nv==3) then
                        di=1;dj=1
                    endif
                    if (nv==4) then
                        di=1;dj=0
                    endif
                    AuxStag(nv,j,i) = Lat_Stag(j+dj,i+di)
                enddo
            enddo
            enddo

            !write lat staggered
            STAT_CALL = nf90_put_var(Me%ncid, LatStagID, AuxStag)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR110'

            do j = Me%Dims(1)%LB,Me%Dims(1)%UB
            do i = Me%Dims(2)%LB,Me%Dims(2)%UB
                do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                    if (nv==1) then
                        di=0;dj=0
                    endif
                    if (nv==2) then
                        di=0;dj=1
                    endif
                    if (nv==3) then
                        di=1;dj=1
                    endif
                    if (nv==4) then
                        di=1;dj=0
                    endif
                    AuxStag(nv,j,i) = Lon_Stag(j+dj,i+di)
                enddo
            enddo
            enddo

            !write lon staggered
            STAT_CALL = nf90_put_var(Me%ncid, LonStagID, AuxStag)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR120'
            
            !write spheric X staggered
            if (present   (SphericX)) then
            if (associated(SphericX)) then
            
                do j = Me%Dims(1)%LB,Me%Dims(1)%UB
                do i = Me%Dims(2)%LB,Me%Dims(2)%UB
                    do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                        if (nv==1) then
                            di=0;dj=0
                        endif
                        if (nv==2) then
                            di=0;dj=1
                        endif
                        if (nv==3) then
                            di=1;dj=1
                        endif
                        if (nv==4) then
                            di=1;dj=0
                        endif
                        AuxStag(nv,j,i) = SphericX(j+dj,i+di)
                    enddo
                enddo
                enddo
            
                STAT_CALL = nf90_put_var(Me%ncid, SphericXVarID, AuxStag)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR130'
            endif
            endif
            
            !write spheric Y staggered
            if (present   (SphericY)) then
            if (associated(SphericY)) then

                do j = Me%Dims(1)%LB,Me%Dims(1)%UB
                do i = Me%Dims(2)%LB,Me%Dims(2)%UB
                    do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                        if (nv==1) then
                            di=0;dj=0
                        endif
                        if (nv==2) then
                            di=0;dj=1
                        endif
                        if (nv==3) then
                            di=1;dj=1
                        endif
                        if (nv==4) then
                            di=1;dj=0
                        endif
                        AuxStag(nv,j,i) = SphericY(j+dj,i+di)
                    enddo
                enddo
                enddo            
            
                STAT_CALL = nf90_put_var(Me%ncid, SphericYVarID, AuxStag)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon - ModuleNETCDF - ERR140'
            
            endif    
            endif
            

            call NETCDFWriteLineColumn            

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteLatLon

    !--------------------------------------------------------------------------    

    !--------------------------------------------------------------------------

    subroutine NETCDFWriteLatLon1D(NCDFID, Lat, Lon, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:,:), pointer                   :: Lat, Lon
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:), pointer                     :: Lat1D, Lon1D
        real                                            :: FillValue, MissingValue
        !real                                            :: Scale, Aux1, Aux2, Aux3, Dif  
        integer                                         :: i, j !, nv, di, dj
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR10'

            !define matrix lines as variable
            STAT_CALL = nf90_def_var(Me%ncid, trim(Lat_Name), NF90_FLOAT, Me%Dims(2)%ID%Number, Me%Dims(2)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR20'

            !define longitude as variable 
            STAT_CALL = nf90_def_var(Me%ncid, trim(Lon_Name), NF90_FLOAT, Me%Dims(1)%ID%Number, Me%Dims(1)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR30'

            FillValue    = FillValueReal
            MissingValue = FillValueReal
            
            call NETCDFWriteAttributes(Me%Dims(1)%VarID, LongName     = "longitude",    &
                                                         StandardName = "longitude",    &
                                                         Units        = "degrees_east", &
                                                         ValidMin    = -180.,           &
                                                         ValidMax    =  180.,           &
                                                         axis           ="X",           &
                                                         reference    = "geographical coordinates, WGS84 projection")

            call NETCDFWriteAttributes(Me%Dims(2)%VarID, LongName     = "latitude",     &
                                                         StandardName = "latitude",     &
                                                         Units        = "degrees_north",&
                                                         ValidMin    = -90.,            &
                                                         ValidMax    =  90.,            &
                                                         axis           ="Y",           &
                                                         reference    = "geographical coordinates, WGS84 projection")
                                                         
            allocate(Lon1D(Me%Dims(1)%LB:Me%Dims(1)%UB))
            allocate(Lat1D(Me%Dims(2)%LB:Me%Dims(2)%UB))  
            
!            Scale = Lon(2,1) - Lon(1,1)
!            
!            if (Scale > 0.) then
!                Scale = int(log10(100. / Scale))
!                Scale = 10**(Scale +1)
!            else
!                stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR40'
!            endif
!            
!            do j=Me%Dims(1)%LB,Me%Dims(1)%UB
!                Aux1 = Lon(j,1) * Scale
!                Aux2 = FLOAT (INT(Aux1))
!                Dif  = Aux1 - Aux2
!                if (abs(Dif) > 0.5) then
!                    if (Aux2>0) then
!                        Aux3 = (Aux2 + 1.)/Scale
!                    else
!                        Aux3 = (Aux2 - 1.)/Scale
!                    endif                        
!                else                    
!                    Aux3 = (Aux2     )/Scale
!                endif
!                                    
!                Lon1D  (j)= Aux3
!            enddo 
!            
!            Scale = Lat(1,2) - Lat(1,1)
!            
!            if (Scale > 0.) then
!                Scale = int(log10(100. / Scale))
!                Scale = 10**(Scale +1)
!            else
!                stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR40'
!            endif                   
!            
!            do i=Me%Dims(2)%LB,Me%Dims(2)%UB
!                Aux1 = Lat(1,i) * Scale
!                Aux2 = FLOAT (INT(Aux1))
!                Dif  = Aux1 - Aux2
!                if (Dif > 0.5) then
!                    Aux3 = (Aux2 + 1.)/Scale
!                else                    
!                    Aux3 = (Aux2     )/Scale
!                endif
!                                    
!                Lat1D  (i)= Aux3
!            enddo        

            do j=Me%Dims(1)%LB,Me%Dims(1)%UB
                Lon1D  (j)= Lon(j,1)
            enddo   
            
            do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                Lat1D  (i)= Lat(1,i)
            enddo                   

            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR40'

            !write column indexes
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(1)%VarID, Lon1D)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR50'

            !write line indexes
            STAT_CALL = nf90_put_var(Me%ncid, Me%Dims(2)%VarID, Lat1D)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteLatLon1D - ModuleNETCDF - ERR60'

            deallocate(Lon1D)
            deallocate(Lat1D) 

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteLatLon1D

    !--------------------------------------------------------------------------    


    !--------------------------------------------------------------------------

    subroutine NETCDFGetDimensions (NCDFID, JUB, IUB, KUB, nTime, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID
        integer, optional                           :: JUB
        integer, optional                           :: IUB        
        integer, optional                           :: KUB
        integer, optional                           :: nTime
        integer, optional                           :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: dimid
        integer                                     :: STAT_, ready_
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            
            if (present(JUB)) then

                Me%Dims(1)%ID%Name = trim(Column_Name)

                STAT_CALL = nf90_inq_dimid(Me%ncid, Me%Dims(1)%ID%Name, dimid)            
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR10'
                
                STAT_CALL = nf90_inquire_dimension(Me%ncid, dimid, Me%Dims(1)%ID%Name, JUB)     
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR20'
                
                Me%Dims(1)%LB      = 1
                Me%Dims(1)%UB      = JUB

                Me%Dims(5)%ID%Name = trim(HorizCell_Vertices_Names)
                Me%Dims(5)%LB      = 1
                Me%Dims(5)%UB      = 4
                        
            endif

            if (present(IUB)) then 

                Me%Dims(2)%ID%Name = trim(Line_Name)
                
                STAT_CALL = nf90_inq_dimid(Me%ncid, Me%Dims(2)%ID%Name, dimid)            
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR30'
                
                STAT_CALL = nf90_inquire_dimension(Me%ncid, dimid, Me%Dims(2)%ID%Name, IUB)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR40'
                
                Me%Dims(2)%LB      = 1
                Me%Dims(2)%UB      = IUB

                Me%Dims(6)%ID%Name = trim(Layer_Vertices_Names)
                Me%Dims(6)%LB      = 1
                Me%Dims(6)%UB      = 2
                        
            endif

            if (present(KUB)) then 

                Me%Dims(3)%ID%Name = trim(Layer_Name)

                STAT_CALL = nf90_inq_dimid(Me%ncid, Me%Dims(3)%ID%Name, dimid)            
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR50'
                
                STAT_CALL = nf90_inquire_dimension(Me%ncid, dimid, Me%Dims(3)%ID%Name, KUB)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR60'
                
                Me%Dims(3)%LB      = 1
                Me%Dims(3)%UB      = KUB

                Me%Dims(6)%ID%Name = trim(Layer_Vertices_Names)
                Me%Dims(6)%LB      = 1
                Me%Dims(6)%UB      = 2
                        
            endif
            
            if (present(nTime)) then 

                Me%Dims(4)%ID%Name = trim(Time_Name)

                STAT_CALL = nf90_inq_dimid(Me%ncid, Me%Dims(4)%ID%Name, dimid)            
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR70'
                
                STAT_CALL = nf90_inquire_dimension(Me%ncid, dimid, Me%Dims(4)%ID%Name, nTime)
                if(STAT_CALL /= nf90_noerr) stop 'NETCDFGetDimensions - ModuleNETCDF - ERR80'
                
                Me%Dims(4)%LB      = 1
                Me%Dims(4)%UB      = nTime

            endif
            

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end subroutine NETCDFGetDimensions

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine NETCDFReadGrid2D(NCDFID, Lat, Lon, Lat_Stag, Lon_Stag, SphericX, SphericY, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real(8), dimension(:,:), pointer                :: Lat, Lon
        real(8), dimension(:,:), pointer, optional      :: Lat_Stag, Lon_Stag
        real(8), dimension(:,:), pointer, optional      :: SphericX, SphericY       
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real(8), pointer, dimension(:,:  )              :: Aux2D
        real(8), pointer, dimension(:,:,:)              :: Aux3D        
        integer                                         :: VarID, JUB, IUB, numDims, i, j
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL, nv

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !Get the spatial horizontal dimensions        
            call NETCDFGetDimensions (NCDFID, JUB = JUB, IUB = IUB, STAT = STAT_CALL) 
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR10'

            !Check if the grid is defined with 2D matrixes
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Lon_Name),VarID)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR20'

            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR30'
        
            if (numDims /= 2) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR40' 

            allocate(Aux2D(Me%Dims(1)%LB:Me%Dims(1)%UB,Me%Dims(2)%LB:Me%Dims(2)%UB))
        
            !Read longitude
            STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR50'            
            
            do j=Me%Dims(1)%LB,Me%Dims(1)%UB
            do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                Lon(i,j) = Aux2D(j,i)
            enddo
            enddo

            !Read latitude
            STAT_CALL = nf90_inq_varid(Me%ncid,trim(Lat_Name),VarID)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR60'

            STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR70'
            
            do j=Me%Dims(1)%LB,Me%Dims(1)%UB
            do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                Lat(i,j) = Aux2D(j,i)
            enddo
            enddo
                
            deallocate(Aux2D)
            
            if (present(Lat_Stag) .and. present(Lon_Stag)) then

                allocate(Aux3D(Me%Dims(5)%LB:Me%Dims(5)%UB,                             &
                               Me%Dims(1)%LB:Me%Dims(1)%UB,                             &
                               Me%Dims(2)%LB:Me%Dims(2)%UB))
            
                !Read longitude staggered
                STAT_CALL = nf90_inq_varid(Me%ncid,trim(Lon_Stag_Name),VarID)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR80'
                
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR90'            
                
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                    Lon_Stag(i,j) = Aux3D(1,j,i)
                enddo
                enddo
                enddo

                !Upper limit  
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB                              
                    Lon_Stag(Me%Dims(2)%UB+1,j) = Aux3D(4,j,Me%Dims(2)%UB)
                enddo

                !Right limit  
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB                              
                    Lon_Stag(i,Me%Dims(1)%UB+1) = Aux3D(2,Me%Dims(1)%UB,i)
                enddo
                
                !Upper right corner
                Lon_Stag(Me%Dims(2)%UB+1,Me%Dims(1)%UB+1) = Aux3D(3,Me%Dims(1)%UB,Me%Dims(2)%UB)

                !Read latitude staggered
                STAT_CALL = nf90_inq_varid(Me%ncid,trim(Lat_Stag_Name),VarID)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR100'

                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR110'
                
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                    Lat_Stag(i,j) = Aux3D(1,j,i)
                enddo
                enddo
                enddo

                !Upper limit  
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB                              
                    Lat_Stag(Me%Dims(2)%UB+1,j) = Aux3D(4,j,Me%Dims(2)%UB)
                enddo

                !Right limit  
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB                              
                    Lat_Stag(i,Me%Dims(1)%UB+1) = Aux3D(2,Me%Dims(1)%UB,i)
                enddo
                
                !Upper right corner
                Lat_Stag(Me%Dims(2)%UB+1,Me%Dims(1)%UB+1) = Aux3D(3,Me%Dims(1)%UB,Me%Dims(2)%UB) 
                    
                deallocate(Aux3D)            
            endif

            if (present(SphericX) .and. present(SphericY)) then

                allocate(Aux3D(Me%Dims(5)%LB:Me%Dims(5)%UB,                             &
                               Me%Dims(1)%LB:Me%Dims(1)%UB,                             &
                               Me%Dims(2)%LB:Me%Dims(2)%UB))

            
                !Read google maps X staggered
                STAT_CALL = nf90_inq_varid(Me%ncid,trim(gmaps_x_Name),VarID)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR120'
                
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR130'            
                
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                    SphericX(i,j) = Aux3D(1,j,i)
                enddo
                enddo
                enddo

                !Upper limit  
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB                              
                    SphericX(Me%Dims(2)%UB+1,j) = Aux3D(4,j,Me%Dims(2)%UB)
                enddo

                !Right limit  
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB                              
                    SphericX(i,Me%Dims(1)%UB+1) = Aux3D(2,Me%Dims(1)%UB,i)
                enddo
                
                !Upper right corner
                SphericX(Me%Dims(2)%UB+1,Me%Dims(1)%UB+1) = Aux3D(3,Me%Dims(1)%UB,Me%Dims(2)%UB) 
                    
                
                !Read google maps Y staggered
                STAT_CALL = nf90_inq_varid(Me%ncid,trim(gmaps_y_Name),VarID)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR140'
                
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadGrid2D - ModuleNETCDF - ERR150'            
                
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB
                do nv= Me%Dims(5)%LB,Me%Dims(5)%UB
                    SphericY(i,j) = Aux3D(1,j,i)
                enddo
                enddo
                enddo

                !Upper limit  
                do j=Me%Dims(1)%LB,Me%Dims(1)%UB                              
                    SphericY(Me%Dims(2)%UB+1,j) = Aux3D(4,j,Me%Dims(2)%UB)
                enddo

                !Right limit  
                do i=Me%Dims(2)%LB,Me%Dims(2)%UB                              
                    SphericY(i,Me%Dims(1)%UB+1) = Aux3D(2,Me%Dims(1)%UB,i)
                enddo
                
                !Upper right corner
                SphericY(Me%Dims(2)%UB+1,Me%Dims(1)%UB+1) = Aux3D(3,Me%Dims(1)%UB,Me%Dims(2)%UB) 
                    
                deallocate(Aux3D)          
            endif                            
            
            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadGrid2D

    !--------------------------------------------------------------------------    

    subroutine NETCDFWriteVert(NCDFID, Vert, VertCoordinate, SimpleGrid, OffSet, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:  ), pointer                   :: Vert
        logical                                         :: VertCoordinate
        logical, intent(in), optional                   :: SimpleGrid
        real,    intent(in), optional                   :: OffSet
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real                                            :: FillValue, MissingValue
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        logical                                         :: SimpleGrid_

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVert - ModuleNETCDF - ERR00'
            
            !define depth as variable
            STAT_CALL = nf90_def_var(Me%ncid, depth_Name, NF90_FLOAT, Me%Dims(3)%ID%Number, Me%Dims(3)%VarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVert - ModuleNETCDF - ERR01'
            
            MissingValue    = FillValueReal
            FillValue       = FillValueReal
            

            if (present(SimpleGrid)) then
                SimpleGrid_ = SimpleGrid
            else
                SimpleGrid_ = .false. 
            endif                

                

            !Is it a sigma grid?
            if(VertCoordinate)then

                call NETCDFWriteAttributes(Me%Dims(3)%VarID, LongName     = depth_Name,          &
                                                                StandardName = depth_Name,          &
                                                                Units        = null_str,            &
                                                                Positive     = "down",              &
!                                                                FillValue    = MissingValue,        &
                                                                ValidMin     = 0.,                  &
                                                                ValidMax     = 1.,                  & 
                                                                Add_Offset   = OffSet,              &
!                                                                MissingValue = MissingValue,        &
                                                                bounds       = trim(depth_Stag_Name))
               
                !No? Then it must be a z-coordinate grid
            else
                
                !if (SimpleGrid_) then

                    call NETCDFWriteAttributes(Me%Dims(3)%VarID, LongName               = "Depth",      &
                                                                 StandardName           = depth_Name,   &
                                                                 Units                  = "m",          &
                                                                 Positive               = "down",       &
                                                                 axis                   = "Z",          &
                                                                 CoordinateAxisType     = "Height",     &
                                                                 CoordinateZisPositive  = "down",       &    
                                                                 bounds                 = trim(depth_Stag_Name))                    

!                else                    
!                
!                    call NETCDFWriteAttributes(Me%Dims(3)%VarID, LongName     = depth_Name,          &
!                                                                 StandardName = depth_Name,          &
!                                                                 Units        = "m",                 &
!                                                                 Positive     = "down",              &
!!                                                                 FillValue    = FillValue,           &
!                                                                 ValidMin     = -50.,                &
!                                                                 ValidMax     = 10000.,              &
!                                                                 Add_Offset   = OffSet,              &
!!                                                                 MissingValue = MissingValue,        &
!                                                                 bounds       = trim(depth_Stag_Name))
!                endif                                                                                                
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
    !--------------------------------------------------------------------------    

    subroutine NETCDFWriteVertStag(NCDFID, VertStag, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real, dimension(:  ), pointer                   :: VertStag
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real,    dimension(:,:), pointer                :: Aux2D
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: VertStagVarID, nv, k
        integer, dimension(2)                           :: VertStagDimID

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            !enter definition mode
            STAT_CALL = nf90_redef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVertStag - ModuleNETCDF - ERR10'

            VertStagDimID(1) = Me%Dims(6)%ID%Number ! nvv                     
            VertStagDimID(2) = Me%Dims(3)%ID%Number !z
            
            !define depth as variable
            STAT_CALL = nf90_def_var(Me%ncid, depth_Stag_Name, NF90_FLOAT, VertStagDimID, VertStagVarID)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVertStag - ModuleNETCDF - ERR20'


            !exit definition mode
            STAT_CALL = nf90_enddef(ncid = Me%ncid)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVertStag - ModuleNETCDF - ERR30'
            
            allocate(Aux2D(Me%Dims(6)%LB:Me%Dims(6)%UB,Me%Dims(3)%LB:Me%Dims(3)%UB))
            
            do k = Me%Dims(3)%LB,Me%Dims(3)%UB
                do nv = Me%Dims(6)%LB,Me%Dims(6)%UB
                    Aux2D(nv,k)= VertStag (k+nv-1)
                enddo
            enddo
            

            !write Vertical coordinate
            STAT_CALL = nf90_put_var(Me%ncid, VertStagVarID, Aux2D)
            if(STAT_CALL /= nf90_noerr) stop 'NETCDFWriteVertStag - ModuleNETCDF - ERR40'
            
            deallocate(Aux2D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFWriteVertStag

    !--------------------------------------------------------------------------
    logical function NETCDFWithVert(NCDFID, STAT)
    
        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL, VarID

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then                
                
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Depth_Name),VarID)
            if (STAT_CALL == nf90_noerr) then
                NETCDFWithVert = .true.
            else                  
                NETCDFWithVert = .false.
            endif

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_

    end function NETCDFWithVert
                
    !--------------------------------------------------------------------------    

    subroutine NETCDFReadVert(NCDFID, CenterCellDepth, VerticalZ, nInstant,             &
                              ILB, IUB, JLB, JUB, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real,    dimension(:,:,:), pointer, optional    :: CenterCellDepth        
        real,    dimension(:,:,:), pointer, optional    :: VerticalZ
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB, KLB, KUB         
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D
        real, dimension(:,:  ), pointer                 :: Aux2D        
        real, dimension(:    ), pointer                 :: Aux1D                
        integer                                         :: VarID, numDims, i, j, k
        integer                                         :: ILB_, IUB_, JLB_, JUB_, KLB_, KUB_
        integer                                         :: in, jn, kn
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            if (present(KUB)) then 
                KUB_ = KUB
            else
                !Get the spatial vertical dimension        
                call NETCDFGetDimensions (NCDFID, KUB = KUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR10'
            endif
            
            if (present(KLB)) then 
                KLB_ = KLB
            else
                KLB_ = 1
            endif
            
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
            
            in = IUB_ - ILB_ + 1
            jn = JUB_ - JLB_ + 1                    
            kn = KUB_ - KLB_ + 1              

            if (present(CenterCellDepth)) then


                !Check if depth is defined in :
                ! a) 1D (variable in depth)
                ! b) 3D (variable in depth and horizontaly)
                ! c) 4D (variable in depth, horizontaly and in time)
                
                STAT_CALL=nf90_inq_varid(Me%ncid,trim(Depth_Name),VarID)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR20'

                STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR30'
            
                if (numDims /= 1 .and. numDims /= 3 .and. numDims /= 4) stop 'NETCDFReadVert - ModuleNETCDF - ERR40' 
                
                if (numDims == 4 .and. .not. present(nInstant)) stop 'NETCDFReadVert - ModuleNETCDF - ERR50' 
                

                if (numDims == 1) then
                
                    allocate(Aux1D(1:kn))
            
                    !Read 1D depth
                    STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux1D)
                    if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR60'            
                
                    do j=JLB_, JUB_
                    do i=ILB_, IUB_
                        CenterCellDepth(i,j,KLB_:KUB_) = Aux1D(KLB_:KUB_)
                    enddo
                    enddo
                    
                    deallocate(Aux1D)

                elseif (numDims == 3 .or. numDims == 4) then
                
                    allocate(Aux3D (1:jn,1:in,1:kn))
                                    
                    if (numDims == 3) then                                    
            
                        !Read 3D depth
                        STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                   &
                            start = (/ JLB_, ILB_, KLB_/),                              &
                            count = (/   jn,   in,   kn/))                             
                        if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR70'            
                    else
                        !Read 4D depth
                        STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                   &
                            start = (/ JLB_, ILB_, KLB_,  ninstant/),                   &
                            count = (/   jn,   in,   kn,  1       /))                             
                        if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR80'            
                    
                    endif

                    do j=JLB_,JUB_
                    do i=ILB_,IUB_
                    do k=KLB_,KUB_
                        CenterCellDepth(i,j,k) = Aux3D(j-JLB_+1,i-ILB_+1,k-KLB_+1)
                    enddo
                    enddo
                    enddo

                    deallocate(Aux3D)
                else
                    stop 'NETCDFReadVert - ModuleNETCDF - ERR90'            
                endif
            endif

            if (present(VerticalZ)) then

                !Check if depth of the cel interfaces are defined in :
                ! a) 1D (variable in depth)
                ! b) 3D (variable in depth and horizontaly)
                ! c) 4D (variable in depth, horizontaly and in time)
                
                STAT_CALL=nf90_inq_varid(Me%ncid,trim(depth_Stag_Name),VarID)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR100'

                STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR110'
            
            
                allocate(Aux2D(2, 1:kn))
        
                !Read 1D depth
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D)
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadVert - ModuleNETCDF - ERR140'            
            
                do j=JLB_,JUB_
                do i=ILB_,IUB_
                do k=KLB_,KUB_
                    VerticalZ(i,j,k)                = Aux2D(2,k)
                    if (k==KLB_) VerticalZ(i,j,k-1) = Aux2D(1,k)
                enddo
                enddo
                enddo
                
                deallocate(Aux2D)
            else
                stop 'NETCDFReadVert - ModuleNETCDF - ERR150'
            endif


            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadVert

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------    

    subroutine NETCDFReadDataI4_3D(NCDFID, Array3D, Name, nInstant,     &
                                   ILB, IUB, JLB, JUB, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        integer(4), dimension(:,:,:), pointer           :: Array3D        
        character(len = *)                              :: Name
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB, KLB, KUB 
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D        
        integer                                         :: VarID, numDims, i, j, k
        integer                                         :: ILB_, IUB_, JLB_, JUB_, KLB_, KUB_         
        integer                                         :: in, jn, kn        
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(KUB)) then 
                KUB_ = KUB
            else
                !Get the spatial vertical dimension        
                call NETCDFGetDimensions (NCDFID, KUB = KUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR10'
            endif
            
            if (present(KLB)) then 
                KLB_ = KLB
            else
                KLB_ = 1
            endif
            
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
                        
            !Check if depth is defined in :
            ! b) 3D (variable in depth and horizontaly)
            ! c) 4D (variable in depth, horizontaly and in time)
            
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Name),VarID)
            if (STAT_CALL /= nf90_noerr) then
                write(*,*) "Property ", trim(Name)," not found in NetCDF file ",trim(Me%FileName)
                stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR40'
            endif

            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR50'
        
            if (numDims /= 3 .and. numDims /= 4) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR60' 
            
            if (numDims == 4 .and. .not. present(nInstant)) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR70' 
            
           
            jn = 1-JLB_+JUB_
            in = 1-ILB_+IUB_
            kn = 1-KLB_+KUB_
            
            allocate(Aux3D (1:jn,1:in,1:kn))
                            
            if (numDims == 3) then                                    
    
                !Read 3D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                           &
                    start = (/ JLB_, ILB_, KLB_/),                                      &
                    count = (/   jn,   in,   kn/))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR80'                                
            else
                !Read 4D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                           &
                    start = (/ JLB_, ILB_, KLB_,  ninstant /),                          &
                    count = (/   jn,   in,   kn, 1       /))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_3D - ModuleNETCDF - ERR90'            
            
            endif
            
            do j=JLB_,JUB_
            do i=ILB_,IUB_
            do k=KLB_,KUB_
                Array3D(i,j,k) = Aux3D(j-JLB_+1,i-ILB_+1,k-KLB_+1)
            enddo
            enddo
            enddo 

            deallocate(Aux3D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadDataI4_3D
    !--------------------------------------------------------------------------    


    !--------------------------------------------------------------------------    

    subroutine NETCDFReadDataR4_3D(NCDFID, Array3D, Name, nInstant, &
                                   ILB, IUB, JLB, JUB, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real(4), dimension(:,:,:), pointer              :: Array3D        
        character(len = *)                              :: Name
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB, KLB, KUB 
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D        
        integer                                         :: VarID, numDims, i, j, k
        integer                                         :: ILB_, IUB_, JLB_, JUB_, KLB_, KUB_         
        integer                                         :: in, jn, kn
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(KUB)) then 
                KUB_ = KUB
            else
                !Get the spatial vertical dimension        
                call NETCDFGetDimensions (NCDFID, KUB = KUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR10'
            endif
            
            if (present(KLB)) then 
                KLB_ = KLB
            else
                KLB_ = 1
            endif
            
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
                        
            !Check if depth is defined in :
            ! b) 3D (variable in depth and horizontaly)
            ! c) 4D (variable in depth, horizontaly and in time)
            
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Name),VarID)
            if (STAT_CALL /= nf90_noerr) then
                write(*,*) "Property ", trim(Name)," not found in NetCDF file ",trim(Me%FileName)
                stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR40'
            endif

            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR50'
        
            if (numDims /= 3 .and. numDims /= 4) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR60' 
            
            if (numDims == 4 .and. .not. present(nInstant)) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR70' 
            
            
            jn = 1-JLB_+JUB_
            in = 1-ILB_+IUB_
            kn = 1-KLB_+KUB_
            
            allocate(Aux3D (1:jn,1:in,1:kn))
                            
            if (numDims == 3) then                                    
    
                !Read 3D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                           &
                    start = (/ JLB_, ILB_, KLB_/),                                      &
                    count = (/   jn,   in,   kn/))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR80'                                
            else
                !Read 4D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                           &
                    start = (/ JLB_, ILB_, KLB_,  ninstant /),                          &
                    count = (/   jn,   in,   kn, 1       /))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_3D - ModuleNETCDF - ERR90'            
            
            endif
            
            do j=JLB_,JUB_
            do i=ILB_,IUB_
            do k=KLB_,KUB_
                Array3D(i,j,k) = Aux3D(j-JLB_+1,i-ILB_+1,k-KLB_+1)
            enddo
            enddo
            enddo 

            deallocate(Aux3D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadDataR4_3D
    !--------------------------------------------------------------------------    
!--------------------------------------------------------------------------    

    subroutine NETCDFReadDataR8_3D(NCDFID, Array3D, Name, nInstant,      &
                                   ILB, IUB, JLB, JUB, KLB, KUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real(8), dimension(:,:,:), pointer              :: Array3D        
        character(len = *)                              :: Name
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB, KLB, KUB 
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D        
        integer                                         :: VarID, numDims, i, j, k
        integer                                         :: ILB_, IUB_, JLB_, JUB_, KLB_, KUB_         
        integer                                         :: in, jn, kn        
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(KUB)) then 
                KUB_ = KUB
            else
                !Get the spatial vertical dimension        
                call NETCDFGetDimensions (NCDFID, KUB = KUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR10'
            endif
            
            if (present(KLB)) then 
                KLB_ = KLB
            else
                KLB_ = 1
            endif
            
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
                        
            !Check if depth is defined in :
            ! b) 3D (variable in depth and horizontaly)
            ! c) 4D (variable in depth, horizontaly and in time)
            
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Name),VarID)
            if (STAT_CALL /= nf90_noerr) then
                write(*,*) "Property ", trim(Name)," not found in NetCDF file ",trim(Me%FileName)
                stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR40'
            endif

            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR50'
        
            if (numDims /= 3 .and. numDims /= 4) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR60' 
            
            if (numDims == 4 .and. .not. present(nInstant)) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR70' 
            
           
            jn = 1-JLB_+JUB_
            in = 1-ILB_+IUB_
            kn = 1-KLB_+KUB_
            
            allocate(Aux3D (1:jn,1:in,1:kn))
                            
            if (numDims == 3) then                                    
    
                !Read 3D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                           &
                    start = (/ JLB_, ILB_, KLB_/),                                      &
                    count = (/   jn,   in,   kn/))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR80'                                
            else
                !Read 4D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux3D,                           &
                    start = (/ JLB_, ILB_, KLB_,  ninstant /),                          &
                    count = (/   jn,   in,   kn, 1       /))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_3D - ModuleNETCDF - ERR90'            
            
            endif
            
            do j=JLB_,JUB_
            do i=ILB_,IUB_
            do k=KLB_,KUB_
                Array3D(i,j,k) = Aux3D(j-JLB_+1,i-ILB_+1,k-KLB_+1)
            enddo
            enddo
            enddo 

            deallocate(Aux3D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadDataR8_3D
    !--------------------------------------------------------------------------    


    
    !--------------------------------------------------------------------------    

    subroutine NETCDFReadDataI4_2D(NCDFID, Array2D, Name, nInstant,                     &
                                  ILB, IUB, JLB, JUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        integer(4), dimension(:,:), pointer             :: Array2D        
        character(len = *)                              :: Name
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D        
        integer                                         :: VarID, numDims, i, j
        integer                                         :: ILB_, IUB_, JLB_, JUB_
        integer                                         :: in, jn
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL
        integer                                         :: nInstant_        

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
                        
            ! 1) 2D (variable in horizontaly)
            ! 2) 3D (variable in horizontaly and in time)
            ! 3) 3D (variable in horizontaly and one layer)
            
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Name),VarID)
            if (STAT_CALL /= nf90_noerr) then
                write(*,*) "Property ", trim(Name)," not found in NetCDF file ",trim(Me%FileName)
                stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR40'
            endif

            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR50'
        
            if (numDims /= 2 .and. numDims /= 3) stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR60' 
            
            if (numDims == 3) then
                if (present(nInstant)) then
                    nInstant_ = nInstant
                else                    
                    nInstant_ = 1
                    !stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR70' 
                endif    
            endif
            
            jn = 1-JLB_+JUB_
            in = 1-ILB_+IUB_
            
            allocate(Aux2D (1:jn,1:in))
                            
            if (numDims == 2) then                                    
    
                !Read 2D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D,                           &
                    start = (/ JLB_, ILB_/),                                            &
                    count = (/   jn,   in/))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR80'                                
            else
                !Read 3D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D,                           &
                    start = (/ JLB_, ILB_,  ninstant_ /),                                &
                    count = (/   jn,   in, 1       /))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataI4_2D - ModuleNETCDF - ERR90'            
            
            endif
            
            do j=JLB_,JUB_
            do i=ILB_,IUB_
                Array2D(i,j) = Aux2D(j-JLB_+1,i-ILB_+1)
            enddo
            enddo 

            deallocate(Aux2D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadDataI4_2D
    
    !--------------------------------------------------------------------------    

    subroutine NETCDFReadDataR4_2D(NCDFID, Array2D, Name, nInstant,                     &
                                  ILB, IUB, JLB, JUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real(4), dimension(:,:), pointer                :: Array2D        
        character(len = *)                              :: Name
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D        
        integer                                         :: VarID, numDims, i, j
        integer                                         :: ILB_, IUB_, JLB_, JUB_
        integer                                         :: in, jn
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
                        
            ! 1) 2D (variable in horizontaly)
            ! 2) 3D (variable in horizontaly and in time)
            
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Name),VarID)
            if (STAT_CALL /= nf90_noerr) then
                write(*,*) "Property ", trim(Name)," not found in NetCDF file ",trim(Me%FileName)
                stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR40'
            endif

            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR50'
        
            if (numDims /= 2 .and. numDims /= 3) stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR60' 
            
            if (numDims == 3 .and. .not. present(nInstant)) stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR70' 
            
            
            jn = 1-JLB_+JUB_
            in = 1-ILB_+IUB_
            
            allocate(Aux2D (1:jn,1:in))
                            
            if (numDims == 2) then                                    
    
                !Read 2D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D,                           &
                    start = (/ JLB_, ILB_/),                                            &
                    count = (/   jn,   in/))                             
                if (STAT_CALL /= nf90_noerr) then
                    stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR80'
                endif
            else
                !Read 3D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D,                           &
                    start = (/ JLB_, ILB_,  ninstant /),                                &
                    count = (/   jn,   in, 1       /))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR4_2D - ModuleNETCDF - ERR90'            
            
            endif
            
            do j=JLB_,JUB_
            do i=ILB_,IUB_
                Array2D(i,j) = Aux2D(j-JLB_+1,i-ILB_+1)
            enddo
            enddo 

            deallocate(Aux2D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadDataR4_2D
    !--------------------------------------------------------------------------    

    subroutine NETCDFReadDataR8_2D(NCDFID, Array2D, Name, nInstant,                     &
                                  ILB, IUB, JLB, JUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: NCDFID
        real(8), dimension(:,:), pointer                :: Array2D        
        character(len = *)                              :: Name
        integer, optional                               :: nInstant
        integer, optional                               :: ILB, IUB, JLB, JUB
        integer, optional                               :: STAT
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D        
        integer                                         :: VarID, numDims, i, j
        integer                                         :: ILB_, IUB_, JLB_, JUB_
        integer                                         :: in, jn
        integer                                         :: STAT_, ready_
        integer                                         :: STAT_CALL

        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready (NCDFID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            if (present(JUB)) then 
                JUB_ = JUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, JUB = JUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR20'
            endif
            
            if (present(JLB)) then 
                JLB_ = JLB
            else
                JLB_ = 1
            endif
                    
            if (present(IUB)) then 
                IUB_ = IUB
            else
                !Get one of the horizontal spatial dimension        
                call NETCDFGetDimensions (NCDFID, IUB = IUB_, STAT = STAT_CALL) 
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR30'
            endif
            
            if (present(ILB)) then 
                ILB_ = ILB
            else
                ILB_ = 1
            endif                    
                        
            ! 1) 2D (variable in horizontaly)
            ! 2) 3D (variable in horizontaly and in time)
            
            STAT_CALL=nf90_inq_varid(Me%ncid,trim(Name),VarID)
            if (STAT_CALL /= nf90_noerr) then
                write(*,*) "Property ", trim(Name)," not found in NetCDF file ",trim(Me%FileName)
                stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR40'
            endif
            
            STAT_CALL = nf90_inquire_variable(Me%ncid, VarID, ndims = numDims)
            if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR50'
        
            if (numDims /= 2 .and. numDims /= 3) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR60' 
            
            if (numDims == 3 .and. .not. present(nInstant)) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR70' 
            
            
            jn = 1-JLB_+JUB_
            in = 1-ILB_+IUB_
            
            allocate(Aux2D (1:jn,1:in))
                            
            if (numDims == 2) then                                    
    
                !Read 2D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D,                           &
                    start = (/ JLB_, ILB_/),                                            &
                    count = (/   jn,   in/))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR80'                                
            else
                !Read 3D Field
                STAT_CALL = NF90_GET_VAR(Me%ncid,VarID,Aux2D,                           &
                    start = (/ JLB_, ILB_,  ninstant /),                                &
                    count = (/   jn,   in, 1       /))                             
                if (STAT_CALL /= nf90_noerr) stop 'NETCDFReadDataR8_2D - ModuleNETCDF - ERR90'            
            
            endif
            
            do j=JLB_,JUB_
            do i=ILB_,IUB_
                Array2D(i,j) = Aux2D(j-JLB_+1,i-ILB_+1)
            enddo
            enddo 

            deallocate(Aux2D)

            STAT_ = SUCCESS_

        else

            STAT_ = ready_

        endif

        if (present(STAT)) STAT = STAT_


    end subroutine NETCDFReadDataR8_2D
    
    
    !--------------------------------------------------------------------------    


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillNETCDF(NCDFID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: NCDFID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers
        integer                             :: STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(NCDFID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mNETCDF_,  Me%InstanceID)

            if (nUsers == 0) then

                STAT_CALL =nf90_close(Me%ncid)
                if(STAT_CALL /= nf90_noerr) stop 'ModuleNETCDF - KillNETCDF - ERR01' 

                !Deallocates Instance
                call DeallocateInstance ()

                NCDFID = 0
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

    subroutine LocateObjNETCDF (NCDFID)

        !Arguments-------------------------------------------------------------
        integer                                     :: NCDFID

        !Local-----------------------------------------------------------------

        Me => FirstObjNETCDF
        do while (associated (Me))
            if (Me%InstanceID == NCDFID) exit
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









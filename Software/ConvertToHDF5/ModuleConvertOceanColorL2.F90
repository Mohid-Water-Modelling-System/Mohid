!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ConvertOceanColorL2
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Pedro Pina - v4.0
! DESCRIPTION   : Module to convert OceanColorL2 (Modis and SeaWifs) Format files into HDF5 format
! DOWNLOAD      : To download OCL2 files go to http://oceancolor.gsfc.nasa.gov/
!------------------------------------------------------------------------------


Module ModuleConvertOceanColorL2

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleTime
    use ModuleGridData
    use ModuleDrawing


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

  
    public  ::  StartConvertOceanColorL2
    private ::  Construct_FieldsList
    private ::        Construct_Field
    private ::        Add_Field
    !Types---------------------------------------------------------------------
    
    private :: T_Attribute
    type       T_Attribute                       
        character(len=256)                  :: Name
        integer                             :: Att_Type
        integer                             :: Length
        character(len=256)                  :: String = ' '
        real                                :: ValueReal
        integer                             :: ValueInt
    end type T_Attribute
    
    private :: T_OceanColorL2
    type       T_OceanColorL2
        
        character(len=PathLength)                :: L2FileName
        character(len=StringLength)              :: ParameterName
        character(len=StringLength)              :: Units  
        integer                                  :: IDNumber
        integer, dimension(:,:,:),  pointer      :: OpenPoints3D
        integer,dimension(:),pointer             :: ProdFound 
        type(T_Time)                             :: Date
        real, dimension(:,:,:),    pointer       :: Scalar
        type(T_OceanColorL2),      pointer       :: Next
        type(T_OceanColorL2),      pointer       :: Prev
        type(T_Attribute), dimension(:), pointer :: Attributes
        character(len=PathLength)                :: GridFileName
        type (T_Size2D)                          :: WorkSize
        type (T_Size2D)                          :: Size
        real, dimension(:, :), pointer           :: XX_IE, YY_IE
        real, dimension(:, :), pointer           :: LongitudeConn, LatitudeConn
        integer, dimension(:, :), pointer        :: DefineCellsMap
    end type  T_OceanColorL2
    
    

   


    private :: T_OceanColorL2Format
    type       T_OceanColorL2Format
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjModisL3MappedGrid = 0
        integer                                 :: ObjNewGrid           = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjBathymetry        = 0
! Header parameters
        real                                    :: LatStep  
        real                                    :: LongStep 
        real                                    :: OrigLat  
        real                                    :: OrigLong 

        
     
        character(len=PathLength)               :: OutputFileName
        logical                                 :: FirstFile            = .true. 
        integer                                 :: NoExtraction         
        integer                                 :: Unit
        real                                    :: MaxLong
        real                                    :: MinLong
        real                                    :: MaxLat
        real                                    :: MinLat
        integer                                 :: FieldsNumber
        real, dimension(:,:),  pointer          :: Bathymetry
        real , dimension(:),  pointer           :: XX
        real , dimension(:),  pointer           :: YY
        type(T_OceanColorL2),      pointer           :: FirstField
        type(T_OceanColorL2),      pointer           :: LastField
        

        
        integer                                 :: CHLOR_A
        integer                                 :: SST
        integer                                 :: K_490
        integer                                 :: TAU_865
        integer                                 :: TAU_869
        integer                                 :: EPS_78
        integer                                 :: ANGSTROM_510
        integer                                 :: ANGSTROM_531
        integer                                 :: NLW_412
        integer                                 :: NLW_443
        integer                                 :: NLW_490
        integer                                 :: NLW_510
        integer                                 :: NLW_555
        integer                                 :: NLW_670
        integer                                 :: NLW_488
        integer                                 :: NLW_531
        integer                                 :: NLW_551
        integer                                 :: NLW_667
		integer                                 :: CALCITE
		integer                                 :: TSM_CLARK
		integer                                 :: POC_CLARK
		integer                                 :: CHL_CLARK 
        
        integer                                 :: CHLOR_A_      = 1
        integer                                 :: SST_          = 2
        integer                                 :: K_490_        = 3
        integer                                 :: TAU_865_      = 4
        integer                                 :: TAU_869_      = 5
        integer                                 :: EPS_78_       = 6
        integer                                 :: ANGSTROM_510_ = 7
        integer                                 :: ANGSTROM_531_ = 8
        integer                                 :: NLW_412_      = 9
        integer                                 :: NLW_443_      = 10
        integer                                 :: NLW_490_      = 11
        integer                                 :: NLW_510_      = 12
        integer                                 :: NLW_555_      = 13
        integer                                 :: NLW_670_      = 14
        integer                                 :: NLW_488_      = 15
        integer                                 :: NLW_531_      = 16
        integer                                 :: NLW_551_      = 17
        integer                                 :: NLW_667_      = 18
        integer                                 :: CALCITE_      = 19
		integer                                 :: TSM_CLARK_    = 20
		integer                                 :: POC_CLARK_    = 21
		integer                                 :: CHL_CLARK_    = 22  
        

        
        character(len=40)                       :: Char_chlor_a       ='chlor_a'
        character(len=40)                       :: Char_sst           ='sst'
        character(len=40)                       :: Char_K_490         ='K_490'
        character(len=40)                       :: Char_tau_865       ='tau_865'
        character(len=40)                       :: Char_tau_869       ='tau_869'
        character(len=40)                       :: Char_angstrom_510  ='angstrom_510'
        character(len=40)                       :: Char_angstrom_531  ='angstrom_531'
        character(len=40)                       :: Char_eps_78        ='eps_78'
        character(len=40)                       :: Char_nLw_412       ='nLw_412'
        character(len=40)                       :: Char_nLw_443       ='nLw_443'
        character(len=40)                       :: Char_nLw_490       ='nLw_490'
        character(len=40)                       :: Char_nLw_510       ='nLw_510'
        character(len=40)                       :: Char_nLw_555       ='nLw_555'
        character(len=40)                       :: Char_nLw_670       ='nLw_670'
        character(len=40)                       :: Char_nLw_488       ='nLw_488'
        character(len=40)                       :: Char_nLw_531       ='nLw_531'
        character(len=40)                       :: Char_nLw_551       ='nLw_551'
        character(len=40)                       :: Char_nLw_667       ='nLw_667'
		character(len=40)                       :: Char_calcite       ='calcite'
		character(len=40)                       :: Char_tsm_clark     ='tsm_clark'
		character(len=40)                       :: Char_poc_clark     ='poc_clark '
		character(len=40)                       :: Char_chl_clark     ='chl_clark'
        character(len=40)                       :: Char_l2_flags      ='l2_flags'

        character(len=40)                       :: Units_chla         ='mg/m^3'
        character(len=40)                       :: Units_sst          ='ºC' 
        character(len=40)                       :: Units_nLw          ='mW cm^-2 um^-1 sr^-1' 
        character(len=40)                       :: Units_K            ='m-1'
        character(len=40)                       :: Units_Calcite      ='moles/m^3'
        character(len=40)                       :: Units_TSM          ='mg/l'
		character(len=40)                       :: Units_POC          ='mg/l'
        character(len=40)                       :: Units_NoDim        ='dimensionless'
        
        integer                                 :: nProd
        character(len=40),dimension(:),pointer  :: ProdNames
        integer,dimension(:),pointer            :: ProdID
        character(len=40),dimension(:),pointer  :: ProdUnits
             
       
    end type  T_OceanColorL2Format


    type(T_OceanColorL2Format),         pointer     :: Me

    !--------------------------------------------------------------------------

      INTERFACE TO SUBROUTINE ReadDataL2 [C,ALIAS:'_ReadDataL2'] &
      (Mini,Maxi, ProdId)
      INTEGER Mini [VALUE]
      INTEGER Maxi [VALUE]
      INTEGER ProdId [VALUE]
      END

      INTERFACE TO SUBROUTINE ReadFirstL2 [C,ALIAS:'_ReadFirstL2'] ()
      
      END


      INTERFACE TO SUBROUTINE closeFile [C,ALIAS:'_closeFile'] ()
      
      END 
 
      INTERFACE TO SUBROUTINE GetMetaData [C,ALIAS:'_GetMetaData'] &
                              (ullat, ullon, urlat, urlon, &
                               lllat, lllon, lrlat, lrlon,  &
                               northlat, southlat, westlon,eastlon,&
                               columns,lines,nProd,geointerp)

        REAL*4  ullat     [REFERENCE]
        REAL*4  ullon     [REFERENCE]
        REAL*4  lllat     [REFERENCE]
        REAL*4  southlat  [REFERENCE]
        REAL*4  urlat     [REFERENCE]
        REAL*4  urlon     [REFERENCE]
        REAL*4  lllon     [REFERENCE]
        REAL*4  lrlat     [REFERENCE]
        REAL*4  lrlon     [REFERENCE]
        REAL*4  northlat  [REFERENCE]
        REAL*4  westlon   [REFERENCE]
        REAL*4  eastlon   [REFERENCE]
        INTEGER lines     [REFERENCE]
        INTEGER columns   [REFERENCE]
        INTEGER nProd     [REFERENCE]
        INTEGER geointerp [REFERENCE]

      END 

      INTERFACE TO SUBROUTINE GetLatLon [C,ALIAS:'_GetLatLon'] &
                              (i,j,lat,lon)
        REAL*4    lat     [REFERENCE]
        REAL*4    lon     [REFERENCE]
        INTEGER i [VALUE]
        INTEGER j [VALUE]
      END

      INTERFACE TO SUBROUTINE GetDataL2 [C,ALIAS:'_GetData'] & 
         (i,j,idProd, sdata,flag)

        REAL*4  sdata [REFERENCE]
        INTEGER flag  [REFERENCE]
        INTEGER i [VALUE]
        INTEGER j [VALUE]
        INTEGER idProd [VALUE]

      END
    
      
     INTERFACE TO SUBROUTINE OpenL2Seadas [ALIAS:'_OpenL2Seadas@8'] & 
         (filename)
       CHARACTER*(*) filename
     END

     INTERFACE TO SUBROUTINE GetProdName [ALIAS:'_GetProdName@1'] &
        (ProductName) 
       CHARACTER*(*) ProductName
     END 
    !--------------------------------------------------------------------------
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartConvertOceanColorL2(EnterDataID, ClientNumber,  STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !External----------------------------------------------------------------
        
        
        !Local-------------------------------------------------------------------
        integer                                         :: nUsers, i
        integer                                         :: STAT_CALL
    
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)
             

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call Construct_FieldsList(ClientNumber)
        call OutputFields

        Call KillOceanColorL2

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'StartConvertOceanColorL2 - ModuleConvertOceanColorL2 - ERR01' 

        STAT = SUCCESS_

        deallocate(Me)
        nullify   (Me)



    end subroutine StartConvertOceanColorL2
  ! ------------------------------------------------------------------------------
   
    subroutine Construct_FieldsList(ClientNumber)

        ! Arguments -----------------------------------------------------------

        integer,           intent(IN )              :: ClientNumber

        !External----------------------------------------------------------------
        integer                         :: STAT_CALL
        logical                         :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_OceanColorL2), pointer      :: NewField
        integer                              :: iflag
        !------------------------------------------------------------------------

       
        me%nProd = 0
        allocate (me%ProdNames(20))
        allocate (me%ProdUnits(20))

        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_FILE',                      &
                     ClientModule = 'ConvertOceanColorL2',              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of Output file'
            stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR02'
        end if
       
        call GetData(Me%NoExtraction,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'NO_EXTRACTION',                    &
                     ClientModule = 'ConvertOceanColorL2',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR03'
        
        if (iflag == 0)then
            write(*,*)'Must specify NoExtraction'
            stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR09'
        end if
       if (me%NoExtraction.eq.0) then 
        call GetData(Me%MaxLong,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MAX_LONG',                         &
                     ClientModule = 'ConvertOceanColorL2',              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR02'
        
        if (iflag == 0)then
            write(*,*)'Must specify MAX_LONG'
            stop 'Construct_FieldsList - ModuleEUCenterFormat - ERR06'
        end if


        call GetData(Me%MinLong,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MIN_LONG',                         &
                     ClientModule = 'ConvertOceanColorL2',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR03'
        
        if (iflag == 0)then
            write(*,*)'Must specify MIN_LONG'
            stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR07'
        end if
        
        
        call GetData(Me%MaxLat,                                         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MAX_LAT',                          &
                     ClientModule = 'ConvertOceanColorL2',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR04'
        
        if (iflag == 0)then
            write(*,*)'Must specify MAX_LAT'
            stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR08'
        end if


        call GetData(Me%MinLat,                                         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MIN_LAT',                          &
                     ClientModule = 'ConvertOceanColorL2',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR03'
        
        if (iflag == 0)then
            write(*,*)'Must specify MIN_LAT'
            stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR09'
        end if

     endif


          call GetData(Me%CHLOR_A, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='CHLOR_A', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR02 .'

       if (Me%CHLOR_A.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_CHLOR_A
       endif
 
        call GetData(Me%SST, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='SST', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR03 .'

       if (Me%SST.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd)  = me%Char_SST
         me%ProdUnits(me%NProd)  = me%Units_sst
       endif
 
        call GetData(Me%K_490, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='K_490', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR04 .'

       if (Me%K_490.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_K_490
         me%ProdUnits(me%NProd)  = me%Units_K
       endif
 
        call GetData(Me%TAU_865, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='TAU_865', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR05 .'

       if (Me%TAU_865.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_TAU_865
         me%ProdUnits(me%NProd)  = me%Units_NoDim
       endif
 
        call GetData(Me%TAU_869, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='TAU_869', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR06 .'

       if (Me%TAU_869.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_TAU_869
         me%ProdUnits(me%NProd)  = me%Units_NoDim
       endif
 
        call GetData(Me%EPS_78, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='EPS_78', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR07 .'

       if (Me%EPS_78.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_EPS_78
         me%ProdUnits(me%NProd) = me%Units_NoDim
       endif
 
        call GetData(Me%ANGSTROM_510, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='ANGSTROM_510', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR08 .'

       if (Me%ANGSTROM_510.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_ANGSTROM_510
         me%ProdUnits(me%NProd) = me%Units_NoDim
       endif
 
        call GetData(Me%ANGSTROM_531, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='ANGSTROM_531', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR09 .'

       if (Me%ANGSTROM_531.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_ANGSTROM_531
         me%ProdUnits(me%NProd) = me%Units_NoDim
       endif
 
        call GetData(Me%NLW_412, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_412', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR010 .'

       if (Me%NLW_412.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_412
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_443, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_443', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR011 .'

       if (Me%NLW_443.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_443
         me%ProdUnits(me%NProd) = me%Units_nLw
       endif
 
        call GetData(Me%NLW_490, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_490', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR012 .'

       if (Me%NLW_490.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_490
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_510, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_510', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR013 .'

       if (Me%NLW_510.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_510
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_555, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_555', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR014 .'

       if (Me%NLW_555.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_555
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_670, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_670', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR015 .'

       if (Me%NLW_670.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_670
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_488, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_488', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR016 .'

       if (Me%NLW_488.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_488
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_531, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_531', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR017 .'

       if (Me%NLW_531.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_531
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_551, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_551', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR018 .'

       if (Me%NLW_551.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_551
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif
 
        call GetData(Me%NLW_667, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='NLW_667', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR019 .'

       if (Me%NLW_667.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_NLW_667
         me%ProdUnits(me%NProd)  = me%Units_nLw
       endif

	   
        call GetData(Me%CALCITE, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='CALCITE', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR020 .'

       if (Me%CALCITE.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_CALCITE
         me%ProdUnits(me%NProd)  = me%Units_calcite
       endif

	   
	   call GetData(Me%TSM_CLARK, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='TSM_CLARK', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR021 .'

       if (Me%TSM_CLARK.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_TSM_CLARK
         me%ProdUnits(me%NProd)  = me%Units_TSM
       endif

	   call GetData(Me%POC_CLARK, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='POC_CLARK', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR022 .'

       if (Me%POC_CLARK.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_POC_CLARK
         me%ProdUnits(me%NProd)  = me%Units_POC
       endif

	   call GetData(Me%CHL_CLARK, &
               Me%ObjEnterData, iflag, &
               SearchType   = FromBlock,  keyword='CHL_CLARK', &
               default    = 0 , &
               ClientModule = 'ConvertOceanColorL2', &
               STAT       = STAT_CALL)
 
       if (STAT_CALL .NE. SUCCESS_) &
           Stop 'Construct_FieldsList - ModuleConvertOceanColorL2 ERR022 .'

       if (Me%POC_CLARK.eq.1) then
         me%NProd = me%NProd + 1   
         me%ProdNames(me%NProd) = me%Char_CHL_CLARK
         me%ProdUnits(me%NProd)  = me%Units_chla
       endif
       
	   allocate (me%ProdId(me%NProd))
        
       
       

        
do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData,                                &
                                        ClientNumber           = ClientNumber,         &
                                        block_begin            = '<<begin_field>>',    &
                                        block_end              = '<<end_field>>',      &
                                        BlockInBlockFound      = BlockFound,           &
                                        STAT                   = STAT_CALL)
cd1 :       if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :           if (BlockFound) then                                                  
                    ! Construct a New Field
                    Call Construct_Field(NewField)

                    ! Add new Property to the Fields List 
                    Call Add_Field(NewField)
               
                    call ReadFileHDF(NewField)
               
                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR02'
            else cd1
                stop 'Construct_FieldsList - ModuleConvertOceanColorL2 - ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_FieldsList

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new Field       

    subroutine Construct_Field(NewField)

        !Arguments-------------------------------------------------------------
 
         type(T_OceanColorL2), pointer       :: NewField
  
        !Local ----------------------------------------------------------------
         integer                               :: STAT_CALL
         integer                               :: iflag


        !----------------------------------------------------------------------
             
        allocate (NewField)

        nullify(NewField%Scalar        )
        nullify(NewField%Prev,NewField%Next)


         call GetData(NewField%L2FileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'INPUTFILE_L2',                     &
                     ClientModule = 'ConvertOceanColorL2',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertOceanColorL2 - ERR03'

        if (iflag == 0)then
            write(*,*)'Must specify name of L2 File to convert'
            stop 'Construct_Field - ModuleConvertOceanColorL2 - ERR04'
        end if


    end subroutine Construct_Field

    !--------------------------------------------------------------------------
    

    ! This subroutine adds a new Field to the Fields List  

    subroutine Add_Field(NewField)

        !Arguments-------------------------------------------------------------
        type(T_OceanColorL2),           pointer     :: NewField

        !----------------------------------------------------------------------

        ! Add to the WaterProperty List a new property
        if (.not.associated(Me%FirstField)) then
            Me%FieldsNumber = 1
            Me%FirstField    => NewField
            Me%LastField     => NewField
        else
            NewField%Prev                     => Me%LastField
            Me%LastField%Next => NewField
            Me%LastField      => NewField
            Me%FieldsNumber   =  Me%FieldsNumber + 1
        end if 

        !----------------------------------------------------------------------

    end subroutine Add_Field 

        !----------------------------------------------------------------------                     

    subroutine ReadFileHDF (NewField)

    !Types -------------------------------------------------------------
        type(T_OceanColorL2),  pointer                    :: NewField
        type(T_Polygon),  pointer                    :: OriginalArea
        type(T_PointF),   pointer                    :: AuxPoint
        type T_selectedArea
           type(T_PointF)                            :: pointF
           type(T_Point)                             :: point
        endType

        type(T_selectedArea), dimension(:),  pointer :: UserArea


        
    !----------------------------------------------------------------------

    real, dimension(:, :), pointer                :: XX_IE,YY_IE, LatitudeConn, LongitudeConn
    integer, dimension(:, :), pointer             :: DefineCellsMap
    character(len=256)                            :: L2FileName
    character(len=40)                             :: auxprod
    character(len=40), dimension(16)              :: FlagNames
    real, dimension(6)              :: AuxTime
    character(len=256)              :: CharTime, auxchar, ctime, cdate
    REAL*4  lat    
    REAL*4  lon  
    REAL*4, dimension(:,:), pointer :: LatArray
    REAL*4, dimension(:,:), pointer :: LonArray
    REAL*4, dimension(:,:), pointer :: GridData
	REAL*4  sdata 
	REAL*4  aux
	INTEGER flag  
    INTEGER i 
	INTEGER j 
	REAL*4  ullat, ullon, urlat, urlon
    REAL*4  lllat, lllon, lrlat, lrlon
	REAL*4  northlat,	southlat, westlon,eastlon
	INTEGER lines,columns
    LOGICAL outpoint(4), pointin
    integer npoint, line, column, maxj,minj,maxi,mini,limite,nProd,prodID
    real, dimension(:,:), pointer           :: GridLat
    real, dimension(:,:), pointer           :: GridLong
    integer*4 status,vg_ref,vgid, file_info_status, sds_id, FileID,Hclose
    integer*4 vfdtch, Sfend, vfend
    integer*4 hopen, sfstart, n_datasets, n_file_attributes, sffinfo,Vstart
    integer*4 vfind,vfatch, Vntagrefs, tag, ref, Vgettagref, SDreftoindex
    integer*4 sfselect,rank, dtype, nattrs,sfginfo,sfendacc,vfstart,vntrc,vfgttr
    integer*4 sfref2index, test,sd_id
    character(40) buffer
    integer, dimension(3)                   :: dims
    character*1 access 
    integer STAT_CALL
    integer geointerp, geoCorrection
    integer iub,ilb,jub,jlb,prodFound,contaId
    ! ----------------------------------------------------------------------

 FlagNames(1) = 'EPSILON1'
 FlagNames(2) = 'LAND1'
 FlagNames(3) = 'ANCIL1'
 FlagNames(4) = 'SUNGLINT1'
 FlagNames(5) = 'HIGHLT1'
 FlagNames(6) = 'SATZEN1'
 FlagNames(7) = 'COASTZ1'
 FlagNames(8) = 'NEGLW1'
 FlagNames(9) = 'STRAYLIGHT1'
 FlagNames(10)= 'CLDICE1'
 FlagNames(11)= 'COCCOLITH1'
 FlagNames(12)= 'TRUBIDW1'
 FlagNames(13)= 'SOLZEN1'
 FlagNames(14)= 'HIGHTAU1'
 FlagNames(15)= 'LOWLW1'
 FlagNames(16)= 'CHLOR1'


    !----------------------------------------------------------------------

    L2FileName    = NewField%L2FileName
    
    
    allocate(OriginalArea)
    allocate(OriginalArea%VerticesF(1:5))
    allocate (NewField%ProdFound(me%NProd))
    
    nullify (UserArea)
    allocate(UserArea(1:5))

    nullify (Auxpoint)
    allocate(Auxpoint)

    FileID      = hopen (L2FileName, 1, 0)

    write(*,*) 'Processing  file: '// trim(adjustl(L2FileName))
    write(*,*)

    if (FileID.eq. -1) then
      write(*,*)'Error File not Found!!'
    endif

    sd_id      = sfstart   (L2FileName, 1)

    file_info_status = sffinfo(FileID, n_datasets, n_file_attributes)
    
    status = vfstart(FileID);
    vg_ref = vfind(FileID, 'Geophysical Data') 
    vgid = vfatch(FileID, vg_ref, 'r')
   
    
       NewField%ProdFound=0
       
       do j=1, me%Nprod
        prodFound=0
        i=0
        test=0
        contaId = 0
        do while (test.ne.-1)
          
          test = vfgttr(vgid, i, tag, ref)  
          sds_id = sfselect(sd_id, sfref2index(sd_id, ref))
          status = sfginfo(sds_id, buffer, rank, dims, dtype, nattrs)
          status = sfendacc(sds_id)
        
          if (test.ne.-1) then
            if(buffer.eq.me%ProdNames(j)) then
              me%ProdID(j)=contaId
              prodFound=1
            endif
           endif
        
           if (buffer.ne.me%Char_l2_flags) then
             contaId=contaId +1
           endif
           i=i+1
        enddo

        if (prodFound.eq.0)then
         write(*,*) ' Product '// me%ProdNames(j) // ' not found!'
         NewField%ProdFound(j)=0
         me%ProdID(j)=-1
        else
         NewField%ProdFound(j)=1
        endif
         
    enddo

    status = Sfend(sd_id);
    status = vfdtch(vgid);
    status = vfend(FileID);

    status = Hclose(FileID);

    
 ! Geo File Get Georeference   ---------------------------------------------
  
      call OpenL2Seadas (trim(adjustl(L2FileName))//''C)
	  call GetMetaData( ullat, ullon, urlat, urlon, &
                        lllat, lllon, lrlat, lrlon,	&
                        northlat,	southlat, westlon, eastlon, columns, &
                        lines,nProd,geointerp)
	  
      iub = lines - 2
      ilb = 2
      jub = columns -4
      jlb = 4
      
      if (ullat.lt.lllat) then
       aux=lllat
       lllat=ullat
       ullat=aux
      endif

      if (urlat.lt.lrlat) then
       aux=lrlat
       lrlat=urlat
       urlat=aux
      endif
	  
      if (urlon.lt.ullon) then
       aux=urlon
       urlon=ullon
       ullon=aux
      endif

      if (lrlon.lt.lllon) then
       aux=lllon
       lllon=lrlon
       lrlon=aux
      endif
    
      call ReadFirstL2()

! Get date ----------------------------------------------------------------
   
   call ReadDate (NewField, L2FileName, AuxTime, cdate, ctime)

! -------------- Get Lat and Lon --------------------------

      nullify (LatArray)
	  nullify (LonArray)
	  allocate(LatArray(ilb:iub,jlb:jub))
      allocate(LonArray(ilb:iub,jlb:jub))
    
      do i=ilb,iub
	   do j=jlb,jub
        
        call GetLatLon(j,i,LatArray(i,j),LonArray(i,j))

	   enddo
      enddo


if (Me%NoExtraction.eq.0) then
! Get Area ------------------------------------------------
   
   OriginalArea%Count = 5


                OriginalArea%VerticesF(1)%X = LonArray(ilb,jlb) 
                OriginalArea%VerticesF(2)%X = LonArray(iub,jlb) 
                OriginalArea%VerticesF(3)%X = LonArray(iub,jub) 
                OriginalArea%VerticesF(4)%X = LonArray(ilb,jub)
                OriginalArea%VerticesF(5)%X = LonArray(ilb,jlb) 

                OriginalArea%VerticesF(1)%Y = LatArray(ilb,jlb)
                OriginalArea%VerticesF(2)%Y = LatArray(iub,jlb)
                OriginalArea%VerticesF(3)%Y = LatArray(iub,jub)  
                OriginalArea%VerticesF(4)%Y = LatArray(ilb,jub)
                OriginalArea%VerticesF(5)%Y = LatArray(ilb,jlb)

! --- Confirm that the selected area is under the original one

   call SetLimits(OriginalArea)

   UserArea(1)%PointF%X = Me%MinLong
   UserArea(1)%PointF%Y = Me%MinLat

   UserArea(2)%PointF%X = Me%MinLong
   UserArea(2)%PointF%Y = Me%MaxLat

   UserArea(3)%PointF%X = Me%MaxLong 
   UserArea(3)%PointF%Y = Me%MinLat

   UserArea(4)%PointF%X = Me%MaxLong
   UserArea(4)%PointF%Y = Me%MaxLat

do npoint=1, 4

   Auxpoint%X=UserArea(npoint)%PointF%X
   Auxpoint%Y=UserArea(npoint)%PointF%Y
   outpoint(npoint)=.false.
   
   pointin = IsPointInsidePolygon(Auxpoint, OriginalArea)
   
   if (pointin.eq..false.) then    
       outpoint(npoint)=.true.
   endif

 enddo 


  if (outpoint(1).AND.outpoint(2).AND.outpoint(3).AND.outpoint(4)) then
   Stop "All User Points are outside the image"
  endif

  do npoint=1,4
    if (outpoint(npoint))then
       
         Auxpoint%X=UserArea(npoint)%PointF%X
                     

        if ((Auxpoint%X.gt.urlon).or.(Auxpoint%X.gt.lrlon)) then
            UserArea(npoint)%Point%J= jub
        else if ((Auxpoint%X.lt.ullon).or.(Auxpoint%X.lt.lllon)) then
            UserArea(npoint)%Point%J= jlb
        else
                        !isto é para ver qual é a coluna tento em conta que o j ta dentro dos limites
                        Auxpoint%Y=(min(ullat,urlat)+max(lllat,lrlat))/2
                        
                        do column = jlb,jub-1
 
                       OriginalArea%VerticesF(1)%X = LonArray(ilb, column) 
                       OriginalArea%VerticesF(2)%X = LonArray(ilb, column +1) 
                       OriginalArea%VerticesF(3)%X = LonArray(iub, column +1) 
                       OriginalArea%VerticesF(4)%X = LonArray(iub, column)
                       OriginalArea%VerticesF(5)%X = LonArray(ilb, column)            
   
                       OriginalArea%VerticesF(1)%Y = LatArray(ilb, column) 
                       OriginalArea%VerticesF(2)%Y = LatArray(ilb, column +1)
                       OriginalArea%VerticesF(3)%Y = LatArray(iub, column +1)  
                       OriginalArea%VerticesF(4)%Y = LatArray(iub, column)
                       OriginalArea%VerticesF(5)%Y = LatArray(ilb, column)   
                 
                        call SetLimits(OriginalArea)

                           
                         if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
              
                                        UserArea(npoint)%Point%J= column
                                        exit

                         endif            

                 enddo
            
        endif

       !a linha 1 é a mais a norte!!
        Auxpoint%Y=UserArea(npoint)%PointF%Y 
        if ((Auxpoint%Y.gt.urlat).or.(Auxpoint%Y.gt.ullat)) then
            UserArea(npoint)%Point%I= ilb
        else if ((Auxpoint%Y.lt.lllat).or.(Auxpoint%Y.lt.lrlat)) then
            UserArea(npoint)%Point%I= iub
        else

        !isto é para ver qual é a linha tendo em conta que o i ta dentro dos limites
                        Auxpoint%X=(min(urlon,lrlon)+max(lllon,ullon))/2

             
             do line = ilb, iub-1
 
                OriginalArea%VerticesF(1)%X = LonArray(line  ,jlb) 
                OriginalArea%VerticesF(2)%X = LonArray(line+1,jlb) 
                OriginalArea%VerticesF(3)%X = LonArray(line+1,jub) 
                OriginalArea%VerticesF(4)%X = LonArray(line  ,jub)
                OriginalArea%VerticesF(5)%X = LonArray(line  ,jlb) 

                OriginalArea%VerticesF(1)%Y = LatArray(line  ,jlb) 
                OriginalArea%VerticesF(2)%Y = LatArray(line+1,jlb)
                OriginalArea%VerticesF(3)%Y = LatArray(line+1,jub)  
                OriginalArea%VerticesF(4)%Y = LatArray(line  ,jub)
                OriginalArea%VerticesF(5)%Y = LatArray(line  ,jlb) 


                call SetLimits(OriginalArea)

                if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
                        UserArea(npoint)%Point%I= line
                        exit
                endif

              enddo

        endif
     endif   
  enddo


 !Get i,j for each user point ------------------------------------------------------
 
  OriginalArea%Count = 5

 
 do line = ilb, iub-1
 
                OriginalArea%VerticesF(1)%X = LonArray(line  ,jlb) 
                OriginalArea%VerticesF(2)%X = LonArray(line+1,jlb) 
                OriginalArea%VerticesF(3)%X = LonArray(line+1,jub) 
                OriginalArea%VerticesF(4)%X = LonArray(line  ,jub)
                OriginalArea%VerticesF(5)%X = LonArray(line  ,jlb) 

                OriginalArea%VerticesF(1)%Y = LatArray(line  ,jlb) 
                OriginalArea%VerticesF(2)%Y = LatArray(line+1,jlb)
                OriginalArea%VerticesF(3)%Y = LatArray(line+1,jub)  
                OriginalArea%VerticesF(4)%Y = LatArray(line  ,jub)
                OriginalArea%VerticesF(5)%Y = LatArray(line  ,jlb)

                call SetLimits(OriginalArea)

   
   do npoint=1, 4
    if (.not.outpoint(npoint))then

      Auxpoint%X=UserArea(npoint)%PointF%X
      Auxpoint%Y=UserArea(npoint)%PointF%Y

       if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
              
                
                UserArea(npoint)%Point%I= line
       
           
                do column = jlb, jub-1



                       OriginalArea%VerticesF(1)%X = LonArray(max(ilb,line - 100), column) 
                       OriginalArea%VerticesF(2)%X = LonArray(max(ilb,line - 100), column +1) 
                       OriginalArea%VerticesF(3)%X = LonArray(min(iub,line + 100), column +1) 
                       OriginalArea%VerticesF(4)%X = LonArray(min(iub,line + 100), column)
                       OriginalArea%VerticesF(5)%X = LonArray(max(ilb,line - 100), column)            
   
                       OriginalArea%VerticesF(1)%Y = LatArray(max(ilb,line - 100), column) 
                       OriginalArea%VerticesF(2)%Y = LatArray(max(ilb,line - 100), column +1)
                       OriginalArea%VerticesF(3)%Y = LatArray(min(iub,line + 100), column +1)  
                       OriginalArea%VerticesF(4)%Y = LatArray(min(iub,line + 100), column)
                       OriginalArea%VerticesF(5)%Y = LatArray(max(ilb,line - 100), column) 
                 
                        call SetLimits(OriginalArea)

                           
                         if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
              
                                        UserArea(npoint)%Point%J= column
                                        exit

                         endif
                          
                            

                 enddo
       
       endif
     endif
    enddo 

 enddo

 MaxJ = max(UserArea(1)%Point%J, &
            UserArea(2)%Point%J, &
            UserArea(3)%Point%J, &
            UserArea(4)%Point%J)
 
 MinJ = min(UserArea(1)%Point%J, &
            UserArea(2)%Point%J, &
            UserArea(3)%Point%J, &
            UserArea(4)%Point%J)


 MaxI = max(UserArea(1)%Point%I, &
            UserArea(2)%Point%I, &
            UserArea(3)%Point%I, &
            UserArea(4)%Point%I)

 MinI = min(UserArea(1)%Point%I, &
            UserArea(2)%Point%I, &
            UserArea(3)%Point%I, &
            UserArea(4)%Point%I)

 if (min(MaxJ,MinJ,MaxI,MinI).lt.0) then
   write(*,*)  'Atention! Error reading line and column'
   stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR01'
 endif
else

    MaxJ=columns - 5
    MinJ=5
    MinI=3
    MaxI=lines - 3

endif !no Extraction!!!
! Generate Grid ---------------------------------------------------------

if (mini.le.2) then
 minI =3
endif

if (minj.le.4) then
 minJ =5
endif

if (maxj.ge.(columns-4)) then
 maxj =columns -5
endif

if (maxi.eq.(lines-2)) then
 maxi = line - 3
endif
     
     NewField%Size%ILB  = 0
     NewField%Size%IUB  = maxI - minI + 2
     NewField%Size%JLB  = 0
     NewField%Size%JUB  = maxJ -minJ + 2

     NewField%WorkSize%ILB = 1
     NewField%WorkSize%IUB = maxI - minI + 1
     NewField%WorkSize%JLB = 1
     NewField%WorkSize%JUB = maxJ - minJ + 1


Allocate (GridLat                (NewField%Size%ILB:NewField%Size%IUB, NewField%Size%JLB:NewField%Size%JUB))
Allocate (GridLong               (NewField%Size%ILB:NewField%Size%IUB, NewField%Size%JLB:NewField%Size%JUB))
     
Allocate (NewField%Scalar        (me%Nprod, NewField%Size%ILB:NewField%Size%IUB , NewField%Size%JLB:NewField%Size%JUB ))
Allocate (NewField%OpenPoints3D  (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB,1))
Allocate (NewField%XX_IE         (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
Allocate (NewField%YY_IE         (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
        
Allocate(NewField%LatitudeConn   (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
Allocate(NewField%LongitudeConn  (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))

Allocate(NewField%DefineCellsMap (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))

Allocate(Me%Bathymetry           (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))

Allocate (XX_IE                  (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
Allocate (YY_IE                  (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
        
Allocate(LatitudeConn            (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
Allocate(LongitudeConn           (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))

Allocate(DefineCellsMap          (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))

Allocate(GridData                (NewField%Size%ILB : NewField%Size%IUB, NewField%Size%JLB : NewField%Size%JUB))
        
        GridLat  = null_real 
        GridLong = null_real
        NewField%Scalar   = null_real
        NewField%XX_IE    = null_real
        NewField%YY_IE    = null_real
        Me%Bathymetry     = null_real

        NewField%LatitudeConn  = null_real
        NewField%LongitudeConn = null_real
        
        NewField%DefineCellsMap(:,:)  = 1


Line  = 1
Column =1


do i = minI , maxI+1
 do j= minJ , maxJ+1 
   
   GridLat   (Line,Column) = (LatArray(i-1,j)+LatArray(i,j))/2
   GridLong  (Line,Column) = (LonArray(i,j-1)+LonArray(i,j))/2
   
   Column = Column + 1 
 
 enddo
  Column = 1
  Line   = Line + 1
enddo

deallocate (LatArray)
deallocate (LonArray)
nullify(LatArray)
nullify(LonArray)


! -------- Reads Geophysical Data ----------------------------------------  
     

    do Nprod=1,me%Nprod
     if (NewField%ProdFound(Nprod).eq.1) then
      prodID = me%prodID(Nprod)
     

        call ReadDataL2(minI,maxI, prodID)
	        column=1
            line=1
	     
         
         do i=minI , maxI
          do j = minJ , maxJ
              call GetDataL2(j,i,prodID, NewField%Scalar(Nprod,line,column),flag)
              column = column + 1
              NewField%OpenPoints3D(line,column,1) = 1
              Me%Bathymetry(line,column)=NewField%Scalar(Nprod,line,column-1)
             
             ! if (NewField%Scalar(Nprod,line,column).le.-1) then
              !   NewField%OpenPoints3D(line,column,1) = 0
              !endif 
	      enddo
	          column = 1
              line   = line + 1
         enddo
      endif
    enddo
    
      

! -------- Writes Grid ----------------------------------------------------
    
      
      
      NewField%Gridfilename = trim(adjustl(NewField%L2FileName))//'.dat'

      call WriteGridData(FileName       = NewField%Gridfilename,             &
                         COMENT1        = 'OC L2 file'//trim(adjustl(L2Filename)), &
                         COMENT2        = trim(adjustl(cdate))//' '//trim(adjustl(ctime))//' '//'Chl a mg/m^3',     &
                         ConnectionX    = GridLong,                                                                 &
                         ConnectionY    = GridLat,                                                                  &
                         WorkSize       = NewField%WorkSize,                                                                 &
                         CoordType      = 4,                                                                        &
                         Xorig          = gridlong(1,NewField%Size%jub),                                  &
                         Yorig          = gridlat(1,1),                                   &
                         Zone           = 29,                                                                       &
                         GRID_ANGLE     = 0.,                                                                       &
                         Latitude       = gridlat(1,1),                                  &
                         Longitude      = gridlong(1,NewField%Size%jub),                                  &
                         FillValue      = -99.,                                                                     &
                         Overwrite      =.true.,                                                                   &
                         GridData2D_Real= Me%Bathymetry,                                                                     &
                         STAT           = STAT_CALL) 
        if(STAT_CALL .ne. SUCCESS_)stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR02'
     
     call ConstructHorizontalGrid(Me%ObjHorizontalGrid, NewField%Gridfilename, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR03'
    

     call GetHorizontalGrid(Me%ObjHorizontalGrid, XX_IE , YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR04'

     
     call GetGridLatitudeLongitude(Me%ObjHorizontalGrid, LatitudeConn, LongitudeConn,STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR05'

     call GetDefineCellsMap(Me%ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR06'
       
      NewField%XX_IE          = XX_IE
      NewField%YY_IE          = YY_IE
      NewField%LatitudeConn   = LatitudeConn
      NewField%LongitudeConn  = LongitudeConn
      NewField%DefineCellsMap = DefineCellsMap 
     
     call UnGetHorizontalGrid (Me%ObjHorizontalGrid, XX_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR07'

     call UnGetHorizontalGrid (Me%ObjHorizontalGrid, YY_IE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR08'

     call UnGetHorizontalGrid (Me%ObjHorizontalGrid, LatitudeConn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR09'

     call UnGetHorizontalGrid (Me%ObjHorizontalGrid, LongitudeConn, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR010'

     call UnGetHorizontalGrid (Me%ObjHorizontalGrid, DefineCellsMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR011'
         
     call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ReadFileHDF - ModuleConvertOceanColorL2 - ERR012'

     deallocate (GridData)
     deallocate (GridLat)
     deallocate (GridLong)
     
     
     call closeFile
  
 end subroutine ReadFileHDF

 !------------------------------------------------------------------------


 subroutine ReadDate (NewField, Filename, AuxTime, cdate, ctime)

    !Arguments-----------------------------------------------------------
    type(T_Attribute), dimension(:), pointer     :: Attributes
    type(T_Time)                                 :: Date
    type(T_OceanColorL2),  pointer               :: NewField
    
    character(len=256), intent (IN):: Filename
    real, dimension(6),intent (OUT):: AuxTime
    character(len=256),intent (OUT):: cdate, ctime
    character(len=256)             :: chartime, charday, charyear, buf,aux
    character(len=4)               :: char_i , cmonth, cday                   
    integer(4)                     :: FileID,  n_file_attributes, file_attr_index, file_attr_status
    integer(4)                     :: read_attr_status, dayj, AuxDayj, count,sd_id,sfstart
    logical                        :: FileFound
    integer(4)                     :: sffattr, sfgainfo, attr_index, dtype, day, month, i
    integer(4)                     :: status,sfrcatt,Sfend,hclose,vfend,file_info_status
    integer(4)                     :: n_datasets,hopen,sffinfo
    !---------------------------------------------------------------------
 

    FileID      = hopen(FileName, 1, 0)


    if (FileID.eq.-1) then
    
       write(*,*) 'Subroutine ReadDate ModuleConvertOceanColorL2 ERR01'
       stop 'File not found'
       
    endif
            
    sd_id      = sfstart   (Filename, 1)
    
    aux = 'Start Time'
    do i=0,200
      status = sfgainfo(sd_id, i, buf, dtype, count)
      
      if (buf.eq.aux) then
         status = sfrcatt(sd_id, i, charTime)
         exit
      endif

      if (status.ne.0) then
        exit
      endif 
    
    enddo 
  
      
      read(CharTime(1:4), '(i4)') i
      AuxTime(1) = i

      read(CharTime(8:9), '(i4)') i
      AuxTime(4) = i

      read(CharTime(10:11), '(i4)') i
      AuxTime(5) = i

      read(CharTime(12:13), '(i4)') i
      AuxTime(6) = i
    
      read(CharTime(5:7), '(i4)') i
      dayj = i
      
      call GetDateFromJDay (AuxTime(1),AuxTime,1, real(dayj))
         
    ctime = CharTime(8:9)//'-'//CharTime(10:11)//'-'//CharTime(12:13)
    write(cmonth, '(i4)')int(AuxTime(2))
    write(cday, '(i4)')  int(AuxTime(3))
    cdate = CharTime(1:4)//'-'//trim(adjustl(cmonth))//'-'//trim(adjustl(cday))

    call SetDate(NewField%Date, Year = AuxTime(1), Month  = AuxTime(2), Day    = AuxTime(3), &
                  Hour = AuxTime(4), Minute = AuxTime(5), Second = AuxTime(6))

    status = Sfend(sd_id);
    status = Hclose(FileID);

   end subroutine
!------------------------------------------------------------------------
 subroutine GetDateFromJDay (Year, AuxTime, StartInZero, JulianDay )
   
    real              , intent (IN)   ::  Year
	integer           , intent (IN)   ::  StartInZero
    real, dimension(6), intent (OUT)  ::  AuxTime
	real              , intent (IN)   ::  JulianDay

	integer                         ::  iMonth,iDay, JulianDayInt,Month
	integer                         ::  Day,Hour,Minute,Second
	real                            ::  Rest
    integer, Dimension(12)          ::  NDay

    Data NDay/31,28,31,30,31,30,31,31,30,31,30,31/

    !if StartInZero = 0 then 12:00 1st January is Julian 0.5
    !if StartInZero = 1 then 12:00 1st January is Julian 1.5
        

    JulianDayInt= int(JulianDay)
    
    !Leap year
    !lp = ((y % 4 = 0) AND (y % 100 <> 0)) OR (y % 400 = 0) 
    !lp means leap year. 
    !y is the year you want to calculate. 
    !% means take the remainder of. 
    !Corrected Frank 8-1999
    !Recorrected 10-2000
        if ( (mod(Year,4.) == 0 .and. mod(Year,100.) /= 0) .or.         &
             (mod(Year,400.) == 0))  then
              NDay(2)=29
        else
              NDay(2)=28
        endif
    
    Day = 0
    Month = 1
    
    do iDay = StartInZero,JulianDayInt
    
        Day = Day + 1
        
        If (Day.gt.NDay(Month)) Then
            Month = Month + 1
            Day = 1
        End If
        
    enddo
    
	AuxTime(2)=real(Month)
	AuxTime(3)=real(day)

    !JulianToDate = CDate(Str(Year) + "-" + Str(Month) + "-" + Str(Day) + " " + _
     !                    Str(Hour) + ":" + Str(Minute) + ":" + Str(Second))

   
 end subroutine GetDateFromJDay
! ------------------------------------------------------------------------

 subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertOceanColorL2 - ERR01'
                
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertOceanColorL2 - ERR07'

    end subroutine Open_HDF5_OutPut_File
 
   !------------------------------------------------------------------------
    

    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        real,    dimension(:,:), pointer                :: Scalar2D
        integer                                         :: STAT_CALL, OutputNumber, nprod
        type (T_OceanColorL2), pointer                  :: NewField
        type(T_Time)                                    :: CurrentDate

        !Begin-----------------------------------------------------------------
        
             
        call Open_HDF5_OutPut_File

        OutputNumber = 1
        CurrentDate  = Me%FirstField%Date
        
        call ExtractDate   (CurrentDate,                                &
                            AuxTime(1), AuxTime(2), AuxTime(3),         &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR01'


        call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                             "Time", "YYYY/MM/DD HH:MM:SS",             &
                             Array1D = TimePtr,                         &
                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR02'


        call HDF5SetLimits  (Me%ObjHDF5, Me%FirstField%WorkSize%ILB, Me%FirstField%WorkSize%IUB - 1,&
                             Me%FirstField%WorkSize%JLB, Me%FirstField%WorkSize%JUB - 1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR03'

        
        Allocate (Scalar2D(Me%FirstField%Size%ILB:Me%FirstField%Size%IUB , Me%FirstField%Size%JLB:Me%FirstField%Size%JUB ))
        Scalar2D(:,:) = Me%FirstField%Scalar(1,:,:)

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D = Scalar2D , STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertModisL3Mapped - ERR04'            

        Deallocate (Scalar2D)
        
        NewField => Me%FirstField

        do while(associated(NewField))

        Write(*,*) 'Writing to HDF5 file: '// trim(NewField%L2FileName)

         
            

            if(NewField%Date .ne. CurrentDate)then

               

                CurrentDate = NewField%Date

                OutputNumber = OutputNumber + 1


                call ExtractDate   (CurrentDate,                                &
                                    AuxTime(1), AuxTime(2), AuxTime(3),         &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime
            
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR05'
                

                call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",             &
                                     Array1D = TimePtr,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR06'

            
            end if


            if(NewField%Date .eq. CurrentDate)then


                 call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB,&
                               NewField%WorkSize%JLB, NewField%WorkSize%JUB, STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR09'

                call HDF5WriteData   (me%ObjHDF5, "/Grid/ConnectionX", "ConnectionX", "m",          &
                                      Array2D = NewField%XX_IE,                                  &
                                      OutputNumber = OutputNumber,                               &
                                      STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'OutputFields - ModuleConvertOceanColorL2 - ERR011'


                call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB,&
                               NewField%WorkSize%JLB, NewField%WorkSize%JUB, STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR09'

               
                call HDF5WriteData   (me%ObjHDF5, "/Grid/ConnectionY", "ConnectionY", "m",          &
                                      Array2D = NewField%YY_IE,                                  &
                                      OutputNumber = OutputNumber,                               &
                                      STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'OutputFields - ModuleConvertOceanColorL2 - ERR012'

                 call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB,&
                               NewField%WorkSize%JLB, NewField%WorkSize%JUB, STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR09'

                
                
                call HDF5WriteData   (me%ObjHDF5, "/Grid/Longitude", "Longitude", "º",              &
                                      Array2D = NewField%LongitudeConn,                          &
                                      OutputNumber = OutputNumber,                               &
                                      STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'OutputFields - ModuleConvertOceanColorL2 - ERR013'

                 
                 call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB,&
                               NewField%WorkSize%JLB, NewField%WorkSize%JUB, STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR09'

                
                call HDF5WriteData   (me%ObjHDF5, "/Grid/Latitude", "Latitude", "º",                &
                                      Array2D = NewField%LatitudeConn,                           &
                                      OutputNumber = OutputNumber,                               &
                                      STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'OutputFields - ModuleConvertOceanColorL2 - ERR014'

               
                 call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB,&
                               NewField%WorkSize%JLB, NewField%WorkSize%JUB, STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR09'

               
                 call HDF5WriteData   (me%ObjHDF5, "/Grid/Define_Cells", "Define Cells", "-",    &
                                      Array2D = NewField%DefineCellsMap,                     &
                                      OutputNumber = OutputNumber,                           &
                                      STAT = STAT_CALL)
                   if (STAT_CALL /= SUCCESS_) stop 'OutputFields - ModuleConvertOceanColorL2 - ERR015'
   
        
                 Allocate (Scalar2D(NewField%Size%ILB:NewField%Size%IUB , NewField%Size%JLB:NewField%Size%JUB ))


                 call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB-1,&
                                      NewField%WorkSize%JLB, NewField%WorkSize%JUB-1, 1 , 1, STAT = STAT_CALL)
                  if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR016'

                 call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",       &
                 
                                     "-", Array3D = NewField%OpenPoints3D,               &
                                      OutputNumber = OutPutNumber, STAT = STAT_CALL)

                  do nprod=1,me%nprod
                    if (NewField%ProdFound(nprod).eq.1) then

                       Scalar2D(:,:) = NewField%Scalar(nprod,:,:)
                
                       call HDF5SetLimits  (Me%ObjHDF5, NewField%WorkSize%ILB, NewField%WorkSize%IUB-1,&
                                            NewField%WorkSize%JLB, NewField%WorkSize%JUB-1, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR017'
                
                    call HDF5WriteData(Me%ObjHDF5,                                      &
                                   "/Results/"//trim(adjustl(me%ProdNames(nprod))), &
                                   trim(adjustl(me%ProdNames(nprod))),              &
                                   trim(adjustl(me%ProdUnits(nprod))),              &
                                   Array2D      = Scalar2D,                         &
                                   OutputNumber = OutPutNumber,                     &
                                   STAT         = STAT_CALL)
                     if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR018'
             
             endif
            
            enddo

            endif
           
            deallocate(Scalar2D)

            NewField => NewField%Next

        end do

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertOceanColorL2 - ERR021'
         
        
   
    end subroutine OutputFields    
    
    subroutine KillOceanColorL2
        
        !Local-----------------------------------------------------------------
        integer                             :: STAT_CALL
        type (T_OceanColorL2), pointer      :: NewField
        
        !Begin-----------------------------------------------------------------

               
        
        deallocate (me%Bathymetry)
        deallocate (me%ProdNames)
        deallocate (me%ProdUnits)
        deallocate (me%ProdId)
         

           NewField => Me%FirstField

           do while(associated(NewField))

              deallocate (NewField%Scalar)
              deallocate (NewField%OpenPoints3D)
              deallocate (NewField%ProdFound)
                        
              deallocate (NewField%XX_IE)        
              deallocate (NewField%YY_IE)        
        
              deallocate(NewField%LatitudeConn)  
              deallocate(NewField%LongitudeConn) 

              deallocate(NewField%DefineCellsMap)
           
              NewField => NewField%Next

           end do


    end subroutine KillOceanColorL2


    !------------------------------------------------------------------------
 
end module ModuleConvertOceanColorL2










!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ConvertModisL2
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Pedro Pina - v4.0
! DESCRIPTION   : Module to convert ModisL3Mapped Format files into HDF5 format
! DOWNLOAD      : To download ModisL3Mapped files go to http://oceancolor.gsfc.nasa.gov/cgi/level3.pl
!------------------------------------------------------------------------------


Module ModuleConvertModisL2

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

  
    public  ::  StartConvertModisL2
    private ::  Construct_FieldsList
    private ::        Construct_Field
    private ::        Add_Field
    private ::  ReadHDF4
    private ::  KillModisL2
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
    
    private :: T_ModisL2
    type       T_ModisL2
        character(len=PathLength)                :: OutputFileName
        character(len=PathLength)                :: L2FileName
        character(len=PathLength)                :: GeoFileName
        character(len=StringLength)              :: ParameterName
        character(len=StringLength)              :: Units  
        integer                                  :: IDNumber
        type(T_Time)                             :: Date
        real, dimension(:,:),       pointer      :: Scalar
        type(T_ModisL2),      pointer      :: Next
        type(T_ModisL2),      pointer      :: Prev
        type(T_Attribute), dimension(:), pointer :: Attributes
    end type  T_ModisL2
    
    

   


    private :: T_ModisL2Format
    type       T_ModisL2Format
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjModisL3MappedGrid = 0
        integer                                 :: ObjNewGrid           = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
! Header parameters
        real                                    :: LatStep  
        real                                    :: LongStep 
        real                                    :: OrigLat  
        real                                    :: OrigLong 

        
     
        
        character(len=PathLength)               :: GridFileName
        logical                                 :: FirstFile            = .true. 
        integer                                 :: Unit
        real                                    :: MaxLong
        real                                    :: MinLong
        real                                    :: MaxLat
        real                                    :: MinLat
        integer                                 :: FieldsNumber
        real, dimension(:,:),  pointer          :: Bathymetry
        integer, dimension(:,:), pointer        :: Frequency
        real , dimension(:),  pointer           :: XX
        real , dimension(:),  pointer           :: YY
        type(T_ModisL2),      pointer           :: FirstField
        type(T_ModisL2),      pointer           :: LastField

                
       
    end type  T_ModisL2Format


    type(T_ModisL2Format),         pointer     :: Me

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartConvertModisL2(EnterDataID, ClientNumber,  STAT)

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

        Call KillModisL2

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'StartConvertModisL2 - ModuleConvertModisL2 - ERR01' 

        STAT = SUCCESS_

        deallocate(Me)
        nullify   (Me)



    end subroutine StartConvertModisL2
!---------------------
    subroutine str

 

    integer                 :: i
    character(len=4)        :: char_i

    char_i= '2004'

    read(char_i, '(i4)')i

    write(*,*) char_i
    write(*,*) i

    i = i + 1

    write(char_i, '(i4)')i

    write(*,*) char_i

    end subroutine
  ! ------------------------------------------------------------------------------
   
    subroutine Construct_FieldsList(ClientNumber)

        ! Arguments -----------------------------------------------------------

        integer,           intent(IN )              :: ClientNumber

        !External----------------------------------------------------------------
        integer                         :: STAT_CALL
        logical                         :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_ModisL2), pointer      :: NewField
        integer                              :: iflag
        !------------------------------------------------------------------------

       
       

        call GetData(Me%MaxLong,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MAX_LONG',                          &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR02'
        
        if (iflag == 0)then
            write(*,*)'Must specify MaxLong'
            stop 'Construct_FieldsList - ModuleEUCenterFormat - ERR06'
        end if


        call GetData(Me%MinLong,                                        &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MIN_LONG',                         &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR03'
        
        if (iflag == 0)then
            write(*,*)'Must specify MinLong'
            stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR07'
        end if
        
        
        call GetData(Me%MaxLat,                                         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MAX_LAT',                          &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR04'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR08'
        end if


        call GetData(Me%MinLat,                                         &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'MIN_LAT',                          &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleEUCenterFormat - ERR03'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR09'
        end if

        
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
                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR02'
            else cd1
                stop 'Construct_FieldsList - ModuleConvertModisL2 - ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_FieldsList

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new Field       

    subroutine Construct_Field(NewField)

        !Arguments-------------------------------------------------------------
 
         type(T_ModisL2), pointer       :: NewField
  
        !Local ----------------------------------------------------------------
         integer                               :: STAT_CALL
         integer                               :: iflag


        !----------------------------------------------------------------------
             
        allocate (NewField)

        nullify(NewField%Scalar        )
        nullify(NewField%Prev,NewField%Next)

        call GetData(NewField%OutputFileName,                           &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'OUTPUT_FILE',                      &
                     ClientModule = 'ConvertModisL2',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertModisL2 - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of Output file'
            stop 'Construct_Field - ModuleConvertModisL2 - ERR02'
        end if

         call GetData(NewField%L2FileName,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'INPUTFILE_L2',                     &
                     ClientModule = 'ConvertModisL2',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertModisL2 - ERR03'

        if (iflag == 0)then
            write(*,*)'Must specify name of L2 File to convert'
            stop 'Construct_Field - ModuleConvertModisL2 - ERR04'
        end if
         
         call GetData(NewField%GeoFileName,                             &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'INPUTFILE_GEO',                    &
                     ClientModule = 'ConvertModisL2',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertModisL2 - ERR05'

        if (iflag == 0)then
            write(*,*)'Must specify name of GEO file to convert'
            stop 'Construct_Field - ModuleConvertModisL2 - ERR06'
        end if

        call GetData(NewField%ParameterName,                            &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'PARAMETER',                        &
                     ClientModule = 'ConvertModisL2',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertModisL2 - ERR03'

        if (iflag == 0)then
            write(*,*)'Must specify parameter'
            stop 'Construct_Field - ModuleConvertModisL2 - ERR04'
        end if

        call GetData(NewField%Units,                                    &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'UNITS',                            &
                     ClientModule = 'ConvertModisL2',                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertModisL2 - ERR031'

        

        call ReadHDF4 (NewField)



    end subroutine Construct_Field

    !--------------------------------------------------------------------------
    

    ! This subroutine adds a new Field to the Fields List  

    subroutine Add_Field(NewField)

        !Arguments-------------------------------------------------------------
        type(T_ModisL2),           pointer     :: NewField

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

    subroutine ReadHDF4 (NewField)

    !Types -------------------------------------------------------------
        type(T_ModisL2),  pointer                    :: NewField
        type(T_Polygon),  pointer                    :: OriginalArea
        type(T_PointF),   pointer                    :: AuxPoint
        type T_selectedArea
           type(T_PointF)                            :: pointF
           type(T_Point)                             :: point
        endType

        type(T_selectedArea), dimension(:),  pointer :: UserArea

        type (T_Size2D)                              :: WorkSize
        type (T_Size2D)                              :: Size

		
    !----------------------------------------------------------------------
    integer(4)                              :: hopen, sfstart, istat, index, sfselect, sfgainfo
    integer(4)                              :: read_attr_status, file_attr_index, sfrcatt, sfrnatt
    integer(4)                              ::  sffinfo, n_file_attributes
    character(256), dimension(:), pointer   :: file_attributes, file_attr_names
    integer(4)                              :: FileID, Access, sds_id, sfginfo, sfrdata
    
    integer, dimension(3)                   :: dimsizes, limits
    integer(4), parameter                   :: LABLEN=64
    integer(4)                              :: rank, data_type, num_attrs
    integer(4), dimension(2)                :: start, stride
    character*1, dimension(:,:), allocatable:: Array2DChar
    
    
    real, dimension(:,:), pointer           :: auxArray2DReal
	real, dimension(:,:), pointer           :: auxArray2DReal2
    integer(2), dimension(:,:), pointer        :: auxArray2DInt
    
    real, dimension(:,:), pointer           :: Latitude
    real, dimension(:,:), pointer           :: Longitude
    real, dimension(:,:), pointer           :: Chla
    real, dimension(:,:), pointer           :: sst
    real, dimension(:,:), pointer           :: GridLat
    real, dimension(:,:), pointer           :: GridLong

    character(len=256)                      :: L2FileName, GeoFileName, ParameterName
     character(len=256)                     :: CharTime, auxchar, ctime, cdate
    character(len=LABLEN)                   :: sds_name
    integer                                 :: STAT_CALL, sfrcdata
    integer                                 :: file_attr_status, file_attr_type, file_attr_length
    real                                    ::  aux
    integer                                 :: nlines, ncolumns, AuxIMax, AuxJMax, AuxIMin, AuxJMin
    integer                                 :: i, j , dayj, AuxDayj, HDF5_CREATE
    real                                    :: base, slope, intercept
    real                                    :: month, day ,Year, maxOrigLat, minOrigLat
    real, dimension(6)                      :: AuxTime
    integer, dimension(255)                 :: array
    integer                                 :: num, paramIndex, sfn2index, npoint, npoint2
    logical                                 :: FileFound
    integer                                 :: line,column, lines, columns
    integer                                 :: maxI, maxJ, minI, minJ
    character(len=4)                        :: char_i, char_h
    !----------------------------------------------------------------------

    L2FileName    = NewField%L2FileName; L2FileName   = trim(adjustl(L2FileName))
    GeoFilename   = NewField%GeoFileName;GeoFilename  = trim(adjustl(GeoFilename))     
    
    allocate(OriginalArea)
    allocate(OriginalArea%VerticesF(1:5))
    
    nullify (UserArea)
    allocate(UserArea(1:5))

    nullify (Auxpoint)
    allocate(Auxpoint)
    
 ! Get date ----------------------------------------------------------------

   call ReadDate (L2FileName, AuxTime, cdate, ctime)

 ! Geo File Get Georeference   ---------------------------------------------
    FileFound=LookforFile (GeoFilename, FileID, n_file_attributes)
    
    if(FileFound.eq..false.) then 
                write(*,*)  
                write(*,*) 'Error! File Not Found: ', GeoFilename
                stop 'ReadHDF4 - ModuleConvertModisL2 - ERR01'
    endif

     !-- Longitude --------------------------------------------
     sds_name    = "Longitude"
     paramIndex  = sfn2index(FileID, sds_name)
     sds_id      = sfselect  (FileID, paramIndex)

     istat       = sfginfo   (sds_id, sds_name, rank, dimsizes, data_type, num_attrs)
    
     start       = 0
     stride      = 1

     !allocate to correct dimensions !ATTENTION remeber that HDF changes lines with columns
    allocate(Latitude(1:dimsizes(2),1:dimsizes(1)))
    allocate(Longitude(1:dimsizes(2),1:dimsizes(1)))
    allocate(AuxArray2Dreal(1:dimsizes(1),1:dimsizes(2)))
     
     !reads data
     istat       = sfrdata   (sds_id, start, stride, dimsizes, AuxArray2Dreal)

    do i=1 , dimsizes(2)
     do j = 1 , dimsizes(1)
       Longitude(i,j) = AuxArray2Dreal(j,i)
     enddo
    enddo
  
  !-- Latitude --------------------------------------------
     sds_name    = "Latitude"
     paramIndex  = sfn2index(FileID, sds_name)
     sds_id      = sfselect  (FileID, paramIndex)

     istat       = sfginfo   (sds_id, sds_name, rank, dimsizes, data_type, num_attrs)
    
     start       = 0
     stride      = 1
     
     
     !reads data
     istat       = sfrdata   (sds_id, start, stride, dimsizes, AuxArray2Dreal)

    do i=1 , dimsizes(2)
     do j = 1 , dimsizes(1)
       Latitude(i,j) = AuxArray2Dreal(j,i)
     enddo
    enddo           
     
    
 ! Get Area ------------------------------------------------
   
   OriginalArea%Count = 5

   OriginalArea%VerticesF(1)%X = longitude(1,1) 
   OriginalArea%VerticesF(2)%X = longitude(dimsizes(2),1) 
   OriginalArea%VerticesF(3)%X = longitude(dimsizes(2),dimsizes(1)) 
   OriginalArea%VerticesF(4)%X = longitude(1,dimsizes(1))
   OriginalArea%VerticesF(5)%X = longitude(1,1) 

   
   OriginalArea%VerticesF(1)%Y = latitude(1,1) 
   OriginalArea%VerticesF(2)%Y = latitude(dimsizes(2),1) 
   OriginalArea%VerticesF(3)%Y = latitude(dimsizes(2),dimsizes(1)) 
   OriginalArea%VerticesF(4)%Y = latitude(1,dimsizes(1)) 
   OriginalArea%VerticesF(5)%Y = latitude(1,1)

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

   if (.Not.IsPointInsidePolygon(Auxpoint, OriginalArea)) then
                write(*,*)  
                write(*,*) 'Error! Point Lat: ', Auxpoint%Y, ', Long: ',Auxpoint%X 
                write(*,*) 'Is Outside the original image:',GeoFilename 
                stop 'ReadHDF4 - ModuleConvertModisL2 - ERR01aa'
   endif
 
 enddo 
 
 !Get i,j for each point ------------------------------------------------------
 
  OriginalArea%Count = 5

 
 do line = 1, (dimsizes(2) -1)
 
   OriginalArea%VerticesF(1)%X = longitude(line,1) 
   OriginalArea%VerticesF(2)%X = longitude(line+1,1) 
   OriginalArea%VerticesF(3)%X = longitude(line+1,dimsizes(1)) 
   OriginalArea%VerticesF(4)%X = longitude(line,dimsizes(1))
   OriginalArea%VerticesF(5)%X = longitude(line,1) 

   
   OriginalArea%VerticesF(1)%Y = latitude(line,1) 
   OriginalArea%VerticesF(2)%Y = latitude(line+1,1) 
   OriginalArea%VerticesF(3)%Y = latitude(line+1,dimsizes(1))  
   OriginalArea%VerticesF(4)%Y = latitude(line,dimsizes(1)) 
   OriginalArea%VerticesF(5)%Y = latitude(line,1)


   call SetLimits(OriginalArea)

   
   do npoint=1, 4

      Auxpoint%X=UserArea(npoint)%PointF%X
      Auxpoint%Y=UserArea(npoint)%PointF%Y

       if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
              
                UserArea(npoint)%Point%I= line
       
           
                do column = 1, (dimsizes(1) -1)
 
                       OriginalArea%VerticesF(1)%X = longitude(line,column) 
                       OriginalArea%VerticesF(2)%X = longitude(line,column +1) 
                       OriginalArea%VerticesF(3)%X = longitude(line + 1, column +1) 
                       OriginalArea%VerticesF(4)%X = longitude(line + 1, column)
                       OriginalArea%VerticesF(5)%X = longitude(line, column)            
 

  
                       OriginalArea%VerticesF(1)%Y = latitude(line,column) 
                       OriginalArea%VerticesF(2)%Y = latitude(1,column +1) 
                       OriginalArea%VerticesF(3)%Y = latitude(line + 1, column +1)  
                       OriginalArea%VerticesF(4)%Y = latitude(line + 1, column) 
                       OriginalArea%VerticesF(5)%Y = latitude(line,column)
                 
                        call SetLimits(OriginalArea)

                         
                           
                         if (IsPointInsidePolygon(Auxpoint, OriginalArea)) then
              
                                        UserArea(npoint)%Point%J= column
                                         exit

                         else
                           stop ' Erro'
                         endif
                            

                 enddo
       
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

! Generate Grid ---------------------------------------------------------

! como o numero de linhas e colunas é igual nas lat, longs e campo geofisicos vou assumir que
! a lat e long soa no centro de uma malha malha e nao os cantos


     Size%ILB  = 0
     Size%IUB  = maxI - minI + 2
     Size%JLB  = 0
     Size%JUB  = maxJ -minJ + 2

     WorkSize%ILB = 1
     WorkSize%IUB = maxI - minI + 1
     WorkSize%JLB = 1
     WorkSize%JUB = maxJ - minJ + 1


Allocate (GridLat  (Size%ILB:Size%IUB, Size%JLB:Size%JUB))
Allocate (GridLong (Size%ILB:Size%IUB, Size%JLB:Size%JUB))

Lines  = 1
Columns =1
do i = minI , maxI + 1
 do j= minJ , maxJ + 1
   
   
   GridLat   (Lines,Columns) = (Latitude(i-1,j)+Latitude(i,j))/2
   GridLong  (Lines,Columns) = (Longitude(i,j-1)+Longitude(i,j))/2
   
   Columns = Columns + 1 
 
 enddo
  Columns = 1
  Lines   = Lines + 1
enddo

deallocate (Latitude)
deallocate (Longitude)
deallocate (AuxArray2Dreal)

   
! -------- Reads Geophysical Data ----------------------------------------

     FileFound=LookforFile (L2Filename, FileID, n_file_attributes)
    
    if(FileFound.eq..false.) then 
                write(*,*)  
                write(*,*) 'Error! File Not Found: ', GeoFilename
                stop 'ReadHDF4 - ModuleConvertModisL2 - ERR01'
    endif
	
	 allocate(Chla(Size%ILB:Size%IUB, Size%JLB:Size%JUB))
     allocate(sst (Size%ILB:Size%IUB, Size%JLB:Size%JUB))


     sds_name    = "chlor_a"
     paramIndex  = sfn2index (FileID, sds_name)
	 sds_id      = sfselect  (FileID, paramIndex)
     
	 
     istat       = sfginfo   (sds_id, sds_name, rank, dimsizes, data_type, num_attrs)
    
	 start = 0
     
     stride = 1

     !allocate to correct dimensions !ATTENTION remember that HDF changes lines with columns
      
      allocate(auxArray2DReal2(1:dimsizes(1),1:dimsizes(2)))
      allocate(auxArray2DInt(1:dimsizes(1),1:dimsizes(2)))
      allocate(Array2DChar(1:dimsizes(1),1:dimsizes(2)))

     !reads data
      istat       = sfrdata   (sds_id, start, stride, dimsizes, AuxArray2Dreal2)
 
       Lines   =1
	   Columns =1
    
	   do i=minI , maxI
          do j = minJ , maxJ
              Chla(Lines,Columns) = AuxArray2Dreal2(j,i)
              Columns = Columns + 1 
	      enddo
	          Columns = 1
              Lines   = Lines + 1
        enddo

! sst ------

     sds_name    = "sst"
     paramIndex  = sfn2index(FileID, sds_name)
     sds_id      = sfselect  (FileID, paramIndex)

     istat       = sfginfo   (sds_id, sds_name, rank, dimsizes, data_type, num_attrs)
    
     !reads data
      istat       = sfrdata   (sds_id, start, stride, dimsizes, auxArray2DInt)


      
   
 
       Lines   =1
	   Columns =1
    
	   do i=minI , maxI
          do j = minJ , maxJ
              sst(Lines,Columns) = auxArray2DInt(j,i) * 0.005
              Columns = Columns + 1 
	       enddo
	          Columns = 1
              Lines   = Lines + 1
         enddo


    call WriteGridData  (FileName       = 'chla_'//trim(adjustl(cdate))//'_'//trim(adjustl(ctime))//'.dat',                  &
                         COMENT1        = 'ModisL2 file'//trim(adjustl(L2Filename))//' Modis Geo File: '//trim(adjustl(GeoFileName)), &
                         COMENT2        = trim(adjustl(cdate))//' '//trim(adjustl(ctime))//' '//'Chl a mg/m^3',     &
                         ConnectionX    = GridLong,                                                                 &
                         ConnectionY    = GridLat,                                                                  &
                         WorkSize       = WorkSize,                                                                 &
                         CoordType      = 4,                                                                        &
                         Xorig          = gridlong(1,size%jub),                                  &
                         Yorig          = gridlat(1,1),                                   &
                         Zone           = 29,                                                                       &
                         GRID_ANGLE     = 0.,                                                                       &
                         Latitude       = gridlat(1,1),                                  &
                         Longitude      = gridlong(1,size%jub),                                  &
                         FillValue      = -99.,                                                                     &
                         Overwrite      = .true.,                                                                   &
                         GridData2D_Real= Chla,                                                                     &
                         STAT           = STAT_CALL) 
     if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleConvertModisL2 - ERR066'

	 call WriteGridData  (FileName       = 'sst_'//trim(adjustl(cdate))//'_'//trim(adjustl(ctime))//'.dat',                                                          &
                         COMENT1        = 'ModisL2 file'//trim(L2Filename)//' Modis Geo File: '//trim(GeoFileName), &
                         COMENT2        = trim(adjustl(cdate))//' '//trim(adjustl(ctime))//' '//'Chl a mg/m^3',                                                                       &
                         ConnectionX    = GridLong,                                                                 &
                         ConnectionY    = GridLat,                                                                  &
                         WorkSize       = WorkSize,                                                                 &
                         CoordType      = 4,                                                                        &
                         Xorig          = gridlong(1,size%jub),                                  &
                         Yorig          = gridlat(1,1),                                   &
                         Zone           = 29,                                                                       &
                         GRID_ANGLE     = 0.,                                                                       &
                         Latitude       = gridlat(1,1),                                  &
                         Longitude      = gridlong(1,size%jub),                                  &
                         FillValue      = -99.,                                                                     &
                         Overwrite      = .true.,                                                                   &
                         GridData2D_Real= sst,                                                                     &
                         STAT           = STAT_CALL) 
     if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleConvertModisL2 - ERR066'
 
    deallocate(Chla)
    deallocate(sst)
	deallocate(auxArray2DReal2)
    deallocate(auxArray2DInt)
    deallocate (GridLat)
    deallocate(GridLong)
    end subroutine ReadHDF4

!-------------------------------------------------------------------------

 subroutine ReadDate (Filename, AuxTime, cdate, ctime)

    !Arguments-----------------------------------------------------------
    type(T_Attribute), dimension(:), pointer     :: Attributes
	type(T_Time)                                 :: Date
	
	character(len=256), intent (IN):: Filename
	real, dimension(6),intent (OUT):: AuxTime
	character(len=256),intent (OUT):: cdate, ctime
	character(len=256)             :: chartime 
	character(len=4)               :: char_i , cmonth, cday                   
    integer(4)                     :: FileID,  n_file_attributes, file_attr_index, file_attr_status
	integer(4)                     :: read_attr_status, dayj, AuxDayj
	logical                        :: FileFound
    integer(4)                     :: sfrnatt, sfrcatt, sfgainfo, sffinfo, day, month, i
	
    !---------------------------------------------------------------------
 

	 FileFound=LookforFile (Filename, FileID, n_file_attributes)
    
	  
     
	 if(FileFound.eq..false.) then 
                write(*,*)  
                write(*,*) 'Error! File Not Found: ', Filename
                stop 'ReadDate - ModuleConvertModisL2 - ERR01a'
      endif
    
	  nullify  (Attributes)


      allocate (Attributes(0:n_file_attributes-1))

       
	!read all file attributes
    do file_attr_index = 0, n_file_attributes-1
	  
	     
        
        !get file attributes
        file_attr_status = sfgainfo(FileID, file_attr_index,                &
                                    Attributes(file_attr_index)%name,       &
                                    Attributes(file_attr_index)%att_type,   &
                                    Attributes(file_attr_index)%length)

        if(Attributes(file_attr_index)%att_type == 4)then
             
            !reads string type attributes
            read_attr_status = sfrcatt (FileID, file_attr_index, Attributes(file_attr_index)%String)
             
		end if

    enddo

      
    CharTime=trim(adjustl(attributes(20)%String))
      
      
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
      
	  

   
   do day =1, 31
    do month =1, 12
     
    AuxTime(2) =  Month
    AuxTime(3) =  Day

      call SetDate( Date,  Year = AuxTime(1), Month  = AuxTime(2), Day    = AuxTime(3), &
                                   Hour = AuxTime(4), Minute = AuxTime(5), Second = AuxTime(6)) 
    
	 call JulianDay( Date, AuxDayj)
    
	 if (AuxDayj.eq.DayJ) then 
        exit
     endif
   
    enddo
     
    if (AuxDayj.eq.DayJ) exit
    
   enddo
   
    
    ctime = CharTime(8:9)//'-'//CharTime(10:11)//'-'//CharTime(12:13)
    write(cmonth, '(i4)')int(AuxTime(2))
	write(cday, '(i4)')  int(AuxTime(3))
	cdate = CharTime(1:4)//'-'//trim(adjustl(cmonth))//'-'//trim(adjustl(cday))

   deallocate (Attributes)

   end subroutine

!-------------------------------------------------------------------------

    Function LookforFile (Filename, FileID, n_file_attributes)
    logical :: LookforFile

    !Arguments-----------------------------------------------------------
    character(len=256), intent (IN):: Filename                      
    integer(4), intent (OUT)       :: FileID,  n_file_attributes
    integer(4), parameter          :: DFACC_READ_=1, DFACC_WRITE_=2, DFACC_CREATE_=4  !File Access
    integer(4)                     :: file_info_status, n_datasets, hopen, sfstart, sffinfo
    !---------------------------------------------------------------------

    FileID      = hopen     (FileName, DFACC_READ_, 0)
    FileID      = sfstart   (FileName, DFACC_READ_)

    file_info_status = sffinfo(FileID, n_datasets, n_file_attributes)

    if (file_info_status.eq.-1) then

       LookforFile=.false.
        
    else
       
       LookforFile=.true.
       
    endif
    end function
 !-------------------------------------------------------------------------
    
    subroutine KillModisL2
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

               
        
        
    
    end subroutine KillModisL2


    !------------------------------------------------------------------------
 
end module ModuleConvertModisL2










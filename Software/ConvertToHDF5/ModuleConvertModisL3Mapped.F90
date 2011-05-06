!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : ConvertModisL3Mapped
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : July 2003
! REVISION      : Pedro Pina - v4.0
! DESCRIPTION   : Module to convert ModisL3Mapped Format files into HDF5 format
! DOWNLOAD      : To download ModisL3Mapped files go to http://oceancolor.gsfc.nasa.gov/cgi/level3.pl
!------------------------------------------------------------------------------


Module ModuleConvertModisL3Mapped

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleTime
    use ModuleGridData


    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    public  ::  StartConvertModisL3Mapped
    private ::  AllocateVariables
    private ::  Open_HDF5_OutPut_File
    private ::  Construct_FieldsList
    private ::        Construct_Field
    private ::        Add_Field
    private ::  ReadHDF4
    private ::  OutputFields
    private ::  KillModisL3MappedFormat

    
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
    
    private :: T_ModisL3Mapped
    type       T_ModisL3Mapped
        character(len=PathLength)                :: FileName
        character(len=StringLength)              :: ParameterName
        character(len=StringLength)              :: Units  
        integer                                  :: IDNumber
        type(T_Time)                             :: Date
        real, dimension(:,:),       pointer      :: Scalar
        integer, dimension(:,:,:),  pointer      :: OpenPoints3D
        type(T_ModisL3Mapped),      pointer      :: Next
        type(T_ModisL3Mapped),      pointer      :: Prev
        type(T_Attribute), dimension(:), pointer :: Attributes
    end type  T_ModisL3Mapped
    
    

   


    private :: T_ModisL3MappedFormat
    type       T_ModisL3MappedFormat
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
        real                                    :: MaxOrigLat

        
     
        character(len=PathLength)               :: OutputFileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: InputGrid
        integer                                 :: InputType
        logical                                 :: FirstFile            = .true. 
        integer                                 :: Unit
        real                                    :: MaxLong
        real                                    :: MinLong
        real                                    :: MaxLat
        real                                    :: MinLat
        integer                                 :: FieldsNumber
        real, dimension(:,:),  pointer          :: Bathymetry
        integer, dimension(:,:), pointer        :: Frequency

        type (T_Size2D)                         :: WorkSize
        type (T_Size2D)                         :: Size
        type (T_Size2D)                         :: AuxSize 

        real , dimension(:),  pointer           :: XX
        real , dimension(:),  pointer           :: YY
        type(T_ModisL3Mapped),      pointer     :: FirstField
        type(T_ModisL3Mapped),      pointer     :: LastField

                
       
    end type  T_ModisL3MappedFormat


    type(T_ModisL3MappedFormat),         pointer     :: Me

    integer, parameter :: Grid_              = 1
    integer, parameter :: Area_              = 2

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartConvertModisL3Mapped(EnterDataID, ClientNumber,  STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !External----------------------------------------------------------------
        
        
        !Local-------------------------------------------------------------------
        integer                                         :: nUsers
    
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions (ClientNumber)
        
        call OutputFields

        Call KillModisL3MappedFormat

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'ConvertModisL3Mapped - ModuleConvertModisL3Mapped - ERR01' 

        STAT = SUCCESS_

        deallocate(Me)
        nullify   (Me)



    end subroutine StartConvertModisL3Mapped

    !------------------------------------------------------------------------

    

    
    subroutine ReadOptions (ClientNumber)

        ! Arguments -----------------------------------------------------------

        integer,           intent(IN )              :: ClientNumber
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
       
        

        
      call Construct_FieldsList (ClientNumber)
       


    end subroutine ReadOptions 

    !------------------------------------------------------------------------

    subroutine AllocateVariables


       


    end subroutine AllocateVariables

    !------------------------------------------------------------------------
   
    subroutine Construct_FieldsList(ClientNumber)

        ! Arguments -----------------------------------------------------------

        integer,           intent(IN )              :: ClientNumber

        !External----------------------------------------------------------------
        integer                         :: STAT_CALL
        logical                         :: BlockFound

        !Local-------------------------------------------------------------------
        type (T_ModisL3Mapped), pointer      :: NewField
        integer                              :: iflag
        !------------------------------------------------------------------------

        call GetData(Me%InputType,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'INPUTTYPE',                        &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify input Type'
            stop 'ReadOptions - ModuleConvertModisL3Mapped - ERR05'
        end if
        
        if (Me%InputType==Grid_) then
          
                call GetData(Me%InputGrid,                                     &
                            Me%ObjEnterData, iflag,                            &
                            SearchType   = FromBlock,                          &
                            keyword      = 'INPUTGRID',                        &
                            ClientModule = 'ConvertModisL3Mapped',             &
                            STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR01'
        
                if (iflag == 0)then
                            write(*,*)'Must specify input Grid'
                            stop 'ReadOptions - ModuleConvertModisL3Mapped - ERR0345'
                end if

         elseif (Me%InputType==area_) then

                    call GetData(Me%MaxLong,                                        &
                                Me%ObjEnterData, iflag,                            &
                                SearchType   = FromBlock,                          &
                                keyword      = 'MAX_LONG',                          &
                                ClientModule = 'ConvertModisL3Mapped',             &
                                STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR02'
        
                    if (iflag == 0)then
                                write(*,*)'Must specify MaxLong'
                                stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR06'
                    end if


                    call GetData(Me%MinLong,                                        &
                                Me%ObjEnterData, iflag,                            &
                                SearchType   = FromBlock,                          &
                                keyword      = 'MIN_LONG',                         &
                                ClientModule = 'ConvertModisL3Mapped',             &
                                STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR03'
        
                    if (iflag == 0)then
                        write(*,*)'Must specify MinLong'
                        stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR07'
                    end if
        
        
                    call GetData(Me%MaxLat,                                         &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'MAX_LAT',                          &
                                 ClientModule = 'ConvertModisL3Mapped',             &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR04'
        
                    if (iflag == 0)then
                        write(*,*)'Must specify name of file to convert'
                        stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR08'
                    end if


                    call GetData(Me%MinLat,                                         &
                                 Me%ObjEnterData, iflag,                            &
                                 SearchType   = FromBlock,                          &
                                 keyword      = 'MIN_LAT',                          &
                                 ClientModule = 'ConvertModisL3Mapped',             &
                                 STAT         = STAT_CALL)        
                    if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleConvertModisL3Mapped - ERR03'
        
                    if (iflag == 0)then
                        write(*,*)'Must specify name of file to convert'
                        stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR09'
                    end if
        
                    call GetData(Me%GridFileName,                                   &
                                Me%ObjEnterData, iflag,                            &
                                SearchType   = FromBlock,                          &
                                keyword      = 'GRIDFILENAME',                     &
                                ClientModule = 'ConvertModisL3Mapped',             &
                                STAT         = STAT_CALL)        
                        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR01a'
        
                    if (iflag == 0)then
                            write(*,*)'Must specify name for output grid file'
                            stop 'ReadOptions - ModuleConvertModisL3Mapped - ERR05a'
                            end if
        else
 
                    write(*,*)'Must specify correct input type number (1-Grid, 2-Area)'
                    stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR09'
        endif
        
        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR01'
        
        if (iflag == 0)then
            write(*,*)'Must specify name of output file'
            stop 'ReadOptions - ModuleConvertModisL3Mapped - ERR05'
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
                        stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR01'

                    exit do1    !No more blocks
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR02'
            else cd1
                stop 'Construct_FieldsList - ModuleConvertModisL3Mapped - ERR03'
            end if cd1
        end do do1

        !------------------------------------------------------------------------

    end subroutine Construct_FieldsList

    !--------------------------------------------------------------------------
    !This subroutine reads all the information needed to construct a new Field       

    subroutine Construct_Field(NewField)

        !Arguments-------------------------------------------------------------
 
         type(T_ModisL3Mapped), pointer       :: NewField
  
        !Local ----------------------------------------------------------------
         integer                               :: STAT_CALL
         integer                               :: iflag


        !----------------------------------------------------------------------
             
        allocate (NewField)

        nullify(NewField%Scalar        )
        nullify(NewField%OpenPoints3D  )
        nullify(NewField%Prev,NewField%Next)

         call GetData(NewField%FileName,                                &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlockInBlock,                   &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'ConvertModisL3Mapped',             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'Construct_Field - ModuleConvertModisL3Mapped - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'Construct_Field - ModuleConvertModisL3Mapped - ERR02'
        end if
        

        call ReadHDF4 (NewField)



    end subroutine Construct_Field

    !--------------------------------------------------------------------------
    

    ! This subroutine adds a new Field to the Fields List  

    subroutine Add_Field(NewField)

        !Arguments-------------------------------------------------------------
        type(T_ModisL3Mapped),           pointer     :: NewField

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

    !Arguments-------------------------------------------------------------
        type(T_ModisL3Mapped), pointer       :: NewField
        
    !----------------------------------------------------------------------
    integer(4)                              :: hopen, sfstart, istat, index, sfselect, sfgainfo
    integer(4)                              :: read_attr_status, file_attr_index, sfrcatt, sfrnatt
    integer(4)                              :: file_info_status, sffinfo, n_datasets, n_file_attributes
    integer(4)                              :: FileID, sds_id, sfginfo
    integer(4), parameter                   :: DFACC_READ_=1, DFACC_WRITE_=2, DFACC_CREATE_=4  !File Access
    integer, dimension(3)                   :: dimsizes
    integer(4), parameter                   :: LABLEN=64
    integer(4)                              :: rank, data_type, num_attrs
    integer(4), dimension(2)                :: start, stride
    character, dimension(:,:), pointer      :: Array2DChar
    integer, dimension(:,:), pointer        :: Array2DInt

    

    character(len=256)                      :: FileName
    character(len=LABLEN)                   :: sds_name
    integer                                 :: sfrcdata
    integer                                 :: file_attr_status

    integer                                 :: i, j , dayj
    real                                    :: base, slope, intercept
    real                                    :: month, day ,Year
    real, dimension(6)                      :: AuxTime
    logical                                 :: Calc_sst    
    type(T_Time)                            :: AuxDate
    !----------------------------------------------------------------------

    
    
    FileName    = NewField%FileName; FileName    = trim(adjustl(FileName))
    FileID      = hopen     (FileName, DFACC_READ_, 0)
    FileID      = sfstart   (FileName, DFACC_READ_)

    file_info_status = sffinfo(FileID, n_datasets, n_file_attributes)

    write(*,*) 'Processing : ', FileName
    
    if (file_info_status.eq.-1) then
        write(*,*)  
        write(*,*) 'File not found. ', FileName
        stop 'ReadHDF4 - ModuleConvertModisL3Mapped - ERR01'
    endif
    !allocate at
    nullify  (NewField%Attributes)
!    nullify  (Array2D)

    allocate (NewField%Attributes(0:n_file_attributes-1))

    !read all file attributes
    do file_attr_index = 0, n_file_attributes-1
        
        !get file attributes
        file_attr_status = sfgainfo(FileID, file_attr_index,                &
                                    NewField%Attributes(file_attr_index)%name,       &
                                    NewField%Attributes(file_attr_index)%att_type,   &
                                    NewField%Attributes(file_attr_index)%length)

        if(NewField%Attributes(file_attr_index)%att_type == 4)then
            
            !reads string type attributes
            read_attr_status = sfrcatt (FileID, file_attr_index, NewField%Attributes(file_attr_index)%String)
                                        
        elseif (NewField%Attributes(file_attr_index)%att_type == 5)then
            
            !reads numerical type attributes
            read_attr_status = sfrnatt (FileID, file_attr_index, NewField%Attributes(file_attr_index)%ValueReal)

         elseif (NewField%Attributes(file_attr_index)%att_type == 22)then
            
            !reads numerical type attributes
            read_attr_status = sfrnatt (FileID, file_attr_index, NewField%Attributes(file_attr_index)%ValueInt)
        
        end if

    enddo
    
    !first data set in file
    index       = 0 

    !select data set and get information
    sds_id      = sfselect  (FileID, index)
    sds_name    = "l3m_data"
    istat       = sfginfo   (sds_id, sds_name, rank, dimsizes, data_type, num_attrs)

    
    NewField%ParameterName = trim(adjustl(NewField%Attributes(49)%string(1:5)))
   

    if (NewField%ParameterName(1:5).eq.'Sea S') then
     Calc_sst = .true.
     NewField%ParameterName='sst'
    else
     Calc_sst= .false.
    endif
   
      NewField%Units         = trim(adjustl(NewField%Attributes(51)%string))

    !allocate to correct dimensions
    allocate(Array2Dint(1:dimsizes(2),1:dimsizes(1)))
    allocate(Array2Dchar(1:dimsizes(1),1:dimsizes(2)))
    
    start       = 0
    stride      = 1

    !reads data
    istat       = sfrcdata   (sds_id, start, stride, dimsizes, Array2Dchar)

    do i=1 , dimsizes(2)
     do j = 1 , dimsizes(1)
       Array2Dint(i,j) = IChar(Array2Dchar(j,i))
     enddo
    enddo

    
      me%LatStep   = NewField%Attributes(42)%valueReal
      me%LongStep  = NewField%Attributes(43)%valueReal

      me%OrigLat   = NewField%Attributes(44)%valueReal
      me%OrigLong  = NewField%Attributes(45)%valueReal

   
      Me%MaxOrigLat   = me%OrigLat +  dimsizes(2) * me%LatStep ! Max North Coordenate

     
     
    !select information to extract

   if (me%firstFile) then
     
     if(Me%InputType==Grid_) then
       Call DefineGrid
     else 
       Call DefineArea
     endif
     
     me%firstFile=.false.
   endif

     !checks for any changes on grid (not probable but just in case ...)
     if      ((me%LatStep.ne.NewField%Attributes(42)%valueReal).or.&
            (me%LongStep.ne.NewField%Attributes(43)%valueReal).or.&
            (me%OrigLat.ne.NewField%Attributes(44)%valueReal).or.&
            (me%OrigLong.ne.NewField%Attributes(45)%valueReal)) then
            stop 'ReadHDF4 - ModuleConvertModisL3Mapped - ERR0fd'
    
     endif
     
     
     Allocate (Me%Frequency(me%Size%ILB:me%Size%IUB, me%Size%JLB:me%Size%JUB))
     Me%Frequency  = 0
    

   if (Calc_sst) then
      slope     = NewField%Attributes(54)%valueReal
      intercept = NewField%Attributes(55)%valueReal
   else 
       base      = NewField%Attributes(54)%valueReal
       slope     = NewField%Attributes(55)%valueReal
       intercept = NewField%Attributes(56)%valueReal
   endif
  

       dayj      = NewField%Attributes(27)%valueInt
       year      = NewField%Attributes(29)%valueInt
    
     
     Allocate (NewField%Scalar(me%Size%ILB:me%Size%IUB , me%Size%JLB:me%Size%JUB ))
     Allocate (NewField%OpenPoints3D (me%Size%ILB:me%Size%IUB,me%Size%JLB:me%Size%JUB,1))

     !Extract and apply scalling equation
    
     !do i = me%WorkSize%ILB, me%WorkSize%IUB
     ! do j = me%WorkSize%JLB, me%WorkSize%JUB
     !   if (Array2DInt(Me%AuxSize%IUB - i +1 ,j + Me%AuxSize%JLB - 1).eq.255) then
     !      NewField%OpenPoints3D(i,j,1) = 0
     !      NewField%Scalar(i,j)         = 0.0
     !   else
     !         Me%Frequency(i,j) = Me%Frequency(i,j) + 1
     !         NewField%OpenPoints3D(i,j,1) = 1.0
              !attention in Mohid i=1 is on the SW point and in Modis =1 is in NW point
     !         if (Calc_sst) then
     !              NewField%Scalar(i,j)  =  (Slope * Array2DInt(Me%AuxSize%IUB - i +1 ,j + Me%AuxSize%JLB - 1)) + Intercept
     !         else 
     !              NewField%Scalar(i,j)  =  Base ** ((Slope * Array2DInt(Me%AuxSize%IUB - i +1 ,j + Me%AuxSize%JLB - 1)) + Intercept)
     !         endif 
     !   endif
     ! enddo
     !enddo 

     do i = me%WorkSize%ILB, me%WorkSize%IUB
     do j = me%WorkSize%JLB, me%WorkSize%JUB
         
         if (Calc_sst) then
           NewField%Scalar(i,j)= Slope * float(Array2DInt((Me%AuxSize%IUB - i +1) ,(j + Me%AuxSize%JLB - 1))) + Intercept
         else
           NewField%Scalar(i,j)  =  Base ** ((Slope * float(Array2DInt(Me%AuxSize%IUB - i +1 ,j + Me%AuxSize%JLB - 1))) + Intercept)
         endif
          
           NewField%OpenPoints3D(i,j,1) = 1
     enddo
     enddo
    
     do i = me%WorkSize%ILB, me%WorkSize%IUB
     do j = me%WorkSize%JLB, me%WorkSize%JUB

       if (Array2DInt(Me%AuxSize%IUB - i +1 ,j + Me%AuxSize%JLB - 1).eq.255) then
         NewField%OpenPoints3D(i,j,1) = 0
         NewField%Scalar(i,j)         = 0.0
       else
         Me%Frequency(i,j) = Me%Frequency(i,j) + 1
       endif

     enddo
     enddo
     
    !do i=Me%WorkSize%ILB, Me%WorkSize%IUB
    ! do j=Me%WorkSize%JLB, Me%WorkSize%JUB 
    !    if (me%bathymetry(i,j).ne.-99) then
     !       NewField%OpenPoints3D(i,j,1) = 1
     !   else
     !       NewField%OpenPoints3D(i,j,1) = 0
     !   endif
     !enddo
    !enddo

    call JulianDayToMonthDay(1980, DayJ, AuxDate)
    call ExtractDate(AuxDate, Month = Month, Day = Day)

    AuxTime(1) = Year
    AuxTime(2) = Month
    AuxTime(3) = Day
    AuxTime(4) = 0.0
    AuxTime(5) = 0.0
    AuxTime(6) = 0.0

    call SetDate(NewField%Date, Year = AuxTime(1), Month  = AuxTime(2), Day    = AuxTime(3), &
                                   Hour = AuxTime(4), Minute = AuxTime(5), Second = AuxTime(6)) 

    deallocate(Array2Dchar)
    deallocate(Array2Dint)

    end subroutine ReadHDF4

    !-------------------------------------------------------------------------

    Subroutine DefineGrid 

    ! Arguments -------------------------------------------------------------

    integer                                          :: STAT_CALL
    real                                             :: aux1,aux2 

    !------------------------------------------------------------------------

     call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%InputGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGrid - ModuleConvertModisL3Mapped - ERR01'

     call ConstructGridData          (GridDataID       = Me%ObjBathymetry,       &
                                     HorizontalGridID = Me%ObjHorizontalGrid,    &
                                     FileName         = Me%InputGrid,            &
                                     STAT             = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_) stop 'DefineGrid - ModuleConvertModisL3Mapped - ERR05'

     call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGrid - ModuleConvertModisL3Mapped - ERR02'

     call GetGridOrigin(Me%ObjHorizontalGrid, Me%MinLong, Me%MinLat, STAT = STAT_CALL)

     call GetGridData(Me%ObjBathymetry, Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'DefineGrid - ModuleConvertModisL3Mapped - ERR03'

     call GetHorizontalGrid(Me%ObjHorizontalGrid, XX = Me%XX, YY = Me%YY, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DefineGrid - ModuleConvertModisL3Mapped - ERR04'
       
     aux1=Abs((Me%XX(2) - Me%XX(1))/me%LongStep-1)
     aux2=Abs((Me%YY(2) - Me%YY(1))/me%LatStep-1)

     ! tests grid consistency
        if ((aux1.gt.0.00001).or.(aux2.gt.0.00001)) then
               write(*,*)  
               write(*,*) 'Grid inconsistency ' 
               stop 'ReadHDF4 - ModuleConvertModisL3Mapped - ERR01'

        endif   
     ! 

       Me%MaxLat  = Me%MinLat + Me%YY(me%Size%IUB)
       Me%MaxLong = Me%MinLong + Me%XX(me%Size%JUB)

       !i=1 no topo superior esquerdo

       Me%AuxSize%ILB  = (Me%maxOrigLat - Me%MaxLat) / Me%LatStep
       Me%AuxSize%IUB  = (Me%maxOrigLat - Me%MinLat) / Me%LatStep

       Me%AuxSize%JUB   = int(Abs(Me%MaxLong - me%OrigLong) / me%LongStep)
       Me%AuxSize%JLB   = int(Abs(Me%MinLong - me%OrigLong) / me%LongStep)

    end subroutine DefineGrid

!  -------------------------------------------------------------------------
  
  Subroutine DefineArea 

    ! Arguments -------------------------------------------------------------

    !Local-------------------------------------------------------------------

   
    integer                      ::  nlines, ncolumns, i, j   
    integer                      ::  aux, STAT_CALL               

    !------------------------------------------------------------------------
    
   
      !i=1 no topo superior esquerdo
    
       Me%AuxSize%IUB  = (Me%maxOrigLat - Me%MaxLat) / Me%LatStep
    
       Me%AuxSize%ILB  = (Me%maxOrigLat - Me%MinLat) / Me%LatStep
     

     if (Me%AuxSize%IUB.lt.Me%AuxSize%ILB) then
       aux = Me%AuxSize%IUB
       Me%AuxSize%IUB = Me%AuxSize%ILB
       Me%AuxSize%ILB = aux
      endif

     Me%AuxSize%JUB   = int(Abs(Me%MaxLong - me%OrigLong) / me%LongStep)
     Me%AuxSize%JLB   = int(Abs(Me%MinLong - me%OrigLong) / me%LongStep)

     if (Me%AuxSize%JUB.lt.Me%AuxSize%JLB) then
       aux = Me%AuxSize%JUB
       Me%AuxSize%JUB = Me%AuxSize%JLB
       Me%AuxSize%JLB = aux
      endif

     nlines    = Abs (Me%AuxSize%IUB - Me%AuxSize%ILB)
     ncolumns  = Abs (Me%AuxSize%JUB - Me%AuxSize%JLB)

     Me%Size%ILB  = 0
     Me%Size%IUB  = nlines + 1
     Me%Size%JLB  = 0
     Me%Size%JUB  = ncolumns + 1

     Me%WorkSize%ILB = 1
     Me%WorkSize%IUB = nlines
     Me%WorkSize%JLB = 1
     Me%WorkSize%JUB = ncolumns

     !make grid

    Allocate (me%XX(me%Size%JLB:me%Size%JUB))
    Allocate (me%YY(me%Size%ILB:me%Size%IUB))
    

    me%YY(1) = 0.
    do i = 2, me%Size%IUB
      me%YY(i) =   i * Me%LatStep
    enddo

    me%XX(1) = 0.
    do j = 2, me%Size%JUB
      me%XX(j) =   j * Me%LongStep
    enddo

    Allocate (Me%Bathymetry(me%Size%ILB:me%Size%IUB, me%Size%JLB:me%Size%JUB))
    
          Me%Bathymetry = -5.0
          
    
     
    call WriteGridData  (FileName       = Me%GridFileName,                              &
                         COMENT1        = 'Grid Data created from ModisL3Mapped file',  &
                         COMENT2        = trim(Me%OutputFilename),                      &
                         XX             = Me%XX,                                        &
                         YY             = Me%YY,                                        &
                         WorkSize       = me%WorkSize,                                  &
                         CoordType      = 4,                                            &
                         Xorig          = Me%OrigLong + (Me%AuxSize%JLB-2) * Me%LongStep,      &
                         Yorig          = Me%maxOrigLat  - (Me%AuxSize%IUB+1) * Me%LatStep,           &
                         Zone           = 29,                                           &
                         GRID_ANGLE     = 0.,                                           &
                         Latitude       = Me%maxOrigLat  - (Me%AuxSize%IUB +1) * Me%LatStep,           &
                         Longitude      = Me%OrigLong + (Me%AuxSize%JLB -2) * Me%LongStep,      &
                         FillValue      = -99.,                                         &
                         Overwrite      = .true.,                                       &
                         GridData2D_Real= Me%Bathymetry,                                &
                         STAT           = STAT_CALL) 
      if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleConvertModisL3Mapped - ERR066'
 
      call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleConvertModisL3Mapped - ERR0iu2'
   
   
   
    end subroutine DefineArea

    !-------------------------------------------------------------------------
 
    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertModisL3Mapped - ERR01'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR02'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",       &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertModisL3Mapped - ERR03'            

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertModisL3Mapped - ERR04'

        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR041'
        
        call HDF5WriteData   (Me%ObjHDF5, "/Frequency", "Frequency", "-",       &
                              Array2D =  Me%Frequency, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertModisL3Mapped - ERR05'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleConvertModisL3Mapped - ERR07'

    end subroutine Open_HDF5_OutPut_File
 
   !------------------------------------------------------------------------


    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, OutputNumber
        type (T_ModisL3Mapped), pointer                 :: NewField
        type(T_Time)                                    :: CurrentDate

        !Begin-----------------------------------------------------------------
        
        call Open_HDF5_OutPut_File

        OutputNumber = 1
        CurrentDate  = Me%FirstField%Date
        
        write(*,*)'Writing field number: ', OutputNumber

        call ExtractDate   (CurrentDate,                                &
                            AuxTime(1), AuxTime(2), AuxTime(3),         &
                            AuxTime(4), AuxTime(5), AuxTime(6))
        TimePtr => AuxTime

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR01'


        call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                             "Time", "YYYY/MM/DD HH:MM:SS",             &
                             Array1D = TimePtr,                         &
                             OutputNumber = OutPutNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR02'

        NewField => Me%FirstField

        do while(associated(NewField))

            if(NewField%Date .gt. CurrentDate)then

                CurrentDate = NewField%Date

                OutputNumber = OutputNumber + 1

                write(*,*)'Writing field number: ', OutputNumber

                call ExtractDate   (CurrentDate,                                &
                                    AuxTime(1), AuxTime(2), AuxTime(3),         &
                                    AuxTime(4), AuxTime(5), AuxTime(6))
                TimePtr => AuxTime
            
                call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR03'
                

                call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                     "Time", "YYYY/MM/DD HH:MM:SS",             &
                                     Array1D = TimePtr,                         &
                                     OutputNumber = OutPutNumber, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR04'

            end if


            if(NewField%Date .eq. CurrentDate)then

                call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                                     Me%WorkSize%JLB, Me%WorkSize%JUB, 1 , 1, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR05'

                call HDF5WriteData  (Me%ObjHDF5, "/Grid/OpenPoints", "OpenPoints",       &
                                     "-", Array3D = NewField%OpenPoints3D,               &
                                      OutputNumber = OutPutNumber, STAT = STAT_CALL)

                call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                                     Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR05'
                
                call HDF5WriteData(Me%ObjHDF5,                                      &
                                   "/Results",                                      &
                                   NewField%ParameterName,                          &
                                   NewField%Units,                                  &
                                   Array2D      = NewField%Scalar,                  &
                                   OutputNumber = OutPutNumber,                     &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR06'
                
                
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR06'



            end if

            NewField => NewField%Next

        end do

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleConvertModisL3Mapped - ERR07'

    end subroutine OutputFields    
    !---------------------------------------------------------------------------

    subroutine KillModisL3MappedFormat
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        if (me%inputType==Grid_) then
          call UngetGridData(Me%ObjBathymetry, Me%Bathymetry, STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleInterpolateGrids - ERR02'
        endif
        
        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillModisL3MappedFormat - ModuleConvertModisL3Mapped - ERR01'
        
        
        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillModisL3MappedFormat - ModuleConvertModisL3Mapped - ERR02'

        deallocate (Me%Bathymetry)
        deallocate (Me%Frequency)
        deallocate (Me%XX)
        deallocate (Me%YY)
    
    end subroutine KillModisL3MappedFormat


    !------------------------------------------------------------------------
 
end module ModuleConvertModisL3Mapped










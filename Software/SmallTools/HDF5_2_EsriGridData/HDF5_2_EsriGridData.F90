!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Small Tools
! MODULE        : HDF5ToASCIIandBIN
! URL           : http://www.mohid.com
! AFFILIATION   : Hidromod Lda
! DATE          : February 2009
! REVISION      : Paulo Leitao - v1.0
! DESCRIPTION   : Module to convert matrixes from HDF5 to ESRI GRID DATA
!
!------------------------------------------------------------------------------


Program HDF5_2_EsriGridData

    use ModuleGlobalData
    use ModuleTime
    use ModuleFunctions
    use ModuleHDF5
    use ModuleEnterData
    use ModuleHorizontalGrid

    implicit none



    !Subroutines---------------------------------------------------------------

    !Constructor

    !Parameters----------------------------------------------------------------
    integer, parameter  :: HDF5_2_ESRI = 1, ESRI_2_HDF5 = 2

    

    type       T_HDF5_2_EsriGridData
        integer                                             :: ObjEnterData         = 0
        integer                                             :: ClientNumber
        integer                                             :: ObjHDF5              = 0
        integer                                             :: ObjHorizontalGrid    = 0
        character(len=PathLength), dimension(:), pointer    :: FieldName
        character(len=PathLength), dimension(:), pointer    :: OutputESRI
        
        integer,                   dimension(:), pointer    :: PropESRI
        type(T_Time),              dimension(:), pointer    :: TimeESRI
        
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: InPutFileName
        character(len=PathLength)                           :: OutPutFileName

        integer                                             :: FieldNumber         = FillValueInt
        integer                                             :: PropNumber          = FillValueInt        
        integer                                             :: InstantNumber       = FillValueInt                
        
        integer                                             :: Conversion          = FillValueInt

        type(T_Size2D)                                      :: WorkSize, Size
        real                                                :: FillValue
        integer, dimension(4)                               :: Window
        logical                                             :: WindowON
        integer                                             :: Layer, KUB
        character(30)                                       :: DX
        character(30), dimension(2)                         :: Origin        
    end type  T_HDF5_2_EsriGridData

    type(T_HDF5_2_EsriGridData), pointer              :: Me


    !--------------------------------------------------------------------------
    
        nullify (Me)
        allocate(Me) 
        
        call ConstructHDF5_2_EsriGridData
        
        if      (Me%Conversion == HDF5_2_ESRI) then
        
            call ModifyHDF5_2_EsriGridData
            
        else if (Me%Conversion == ESRI_2_HDF5) then
        
            call ModifyEsriGridData_2_HDF5    
            
        endif
        
        !call KillHDF5_2_EsriGridData
    
    
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructHDF5_2_EsriGridData

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: STAT_CALL
        !------------------------------------------------------------------------

        nullify (Me)
        allocate(Me)

        call ConstructEnterData (Me%ObjEnterData, "HDF5_2_EsriGridData.dat", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)  stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR10'

        call ReadGlobalOptions
        
        if      (Me%Conversion == HDF5_2_ESRI) then
        
            call ReadHDF5_2_ESRI_Options
            
        else if (Me%Conversion == ESRI_2_HDF5) then
        
            call ReadESRI_2_HDF5_Options
            
        endif
                
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR20'

    end subroutine ConstructHDF5_2_EsriGridData
 
    !--------------------------------------------------------------------------

    subroutine ReadHDF5Fields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL, Line, FirstLine, LastLine
        character(len=StringLength)                 :: AuxChar
        logical                                     :: BlockFound
        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBuffer (Me%ObjEnterData,                                   &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = '<begin_field>',                &
                                    block_end         = '<end_field>',                  &
                                    BlockFound        = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :       if (.not. BlockFound) then                                                  
                stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR10'
            end if cd2
        else cd1
            stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR20'
        end if cd1


        Me%FieldNumber = LastLine - FirstLine - 1

        allocate(Me%FieldName  (Me%FieldNumber))  
        allocate(Me%OutputESRI (Me%FieldNumber))  
  

d1:     do l= 1, Me%FieldNumber

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR30'


            Me%FieldName(l) = AuxChar

         enddo d1

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR40'


        call ExtractBlockFromBuffer  (Me%ObjEnterData,                                  &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = '<begin_esri>',                 &
                                    block_end         = '<end_esri>',                   &
                                    BlockFound        = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)
 

d2:     do l= 1, Me%FieldNumber

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR50'


            Me%OutputESRI(l) = AuxChar

         enddo d2

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR60'

    end subroutine ReadHDF5Fields 

    !--------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag
        integer, dimension(4)                       :: Iaux

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'INPUT_GRID_FILENAME',                              &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR20'

        Iaux (:) = -99

        call GetData(Me%WindowON,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'WINDOW_ON',                                        &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR30'
        
        if (Me%WindowON) then        
                    call GetData(Me%WindowON,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'WINDOW',                                           &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR30'

        endif
        
        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.,                                               &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR40'
        


        call GetData(Me%Conversion,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'CONVERSION_TYPE',                                  &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     default      = HDF5_2_ESRI,                                        &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR50'
        if (Me%Conversion /= HDF5_2_ESRI .and. Me%Conversion /= ESRI_2_HDF5)            &
                                   stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR60'

 
    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------
    
    subroutine ReadHDF5_2_ESRI_Options

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%InPutFileName,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'INPUTFILENAME',                                    &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR10'
        if (iflag     == 0)        stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR20'

        call GetData(Me%Layer,                                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'LAYER',                                            &
                     default      = 1,                                                  &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR30'
        
        call GetData(Me%KUB,                                                            &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'KUB',                                              &
                     default      = 1,                                                  &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR40'

        call GetData(Me%DX,                                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DX',                                               &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR50'

        if (iflag == 0) then
            stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR60'
        endif

        call GetData(Me%Origin,                                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'ORIGIN',                                           &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR70'

        if (iflag == 0) then
            stop 'ReadHDF5_2_ESRI_Options - HDF5_2_EsriGridData - ERR80'
        endif

    end subroutine ReadHDF5_2_ESRI_Options

    !--------------------------------------------------------------------------
    
    subroutine ReadESRI_2_HDF5_Options

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag

        !Begin-----------------------------------------------------------------

        call GetData(Me%OutPutFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadESRI_2_HDF5_Options - HDF5_2_EsriGridData - ERR10'
        if (iflag     == 0)        stop 'ReadESRI_2_HDF5_Options - HDF5_2_EsriGridData - ERR20'

    end subroutine ReadESRI_2_HDF5_Options

    !--------------------------------------------------------------------------


    subroutine ReadESRIFields

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL, Line, FirstLine, LastLine
        character(len=StringLength)                 :: AuxChar
        real, dimension(6)                          :: AuxVector
        logical                                     :: BlockFound
       
        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBuffer (Me%ObjEnterData,                                   &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = '<begin_properties>',           &
                                    block_end         = '<end_properties>',             &
                                    BlockFound        = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :       if (.not. BlockFound) then                                                  
                stop 'ReadESRIFields - ModuleHDF5ToASCIIandBIN - ERR10'
            end if cd2
        else cd1
            stop 'ReadESRIFields - ModuleHDF5ToASCIIandBIN - ERR20'
        end if cd1


        Me%PropNumber = LastLine - FirstLine - 1

        allocate(Me%PropESRI  (Me%PropNumber))  

d1:     do l= 1, Me%PropNumber

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadESRIFields - HDF5_2_EsriGridData - ERR30'

            if (.not. CheckPropertyName (AuxChar, Me%PropESRI(l))) then
            
                stop 'ReadESRIFields - HDF5_2_EsriGridData - ERR40'
            
            endif

         enddo d1

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5Fields - HDF5_2_EsriGridData - ERR50'


        call ExtractBlockFromBuffer  (Me%ObjEnterData,                                  &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = '<begin_instants>',             &
                                    block_end         = '<end_instants>',               &
                                    BlockFound        = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)


cd3 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd4 :       if (.not. BlockFound) then                                                  
                stop 'ReadESRIFields - ModuleHDF5ToASCIIandBIN - ERR60'
            end if cd4
        else   cd3
            stop 'ReadESRIFields - ModuleHDF5ToASCIIandBIN - ERR70'
        end if cd3

        Me%InstantNumber = LastLine - FirstLine - 1

        allocate(Me%TimeESRI  (Me%InstantNumber)) 
 
d2:     do l= 1, Me%InstantNumber

            line = FirstLine + l

            call GetData(AuxVector, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadESRIFields - HDF5_2_EsriGridData - ERR80'
            
            call SetDate (Me%TimeESRI(l), Year  = AuxVector(1), Month = AuxVector(2),   &
                                          Day   = AuxVector(3), Hour  = AuxVector(4),   &
                                          Minute= AuxVector(5), Second= AuxVector(6))

         enddo d2

        call RewindBuffer(Me%ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadESRIFields - HDF5_2_EsriGridData - ERR90'

    end subroutine ReadESRIFields 

    !--------------------------------------------------------------------------

    subroutine ModifyHDF5_2_EsriGridData

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D                      
        character(len=50000)                            :: Line
        character(len=PathLength)                       :: VGroup, Field, AuxChar
        integer                                         :: l, i, j, k, STAT_CALL, Unit
        integer                                         :: ILB, IUB, JLB, JUB, HDF5_READ, a, ba, bt
        logical                                         :: Found2Blanks
        
        !------------------------------------------------------------------------

        call ReadHDF5Fields

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%InPutFileName, HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR30'

        allocate(Aux3D     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))  

        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                           Me%WorkSize%JLB,Me%WorkSize%JUB, 1, Me%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR40'


        if (Me%WindowON) then

            ILB = Me%Window(1); IUB = Me%Window(2);        
            JLB = Me%Window(3); JUB = Me%Window(4);                    

        else

            ILB = Me%WorkSize%ILB; IUB = Me%WorkSize%IUB;        
            JLB = Me%WorkSize%JLB; JUB = Me%WorkSize%JUB;                    
        
        endif
        
        k        = Me%Layer


d11:    do l = 1, Me%FieldNumber

                i = scan(Me%FieldName(l),'/',Back = .true.)

                AuxChar = trim(Me%FieldName(l))
                
                VGroup  = AuxChar(1:i-1)
                
                j = len_trim(Me%FieldName(l))
                
                Field  = AuxChar(i+1:j)
                
                call HDF5ReadData(Me%ObjHDF5,                                           &
                                   trim(Vgroup),                                        &
                                   trim(Field),                                         &
                                   Array3D      = Aux3D,                                &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR50'
                
                call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'OutputMohidBin - ModuleTecnoceanAscii - ERR10'

                open(Unit   = Unit,                                                     &
                     File   = trim(Me%OutputESRI(l)),                                   &
                     Form   = 'FORMATTED',                                              &
                     STATUS = 'UNKNOWN',                                                &
                     Action = 'WRITE',                                                  &
                     IOSTAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'OutputMohidBin - ModuleTecnoceanAscii - ERR20'
                

                write(Unit,'(A14,I6)'   ) 'ncols         ', JUB - JLB + 1
                write(Unit,'(A14,I6)'   ) 'nrows         ', IUB - ILB + 1
                write(Unit,'(A14,A)')     'xllcorner     ', trim(Me%Origin(1))
                write(Unit,'(A14,A)')     'yllcorner     ', trim(Me%Origin(2))
                write(Unit,'(A14,A)')     'cellsize      ', trim(Me%DX)
                write(Unit,'(A14,f12.6)') 'nodata_value  ', Me%FillValue
               
                do i = IUB, ILB, -1
                    do j = JLB, JUB
                        if (abs(Aux3D(i,j,k)) > abs(Me%FillValue)) Aux3D(i,j,k) = Me%FillValue
                    enddo
                    write(Line,'(2000(f12.6,1x))') Aux3D(i,JLB:JUB,k)
                    Line = adjustl(Line)
                    Found2Blanks = .true.

                    bt = len_Trim(Line)
                    ba = bt
                    do a = 1, bt-1
                        do while (Line(a:a+1)=='  ') 
                            Line(a:ba)=Line(a+1:ba+1)
                            ba = ba - 1
                            if (ba == a) exit
                        enddo
                        if (ba == a) exit
                    enddo
                    write(Unit,'(A)') trim(Line)
                enddo
                
                call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'OutputMohidBin - ModuleTecnoceanAscii - ERR10'

        enddo d11
        

        deallocate(Aux3D)  

    end subroutine ModifyHDF5_2_EsriGridData
 
    !--------------------------------------------------------------------------    

    subroutine ModifyEsriGridData_2_HDF5

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real,    dimension(:,:,:), pointer    :: Aux3D                      
        integer, dimension(:,:,:), pointer    :: WaterPoints    
        real,    dimension(:,:  ), pointer    :: Bathymetry
        character(len=PathLength)             :: Field, AuxChar, FileName
        character(len=StringLength)           :: PropName
        integer                               :: l, i, j, k, STAT_CALL, Unit, it, ncols, nrows
        integer                               :: ILB, IUB, JLB, JUB, HDF5_CREATE, istart, iend
        real,    dimension(6    ), target     :: AuxTime
        real,    dimension(:    ), pointer    :: TimePtr
        real                                  :: Xorig, Yorig, xllcorner, yllcorner, nodata_value     
        logical                               :: FirstTime = .true.                                           
        
        !------------------------------------------------------------------------

        call ReadESRIFields

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        Me%ObjHDF5 = 0

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutPutFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR10'
        
        call GetGridOrigin(Me%ObjHorizontalGrid, Xorig, Yorig, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR20'        
        

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR30'
        
        Me%KUB = 1
        
        allocate(Aux3D      (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))  
        allocate(WaterPoints(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))
        allocate(Bathymetry (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB          ))



!        if (Me%WindowON) then

!            ILB = Me%Window(1); IUB = Me%Window(2);        
!            JLB = Me%Window(3); JUB = Me%Window(4);                    

!        else

            ILB = Me%WorkSize%ILB; IUB = Me%WorkSize%IUB;        
            JLB = Me%WorkSize%JLB; JUB = Me%WorkSize%JUB;                    
        
!        endif

        Me%Layer = 1
        k        = Me%Layer

        call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR40'

        do it = 1, Me%InstantNumber
        
            !Writes current time
            call ExtractDate   (Me%TimeESRI(it), AuxTime(1), AuxTime(2), AuxTime(3),  &
                                AuxTime(4), AuxTime(5), AuxTime(6))
            TimePtr => AuxTime

            call HDF5WriteData  (Me%ObjHDF5, "/Time", "Time", "YYYY/MM/DD HH:MM:SS", &
                                 Array1D = TimePtr, OutputNumber = it, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR50'        
        
        enddo
        
        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                           Me%WorkSize%JLB,Me%WorkSize%JUB, 1, Me%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR60'
      

do1:    do it = 1, Me%InstantNumber
do2:    do l  = 1, Me%PropNumber

            call ExtractDate   (Me%TimeESRI(it), AuxTime(1), AuxTime(2), AuxTime(3),  &
                                AuxTime(4), AuxTime(5), AuxTime(6))
                                
            AuxChar(:)        = " "
            
            write(AuxChar(1:4  ),'(I4)') int(AuxTime(1))
            write(AuxChar(5:5  ),'(A1)') "_"
            write(AuxChar(6:7  ),'(I2)') int(AuxTime(2))
            if   (AuxChar(6:6)==" ") AuxChar(6:6) = "0"
            write(AuxChar(8:8  ),'(A1)') "_"
            write(AuxChar(9:10 ),'(I2)') int(AuxTime(3))
            if   (AuxChar(9:9)==" ") AuxChar(9:9) = "0"
            write(AuxChar(11:11),'(A1)') "_"
            write(AuxChar(12:13),'(I2)') int(AuxTime(4))
            if   (AuxChar(12:12)==" ") AuxChar(12:12) = "0"
            write(AuxChar(14:14),'(A1)') "_"
            write(AuxChar(15:16),'(I2)') int(AuxTime(5))
            if   (AuxChar(15:15)==" ") AuxChar(15:15) = "0"
            write(AuxChar(17:17),'(A1)') "_"
            write(AuxChar(18:19),'(I2)') int(AuxTime(6))
            if   (AuxChar(18:18)==" ") AuxChar(18:18) = "0"
            write(AuxChar(20:20),'(A1)') "_"
            
            PropName = GetPropertyName(Me%PropESRI(l))

            FileName = trim(AuxChar)//trim(PropName)//".dat"

            open(Unit   = Unit,                                                     &
                 File   = trim(FileName),                                           &
                 Form   = 'FORMATTED',                                              &
                 STATUS = 'UNKNOWN',                                                &
                 Action = 'READ',                                                   &
                 IOSTAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR70'
            

            read(Unit,'(A100)') AuxChar

            istart = scan(AuxChar, " " )
            iend   = Len_Trim(AuxChar) 
            
            read(AuxChar(istart:iend),*) ncols 
            
            if (ncols /= JUB - JLB + 1) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR80'
            
            read(Unit,'(A100)') AuxChar

            istart = scan(AuxChar, " " )
            iend   = Len_Trim(AuxChar) 
            
            read(AuxChar(istart:iend),*) nrows 
            
            if (nrows /= IUB - ILB + 1) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR90'

            read(Unit,'(A100)') AuxChar

            istart = scan(AuxChar, " " )
            iend   = Len_Trim(AuxChar) 
            
            read(AuxChar(istart:iend),*) xllcorner 
            
            if (xllcorner /= Xorig) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR100'
            
            read(Unit,'(A100)') AuxChar
            
            istart = scan(AuxChar, " " )
            iend   = Len_Trim(AuxChar) 
            
            read(AuxChar(istart:iend),*) yllcorner 
            
            if (yllcorner /= Yorig) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR110'          
              
            read(Unit,'(A100)') AuxChar

            read(Unit,'(A100)') AuxChar

            istart = scan(AuxChar, " " )
            iend   = Len_Trim(AuxChar) 
            
            read(AuxChar(istart:iend),*) nodata_value
                        
            do i = IUB, ILB, -1
                read(Unit,*) (Aux3D(i,j,k),j=JLB,JUB)
            enddo
            
            
            call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR120'
            
            call HDF5WriteData(Me%ObjHDF5,                                              &
                               GroupName    = "/Results/"//PropName,                    &
                               Name         = PropName,                                 &
                               Units        = '-',                                      &       
                               Array3D      = Aux3D,                                    &
                               OutputNumber = it,                                       &
                               STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR130'
            
            call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR140'
            
            
            if (FirstTime) then
                WaterPoints(:,:,:) = 0
                
                do j = JLB, JUB
                do i = ILB, IUB
                    if (Aux3D(i,j,k) /= nodata_value) WaterPoints(i,j,k) = 1
                enddo
                enddo
                
                Bathymetry(:,:) = WaterPoints(:,:,k)

                call HDF5WriteData(Me%ObjHDF5,                                          &
                                   GroupName    = "/Grid",                              &
                                   Name         = "WaterPoints3D",                      &
                                   Units        = '-',                                  &       
                                   Array3D      = WaterPoints,                          &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR150'

                deallocate(WaterPoints)



                call HDF5WriteData(Me%ObjHDF5,                                          &
                                   GroupName    = "/Grid",                              &
                                   Name         = "Bathymetry",                         &
                                   Units        = '-',                                  &       
                                   Array2D      = Bathymetry,                           &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR160'

                deallocate(Bathymetry )        
                
                FirstTime = .false.
            endif            
                
        enddo do2
        enddo do1
        

        deallocate(Aux3D)  

    end subroutine ModifyEsriGridData_2_HDF5
 
    !--------------------------------------------------------------------------
    End Program HDF5_2_EsriGridData
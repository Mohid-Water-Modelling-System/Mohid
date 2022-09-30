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
    use ModuleGridData
    use ModuleDrawing 

    implicit none



    !Subroutines---------------------------------------------------------------

    !Constructor

    !Parameters----------------------------------------------------------------
    integer, parameter  :: HDF5_2_ESRI = 1, ESRI_2_HDF5 = 2
    integer, parameter  :: InstField_  = 1, AverageField_ = 2, SumField_ = 3

    

    type       T_HDF5_2_EsriGridData
        integer                                             :: ObjEnterData         = 0
        integer                                             :: ClientNumber
        integer                                             :: ObjHDF5              = 0
        integer                                             :: ObjHorizontalGrid    = 0
        integer                                             :: ObjHorizontalGridOut = 0        
        character(len=PathLength), dimension(:), pointer    :: FieldName
        character(len=PathLength), dimension(:), pointer    :: OutputESRI
        integer,               dimension(:,:,:), pointer    :: InPutMap
        integer,               dimension(:,:  ), pointer    :: OutPutMap     
        logical                                             :: InterpolOut
        
        character(len=PathLength)                           :: LandMaskFile
        type (T_Polygon), pointer                           :: LandMask
        logical                                             :: LandMaskON
        
        
        integer,                   dimension(:), pointer    :: PropESRI
        type(T_Time),              dimension(:), pointer    :: TimeESRI
        
        character(len=PathLength)                           :: GridFileName
        character(len=PathLength)                           :: GridFileNameOut        
        character(len=PathLength)                           :: InPutFileName
        character(len=PathLength)                           :: OutPutFileName

        integer                                             :: FieldNumber         = FillValueInt
        integer                                             :: PropNumber          = FillValueInt        
        integer                                             :: InstantNumber       = FillValueInt                
        
        integer                                             :: Conversion          = FillValueInt

        type(T_Size2D)                                      :: WorkSize, Size
        type(T_Size2D)                                      :: WorkSizeOut, SizeOut
        real                                                :: FillValue
        real                                                :: MinValue
        real                                                :: MaxValue
        real                                                :: MultiplyFactor
        real                                                :: AddFactor        
        integer, dimension(4)                               :: Window
        logical                                             :: WindowON
        integer                                             :: Layer, KUB
        character(30)                                       :: DX
        character(30), dimension(2)                         :: Origin  
        logical                                             :: Dim2D      
        logical                                             :: NoNegativeValues    
        logical                                             :: TransferToOutGrid
        logical                                             :: ExportXYZ    
        logical                                             :: ExportXY_Vector
        logical                                             :: VectorON
        character(len=StringLength)                         :: Vector_X, Vector_Y
        
        type(T_Time)                                        :: StartTime, EndTime
        integer                                             :: ComputeOption
        
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
        
        if (Me%TransferToOutGrid) then
            call ConstructHorizontalGrid(Me%ObjHorizontalGridOut, Me%GridFileNameOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR30'

            call GetHorizontalGridSize(Me%ObjHorizontalGridOut, Me%SizeOut, Me%WorkSizeOut, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR40'

            !call ConstructFatherGridLocation(Me%ObjHorizontalGridOut, Me%ObjHorizontalGrid, &
            !                                 OkCross = .false., OkZ = .true.,               &
            !                                 OkU = .false., OkV = .false., STAT = STAT_CALL)
            !if (STAT_CALL /= SUCCESS_) stop 'ConstructHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR50'
            
        endif            
        

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
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR20'

        Iaux (:) = -99

        call GetData(Me%WindowON,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'WINDOW_ON',                                        &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR30'
        
        if (Me%WindowON) then        
                    call GetData(Me%WindowON,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'WINDOW',                                           &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR30'

        endif
        
        call GetData(Me%FillValue,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.,                                               &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR40'

        call GetData(Me%Conversion,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'CONVERSION_TYPE',                                  &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     default      = HDF5_2_ESRI,                                        &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR50'
        if (Me%Conversion /= HDF5_2_ESRI .and. Me%Conversion /= ESRI_2_HDF5)            &
                                   stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR60'


        call GetData(Me%TransferToOutGrid,                                              &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'TRANSFER_TO_OUT_GRID',                             &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR70'

        if (Me%TransferToOutGrid) then

            call GetData(Me%GridFileNameOut,                                                &
                         Me%ObjEnterData, iflag,                                            &
                         SearchType   = FromFile,                                           &
                         keyword      = 'OUTPUT_GRID_FILENAME',                             &
                         ClientModule = 'HDF5_2_EsriGridData',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR80'
            if (iflag     == 0)        stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR90'
 

            call GetData(Me%InterpolOut,                                                &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'INTERPOL_OUTPUT',                              &
                         default      = .false.,                                        &  
                         ClientModule = 'HDF5_2_EsriGridData',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR100'
            
            call GetData(Me%LandMaskFile,                                               &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'LAND_MASK_FILE',                               &
                         ClientModule = 'HDF5_2_EsriGridData',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR110'
                
            if (iflag == 0) then
                Me%LandMaskON  = .false.
            else
                Me%LandMaskON  = .true.
                call New(Polygons = Me%LandMask, PolygonsFileName = Me%LandMaskFile)
            endif                
 
        endif
        
        call GetData(Me%ExportXYZ,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'EXPORT_XYZ',                                       &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR120'
        
 
        call GetData(Me%VectorON,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'VECTOR_ON',                                        &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR130'
        
        if (Me%VectorON) then

            call GetData(Me%Vector_X,                                                   &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'VECTOR_X',                                     &
                         ClientModule = 'HDF5_2_EsriGridData',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR140'        
            if (iflag     == 0       ) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR150'        

            call GetData(Me%Vector_Y,                                                   &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromFile,                                       &
                         keyword      = 'VECTOR_Y',                                     &
                         ClientModule = 'HDF5_2_EsriGridData',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR160'   
            if (iflag     == 0       ) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR170'        
            
        endif
        
        call GetData(Me%ComputeOption,                                                  &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'COMPUTE_OPTION',                                   &
                     default      = InstField_,                                         &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR180'    
        
        
        if (Me%ComputeOption > InstField_) then
            
            !Reads Begin Time
            call GetData(value          = Me%StartTime,                                 &
                         EnterDataID    = Me%ObjEnterData,                              &
                         flag           = iflag,                                        & 
                         keyword        = 'START',                                      &
                         SearchType     = FromFile,                                     &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR190'

            !Reads End Time
            call GetData(value          = Me%EndTime,                                   &
                         EnterDataID    = Me%ObjEnterData,                              &
                         flag           = iflag,                                        & 
                         keyword        = 'END',                                        &
                         SearchType     = FromFile,                                     &
                         STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR200'    
        
        endif          
        
        
        call GetData(Me%MinValue,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MIN_VALUE',                                        &
                     default      = FillValueReal,                                      &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR210'
        
        call GetData(Me%MaxValue,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MAX_VALUE',                                        &
                     default      = -FillValueReal,                                     &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR220'
        
        call GetData(Me%AddFactor,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'ADD_FACTOR',                                       &
                     default      = 0.,                                                 &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR230'
        
        call GetData(Me%MultiplyFactor,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MULTIPLY_FACTOR',                                  &
                     default      = 1.,                                                 &
                     ClientModule = 'HDF5_2_EsriGridData',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - HDF5_2_EsriGridData - ERR240'
        
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


        call GetData(Me%Dim2D,                                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DIM_2D',                                           &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     Default      = .true.,                                             &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadESRI_2_HDF5_Options - HDF5_2_EsriGridData - ERR30'
        

        call GetData(Me%NoNegativeValues,                                               &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NO_NEGATVIE_VALUES',                               &
                     ClientModule = 'HDF5ToASCIIandBIN',                                &
                     Default      = .false.,                                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadESRI_2_HDF5_Options - HDF5_2_EsriGridData - ERR40'

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
        real, dimension(:,:,:), pointer                 :: Aux3D, Aux3DOut
        real, dimension(:,:), pointer                   :: Aux2DOut, Aux2D
        character(len=50000)                            :: Line
        character(len=PathLength)                       :: VGroup, Field, AuxChar
        integer                                         :: l, i, j, k, STAT_CALL, Unit
        integer                                         :: ILB, IUB, JLB, JUB, HDF5_READ, a, ba, bt
        integer                                         :: ILBout, IUBout, JLBout, JUBout, iout, jout
        logical                                         :: Found2Blanks
        real,  dimension(:,:), pointer                  :: CoordX, CoordY
        real,  dimension(:,:), pointer                  :: CoordXout, CoordYout
        integer, dimension(:,:), pointer                :: InPutMap2D
        logical                                         :: Exist
        type (T_PointF),                pointer         :: Point
        
        !------------------------------------------------------------------------

        call ReadHDF5Fields

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%InPutFileName, HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR10'

        allocate(Aux3D     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))  
        allocate(Aux2D     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB          ))  
        
        allocate(Point)
        
        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                           Me%WorkSize%JLB,Me%WorkSize%JUB, 1, Me%KUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR20'


        if (Me%WindowON) then

            ILB = Me%Window(1); IUB = Me%Window(2);        
            JLB = Me%Window(3); JUB = Me%Window(4);                    

        else

            ILB = Me%WorkSize%ILB; IUB = Me%WorkSize%IUB;        
            JLB = Me%WorkSize%JLB; JUB = Me%WorkSize%JUB;                    
        
        endif
        
        if (Me%TransferToOutGrid) then
            
            allocate(Aux3DOut(Me%SizeOut%ILB:Me%SizeOut%IUB, Me%SizeOut%JLB:Me%SizeOut%JUB,1: Me%KUB))
            allocate(Aux2DOut(Me%SizeOut%ILB:Me%SizeOut%IUB, Me%SizeOut%JLB:Me%SizeOut%JUB))
            
            allocate(Me%InPutMap (Me%Size%ILB   :Me%Size%IUB   , Me%Size%JLB   :Me%Size%JUB,1: Me%KUB))
            allocate(   InPutMap2D (Me%Size%ILB   :Me%Size%IUB   , Me%Size%JLB   :Me%Size%JUB   ))            
            allocate(Me%OutPutMap(Me%SizeOut%ILB:Me%SizeOut%IUB, Me%SizeOut%JLB:Me%SizeOut%JUB))
            
            
            ILBout = Me%WorkSizeOut%ILB
            IUBout = Me%WorkSizeOut%IUB
            JLBout = Me%WorkSizeOut%JLB
            JUBout = Me%WorkSizeOut%JUB            

            call GetZCoordinates(Me%ObjHorizontalGrid, CoordX, CoordY)
            call GetZCoordinates(Me%ObjHorizontalGridOut, CoordXout, CoordYout)
            
            if (Me%LandMaskON) then
            
                do j = Me%WorkSize%JLB, Me%WorkSize%JUB
                do i = Me%WorkSize%ILB, Me%WorkSize%IUB

                    Point%X = CoordX(i, j)
                    Point%Y = CoordY(i, j)
                    if (IsVisible(Me%LandMask, Point)) then
                        Me%InPutMap(i, j,:) = 0
        else
                        Me%InPutMap(i, j,:) = 1
                    endif
                
                enddo
                enddo
            
                do j = Me%WorkSizeOut%JLB, Me%WorkSizeOut%JUB
                do i = Me%WorkSizeOut%ILB, Me%WorkSizeOut%IUB
                
                    Point%X = CoordXout(i, j)
                    Point%Y = CoordYout(i, j)
                    if (IsVisible(Me%LandMask, Point)) then
                        Me%OutPutMap(i, j) = 0
                    else
                        Me%OutPutMap(i, j) = 1
                    endif
                
                enddo
                enddo
                
            endif
        else
            Aux2DOut => Aux2D
            Aux3DOut => Aux3D
            Me%ObjHorizontalGridOut = Me%ObjHorizontalGrid
            ILBout = ILB
            IUBout = IUB
            JLBout = JLB
            JUBout = JUB
            
        endif            
                
        
        k        = Me%Layer


d11:    do l = 1, Me%FieldNumber

                i = scan(Me%FieldName(l),'/',Back = .true.)

                AuxChar = trim(Me%FieldName(l))
                
                VGroup  = AuxChar(1:i-1)
                
                j = len_trim(Me%FieldName(l))
                
                Field  = AuxChar(i+1:j)

                if (Me%ComputeOption == InstField_) then
                    
                    call GetHDF5DataSetExist (Me%ObjHDF5, DataSetName =trim(Me%FieldName(l)),&
                                              Exist = Exist, STAT= STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR30'

                    if (.not.Exist) then
                        write(*,*) 'The field'
                        write(*,*) trim(Me%FieldName(l))
                        write(*,*) 'does not exist'
                        stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR40'
                    endif      
                
                    call HDF5ReadData(Me%ObjHDF5,                                           &
                                       trim(Vgroup),                                        &
                                       trim(Field),                                         &
                                       Array3D      = Aux3D,                                &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR50'
                
                else
                    
                    call ReadHDF5_FieldInTime (Aux3D, Vgroup, Field)                    
                    
                endif
                
                where (Aux3D < Me%MinValue) Aux3D = Me%FillValue
                
                where (Aux3D > Me%MaxValue) Aux3D = Me%FillValue
                
                
                
                where (Aux3D > Me%FillValue) Aux3D = Aux3D * Me%MultiplyFactor + Me%AddFactor
                    
                call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR60'

                open(Unit   = Unit,                                                     &
                     File   = trim(Me%OutputESRI(l)),                                   &
                     Form   = 'FORMATTED',                                              &
                     STATUS = 'UNKNOWN',                                                &
                     Action = 'WRITE',                                                  &
                     IOSTAT = STAT_CALL) 
                if (STAT_CALL /= SUCCESS_) stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR70'
                

                write(Unit,'(A14,I6)'   ) 'ncols         ', JUBout - JLBout + 1
                write(Unit,'(A14,I6)'   ) 'nrows         ', IUBout - ILBout + 1
                write(Unit,'(A14,A30)')     'xllcorner     ', trim(adjustl(Me%Origin(1)))
                write(Unit,'(A14,A30)')     'yllcorner     ', trim(adjustl(Me%Origin(2)))
                write(Unit,'(A14,A30)')     'cellsize      ', trim(adjustl(Me%DX))
                write(Unit,'(A14,f12.6)') 'nodata_value  ', Me%FillValue
                
                if (Me%TransferToOutGrid) then
                
                    Aux2DOut(:,:) = Me%FillValue
                    Aux3DOut(:,:,:) = Me%FillValue
     
                Aux2D   (:,:  ) = Aux3D(:,:,k)
                if (Me%LandMaskON) then
                    InPutMap2D   (:,:)  = Me%InPutMap(:,:,k)                
                else
                    InPutMap2D   (:,:)  = 1
                    Me%OutPutMap (:,:)  = 1
                endif
                
                where (Aux2D == Me%FillValue) InPutMap2D = 0
                    
                if (Me%InterpolOut) then

                    !if (Me%Regular) then
                    !
                    !    call InterpolRegularGrid(HorizontalGridSonID      = Me%ObjHorizontalGridOut,&
                    !                             HorizontalGridFatherID   = Me%ObjHorizontalGrid,   &
                    !                             Field2DFather            = Aux2D,                  & 
                    !                             Field2DSon               = Aux2DOut,               &
                    !                             ComputeFather            = Me%InPutMap,            &
                    !                             KUBFather                = k,                      &
                    !                             STAT                     = STAT_CALL)
                    !    if (STAT_CALL /= SUCCESS_)stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR80'
                    !
                    !else
                    

                        
                    do j = Me%WorkSizeOut%JLB, Me%WorkSizeOut%JUB
                    do i = Me%WorkSizeOut%ILB, Me%WorkSizeOut%IUB

                        Aux2DOut(i,j) = InterpolXYPoint(HorizontalGridID = Me%ObjHorizontalGrid,    &
                                             Field2DFather            = Aux2D,                  & 
                                                        ComputeFather    = InPutMap2D,      &
                                                        XInput           = CoordXout(i, j), &
                                                        YInput           = CoordYout(i, j), &
                                             STAT                     = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_)stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR80'
                    
                    enddo
                    enddo
                
                    !endif
                    
                    Aux3DOut(:,:,k) = Aux2DOut(:,:)
                    
                    do j = Me%WorkSizeOut%JLB, Me%WorkSizeOut%JUB
                    do i = Me%WorkSizeOut%ILB, Me%WorkSizeOut%IUB
                        
                        if (Me%OutPutMap(i, j) == 0) then
                            Aux2DOut(i,j  ) = Me%FillValue
                            Aux3DOut(i,j,:) = Me%FillValue
                        endif
                
                    enddo
                    enddo                    
                    
                else
                        
                    do i = ILB, IUB
                        do j = JLB, JUB
                            call GetXYCellZ(Me%ObjHorizontalGridOut, CoordX(i, j), CoordY(i, j), iout,jout)
                            if (iout > 0 .and. jout>0) then
                            Aux2DOut(iout,jout)   = Aux3D(i,j,k)
                            Aux3DOut(iout,jout,k) = Aux3D(i,j,k)
                            endif
                        enddo
                    enddo
                    
                endif
               
            endif
                
            call WriteHDF5_To_GridData(Aux2DOut, Me%OutputESRI(l))            
                
                
                if (Me%ExportXYZ) then
                
                call Export_To_XYZ(Aux3DOut, k, ILBout, IUBout, JLBout, JUBout, Me%OutputESRI(l))
                
                endif 
                
                do i = IUBout, ILBout, -1
                    do j = JLBout, JUBout
                        if (Aux3DOut(i,j,k) <=Me%FillValue) Aux3DOut(i,j,k) = Me%FillValue
                    enddo
                    write(Line,'(4000(e12.6,1x))') Aux3Dout(i,JLBout:JUBout,k)
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
            if (STAT_CALL /= SUCCESS_) stop 'ModifyHDF5_2_EsriGridData - HDF5_2_EsriGridData - ERR90'

        enddo d11
        
        deallocate(Aux2D)  
        deallocate(Aux3D)  

        if (Me%TransferToOutGrid) then
            deallocate(Aux3DOut)
            deallocate(Aux2DOut)
            deallocate(InPutMap2D)
        endif            
        

    end subroutine ModifyHDF5_2_EsriGridData
 
    !--------------------------------------------------------------------------    

    subroutine ModifyEsriGridData_2_HDF5

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real,    dimension(:,:,:), pointer    :: Aux3D
        real,    dimension(:,:  ), pointer    :: Aux2D        
        integer, dimension(:,:,:), pointer    :: WaterPoints    
        integer, dimension(:,:) ,  pointer    :: WaterPoints2D        
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
        

        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR30'

        if (Me%Dim2D) then
        
            allocate(Aux2D        (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))  
            allocate(WaterPoints2D(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
        
        else

            Me%KUB = 1
            allocate(Aux3D      (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))  
            allocate(WaterPoints(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))
        
        endif

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
        
        if (Me%Dim2D) then        
        
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                               Me%WorkSize%JLB,Me%WorkSize%JUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR60'
      
        else
        
            call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,                &
                               Me%WorkSize%JLB,Me%WorkSize%JUB, 1, Me%KUB, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR65'
        endif
        
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

            FileName = trim(AuxChar)//trim(PropName)//".asc"

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
            
            if (Me%Dim2D) then
                do i = IUB, ILB, -1
                    read(Unit,*) (Aux2D(i,j),j=JLB,JUB)
                    if (Me%NoNegativeValues) then
                        do j=JLB,JUB
                            if (Aux2D(i,j)<0.) then
                                Aux2D(i,j) = 0.
                            endif
                        enddo  
                    endif                      
                enddo
                
                call HDF5WriteData(Me%ObjHDF5,                                              &
                                   GroupName    = "/Results/"//PropName,                    &
                                   Name         = PropName,                                 &
                                   Units        = '-',                                      &       
                                   Array2D      = Aux2D,                                    &
                                   OutputNumber = it,                                       &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR125'

                
            else                        
                do i = IUB, ILB, -1
                    read(Unit,*) (Aux3D(i,j,k),j=JLB,JUB)
                    if (Me%NoNegativeValues) then
                        do j=JLB,JUB
                            if (Aux3D(i,j,k)<0.) then
                                Aux3D(i,j,k) = 0.
                            endif
                        enddo  
                    endif                      
                enddo
                
                call HDF5WriteData(Me%ObjHDF5,                                              &
                                   GroupName    = "/Results/"//PropName,                    &
                                   Name         = PropName,                                 &
                                   Units        = '-',                                      &       
                                   Array3D      = Aux3D,                                    &
                                   OutputNumber = it,                                       &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR130'
                
            endif
            
            
            call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR120'

            
            
            if (FirstTime) then
            
                if (Me%Dim2D) then            
                    WaterPoints2D(:,:) = 0
                    
                    do j = JLB, JUB
                    do i = ILB, IUB
                        if (Aux2D(i,j) /= nodata_value) then
                            WaterPoints2D(i,j) = 1
                        endif                        
                    enddo
                    enddo
                    
                    Bathymetry(:,:) = WaterPoints2D(:,:)

                    call HDF5WriteData(Me%ObjHDF5,                                          &
                                       GroupName    = "/Grid",                              &
                                       Name         = "WaterPoints2D",                      &
                                       Units        = '-',                                  &       
                                       Array2D      = WaterPoints2D,                        &
                                       STAT         = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'ModifyEsriGridData_2_HDF5 - HDF5_2_EsriGridData - ERR150'

                    deallocate(WaterPoints2D)

                else
                    WaterPoints(:,:,:) = 0
                    
                    do j = JLB, JUB
                    do i = ILB, IUB
                        if (Aux3D(i,j,k) /= nodata_value) then
                            WaterPoints(i,j,k) = 1
                        endif                        
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

                endif

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
        
        if (me%Dim2D) then
            deallocate(Aux2D)
        else
            deallocate(Aux3D)
        endif                        

    end subroutine ModifyEsriGridData_2_HDF5
 
    !--------------------------------------------------------------------------

    subroutine Export_To_XYZ(Aux3D, k, ILB, IUB, JLB, JUB, Filename)

        !Arguments---------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D
        integer                                         :: i, j, STAT_CALL, Unit, SplitByExtension
        integer                                         :: k, ILB, IUB, JLB, JUB
        real,  dimension(:,:), pointer                  :: CoordX, CoordY
        character(len=PathLength)                       :: Filename
        
        !------------------------------------------------------------------------

        write(*,*)"Writing MOHID .xyz file..."

        call GetZCoordinates(Me%ObjHorizontalGridOut, CoordX, CoordY)
              
        SplitByExtension = scan(trim(Filename),'.',Back = .true.)
                    
        open(Unit   = Unit,                                                 &
         File   = trim(Filename(1:SplitByExtension-1)//'.xyz'),             &
         Form   = 'FORMATTED',                                              &
         STATUS = 'UNKNOWN',                                                &
         Action = 'WRITE',                                                  &
         IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Export_To_XYZ - HDF5_2_EsriGridData - ERR10'
                    
        write(Unit,*)'<begin_xyz>'
        
        do i = ILB, IUB
            do j = JLB, JUB
                if (Aux3D(i,j,k) <= FillValueReal/1e4 .and. Aux3D(i,j,k) /= Me%FillValue) then
                    write(Unit,*)CoordX(i, j), CoordY(i, j), Me%FillValue
                else
                    write(Unit,*)CoordX(i, j), CoordY(i, j), Aux3D(i,j,k)
                endif
            enddo
        enddo
        
        write(Unit,*)'<end_xyz>'
                        
        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Export_To_XYZ - HDF5_2_EsriGridData - ERR20'
        
    end subroutine Export_To_XYZ
 
    !--------------------------------------------------------------------------
    
    subroutine WriteHDF5_To_GridData (Aux2DOut, Filename)
        
        !External--------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2DOut
        integer                                         :: STAT_CALL, SplitByExtension
        character(len=PathLength)                       :: Filename
        
        !Begin-----------------------------------------------------------------

        SplitByExtension = scan(trim(Filename),'.',Back = .true.)
        
        write(*,*)"Writing griddata file..."
        
        call WriteGridData(FileName         = trim(Filename(1:SplitByExtension-1)//'.dat'),        &
               COMENT1          = 'File generated by',                                             &
               COMENT2          = 'HDF5_2_EsriGridData',                                           &
               HorizontalGridID = Me%ObjHorizontalGridOut,                                         &
               FillValue        = Me%FillValue,                                                    &
               Overwrite        = .true.,                                                          &
               GridData2D_Real  = Aux2DOut,                                                        &
               STAT             = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'WriteHDF5_To_GridData - HDF5_2_EsriGridData - ERR10'
        
    end subroutine WriteHDF5_To_GridData
    
    !--------------------------------------------------------------------------   
    
    subroutine ReadHDF5_FieldInTime (Aux3D, Vgroup, Field)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3D
        character (len=*)                               :: Vgroup, Field
        
        !Local-----------------------------------------------------------------
        real, dimension(:,:,:), pointer                 :: Aux3DInst        
        real                                            :: dt, dt_sum, dt_old
        integer                                         :: i, j, k, n
        integer                                         :: NumberOfInstants, STAT_CALL
        type (T_Time)                                   :: CurrentTime, NextTime
        logical                                         :: Exist
        
        !Begin-----------------------------------------------------------------
        
        call GetHDF5GroupExist (HDF5ID      = Me%ObjHDF5,                               &
                                GroupName   = trim(Vgroup),                             &
                                Exist       = Exist,                                    &
                                STAT        = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_FieldInTime - HDF5_2_EsriGridData - ERR10'
                    
        if (.not.Exist) then
            write(*,*) 'The Vgroup'
            write(*,*) trim(Vgroup)
            write(*,*) 'does not exist'
            stop 'ReadHDF5_FieldInTime - HDF5_2_EsriGridData - ERR20'
        endif                     
        
    
        allocate(Aux3DInst     (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB,1: Me%KUB))
        
        Aux3D    (:,:,:) = 0.
        Aux3DInst(:,:,:) = 0.
    
        call GetHDF5GroupNumberOfItems(Me%ObjHDF5, "/Time",                             &
                                       NumberOfInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_FieldInTime - HDF5_2_EsriGridData - ERR030'

        dt_sum = 0.
        dt_old = 0
        
        do n=1, NumberOfInstants
            
            CurrentTime = HDF5TimeInstant(n  , Me%ObjHDF5)
            
            if (n < NumberOfInstants) then
                NextTime    = HDF5TimeInstant(n+1, Me%ObjHDF5)
            else
                NextTime    = CurrentTime
            endif
            
            if (NextTime > Me%EndTime) then
                NextTime = CurrentTime
            endif
            
            if (CurrentTime >= Me%StartTime .and. CurrentTime <= Me%EndTime) then
                
                call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,        &
                                   Me%WorkSize%JLB,Me%WorkSize%JUB, 1, Me%KUB, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'ReadHDF5_FieldInTime - HDF5_2_EsriGridData - ERR40'
                
                call HDF5ReadData(HDF5ID         = Me%ObjHDF5,                          &
                                  GroupName      = trim(Vgroup),                        &
                                  Name           = trim(Field),                         &
                                  Array3D        = Aux3DInst,                           &
                                  OutputNumber   = n,                                   &
                                  STAT           = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'ReadHDF5_FieldInTime - HDF5_2_EsriGridData - ERR50'
                
                dt      = dt_old + (NextTime - CurrentTime)/2
                dt_sum  = dt_sum + dt
                
                if (Me%ComputeOption == AverageField_) then
                                
                    do i=Me%Size%ILB,Me%Size%IUB
                    do j=Me%Size%JLB,Me%Size%JUB
                    do k=1          ,Me%KUB
                        Aux3D(i, j, k) =  Aux3D(i, j, k) + Aux3DInst(i, j, k) * dt 
                    enddo
                    enddo
                    enddo
                    
                elseif (Me%ComputeOption  == SumField_) then

                    do i=Me%Size%ILB,Me%Size%IUB
                    do j=Me%Size%JLB,Me%Size%JUB
                    do k=1          ,Me%KUB
                        Aux3D(i, j, k) =  Aux3D(i, j, k) + Aux3DInst(i, j, k)
                    enddo
                    enddo
                    enddo                    
                
                endif
            endif
            
        enddo
        
        if (Me%ComputeOption == AverageField_) then
            if (dt_sum > 0.) then
            
                do i=Me%Size%ILB,Me%Size%IUB
                do j=Me%Size%JLB,Me%Size%JUB
                do k=1          ,Me%KUB
                    Aux3D(i, j, k) = Aux3D(i, j, k) / dt_sum
                enddo
                enddo
                enddo            
            else
                
                stop 'ReadHDF5_FieldInTime - HDF5_2_EsriGridData - ERR60'             
                
            endif
        endif        
        
        deallocate(Aux3DInst)
        
        
    end subroutine ReadHDF5_FieldInTime
    
   !----------------------------------------------------------------------------


    type(T_Time) function HDF5TimeInstant(Instant, ObjHDF5)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        integer                                 :: ObjHDF5

        !Local-----------------------------------------------------------------
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------

        call HDF5SetLimits  (ObjHDF5, 1, 6, STAT = STAT_CALL)

        allocate(TimeVector(6))

        call HDF5ReadData   (HDF5ID         = ObjHDF5,                                  &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - HDF5_2_EsriGridData - ERR10'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))


        deallocate(TimeVector)

    end function HDF5TimeInstant
    
    
    End Program HDF5_2_EsriGridData
    
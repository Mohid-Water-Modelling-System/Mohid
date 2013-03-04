!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : SWAN
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group & MeteoGalicia
! DATE          : September 2003
! REVISION      : Joao Ribeiro - v1.0
! DESCRIPTION   : Module to convert SWAN NonStationary run table format
!                 files into HDF5 format. For reading into Mohid module
!                  HydrodynamicFile
!
!------------------------------------------------------------------------------


Module ModuleReadSWANNonStationary

    use ModuleTime
    use ModuleGlobalData
    use ModuleFunctions
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleHorizontalGrid
    use ModuleGridData
    use ModuleHorizontalMap

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertSWANNonStationary
    private ::      ReadGlobalOptions
    private ::      ReadProperties
    private ::      ConstructBathymetry
    private ::      Open_HDF5_OutPut_File
    private ::      OutputFields
    private ::      KillSWAN

    !Parameters----------------------------------------------------------------

    character(LEN = StringLength), parameter    :: prop_block_begin   = '<<beginproperty>>'
    character(LEN = StringLength), parameter    :: prop_block_end     = '<<endproperty>>'

    character(LEN = StringLength), parameter    :: units_block_begin  = '<<beginunits>>'
    character(LEN = StringLength), parameter    :: units_block_end    = '<<endunits>>'

    character(LEN = StringLength), parameter    :: convert_block_begin= '<<beginconvert>>'
    character(LEN = StringLength), parameter    :: convert_block_end  = '<<endconvert>>'

    integer,                       parameter    :: NoConversion_          = 0
    
    !Types---------------------------------------------------------------------
    
    type       T_OutPut                                 
         type (T_Time), pointer, dimension(:)   :: OutTime
         integer                                :: Number
         integer                                :: NextOutPut
         logical                                :: ON
         logical                                :: Yes                  = .false.
         integer                                :: TotalOutputs
    end type T_OutPut                                   
    
    private :: T_SWAN
    type       T_SWAN
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ClientNumber
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: ObjBathymetry        = 0
        integer                                 :: ObjHorizontalMap     = 0
        integer                                 :: ObjTime              = 0
        integer                                 :: Unit
        integer                                 :: NauticalToCartesianDegrees
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: SwanFileName
        character(len=PathLength)               :: OutputFileName
        character(len=StringLength)             :: Generic4DProperty
        character(len=StringLength)             :: Generic4DUnits
        
        integer                                 :: NumberDates          = FillValueInt
        integer                                 :: NumberFields         = FillValueInt
        integer                                 :: NumberProps          = FillValueInt
        integer                                 :: NumberUnits          = FillValueInt
        integer                                 :: NumberCovertions     = FillValueInt

        integer                                 :: DirectionReferential = FillValueInt
        
        character(len=StringLength), dimension(:), pointer :: PropsName
        character(len=StringLength), dimension(:), pointer :: PropsUnits
        integer,                     dimension(:), pointer :: ConvertProp

        type(T_Time)                            :: BeginTime, EndTime
        real,   dimension(:,:,:,:),     pointer :: Fields 
        real,   dimension(:),           pointer :: PropVector

        integer                                 :: Clientumber

        logical                                 :: WriteVelModulus = .false., WriteWindModulus = .false.
        real                                    :: FillValue
        real                                    :: Generic4DValue
        integer                                 :: Generic4DPropertyIndeX
        integer                                 :: Generic4D
        integer                                 :: ReadType
        integer, dimension(:,:  ),  pointer     :: WaterPoints2D
        type(T_OutPut)                          :: OutPut
        type(T_Size2D)                          :: WorkSize, Size
    end type  T_SWAN

    type(T_SWAN), pointer              :: Me

    !--------------------------------------------------------------------------
    
    contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertSWANNonStationary(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !Local-------------------------------------------------------------------
!        integer                                         :: l
        integer                                         :: p
        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        Me%ClientNumber = ClientNumber

        call ReadGlobalOptions

        call ConstructBathymetry

        call ReadProperties

        call ReadConversion

        call Open_HDF5_OutPut_File

        call ReadFieldFromFile
            
        Me%OutPut%NextOutPut = 1.0

d1:     do p=1, Me%NumberProps

           call OutputFields (p)

                !if (Me%WriteVelModulus) then
                    !call WriteVelocityModulus(VelocityU_, VelocityV_, VelocityModulus_)
                !endif

                !if (Me%WriteWindModulus) then
                    !call WriteVelocityModulus(WindVelocityX_, WindVelocityY_, WindModulos_)
                !endif

        enddo d1

        call KillSWAN

        STAT = SUCCESS_

    end subroutine ConvertSWANNonStationary

    !------------------------------------------------------------------------

    subroutine ReadGlobalOptions

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag, iflag1

        !Begin-----------------------------------------------------------------

        call GetData(Me%GridFileName,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUT_GRID_FILENAME',                              &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR10'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR20'

        call GetData(Me%OutputFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR30'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR40'

        call GetData(Me%SwanFileName,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'INPUT_SWAN_FILENAME',                                   &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR50'
        if (iflag     == 0)        stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR60'
       
        call GetData(Me%FillValue,                                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'FILL_VALUE',                                       &
                     default      = -99.999900,                                         &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR70'
        
        call GetData(Me%Generic4D,                                            &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'Generic4D',                              &
                     default      = 0,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR80'
        
        call GetData(Me%Generic4DProperty,                                    &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'Generic4D_Property',                     &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR90'
        if (iflag     == 0 .AND. Me%Generic4D == 1 )                          &
                stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR100'

        call GetData(Me%BeginTime,                                      &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'START',                            &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR110'

        call GetData(Me%EndTime,                                        &
                     Me%ObjEnterData, iflag1,                           &
                     SearchType   = FromBlock,                          &
                     keyword      = 'END',                              &
                     ClientModule = 'ConvertToHDF5',                    &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR120'

        if (iflag==0 .and. iflag1==0)                                   &
                stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR130'
                
        call GetOutPutTime(Me%ObjEnterData,                                   &
                     CurrentTime      = Me%BeginTime,                         &
                     EndTime          = Me%EndTime,                           &
                     keyword          = 'OUTPUT_TIME',                        &
                     SearchType       = FromFile,                             &
                     OutPutsTime      = Me%OutPut%OutTime,                    &
                     OutPutsOn        = Me%OutPut%Yes,                        &
                     OutPutsNumber    = Me%OutPut%TotalOutputs,               &
                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                            &
                stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR140' 

        call GetData(Me%NauticalToCartesianDegrees,                           &
                     Me%ObjEnterData, iflag,                                  &
                     SearchType   = FromBlock,                                &
                     keyword      = 'NAUTICAL_TO_CARTESIAN',                  &
                     default      = 0,                                        &
                     ClientModule = 'SWAN',                                   &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_)                                            &
         stop 'ReadGlobalOptions - ModuleReadSWANNonStationary - ERR150'
        
    end subroutine ReadGlobalOptions

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadProperties

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL, Line, FirstLine, LastLine
        character(len=StringLength)                 :: AuxChar
        logical                                     :: BlockFound
        !Begin-----------------------------------------------------------------


        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = prop_block_begin,               &
                                    block_end         = prop_block_end,                 &
                                    BlockInBlockFound = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd1 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :       if (.not. BlockFound) then                                                  
                stop 'ReadProperties - ModuleReadSWANNonStationary - ERR10'
            end if cd2
        else cd1
            stop 'ReadProperties - ModuleReadSWANNonStationary - ERR20'
        end if cd1


        Me%NumberProps = LastLine - FirstLine - 1

        allocate(Me%PropsName (Me%NumberProps))  
  
        allocate(Me%PropVector(Me%NumberProps))  


d1:     do l= 1, Me%NumberProps

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR30'

            Me%PropsName(l) = AuxChar

            if (.not. CheckPropertyName (AuxChar)) then
                write(*,*) 'The name ',trim(AuxChar),' is not valid name for the MOHID system'
            endif

        enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR40'

        call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = units_block_begin,              &
                                    block_end         = units_block_end,                &
                                    BlockInBlockFound = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd3 :   if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd4 :       if (.not. BlockFound) then                                                  
                stop 'ReadProperties - ModuleReadSWANNonStationary - ERR50'
            end if cd4
        else cd3
            stop 'ReadProperties - ModuleReadSWANNonStationary - ERR60'
        end if cd3

        Me%NumberUnits = LastLine - FirstLine - 1

        allocate(Me%PropsUnits(Me%NumberUnits))

d2:     do l= 1, Me%NumberUnits

            line = FirstLine + l

            call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR70'

            Me%PropsUnits(l) = AuxChar

        enddo d2

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadProperties - ModuleReadSWANNonStationary - ERR80'

    end subroutine ReadProperties 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadConversion

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: l, iflag, STAT_CALL
        integer                                     :: ConversionType, j
        character(len=StringLength)                 :: PropName, AuxChar
        logical                                     :: BlockFound
        integer                                     :: Line, FirstLine, LastLine
        !Begin-----------------------------------------------------------------

        if (Me%NauticalToCartesianDegrees) then

             allocate(Me%ConvertProp(Me%NumberProps))  
  
             Me%ConvertProp(:) = NoConversion_
        
             call ExtractBlockFromBlock (Me%ObjEnterData,                                    &
                                    ClientNumber      = Me%ClientNumber,                &
                                    block_begin       = convert_block_begin,            &
                                    block_end         = convert_block_end,              &
                                    BlockInBlockFound = BlockFound,                     &
                                    FirstLine         = FirstLine,                      &
                                    LastLine          = LastLine,                       &
                                    STAT              = STAT_CALL)

cd1 :        if      (STAT_CALL .EQ. SUCCESS_     ) then    
cd2 :            if (.not. BlockFound) then                                                  
                     stop 'ReadConversion - ModuleReadSWANNonStationary - ERR10'
                 end if cd2
             else cd1
                 stop 'ReadConversion - ModuleReadSWANNonStationary - ERR20'
             end if cd1

             Me%NumberCovertions = LastLine - FirstLine - 1

d1:          do l= 1, Me%NumberCovertions

                 line = FirstLine + l

                 call GetData(AuxChar, EnterDataID = Me%ObjEnterData, flag = iflag,          &
                         SearchType = FromBlock, Buffer_Line = line, STAT = STAT_CALL)

                 if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleReadSWANNonStationary - ERR30'

                 do j = 1, Me%NumberProps
                     if (trim(Me%PropsName(j))==trim(AuxChar)) then
                         Me%ConvertProp(j) = 1
!                         write(*,*)  trim(AuxChar)
                         exit
                     endif
                 enddo

             enddo d1

        call RewindBlock(Me%ObjEnterData, Me%ClientNumber, STAT = STAT_CALL)
        
        if (STAT_CALL /= SUCCESS_) stop 'ReadConversion - ModuleReadSWANNonStationary - ERR40'
        endif
        
    end subroutine ReadConversion 

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine ReadFieldFromFile 

        !Arguments-------------------------------------------------------------

        !Local----------------------------------------------------------------
        integer                                     :: i, j, p, STAT_CALL,k
        !Begin-----------------------------------------------------------------
        
        allocate(Me%Fields(Me%OutPut%TotalOutputs,Me%NumberProps,                      &
        Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))

        call UnitsManager(Me%Unit, OPEN_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleReadSWANNonStationary - ERR10'

        open(Unit   = Me%Unit,                                                          &
             File   = trim(Me%SwanFileName),                                                  &
             Form   = 'FORMATTED',                                                      &
             STATUS = 'OLD',                                                            &
             Action = 'READ',                                                           &
             IOSTAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleReadSWANNonStationary - ERR20'
        
        do k=1,Me%OutPut%TotalOutputs
        
             do i=Me%WorkSize%ILB, Me%WorkSize%IUB
             do j=Me%WorkSize%JLB, Me%WorkSize%JUB

                 read(Me%Unit,*) Me%PropVector

                 do p = 1, Me%NumberProps
 
                      if(trim(Me%Generic4DProperty) == trim(Me%PropsName(p))) Then 
                
                           Me%Generic4DPropertyIndeX=p
                
                      end if
                      
                      Me%Fields(k, p, i, j) = Me%PropVector(p)
                
                      if (Me%Fields(k, p, i, j).eq.-9.0 .OR. Me%Fields(k, p, i, j).eq.-99.0 &
                          .OR. Me%Fields(k, p, i, j).eq.-999.0) Then
                       
                           Me%Fields(k, p, i, j) = FillValueInt
                       
                      end if
                
                      if(Me%Generic4D==1 .AND. p==Me%Generic4DPropertyIndeX) Then
                      Me%Generic4DValue=Me%PropVector(p)
                      Me%Generic4DUnits=trim(Me%PropsUnits(p))
                      end if

                 enddo

             enddo
             enddo
             
        enddo

        call UnitsManager(Me%Unit, CLOSE_FILE, STAT = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'ReadFieldFromFile - ModuleReadSWANNonStationary - ERR30'

    end subroutine ReadFieldFromFile
    
    !----------------------------------------------------------------------------
   
    !------------------------------------------------------------------------
    
    subroutine ConstructBathymetry
        
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        
        !Begin-----------------------------------------------------------------

        !Me%GridFileName="NewGrid.dat_.new"
        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR10'

        call GetHorizontalGridSize(Me%ObjHorizontalGrid, Me%Size, Me%WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR20'

        call ConstructGridData(Me%ObjBathymetry, Me%ObjHorizontalGrid, FileName = Me%GridFileName,&
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR30'


        call ConstructHorizontalMap(Me%ObjHorizontalMap, Me%ObjBathymetry, Me%ObjHorizontalGrid, &
                               STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR40'

        call GetWaterPoints2D   (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'ConstructBathymetry - ModuleReadSWANNonStationary - ERR50'

    end subroutine ConstructBathymetry
    
    !------------------------------------------------------------------------
    
    !------------------------------------------------------------------------
    
    subroutine OutputFields(p)

        !Arguments-------------------------------------------------------------
        integer                                         :: p
        !Local-----------------------------------------------------------------
        real, dimension(:,:), pointer                   :: Aux2D, Aux2DX, Aux2DY
        real, dimension(:,:), pointer                   :: Aux2DXY
        real,    dimension(:), pointer                  :: RealArray1D
        REAL, PARAMETER                                 :: RADE= 180./pi
        real                                            :: Angle
        integer                                         :: STAT_CALL, i, ii, jj
        integer                                         :: kk, k
        character(len=StringLength)                     :: FieldName
        !Begin-----------------------------------------------------------------

        if (p==1) i = Me%OutPut%NextOutPut

        if (i==1 .and. p==1) then
            allocate(Aux2D (Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Aux2DX(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(Aux2DY(Me%Size%ILB:Me%Size%IUB, Me%Size%JLB:Me%Size%JUB))
            allocate(RealArray1D(1:1))
        endif
            
i2:     if (p==1 .AND. Me%Generic4D==1) then
                
            RealArray1D(1) = Me%Generic4DValue
                
            call HDF5SetLimits  (Me%ObjHDF5, 1, 1, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR30'
                
            call HDF5WriteData  (Me%ObjHDF5, "/Generic4D",                           &
                                 "Generic4D", Me%Generic4DUnits,                     &
                                 Array1D = RealArray1D,                              &
                                 OutputNumber = k, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR40'
                
        endif i2

        call HDF5SetLimits(Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,            &
                           Me%WorkSize%JLB, Me%WorkSize%JUB,                        &
                           STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR60'
            
d1:      do kk=1,Me%OutPut%TotalOutputs

             Aux2D(:,:) = Me%Fields(kk, p,:,:)
             
             if (Me%ConvertProp(p) == Me%NauticalToCartesianDegrees) then
             
!                  write(*,*)  trim(Me%PropsName(p))
                  
                  do jj = Me%WorkSize%JLB, Me%WorkSize%JUB
                  do ii = Me%WorkSize%ILB, Me%WorkSize%IUB
                 
                     if (Aux2D(ii,jj) /= FillValueInt)  then
                          Aux2D(ii,jj) = MOD(630.-(Aux2D(ii,jj)*(pi/180))*RADE,360.)
                     endif
                
                  enddo
                  enddo
                  
             endif
             
             call HDF5WriteData(Me%ObjHDF5,                                         &
                           "/Results/"//trim(Me%PropsName(p)),                      &
                           trim(Me%PropsName(p)),                                   &
                           trim(Me%PropsUnits(p)),                                  &
                           Array2D      = Aux2D,                                    &
                           OutputNumber = kk,                                       &
                           STAT         = STAT_CALL)
             if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR70'
            
       
            
ip:     if (trim(Me%PropsName(p)) == GetPropertyName(MeanWaveDirection_)) then

            Aux2DX(:,:) = 0.
            Aux2DY(:,:) = 0.

d2:         do k=1,2

                do jj= Me%WorkSize%JLB,Me%WorkSize%JUB  
                do ii= Me%WorkSize%ILB,Me%WorkSize%IUB  

                    if (Me%WaterPoints2D(ii,jj) == WaterPoint) then

                        Angle = Aux2D(ii,jj)

                        if (k==1) then
                            Aux2DX(ii,jj) = cos(Angle * Pi / 180.) 
                        else
                            Aux2DY(ii,jj) = sin(Angle * Pi / 180.) 
                        endif

                    endif

                enddo
                enddo
           
                if (k==1) then
                    FieldName = trim(Me%PropsName(p))//'_x'
                    Aux2DXY => Aux2DX
                else
                    FieldName = trim(Me%PropsName(p))//'_y'
                    Aux2DXY => Aux2DY
                endif

                call HDF5WriteData(Me%ObjHDF5,                                      &
                                   "/Results/"//FieldName,                          &
                                   FieldName,                                       &
                                   '-',                                             &
                                   Array2D      = Aux2DXY,                          &
                                   OutputNumber = kk,                               &
                                   STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR80'

            enddo d2

        endif ip
        
        enddo d1

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleReadSWANNonStationary - ERR90'

        if (p == Me%NumberProps) then

            Me%OutPut%NextOutPut = i + 1

        endif

        if (i==Me%NumberDates .and. p==Me%NumberProps) then
            deallocate(Aux2D)
            deallocate(Aux2DX)
            deallocate(Aux2DY)
            deallocate(RealArray1D)
            nullify(Aux2D, Aux2DX, Aux2DY, Aux2DXY)
        endif

    end subroutine OutputFields

   !----------------------------------------------------------------------

   !------------------------------------------------------------------------
    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        real,    dimension(6), target               :: AuxTime
        real,    dimension(:), pointer              :: TimePtr
        real,    dimension(:,:), pointer            :: Bathymetry
        integer                                     :: STAT_CALL,i
        integer                                     :: HDF5_CREATE

        !----------------------------------------------------------------------

        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        call GetGridData        (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR10'

      
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR30'
        
        
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR40'

            
        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "m",                   &
                              Array2D =  Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR50'            


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR60'            
   
       
        call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,              &
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR70'            

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints2D", "-",                &
                              Array2D = Me%WaterPoints2D,  STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR80'


        call UnGetGridData      (Me%ObjBathymetry, Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR90'


        do i=1,Me%Output%TotalOutputs

             call ExtractDate   (Me%Output%OutTime(i),                               &
                                AuxTime(1), AuxTime(2), AuxTime(3),                 &
                                AuxTime(4), AuxTime(5), AuxTime(6))

             TimePtr => AuxTime

             call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
             if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR100'


             call HDF5WriteData  (Me%ObjHDF5, "/Time",                               &
                                 "Time", "YYYY/MM/DD HH:MM:SS",                     &
                                 Array1D = TimePtr,                                 &
                                 OutputNumber = i, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR110'
        
        enddo
        
        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleReadSWANNonStationary - ERR120'

    end subroutine Open_HDF5_OutPut_File

    !--------------------------------------------------------------------------

   
    !--------------------------------------------------------------------------

    
    subroutine KillSWAN
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------

        call UnGetHorizontalMap (Me%ObjHorizontalMap, Me%WaterPoints2D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR10'

        call KillHorizontalMap(Me%ObjHorizontalMap, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR20'

        call KillGridData(Me%ObjBathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR30'

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillSWAN - ModuleReadSWANNonStationary - ERR40'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillSWAN - ModuleReadSWANNonStationary - ERR60'

        deallocate(Me%OutPut%OutTime )
        nullify   (Me%OutPut%OutTime )

        deallocate(Me%PropsName )
        nullify   (Me%PropsName )

        deallocate(Me%PropsUnits)
        nullify   (Me%PropsUnits)

        deallocate(Me%ConvertProp)
        nullify   (Me%ConvertProp)

        deallocate(Me%Fields    )
        nullify   (Me%Fields    )

        deallocate(Me%PropVector)
        nullify   (Me%PropVector)

        deallocate(Me)
        nullify   (Me)
    
    end subroutine KillSWAN

    !--------------------------------------------------------------------------
 
end module ModuleReadSWANNonStationary

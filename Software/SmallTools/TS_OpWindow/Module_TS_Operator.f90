!------------------------------------------------------------------------------
!        HIDROMOD : Modelação em Engenharia
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Time Serie
! URL           : http://www.mohid.com
! AFFILIATION   : HIDROMOD
! DATE          : Nov2012
! REVISION      : Paulo Leitão - v1.0
! DESCRIPTION   : Module to analyse (filter, interpolate, identify patterns, compare) Time Series
!
!------------------------------------------------------------------------------
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License 
!version 2, as published by the Free Software Foundation.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
!
!------------------------------------------------------------------------------

Module Module_TS_Operator

    use interpreter
    use ModuleGlobalData
    use ModuleFunctions
    use ModuleEnterData
    use ModuleTime
    use ModuleTimeSerie
    use Module_TS_Synch
    
    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: Construct_TS_Operator
    private ::      AllocateInstance

    !Selector
    public  :: Get_TS_Operator_Values
    public  :: Get_TS_OperatorID
    public  :: UnGet_TS_Operator
                     
    
    !Modifier
    public  :: Modify_TS_Operator

    !Destructor
    public  :: Kill_TS_Operator                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObj_TS_Operator 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGet_TS_Operator1D_R8
    interface  UnGet_TS_Operator
        module procedure UnGet_TS_Operator1D_R8
    end interface  UnGet_TS_Operator


    !Parameter-----------------------------------------------------------------
    
    real(8), parameter  :: Pi_ = 3.1415926535897932384626433832795
    !Input / Output
    integer, parameter  :: FileOpen = 1, FileClose = 0
    
    !Interpolation methods
    integer, parameter  :: LinearTS_ = 1, BackwardTS_ = 2, CumulativeTS_ = 3
    
    !Types---------------------------------------------------------------------
    
    type T_TS_Input
    
        integer                                                 :: Obj_TS_Synch     =  0
        real(8),                    dimension(:),   pointer     :: Values         => null()
	    character(len=StringLength)                             :: Name           =  null_str         
        
    end type T_TS_Input
    
    type T_TS_Operation
        character(len=StringLength)                             :: Name             = null_str
        real(8),                    dimension(:),   pointer     :: TimeSerie        => null()
	    real                                                    :: DT_Synch         = null_real
        real                                                    :: FillValue        = null_real
        type (T_Time)                                           :: BeginTime
        type (T_Time)                                           :: EndTime
	    integer                                                 :: NValues          = null_int    
        
	    real(8)                                                 :: CoordX           = null_real
	    real(8)                                                 :: CoordY           = null_real
        
        
        integer                                                 :: N_TS_Input
        type (T_TS_Input),  dimension(:),   pointer             :: TS_Input         => null()
        
        character(len=1000        ), dimension(:),   pointer    :: Expression       => null()
        
        integer                                                 :: NExpressions     =  null_int
        
        real(8),                    dimension(:,:),  pointer     :: AuxTS            => null()        
        
        character(len=PathLength  )                             :: OutputFileName   =  null_str         
        character(len=PathLength  )                             :: OutputFileNamev2 =  null_str  
        logical                                                 :: Out_V2_On        = .false.  
        integer                                                 :: iOut             =  null_int
        integer                                                 :: iOut_v2          =  null_int        
        
        !Operation window
        logical                                                 :: Operation_Window = .false. 
        real(8)                                                 :: Max_Window       = null_real
        real(8)                                                 :: Min_Window       = null_real        
        real(8)                                                 :: DT_Window        = null_real        
        
        character(len=PathLength  )                             :: OpWindowFileName      =  null_str         
        character(len=PathLength  )                             :: OpWindowFileNamev2    =  null_str          
        character(len=PathLength  )                             :: OpWindowFileNameStart =  null_str
        character(len=PathLength  )                             :: OpWindowFileNameEnd   =  null_str
        
        integer                                                 :: iW               =  null_int
        integer                                                 :: iW_v2            =  null_int
        integer                                                 :: iW_v3            =  null_int        
        integer                                                 :: iW_v4            =  null_int                
        
    end type T_TS_Operation    

    type T_TS_Operator
    
        integer                                                 :: InstanceID       =  0
        
                                               
	    integer                                                 :: ObjEnterData     = 0
        
        integer                                                 :: ExtractType      = null_int
        
        integer                                                 :: ExtractTypeExp   = null_int
        
        integer                                                 :: ExtractTypeOp    = null_int        
        
        integer                                                 :: ExtractTypeInput = null_int
        
        integer                                                 :: ClientNumber         = null_int
        
        type (T_TS_Operation),  dimension(:),   pointer         :: TS_Operation     => null()        
        
        integer                                                 :: N_Operations     = null_int
        
        character(len=StringLength  )                           :: BlockBeginOperator
        character(len=StringLength  )                           :: BlockEndOperator
        
        character(len=StringLength  )                           :: BlockBeginExpression
        character(len=StringLength  )                           :: BlockEndExpression
        
        character(len=StringLength  )                           :: BlockBeginInputTS
        character(len=StringLength  )                           :: BlockEndInputTS


        type(T_TS_Operator), pointer                     :: Next

    end type T_TS_Operator    

    !Global Variables
    type (T_TS_Operator), pointer                        :: FirstObj_TS_Operator
    type (T_TS_Operator), pointer                        :: Me    


    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Construct_TS_Operator(Obj_TS_OperatorID, EnterDataID, ExtractType,       &
                                     ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: Obj_TS_OperatorID 
        integer      ,          optional, intent(IN )   :: EnterDataID
        integer      ,          optional, intent(IN )   :: ExtractType      
        integer      ,          optional, intent(IN )   :: ClientNumber              
        integer      ,          optional, intent(OUT)   :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(m_TS_Operator_)) then
            nullify (FirstObj_TS_Operator)
            call RegisterModule (m_TS_Operator_) 
        endif

        call Ready(Obj_TS_OperatorID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            if (present(EnterDataID)) then
                Me%ObjEnterData      = AssociateInstance (mENTERDATA_,      EnterDataID     )            
            endif                    
                
            if (present(ExtractType))  then
                Me%ExtractType  =  ExtractType
            else 
                Me%ExtractType  = FromFile_
            endif     
            
            if (present(ClientNumber)) then
                Me%ClientNumber     = ClientNumber            
            else
                Me%ClientNumber     = FillValueInt    
            endif
            
            call ReadKeywords
            
            call ConstructOperatorTS
            
            call ConstructOpWindowTS
            
            
            !Returns ID
            Obj_TS_OperatorID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'Module_TS_Operator - Construct_TS_Operator - ERR140' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine Construct_TS_Operator
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Arguments----------------------------------------------------------

        !Local--------------------------------------------------------------
        integer                 :: STAT_CALL, NumberOfBlocks, i
        logical                 :: BlockFound
        !Begin--------------------------------------------------------------
        
        if (Me%ExtractType == FromFile_) then
        
            Me%ObjEnterData = 0
        
            call ConstructEnterData(Me%ObjEnterData, "TS_Operator.dat", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'Module_TS_Operator - ReadKeywords - ERR10'                
            endif
            
        endif
        
        call DefineOperatorBlock(ExtractType    = Me%ExtractType,                       &
                                 block_begin    = Me%BlockBeginOperator,                &
                                 block_end      = Me%BlockEndOperator,                  &
                                 ExtractTypeOp  = Me%ExtractTypeOp)        
        
        call DefineExpressionBlock(ExtractType      = Me%ExtractType,                   &
                                   block_begin      = Me%BlockBeginExpression,          &
                                   block_end        = Me%BlockEndExpression,            &
                                   ExtractTypeExp   = Me%ExtractTypeExp) 
        
        call DefineInputTS_Block  (ExtractType      = Me%ExtractType,                   &
                                   block_begin      = Me%BlockBeginInputTS,             &
                                   block_end        = Me%BlockEndInputTS,               &
                                   ExtractTypeInput = Me%ExtractTypeInput) 
                
        
        call GetNumberOfBlocks (EnterDataID     = Me%ObjEnterData,                      & 
                                BlockBegin      = Me%BlockBeginOperator,                &
                                BlockEnd        = Me%BlockEndOperator,                  &
                                SearchType      = Me%ExtractType,                       &
                                NumberOfBlocks  = NumberOfBlocks,                       &
                                ClientNumber    = Me%ClientNumber,                      &
                                STAT            = STAT_CALL) 
        
        Me%N_Operations = NumberOfBlocks
        
        allocate(Me%TS_Operation(Me%N_Operations))
        
        
        
        do i = 1, Me%N_Operations
            
            call ReadBlock (ExtractType     = Me%ExtractTypeOp,                         &
                            ClientNumber    = Me%ClientNumber,                          &
                            block_begin     = Me%BlockBeginOperator,                    &
                            block_end       = Me%BlockEndOperator,                      &
                            BlockFound      = BlockFound)

            if (BlockFound) then
                call Read_TS_Operator(Me%ExtractTypeOp, i)
                
                call Read_TS_OpWindow(Me%ExtractTypeOp, i)
                
                call Read_Operator_Expressions(Me%ExtractTypeExp, i)
                
                call Read_Operator_InputTS(Me%ExtractTypeInput, i)
            else
                stop 'Module_TS_Operator - ReadKeywords - ERR20'
            endif
            
        enddo
           
        
        if (Me%ExtractType == FromBlock_) then
        
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - ReadKeywords - ERR240'

        endif
    
    end subroutine ReadKeywords
    
    !-------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------

    subroutine Read_TS_Operator(ExtractType, i)

        !Arguments----------------------------------------------------------
        integer                 :: ExtractType, i
        !Local--------------------------------------------------------------
        integer                 :: flag, STAT_CALL
        integer                 :: EnterAux = 0, ExtractAux
        logical                 :: Exist
        
        !Begin--------------------------------------------------------------
        

        
        call GetData(Me%TS_Operation(i)%DT_Synch,                                       &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='DT_SHYNCRONISATION',                                &
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - ReadKeywords - ERR20'
        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR30'
        endif
        
        !Verifies if file exists
        inquire(FILE = "DefaultInput.dat", EXIST = Exist)
        if (Exist) then
            
            EnterAux    = 0
            ExtractAux  = FromFile
            
            call ConstructEnterData(EnterAux, "DefaultInput.dat", STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR40'            
            

        else
            
            EnterAux    = Me%ObjEnterData
            ExtractAux  = ExtractType            
            
        endif
            
        !Reads Begin Time
        call GetData(Me%TS_Operation(i)%BeginTime,                                      &
                        EnterAux,                                                       &
                        flag,                                                           &
                        SearchType   = ExtractAux,                                      &
                        keyword      ='START',                                          &
                        ClientModule ='Module_TS_Operator',                             &
                        STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR50'
        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR60'
        endif        
        
        !Reads End Time
        call GetData(Me%TS_Operation(i)%EndTime,                                        &
                        EnterAux,                                                       &
                        flag,                                                           &
                        SearchType   = ExtractAux,                                      &
                        keyword      ='END',                                            &
                        ClientModule ='Module_TS_Operator',                             &
                        STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR70'

        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR80'
        endif 
        
        if (Exist) then
            
            call KillEnterData(EnterAux, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR100'
            
        endif        
        
        
	    Me%TS_Operation(i)%NValues = (Me%TS_Operation(i)%EndTime - Me%TS_Operation(i)%BeginTime) / &
                                     Me%TS_Operation(i)%DT_Synch + 1
        
        allocate(Me%TS_Operation(i)%TimeSerie(1:Me%TS_Operation(i)%NValues))
        
        call GetData(Me%TS_Operation(i)%FillValue,                                      &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='FILL_VALUE',                                        &
                     default      = FillValueReal,                                      &            
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR150'
        endif
        
        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR160'
        endif              

        call GetData(Me%TS_Operation(i)%OutputFileName,                                 &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='OUTPUT_FILE',                                       &
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR100'
        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR170'
        endif             
        
        call GetData(Me%TS_Operation(i)%OutputFileNameV2,                               &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='OUTPUT_FILE_V2',                                    &
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR180'
        if (flag == 0) then
            Me%TS_Operation(i)%Out_V2_On = .false.
        else
            Me%TS_Operation(i)%Out_V2_On = .true.
        endif                  
        
        call GetData(Me%TS_Operation(i)%CoordX,                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='COORD_X',                                           &
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR190'
        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR130'
        endif              

        call GetData(Me%TS_Operation(i)%CoordY,                                         &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='COORD_Y',                                           &
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_Operator - ERR200'
        if (flag == 0) then
            stop 'Module_TS_Operator - Read_TS_Operator - ERR210'
        endif              
        
    end subroutine Read_TS_Operator

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine Read_TS_OpWindow(ExtractType, i)

        !Arguments----------------------------------------------------------
        integer                 :: ExtractType, i
        !Local--------------------------------------------------------------
        integer                 :: flag, STAT_CALL
        
        !Begin--------------------------------------------------------------
        
        call GetData(Me%TS_Operation(i)%Operation_Window,                               &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = ExtractType,                                        &
                     keyword      ='OPERATION_WINDOW',                                  &
                     default      = .false.,                                            &
                     ClientModule ='Module_TS_Operator',                                &
                     STAT         = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR10'
        
        
        if (Me%TS_Operation(i)%Operation_Window) then
        
            call GetData(Me%TS_Operation(i)%OpWindowFileName,                           &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='OPERATION_WINDOW_FILE',                         &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR20'
            if (flag == 0) then
                stop 'Module_TS_Operator - Read_TS_OpWindow - ERR30'
            endif             
        
            call GetData(Me%TS_Operation(i)%OpWindowFileNameV2,                         &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='OPERATION_WINDOW_FILE_V2',                      &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR40'
            
            call GetData(Me%TS_Operation(i)%OpWindowFileNameStart,                      &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='OPERATION_WINDOW_START',                        &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR50'
            
            call GetData(Me%TS_Operation(i)%OpWindowFileNameEnd,                        &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='OPERATION_WINDOW_END',                          &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR60'
        
            call GetData(Me%TS_Operation(i)%Max_Window,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='MAX_WINDOW',                                    &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR70'
            if (flag == 0) then
                stop 'Module_TS_Operator - Read_TS_OpWindow - ERR80'
            endif              

            call GetData(Me%TS_Operation(i)%Min_Window,                                 &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='MIN_WINDOW',                                    &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR100'
            if (flag == 0) then
                stop 'Module_TS_Operator - Read_TS_OpWindow - ERR110'
            endif  
            
            call GetData(Me%TS_Operation(i)%DT_Window,                                  &
                         Me%ObjEnterData,                                               &
                         flag,                                                          &
                         SearchType   = ExtractType,                                    &
                         keyword      ='DT_WINDOW',                                     &
                         ClientModule ='Module_TS_Operator',                            &
                         STAT         = STAT_CALL) 
            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - Read_TS_OpWindow - ERR120'
            if (flag == 0) then
                stop 'Module_TS_Operator - Read_TS_OpWindow - ERR130'
            endif              
            
        endif
        
    end subroutine Read_TS_OpWindow

    !--------------------------------------------------------------------------
    

    subroutine DefineOperatorBlock(ExtractType, block_begin, block_end, ExtractTypeOp)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        character(len = * )                         :: block_begin
        character(len = * )                         :: block_end
        integer                                     :: ExtractTypeOp
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------
        if  (ExtractType == FromFile_  ) then

            block_begin = "<begin_TS_operator>"
            block_end   = "<end_TS_operator>"
            
            ExtractTypeOp = FromBlock_

        elseif  (ExtractType == FromBlock_  ) then

            block_begin = "<<begin_TS_operator>>"
            block_end   = "<<end_TS_operator>>"
            
            ExtractTypeOp = FromBlockInBlock_

        elseif  (ExtractType == FromBlockInBlock_  ) then

            block_begin = "<<<begin_TS_operator>>>"
            block_end   = "<<<end_TS_operator>>>"
            
            ExtractTypeOp = FromBlockInBlockInBlock_
            
        endif
    
    end subroutine DefineOperatorBlock
    
    !-------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------

    subroutine DefineExpressionBlock(ExtractType, block_begin, block_end, ExtractTypeExp)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        character(len = * )                         :: block_begin
        character(len = * )                         :: block_end
        integer                                     :: ExtractTypeExp
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        if  (ExtractType == FromFile_  ) then

            block_begin = "<<begin_TS_expressions>>"
            block_end   = "<<end_TS_expressions>>"
            
            ExtractTypeExp = FromBlockInBlock_
            

        elseif  (ExtractType == FromBlock_  ) then

            block_begin = "<<<begin_TS_expressions>>>"
            block_end   = "<<<end_TS_expressions>>>"
            
            ExtractTypeExp = FromBlockInBlockInBlock_

        elseif  (ExtractType == FromBlockInBlock_  ) then

            stop 'Module_TS_Operator - DefineExpressionBlock - ERR10'
            
        endif
    
    end subroutine DefineExpressionBlock
    
    !-------------------------------------------------------------------------
    
!--------------------------------------------------------------------------

    subroutine DefineInputTS_Block(ExtractType, block_begin, block_end, ExtractTypeInput)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType
        character(len = * )                         :: block_begin
        character(len = * )                         :: block_end
        integer                                     :: ExtractTypeInput
        !Local-----------------------------------------------------------------

        !Begin-----------------------------------------------------------------

        if  (ExtractType == FromFile_  ) then

            block_begin = "<<begin_TS_input>>"
            block_end   = "<<end_TS_input>>"
            
            ExtractTypeInput = FromBlockInBlock_
            

        elseif  (ExtractType == FromBlock_  ) then

            block_begin = "<<<begin_TS_input>>>"
            block_end   = "<<<end_TS_input>>>"
            
            ExtractTypeInput = FromBlockInBlockInBlock_

        elseif  (ExtractType == FromBlockInBlock_  ) then

            stop 'Module_TS_Operator - DefineExpressionBlock - ERR10'
            
        endif
    
    end subroutine DefineInputTS_Block
    
    !-------------------------------------------------------------------------    
    subroutine Read_Operator_Expressions(ExtractType, i)    
    

        !Arguments----------------------------------------------------------
        integer                 :: ExtractType, i
        !Local--------------------------------------------------------------
        integer                 :: flag, STAT_CALL, j, FirstLine, LastLine, jaux
        logical                 :: BlockFound
        
        !Begin--------------------------------------------------------------
        

        call ReadBlock (ExtractType     = ExtractType,                                  &
                        ClientNumber    = Me%ClientNumber,                              &
                        block_begin     = Me%BlockBeginExpression,                      &
                        block_end       = Me%BlockEndExpression,                        &
                        BlockFound      = BlockFound,                                   &
                        FirstLine       = FirstLine,                                    &
                        LastLine        = LastLine)
        
        if (BlockFound) then
            
            jaux = LastLine - FirstLine - 1
            
            Me%TS_Operation(i)%NExpressions = jaux
            
            allocate(Me%TS_Operation(i)%Expression(jaux))
            
            if (jaux > 1) then
                
                allocate(Me%TS_Operation(i)%AuxTS(1:jaux-1, 1: Me%TS_Operation(i)%NValues))
            
            endif

            do j = 1, Me%TS_Operation(i)%NExpressions

                call GetData(Me%TS_Operation(i)%Expression(j),                          &
                             Me%ObjEnterData,   flag,                                   &
                             Buffer_Line  = FirstLine + j,                              &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - Read_Operator_Expressions - ERR10'
                endif

            enddo            
            
        else
            
            stop 'Module_TS_Operator - Read_Operator_Expressions - ERR20'
            
        endif
    


    end subroutine Read_Operator_Expressions

    !-------------------------------------------------------------------------
    
    subroutine Read_Operator_InputTS(ExtractType, i)    
    

        !Arguments----------------------------------------------------------
        integer                 :: ExtractType, i
        !Local--------------------------------------------------------------
        integer                 :: STAT_CALL, j
        integer                 :: NumberOfBlocks
        logical                 :: BlockFound
        
        !Begin--------------------------------------------------------------
        
        
        call GetNumberOfBlocks (EnterDataID     = Me%ObjEnterData,                      & 
                                BlockBegin      = Me%BlockBeginInputTS,                 &
                                BlockEnd        = Me%BlockEndInputTS,                   &
                                SearchType      = Me%ExtractTypeOp,                     &
                                NumberOfBlocks  = NumberOfBlocks,                       &
                                ClientNumber    = Me%ClientNumber,                      &
                                STAT            = STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) then
            stop 'Module_TS_Operator - Read_Operator_InputTS - ERR10'
        endif        
        
        Me%TS_Operation(i)%N_TS_Input = NumberOfBlocks
        
        allocate(Me%TS_Operation(i)%TS_Input(1:NumberOfBlocks))
        
        do j = 1, Me%TS_Operation(i)%N_TS_Input 

            call ReadBlock (ExtractType     = ExtractType,                              &
                            ClientNumber    = Me%ClientNumber,                          &
                            block_begin     = Me%BlockBeginInputTS,                     &
                            block_end       = Me%BlockEndInputTS,                       &
                            BlockFound      = BlockFound)
        
            if (BlockFound) then
            
                call Construct_TS_Synch(Obj_TS_SynchID      = Me%TS_Operation(i)%TS_Input(j)%Obj_TS_Synch,  & 
                                        EnterDataID         = Me%ObjEnterData,                              &
                                        ExtractType         = ExtractType,                                  &
                                        BeginTime           = Me%TS_Operation(i)%BeginTime,                 &
                                        EndTime             = Me%TS_Operation(i)%EndTime  ,                 &
                                        DT_Synch            = Me%TS_Operation(i)%DT_Synch ,                 &
                                        FillValue           = Me%TS_Operation(i)%FillValue,                 &
                                        STAT                = STAT_CALL)

                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - Read_Operator_InputTS - ERR10'
                endif        
                
                call Get_TS_SyncName(Obj_TS_SynchID = Me%TS_Operation(i)%TS_Input(j)%Obj_TS_Synch,      & 
                                     Name           = Me%TS_Operation(i)%TS_Input(j)%Name,              &
                                     STAT           = STAT_CALL)
                
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - Read_Operator_InputTS - ERR20'
                endif        
                
            else
            
                stop 'Module_TS_Operator - Read_Operator_InputTS - ERR20'
            
            endif
    
        enddo

    end subroutine Read_Operator_InputTS

    !-------------------------------------------------------------------------    

    subroutine ConstructOperatorTS
    
        !Local----------------------------------------------------------------
        integer                         :: STAT_CALL, i

        !Begin----------------------------------------------------------------
        
        do i = 1, Me%N_Operations
            
            !Open Output files
            call UnitsManager(Me%TS_Operation(i)%iOut, FileOpen, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'Module_TS_Operator - ConstructOperatorTS - ERR10'
            endif
        
            open(unit = Me%TS_Operation(i)%iOut, file = trim(Me%TS_Operation(i)%OutputFileName), &
                 form = 'FORMATTED', status = 'UNKNOWN')      
            
            
            if (Me%TS_Operation(i)%Out_V2_On) then

                !Open Output files
                call UnitsManager(Me%TS_Operation(i)%iOut_v2, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - ConstructOperatorTS - ERR20'
                endif
        
                open(unit = Me%TS_Operation(i)%iOut_v2, file = trim(Me%TS_Operation(i)%OutputFileNameV2), &
                     form = 'FORMATTED', status = 'UNKNOWN')                 
            endif
            
        enddo    
    
    end subroutine ConstructOperatorTS
    
    !-------------------------------------------------------------------------    
    


    subroutine ConstructOpWindowTS
    
        !Local----------------------------------------------------------------
        integer                         :: STAT_CALL, i

        !Begin----------------------------------------------------------------
        
        do i = 1, Me%N_Operations
            
            if (Me%TS_Operation(i)%Operation_Window) then
            
                !Open Output files
                call UnitsManager(Me%TS_Operation(i)%iW, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - ConstructOpWindowTS - ERR10'
                endif
        
                open(unit = Me%TS_Operation(i)%iW, file = trim(Me%TS_Operation(i)%OpWindowFileName), &
                     form = 'FORMATTED', status = 'UNKNOWN')      
            
            
                !Open Output files
                call UnitsManager(Me%TS_Operation(i)%iW_v2, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - ConstructOpWindowTS - ERR20'
                endif
        
                open(unit = Me%TS_Operation(i)%iW_v2, file = trim(Me%TS_Operation(i)%OpWindowFileNameV2), &
                        form = 'FORMATTED', status = 'UNKNOWN')                 
                    
                call UnitsManager(Me%TS_Operation(i)%iW_v3, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - ConstructOpWindowTS - ERR30'
                endif
        
                open(unit = Me%TS_Operation(i)%iW_v3, file = trim(Me%TS_Operation(i)%OpWindowFileNameStart), &
                        form = 'FORMATTED', status = 'UNKNOWN')                 
                    
                call UnitsManager(Me%TS_Operation(i)%iW_v4, FileOpen, STAT = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) then
                    stop 'Module_TS_Operator - ConstructOpWindowTS - ERR40'
                endif
        
                open(unit = Me%TS_Operation(i)%iW_v4, file = trim(Me%TS_Operation(i)%OpWindowFileNameEnd), &
                        form = 'FORMATTED', status = 'UNKNOWN')                 
                    
                
            endif
            
        enddo    
    
    end subroutine ConstructOpWindowTS
    
    !-------------------------------------------------------------------------        
    
    subroutine ReadBlock(ExtractType, ClientNumber, block_begin, block_end, BlockFound, &
                         FirstLine, LastLine)

        !Arguments-------------------------------------------------------------
        integer                                     :: ExtractType, ClientNumber
        character(len=*)                            :: block_begin, block_end
        logical                                     :: BlockFound
        integer, optional                           :: FirstLine, LastLine

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Begin-----------------------------------------------------------------


        if      (ExtractType == FromBlock_   ) then

            call ExtractBlockFromBuffer(EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = ClientNumber,             &
                                        block_begin         = block_begin,              &
                                        block_end           = block_end,                &
                                        BlockFound          = BlockFound,               &
                                        FirstLine           = FirstLine,                &
                                        LastLine            = LastLine,                 &
                                        STAT                = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - ReadBlock - ERR10'
            if (.not. BlockFound     ) then
                stop 'Module_TS_Operator - ReadBlock - ERR20'
            endif


        elseif  (ExtractType == FromBlockInBlock_  ) then

            call ExtractBlockFromBlock (EnterDataID         = Me%ObjEnterData,          &
                                        ClientNumber        = ClientNumber,             &
                                        block_begin         = block_begin,              &
                                        block_end           = block_end,                &
                                        BlockInBlockFound   = BlockFound,               &
                                        FirstLine           = FirstLine,                &
                                        LastLine            = LastLine,                 &
                                        STAT                = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - ReadBlock - ERR30'
            if (.not. BlockFound     ) stop 'Module_TS_Operator - ReadBlock - ERR40'

        elseif  (ExtractType == FromBlockInBlockInBlock_  ) then

            call ExtractBlockFromBlockFromBlock(EnterDataID              = Me%ObjEnterData, &
                                                ClientNumber             = ClientNumber,    &
                                                block_begin              = block_begin,     &
                                                block_end                = block_end,       &
                                                BlockInBlockInBlockFound = BlockFound,      &
                                                FirstLine                = FirstLine,       &
                                                LastLine                 = LastLine,        &
                                                STAT                     = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'Module_TS_Operator - ReadBlock - ERR50'
            if (.not. BlockFound     ) stop 'Module_TS_Operator - ReadBlock - ERR60'


        else

            stop 'Module_TS_Operator - ReadBlock - ERR70'
            
        endif
        
    end subroutine ReadBlock        
    
    !-------------------------------------------------------------------------    
    
    !
    !character(len=PathLength) function AddString2FileName(Filename, AddString)
    !
    !    !Arguments------------------------------------------------------------    
    !    character(len=*)                :: Filename, AddString
    !    
    !    !Local----------------------------------------------------------------
    !    integer                         :: n, i, k
    !
    !    !Begin----------------------------------------------------------------    
    !    
    !    n = len_trim(Filename)
    !    
    !    k = FillValueInt
    !    
    !    do i=n,1,-1
    !        if (Filename(i:i) == "/" .or. Filename(i:i) == "\") then
    !            k = i
    !            exit
    !        endif                
    !    enddo
    !    
    !    if (k > FillValueInt) then
    !        AddString2FileName = Filename(1:k)//trim(AddString)//Filename(k+1:n)
    !    else
    !        AddString2FileName = trim(AddString)//trim(Filename)
    !    endif            
    !    
    !
    !end function AddString2FileName

!-------------------------------------------------------------------------    

    


    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_TS_Operator), pointer                         :: NewObj_TS_Operator
        type (T_TS_Operator), pointer                         :: PreviousObj_TS_Operator


        !Allocates new instance
        allocate (NewObj_TS_Operator)
        nullify  (NewObj_TS_Operator%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObj_TS_Operator)) then
            FirstObj_TS_Operator         => NewObj_TS_Operator
            Me                    => NewObj_TS_Operator
        else
            PreviousObj_TS_Operator      => FirstObj_TS_Operator
            Me                    => FirstObj_TS_Operator%Next
            do while (associated(Me))
                PreviousObj_TS_Operator  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObj_TS_Operator
            PreviousObj_TS_Operator%Next => NewObj_TS_Operator
        endif

        Me%InstanceID = RegisterNewInstance (m_TS_Operator_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine Get_TS_Operator_Values (Obj_TS_OperatorID, i, ValuesTS, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_OperatorID
        integer                                         :: i    
        real(8), dimension(:),  pointer                 :: ValuesTS
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_OperatorID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(m_TS_Operator_, Me%InstanceID)

            ValuesTS => Me%TS_Operation(i)%TimeSerie

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_Operator_Values
    
    !--------------------------------------------------------------------------
    
    subroutine Get_TS_OperatorID (Obj_TS_OperatorID, ID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_OperatorID
        real                                            :: ID
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_OperatorID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            ID = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Get_TS_OperatorID

    !--------------------------------------------------------------------------


    subroutine UnGet_TS_Operator1D_R8(Obj_TS_OperatorID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: Obj_TS_OperatorID
        real(8), dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_OperatorID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(m_TS_Operator_, Me%InstanceID,  "UnGet_TS_Operator3D_R8")


            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGet_TS_Operator1D_R8

    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine Modify_TS_Operator(Obj_TS_OperatorID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: Obj_TS_OperatorID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_, i

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_OperatorID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            
            do i = 1, Me%N_Operations
                call ModifyTimeSeriesOperation(i)
            enddo     
                        
            
            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine Modify_TS_Operator
    
    !--------------------------------------------------------------------------
    
    subroutine ModifyTimeSeriesOperation(i)

        !Arguments-------------------------------------------------------------
        integer                                     :: i
        !Local-----------------------------------------------------------------
        real                                        :: Year, Month, Day, hour, minute, second
        real(8)                                     :: NewValue, Tx
        integer                                     :: j, NValues, n, e, e1, STAT_CALL
        integer                                     :: Nvar
        logical                                     :: FillValue
        character(len = 1000)                       :: func, Expression
        character(len = 5)                          :: statusflag
        character(len = 10)                         :: Name
        character(len = 10),  dimension(:), allocatable :: variables
        real(8),       dimension(:,:), allocatable :: variablesvalues 
        logical                                     :: AssumeFillValue, Conditional
        integer, dimension(5)                       :: ConditionTS
        real(8), dimension(5)                       :: ConditionValues
        
        
        !----------------------------------------------------------------------
        
        NVar =  Me%TS_Operation(i)%N_TS_Input + Me%TS_Operation(i)%NExpressions
        
        allocate(variables(1:Nvar))
        allocate(variablesvalues(1:Nvar,1:Me%TS_Operation(i)%NValues))
        
        !Synchronise all input times series
        
        do j = 1, Me%TS_Operation(i)%N_TS_Input 

            call Modify_TS_Synch(Obj_TS_SynchID = Me%TS_Operation(i)%TS_Input(j)%Obj_TS_Synch,& 
                                 STAT           = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                stop 'Module_TS_Operator - ModifyTimeSeriesOperation - ERR10'
            endif        
                                                                                                    

            call Get_TS_Synch_nValues(Obj_TS_SynchID = Me%TS_Operation(i)%TS_Input(j)%Obj_TS_Synch, & 
                                      nValues        = nValues,                                     &
                                      STAT           = STAT_CALL)

            if (STAT_CALL /= SUCCESS_) then
                stop 'Module_TS_Operator - ModifyTimeSeriesOperation - ERR30'
            endif 
                
            if (nValues /= Me%TS_Operation(i)%NValues) then
                stop 'Module_TS_Operator - ModifyTimeSeriesOperation - ERR40'
            endif
            
            variables(j) = trim(Me%TS_Operation(i)%TS_Input(j)%Name)
            
        enddo
        
        do e =1, Me%TS_Operation(i)%NExpressions
                    
            Expression = Me%TS_Operation(i)%Expression(e)
                
            call GetExpression(Expression, Name, Func, Conditional)
                
            e1 = Me%TS_Operation(i)%N_TS_Input + e
            
            if (Conditional) then
                
                call init_condition (func, variables(1:e1-1), ConditionTS, ConditionValues, statusflag)

            else
                    
                call init (func, variables(1:e1-1), statusflag)

            endif
            
            if (statusflag /= "ok") then
                write(*,*) "Something wrong with the follow expression ", trim(Expression)
            endif                        
                
            variables(e1) = Name
            
            if (e1 == Nvar) then
            
                Me%TS_Operation(i)%Name         = Name
                
            endif
            
        
            do n =1, Me%TS_Operation(i)%NValues
            
                 Me%TS_Operation(i)%TimeSerie(n) = Me%TS_Operation(i)%FillValue
             
             
                do j = 1, Me%TS_Operation(i)%N_TS_Input 
                
                    call Get_TS_Synch_ValueX(Obj_TS_SynchID = Me%TS_Operation(i)%TS_Input(j)%Obj_TS_Synch,  & 
                                            n               = n,                                            &
                                            ValueX          = Tx,                                           &
                                            FillValue       = FillValue,                                    &
                                            STAT            = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) then
                        stop 'Module_TS_Operator - ModifyTimeSeriesOperation - ERR20'
                    endif  
                

                    if (FillValue) then
                        AssumeFillValue     = .true.
                        exit
                    else
                        AssumeFillValue     = .false.
                        variablesvalues(j,n)  = Tx
                    endif
                
                enddo
            
                if (AssumeFillValue) then
                    NewValue = Me%TS_Operation(i)%FillValue
            
                else
                    
                    if (Conditional) then
                        variablesvalues(e1,n) = compute_condition(ConditionTS, ConditionValues, variablesvalues(1:e1-1,n))
                    else                    
                        variablesvalues(e1,n) = evaluate (variablesvalues(1:e1-1,n))
                    endif
                    

                    NewValue = variablesvalues(e1,n)
                
                endif
            
                Me%TS_Operation(i)%TimeSerie(n) =  NewValue
            

            enddo     
            
            if (.not. Conditional) then
                call destroyfunc()
            endif

        
        enddo
        
       
        write(Me%TS_Operation(i)%iOut,'(A25,A6)') "NAME                    : ", trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iOut,*         ) "LOCALIZATION_I          : -999999"
        write(Me%TS_Operation(i)%iOut,*         ) "LOCALIZATION_J          : -999999"
        write(Me%TS_Operation(i)%iOut,*         ) "LOCALIZATION_K          : -999999"
        
        call ExtractDate(Me%TS_Operation(i)%BeginTime, Year, Month, Day, hour, minute, second)
        
        write(Me%TS_Operation(i)%iOut,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%TS_Operation(i)%iOut,*) "TIME_UNITS              : SECONDS"
        write(Me%TS_Operation(i)%iOut,*) "COORD_X                 : ", Me%TS_Operation(i)%CoordX
        write(Me%TS_Operation(i)%iOut,*) "COORD_Y                 : ", Me%TS_Operation(i)%CoordY
        write(Me%TS_Operation(i)%iOut,*) "FILL_VALUE              : ", Me%TS_Operation(i)%FillValue
        
        write(Me%TS_Operation(i)%iOut,'(A5,A6)') "Time ", variables(e1)
        
        write(Me%TS_Operation(i)%iOut,*) "<BeginTimeSerie>"         
        
        do n=1, Me%TS_Operation(i)%NValues        
            write(Me%TS_Operation(i)%iOut,'(f16.1,e12.4)') real(n-1)*Me%TS_Operation(i)%DT_Synch, variablesvalues(e1,n)
        enddo
        
        write(Me%TS_Operation(i)%iOut,*) "<EndTimeSerie>"      
        

        !Open Output files
        call UnitsManager(Me%TS_Operation(i)%iOut, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Operator - ModifyTimeSeriesOperation - ERR40"
        
        if (Me%TS_Operation(i)%Out_V2_On) then        
        
            write(Me%TS_Operation(i)%iOut_v2,'(A25,A6)') "NAME                    : ", trim(Me%TS_Operation(i)%Name)
            write(Me%TS_Operation(i)%iOut_v2,*         ) "LOCALIZATION_I          : -999999"
            write(Me%TS_Operation(i)%iOut_v2,*         ) "LOCALIZATION_J          : -999999"
            write(Me%TS_Operation(i)%iOut_v2,*         ) "LOCALIZATION_K          : -999999"
        
            call ExtractDate(Me%TS_Operation(i)%BeginTime, Year, Month, Day, hour, minute, second)
        
            write(Me%TS_Operation(i)%iOut_v2,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
            write(Me%TS_Operation(i)%iOut_v2,*) "TIME_UNITS              : SECONDS"
            write(Me%TS_Operation(i)%iOut_v2,*) "COORD_X                 : ", Me%TS_Operation(i)%CoordX
            write(Me%TS_Operation(i)%iOut_v2,*) "COORD_Y                 : ", Me%TS_Operation(i)%CoordY
            write(Me%TS_Operation(i)%iOut_v2,*) "FILL_VALUE              : ", Me%TS_Operation(i)%FillValue
        
            write(Me%TS_Operation(i)%iOut_v2,'(A5,100A6)') "Time ", variables
        
            write(Me%TS_Operation(i)%iOut_v2,*) "<BeginTimeSerie>"         
        
            do n=1, Me%TS_Operation(i)%NValues        
                write(Me%TS_Operation(i)%iOut_v2,'(f16.1,100e12.4)') real(n-1)*Me%TS_Operation(i)%DT_Synch, variablesvalues(1:e1,n)
            enddo
        
            write(Me%TS_Operation(i)%iOut_v2,*) "<EndTimeSerie>"      
        

            !Open Output files
            call UnitsManager(Me%TS_Operation(i)%iOut_v2, FileClose, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop "Module_TS_Operator - ModifyTimeSeriesOperation - ERR40"            

        endif
        
        if (Me%TS_Operation(i)%Operation_Window) then
            call Modify_Operation_Window(i, variablesvalues(e1,:))
        endif
        
        
        deallocate(variables)
        deallocate(variablesvalues)        
        
    end subroutine ModifyTimeSeriesOperation

    !--------------------------------------------------------------------------
    
    subroutine Modify_Operation_Window(i, TS_Values)

        !Arguments-------------------------------------------------------------
        integer                                     :: i
        real(8),       dimension(:)                 :: TS_Values         
        !Local-----------------------------------------------------------------
        real                                        :: DT_Aux, DT_Synch
        real                                        :: Year, Month, Day, hour, minute, second
        integer                                     :: Nvar, n, STAT_CALL, is, aux, NValues
        integer,       dimension(:), allocatable    :: WindowTS, AuxTS
        logical                                     :: StartC
        
        !----------------------------------------------------------------------
        
        Nvar = size(TS_Values)
        
        allocate(WindowTS(1:Nvar), AuxTS(1:Nvar))
        
        WindowTS(:) = 0
        
        do n =1, Nvar
            
            if (TS_Values(n) >= Me%TS_Operation(i)%Min_Window .and. TS_Values(n) <= Me%TS_Operation(i)%Max_Window) then
                AuxTS(n) = 1
            else
                AuxTS(n) = 0
            endif

        enddo
        
        DT_Aux = 0
        StartC = .true.
        do n =1, Nvar
            
            if (AuxTS(n) == 1) then
                if (StartC) then
                    is = n
                    StartC = .false.
                endif
                DT_Aux = DT_Aux + Me%TS_Operation(i)%DT_Synch
                if (DT_Aux >= Me%TS_Operation(i)%DT_Window) then
                    WindowTS(is:n) = 1   
                endif
            else
                DT_Aux = 0.  
                StartC = .true. 
            endif

        enddo        
       
       
        write(Me%TS_Operation(i)%iW,'(A25,A50)') "NAME                    : ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW,*          ) "LOCALIZATION_I          : -999999"
        write(Me%TS_Operation(i)%iW,*          ) "LOCALIZATION_J          : -999999"
        write(Me%TS_Operation(i)%iW,*          ) "LOCALIZATION_K          : -999999"
        
        call ExtractDate(Me%TS_Operation(i)%BeginTime, Year, Month, Day, hour, minute, second)
        
        write(Me%TS_Operation(i)%iW,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%TS_Operation(i)%iW,*) "TIME_UNITS              : SECONDS"
        write(Me%TS_Operation(i)%iW,*) "COORD_X                 : ", Me%TS_Operation(i)%CoordX
        write(Me%TS_Operation(i)%iW,*) "COORD_Y                 : ", Me%TS_Operation(i)%CoordY
        write(Me%TS_Operation(i)%iW,*) "FILL_VALUE              : -99"
        
        write(Me%TS_Operation(i)%iW,'(A5,A60)') "Time ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        
        write(Me%TS_Operation(i)%iW,*) "<BeginTimeSerie>"         
        

        do n=1, Me%TS_Operation(i)%NValues
            write(Me%TS_Operation(i)%iW,'(f16.1,i4)') real(n-1)*Me%TS_Operation(i)%DT_Synch, WindowTS(n)
        enddo
                
        
        write(Me%TS_Operation(i)%iW,*) "<EndTimeSerie>"      
        

        !Open Output files
        call UnitsManager(Me%TS_Operation(i)%iW, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Operator - ModifyTimeSeriesOperation - ERR40"

        call ExtractDate(Me%TS_Operation(i)%BeginTime, Year, Month, Day, hour, minute, second)
           
        write(Me%TS_Operation(i)%iW_v2,'(A25,A50)') "NAME                    : ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW_v2,*          ) "LOCALIZATION_I          : -999999"
        write(Me%TS_Operation(i)%iW_v2,*          ) "LOCALIZATION_J          : -999999"
        write(Me%TS_Operation(i)%iW_v2,*          ) "LOCALIZATION_K          : -999999"
        write(Me%TS_Operation(i)%iW_v2,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%TS_Operation(i)%iW_v2,*) "TIME_UNITS              : SECONDS"
        write(Me%TS_Operation(i)%iW_v2,*) "COORD_X                 : ", Me%TS_Operation(i)%CoordX
        write(Me%TS_Operation(i)%iW_v2,*) "COORD_Y                 : ", Me%TS_Operation(i)%CoordY
        write(Me%TS_Operation(i)%iW_v2,*) "FILL_VALUE              : -99"
        write(Me%TS_Operation(i)%iW_v2,'(A5,A60)') "Time ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW_v2,*) "<BeginTimeSerie>"         
        
        write(Me%TS_Operation(i)%iW_v3,'(A25,A50)') "NAME                    : ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW_v3,*          ) "LOCALIZATION_I          : -999999"
        write(Me%TS_Operation(i)%iW_v3,*          ) "LOCALIZATION_J          : -999999"
        write(Me%TS_Operation(i)%iW_v3,*          ) "LOCALIZATION_K          : -999999"
        write(Me%TS_Operation(i)%iW_v3,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%TS_Operation(i)%iW_v3,*) "TIME_UNITS              : SECONDS"
        write(Me%TS_Operation(i)%iW_v3,*) "COORD_X                 : ", Me%TS_Operation(i)%CoordX
        write(Me%TS_Operation(i)%iW_v3,*) "COORD_Y                 : ", Me%TS_Operation(i)%CoordY
        write(Me%TS_Operation(i)%iW_v3,*) "FILL_VALUE              : -99"
        write(Me%TS_Operation(i)%iW_v3,'(A5,A60)') "Time ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW_v3,*) "<BeginTimeSerie>"         
        
        write(Me%TS_Operation(i)%iW_v4,'(A25,A50)') "NAME                    : ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW_v4,*          ) "LOCALIZATION_I          : -999999"
        write(Me%TS_Operation(i)%iW_v4,*          ) "LOCALIZATION_J          : -999999"
        write(Me%TS_Operation(i)%iW_v4,*          ) "LOCALIZATION_K          : -999999"
        write(Me%TS_Operation(i)%iW_v4,'(A26,5F6.0,1f8.2)') "SERIE_INITIAL_DATA      : ", Year, Month, Day, hour, minute, second
        write(Me%TS_Operation(i)%iW_v4,*) "TIME_UNITS              : SECONDS"
        write(Me%TS_Operation(i)%iW_v4,*) "COORD_X                 : ", Me%TS_Operation(i)%CoordX
        write(Me%TS_Operation(i)%iW_v4,*) "COORD_Y                 : ", Me%TS_Operation(i)%CoordY
        write(Me%TS_Operation(i)%iW_v4,*) "FILL_VALUE              : -99"
        write(Me%TS_Operation(i)%iW_v4,'(A5,A60)') "Time ", "Operation_Window_"//trim(Me%TS_Operation(i)%Name)
        write(Me%TS_Operation(i)%iW_v4,*) "<BeginTimeSerie>"         
        
        
            
        if (sum(WindowTS) == 0) then
            write(Me%TS_Operation(i)%iW_v2,'(f16.1,i4)') 0., -99
            write(Me%TS_Operation(i)%iW_v2,'(f16.1,i4)') real(Nvar-1)*Me%TS_Operation(i)%DT_Synch, -99
            
            write(Me%TS_Operation(i)%iW_v3,'(f16.1,i4)') 0., -99
            write(Me%TS_Operation(i)%iW_v3,'(f16.1,i4)') real(Nvar-1)*Me%TS_Operation(i)%DT_Synch, -99
            
            write(Me%TS_Operation(i)%iW_v4,'(f16.1,i4)') 0., -99
            write(Me%TS_Operation(i)%iW_v4,'(f16.1,i4)') real(Nvar-1)*Me%TS_Operation(i)%DT_Synch, -99
                
        else
            aux      = -99
            NValues  = Me%TS_Operation(i)%NValues
            DT_Synch = Me%TS_Operation(i)%DT_Synch 
                
            do n=1,NValues-1
                if (WindowTS(n+1) == 1 .and. WindowTS(n) == 0) then
                    aux = 0
                    write(Me%TS_Operation(i)%iW_v2,'(f16.1,i4)') real(n)*DT_Synch, aux                        
                    write(Me%TS_Operation(i)%iW_v3,'(f16.1,i4)') real(n)*DT_Synch, aux                                            
                endif
                if (WindowTS(n+1) == 0 .and. WindowTS(n) == 1) then
                    aux = 1
                    write(Me%TS_Operation(i)%iW_v2,'(f16.1,i4)') real(n-1)*DT_Synch, aux
                    write(Me%TS_Operation(i)%iW_v4,'(f16.1,i4)') real(n-1)*DT_Synch, aux                    
                endif                    
            enddo
                
        endif
        
        write(Me%TS_Operation(i)%iW_v2,*) "<EndTimeSerie>"
        write(Me%TS_Operation(i)%iW_v3,*) "<EndTimeSerie>"
        write(Me%TS_Operation(i)%iW_v4,*) "<EndTimeSerie>"        
        
        !Open Output files
        call UnitsManager(Me%TS_Operation(i)%iW_v2, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Operator - ModifyTimeSeriesOperation - ERR40"            

        call UnitsManager(Me%TS_Operation(i)%iW_v3, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Operator - ModifyTimeSeriesOperation - ERR50"
        
        call UnitsManager(Me%TS_Operation(i)%iW_v4, FileClose, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop "Module_TS_Operator - ModifyTimeSeriesOperation - ERR60"        
        
        deallocate(WindowTS, AuxTS)
         
    end subroutine Modify_Operation_Window
    
    !--------------------------------------------------------------------------      
    
    !--------------------------------------------------------------------------    
    
    subroutine init_condition(Func, variables, ConditionTS, ConditionValues, statusflag)
    
        !Arguments-----------------------------------------------------------------
        character(len=*),               intent (IN ) :: Func
        character(len=*), dimension(:), intent (IN ) :: variables        
        integer,          dimension(:), intent (OUT) :: ConditionTS
        real(8),          dimension(:), intent (OUT) :: ConditionValues        
        character(len=*)              , intent (OUT) :: statusflag
    
        !Local-----------------------------------------------------------------    
        integer                                      :: l, p1, p2, v1, v2, v3, v4
        character(len=StringLength)                  :: AuxChar
        !Begin-----------------------------------------------------------------
        statusflag = "ok"
        
        l = len(Func)
        
        p1 = scan(Func,"(")
        if (p1 == 0) then
             statusflag = "err_1"
        endif

        p2 = scan(Func,")")
        if (p2 == 0) then
             statusflag = "err_2"
        endif 
        
        
        v1 = scan(Func,",")
        if (v1 == 0) then
             statusflag = "err_3"
        endif 

        v2 = v1 + scan(Func(v1+1:p2),",")
        if (v2 == 0) then
             statusflag = "err_4"
        endif 
        
        v3 = v2 + scan(Func(v2+1:p2),",")        
        if (v3 == 0) then
             statusflag = "err_5"
        endif 
        
        v4 = v3 + scan(Func(v3+1:p2),",")                
        if (v4 == 0) then
             statusflag = "err_6"
        endif     
        
        ConditionValues(:) = FillValueReal
        
        ConditionTS(1) = variableID (Func(p1+1:v1-1), variables)
        
        if (ConditionTS(1) == 0) then
            !stop "Module_TS_Operator - init_condition - ERR10"
            AuxChar = Func(p1+1:v1-1)
            read(AuxChar,*,err=10) ConditionValues(1)
        endif
        
        ConditionTS(2) = conditionID (Func(v1+1:v2-1))        
        if (ConditionTS(2) == 0) then
            stop "Module_TS_Operator - init_condition - ERR20"
        endif        
        
        ConditionTS(3) = variableID (Func(v2+1:v3-1), variables)
        if (ConditionTS(3) == 0) then
            !stop "Module_TS_Operator - init_condition - ERR30"
            AuxChar = Func(v2+1:v3-1)
            read(AuxChar,*,err=10) ConditionValues(3)            
        endif        
        
        ConditionTS(4) = variableID (Func(v3+1:v4-1), variables)
        if (ConditionTS(4) == 0) then
            !stop "Module_TS_Operator - init_condition - ERR40"
            AuxChar = Func(v3+1:v4-1)
            read(AuxChar,*,err=10) ConditionValues(4)                        
        endif
        
        ConditionTS(5) = variableID (Func(v4+1:p2-1), variables)       
        if (ConditionTS(5) == 0) then
            !stop "Module_TS_Operator - init_condition - ERR50"
            AuxChar = Func(v4+1:p2-1)
            read(AuxChar,*,err=10) ConditionValues(5)                                    
        endif   
        
        return 
10      write(*,*) 'Not valid variable or value =',  trim(AuxChar)
        write(*,*) 'see expression =',  trim(Func)
        stop 
    
    end subroutine init_condition
    
    !--------------------------------------------------------------------------    
    
    
    integer function variableID(varX, variables)
    
        !Arguments-----------------------------------------------------------------
        character(len=*),               intent (IN ) :: varX
        character(len=*), dimension(:), intent (IN ) :: variables        
    
        !Local-----------------------------------------------------------------    
        integer                 :: l, i
    
        !Begin-----------------------------------------------------------------
        l = size(variables)
        
        !By default not found
        variableID = 0
        
        do i=1, l
            if (trim(adjustl(varX)) == trim(adjustl(variables(i)))) then
                variableID = i
                exit
            endif
        enddo
        
    
    end function variableID
    
    !--------------------------------------------------------------------------     
    
    
    integer function conditionID(condX)
    
        !Arguments-----------------------------------------------------------------
        character(len=*),               intent (IN ) :: condX
    
        !Local-----------------------------------------------------------------    

    
        !Begin-----------------------------------------------------------------
        
        if      (trim(adjustl(condX)) == "<" ) then
            conditionID = 1
        elseif  (trim(adjustl(condX)) == "<=") then
            conditionID = 2
        elseif  (trim(adjustl(condX)) == "==") then
            conditionID = 3            
        elseif  (trim(adjustl(condX)) == "/=") then
            conditionID = 4            
        elseif  (trim(adjustl(condX)) == ">=") then
            conditionID = 5
        elseif  (trim(adjustl(condX)) == ">" ) then
            conditionID = 6           
        else
            !not found
            conditionID = 0
        endif
    
    end function conditionID
    
    !--------------------------------------------------------------------------         
    
    
    logical function conditionViaID(conditionID, a, b)
    
        !Arguments-----------------------------------------------------------------
        integer,                        intent (IN ) :: conditionID
        real(8),                        intent (IN ) :: a, b
        !Local-----------------------------------------------------------------    
    
        !Begin-----------------------------------------------------------------
        
        if      (conditionID == 1) then
            if (a < b) then
                conditionViaID = .true.
            else
                conditionViaID = .false.                
            endif
        elseif (conditionID == 2) then
            if (a <= b) then
                conditionViaID = .true.
            else
                conditionViaID = .false.                
            endif   
        elseif (conditionID == 3) then
            if (a == b) then
                conditionViaID = .true.
            else
                conditionViaID = .false.                
            endif   
        elseif (conditionID == 4) then
            if (a /= b) then
                conditionViaID = .true.
            else
                conditionViaID = .false.                
            endif   
        elseif (conditionID == 5) then
            if (a >= b) then
                conditionViaID = .true.
            else
                conditionViaID = .false.                
            endif   
        elseif (conditionID == 6) then
            if (a >  b) then
                conditionViaID = .true.
            else
                conditionViaID = .false.                
            endif  
        else
            write(*,*) "Logical operator no valid"
            stop "Module_TS_Operator - conditionViaID - ERR10"
        endif
    
    end function conditionViaID
    
    !--------------------------------------------------------------------------         
        
        
        
    
    real(8) function compute_condition(ConditionTS, ConditionValues, variablesvalue)
    
        !Arguments-----------------------------------------------------------------
        integer,          dimension(5), intent (IN ) :: ConditionTS
        real(8),          dimension(5), intent (IN ) :: ConditionValues        
        real(8),          dimension(:), intent (IN ) :: variablesvalue        
        
    
        !Local-----------------------------------------------------------------    
        logical                                     :: condition
        integer                                     :: IDa, IDb, IDc, IDd, IDcondition
        real(8)                                     :: a, b, c, d
    
        !Begin-----------------------------------------------------------------
        
        !if (a condition b) then c else d endif
        IDa         = ConditionTS(1)
        IDcondition = ConditionTS(2)        
        IDb         = ConditionTS(3)
        IDc         = ConditionTS(4)
        IDd         = ConditionTS(5)
        
        if (IDa == 0) then
            a = ConditionValues(1)
        else
            a = variablesvalue(IDa)
        endif
        
        if (IDb == 0) then
            b = ConditionValues(3)
        else
            b = variablesvalue(IDb)
        endif

        if (IDc == 0) then
            c = ConditionValues(4)
        else
            c = variablesvalue(IDc)
        endif
        
        if (IDd == 0) then
            d = ConditionValues(5)
        else
            d = variablesvalue(IDd)
        endif
                
        
        condition = conditionViaID(IDcondition, a, b)
        
        if (condition) then
            compute_condition = c
        else
            compute_condition = d
        endif

    
    end function compute_condition    
    
    !--------------------------------------------------------------------------        
    
    subroutine GetExpression(Expression, Name, Func, Conditional)
    
        !Arguments-----------------------------------------------------------------
        character(len=*), intent (IN)   :: Expression
        character(len=*), intent (OUT)  :: Name    
        character(len=*), intent (OUT)  :: Func       
        logical         , intent (OUT)  :: Conditional
    
        !Local-----------------------------------------------------------------    
        integer                         :: l, s, i
    
        !Begin-----------------------------------------------------------------
    
        l = scan(trim(Expression),'=')
        if (l==0) then
            write(*,*) "In expression ",trim(Expression) 
            write(*,*) "is missing the equal sign '=' "
            stop "Module_TS_Operator - GetExpression - ERR10"
        endif    
    
        s = len(Expression)    
    
        i = scan(trim(Expression),'if') 
    
        !conditional relation (if (a,<=,b,c,d)
        !a and b input and c true, d false
        if (i > 0) then
            !Func = (a,<=,b,c,d)
            Func = trim(Expression(i+2:s))        
            Conditional = .true.
        else !algebric function
            !Func = a*b + c*d
            Func = trim(Expression(l+1:s))        
            Conditional = .false.
        endif
    
        Name = trim(Expression(1:l-1))

    
    
    end subroutine GetExpression    
    
    !--------------------------------------------------------------------------        

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine Kill_TS_Operator(Obj_TS_OperatorID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: Obj_TS_OperatorID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Obj_TS_OperatorID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(m_TS_Operator_,  Me%InstanceID)

            if (nUsers == 0) then
            
                call KillVariablesAndFiles
            
                !Deallocates Instance
                call DeallocateInstance ()

                Obj_TS_OperatorID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine Kill_TS_Operator
        

    !------------------------------------------------------------------------
    


    subroutine KillVariablesAndFiles
    
    
        !Local--------------------------------------------------------------------------
        integer                         :: i
        !Begin--------------------------------------------------------------------------
    
    
        do i=1, Me%N_Operations
            deallocate(Me%TS_Operation(i)%TimeSerie)
        enddo
        
        deallocate(Me%TS_Operation)
    
    end subroutine KillVariablesAndFiles    
    
    !--------------------------------------------------------------------------    
    
   !--------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_TS_Operator), pointer          :: AuxObj_TS_Operator
        type (T_TS_Operator), pointer          :: PreviousObj_TS_Operator

        !Updates pointers
        if (Me%InstanceID == FirstObj_TS_Operator%InstanceID) then
            FirstObj_TS_Operator => FirstObj_TS_Operator%Next
        else
            PreviousObj_TS_Operator => FirstObj_TS_Operator
            AuxObj_TS_Operator      => FirstObj_TS_Operator%Next
            do while (AuxObj_TS_Operator%InstanceID /= Me%InstanceID)
                PreviousObj_TS_Operator => AuxObj_TS_Operator
                AuxObj_TS_Operator      => AuxObj_TS_Operator%Next
            enddo

            !Now update linked list
            PreviousObj_TS_Operator%Next => AuxObj_TS_Operator%Next

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

    subroutine Ready (Obj_TS_Operator_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: Obj_TS_Operator_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (Obj_TS_Operator_ID > 0) then
            call LocateObj_TS_Operator (Obj_TS_Operator_ID)
            ready_ = VerifyReadLock (m_TS_Operator_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObj_TS_Operator (Obj_TS_OperatorID)

        !Arguments-------------------------------------------------------------
        integer                                     :: Obj_TS_OperatorID

        !Local-----------------------------------------------------------------

        Me => FirstObj_TS_Operator
        do while (associated (Me))
            if (Me%InstanceID == Obj_TS_OperatorID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'Module_TS_Operator - LocateObj_TS_Operator - ERR01'

    end subroutine LocateObj_TS_Operator

    !--------------------------------------------------------------------------

    end module Module_TS_Operator
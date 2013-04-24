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

Module ModuleNetworkStatistics

    use HDF5
    use ModuleGlobalData
    use ModuleFunctions
    use ModuleEnterData
    use ModuleTime
    use ModuleHDF5
    !use nr; use nrtype; use nrutil

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConstructNetworkStatistics
    private ::      AllocateInstance

    !Selector
    public  :: GetNetworkStatisticsPointer
    public  :: GetNetworkStatisticsInteger
    public  :: UnGetNetworkStatistics
                     
    
    !Modifier
    public  :: ModifyNetworkStatistics

    !Destructor
    public  :: KillNetworkStatistics                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjNetworkStatistics 
    
    !Interfaces----------------------------------------------------------------
    private :: UnGetNetworkStatistics3D_I
    private :: UnGetNetworkStatistics3D_R8
    interface  UnGetNetworkStatistics
        module procedure UnGetNetworkStatistics3D_I
        module procedure UnGetNetworkStatistics3D_R8
    end interface  UnGetNetworkStatistics


    !Parameter-----------------------------------------------------------------
    
    real(8), parameter  :: Pi_ = 3.1415926535897932384626433832795
    !Input / Output
    integer, parameter  :: FileOpen = 1, FileClose = 0
    

    integer, parameter  :: Day_ = 1, Week_ = 2, Month_ = 3    
    
    integer, parameter  :: WeekEnd_ = 1, WeekWay_ = 2
    


    
    !Types---------------------------------------------------------------------
    type T_NetworkStatistics
    
        integer                                                 :: InstanceID
    
        character(len=PathLength)                               :: InputFile, OutputFile
        
        type(T_Time)                                            :: StartTime, EndTime
        type(T_Time), dimension(:), pointer                     :: InputTime, OutputTime
        
        integer                                                 :: NInstants, NGroups, NCopyGroups
        character(len=StringLength), dimension(:), pointer      :: AnalysisGroups, AnalysisDataSets
        character(len=StringLength), dimension(:), pointer      :: CopyGroups

        real, dimension(:), pointer                             :: ValuesToAnalyse
        integer                                                 :: DTAnalysis
        real, dimension(:,:), pointer                           :: PeakDemandStart, PeakDemandEnd
        real, dimension(:,:), pointer                           :: LowDemandStart,  LowDemandEnd

	    integer                                                 :: ObjEnterData          = 0
	    integer                                                 :: ObjTime               = 0
	    integer                                                 :: ObjHDF5_In            = 0
	    integer                                                 :: ObjHDF5_Out           = 0	    
	    
        type(T_NetworkStatistics), pointer                     :: Next

    end type T_NetworkStatistics    

    !Global Variables
    type (T_NetworkStatistics), pointer                        :: FirstObjNetworkStatistics
    type (T_NetworkStatistics), pointer                        :: Me    


    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConstructNetworkStatistics(ObjNetworkStatisticsID, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjNetworkStatisticsID 
        integer, optional, intent(OUT)                  :: STAT     

        !External----------------------------------------------------------------
        integer                                         :: ready_         

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_, STAT_CALL

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mNetworkStatistics_)) then
            nullify (FirstObjNetworkStatistics)
            call RegisterModule (mNetworkStatistics_) 
        endif

        call Ready(ObjNetworkStatisticsID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance

            call ReadKeywords
            
            call OpenHDF5Files
            
            call ConstructDemandPatterns            
            
            call ConstructInputTime
            
            call ConstructVGroupNames
            
            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL)
            if(STAT_CALL /= SUCCESS_) stop 'ModuleNetworkStatistics - ReadKeywords - ERR330'            
            
            !Returns ID
            ObjNetworkStatisticsID          = Me%InstanceID

            STAT_ = SUCCESS_

        else cd0
            
            stop 'ModuleNetworkStatistics - ConstructNetworkStatistics - ERR01' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructNetworkStatistics
 
    !--------------------------------------------------------------------------

    subroutine ReadKeywords

        !Local--------------------------------------------------------------
        integer                 :: status, flag, STAT_CALL
        
        !Begin--------------------------------------------------------------

    
        Me%ObjEnterData = 0
        
        call ConstructEnterData(Me%ObjEnterData, "NetworkStatistics.dat", STAT = STAT_CALL)
        if(STAT_CALL /= SUCCESS_) stop 'ModuleNetworkStatistics - ReadKeywords - ERR10'
        
        
        
        call GetData(Me%DTAnalysis,                                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='DT_ANALYSIS',                                       &
                     Default      = Day_,                                               &
                     ClientModule ='ModuleNetworkStatistics',                           &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleNetworkStatistics - ReadKeywords - ERR20'
        
        call GetData(Me%InputFile,                                                      &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='INPUT_FILE',                                        &
                     ClientModule ='ModuleNetworkStatistics',                           &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleNetworkStatistics - ReadKeywords - ERR30'
        if (flag == 0) then
            write(*,*) 'Needs the input file'
            stop 'ModuleNetworkStatistics - ReadKeywords - ERR40'
        endif
        
        call GetData(Me%OutputFile,                                                     &
                     Me%ObjEnterData,                                                   &
                     flag,                                                              &
                     SearchType   = FromFile,                                           &
                     keyword      ='OUTPUT_FILE',                                       &
                     ClientModule ='ModuleNetworkStatistics',                           &
                     STAT         = STAT_CALL)        
        if(STAT_CALL /= SUCCESS_) stop 'ModuleNetworkStatistics - ReadKeywords - ERR50'
        if (flag == 0) then
            write(*,*) 'Needs the output file'
            stop 'ModuleNetworkStatistics - ReadKeywords - ERR60'
        endif

        


    
    end subroutine ReadKeywords
    
    !-------------------------------------------------------------------------
 
    
    subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_NetworkStatistics), pointer                         :: NewObjNetworkStatistics
        type (T_NetworkStatistics), pointer                         :: PreviousObjNetworkStatistics


        !Allocates new instance
        allocate (NewObjNetworkStatistics)
        nullify  (NewObjNetworkStatistics%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjNetworkStatistics)) then
            FirstObjNetworkStatistics         => NewObjNetworkStatistics
            Me                    => NewObjNetworkStatistics
        else
            PreviousObjNetworkStatistics      => FirstObjNetworkStatistics
            Me                    => FirstObjNetworkStatistics%Next
            do while (associated(Me))
                PreviousObjNetworkStatistics  => Me
                Me                => Me%Next
            enddo
            Me                    => NewObjNetworkStatistics
            PreviousObjNetworkStatistics%Next => NewObjNetworkStatistics
        endif

        Me%InstanceID = RegisterNewInstance (mNetworkStatistics_)


    end subroutine AllocateInstance

    !--------------------------------------------------------------------------

    subroutine OpenHDF5Files

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        logical                                     :: exist
        integer                                     :: HDF5_READ, HDF5_CREATE

      
        !Begin-----------------------------------------------------------------

        !Verifies if file exists
        inquire(FILE = Me%InputFile, EXIST = exist)
        if (.not. exist) then
            write(*,*)'HDF5 file does not exist:'//trim(Me%InputFile)
            stop 'OpenHDF5Files - ModuleNetworkStatistics - ERR10'
        endif

        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        !Open HDF5 file
        call ConstructHDF5 (Me%ObjHDF5_In, trim(Me%InputFile),                          &
                            HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)                                                      &
        stop 'OpenHDF5Files - ModuleNetworkStatistics - ERR20'
        
        !Obtain start and end times of HDF5 file
        !(obtain number of instants) 
        call GetHDF5GroupNumberOfItems(Me%ObjHDF5_In, "/Time", Me%NInstants, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5Files - ModuleNetworkStatistics - ERR30'
        
        !Create HDF5 file
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)

        !Opens HDF File
        call ConstructHDF5      (Me%ObjHDF5_Out, trim(Me%OutputFile),                   &
                                 HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'OpenHDF5Files - ModuleNetworkStatistics - ERR40'
        
        
    end subroutine OpenHDF5Files

    !--------------------------------------------------------------------------

    subroutine ConstructInputTime

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, i
        real,   dimension(:), pointer               :: Aux6
        character(StringLength)                     :: obj_name
        integer                                     :: obj_type
        integer(HID_T)                              :: FileID_In
      
        !Begin-----------------------------------------------------------------
        
        allocate(Me%InputTime(Me%Ninstants))

        allocate(Aux6(6))        
        
        call GetHDF5FileID (Me%ObjHDF5_In, FileID_In,   STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CopyNetwork - ModuleNetworkStatistics - ERR10'

        do i = 1, Me%Ninstants

            !Gets information about the group
            call h5gget_obj_info_idx_f(FileID_In, "/Time", i-1, obj_name, obj_type,  & 
                                       STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'ConstructInputTime - ModuleNetworkStatistics - ERR10'
            endif
            
            
            call HDF5SetLimits  (Me%ObjHDF5_In, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ConstructInputTime - ModuleNetworkStatistics - ERR20'
            

            call HDF5ReadData   (HDF5ID         = Me%ObjHDF5_In,                &
                                 GroupName      = "/Time",                      &
                                 Name           = trim(obj_Name),               &
                                 Array1D        = Aux6,                         &
                                 STAT           = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)                                          &
                stop 'ConstructInputTime - ModuleNetworkStatistics - ERR30'
                
            call SetDate(Me%InputTime(i), Aux6(1), Aux6(2), Aux6(3), Aux6(4), Aux6(5), Aux6(6))
            
            if (Me%FirstOutPutInst < -99 .and. ) then
            
            
        enddo                            
            
        deallocate(Aux6)
        
    end subroutine ConstructInputTime



    !--------------------------------------------------------------------------


    subroutine ConstructDemandPatterns

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: GroupExist 
        integer                                     :: i, ClientNumber, line, iflag
        integer                                     :: FirstLine, LastLine, STAT_CALL
        logical                                     :: BlockFound
      
        !Begin-----------------------------------------------------------------
        
        real, dimension(:,:), pointer                           :: PeakDemandStart, PeakDemandEnd
        real, dimension(:,:), pointer                           :: LowDemandStart,  LowDemandEnd

        

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginDemandPattern>',           &
                                    block_end       = '<EndDemandPattern>',             &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS:     if (STAT_CALL == SUCCESS_) then

BF:         if (BlockFound) then
                 
                Me%NPatterns = LastLine - FirstLine - 1
                
                allocate(Me%AnalysisGroups(Me%NGroups))

                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(Me%AnalysisGroups(i), EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR10'
                    if (iflag == 0) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR20'
                enddo
            
            else BF
                Me%NGroups = 0.
            endif BF
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS
            
            stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR30'
        
        endif IS           
        


        do i=1, Me%NGroups

            call GetHDF5GroupExist (Me%ObjHDF5_In, Me%AnalysisGroups(i), GroupExist)

            !check if file contains parameter required
            if (.NOT. GroupExist) then  
                write(*,*)'HDF5 file do not contain parameter required:'            &
                           //trim(Me%InputFile)
                write(*,*)'Parameter required:'//trim(Me%AnalysisGroups(i))
                stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR40'
            end if
            
        enddo        


        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginGroupsCopy>',              &
                                    block_end       = '<EndGroupsCopy>',                &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS1:    if (STAT_CALL == SUCCESS_) then

BF1:        if (BlockFound) then
                 
                Me%NCopyGroups = LastLine - FirstLine - 1
                
                allocate(Me%CopyGroups(Me%NCopyGroups))

                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(Me%CopyGroups(i), EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR50'
                    if (iflag == 0) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR60'
                enddo
            
            else BF1
                Me%NCopyGroups = 0.
            endif BF1
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS1
            
            stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR70'
        
        endif IS1           
        


        do i=1, Me%NCopyGroups

            call GetHDF5GroupExist (Me%ObjHDF5_In, Me%CopyGroups(i), GroupExist)

            !check if file contains parameter required
            if (.NOT. GroupExist) then  
                write(*,*)'HDF5 file do not contain parameter required:'            &
                           //trim(Me%InputFile)
                write(*,*)'Parameter required:'//trim(Me%CopyGroups(i))
                stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR80'
            end if
            
        enddo        
        
    end subroutine ConstructVGroupNames


    !--------------------------------------------------------------------------


    subroutine ConstructVGroupNames

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        logical                                     :: GroupExist 
        integer                                     :: i, ClientNumber, line, iflag
        integer                                     :: FirstLine, LastLine, STAT_CALL
        logical                                     :: BlockFound
      
        !Begin-----------------------------------------------------------------

        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginGroupsStat>',              &
                                    block_end       = '<EndGroupsStat>',                &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS:     if (STAT_CALL == SUCCESS_) then

BF:         if (BlockFound) then
                 
                Me%NGroups = LastLine - FirstLine - 1
                
                allocate(Me%AnalysisGroups(Me%NGroups))

                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(Me%AnalysisGroups(i), EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR10'
                    if (iflag == 0) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR20'
                enddo
            
            else BF
                Me%NGroups = 0.
            endif BF
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS
            
            stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR30'
        
        endif IS           
        


        do i=1, Me%NGroups

            call GetHDF5GroupExist (Me%ObjHDF5_In, Me%AnalysisGroups(i), GroupExist)

            !check if file contains parameter required
            if (.NOT. GroupExist) then  
                write(*,*)'HDF5 file do not contain parameter required:'            &
                           //trim(Me%InputFile)
                write(*,*)'Parameter required:'//trim(Me%AnalysisGroups(i))
                stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR40'
            end if
            
        enddo        


        call ExtractBlockFromBuffer(Me%ObjEnterData,                                    &
                                    ClientNumber    = ClientNumber,                     &
                                    block_begin     = '<BeginGroupsCopy>',              &
                                    block_end       = '<EndGroupsCopy>',                &
                                    BlockFound      = BlockFound,                       &
                                    FirstLine       = FirstLine,                        &
                                    LastLine        = LastLine,                         &
                                    STAT            = STAT_CALL)
IS1:    if (STAT_CALL == SUCCESS_) then

BF1:        if (BlockFound) then
                 
                Me%NCopyGroups = LastLine - FirstLine - 1
                
                allocate(Me%CopyGroups(Me%NCopyGroups))

                i=0
                do line=FirstLine +1, LastLine-1
                    i = i + 1
                    call GetData(Me%CopyGroups(i), EnterDataID = Me%ObjEnterData, flag = iflag, &
                                 Buffer_Line = line, STAT = STAT_CALL) 
                    if (STAT_CALL /= SUCCESS_) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR50'
                    if (iflag == 0) stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR60'
                enddo
            
            else BF1
                Me%NCopyGroups = 0.
            endif BF1
             
            call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL)                  
        
        else IS1
            
            stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR70'
        
        endif IS1           
        


        do i=1, Me%NCopyGroups

            call GetHDF5GroupExist (Me%ObjHDF5_In, Me%CopyGroups(i), GroupExist)

            !check if file contains parameter required
            if (.NOT. GroupExist) then  
                write(*,*)'HDF5 file do not contain parameter required:'            &
                           //trim(Me%InputFile)
                write(*,*)'Parameter required:'//trim(Me%CopyGroups(i))
                stop 'ConstructVGroupNames - ModuleNetworkStatistics - ERR80'
            end if
            
        enddo        
        
    end subroutine ConstructVGroupNames

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    
    !--------------------------------------------------------------------------
    subroutine GetNetworkStatisticsPointer (ObjNetworkStatisticsID, Matrix, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjNetworkStatisticsID
        real(8), dimension(:, :, :),  pointer           :: Matrix
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjNetworkStatisticsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mNetworkStatistics_, Me%InstanceID)

            !Matrix => Me%Matrix

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetNetworkStatisticsPointer
    
    !--------------------------------------------------------------------------
    
    subroutine GetNetworkStatisticsInteger (ObjNetworkStatisticsID, Int, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: ObjNetworkStatisticsID
        real                                            :: Int
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjNetworkStatisticsID, ready_)

        if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                            &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            Int = Me%InstanceID

            STAT_ = SUCCESS_

        else 
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine GetNetworkStatisticsInteger

    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ModifyNetworkStatistics(ObjNetworkStatisticsID, STAT)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjNetworkStatisticsID
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjNetworkStatisticsID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then
        
            call CopyNetwork
            
            !call ModifyNetworkStatistics

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine ModifyNetworkStatistics
    
    !--------------------------------------------------------------------------

    subroutine CopyNetwork

        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        integer(HID_T)                                          :: FileID_In, gr_id
        integer                                                 :: STAT_CALL, i

        !Begin-----------------------------------------------------------------
    
        call GetHDF5FileID (Me%ObjHDF5_In, FileID_In,   STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CopyNetwork - ModuleNetworkStatistics - ERR10'
        
        do i=1, Me%NCopyGroups
            call h5gopen_f (FileID_In, trim(adjustl(Me%CopyGroups(i))), gr_id, STAT_CALL)
            call CopyGroup (gr_id,     trim(adjustl(Me%CopyGroups(i))))
        enddo            
        
    end subroutine CopyNetwork        
    
   !--------------------------------------------------------------------------

    recursive subroutine CopyGroup (ID, GroupName)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: GroupName
        integer(HID_T)                              :: ID
        
        !Local-----------------------------------------------------------------
        character(StringLength)                     :: obj_name
        integer                                     :: obj_type
        integer(HID_T)                              :: gr_id, dset_id
        integer(HID_T)                              :: datatype_id, class_id, size        
        integer                                     :: STAT_CALL
        character(StringLength)                     :: NewGroupName, LastGroupName
        integer                                     :: ItensNumber
        integer                                     :: i, imax, jmax
        character(len=StringLength)                 :: Name
        logical                                     :: TimeIndependent = .false.
        real(4),  dimension(:), pointer             :: ArrayReal1D
        integer,  dimension(:), pointer             :: ArrayInt1D
        real(4),  dimension(:,:), pointer           :: ArrayReal2D
        integer,  dimension(:,:), pointer           :: ArrayInt2D
        
        character(len=StringLength)                 :: Units
        integer                                     :: Rank
        integer,dimension(7)                        :: Dimensions        
        
        !Begin-----------------------------------------------------------------

        call HDF5CreateGroup  (Me%ObjHDF5_Out, GroupName, STAT = STAT_CALL)    
        if (STAT_CALL /= SUCCESS_) stop 'InquireSubGroup - ModuleHDF5Extractor - ERR10'

        ItensNumber = 0

        !Get the number of members in the Group
        call h5gn_members_f(ID, GroupName, ItensNumber, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) return
  
        do i = 1, ItensNumber

            !Gets information about the group
            call h5gget_obj_info_idx_f(ID, GroupName, i-1, obj_name, obj_type,  & 
                                       STAT_CALL)
            if (STAT_CALL /= SUCCESS_) exit

            if (obj_type == H5G_DATASET_F) then

                !Get item specifics
                call GetHDF5GroupID(Me%ObjHDF5_In, GroupName, i,                    &
                                    obj_name,                                       &
                                    Rank       = Rank,                              &
                                    Dimensions = Dimensions,                        &
                                    !Units      = Units,                             &
                                    STAT       = STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'InquireSubGroup - ModuleHDF5Extractor - ERR10'
                
                !Dummy value
                Units='-'
                
                if (Rank > 2) then
                    write(*,*) 'In network type files the fields are always assumed 1D or 2D' 
                    stop 'InquireSubGroup - ModuleHDF5Extractor - ERR10'
                endif

                !Get data type (integer or real)
                !(for time dependent itens assumed that data type equal for all fields)
                !Opens data set
                call h5dopen_f     (ID, trim(adjustl(obj_name)), dset_id, STAT_CALL)
                !Gets datatype
                call h5dget_type_f (dset_id, datatype_id,   STAT_CALL)
                !call h5tget_size_f (datatype_id, size,      STAT_CALL)
                call h5tget_class_f(datatype_id, class_id,  STAT_CALL) 
                
                
                imax = Dimensions(1)
                if (Rank == 1) then
                    
                    call HDF5SetLimits  (Me%ObjHDF5_In, 1, imax, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR10'

                    call HDF5SetLimits  (Me%ObjHDF5_Out, 1, imax, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR10'

                else
                    jmax = Dimensions(2)
                    
                    call HDF5SetLimits  (Me%ObjHDF5_In, 1, imax, 1, jmax, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR10'

                    call HDF5SetLimits  (Me%ObjHDF5_Out, 1, imax, 1, jmax, STAT = STAT_CALL)
                    if (STAT_CALL /= SUCCESS_) stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR10'
                    
                endif
                if     (class_id == H5T_FLOAT_F  ) then
                
                    if (Rank==1) then
                        allocate(ArrayReal1D(1:imax))
                        
                        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5_In,                &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Array1D        = ArrayReal1D,                  &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'

                        call HDF5WriteData  (HDF5ID         = Me%ObjHDF5_Out,               &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Units          = Units,                        &
                                             Array1D        = ArrayReal1D,                  &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'
                        
                        deallocate(ArrayReal1D)
                    
                    else

                        allocate(ArrayReal2D(1:imax,1:jmax))
                        
                        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5_In,                &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Array2D        = ArrayReal2D,                  &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'

                        call HDF5WriteData  (HDF5ID         = Me%ObjHDF5_Out,               &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Units          = Units,                        &
                                             Array2D        = ArrayReal2D,                  &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'
                    
                        deallocate(ArrayReal2D)
                    
                    endif

                    
                elseif (class_id == H5T_INTEGER_F) then

                    if (Rank==1) then
                        allocate(ArrayInt1D(1:imax))
                       
                        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5_In,                &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Array1D        = ArrayInt1D,                   &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'

                        call HDF5WriteData  (HDF5ID         = Me%ObjHDF5_Out,               &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Units          = Units,                        &                                             
                                             Array1D        = ArrayInt1D,                   &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'
                        
                        deallocate(ArrayInt1D)
                    else
                        allocate(ArrayInt2D(1:imax,1:jmax))
                       
                        call HDF5ReadData   (HDF5ID         = Me%ObjHDF5_In,                &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Array2D        = ArrayInt2D,                   &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'

                        call HDF5WriteData  (HDF5ID         = Me%ObjHDF5_Out,               &
                                             GroupName      = "/"//trim(GroupName),         &
                                             Name           = trim(obj_Name),               &
                                             Units          = Units,                        &                                             
                                             Array2D        = ArrayInt2D,                   &
                                             STAT           = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_)                                          &
                        stop 'HDF5TimeInstant - ModuleHDF5Extractor - ERR20'
                        
                        deallocate(ArrayInt2D)
                        
                    endif                                            
                                            
                else
                    stop 'InquireSubGroup - ModuleHDF5Extractor - ERR20'
                end if
                
             elseif (obj_type == H5G_GROUP_F) then
             
                LastGroupName = GroupName

                if (GroupName == "/") then
                    NewGroupName = GroupName//trim(adjustl(obj_name))
                else
                    NewGroupName = GroupName//"/"//trim(adjustl(obj_name))
                endif
                call h5gopen_f (ID, trim(adjustl(NewGroupName)), gr_id,     &    
                                STAT_CALL)
                call CopyGroup (gr_id, trim(adjustl(NewGroupName)))
                call h5gclose_f (gr_id, STAT_CALL)

            endif

        enddo

    end subroutine CopyGroup

    !--------------------------------------------------------------------------
    
    
    !--------------------------------------------------------------------------   

    function FreqAnalysis(SortArray,SizeArray, Percentil)    

        !Arguments-----------------------------
        real :: Percentil, FreqAnalysis
        integer :: SizeArray
        real, dimension(:) :: SortArray
        !Local---------------------------------
        real       :: Aux, Raux
        integer    :: Iaux
        
        !Begin---------------------------------
        if (SizeArray==1) then
            FreqAnalysis=SortArray(1)
        else
            Aux = real(SizeArray-1)*Percentil
            Iaux = int(Aux)
            Raux = Aux-real(Iaux)
            if (Percentil == 0 .or. Iaux == 0) then
                FreqAnalysis = SortArray(1)
            else
                if (Raux>0.) then
                    FreqAnalysis = SortArray(Iaux+1)*Raux+SortArray(Iaux)*(1.-Raux)
                else
                    FreqAnalysis = SortArray(Iaux)
                endif
            endif
        endif
    
    
    end function FreqAnalysis
    
    
    logical function WeekWayDay (Date)
    
        !Arguments-----------------------------    
        type (T_Time) :: Date
           
        !Local---------------------------------
        type (T_Time) :: AuxDate
        real(8)       :: AuxSeconds, Weeks, AuxWeek 
        !Begin---------------------------------

        call SetDate (AuxDate, 2013.,2.,11.,0.,0.,0.)
        
        AuxSeconds = Date-AuxDate
        
        if (AuxSeconds >= 0.) then
            Weeks    = AuxSeconds / (168.*3600.)
            AuxWeek  = Weeks - int(Weeks)
            if (AuxWeek <5./7.) then
                WeekWayDay = .true.
            else
                WeekWayDay = .false.
            endif    
        else
            AuxSeconds = - AuxSeconds
            Weeks    = AuxSeconds / (168.*3600.)
            AuxWeek  = Weeks - int(Weeks)
            if (AuxWeek >2./7.) then
                WeekWayDay = .true.
            else
                WeekWayDay = .false.
            endif             
        endif 
    
    end function WeekWayDay 
    
    !----------------------------------------------------------------

    logical function WeekEndDay (Date)
    
        !Arguments-----------------------------    
        type (T_Time) :: Date
           
        !Local---------------------------------
        !Begin---------------------------------

        if (WeekWayDay(Date)) then
            
            WeekEndDay = .false.
            
        else
        
            WeekEndDay = .true.        
        
        endif
    
    end function WeekEndDay 
    
    

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    subroutine KillNetworkStatistics(ObjNetworkStatisticsID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjNetworkStatisticsID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjNetworkStatisticsID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mNetworkStatistics_,  Me%InstanceID)

            if (nUsers == 0) then
            
                call KillVariablesAndFiles
            
                !Deallocates Instance
                call DeallocateInstance ()

                ObjNetworkStatisticsID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_
           

        !------------------------------------------------------------------------

    end subroutine KillNetworkStatistics
        

    !------------------------------------------------------------------------
    



    subroutine KillVariablesAndFiles
    
    
        !Local--------------------------------------------------------------------------
        integer         :: STAT_CALL
        
        !Begin--------------------------------------------------------------------------

        !Kill HDF5 file
        call KillHDF5 (Me%ObjHDF5_In, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillVariablesAndFiles - ModuleNetworkStatistics - ERR10'

        !Kill HDF5 File
        call KillHDF5 (Me%ObjHDF5_Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'KillVariablesAndFiles - ModuleNetworkStatistics - ERR20'

    end subroutine KillVariablesAndFiles    
    
    !--------------------------------------------------------------------------    
    
   !--------------------------------------------------------------------------
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_NetworkStatistics), pointer          :: AuxObjNetworkStatistics
        type (T_NetworkStatistics), pointer          :: PreviousObjNetworkStatistics

        !Updates pointers
        if (Me%InstanceID == FirstObjNetworkStatistics%InstanceID) then
            FirstObjNetworkStatistics => FirstObjNetworkStatistics%Next
        else
            PreviousObjNetworkStatistics => FirstObjNetworkStatistics
            AuxObjNetworkStatistics      => FirstObjNetworkStatistics%Next
            do while (AuxObjNetworkStatistics%InstanceID /= Me%InstanceID)
                PreviousObjNetworkStatistics => AuxObjNetworkStatistics
                AuxObjNetworkStatistics      => AuxObjNetworkStatistics%Next
            enddo

            !Now update linked list
            PreviousObjNetworkStatistics%Next => AuxObjNetworkStatistics%Next

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

    subroutine Ready (ObjNetworkStatistics_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjNetworkStatistics_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjNetworkStatistics_ID > 0) then
            call LocateObjNetworkStatistics (ObjNetworkStatistics_ID)
            ready_ = VerifyReadLock (mNetworkStatistics_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjNetworkStatistics (ObjNetworkStatisticsID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjNetworkStatisticsID

        !Local-----------------------------------------------------------------

        Me => FirstObjNetworkStatistics
        do while (associated (Me))
            if (Me%InstanceID == ObjNetworkStatisticsID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleNetworkStatistics - LocateObjNetworkStatistics - ERR01'

    end subroutine LocateObjNetworkStatistics

    !--------------------------------------------------------------------------

    end module ModuleNetworkStatistics
!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : GlueHDF5Files
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : October 2003
! REVISION      : Paulo Leitao - v4.0
! DESCRIPTION   : Module to convert Glue HDF5 Files. The files must have the same VGroups
!
!------------------------------------------------------------------------------


Module ModuleGlueHDF5Files

    use HDF5
    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
!    NIX
#ifndef _USE_NIX
    use DFWIN
#endif

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: StartGlueHDF5Files
    private ::      ReadOptions
    private ::      GlueProcess
    private ::      KillGlueHDF5Files

    interface ReadInterface
        module procedure ReadInterface1DR4
        module procedure ReadInterface1DR8
        module procedure ReadInterface2DR4
        module procedure ReadInterface2DR8
        module procedure ReadInterface3DR4
        module procedure ReadInterface3DR8
        module procedure ReadInterface1DI4
        module procedure ReadInterface2DI4
        module procedure ReadInterface3DI4
    end interface

    interface WriteInterface
        module procedure WriteInterface1DR4
        module procedure WriteInterface1DR8
    end interface



    
    
    !Types---------------------------------------------------------------------
    
    private :: T_TimeMap
    type       T_TimeMap
        integer,        dimension(:), pointer            :: FirstInstantBest   
        integer,        dimension(:), pointer            :: LastInstantBest    
        integer,        dimension(:), pointer            :: ObjHDF5_ID        
        integer(HID_T), dimension(:), pointer            :: File_ID
        type (T_Time),  dimension(:), pointer            :: FirstDateBest      
        type (T_Time),  dimension(:), pointer            :: LastDateBest       
        integer                                          :: PresentFile
        logical                                          :: BestTimeSerieON     = .true.
    end type  T_TimeMap
    
    
    private :: T_GlueHDF5Files
    type       T_GlueHDF5Files
        integer                                          :: ObjEnterData         = 0
        integer                                          :: ObjHDF5_In           = 0
        integer                                          :: ObjHDF5_Out          = 0
        character(len=PathLength), dimension(:), pointer :: FileNameIn           
        integer, dimension(:), pointer                   :: FirstInstant         
        type (T_Time)                                    :: LastInstant
        character(len=PathLength)                        :: FileNameOut
        integer                                          :: FileNameInNumber
        logical                                          :: File_3D
        logical                                          :: Vert3D, Open3D
        logical                                          :: GlueInTime
        character(len=PathLength)                        :: BaseGroup
        character(len=PathLength)                        :: TimeGroup
        type (T_TimeMap)                                 :: TimeMap
        logical                                          :: CheckHDF5_File      = .false.   
    end type  T_GlueHDF5Files

    type(T_GlueHDF5Files), pointer                       :: Me                  

    !--------------------------------------------------------------------------
    
    contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine StartGlueHDF5Files(EnterDataID, ClientNumber, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID, ClientNumber
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_
        
        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions(ClientNumber)
        
        call GlueOnlyCommonVGroups

        call GlueProcess

        call KillGlueHDF5Files


        STAT = SUCCESS_


    end subroutine StartGlueHDF5Files

    !------------------------------------------------------------------------

    subroutine ReadOptions(ClientNumber)
    
        !Arguments
        integer                                     :: ClientNumber

        !Local-----------------------------------------------------------------
        character(len=PathLength), dimension(:), pointer :: FileName1DAux
        character(len=PathLength)                   :: FileNameAux
        logical                                     :: BlockFound
        integer                                     :: STAT_CALL
        integer                                     :: iflag, FirstLine, LastLine, i, iaux1, iaux2

        !Begin-----------------------------------------------------------------
       

 
        call GetData(Me%FileNameOut,                                                    &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'OUTPUTFILENAME',                                   &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR10'


        call GetData(Me%File_3D,                                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = '3D_FILE',                                          &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR20'

        call GetData(Me%Vert3D,                                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = '3D_VERT',                                          &
                     default      = .false.,                                            &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR30'
        
        if (Me%File_3D) Me%Vert3D =.true.

        call GetData(Me%Open3D,                                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = '3D_OPEN',                                          &
                     default      = .false.,                                            &                     
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR40'

!ROD: next line was commented since it is inconsistent with previous keyword
!        if (Me%File_3D) Me%Open3D =.true.        

        call GetData(Me%BaseGroup,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'BASE_GROUP',                                       &
                     Default      = 'Results',                                          &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR50'

        call GetData(Me%TimeGroup,                                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'TIME_GROUP',                                       &
                     Default      = 'Time',                                             &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR60'
        
        !Check if the user wants to glue files in time (default)
        !if false the time is the same in all the files and will add the files results
        call GetData(Me%GlueInTime,                                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'GLUE_IN_TIME',                                     &
                     Default      = .true.,                                             &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR70'
        
        if (Me%GlueInTime) then
            call GetData(Me%TimeMap%BestTimeSerieON,                                    &
                         Me%ObjEnterData, iflag,                                        &
                         SearchType   = FromBlock,                                      &
                         keyword      = 'BEST_TIME_SERIE',                              &
                         Default      = .false.,                                        &
                         ClientModule = 'ModuleGlueHDF5Files',                          &
                         STAT         = STAT_CALL)        
            if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR80'
        else
            Me%TimeMap%BestTimeSerieON = .false.
        endif            
        
        call GetData(Me%CheckHDF5_File,                                                 &
                     Me%ObjEnterData , iflag,                                           &
                     SearchType   = FromBlock,                                          &
                     keyword      = 'CHECK_HDF5_FILE',                                  &
                     default      = .true.,                                             &
                     ClientModule = 'ModuleGlueHDF5Files',                              &
                     STAT         = STAT_CALL)                                      
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR90'
        
        
do1 :   do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,                   &
                                        '<<begin_list>>', '<<end_list>>', BlockFound,   &
                                        FirstLine = FirstLine, LastLine = LastLine,     &
                                        STAT = STAT_CALL)

if1 :       if(STAT_CALL .EQ. SUCCESS_) then    
if2 :           if (BlockFound) then

                    iAux1 = LastLine - FirstLine - 1

                    allocate (FileName1DAux  (iAux1))
                    
                    iAux2 = 0

                    do i = 1, iAux1

                        call GetData(FileNameAux, Me%ObjEnterData,  iflag,              & 
                                     Buffer_Line  = FirstLine + i,                      &
                                     STAT         = STAT_CALL)
                        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleGlueHDF5Files - ERR100'
                        
                        if (GetHDF5FileOkToRead(FileNameAux)) then
                            iAux2                = iAux2 + 1
                            FileName1DAux(iAux2) = FileNameAux
                        endif

                    enddo
                    
                    Me%FileNameInNumber = iAux2
                                        
                    allocate (Me%FileNameIn  (Me%FileNameInNumber))
                    allocate (Me%FirstInstant(Me%FileNameInNumber))
                    
                    Me%FileNameIn(1:Me%FileNameInNumber) = FileName1DAux(1:Me%FileNameInNumber)
                    
                    deallocate(FileName1DAux)
                    
                    exit do1 
                        
                end if if2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then if1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                if(STAT_CALL .ne. SUCCESS_)stop 'ReadOptions - ModuleGlueHDF5Files - ERR110'
                    
            end if if1
        end do do1


    end subroutine ReadOptions

    !--------------------------------------------------------------------------
    

    subroutine GlueProcess

        !Local-----------------------------------------------------------------
        integer                                     :: i, HDF5_READWRITE, STAT_CALL, iflag
        
!       NIX   
#ifdef _USE_NIX     
        integer                                     :: system
        character(PathLength)                       :: aux
#endif

        !Begin-----------------------------------------------------------------

!         NIX
#ifdef _USE_NIX     
        write(aux, *) 'cp ',trim(Me%FileNameIn(1)), ' ',trim(Me%FileNameOut)      
        iflag = system (aux)
        if (iflag /= 0) stop 'GlueProcess - ConvertHDF5Files - ERR01'
#else
        iflag = copyfile(trim(Me%FileNameIn(1))//""C,trim(Me%FileNameOut)//""C,.false.)
        if (iflag /= 1) stop 'GlueProcess - ConvertHDF5Files - ERR01'
#endif
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
        

        call ConstructHDF5 (Me%ObjHDF5_Out, Me%FileNameOut,                            &
                            Access = HDF5_READWRITE, STAT = STAT_CALL)

        
        if (Me%TimeMap%BestTimeSerieON) then
        
            call BestTimeSerieGlue
        
        else
        
            call InquireFile(Me%FileNameIn(1))        
        


            !only verify group compatibility if glueing in time
            !if adding results groups need to verify if time is the same
            if (Me%GlueInTime) then
                write (*,*)
                write (*,*) 'Glueing HDF files...'
            
                do i=2, Me%FileNameInNumber
                
                    call CheckVGCompatibility(i)

                enddo
            else
                write (*,*)
                write (*,*) 'Merging HDF files...'    
            
                do i=2, Me%FileNameInNumber
                
                    call InquireFile(Me%FileNameIn(i))
                
                    !check if times are the same and bathymetries have same dimension
                    call CheckTimeAndBathymetry(i)
                
                    !add the new groups not existing in the output
                    call CheckGroupExistence(i)
                
                    !neded for glue. first instant is always used because no glueintime occurs
                    Me%FirstInstant(i) = 1
                
                enddo        
            endif

            do i=2, Me%FileNameInNumber
        
                call GlueFileIn(i)

            enddo

            if (Me%GlueInTime) then 
                write (*,*)  
                write (*,*) 'Finished Glueing HDF files!'      
            else
                write (*,*)
                write (*,*) 'Finished Merging HDF files!'
            endif
            
        endif            

    end subroutine GlueProcess

    !--------------------------------------------------------------------------
    

    
    subroutine GlueOnlyCommonVGroups
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j
        logical                                     :: NoDelete

        !Begin-----------------------------------------------------------------
        
        NoDelete = .true.
        
        do While (NoDelete)
            
            NoDelete = .false. 

            do i=1, Me%FileNameInNumber
            do j=1, Me%FileNameInNumber   

                if (i /= j) then    
                    call DeleteVGNotCommon(i, j, NoDelete)
                endif                    

            enddo    
            enddo        
        enddo
        
    end subroutine GlueOnlyCommonVGroups        
        
                
    !--------------------------------------------------------------------------
    subroutine DeleteVGNotCommon(i, j, NoDelete)
    
        !Arguments-------------------------------------------------------------
        integer                                     :: i, j
        logical                                     :: NoDelete

        !Local-----------------------------------------------------------------
        integer                                     :: HDF5_1, HDF5_2
        integer                                     :: HDF5_READWRITE, STAT_CALL, ID1, ID2

        !Begin-----------------------------------------------------------------
       
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READWRITE = HDF5_READWRITE)
        
        HDF5_1 = 0
        HDF5_2 = 0

        call ConstructHDF5 (HDF5_1, Me%FileNameIn(i),                                   &
                            Access = HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DeleteVGNotCommon - ModuleGlueHDF5Files - ERR10'        
       
        call GetHDF5FileID (HDF5_1, ID1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DeleteVGNotCommon - ModuleGlueHDF5Files - ERR20'        
        
        call ConstructHDF5 (HDF5_2, Me%FileNameIn(j),                                   &
                            Access = HDF5_READWRITE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DeleteVGNotCommon - ModuleGlueHDF5Files - ERR30'        
       
        call GetHDF5FileID (HDF5_2, ID2, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DeleteVGNotCommon - ModuleGlueHDF5Files - ERR40'        
        
        call DeleteNoMatchingSubGroups (ID1, ID2, "/", NoDelete)

        call KillHDF5(HDF5_1, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DeleteVGNotCommon - ModuleGlueHDF5Files - ERR50'

        call KillHDF5(HDF5_2, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'DeleteVGNotCommon - ModuleGlueHDF5Files - ERR60'

        
    end subroutine DeleteVGNotCommon

    
    !--------------------------------------------------------------------------

    recursive subroutine DeleteNoMatchingSubGroups (IDOut, IDIn, GroupNameOut, NoDelete)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn
        character(len=*)                            :: GroupNameOut
        logical                                     :: NoDelete

        !Local-----------------------------------------------------------------
        logical                                     :: FoundOK        
        integer                                     :: nmembersOut, nmembersIn
        character(StringLength)                     :: obj_nameIn, obj_nameOut
        integer                                     :: obj_type, idxOut, idxIN
        integer                                     :: STAT

        !Begin-----------------------------------------------------------------
        
        !Turns Error printing of        
        call h5eset_auto_f  (0, STAT)        

        !Get the number of members in the Group        
        call h5gn_members_f(IDOut, GroupNameOut, nmembersOut, STAT)
        
        !Turns Error printing on
        call h5eset_auto_f  (1, STAT)        
        
        if (STAT /= SUCCESS_) then
            return
            !stop 'DeleteNoMatchingSubGroups - ModuleGlueHDF5Files - ERR10'
        endif            
    

        do idxOut = 1, nmembersOut

            !Gets information about the group
            call h5gget_obj_info_idx_f(IDOut, GroupNameOut, idxOut-1, obj_nameOut, obj_type, STAT)
            !if (STAT /= SUCCESS_) then
            !    stop 'DeleteNoMatchingSubGroups - ModuleGlueHDF5Files - ERR30'
            !endif                            
            
            if  (obj_type ==H5G_GROUP_F) then
                
                call h5gn_members_f(IDIn, GroupNameOut, nmembersIn, STAT)
                if (STAT /= SUCCESS_) then
                    return
                    !stop 'DeleteNoMatchingSubGroups - ModuleGlueHDF5Files - ERR20'
                endif                   
            
                FoundOK = .false.
                do idxIN = 1, nmembersIN
                    call h5gget_obj_info_idx_f(IDIn, GroupNameOut, idxIN-1, obj_nameIn, obj_type, STAT)
                    if (STAT /= SUCCESS_) then
                        stop 'DeleteNoMatchingSubGroups - ModuleGlueHDF5Files - ERR10'
                    endif      
                
                    if (trim(obj_nameOut) == trim(obj_nameIn))  then
                        FoundOK = .true.
                        exit
                    endif    
                enddo
                    
                    
                if (FoundOK) then                    
                    call DeleteNoMatchingSubGroups (IDOut, IDIn, trim(obj_nameOut), NoDelete)                    
                else                    
                    call DeleteVGroup(ID = IDOut,  ParentGroupName = GroupNameOut, GroupName = obj_nameOut )         
                    NoDelete = .true.
                endif

            endif
            
        enddo


    end subroutine DeleteNoMatchingSubGroups
    
    
    !--------------------------------------------------------------------------    
    
    subroutine DeleteVGroup(ID, ParentGroupName, GroupName)
    
        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: ID
        character(len=*)                            :: ParentGroupName, GroupName

        !Local-----------------------------------------------------------------
        integer                                     :: STAT
        integer(HID_T)                              :: gr_id

        !Begin-----------------------------------------------------------------        
        
        call h5gopen_f  (ID, trim(ParentGroupName), gr_id, STAT)
        if (STAT /= SUCCESS_) then
            stop 'DeleteVGroup - ModuleGlueHDF5Files - ERR10'
        endif                       
        
        call h5ldelete_f(loc_id     =  gr_id,                                           &
                            name       =  trim(GroupName),                              &
                            hdferr     =  STAT)
        !if (STAT /= SUCCESS_) then
        !    stop 'DeleteVGroup - ModuleGlueHDF5Files - ERR20'
        !endif  
        
        call h5gclose_f       (gr_id, STAT)        
        if (STAT /= SUCCESS_) then
            stop 'DeleteVGroup - ModuleGlueHDF5Files - ERR30'
        endif          
    
    end subroutine DeleteVGroup
    
   !--------------------------------------------------------------------------        
    
    subroutine BestTimeSerieGlue
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, STAT_CALL
        integer                                     :: HDF5_READ
        integer(HID_T)                              :: IDOut, gr_id
        logical                                     :: CheckOK


        !Begin-----------------------------------------------------------------
    
        call GetHDF5FileID (Me%ObjHDF5_Out, IDOut, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR10'  
        
        call h5gopen_f       (IDOut, "/", gr_id, STAT_CALL)        
                          
        call h5ldelete_f(loc_id     =  gr_id,                                            &
                         name       = "/Results",                                        &
                         hdferr     =  STAT_CALL)
        if (STAT_CALL /= 0) stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR20'   
                            
        call h5ldelete_f(loc_id     =  gr_id,                                            &
                         name       = "/Time",                                           &
                         hdferr     =  STAT_CALL)
        if (STAT_CALL /= 0) stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR30'  
        
        
        call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Grid/VerticalZ", Me%Vert3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR40'
        endif
        
        if (Me%Vert3D) then
            call h5ldelete_f(loc_id     =  gr_id,                                       &
                             name       =  "/Grid/VerticalZ",                           &
                             hdferr     =  STAT_CALL)
            if (STAT_CALL /= 0) stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR50'           
        endif

        call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Grid/OpenPoints", Me%Open3D, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) then
            stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR60'
        endif          
        
        if (Me%Open3D) then
            call h5ldelete_f(loc_id     =  gr_id,                                       &
                             name       =  "/Grid/OpenPoints",                          &
                             hdferr     =  STAT_CALL)
            if (STAT_CALL /= 0) stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR70'           
        endif        

        
        allocate(Me%TimeMap%ObjHDF5_ID(1:Me%FileNameInNumber))
        allocate(Me%TimeMap%File_ID   (1:Me%FileNameInNumber))
        
        Me%TimeMap%ObjHDF5_ID(1:Me%FileNameInNumber) = 0
        Me%TimeMap%File_ID   (1:Me%FileNameInNumber) = 0
        
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)
        
        do i=1, Me%FileNameInNumber

            call InquireFile(Me%FileNameIn(i))          

            call ConstructHDF5 (Me%TimeMap%ObjHDF5_ID(i), Me%FileNameIn(i),             &
                                Access = HDF5_READ, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR80'
            endif 
            
            call GetHDF5FileID (HDF5ID = Me%TimeMap%ObjHDF5_ID(i), FileID = Me%TimeMap%File_ID(i), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR90'
            endif 


       enddo
                                
        
        do i=2, Me%FileNameInNumber

            CheckOK = .true.

            call CompareSubGroups (Me%TimeMap%File_ID(i), Me%TimeMap%File_ID(i-1), "/", "/", CheckOK)

            if (.not.CheckOK) then
                write(*,*) trim(Me%FileNameIn(i))//" is not a compatible file"
                stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR100'
            endif

        enddo                            

        
        call AddNewGroupsToOutput(IDin = Me%TimeMap%File_ID(1), GroupName = "/")        
        
        
        call ConstructTimeMap
        
        do i=1, Me%FileNameInNumber
        
            !neded for glue. first instant is always used because no glueintime occurs
            Me%FirstInstant(i) = 1        
        
        enddo        
        
        do i=1, Me%FileNameInNumber
           
            call KillHDF5(Me%TimeMap%ObjHDF5_ID(i), STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                stop 'BestTimeSerieGlue - ModuleGlueHDF5Files - ERR110'
            endif      
        
        enddo
        
        do i=1, Me%FileNameInNumber
           
            call GlueFileIn(i)

        enddo
        
       
        
    
    end subroutine BestTimeSerieGlue
    
    !--------------------------------------------------------------------------  
    
    subroutine ConstructTimeMap    
    
        !Arguments-------------------------------------------------------------    
    
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, STAT_CALL, nItems, PrevInTimeFile
        type (T_Time)                               :: DateAux


        !Begin-----------------------------------------------------------------    
        
        allocate(Me%TimeMap%FirstInstantBest(1:Me%FileNameInNumber))        
        allocate(Me%TimeMap%LastInstantBest (1:Me%FileNameInNumber))                

        allocate(Me%TimeMap%FirstDateBest   (1:Me%FileNameInNumber))        
        allocate(Me%TimeMap%LastDateBest    (1:Me%FileNameInNumber))     
        
        PrevInTimeFile = Me%FileNameInNumber
        
        do i = Me%FileNameInNumber, 1, -1
        
            Me%TimeMap%FirstInstantBest(i) = 0
            Me%TimeMap%LastInstantBest (i) = 0
        

            call GetHDF5GroupNumberOfItems (HDF5ID      = Me%TimeMap%ObjHDF5_ID(i),         &
                                            GroupName   = "/Time",                          &
                                            nItems      = nItems,                           &
                                            STAT        = STAT_CALL)        
            
            if (i == Me%FileNameInNumber) then
            
                Me%TimeMap%FirstInstantBest(i) = 1
                Me%TimeMap%LastInstantBest (i) = nItems  
                
                Me%TimeMap%FirstDateBest   (i) = HDF5TimeInstant(HDF5ID = Me%TimeMap%ObjHDF5_ID(i), Instant = 1)
                Me%TimeMap%LastDateBest    (i) = HDF5TimeInstant(HDF5ID = Me%TimeMap%ObjHDF5_ID(i), Instant = nItems)    
                
            else
                do j = nItems, 1, -1
                    
                    DateAux = HDF5TimeInstant(HDF5ID = Me%TimeMap%ObjHDF5_ID(i), Instant = j)
                    
                    if (DateAux < Me%TimeMap%FirstDateBest   (PrevInTimeFile)) then
                                            
                        Me%TimeMap%FirstInstantBest(i) = 1
                        Me%TimeMap%LastInstantBest (i) = j  
                        
                        Me%TimeMap%FirstDateBest   (i) = HDF5TimeInstant(HDF5ID = Me%TimeMap%ObjHDF5_ID(i), Instant = 1)
                        Me%TimeMap%LastDateBest    (i) = DateAux
                        
                        PrevInTimeFile = i
                        
                        exit

                    endif
                    
                enddo                    
            endif
        
        
        enddo
        
    end subroutine ConstructTimeMap            
        
    !--------------------------------------------------------------------------  
    
   !----------------------------------------------------------------------------


    type(T_Time) function HDF5TimeInstant(HDF5ID, Instant)

        !Arguments-------------------------------------------------------------
        integer                                 :: Instant
        integer                                 :: HDF5ID
        

        !Local-----------------------------------------------------------------
        real,    dimension(:), pointer          :: TimeVector
        integer                                 :: STAT_CALL

        !Begin-----------------------------------------------------------------
        
        allocate(TimeVector(6))

        
        call HDF5SetLimits  (HDF5ID, 1, 6, STAT = STAT_CALL)        

        call HDF5ReadWindow (HDF5ID         = HDF5ID,                                   &
                             GroupName      = "/Time",                                  &
                             Name           = "Time",                                   &
                             Array1D        = TimeVector,                               &
                             OutputNumber   = Instant,                                  &
                             STAT           = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'HDF5TimeInstant - ModuleGlueHDF5Files - ERR10'

        call SetDate(HDF5TimeInstant, Year     = TimeVector(1), Month  = TimeVector(2), &
                                      Day      = TimeVector(3), Hour   = TimeVector(4), &
                                      Minute   = TimeVector(5), Second = TimeVector(6))

                                     
        deallocate(TimeVector)

    end function HDF5TimeInstant

    
    !--------------------------------------------------------------------------    
    
    
    subroutine InquireFile(FileName)
    
        !Argument--------------------------------------------------------------
        character (*)                               :: Filename
        !Local-----------------------------------------------------------------
        logical                                     :: exist
        !Begin-----------------------------------------------------------------
    
        inquire(file = FileName, exist = exist)
        if (.not. exist) then
            write(*,*)
            write(*,*)'HDF file in list does not exist'
            write(*,*)trim(FileName)
            stop 'InquireFile - ModuleGlueHDF5Files - ERR01'
        endif
    
    end subroutine InquireFile

    !--------------------------------------------------------------------------
    subroutine CheckVGCompatibility(i)

        !Local-----------------------------------------------------------------
        integer                                     :: i, HDF5_READ, STAT_CALL, IDIn, IDOut
        logical                                     :: CheckOK

        !Begin-----------------------------------------------------------------
       
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (Me%ObjHDF5_In, Me%FileNameIn(i),                             &
                            Access = HDF5_READ, STAT = STAT_CALL)
       
        call GetHDF5FileID (Me%ObjHDF5_In, IDIn, STAT = STAT_CALL)

        call GetHDF5FileID (Me%ObjHDF5_Out, IDOut, STAT = STAT_CALL)

        CheckOK = .true.

        call CompareSubGroups (IDOut, IDIn, "/", "/", CheckOK)

        if (.not.CheckOK) then
            write(*,*) trim(Me%FileNameIn(i))//" is not a compatible file"
            stop 'CheckVGCompatibility - ModuleGlueHDF5Files - ERR10'
        endif
        
        call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Grid/VerticalZ", Me%Vert3D, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'CheckVGCompatibility - ModuleGlueHDF5Files - ERR20'


        call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Grid/OpenPoints", Me%Open3D, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'CheckVGCompatibility - ModuleGlueHDF5Files - ERR30'

        call KillHDF5(Me%ObjHDF5_In, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckVGCompatibility - ModuleGlueHDF5Files - ERR40'

        


    end subroutine CheckVGCompatibility

    !--------------------------------------------------------------------------
    
    subroutine CheckTimeAndBathymetry(i)

        !Local-----------------------------------------------------------------
        integer                                     :: i, HDF5_READ, STAT_CALL, IDIn, IDOut
        logical                                     :: CheckOK, OutTimeExist, InTimeExist
        logical                                     :: OutBathExist, InBathExist

        !Begin-----------------------------------------------------------------
       
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (Me%ObjHDF5_In, Me%FileNameIn(i),                             &
                            Access = HDF5_READ, STAT = STAT_CALL)
       
        call GetHDF5FileID (Me%ObjHDF5_In, IDIn, STAT = STAT_CALL)

        call GetHDF5FileID (Me%ObjHDF5_Out, IDOut, STAT = STAT_CALL)
        
        
        !Verify if time groups exist
        call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Time", OutTimeExist, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR10'

        call GetHDF5GroupExist(Me%ObjHDF5_In, "/Time", InTimeExist, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR20'

        if (.not. OutTimeExist .or. .not. InTimeExist) then
            write(*,*) "Time group is not available in one of the files"
            stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR030'
        endif

        CheckOK = .true.
        
        call CompareTime (IDOut, IDIn, "/Time", "/Time", CheckOK)

        if (.not.CheckOK) then
            write(*,*)
            write(*,*) trim(Me%FileNameIn(i))//" is not a compatible file"
            write(*,*) 'for merging HDFs because time instants do not match'
            stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR040'
        endif

        write (*,*)
        write (*,*) 'Checked time compatibility for file ', trim(Me%FileNameIn(i))


        !Verify if bathymetry exist
        call GetHDF5DataSetExist(Me%ObjHDF5_Out, "/Grid/Bathymetry", OutBathExist, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR50'

        call GetHDF5DataSetExist(Me%ObjHDF5_In, "/Grid/Bathymetry", InBathExist, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR60'

        if (.not. OutBathExist .or. .not. InBathExist) then
            write(*,*) "Bathymetry group is not available in one of the files"
            stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR070'
        endif

        CheckOK = .true.
        
        call CompareBathymetry (IDOut, IDIn, "/Grid/Bathymetry", "/Grid/Bathymetry", CheckOK)

        if (.not.CheckOK) then
            write(*,*)
            write(*,*) trim(Me%FileNameIn(i))//" is not a compatible file"
            write(*,*) 'for merging HDFs because bathymetries do not match'
            stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR080'
        endif

        write (*,*)
        write (*,*) 'Checked grid compatibility for file ', trim(Me%FileNameIn(i))


        call KillHDF5(Me%ObjHDF5_In, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckTimeAndBathymetry - ModuleGlueHDF5Files - ERR090'


    end subroutine CheckTimeAndBathymetry

    !-------------------------------------------------------------------------


    subroutine CompareTime (IDOut, IDIn, GroupNameOut, GroupNameIn, Check)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn
        character(len=*)                            :: GroupNameOut, GroupNameIn
        logical                                     :: Check


        !Local-----------------------------------------------------------------
        integer                                     :: nmembersOut, nmembersIn
        character(StringLength)                     :: obj_nameIn, obj_nameOut
        integer                                     :: obj_type, idx
        integer(HID_T)                              :: dset_id, gr_id
        integer(HSIZE_T), dimension(7)              :: dimsOut
        integer                                     :: STAT_CALL
        real, allocatable, dimension(:)             :: DataVal
        integer(HSIZE_T), dimension(7)              :: dims
        type (T_Time)                               :: InputTime, OutputTime


        !Get the number of members in the Group
        call h5gn_members_f(IDOut, GroupNameOut, nmembersOut, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR01'
    
        call h5gn_members_f(IDIn, GroupNameIn, nmembersIn, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR02'

        !if the number of time instants is different than file is not compatible
        if (nmembersOut /= nmembersIn) then
            Check = .false.
            return
        endif

        allocate(DataVal(6))

        
        do idx = 1, nmembersIn
            
            !Output
            !Gets information about the group
            call h5gget_obj_info_idx_f(IDOut, GroupNameOut, idx-1, obj_nameOut, obj_type, STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR03'

            !Opens the Group
            call h5gopen_f (IDOut, GroupNameOut, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CompareTime - ModuleHDF5Files - ERR04'

            !Opens data set
            call h5dopen_f (gr_id, trim(adjustl(obj_nameOut)), dset_id,  STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR05'

            call ReadInterface (dset_id, DataVal, dimsOut,  STAT_CALL)
            if (STAT_CALL/=0) stop 'CompareTime - ModuleHDF5Files - ERR06'

            call SetDate  (OutputTime, DataVal(1), DataVal(2), DataVal(3), DataVal(4), DataVal(5), DataVal(6))

            !Closes data set
            call h5dclose_f     (dset_id, STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR07'

            !Closes group
            call h5gclose_f     (gr_id, STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR08'
            
            
            
            !INPUT
            !Gets information about the group
            call h5gget_obj_info_idx_f(IDIn, GroupNameIn, idx-1, obj_nameIn, obj_type,  STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR09'

            !Opens the Group
            call h5gopen_f (IDIn, GroupNameIn, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'CompareTime - ModuleHDF5Files - ERR10'

            !Opens data set
            call h5dopen_f      (gr_id, trim(adjustl(obj_nameIn)), dset_id, STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR11'
            
            call ReadInterface (dset_id, DataVal, dims,  STAT_CALL)            
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR12'
            
            call SetDate  (InputTime, DataVal(1), DataVal(2), DataVal(3), DataVal(4), DataVal(5), DataVal(6))

            !Closes data set
            call h5dclose_f     (dset_id, STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR13'

            !Closes group
            call h5gclose_f     (gr_id, STAT_CALL)
            if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR14'

            !files do not have similar time instants
            if (InputTime /=  OutputTime) then
            
                Check = .false.
                return

            endif

        enddo

        deallocate(DataVal)

    end subroutine CompareTime

    !--------------------------------------------------------------------------

    subroutine CompareBathymetry (IDOut, IDIn, GroupNameOut, GroupNameIn, Check)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn
        character(len=*)                            :: GroupNameOut, GroupNameIn
        logical                                     :: Check


        !Local-----------------------------------------------------------------
        integer(HID_T)                              :: dset_id, gr_id, space_id
        integer(HSIZE_T), dimension(7)              :: dimsOut, dimsIn, maxdims
        integer                                     :: STAT_CALL
        character(StringLength)                     :: ParentGroupName


        !Output
        ParentGroupName = GroupNameOut(1:index(GroupNameOut, "/", .true.)-1)
        
        !Opens the Group
        call h5gopen_f (IDOut, ParentGroupName, gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CompareTime - ModuleHDF5Files - ERR04'

        !Opens data set
        call h5dopen_f (gr_id, GroupNameOut, dset_id,  STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR05'

        !Opens data space
        call h5dget_space_f (dset_id, space_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR06'
        
        !Gets dims
        call h5sget_simple_extent_dims_f  (space_id, dimsOut, maxdims, STAT_CALL) 
        !if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR06.5'
        
        !Closes data set
        call h5dclose_f     (dset_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR07'

        !Closes group
        call h5gclose_f     (gr_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR08'
        
        
        
        !INPUT
        ParentGroupName = GroupNameOut(1:index(GroupNameIn, "/", .true.)-1)
        
        !Opens the Group
        call h5gopen_f (IDIn, ParentGroupName, gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CompareTime - ModuleHDF5Files - ERR10'

        !Opens data set
        call h5dopen_f      (gr_id, GroupNameIn, dset_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR11'
        
        !Opens data space
        call h5dget_space_f (dset_id, space_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR012'
        
        !Gets dims
        call h5sget_simple_extent_dims_f  (space_id, dimsIn, maxdims, STAT_CALL) 
        !if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR012.5'

        !Closes data set
        call h5dclose_f     (dset_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR13'

        !Closes group
        call h5gclose_f     (gr_id, STAT_CALL)
        if (STAT_CALL /= 0) stop 'CompareTime - ModuleHDF5Files - ERR14'

        !files do not have similar dimension so are not compatible
        if (dimsIn(1) /= dimsOut(1)) then 
            Check = .false.
            return
        endif
        if (dimsIn(2) /= dimsOut(2)) then 
            Check = .false.
            return
        endif


    end subroutine CompareBathymetry

    !--------------------------------------------------------------------------

    subroutine CheckGroupExistence(i)
        
        !Arguments-------------------------------------------------------------
        integer                                     :: i
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_READ, IDIn
        !Begin-----------------------------------------------------------------
        
        !check if groups exist. if not add
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (Me%ObjHDF5_In, Me%FileNameIn(i),                             &
                            Access = HDF5_READ, STAT = STAT_CALL)
       
        call GetHDF5FileID (Me%ObjHDF5_In, IDIn, STAT = STAT_CALL)

        call AddNewGroupsToOutput(IDIn, '/'//trim(Me%BaseGroup))
        
        write (*,*) 
        write (*,*) 'Checked if new groups exist in file ', trim(Me%FileNameIn(i))

        
        call KillHDF5(Me%ObjHDF5_In, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'CheckGroupExistence - ModuleGlueHDF5Files - ERR02'

    end subroutine CheckGroupExistence
    
    !--------------------------------------------------------------------------

    recursive subroutine AddNewGroupsToOutput(IDin, GroupName)
        
        !Arguments
        integer                                     :: IDIn
        character(len=*)                            :: GroupName
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL, STAT
        integer                                     :: nmembersIn
        integer                                     :: idx, obj_type
        character(StringLength)                     :: obj_name, NewGroupName, NewGroupNameIn
        integer(HID_T)                              :: gr_idIn
        logical                                     :: Exist
        !Begin-----------------------------------------------------------------
       
        
        call h5gn_members_f(IDIn, GroupName, nmembersIn, STAT_CALL)
        if (STAT_CALL /= 0) stop 'AddNewGroupsToOutput - ModuleHDF5Files - ERR20'
        

        do idx = 1, nmembersIn

            !Gets information about the group
            call h5gget_obj_info_idx_f(IDIn, GroupName, idx-1, obj_name, obj_type,  STAT_CALL)
            if (STAT_CALL /= 0) stop 'AddNewGroupsToOutput - ModuleHDF5Files - ERR30'
            
            !if is a group verify if exists and if there are sub-groups
            if (obj_type ==H5G_GROUP_F) then
            
                !NewGroupName = GroupName//trim(adjustl(obj_name))//"/"
                NewGroupName = GroupName//"/"//trim(adjustl(obj_name))
                
                !check if the group exists in output. if not create it
                call GetHDF5GroupExist(Me%ObjHDF5_Out, NewGroupName, Exist, STAT = STAT_CALL)
                if (STAT_CALL /= 0) stop 'AddNewGroupsToOutput - ModuleGlueHDF5Files - ERR50'
                
                !Create group
                if (.not. Exist) then
                    call HDF5CreateGroup (Me%ObjHDF5_Out, NewGroupName, STAT = STAT_CALL)
                    if (STAT_CALL /= 0) stop 'AddNewGroupsToOutput - ModuleGlueHDF5Files - ERR60'
                    
                    
                endif
                
                !verify if new sub-groups exist
                !if (GroupName == "/") then
                !    NewGroupNameIn = GroupName//trim(adjustl(obj_name))
                !else
                    NewGroupNameIn = GroupName//"/"//trim(adjustl(obj_name))
                !endif                
                call h5gopen_f        (IDIn, trim(adjustl(NewGroupNameIn)), gr_idIn, STAT)
                call AddNewGroupsToOutput (gr_idIn, trim(adjustl(NewGroupNameIn)))
                call h5gclose_f       (gr_idIn, STAT)

            endif
            
        enddo


    end subroutine AddNewGroupsToOutput


    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    subroutine GlueFileIn(i)

        !Local-----------------------------------------------------------------
        integer                                     :: i, HDF5_READ, STAT_CALL
        integer(HID_T)                              :: IDIn, IDOut
        logical                                     :: CheckOK, Exist

        !Begin-----------------------------------------------------------------
       
        !Gets File Access Code
        call GetHDF5FileAccess  (HDF5_READ = HDF5_READ)

        call ConstructHDF5 (Me%ObjHDF5_In, Me%FileNameIn(i),                             &
                            Access = HDF5_READ, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR10'
       
        call GetHDF5FileID (Me%ObjHDF5_In, IDIn, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR20'

        call GetHDF5FileID (Me%ObjHDF5_Out, IDOut, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR30'
        
        CheckOK = .true.
        
        !Only glue time if glueing in time. if not time between files are the same
        if (Me%GlueInTime) then

            if (Me%TimeMap%BestTimeSerieON) then
                call GlueInTimeBest(IDOut        = IDOut,                               &
                                    IDIn         = IDIn,                                &
                                    GroupNameOut = "/Time",                             &
                                    GroupNameIn  = "/Time",                             &
                                    i            = i)  
            else                
            
                call GlueInTime (IDOut, IDIn, '/'//trim(Me%TimeGroup)//'/', '/'         &
                                                 //trim(Me%TimeGroup)//'/',             &
                                                 Me%FirstInstant(i), CheckOK)

            endif                        
            
            if (.not.CheckOK) then
                write(*,*) trim(Me%FileNameIn(i))//" is not a compatible file"
                stop 'GlueFileIn - ModuleGlueHDF5Files - ERR40'
            endif        
        endif
        
        if(Me%Vert3D)then
            if (Me%GlueInTime) then
                call GlueInResults (Me%ObjHDF5_Out, IDOut, IDIn, "/Grid/VerticalZ", Me%FirstInstant(i))
            else
                !only do this once if it does not exist in output (time instants are the same between files)
                !this will work for complete 3D files or mix of 2D and 3D
                call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Grid/VerticalZ", Exist, STAT = STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR50'
                
                if (.not. Exist) call GlueInResults (Me%ObjHDF5_Out, IDOut, IDIn, "/Grid/VerticalZ", Me%FirstInstant(i))

            endif
        endif
        
        if (Me%Open3D) then 
            if (Me%GlueInTime) then          
                call GlueInResults (Me%ObjHDF5_Out, IDOut, IDIn, "/Grid/OpenPoints", Me%FirstInstant(i), 2)
            else
                !only do this once if it does not exist in output (time instants are the same between files)
                !this will work for complete 3D files or mix of 2D and 3D
                call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Grid/OpenPoints", Exist, STAT = STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR60'
                
                if (.not. Exist) call GlueInResults (Me%ObjHDF5_Out, IDOut, IDIn, "/Grid/OpenPoints", Me%FirstInstant(i))            
            
            endif
        endif

        call GlueInResults (Me%ObjHDF5_Out, IDOut, IDIn, '/'//trim(Me%BaseGroup)//'/', Me%FirstInstant(i))

        call GetHDF5GroupExist(Me%ObjHDF5_Out, "/Generic4D/", Exist, STAT = STAT_CALL)
        if (STAT_CALL /= 0) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR70'

        if (Exist) then
            call GlueInResults (Me%ObjHDF5_Out, IDOut, IDIn, "/Generic4D/", Me%FirstInstant(i))
        endif

        if (.not.CheckOK) then
            write(*,*) trim(Me%FileNameIn(i))//" is not a compatible file"
            stop 'GlueFileIn - ModuleGlueHDF5Files - ERR80'
        endif


        call KillHDF5(Me%ObjHDF5_In, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueFileIn - ModuleGlueHDF5Files - ERR90'
       

    end subroutine GlueFileIn


    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    recursive subroutine CompareSubGroups (IDOut, IDIn, GroupNameOut, GroupNameIn, CheckOK)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn
        character(len=*)                            :: GroupNameOut, GroupNameIn
        logical                                     :: CheckOK


        !Local-----------------------------------------------------------------
        integer                                     :: nmembersOut, nmembersIn, nmembers
        character(StringLength)                     :: obj_nameIn, obj_nameOut
        integer                                     :: obj_type, idx
        integer(HID_T)                              :: gr_idIn, gr_idOut, dset_id
        integer(HID_T)                              :: space_id 
        integer                                     :: STAT
        character(StringLength)                     :: NewGroupNameIn, NewGroupNameOut
        integer(HSIZE_T), dimension(7)              :: dimsOut, dimsIn, maxdims


        !Get the number of members in the Group
        call h5gn_members_f(IDOut, GroupNameOut, nmembersOut, STAT)
        if (STAT /= SUCCESS_) return
    
        call h5gn_members_f(IDIn, GroupNameIn, nmembersIn, STAT)
        if (STAT /= SUCCESS_) return
        

        if (trim(GroupNameOut) /= trim(GroupNameIn)) then
            CheckOK = .false.
            return
        endif

        nmembers = min(nmembersOut, nmembersIn)

        do idx = 1, nmembers

            !Gets information about the group
            call h5gget_obj_info_idx_f(IDOut, GroupNameOut, idx-1, obj_nameOut, obj_type, STAT)
            if (STAT /= SUCCESS_) return
!            write(*,*)("+", i=1,Level),trim(adjustl(obj_name))
            call h5gget_obj_info_idx_f(IDIn, GroupNameIn, idx-1, obj_nameIn, obj_type, STAT)
            if (STAT /= SUCCESS_) return

            if     (obj_type == H5G_DATASET_F.and.trim(obj_nameOut)=="Bathymetry") then

                if (trim(obj_nameOut) /= trim(obj_nameIn))  then
                    CheckOK = .false.
                    return
                endif

                !Opens data set
                call h5dopen_f      (IDOut, trim(adjustl(obj_nameOut)), dset_id, STAT)

                !Opens data space
                call h5dget_space_f (dset_id, space_id, STAT)
    
                !Gets dims
                call h5sget_simple_extent_dims_f  (space_id, dimsOut, maxdims, STAT) 

                !Closes data space
                call h5sclose_f     (space_id, STAT)

                !Opens data set
                call h5dopen_f      (IDIn, trim(adjustl(obj_nameIn)), dset_id, STAT)

                !Opens data space
                call h5dget_space_f (dset_id, space_id, STAT)
    
                 !Gets dims
                call h5sget_simple_extent_dims_f  (space_id, dimsIn, maxdims, STAT) 

                !Closes data space
                call h5sclose_f     (space_id, STAT)


                if (dimsIn(1) /= dimsOut(1)) then 
                    CheckOK = .false.
                    return
                endif
                if (dimsIn(2) /= dimsOut(2)) then 
                    CheckOK = .false.
                    return
                endif


            elseif (obj_type ==H5G_GROUP_F) then

                if (trim(obj_nameOut) /= trim(obj_nameIn))  then
                    CheckOK = .false.
                    return
                endif

                !Looks for futher subgroups
                if (GroupNameOut == "/") then
                    NewGroupNameOut = GroupNameOut//trim(adjustl(obj_nameOut))
                else
                    NewGroupNameOut = GroupNameOut//"/"//trim(adjustl(obj_nameOut))
                endif
                if (GroupNameIn == "/") then
                    NewGroupNameIn = GroupNameIn//trim(adjustl(obj_nameIn))
                else
                    NewGroupNameIn = GroupNameIn//"/"//trim(adjustl(obj_nameIn))
                endif

                call h5gopen_f        (IDOut, trim(adjustl(NewGroupNameOut)), gr_idOut, STAT)
                call h5gopen_f        (IDIn, trim(adjustl(NewGroupNameIn)), gr_idIn, STAT)
                call CompareSubGroups (gr_idOut, gr_idIn, trim(adjustl(NewGroupNameOut)),&
                                                          trim(adjustl(NewGroupNameIn)), CheckOK)
                call h5gclose_f       (gr_idOut, STAT)
                call h5gclose_f       (gr_idIn, STAT)
            endif
            
        enddo


    end subroutine CompareSubGroups
    
    
    !--------------------------------------------------------------------------    

    subroutine GlueInTime (IDOut, IDIn, GroupNameOut, GroupNameIn, FirstInstant, Check)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn
        character(len=*)                            :: GroupNameOut, GroupNameIn
        integer                                     :: FirstInstant
        logical                                     :: Check


        !Local-----------------------------------------------------------------
        integer                                     :: nmembersOut, nmembersIn
        character(StringLength)                     :: obj_nameIn, obj_nameOut
        integer                                     :: obj_type, idx, NumType
        integer(HID_T)                              :: dset_id, prp_id, gr_id
        integer(HID_T)                              :: space_id 
        character(StringLength)                     :: Name
        integer(HSIZE_T), dimension(7)              :: dimsOut
        integer                                     :: Rank, STAT_CALL, k
        real, allocatable, dimension(:)             :: DataVal
        logical                                     :: FirstTime
        integer(HSIZE_T), dimension(7)              :: dims
        type (T_Time)                               :: NextTime


        !Get the number of members in the Group
        call h5gn_members_f(IDOut, GroupNameOut, nmembersOut, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR10'
    
        call h5gn_members_f(IDIn, GroupNameIn, nmembersIn, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR20'
        
        !Gets information about the group
        call h5gget_obj_info_idx_f(IDOut, GroupNameOut, nmembersOut-1, obj_nameOut, obj_type, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR30'

        !Opens the Group
        call h5gopen_f (IDOut, GroupNameOut, gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR40'

        !Opens data set
        call h5dopen_f (gr_id, trim(adjustl(obj_nameOut)), dset_id,  STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR50'
        
        allocate(DataVal(6))

        call ReadInterface (dset_id, DataVal, dimsOut,  STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR60'

        call SetDate  (Me%LastInstant, DataVal(1), DataVal(2), DataVal(3), DataVal(4), DataVal(5), DataVal(6))

        !Closes data set
        call h5dclose_f     (dset_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR70'

        !Closes group
        call h5gclose_f     (gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR80'


        FirstTime = .true. 

        k = 0

        do idx = 1, nmembersIn

            !Gets information about the group

            call h5gget_obj_info_idx_f(IDIn, GroupNameIn, idx-1, obj_nameIn, obj_type,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR90'

            !Opens the Group
            call h5gopen_f (IDIn, GroupNameIn, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR100'

            !Opens data set
            call h5dopen_f      (gr_id, trim(adjustl(obj_nameIn)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR110'
            
            call ReadInterface (dset_id, DataVal, dims,  STAT_CALL)            

            !call h5dread_f(dset_id, H5T_NATIVE_REAL, DataVal, dims,   STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR120'
            
            !Closes data set
            call h5dclose_f     (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR30'

            !Closes group
            call h5gclose_f     (gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR40'

            call SetDate  (NextTime, DataVal(1), DataVal(2), DataVal(3), DataVal(4), DataVal(5), DataVal(6))

            if (NextTime <  Me%LastInstant) then
            
                !stop 'GlueInTime - ModuleHDF5Files - ERR15'
                Check = .false.

            elseif (NextTime >  Me%LastInstant) then

                if (FirstTime) then 
                    FirstInstant = idx
                    FirstTime =.false.
                endif

                k = k + 1
                
                call ConstructDSName (trim(Me%TimeGroup), nmembersOut + k, Name)

                dims(1) = 6
                dims(2:7) = 0
                Rank = 1
                NumType = H5T_NATIVE_REAL

                !Opens Group, Creates Dset, etc
                call PrepareWrite (IDOut, Rank, dims, space_id, prp_id, gr_id,      &
                                   dset_id, NumType, GroupNameOut, trim(adjustl(Name)))

                call WriteInterface (dset_id, DataVal, dims,  STAT_CALL)                
                if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR160'
                
                !Closes data set
                call h5dclose_f     (dset_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'GlueInTime - ModuleGlueHDF5Files - ERR170'

                !Closes group
                call h5gclose_f     (gr_id, STAT_CALL)


            endif

        enddo


        deallocate(DataVal)

    end subroutine GlueInTime

    !--------------------------------------------------------------------------
    

    !--------------------------------------------------------------------------    

    subroutine GlueInTimeBest (IDOut, IDIn, GroupNameOut, GroupNameIn, i)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn
        character(len=*)                            :: GroupNameOut, GroupNameIn
        integer                                     :: i


        !Local-----------------------------------------------------------------
        integer                                     :: nmembersOut, nmembersIn
        character(StringLength)                     :: obj_nameIn
        integer                                     :: obj_type, idx, NumType
        integer(HID_T)                              :: dset_id, prp_id, gr_id
        integer(HID_T)                              :: space_id 
        character(StringLength)                     :: Name
        integer                                     :: Rank, STAT_CALL, k
        real, allocatable, dimension(:)             :: DataVal
        integer(HSIZE_T), dimension(7)              :: dims
        type (T_Time)                               :: NextTime
        
        !Begin-----------------------------------------------------------------
        
        Me%TimeMap%PresentFile = i


        !Get the number of members in the Group
        call h5gn_members_f(IDOut, GroupNameOut, nmembersOut, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR10'
        
        call h5gn_members_f(IDIn, GroupNameIn, nmembersIn, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR20'
        
        allocate(DataVal(6))

        k = 0

        do idx = Me%TimeMap%FirstInstantBest(i), Me%TimeMap%LastInstantBest(i)

            !Gets information about the group
        
            if (idx == 0) cycle

            call h5gget_obj_info_idx_f(IDIn, GroupNameIn, idx-1, obj_nameIn, obj_type,  STAT_CALL)
            if (STAT_CALL /= SUCCESS_) then
                write(*,*) "FileName =", Me%FileNameIn(i)            
                write(*,*) "GroupNameIn =", GroupNameIn                         
                write(*,*) "obj_nameIn =", obj_nameIn                
                stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR30'
            endif                

            !Opens the Group
            call h5gopen_f (IDIn, GroupNameIn, gr_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR40'

            !Opens data set
            call h5dopen_f      (gr_id, trim(adjustl(obj_nameIn)), dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR50'
            
            call ReadInterface (dset_id, DataVal, dims,  STAT_CALL)            

            !call h5dread_f(dset_id, H5T_NATIVE_REAL, DataVal, dims,   STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR60'

            !Closes data set
            call h5dclose_f     (dset_id, STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR70'

            !Closes group
            call h5gclose_f     (gr_id, STAT_CALL)

            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR80'


            call SetDate  (NextTime, DataVal(1), DataVal(2), DataVal(3), DataVal(4), DataVal(5), DataVal(6))

            k = k + 1
                
            call ConstructDSName (trim(Me%TimeGroup), nmembersOut + k, Name)

            dims(1) = 6
            dims(2:7) = 0
            Rank = 1
            NumType = H5T_NATIVE_REAL

            !Opens Group, Creates Dset, etc
            call PrepareWrite (IDOut, Rank, dims, space_id, prp_id, gr_id,      &
                                dset_id, NumType, GroupNameOut, trim(adjustl(Name)))

            call WriteInterface (dset_id, DataVal, dims,  STAT_CALL)                
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR90'
            
            !Closes data set
            call h5dclose_f     (dset_id, STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'GlueInTimeBest - ModuleGlueHDF5Files - ERR100'

            !Closes group
            call h5gclose_f     (gr_id, STAT_CALL)

        enddo


        deallocate(DataVal)

    end subroutine GlueInTimeBest

    !--------------------------------------------------------------------------
    

    recursive subroutine GlueInResults (ObjHDF5_Out, IDOut, IDIn, GroupName, FirstInstant, dataType)

        !Arguments-------------------------------------------------------------
        integer(HID_T)                              :: IDOut, IDIn, ObjHDF5_Out
        character(len=*)                            :: GroupName
        integer                                     :: FirstInstant
        integer, optional                           :: dataType



        !Local-----------------------------------------------------------------
        integer                                     :: nmembersOut, nmembersIn
        character(StringLength)                     :: obj_name
        integer                                     :: obj_type, idx
        integer(HID_T)                              :: gr_idIn, gr_idOut
        integer(HID_T)                              :: dset_id, gr_id
        integer(HID_T)                              :: space_id 
        character(StringLength)                     :: NewGroupName
        integer(HSIZE_T), dimension(7)              :: dims, maxdims
        integer, dimension(7)                       :: dims_int
        integer                                     :: Rank, STAT_CALL, k, ia
        real, pointer, dimension(:)                 :: DataVal1D
        real, pointer, dimension(:,:)               :: DataVal2D
        integer, pointer, dimension(:,:)            :: DataInt2D
        real, pointer, dimension(:,:,:)             :: DataVal3D
        integer, pointer, dimension(:,:,:)          :: DataInt3D
        integer(HID_T)                              :: attr_id, type_id
        character(len=StringLength)                 :: Units
        logical                                     :: data_is_integer
        integer                                     :: istart, iend

        !Begin-----------------------------------------------------------------        



        !Get the number of members in the Group
        call h5gn_members_f(IDOut, GroupName, nmembersOut, STAT_CALL)

        if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR10'
    
        call h5gn_members_f(IDIn, GroupName, nmembersIn, STAT_CALL)

        if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR20'
        

        k = 0
        

        do idx = 1, nmembersIn

            !Gets information about the group

            call h5gget_obj_info_idx_f(IDIn, GroupName, idx-1, obj_name, obj_type,  STAT_CALL)

            if (STAT_CALL /= 0) then
                write(*,*) "Groupname =", GroupName
                write(*,*) "obj_name =", obj_name
                stop 'GlueInResults - ModuleHDF5Files - ERR30'
            endif                

            if (obj_type == H5G_DATASET_F) then
            
                if (Me%TimeMap%BestTimeSerieON) then
                    istart  = Me%TimeMap%FirstInstantBest(Me%TimeMap%PresentFile)
                    iend    = Me%TimeMap%LastInstantBest (Me%TimeMap%PresentFile)
                    if (idx < istart .or. idx > iend) then
                        cycle
                    endif                
                endif            


                !Opens the Group
                call h5gopen_f (IDIn, GroupName, gr_id, STAT_CALL)
                if (STAT_CALL /= SUCCESS_) stop 'HDF5ReadDataR4 - ModuleHDF5 - ERR40'

                !Opens data set
                call h5dopen_f      (gr_id, trim(adjustl(obj_name)), dset_id, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR50'

                !Opens data space
                call h5dget_space_f (dset_id, space_id, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR60'

                !Gets dims
                call h5sget_simple_extent_dims_f  (space_id, dims, maxdims, Rank) 
                if (Rank < 0) stop 'GlueInResults - ModuleHDF5Files - ERR70'
                
                if     (Rank==1) then

                    allocate(DataVal1D(1:dims(1)))

                    call ReadInterface(dset_id, DataVal1D, dims,   STAT_CALL)

                    if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR80'
               
                
                elseif (Rank==2) then

                    allocate(DataVal2D(1:dims(1),1:dims(2)))

                    call ReadInterface(dset_id, DataVal2D, dims,   STAT_CALL)
                    
                    if (STAT_CALL /= 0) then
                        
                        data_is_integer = .true.
                        
                        deallocate(DataVal2D)
                        
                        allocate(DataInt2D(1:dims(1),1:dims(2)))
                        
                        call ReadInterface(dset_id, DataInt2D, dims,   STAT_CALL)
                        
                    else
                    
                        data_is_integer = .false. 
                        
                    endif
                        
                    if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR90'


                elseif(Rank == 3) then
                    
                    if ( present(dataType) ) then
                        if ( dataType .eq. 2 ) then
                            allocate(DataInt3D(1:dims(1),1:dims(2),1:dims(3)))
                            call ReadInterface(dset_id, DataInt3D, dims,   STAT_CALL)
                            if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR110'
                        endif
                    endif
                    
                    if ( .not. present(dataType) .or. present(dataType) .and. dataType .eq. 1 ) then
                        allocate(DataVal3D(1:dims(1),1:dims(2),1:dims(3)))
                        call ReadInterface(dset_id, DataVal3D, dims,   STAT_CALL)
                        if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR110'
                    endif
                    

                endif


                
                !Closes data space
                call h5sclose_f     (space_id, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR120'


                !Reads Units
                call h5aopen_name_f     (dset_id, "Units", attr_id, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR130'

                call h5Tcopy_f          (H5T_NATIVE_CHARACTER, type_id, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR140'

                call h5Tset_size_f      (type_id, StringLength, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR150'

                call h5aread_f          (attr_id, type_id, Units, dims, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR160'
                
                call h5aclose_f         (attr_id, STAT_CALL) 
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR170'

                call h5Tclose_f         (type_id, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR180'
                !Closes data set
                call h5dclose_f     (dset_id, STAT_CALL)
                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR190'

                !Closes group
                call h5gclose_f     (gr_id, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR200'


                if (idx >= FirstInstant) then

                    k = k + 1

                    do ia=len_trim(obj_name),1,-1
                        if (obj_name(ia:ia) == '_') exit
                    enddo                  

                    dims_int = dims

                    if (Rank==1) then

                        call HDF5SetLimits  (ObjHDF5_Out, 1, dims_int(1), STAT = STAT_CALL)
                        if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR210'
                    
                        call HDF5WriteData(ObjHDF5_Out, GroupName, trim(adjustl(obj_name(1:ia-1))), & 
                                           Units, Array1D = DataVal1D, OutputNumber = nmembersOut + k, STAT = STAT_CALL)
                        if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR220'

                        deallocate(DataVal1D)

                    elseif (Rank==2) then

                        call HDF5SetLimits  (ObjHDF5_Out, 1, dims_int(1), 1, dims_int(2), STAT = STAT_CALL)
                        if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR220'
                    
                        if (data_is_integer) then
                        
                            call HDF5WriteData(ObjHDF5_Out, GroupName, trim(adjustl(obj_name(1:ia-1))), & 
                                               Units, Array2D = DataInt2D, OutputNumber = nmembersOut + k, STAT = STAT_CALL)
                            if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR225'

                            deallocate(DataInt2D)
                        else                    
                            call HDF5WriteData(ObjHDF5_Out, GroupName, trim(adjustl(obj_name(1:ia-1))), & 
                                               Units, Array2D = DataVal2D, OutputNumber = nmembersOut + k, STAT = STAT_CALL)
                            if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR230'

                            deallocate(DataVal2D)
                        endif

                    elseif(Rank == 3) then

                        call HDF5SetLimits  (ObjHDF5_Out, 1, dims_int(1), 1, dims_int(2),1, dims_int(3),  STAT = STAT_CALL)

                        if ( present(dataType) ) then
                            if ( dataType .eq. 2 ) then
                                call HDF5WriteData(ObjHDF5_Out, GroupName, trim(adjustl(obj_name(1:ia-1))), & 
                                           Units, Array3D = DataInt3D, OutputNumber = nmembersOut + k, STAT = STAT_CALL)
                                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR240'
                                deallocate(DataInt3D)
                            endif
                        endif
                    
                        if ( .not. present(dataType) .or. present(dataType) .and. dataType .eq. 1 ) then
                            call HDF5WriteData(ObjHDF5_Out, GroupName, trim(adjustl(obj_name(1:ia-1))), & 
                                           Units, Array3D = DataVal3D, OutputNumber = nmembersOut + k, STAT = STAT_CALL)
                            if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR240'
                            deallocate(DataVal3D)
                        endif


                    endif

                endif

            elseif (obj_type ==H5G_GROUP_F) then

                NewGroupName = GroupName//trim(adjustl(obj_name))//"/"

                call h5gopen_f        (IDOut, trim(adjustl(NewGroupName)), gr_idOut, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR350'

                call h5gopen_f        (IDIn, trim(adjustl(NewGroupName)), gr_idIn, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR360'

                call GlueInResults    (ObjHDF5_Out, gr_idOut, gr_idIn, trim(adjustl(NewGroupName)), FirstInstant)
                call h5gclose_f       (gr_idOut, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR370'

                call h5gclose_f       (gr_idIn, STAT_CALL)

                if (STAT_CALL /= 0) stop 'GlueInResults - ModuleHDF5Files - ERR1380'

            endif

        enddo

    end subroutine GlueInResults

    subroutine WriteInterface1DR4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(4),   dimension(:)                     :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dWrite_f(dset_id, H5T_NATIVE_REAL, DataVal, dims,  STAT_CALL)

    end subroutine WriteInterface1DR4 


    subroutine WriteInterface1DR8 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(8),   dimension(:)                     :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dWrite_f(dset_id, H5T_NATIVE_DOUBLE, DataVal, dims,  STAT_CALL)

    end subroutine WriteInterface1DR8 



    subroutine ReadInterface1DR4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(4),   dimension(:)                     :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_REAL, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface1DR4 


    subroutine ReadInterface1DR8 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(8),   dimension(:)                     :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface1DR8 



    subroutine ReadInterface2DR4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(4),   dimension(:,:)                   :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_REAL, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface2DR4 


    subroutine ReadInterface2DR8 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(8),   dimension(:,:)                   :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface2DR8 


    subroutine ReadInterface3DR4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(4),   dimension(:,:,:)                 :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_REAL, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface3DR4 


    subroutine ReadInterface3DR8 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        real(8),   dimension(:,:,:)                 :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface3DR8 


    subroutine ReadInterface1DI4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        integer,   dimension(:)                     :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface1DI4 

    subroutine ReadInterface2DI4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        integer,   dimension(:,:)                   :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface2DI4 


    subroutine ReadInterface3DI4 (dset_id, DataVal, dims,  STAT_CALL)

        !Arguments----------------------------------------------------------------
        integer(HID_T)                              :: dset_id
        integer(HSIZE_T), dimension(7)              :: dims
        integer                                     :: STAT_CALL
        integer,   dimension(:,:,:)                 :: DataVal

        !Begin--------------------------------------------------------------------
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, DataVal, dims,  STAT_CALL)

    end subroutine ReadInterface3DI4 



    

    subroutine KillGlueHDF5Files
        
        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers
        
        !Begin-----------------------------------------------------------------
        
        if (Me%TimeMap%BestTimeSerieON) then
        
            deallocate(Me%TimeMap%FirstInstantBest)        
            deallocate(Me%TimeMap%LastInstantBest )                

            deallocate(Me%TimeMap%FirstDateBest   )        
            deallocate(Me%TimeMap%LastDateBest    )   
                    
            deallocate(Me%TimeMap%ObjHDF5_ID      )
            deallocate(Me%TimeMap%File_ID         )     

        endif
        
        if (Me%CheckHDF5_File) then
        
            call GetHDF5AllDataSetsOK (HDF5ID = Me%ObjHDF5_Out, STAT = STAT_CALL)                                      
            if (STAT_CALL /= SUCCESS_) stop 'KillGlueHDF5Files - ModuleGlueHDF5Files - ERR10'
            
        endif              
        
        call KillHDF5(Me%ObjHDF5_Out, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillGlueHDF5Files - ModuleGlueHDF5Files - ERR20'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillGlueHDF5Files - ModuleGlueHDF5Files - ERR30'


        deallocate(Me%FileNameIn)
        deallocate(Me)
        nullify   (Me)

    
    end subroutine KillGlueHDF5Files

    !--------------------------------------------------------------------------

    subroutine ConstructDSName (Name, OutputNumber, AuxChar)

        !Arguments-------------------------------------------------------------
        character(len=*)                            :: Name
        integer                                     :: OutputNumber
        character(len=*)                            :: AuxChar

        !Local-----------------------------------------------------------------
        character(StringLength)                     :: AuxNum

        write(AuxNum, fmt=*)OutputNumber

        if     (OutputNumber < 10     ) then
            AuxChar = trim(adjustl(Name))//"_0000"//trim(adjustl(AuxNum))
        elseif (OutputNumber < 100    ) then
            AuxChar = trim(adjustl(Name))//"_000" //trim(adjustl(AuxNum))
        elseif (OutputNumber < 1000   ) then
            AuxChar = trim(adjustl(Name))//"_00"  //trim(adjustl(AuxNum))
        elseif (OutputNumber < 10000  ) then
            AuxChar = trim(adjustl(Name))//"_0"   //trim(adjustl(AuxNum))
        else
            AuxChar = trim(adjustl(Name))//"_"    //trim(adjustl(AuxNum))
        endif

    end subroutine ConstructDSName

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine PrepareWrite (FileID, Rank, dims, space_id, prp_id, gr_id, dset_id,       &
                             NumType, GroupName, ItemName)

        !Arguments-------------------------------------------------------------
        integer (HID_T)                             :: FileID, Rank, space_id
        integer (HID_T)                             :: prp_id, gr_id, dset_id, NumType
        integer (HSIZE_T), dimension(7)             :: dims
        character(len=*)                            :: GroupName
        character(len=*)                            :: ItemName

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL

        !Creates a simple dataspace
        call h5screate_simple_f(Rank, dims, space_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleGlueHDF5Files - ERR01'

        !Creates a property list
        call h5pcreate_f (H5P_DATASET_CREATE_F, prp_id, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleGlueHDF5Files - ERR02'

        !Sets chunked
        call h5pset_chunk_f(prp_id, Rank, dims, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleGlueHDF5Files - ERR03'

        !Sets the compression
        call h5pset_deflate_f(prp_id, 6, STAT_CALL) 
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleGlueHDF5Files - ERR04'

        !Opens the Group
        call h5gopen_f (FileID, GroupName, gr_id, STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleGlueHDF5Files - ERR05'

        !Creates the dataset with default properties
        call h5dcreate_f(gr_id, ItemName, NumType, space_id, dset_id,  STAT_CALL, prp_id) 
        if (STAT_CALL /= SUCCESS_) stop 'PrepareWrite - ModuleGlueHDF5Files - ERR06'

    end subroutine PrepareWrite



    !--------------------------------------------------------------------------
 
    !--------------------------------------------------------------------------
end module ModuleGlueHDF5Files










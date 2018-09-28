program MohidDomainDecomposition
    
    use ModuleGlobalData
    use ModuleEnterData
    use ModuleHorizontalGrid
    use ModuleGridData  
    use ModuleHorizontalMap
    use ModuleTime
    use ModuleGeometry
    use ModuleMap
    use ModuleFunctions, only: isodd

    implicit none
    
    !Modules
    integer                                     :: ObjHorizontalGrid = 0
    integer                                     :: ObjBathymetry     = 0
    integer                                     :: ObjHorizontalMap  = 0
    integer                                     :: ObjGeometry       = 0
    integer                                     :: ObjMap            = 0

    !Size
    type (T_Size2D)                             :: WorkSize

    !Files
    character(len=PathLength)                   :: BathymetryFile
    character(len=PathLength)                   :: GeometryFile
    character(len=PathLength)                   :: OutputFile             

    !Global variables
    integer                                     :: nSectionsX       = null_int
    integer                                     :: nSectionsY       = null_int
    logical                                     :: HasMaster        = .false.
    logical                                     :: UseCenterMass    = .false.
    integer                                     :: DivisionMethod   = null_int
    
    integer, parameter                          :: Grid             = 1
    integer, parameter                          :: HorizontalStripes= 2
    integer, parameter                          :: VerticalStripes  = 3
    
    !Mapping
    integer, dimension(:,:  ), pointer          :: WaterPoints2D
    integer, dimension(:,:,:), pointer          :: WaterPoints3D
    real,    dimension(:,:  ), pointer          :: Bathymetry
    
    call ConstructMohidDomainDecomposition
    call ModifyMohidDomainDecomposition
    call KillMohidDomainDecomposition

    contains
    
    !--------------------------------------------------------------------------

    subroutine ConstructMohidDomainDecomposition
    
        call StartUpMohid("MohidDomainDecomposition")

        call ReadKeywords
        
        call ConstructDomain

    end subroutine ConstructMohidDomainDecomposition
    
    !--------------------------------------------------------------------------

    subroutine ModifyMohidDomainDecomposition
        
        !Local-----------------------------------------------------------------
        integer                                     :: i, j, n, s, Unit, STAT_CALL
        integer                                     :: ILB, IUB, JLB, JUB
        integer                                     :: dILB, dIUB, dJLB, dJUB
        integer                                     :: nDomains         = 0
        integer                                     :: nInterfaces      = 0
        integer                                     :: nInterfacesSN    = 0
        integer                                     :: ninterfacesWE    = 0
        integer                                     :: MasterDomainID   = 0
        real                                        :: DX, RX, DY, RY
        integer                                     :: start_i, end_i
        real                                        :: nWaterPoints, TotalWaterPoints
        real                                        :: SubdomainPointsEstimation, SubdomainPoints, SubdomainPointsSoFar
        integer, dimension(:),   pointer            :: DomainsID    
        integer, dimension(:),   pointer            :: dimI, dimJ
        integer, dimension(:,:), pointer            :: Ibreaks, Jbreaks, DomainsID2D
        integer                                     :: centermassI, centermassJ
        real                                        :: totalmass, centermassX, centermassY
        real                                        :: maxcourant
        character(len=PathLength)                   :: frmt
        
        ILB = WorkSize%ILB 
        JLB = WorkSize%JLB 
        IUB = WorkSize%IUB
        JUB = WorkSize%JUB
        
        write(*,*)
        write(*,*)"Calculating total number of computation points "
        TotalWaterPoints = sum(WaterPoints3D(ILB:IUB, JLB:JUB, :))     
     
        write(*,*)"     TOTAL  : ", TotalWaterPoints
        
        write(*,*)
        write(*,*)"Calculating 3D grid center of mass "
        
        totalmass   = 0.0
        centermassX = 0.0
        centermassY = 0.0
        
        do i = ILB, IUB
        do j = JLB, JUB
            if(WaterPoints2D(i,j) == WaterPoint)then
                
                totalmass   = totalmass + sum(WaterPoints3D(i,j,:))
                
                centermassX = centermassX + float(j) * sum(WaterPoints3D(i,j,:))
                centermassY = centermassY + float(i) * sum(WaterPoints3D(i,j,:))
                
            endif
        enddo
        enddo
        
        centermassJ = int(centermassX/totalmass)
        centermassI = int(centermassY/totalmass)
        
        write(*,*)"     CENTER OF MASS (I,J) : ", centermassI,  centermassJ
        
        write(*,*)
        write(*,*)"Calculating domain decomposition"

        nullify(dimI, dimJ)

        allocate(dimI(1:2))
        allocate(dimJ(1:2))
        dimI(1) = ILB    
        dimI(2) = IUB   
        dimJ(1) = JLB    
        dimJ(2) = JUB
    
        allocate(Ibreaks(1:nSectionsY,2))
        allocate(Jbreaks(1:nSectionsX,2))
    
        if(HasMaster)then
            MasterDomainID = 1
        else
            MasterDomainID = 0
        endif
    
        nDomains         = nSectionsX    *  nSectionsY
        select case(DivisionMethod)
            case(Grid)
                nInterfacesSN    = (nSectionsX-1)*  nSectionsY
                ninterfacesWE    = nSectionsX    * (nSectionsY-1)
            case(HorizontalStripes)
                nInterfacesSN    = nSectionsY - 1
                ninterfacesWE    = 0
            case(VerticalStripes)
                nInterfacesSN    = 0
                ninterfacesWE    = nSectionsX - 1
            case default
                stop 'ModifyMohidDomainDecomposition - MohidDomainDecomposition - ERR01'
        end select
        
        nInterfaces      = nInterfacesSN + ninterfacesWE

    
        nullify(DomainsID, DomainsID2D)
        allocate(DomainsID(1:nDomains))
        allocate(DomainsID2D(1:nSectionsY, 1:nSectionsX))
    
        do i = MasterDomainID, nDomains + MasterDomainID - 1
            DomainsID(i) = i
        enddo
        
        call UnitsManager(Unit, OPEN_FILE, STAT = STAT_CALL)
    
        open(Unit = Unit, File = trim(outputfile), Status = 'replace')
        
        write(frmt, "(I3, A)") nDomains, "(I3, 1x)"
    
        call WriteDataLine(Unit, "SUBDOMAINS_NUMBER", nDomains)
        write(unit, "(A25, "//trim(frmt)//"))")"SUBDOMAINS              : ", DomainsID
        call WriteDataLine(Unit, "INTERFACES_NUMBER", nInterfaces)
        call WriteDataLine(Unit, "ILB_IUB_ALL",       2, dimI)
        call WriteDataLine(Unit, "ILB_IUB_ALL",       2, dimJ)
    
        n = MasterDomainID
        do i = 1, nSectionsY
        do j = 1, nSectionsX
            DomainsID2D(i, j) = n
            n = n + 1
        enddo
        enddo 
        
        if(UseCenterMass)then
            
            if(isodd(nSectionsX) .or. isodd(nSectionsY)     .or.    &
               .not. (iand (nSectionsX, nSectionsX-1) == 0) .or.    &
               .not. (iand (nSectionsY, nSectionsY-1) == 0))then
                stop 'N_DOMAINS_X and N_DOMAINS_Y must be power of 2 (e.g. 2,4,6,...)'
            endif

                
            if(nSectionsX == 2 .and. nSectionsY == 2)then
                    
                Jbreaks(1, 1) = JLB
                Jbreaks(1, 2) = centermassJ
                Jbreaks(2, 1) = centermassJ+1
                Jbreaks(2, 2) = JUB
                    
                Ibreaks(1, 1) = ILB
                Ibreaks(1, 2) = centermassI
                Ibreaks(2, 1) = centermassI+1
                Ibreaks(2, 2) = IUB

            else 
                
                !compute DX for left section
                DX = centermassJ / (nSectionsX/2.0)
                RX = MOD(float(centermassJ), (nSectionsX/2.0)) 

                start_i = JLB
                end_i   = int(DX) 
    
                do i = 1, nSectionsX/2
                    Jbreaks(i, 1) = start_i
                    Jbreaks(i, 2) = end_i
                    start_i       = end_i + 1
                    end_i         = end_i + int(DX)
                enddo
                
                if(RX > 0.0)then
                    Jbreaks(i, 2) = centermassJ
                endif

                !compute DX for right section
                DX = (JUB-centermassJ) / (nSectionsX/2.0)
                RX = MOD(float(JUB-centermassJ), (nSectionsX/2.0)) 

                start_i = centermassJ
                end_i   = centermassJ + int(DX) 
    
                do i = nSectionsX/2+1, nSectionsX
                    Jbreaks(i, 1) = start_i
                    Jbreaks(i, 2) = end_i
                    start_i       = end_i + 1
                    end_i         = end_i + int(DX)
                enddo
                
                if(RX > 0.0)then
                    Jbreaks(nSectionsX, 2) = JUB
                endif

                
                !compute DY for bottom section
                DY = centermassI / (nSectionsY/2.0)
                RY = MOD(float(centermassI), nSectionsY/2.0) 

                start_i = ILB
                end_i   = int(DY) 
    
                do i = 1, nSectionsY/2
                    Ibreaks(i, 1) = start_i
                    Ibreaks(i, 2) = end_i
                    start_i       = end_i + 1
                    end_i         = end_i + int(DY)
                enddo
                
                if(RY > 0.0)then
                    Ibreaks(i, 2) = centermassI
                endif

                !compute DX for right section
                DY = (IUB-centermassI) / (nSectionsY/2.0)
                RY = MOD(float(IUB-centermassI), nSectionsY/2.0) 

                start_i = centermassI
                end_i   = centermassI + int(DY) 
    
                do i = nSectionsY/2+1, nSectionsY
                    Ibreaks(i, 1) = start_i
                    Ibreaks(i, 2) = end_i
                    start_i       = end_i + 1
                    end_i         = end_i + int(DY)
                enddo
                
                if(RY > 0.0)then
                    Ibreaks(nSectionsY, 2) = IUB
                endif
                    
            endif
            
        else
            
            selectcase(DivisionMethod)
                
                case(Grid)
                
                    DX = JUB / nSectionsX
                    DY = IUB / nSectionsY
                    RX = MOD(JUB, nSectionsX) 
                    RY = MOD(IUB, nSectionsY)
    
                    start_i = 1
                    end_i   = int(DX) 
    
                    do i = 1, nSectionsX
                        Jbreaks(i, 1) = start_i
                        Jbreaks(i, 2) = end_i
                        start_i       = end_i + 1
                        end_i         = end_i + int(DX)
                    enddo
    
                    start_i = 1
                    end_i   = int(DY) 
    
                    if(RX > 0.0)then
                        Jbreaks(nSectionsX, 2) = JUB
                    endif
    
                    do i = 1, nSectionsY
                        Ibreaks(i, 1) = start_i
                        Ibreaks(i, 2) = end_i
                        start_i       = end_i + 1
                        end_i         = end_i + int(DY)
                    enddo
    
                    if(RY > 0.0)then
                        Ibreaks(nSectionsY, 2) = IUB
                    endif
                
                case(HorizontalStripes)
                    
                    SubdomainPointsEstimation = TotalWaterPoints / float(nDomains)
                    
                    Jbreaks(1,1)         = JLB
                    Jbreaks(1,2)         = JUB

                    SubdomainPoints      = 0
                    SubdomainPointsSoFar = 0
                    start_i              = ILB
                    n                    = 1
                    
                    do i = ILB, IUB

                        SubdomainPoints      = SubdomainPoints      + sum(WaterPoints3D(i,:,:))
                        
                        
                        if(SubdomainPoints .ge. SubdomainPointsEstimation)then
                            
                            SubdomainPointsSoFar      = SubdomainPointsSoFar + SubdomainPoints
                            SubdomainPointsEstimation = (TotalWaterPoints - SubdomainPointsSoFar) / float(nDomains - n)
                            
                            Ibreaks(n,1) = start_i
                            Ibreaks(n,2) = i
                            
                            start_i      = i + 1
                            
                            SubdomainPoints = 0
                            n            = n + 1
                            
                        endif
                        
                        Ibreaks(nDomains,1) = Ibreaks(nDomains-1,2)+1
                        Ibreaks(nDomains,2) = IUB

                    enddo
                    
                case(VerticalStripes)
                    
                    SubdomainPointsEstimation = TotalWaterPoints / float(nDomains)
                    
                    Ibreaks(1,1)         = ILB
                    Ibreaks(1,2)         = IUB

                    SubdomainPoints      = 0
                    SubdomainPointsSoFar = 0
                    start_i              = JLB
                    n                    = 1
                    
                    do i = JLB, JUB

                        SubdomainPoints      = SubdomainPoints      + sum(WaterPoints3D(:,i,:))
                        
                        
                        if(SubdomainPoints .ge. SubdomainPointsEstimation)then
                            
                            SubdomainPointsSoFar      = SubdomainPointsSoFar + SubdomainPoints
                            SubdomainPointsEstimation = (TotalWaterPoints - SubdomainPointsSoFar) / float(nDomains - n)
                            
                            Jbreaks(n,1) = start_i
                            Jbreaks(n,2) = i
                            
                            start_i      = i + 1
                            
                            SubdomainPoints = 0
                            n            = n + 1
                            
                        endif
                        
                        Jbreaks(nDomains,1) = Jbreaks(nDomains-1,2)+1
                        Jbreaks(nDomains,2) = JUB

                    enddo

                    
            endselect
            
            
            
            
        endif

100     format(A, i3, 1x, i3)
    
        n = MasterDomainID
        do i = 1, nSectionsY
        do j = 1, nSectionsX
            call WriteDataLine(Unit, '<BeginSubDD>')
            call WriteDataLine(Unit, 'MPI_ID', n)
            write(Unit, 100)'ILB_IUB                 : ', Ibreaks(i,1), Ibreaks(i,2)
            write(Unit, 100)'JLB_JUB                 : ', Jbreaks(j,1), Jbreaks(j,2)
            call WriteDataLine(Unit, '<EndSubDD>')
            n = n + 1
        enddo
        enddo

200     format(i3, 1x, i3)
        !Interfaces
        call WriteDataLine(Unit, '<BeginInterfaceSN>')
        do j = 1, nSectionsX
        do i = 1, nSectionsY - 1
            write(Unit, 200)DomainsID2D(i,j), DomainsID2D(i+1, j)
        enddo
        enddo
        call WriteDataLine(Unit, '<EndInterfaceSN>')
    
        call WriteDataLine(Unit, '<BeginInterfaceWE>')
        do i = 1, nSectionsY
        do j = 1, nSectionsX - 1
            write(Unit, 200)DomainsID2D(i,j), DomainsID2D(i, j+1)
        enddo
        enddo
        call WriteDataLine(Unit, '<EndInterfaceWE>')

        call UnitsManager(Unit, CLOSE_FILE, STAT = STAT_CALL)
        
        write(*,*)
        write(*,*)"Calculating sub-domains number of computation points "
        
        n = MasterDomainID
        !Analyze domains
        do i = 1, nSectionsY
        do j = 1, nSectionsX
            dILB = Ibreaks(i,1)
            dIUB = Ibreaks(i,2)
            dJLB = Jbreaks(j,1)
            dJUB = Jbreaks(j,2)
            
            nWaterPoints = sum(WaterPoints3D(dILB:dIUB, dJLB:dJUB, :))
            !maxcourant   = sqrt(Gravity*maxval(Bathymetry(dILB:dIUB, dJLB:dJUB)))*0.02
                        
            write(*,*)"     DOMAIN : ", n, " --> ", nWaterPoints/TotalWaterPoints*100.0
            !write(*,*)"            : ", n, " --> ", maxcourant
            
            n = n + 1
        enddo
        enddo
    
        deallocate(dimI)
        deallocate(dimJ)
        deallocate(Ibreaks)
        deallocate(Jbreaks)
        deallocate(DomainsID)
        deallocate(DomainsID2D)
    
    end subroutine ModifyMohidDomainDecomposition
   
    !--------------------------------------------------------------------------

    subroutine KillMohidDomainDecomposition

        call ShutdownMohid ("MohidDomainDecomposition", 1.0, 1.0)

    end subroutine KillMohidDomainDecomposition
    
    !--------------------------------------------------------------------------
    
    subroutine ReadKeywords

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile = "InputDomainDecomposition.dat"
        integer                                     :: STAT_CALL, iflag
        integer                                     :: ObjEnterData = 0
        
        write(*,*)
        write(*,*)"Reading options"

        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR01'

        call GetData(BathymetryFile,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'BATHYMETRY_FILE',                  &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR02'

        call GetData(OutputFile,                                        &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'OUTPUT_FILE',                      &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR03'
        
        
        call GetData(nSectionsX,                                        &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'N_DOMAINS_X',                      &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR04'
        
        call GetData(nSectionsY,                                        &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'N_DOMAINS_Y',                      &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR05'
        
        call GetData(HasMaster,                                         &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'HAS_MASTER',                       &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR06'
        
        call GetData(GeometryFile,                                      &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     keyword      = 'GEOMETRY_FILE',                    &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR07'
        
        call GetData(DivisionMethod,                                    &
                     ObjEnterData, iflag   ,                            &
                     SearchType   = FromFile,                           &
                     Default      = Grid,                               &
                     keyword      = 'DIVISION_METHOD',                  &
                     ClientModule = 'MohidDomainDecomposition',         &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR08'
        
        if(DivisionMethod == HorizontalStripes)then
            write(*,*)
            write(*,*)"DIVISION_METHOD is Horizontal Stripes. Setting N_DOMAINS_X to 1."
            write(*,*)
            nSectionsX = 1
        elseif(DivisionMethod == VerticalStripes)then
            write(*,*)
            write(*,*)"DIVISION_METHOD is Vertical Stripes. Setting N_DOMAINS_Y to 1."
            write(*,*)
            nSectionsY = 1            
        endif

        call GetData(UseCenterMass,                                     &
                    ObjEnterData, iflag   ,                             &
                    SearchType   = FromFile,                            &
                    keyword      = 'USE_CENTER_OF_MASS',                &
                    ClientModule = 'MohidDomainDecomposition',          &
                    STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR09'
        
        if((DivisionMethod == HorizontalStripes .or. DivisionMethod == VerticalStripes) .and. UseCenterMass)then
            write(*,*)"USE_CENTER_OF_MASS is only possible when DIVISION_METHOD : 1"
            stop 'ReadKeywords - MohidDomainDecomposition - ERR09'
        endif
 
        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - MohidDomainDecomposition - ERR10'

    end subroutine ReadKeywords
    
    subroutine ConstructDomain
     
        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        type (T_Time  )                             :: CurrentTime
        
        !Begin-----------------------------------------------------------------
    
        write(*,*)"Constructing horizontal grid..."

        !Horizontal Grid
        call ConstructHorizontalGrid(HorizontalGridID = ObjHorizontalGrid,          &
                                     DataFile         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR01'

        write(*,*)"Constructing bathymetry..."
        
        
        !Horizontal Grid Data - Water Column (Bathymetry)
        call ConstructGridData      (GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     FileName         = BathymetryFile,             &
                                     STAT             = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR02'
        
        write(*,*)"Constructing mapping..."
        
        call SetDate(CurrentTime, 2016, 1, 1, 0, 0, 0)
        
        !Horizontal Map
        call ConstructHorizontalMap (HorizontalMapID  = ObjHorizontalMap,           &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     ActualTime       = CurrentTime,                &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR03'
        
        call GetGridData(ObjBathymetry, Bathymetry, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR03a'

        
        call GetWaterPoints2D(HorizontalMapID  = ObjHorizontalMap,                  &
                              WaterPoints2D    = WaterPoints2D,                     &
                              STAT             = STAT_CALL)  

        call GetHorizontalGridSize(ObjHorizontalGrid, WorkSize = WorkSize, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR04'
        
        write(*,*)"Constructing geometry"

        !Geometry - Water Column
        call ConstructGeometry      (GeometryID       = ObjGeometry,                &
                                     GridDataID       = ObjBathymetry,              &
                                     HorizontalGridID = ObjHorizontalGrid,          &
                                     HorizontalMapID  = ObjHorizontalMap,           &
                                     NewDomain        = GeometryFile,               &
                                    STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR05'

        write(*,*)"Constructing 3D map"

        !Map - Water Column            
        call ConstructMap           (Map_ID           = ObjMap,                     &
                                     GeometryID       = ObjGeometry,                &
                                     HorizontalMapID  = ObjHorizontalMap,           &
                                     STAT             = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR06'
        
        call GetWaterPoints3D(ObjMap, WaterPoints3D = WaterPoints3D, STAT = STAT_CALL)  
        if (STAT_CALL /= SUCCESS_) stop 'ConstructDomain - MohidDomainDecomposition - ERR07'

    end subroutine ConstructDomain

end program MohidDomainDecomposition

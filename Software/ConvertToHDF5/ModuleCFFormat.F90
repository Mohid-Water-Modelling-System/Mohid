!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : MM5Format
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : March 2005
! REVISION      : Pedro Montero - v4.0
! DESCRIPTION   : Module to convert CF Surface files into HDF5 format.
!                 Compiler settings: Need netCDF Lib (from MeteoGalicia)
!                 Optimized for CVF6.6. Written in ANSI FORTRAN 95
! HISTORY       : 20060628 Pedro Montero: A file with one z-level is considered
!                                         bidimensional
!
!------------------------------------------------------------------------------
!DataFile
!
!   FILENAME                    : char              -           !Path to CF original file
!   OUTPUT_GRID_FILENAME        : char              -           !Path to grid data file generated from CF file
!
!   <<BeginFields>>
!   air temperature
!   atmospheric pressure
!   wind velocity X
!   ...
!   ... (see below for available fields)
!   <<EndFields>>

!---Available properties in CF output file (copy from MOHID NAME column to input data file)
!
!netcdf 20020924 {
!dimensions:
! lon = 194 ;
! lat = 100 ;
! time = UNLIMITED ; // (4 currently)
!variables:
! float lon(lon) ;
!  lon:long_name = "longitude" ;
!  lon:standard_name = "longitude" ;
!  lon:units = "degrees_east" ;
!  lon:valid_min = 0.f ;
!  lon:valid_max = 0.f ;
!  lon:step = 0.2f ;
! float lat(lat) ;
!  lat:long_name = "latitude" ;
!  lat:standard_name = "latitude" ;
!  lat:units = "degrees_north" ;
!  lat:valid_min = 0.f ;
!  lat:valid_max = 0.f ;
!  lat:step = 0.2f ;
! float time(time) ;
!  time:long_name = "time" ;
!  time:standard_name = "time" ;
!  time:units = "second since 2002-09-24 00:00:00" ;
!  time:valid_min = 0.f ;
!  time:valid_max = 0.f ;
!  time:step = 0.f ;
! float u(time, lat, lon) ;
!  u:long_name = "10 m eastward wind" ;
!  u:standard_name = "eastward_wind" ;
!  u:units = "m s-1" ;
! float v(time, lat, lon) ;
!  v:long_name = "10 m northward wind" ;
!  v:standard_name = "northward_wind" ;
!  v:units = "m s-1" ;
! float rh(time, lat, lon) ;
!  rh:long_name = "2 m relative humidity" ;
!  rh:standard_name = "relative_humidity" ;
!  rh:units = "1" ;
! float t(time, lat, lon) ;
!  t:long_name = "2 m air temperature" ;
!  t:standard_name = "air_temperature" ;
!  t:units = "K" ;
! float e(time, lat, lon) ;
!  e:long_name = "6 hours acumulated evaporation" ;
!  e:standard_name = "water_evaporation_amount" ;
!  e:units = "kg m-2" ;
! float lhf(time, lat, lon) ;
!  lhf:long_name = "latent heat flux" ;
!  lhf:standard_name = "surface_downward_latent_heat_flux" ;
!  lhf:units = "W m-2" ;
!  lhf:comment = "INM Original data is 6 hours acumulated variable.This one is obtained from var/(" ;
! float shf(time, lat, lon) ;
!  shf:long_name = "sensible heat flux" ;
!  shf:standard_name = "surface_downward_sensible_heat_flux" ;
!  shf:units = "W m-2" ;
!  shf:comment = "INM Original data is 6 hours acumulated variable.This one is obtained from var/(" ;
! float mslp(time, lat, lon) ;
!  mslp:long_name = "mean sea level preasure" ;
!  mslp:standard_name = "air_pressure_at_sea_level" ;
!  mslp:units = "Pa" ;
! float sst(time, lat, lon) ;
!  sst:long_name = "sea surface temperature" ;
!  sst:standard_name = "sea_surface_temperature" ;
!  sst:units = "K" ;
! float lwr(time, lat, lon) ;
!  lwr:long_name = "long wave net radiation flux" ;
!  lwr:standard_name = "surface_net_downward_longwave_flux" ;
!  lwr:units = "W m-2" ;
!  lwr:comment = "INM Original data is 6 hours acumulated variable.This one is obtained from var/(" ;
! float swr(time, lat, lon) ;
!  swr:long_name = "short wave net radiation flux" ;
!  swr:standard_name = "surface_net_downward_shortwave_flux" ;
!  swr:units = "W m-2" ;
!  swr:comment = "INM Original data is 6 hours acumulated variable.This one is obtained from var/(" ;
! float tcc(time, lat, lon) ;
!  tcc:long_name = "Total Cloudness" ;
!  tcc:standard_name = "cloud_area_fraction" ;
!  tcc:units = "1" ;
! float tp(time, lat, lon) ;
!  tp:long_name = "6 hours acummulated total precipitation" ;
!  tp:standard_name = "precipitation_amount" ;
!  tp:units = "kg m-2" ;
!

!m s-1   ? baroclinic_eastward_sea_water_velocity
!m s-1   ? baroclinic_northward_sea_water_velocity
!m s-1   ? barotropic_eastward_sea_water_velocity
!m s-1   ? barotropic_northward_sea_water_velocity
!m s-1 50  ? northward_sea_water_velocity
!1   ? ocean_s_coordinate
!1   ? ocean_sigma_coordinate
!m s-1 49  ? eastward_sea_water_velocity
!psu long_name= "water salinity"   standard_name= "sea_water_salinity"
!"celsius deg"  long_name="water temperature" standard_name="sea_water_temperature"

!// global attributes:
!  :Title = "Analisis do INM para a superficie" ;
!  :Conventions = "CF-1.0" ;
!  :netcdf_version_id = "3.5.1" ;
!  :history = "INM2cdf 17-02-2005" ;
!  :date = 0 ;
!  :source = "HIRLAM" ;
!  :institution = "INM" ;
!  :references = "http://www.inm.es" ;
!}



Module ModuleCFFormat

    use ModuleGlobalData
    use ModuleHDF5
    use ModuleEnterData
    use ModuleTime
    use ModuleGridData
    use ModuleHorizontalGrid

    use ncdflib

    implicit none

    private

    !Subroutines---------------------------------------------------------------

    !Constructor
    public  :: ConvertCFFormat

    !Types---------------------------------------------------------------------
    type       T_Date
        type(T_Time)                                        :: Date
        type(T_Date), pointer                               :: Next
    end type  T_Date

    type       T_Field
        character(len=StringLength)                         :: Name
        character(len=StringLength)                         :: Units
        logical                                             :: Convert                  = .false.
        integer                                             :: nDimensions
        integer                                             :: GridLocation
        type(T_Time)                                        :: Date
        real, dimension(:,:  ),     pointer                 :: Values2D
        real, dimension(:,:,:),     pointer                 :: Values3D
        real                                                :: Fillvalue
        integer                                             :: OutputNumber             = 1
        type(T_Size3D)                                      :: Size, WorkSize
        type(T_Field),              pointer                 :: Next
    end type  T_Field

    private :: T_CFFormat
    type       T_CFFormat
        logical                                 :: OutputGridNotSpecified=.false.
        logical                                 :: Proyected=.false.
        logical                                 :: Dim3D=.false.
        logical                                 :: Polcom_data=.false.
        integer                                 :: ObjEnterData         = 0
        integer                                 :: ObjHDF5              = 0
        integer                                 :: ObjHorizontalGrid    = 0
        integer                                 :: Unit
        character(len=PathLength)               :: FileName
        character(len=PathLength)               :: BathymetryFile
        character(len=PathLength)               :: SZZFile
        character(len=PathLength)               :: GridFileName
        character(len=PathLength)               :: OutputFileName
        real, dimension(:,:),       pointer     :: Bathymetry
        real, dimension(:,:),       pointer     :: ConnectionX
        real, dimension(:,:),       pointer     :: ConnectionY
        real, dimension(:,:,:),     pointer     :: SZZ
        integer, dimension(:,:,:),  pointer     :: WaterPoints3D
        real, dimension(:  ),       pointer     :: XX
        real, dimension(:  ),       pointer     :: YY
        real                                    :: Xorig
        real                                    :: Yorig
        type(T_Size3D)                          :: Size, WorkSize
        type(T_Field),              pointer     :: FirstField
        type(T_Date),               pointer     :: FirstDate
    end type  T_CFFormat

    type(T_CFFormat), pointer                   :: Me

        real                                    :: lat1,lon1
        integer                                 :: idump

!--------------------------------------------------------------------------------

   contains

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine ConvertCFFormat(EnterDataID, STAT)

        !Arguments---------------------------------------------------------------
        integer,           intent(IN )                  :: EnterDataID
        integer, optional, intent(OUT)                  :: STAT

        !------------------------------------------------------------------------

        STAT = UNKNOWN_

        nullify (Me)
        allocate(Me)

        Me%ObjEnterData = AssociateInstance (mENTERDATA_, EnterDataID)

        call ReadOptions

        call OpenAndReadCFFile

        call ConstructGrid

        call OutputFields

        call KillCFFormat

        STAT = SUCCESS_

        deallocate(Me)
        nullify   (Me)

    end subroutine ConvertCFFormat


    !------------------------------------------------------------------------


    subroutine ReadOptions



        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: iflag



        !Begin-----------------------------------------------------------------

        call GetData(Me%FileName,                                       &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'FILENAME',                         &
                     ClientModule = 'ModuleCFFormat',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCFFormat - ERR01'

        if (iflag == 0)then
            write(*,*)'Must specify name of file to convert'
            stop 'ReadOptions - ModuleCFFormat - ERR02'
        end if


        call GetData(Me%GridFileName,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUT_GRID_FILENAME',             &
                     ClientModule = 'ModuleCFFormat',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCFFormat - ERR02'

        if (iflag == 0)then
            Me%OutputGridNotSpecified = .true.
            write(*,*)
            write(*,*)'Output grid file (Mohid format) was not specified.'
            write(*,*)'    - Mohid grid will not be constructed.'
            write(*,*)'    - HDF5 file will not be usable to interpolate'
            write(*,*)'      to other grids.'
            write(*,*)
        end if


        call GetData(Me%OutputFileName,                                 &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'OUTPUTFILENAME',                   &
                     ClientModule = 'ModuleCFFormat',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCFFormat - ERR03'



        call GetData(Me%Polcom_data,                              &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'POLCOM_DATA',            &
                     ClientModule = 'ConvertToHDF5',                    &
                     Default      = .false.,                            &
                     STAT         = STAT_CALL)        
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCFPolcomFormat - ERR04'
    
        call GetData(Me%BathymetryFile,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'BATHYMETRYFILE',             &
                     ClientModule = 'ModuleCFPolcomFormat',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCFPolcomFormat - ERR05'

        call GetData(Me%SZZFile,                                   &
                     Me%ObjEnterData, iflag,                            &
                     SearchType   = FromBlock,                          &
                     keyword      = 'SZZFILE',                         &
                     ClientModule = 'ModuleCFPolcomFormat',                  &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadOptions - ModuleCFPolcomFormat - ERR06'



    end subroutine ReadOptions
   !--------------------------------------------------------------------------


    subroutine OpenAndReadCFFile



    type(T_fileCF)                              :: FileCF

    !Others----------------------------------------------------------------
    type(T_Field), pointer                      :: NewField
    logical                                     :: exist
    type(T_Date),pointer                        :: Current

    integer                                     :: WILB,WIUB,WJLB,WJUB,WKLB,WKUB
    integer                                     :: ILB,IUB,JLB,JUB,KLB,KUB
    integer                                     :: i,j,k,n,nv
    integer                                     :: Year, Month, Day, Hour, Minute, Second
    integer                                     :: HourUTC, MinuteUTC
    real                                        :: Factor

    type(T_varCF)                               :: var
    character(len=StringLength)                 :: MohidName


    real                                         :: DIF



   !Begin-----------------------------------------------------------------

     write(*,*)'---------------------------'
     write(*,*)
     write(*,*)'Reading CF output file...'

     nullify(NewField)
     nullify(Me%FirstField)
     nullify(Me%FirstDate )
     nullify(Me%Bathymetry)
     nullify(Me%XX        )
     nullify(Me%YY        )

     !Verifies if file exists
     inquire(file = Me%FileName, exist = exist)
     if (.not. exist) then
        write(*,*)'CF File does not exist'
        stop 'OpenAndReadCFFile - ModuleCFFormat - ERR01'
     endif

   !Read head from CF file

    FileCF%name=Me%FileName
    call NCDF_LE_CAB(FileCF)
    if(FileCF%nDimensions==4 .and. FileCF%zAxis%size>1) Me%Dim3D=.true.

    if (FileCF%projection%pType/=0)  Me%Proyected=.true.

    !Grid Dimensions

    Me%WorkSize%ILB=1
    Me%WorkSize%IUB=FileCF%yAxis%Size
    Me%WorkSize%JLB=1
    Me%WorkSize%JUB=FileCF%xAxis%Size
    if (Me%Dim3D) then
      Me%WorkSize%KLB=1
      Me%WorkSize%KUB=FileCF%zAxis%Size
    else
      Me%WorkSize%KLB=1
      Me%WorkSize%KUB=1
    endif





    Me%Size%ILB     = Me%WorkSize%ILB - 1
    Me%Size%IUB     = Me%WorkSize%IUB + 1
    Me%Size%JLB     = Me%WorkSize%JLB - 1
    Me%Size%JUB     = Me%WorkSize%JUB + 1
    Me%Size%KLB     = Me%WorkSize%KLB - 1
    Me%Size%KUB     = Me%WorkSize%KUB + 1



            WILB            = Me%WorkSize%ILB
            WIUB            = Me%WorkSize%IUB
            WJLB            = Me%WorkSize%JLB
            WJUB            = Me%WorkSize%JUB
            WKLB            = Me%WorkSize%KLB
            WKUB            = Me%WorkSize%KUB

            ILB             = Me%Size%ILB
            IUB             = Me%Size%IUB
            JLB             = Me%Size%JLB
            JUB             = Me%Size%JUB
            KLB             = Me%Size%KLB
            KUB             = Me%Size%KUB

    if (.not. Me%Proyected) then


      !Calculates XX
       allocate(Me%XX(JLB:JUB))
       Me%Xorig=FileCF%xAxis%value(WJLB)-(FileCF%xAxis%value(WJLB+1)-FileCF%xAxis%value(WJLB))/2.
       Me%XX(JLB)=0.
       Me%XX(JUB)=(FileCF%xAxis%value(WJUB)+(FileCF%xAxis%value(WJUB)-FileCF%xAxis%value(WJUB-1))/2)-Me%Xorig
       Me%XX(WJLB)=(FileCF%xAxis%value(WJLB)+FileCF%xAxis%value(WJLB+1))/2.- Me%Xorig

       do j=WJLB+1,WJUB
          Me%XX(j)=(FileCF%xAxis%value(j)+FileCF%xAxis%value(j-1))/2.- Me%Xorig
       enddo



      !Calculates YY
       allocate(Me%YY(ILB:IUB))
       Me%Yorig=FileCF%yAxis%value(WILB)-(FileCF%yAxis%value(WILB+1)-FileCF%yAxis%value(WILB))/2.
       Me%YY(ILB)=0.
       Me%YY(IUB)=(FileCF%yAxis%value(WIUB)+(FileCF%yAxis%value(WIUB)-FileCF%yAxis%value(WIUB-1))/2)-Me%Yorig
       Me%YY(WILB)=(FileCF%yAxis%value(WILB)+FileCF%yAxis%value(WILB+1))/2.- Me%Yorig

       do i=WILB+1,WIUB
          Me%YY(i)=(FileCF%yAxis%value(i)+FileCF%yAxis%value(i-1))/2.- Me%Yorig
       enddo


    else

      nullify(Me%ConnectionX)
      nullify(Me%ConnectionY)
      allocate(Me%ConnectionX(ILB:IUB, JLB:JUB))
      allocate(Me%ConnectionY(ILB:IUB, JLB:JUB))

      !Procuramos a lonxitude
      !FileCF%Nvariables) !Sumase 1 e restase pois a proxeccion contabiliza como var
      do n=(FileCF%nDimensions+1),(FileCF%nDimensions+ FileCF%nVariables) 

        call NCDF_LE_VAR_CAB(FileCF, var,n)
        if (var%standardName=='longitude') then



           call NCDF_LE_VAR(FileCF, var)
           Me%ConnectionX(WILB:WIUB, WJLB:WJUB) = transpose(var%value2D(WJLB:WJUB, WILB:WIUB))

           do i=WILB,WIUB
              Me%ConnectionX(i, WJUB+1) = var%value2D(WJUB, i)+.5*(var%value2D(WJUB, i)-var%value2D(WJUB-1, i))
              Me%ConnectionX(i, WJLB-1) = var%value2D(WJLB, i)+.5*(var%value2D(WJLB+1, i)-var%value2D(WJLB, i))
           enddo

           do j=WJLB,WJUB
              Me%ConnectionX(WIUB+1,j ) = var%value2D(j,WIUB )+.5*(var%value2D( j,WIUB)-var%value2D(j,WIUB-1 ))
              Me%ConnectionX(WILB-1,j ) = var%value2D(j,WILB )+.5*(var%value2D( j,WILB+1)-var%value2D(j,WILB ))
           enddo

           Me%ConnectionX(WIUB+1,WJUB+1) = var%value2D(WJUB,   WIUB)  +.5*(var%value2D(WJUB, WIUB) - &
                                           var%value2D(WJUB-1, WIUB)) +.5*(var%value2D( WJUB,WIUB) - &
                                           var%value2D(WJUB,WIUB-1 ))
           Me%ConnectionX(WILB-1,WJLB-1) = Me%ConnectionX(WILB,WJLB)
           Me%ConnectionX(WIUB+1,WJLB-1) = Me%ConnectionX(WIUB,WJLB)
           Me%ConnectionX(WILB-1,WJUB+1) = Me%ConnectionX(WILB,WJUB)

       elseif (var%standardName=='latitude') then

           call NCDF_LE_VAR(FileCF, var)
           Me%ConnectionY(WILB:WIUB, WJLB:WJUB) = transpose(var%value2D(WJLB:WJUB, WILB:WIUB))



          do i=WILB,WIUB
             Me%ConnectionY(i, WJUB+1) = var%value2D(WJUB, i)+.5*(var%value2D(WJUB, i)-var%value2D(WJUB-1, i))
             Me%ConnectionY(i, WJLB-1) = var%value2D(WJLB, i)+.5*(var%value2D(WJLB+1, i)-var%value2D(WJLB, i))
          enddo

          do j=WJLB,WJUB
             Me%ConnectionY(WIUB+1,j ) = var%value2D(j,WIUB )+.5*(var%value2D( j,WIUB)-var%value2D(j,WIUB-1 ))
             Me%ConnectionY(WILB-1,j ) = var%value2D(j,WILB )+.5*(var%value2D( j,WILB+1)-var%value2D(j,WILB ))
          enddo

          Me%ConnectionY(WIUB+1,WJUB+1) = var%value2D(WJUB,WIUB    )+.5*(var%value2D( WJUB,WIUB)- &
                                          var%value2D(WJUB,WIUB-1 ))+.5*(var%value2D(WJUB, WIUB)- &
                                          var%value2D(WJUB-1, WIUB))
          Me%ConnectionY(WILB-1,WJLB-1) = Me%ConnectionY(WILB,WJLB)
          Me%ConnectionY(WIUB+1,WJLB-1) = Me%ConnectionY(WIUB,WJLB)
          Me%ConnectionY(WILB-1,WJUB+1) = Me%ConnectionY(WILB,WJUB)
      endif


      enddo

    endif


    !Create a zero filled bathymetry


    allocate(Me%Bathymetry(ILB:IUB, JLB:JUB))

   if (.not.Me%Polcom_Data) then
    do j=WJLB,WJUB
     do i=WILB,WIUB
        Me%Bathymetry(i,j) = 0.
     enddo
    enddo

    !Create SZZ
    if (Me%Dim3D) then
       allocate(Me%SZZ(ILB:IUB, JLB:JUB, KLB:KUB))
       do j = Me%WorkSize%JLB, Me%WorkSize%JUB
         do i = Me%WorkSize%ILB, Me%WorkSize%IUB

            Me%SZZ(i,j,Me%WorkSize%KLB:Me%WorkSize%KUB) = FileCF%zAxis%value(Me%WorkSize%KLB:Me%WorkSize%KUB)

         enddo
       enddo
    endif

   endif

! Inicio Polcom Bathymetry Block
   if (Me%Polcom_Data) then
    ! Read Polcom Bathymetry
     
   allocate(Me%WaterPoints3D(Me%Size%ILB:Me%Size%IUB,&
                 Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))  
      Open(30,File=Me%BathymetryFile)
      do i=WILB,WIUB
       do j=WJLB,WJUB

        Read(30,*)  Me%Bathymetry(i,j)

        If (Me%Bathymetry(i,j)==   0.0000000E+00) then
            Me%Bathymetry(i,j)=-999.9
            do k=WKLB,WKUB
                Me%waterpoints3D(i,j,k) = 0
            enddo
        else
            do k=WKLB,WKUB
                Me%waterpoints3D(i,j,k) = 1
            enddo
        endif
      
        enddo
        enddo
        close(30)

        write(*,*) 'ok, batimetria leida bien'

    !Read SZZ array

    if(Me%Dim3D) then
        nullify(Me%SZZ)
        allocate(Me%SZZ(ILB:IUB, JLB:JUB, KLB-1:KUB))

      Open(20,File=Me%SZZfile)

       do i = WILB, WIUB
               read(20,*)
         do j = WJLB, WJUB

           Read(20,110) idump,(Me%SZZ(i,j,k),k=WKLB,WKUB) 

               If (Me%SZZ(i,j,Me%WorkSize%KLB)==0) then

                    Me%SZZ(i,j,Me%WorkSize%KLB)=-999.9
                    Me%SZZ(i,j,Me%WorkSize%KLB-1)=-999.9

               else 
                    Me%SZZ(i,j,Me%WorkSize%KLB-1)=-1
                    Me%SZZ(i,j,Me%WorkSize%KUB)=0.0

                 endif

                  Do k=Me%WorkSize%KLB,Me%WorkSize%KUB-1

                   If (Me%SZZ(i,j,K)==0) then
                     Me%SZZ(i,j,K)=-999.9
                   else
                     DIF=Me%SZZ(i,j,K)-Me%SZZ(i,j,K-1)
                     Me%SZZ(i,j,K)=Me%SZZ(i,j,K)+DIF

                    endif  
                  enddo

 110  FORMAT(i4,32f10.5)

       enddo
       enddo
       close (20)
      endif
     endif

! Fin Polcom Bathymetry Block


    !Loop of dates

    do n=1,FileCF%tAxis%Size

       !Create and add date


        call AddDate(Me%FirstDate,Current)

        call Read_CFDate(FileCF%tAxis%unit,                       &
        Year, Month, Day, Hour, Minute, Second, HourUTC, MinuteUTC, Factor)
   
        
        call SetDate(Current%Date, Year, Month, Day, Hour, Minute, Second)

        Current%Date=Current%Date+FileCF%tAxis%value(n)*Factor
        Current%Date=Current%Date+3600*HourUTC+sign(60*MinuteUTC,HourUTC)
       


       !find out variables and fill fields

       do nv=FileCF%NDimensions+1, (FileCF%Nvariables+FileCF%NDimensions)

          call NCDF_LE_VAR_CAB (FileCF, var, nv)
          if(CheckName(var, MohidName))then

            call AddField(Me%FirstField, NewField)

                    NewField%WorkSize%ILB   = WILB
                    NewField%WorkSize%IUB   = WIUB
                    NewField%WorkSize%JLB   = WJLB
                    NewField%WorkSize%JUB   = WJUB


                    NewField%Size%ILB       = ILB
                    NewField%Size%IUB       = IUB
                    NewField%Size%JLB       = JLB
                    NewField%Size%JUB       = JUB

                    if(var%nDimensions==4) then
                       NewField%WorkSize%KLB   = WKLB
                       NewField%WorkSize%KUB   = WKUB
                       NewField%Size%KLB       = KLB
                       NewField%Size%KUB       = KUB

                    endif

                    NewField%Name        =   trim(MohidName)
                    NewField%Date   = Current%Date
                    NewField%Units  = var%unit
                    NewField%nDimensions=var%nDimensions-1

            if(var%nDimensions==3) then
              allocate(NewField%Values2D(ILB:IUB, JLB:JUB))
              call NCDF_LE_VAR(FileCF,var,nTime=n)
              NewField%Values2D(WILB:WIUB,WJLB:WJUB) = transpose(var%value2D(WJLB:WJUB,WILB:WIUB))

  ! Bloque para meter una mascara correcta
        
            if (.not.Me%Polcom_Data) then

                if (nv==FileCF%Nvariables+FileCF%NDimensions) then
!     write(*,*) var%name, var%nDimensions, var%value2D(WJUB,WIUB)
                    allocate(Me%WaterPoints3D(Me%Size%ILB:Me%Size%IUB,&
                             Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
                    do j=WJLB,WJUB
                    do i=WILB,WIUB
                        if (var%value2D(j,i)== -999) then
                        !                write(*,*) 'si hay valores -999'
                            Me%waterpoints3D(i,j,1) = 0
                        else
                            Me%waterpoints3D(i,j,1) = 1
                        endif
                    enddo
                    enddo
            endif
            endif

! Fin Bloque para meter una mascara Tierra-mar correcta en meteorologicos

            elseif(var%nDimensions==4) then
              allocate(NewField%Values3D(ILB:IUB, JLB:JUB, KLB:KUB))

              call NCDF_LE_VAR(FileCF,var,nTime=n)
              do i=WILB,WIUB
               do j=WJLB,WJUB
                do k=WKLB,WKUB
                   NewField%Values3D(i,j,k) = var%value3D(j,i,k)
                enddo
               enddo
              enddo

            endif

            call transform2Mohid(var,NewField)
        endif

         if(var%nDimensions==4) then
   
              do i=WILB,WIUB
               do j=WJLB,WJUB
                do k=WKLB,WKUB
      if (NewField%Values3D(i,j,k)== NewField%Fillvalue)then
                   NewField%Values3D(i,j,k)=0. 
                   Me%bathymetry(i,j)=-999.9
       Me%waterpoints3D(i,j,k) = 0

                  endif
                enddo
               enddo
              enddo

            endif

        enddo !Variabeis



    enddo !Times


    end subroutine OpenAndReadCFFile



   !--------------------------------------------------------------------------

    logical function CheckName(CFVar, MohidName)

        !Arguments-----------------------------------------------------------
        Type(T_varCF)                   :: CFVar
        character(Len=StringLength)     :: MohidName

        !Begin-----------------------------------------------------------------


        if (trim(CFVar%standardName)/='') then

          select case(trim(CFVar%standardName))

            case('eastward_wind')

                MohidName = GetPropertyName(WindVelocityX_)
                CheckName = .true.

            case('northward_wind')

                MohidName = GetPropertyName(WindVelocityY_)
                CheckName = .true.

            case('air_pressure_at_sea_level')

                MohidName = GetPropertyName(AtmosphericPressure_)
                CheckName = .true.

            case('air_temperature')

                MohidName = GetPropertyName(AirTemperature_)
                CheckName = .true.
            case('relative_humidity')

                MohidName = GetPropertyName(RelativeHumidity_)
                CheckName = .true.
            case('water_evaporation_amount')

                MohidName = GetPropertyName(Evaporation_)
                CheckName = .true.

            case('cloud_area_fraction')

                MohidName = GetPropertyName(CloudCover_)
                CheckName = .true.
            case('precipitation_amount')

                MohidName = GetPropertyName(Precipitation_)
                CheckName = .true.





!            case('surface_temperature')

!                MohidName = GetPropertyName(AirTemperature_)
!                CheckName = .true.


            case('surface_downwelling_longwave_flux_in_air')

               MohidName = GetPropertyName(DownwardLongWaveRadiation_)
               CheckName = .true.
            case('surface_downwelling_shortwave_flux_in_air')

               MohidName = GetPropertyName(SolarRadiation_)
               CheckName = .true.


            case('water_vapor_mixing_ratio')

               MohidName = "2-meter mixing ratio"
               CheckName = .true.

            case('sea_surface_temperature')

               MohidName = "sea water temperature"
               CheckName = .true.

            case('surface_downward_sensible_heat_flux')

               MohidName = GetPropertyName(SensibleHeat_)
               CheckName = .true.

            case('surface_downward_latent_heat_flux')

               MohidName = GetPropertyName(LatentHeat_)
               CheckName = .true.

            case('land_binary_mask')

               MohidName = "land mask"
               CheckName = .true.
! Para modelos oceanográficos. Añadido por Silvia  23/02/06

!            case('sea_water_potential_temperature')
!
!               MohidName = GetPropertyName(Temperature_)
!               CheckName = .true.

            case('sea_water_temperature')

               MohidName = GetPropertyName(Temperature_)
               CheckName = .true.

            case('sea_water_salinity')

               MohidName = GetPropertyName(Salinity_)
               CheckName = .true.

            case('sea_surface_elevation')

               MohidName = GetPropertyName(WaterLevel_)
               CheckName = .true.


           
            case('eastward_sea_water_velocity')

                MohidName = GetPropertyName(VelocityU_)
                CheckName = .true.

            case('northward_sea_water_velocity')

                MohidName = GetPropertyName(VelocityV_)
                CheckName = .true.

            case('barotropic_eastward_sea_water_velocity')

                MohidName = GetPropertyName(BarotropicVelocityU_)
                CheckName = .true.

            case('barotropic_northward_sea_water_velocity')

                MohidName = GetPropertyName(BarotropicVelocityV_)
                CheckName = .true.

            case('sea_surface_height_above_sea_level')

                MohidName = GetPropertyName(WaterLevel_)
                CheckName = .true.


            case default

                CheckName = .false.

          end select

        else
        select case(trim(CFVar%completeName))

            case('6 hours acumulated latent heat flux')

                MohidName = GetPropertyName(LatentHeat_)
                CheckName = .true.
            case('6 hours acumulated latent sensible flux')

                MohidName = GetPropertyName(SensibleHeat_)
                CheckName = .true.
            case default

                CheckName = .false.

          end select


        endif


    end function CheckName

   !------------------------------------------------------------------------
    subroutine Read_CFDate(DateIn,Year, Month, Day, Hour, Minute, Second, HourUTC, MinuteUTC, Factor)

        !Arguments-------------------------------------------------------------
        character(len=*),intent(IN )                    :: DateIn
        integer, intent(OUT )                           :: Year, Month, Day
        integer, intent(OUT )                           :: Hour, Minute, Second
        integer, intent(OUT )                           :: HourUTC, MinuteUTC
        real,    intent(OUT)                            :: Factor

        !Types-----------------------------------------------------------------
        type  :: T_word
            integer              :: beginw
            integer              :: endw
            character(len=80)    :: string
        end type
        !Local-----------------------------------------------------------------
        integer                  :: n
        integer                  :: numWords

        type(T_word), pointer     :: Word(:)
        integer                   :: nword
        integer                   :: nletter

        integer                   :: nseparator
        real                      :: RealSecond




        !Begin-----------------------------------------------------------------

        numWords=1
        do n=1 , len_trim(DateIn)
          if (DateIn(n:n)==" ") numWords=numWords+1
        enddo


        nullify(word)
        allocate(Word(numWords))

        do n=1,numwords
            word(n)%string=""
        enddo


        nletter=1
        nword=1



        do n =1, len_trim(DateIn)
            if (DateIn(n:n)==" ") then
                nword=nword+1
                nletter=1
            else
                word(nword)%string(nletter:nletter)=DateIn(n:n)
                nletter=nletter+1
            endif 
        enddo


        do n=1,3
            read(word(3)%string(1:4),'(i4)')  Year
            read(word(3)%string(6:7),'(i2)')  Month
            read(word(3)%string(9:10),'(i2)') Day
        enddo

        if (numwords>3) then
   
            read(word(4)%string(1:2),'(i2)')  hour
            read(word(4)%string(4:5),'(i2)')  minute
            read(word(4)%string(7:len_trim(word(4)%string)),'(f10.5)') RealSecond
            second=int(RealSecond)
        else 

            hour=0
            minute=0
            second=0

        endif

        if (numwords>4) then

            do n=1,len_trim(word(5)%string) 
                if (word(5)%string(n:n)==":") then
                    nseparator=n
                    read(word(5)%string(1:nseparator-1),'(i3)')  hourUTC
                    read(word(5)%string(nseparator+1:len_trim(word(5)%string)),'(i2)')  minuteUTC
         
                    exit

                else 
                    read(word(5)%string(1:len_trim(word(5)%string)),'(i3)')  hourUTC
                    minuteUTC=0
                endif

            enddo
                    
           
        else 

            hourUTC=0
            minuteUTC=0
 
        endif

        if (word(1)%string(1:1)=="d") then
            factor=3600*24.
        elseif (word(1)%string(1:1)=="h") then
             factor=3600.
        elseif (word(1)%string(1:1)=="m") then
            factor=60.
        elseif (word(1)%string(1:1)=="s") then
            factor=1.
        endif

    end subroutine Read_CFDate

   !
   !------------------------------------------------------------------------

    subroutine transform2Mohid(var,hdfField)

      !Arguments-----------------------------------------------------------
        Type(T_varCF)                   :: var
        Type(T_Field)                   :: hdfField

        !Begin-----------------------------------------------------------------

       if (Me%Polcom_Data) then

        
        hdfField%Values3D=hdfField%Values3D*var%scale_factor+var%add_offset
        hdfField%Values2D=hdfField%Values2D*var%scale_factor+var%add_offset
        hdfField%Fillvalue=var%Fillvalue*var%scale_factor+var%add_offset
       
    end if

        if (var%Unit=="K") then

            hdfField%Values3D=hdfField%Values3D- AbsoluteZero
            hdfField%Values2D=hdfField%Values2D- AbsoluteZero
            hdfField%Units="ºC"
        end if




    end subroutine transform2Mohid

   !
   !------------------------------------------------------------------------


    subroutine AddDate (FirstDate, ObjDate)

        !Arguments-------------------------------------------------------------
        type (T_Date), pointer                   :: FirstDate
        type (T_Date), pointer                   :: ObjDate

        !Local-----------------------------------------------------------------
        type (T_Date), pointer                   :: NewDate
        type (T_Date), pointer                   :: PreviousDate

        !Begin-----------------------------------------------------------------

        !Allocates new instance
        allocate (NewDate)
        nullify  (NewDate%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstDate)) then
            FirstDate         => NewDate
            ObjDate           => NewDate
        else
            PreviousDate      => FirstDate
            ObjDate           => FirstDate%Next
            do while (associated(ObjDate))
                PreviousDate  => ObjDate
                ObjDate       => ObjDate%Next
            enddo
            ObjDate           => NewDate
            PreviousDate%Next => NewDate
        endif


    end subroutine AddDate






     !------------------------------------------------------------------------


    subroutine ConstructGrid

        !Local-----------------------------------------------------------------
        real                        :: Latitude, Longitude
        integer                     :: STAT_CALL
        type(T_Size2D)              :: WorkSize2D

        !Begin-----------------------------------------------------------------

        WorkSize2D%ILB = Me%WorkSize%ILB
        WorkSize2D%IUB = Me%WorkSize%IUB
        WorkSize2D%JLB = Me%WorkSize%JLB
        WorkSize2D%JUB = Me%WorkSize%JUB

    if (.not. Me%Proyected) then

        write(*,*)
        write(*,*)'Constructing grid...'


        call WriteGridData (FileName        = Me%GridFileName,              &
                            XX              = Me%XX,                        &
                            YY              = Me%YY,                        &
                            COMENT1         = " CF Grid based on file :",   &
                            COMENT2         = Me%FileName,                  &
                            WorkSize        = WorkSize2D ,                  &
                            CoordType       = 4,                            &
                            Xorig           = Me%Xorig,                     &
                            Yorig           = Me%Yorig,                     &
                            Zone            = 0,                            &
                            Grid_Angle      = 0.,                           &
                            Latitude        = Latitude,                     &
                            Longitude       = Longitude,                    &
                            GridData2D_Real = Me%Bathymetry,                &
                            Overwrite       = ON,                           &
                            FillValue       = -99.,                         &
                            STAT            = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleCFFormat - ERR01'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleCFFormat - ERR02'
    else

        write(*,*)
        write(*,*)'Constructing grid proyected...'



        call WriteGridData  (FileName       = Me%GridFileName,                              &
                             ConnectionX    = Me%ConnectionX,                               &
                             ConnectionY    = Me%ConnectionY,                               &
                             COMENT1        = 'Grid Data created from CF proyected file',   &
                             COMENT2        = trim(Me%FileName),                            &
                             WorkSize       = WorkSize2D,                                   &
                             CoordType      = 4,                                            &
                             Xorig          = 0.,                                           &
                             Yorig          = 0.,                                           &
                             Zone           = 29,                                           &
                             GRID_ANGLE     = 0.,                                           &
                             Latitude       = Me%ConnectionX(1,1),                          &
                             Longitude      = Me%ConnectionY(1,1),                          &
                             FillValue      = -99.,                                         &
                             Overwrite      = .true.,                                       &
                             GridData2D_Real= Me%Bathymetry,                                &
                             STAT           = STAT_CALL)
        if(STAT_CALL .ne. SUCCESS_)stop 'ConstructGrid - ModuleMM5Format - ERR03'

        call ConstructHorizontalGrid(Me%ObjHorizontalGrid, Me%GridFileName, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructGrid - ModuleCFFormat - ERR04'

    endif



    end subroutine ConstructGrid

        !------------------------------------------------------------------------

    subroutine AddField (FirstField, ObjField)

        !Arguments-------------------------------------------------------------
        type (T_Field), pointer                   :: FirstField
        type (T_Field), pointer                   :: ObjField

        !Local-----------------------------------------------------------------
        type (T_Field), pointer                   :: NewField
        type (T_Field), pointer                   :: PreviousField

        !Begin-----------------------------------------------------------------

        !Allocates new instance
        allocate (NewField)
        nullify  (NewField%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstField)) then
            FirstField         => NewField
            ObjField           => NewField
        else
            PreviousField      => FirstField
            ObjField           => FirstField%Next
            do while (associated(ObjField))
                PreviousField  => ObjField
                ObjField       => ObjField%Next
            enddo
            ObjField           => NewField
            PreviousField%Next => NewField
        endif


    end subroutine AddField

        !------------------------------------------------------------------------

    subroutine OutputFields

        !Local-----------------------------------------------------------------
        real,    dimension(6), target                   :: AuxTime
        real,    dimension(:), pointer                  :: TimePtr
        integer                                         :: STAT_CALL, OutputNumber
        type(T_Field), pointer                          :: Field
        type(T_Date), pointer                           :: CurrentDate

        integer                                         :: i,j,k

        !Begin-----------------------------------------------------------------

        write(*,*)
        write(*,*)'Writing HDF5 file...'



        call Open_HDF5_OutPut_File

        OutputNumber = 1
        CurrentDate  => Me%FirstDate

        do while(associated(CurrentDate))

            call ExtractDate   (CurrentDate%Date,                           &
                                AuxTime(1), AuxTime(2), AuxTime(3),         &
                                AuxTime(4), AuxTime(5), AuxTime(6))

            TimePtr => AuxTime

            call HDF5SetLimits  (Me%ObjHDF5, 1, 6, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR01'


            call HDF5WriteData  (Me%ObjHDF5, "/Time",                       &
                                 "Time", "YYYY/MM/DD HH:MM:SS",             &
                                 Array1D = TimePtr,                         &
                                 OutputNumber = OutPutNumber, STAT = STAT_CALL)
            if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR02'


            Field => Me%FirstField

            do while(associated(Field))


                if(Field%Date == CurrentDate%Date)then

                    Field%OutputNumber = OutputNumber

if0:                if(Field%nDimensions == 2) then

                      call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                       Field%WorkSize%JLB, Field%WorkSize%JUB, STAT = STAT_CALL)
                      if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR05'

                      call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//Field%Name,                         &
                                       Field%Name,                                      &
                                       Field%Units,                                     &
                                       Array2D      = Field%Values2D,                   &
                                       OutputNumber = Field%OutputNumber,               &
                                       STAT         = STAT_CALL)
                      if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR06'

                    elseif(Field%nDimensions == 3) then

                      call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                       Field%WorkSize%JLB, Field%WorkSize%JUB, &
                                       Field%WorkSize%KLB, Field%WorkSize%KUB, STAT = STAT_CALL)
                      if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR05'

                      call HDF5WriteData(Me%ObjHDF5,                                      &
                                       "/Results/"//Field%Name,                         &
                                       Field%Name,                                      &
                                       Field%Units,                                     &
                                       Array3D      = Field%Values3D,                   &
                                       OutputNumber = Field%OutputNumber,               &
                                       STAT         = STAT_CALL)
                      if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR06'

!                  Write Polcom SZZ data into hdf5 file

                      if (Me%Polcom_Data) then

                      call HDF5SetLimits(Me%ObjHDF5, Field%WorkSize%ILB, Field%WorkSize%IUB,&
                                       Field%WorkSize%JLB, Field%WorkSize%JUB, &
                                       Field%WorkSize%KLB, Field%WorkSize%KUB, STAT = STAT_CALL)
                      if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCFPolcomFormat - ERR08'                              
                      
                 Do j=Me%workSize%jLB, Me%workSize%jUB        
                  Do i=Me%WorkSize%ILB, Me%workSize%IUB 
                   Do k=Me%WorkSize%KLB-1,Me%WorkSize%KUB
   
                     If (Me%Bathymetry(i,j)==-999.9) then 
                      Field%Values3D(i,j,k)= -999.9
                     else
                          Field%Values3D(i,j,k)=-Me%SZZ(i,j,k)* Me%Bathymetry(i,j)
                     endif  
                    enddo     
                   enddo
                 enddo
                call HDF5WriteData(Me%ObjHDF5,                                      &
                                      "/Grid/VerticalZ",                                &
                                       "Vertical",                                      &
                                       "m",                                             &
                                       Array3D      = Field%Values3D,                   &
                                       OutputNumber = Field%OutputNumber,               &
                                       STAT         = STAT_CALL)
                      if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleCFPolcomFormat - ERR09'                 
            
                   endif
!                 End write Polcom SZZ data into hdf5 file

     
     endif if0

                end if

                Field => Field%Next

            end do

            OutputNumber = OutputNumber + 1

            CurrentDate => CurrentDate%Next

        end do

        write(*,*)
        write(*,*)'Closing HDF5 file...'

        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'OutputFields - ModuleMM5Format - ERR07'



    end subroutine OutputFields


    !----------------------------------------------------------------------

    subroutine Open_HDF5_OutPut_File

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_CALL
        integer                                     :: HDF5_CREATE
        integer,    dimension(:,:,:), pointer       :: WaterPoints3D

        !----------------------------------------------------------------------



        allocate(WaterPoints3D(Me%Size%ILB:Me%Size%IUB,&
                               Me%Size%JLB:Me%Size%JUB, Me%Size%KLB:Me%Size%KUB))
print *,'flag0.1'
        !WaterPoints3D(:,:,:) = 1
        WaterPoints3D(:,:,:) = Me%waterpoints3D(:,:,:)
        !Gets File Access Code
print *,'flag0.2'
        call GetHDF5FileAccess  (HDF5_CREATE = HDF5_CREATE)
print *,'flag0.3'
        !Opens HDF5 File
        call ConstructHDF5(Me%ObjHDF5, Me%OutputFileName, HDF5_CREATE, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR01'
print *,'flag0.4'



          call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,&
                             Me%WorkSize%JLB, Me%WorkSize%JUB, STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR02'
print *,'flag0.5'

        call HDF5WriteData   (Me%ObjHDF5, "/Grid", "Bathymetry", "-",          &
                              Array2D =  Me%Bathymetry, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR03'


        call WriteHorizontalGrid (Me%ObjHorizontalGrid, Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR04'



          call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB, Me%WorkSize%IUB,    &
                                           Me%WorkSize%JLB, Me%WorkSize%JUB,    &
                                           Me%WorkSize%KLB, Me%WorkSize%KUB,    &
                                           STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR07'


          call HDF5WriteData   (Me%ObjHDF5, "/Grid", "WaterPoints", "-",    &
                              Array3D = WaterPoints3D,  STAT = STAT_CALL)
          if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR08'


        !Writes SZZ - Constant
  
  if (.not.Me%Polcom_Data) then

        if (Me%Dim3D) then
           call HDF5SetLimits  (Me%ObjHDF5, Me%WorkSize%ILB,   Me%WorkSize%IUB,  &
                                            Me%WorkSize%JLB,   Me%WorkSize%JUB,  &
                                            Me%WorkSize%KLB-1, Me%WorkSize%KUB,  &
                                            STAT = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR30'

           call HDF5WriteData  (Me%ObjHDF5, "/Grid/VerticalZ", "Vertical",             &
                            "m", Array3D = Me%SZZ,                                  &
                            STAT = STAT_CALL)
           if (STAT_CALL /= SUCCESS_) stop 'OutPut_Results_HDF - ModuleLevitusFormat - ERR40'
        endif
        endif


        !Writes everything to disk
        call HDF5FlushMemory (Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'Open_HDF5_OutPut_File - ModuleMM5Format - ERR09'

        deallocate(WaterPoints3D)
        nullify   (WaterPoints3D)

    end subroutine Open_HDF5_OutPut_File


    !--------------------------------------------------------------------------

       !--------------------------------------------------------------------------

    subroutine KillCFFormat

        !Local-----------------------------------------------------------------
        integer                     :: STAT_CALL, nUsers

        !Begin-----------------------------------------------------------------

        call KillHorizontalGrid(Me%ObjHorizontalGrid, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCFFormat - ModuleMM5Format - ERR01'


        call KillHDF5(Me%ObjHDF5, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_)stop 'KillCFFormat - ModuleMM5Format - ERR02'

        nUsers = DeassociateInstance(mENTERDATA_, Me%ObjEnterData)
        if (nUsers == 0) stop 'KillCFFormat - ModuleMM5Format - ERR03'


        deallocate(Me%Bathymetry)
        if(.not. Me%proyected) then
           deallocate(Me%XX)
           deallocate(Me%YY)
        else

           deallocate(Me%ConnectionX)
           deallocate(Me%ConnectionY)
        endif

        deallocate(Me%FirstField)
        deallocate(Me)
        nullify   (Me)


    end subroutine KillCFFormat

    !--------------------------------------------------------------------------

end module ModuleCFFormat



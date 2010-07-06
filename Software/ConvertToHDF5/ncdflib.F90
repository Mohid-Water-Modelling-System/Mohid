!****************************************************************************!
!                                                                            !
!  MODULE: ncdflib.f90                                                       !
!                                                                            !
!  PURPOSE:  Modulo que xestiona a saida e entrada do NetCDF.                !
!            A funcion deste modulo é ser unha capa entre a interface        !
!            en fortran e o codigo, de tal xeito que enchendo as estructuras !
!            de dados, sexa facil escrever e ler no convenio CF.             !
!                                                                            !
!            Precisa da proxeccion                                           !
!                                                                            !
!            Autor: Pedro Montero Vilar                                      !
!                                                                            !
!            Data: 16/12/2004                                                !
!            Version : 2.2 (febreiro 2007)                                   !
!            History:                                                        !
!                                                                            !
!            08062006 Corrige: en vez de escoger el primer eje como eje x lo !
!                               selecciona por el nombre                     !
!                                                                            !
!            14062006 Se añade _FillValue, addoffset y scale_factor          !
!****************************************************************************!
!                                                                            !
!   ATENTION: Precisa de compilar coa libreria netcdf                        !
!                                                                            !
!****************************************************************************!
!                                                                            !
!  Subrutinas:                                                               !
!    NCDF_CAB                                                                 !
!                                                                            !
!                                                                            !
!****************************************************************************!

Module ncdflib

#ifdef _USE_NIX
    use netcdf
#else
    Use netcdf90
#endif    
    implicit none
!
!----------------------------------------------------------------------------!
! Subroutines                                                                !
!----------------------------------------------------------------------------!
    private
    
    
    public      ::  NCDF_ESC_CAB                 !Escrebe a cabeceira do ficheiro
    private     ::      NCDF_ESC_EJE_ATT         !Escrebe os atributos do eixo (pertence a cabeceira)
    private     ::      NCDF_ESC_PROJ            !Escrebe unha proxeccion

    public      ::  NCDF_ESC_EJE                 !Escrebe os valores dun eixo (espacial)
    public      ::  NCDF_ESC_TIEMPO              !Escrebe o eixo temporal, un tempo ou todos
    public      ::  NCDF_ESC_VAR_CAB             !Escrebe a cabeceira dunha variable
    private     ::      NCDF_ESC_VAR_ATT         !Escrebe os atributos dunha variable (pertence a cabeceira)
    public      ::  NCDF_ESC_VAR                 !Escrebe os valores dunha variable

    public      ::  NCDF_LE_CAB                  !Le unha cabeceira dun ficheiro
    private     ::      NCDF_LE_EJE              !Le os eixos (a través da leitura da cabeceira)
    public      ::  NCDF_LE_VAR_CAB              !Le a cabeceira duna variable
    public      ::  NCDF_LE_VAR                  !Le unha os valores variable

    private     ::      erro                     !Mensaxe de erro
!----------------------------------------------------------------------------!
! Variabeis                                                                  !
!----------------------------------------------------------------------------!

    


!----------------------------------------------------------------------------!
    integer, parameter          :: lenStr=800             !Tamaño dos strings
    real, parameter             :: valMinD=99999.        !Valor moi grande
    real, parameter             :: valMaxD=-99999.       !Valor moi pequeno
    real, parameter             :: FillValueRef=-123456789.
!----------------------------------------------------------------------------!


    !Proxeccion
    public T_projectionCF 
    Type T_projectionCF
       integer                       :: pType=0             !Tipo: 0 Ningunha
       character(len=lenStr)         :: name                !Nome
       character(len=lenStr)         :: gmName              !Grid mapping nome
       real                          :: trueLat(2), trueLon !latitudes verdadeiras, lonxitudes verdadeiras
       real                          :: ctrLat,ctrLon       !latitude e lonxitude central
       real                          :: sclFactor           !Factor de escala
    endType
!----------------------------------------------------------------------------!

    
    !Eixo
    public T_axisCF
    Type T_axisCF
        character(len=lenStr)         :: name               !Nome do eixo 
        character(len=lenStr)         :: completeName=""    !Nome completo
        character(len=lenStr)         :: standardName=""    !Nome estandar
        character(len=lenStr)         :: unit=""            !unidade
        real                          :: maxVal,minVal      !Valor máximo e mínimo
        real                          :: step               !Paso
        character(len=lenStr)         :: aType=""           !Tipo de eixo
        integer                       :: size               !Tamaño (número de nodos)
        real, pointer                 :: value(:)           !valores
    endType
!----------------------------------------------------------------------------!


    !Ficheiro
    public T_fileCF
    type T_fileCF
        character(len=lenStr)         :: name               !Nome do ficheiro
        type(T_projectionCF)          :: projection         !Proxeccion   
      
        !Atributos globais
        character(len=lenStr)           :: title=""         !Titulo
        character(len=lenStr)           :: convention=""    !Convencion
        character(len=lenStr)           :: version=""       !Version de netcdf
        character(len=lenStr)           :: history=""       !historia
        integer                         :: idate            !
        character(len=lenStr)           :: origin=""        !Orixen
        character(len=lenStr)           :: institution=""   !Institudion
        character(len=lenStr)           :: references=""    !Referencias

        !Dimensions e eixos
        character(len=lenStr),dimension(4) :: nameDim         !Nome da dimension
        type(T_axisCF)                     :: xAxis,yAxis,zAxis,tAxis
                                                              !Eixos
        integer                            :: nVariables      !numero de variables
        integer                            :: nDimensions=-1  !numero de dimensions (espaciais mais a temporal)
        logical                            :: timeDim=.false. !Flag de existencia de dimension temporal
     
    endType
!----------------------------------------------------------------------------!

 

    !Variable
    public T_varCF
    Type T_varCF
        character(len=lenStr)         :: name=""                !Nome da variable
        character(len=lenStr)         :: completeName=""        !Nome longo
        character(len=lenStr)         :: standardName=""        !Nome estandar
        character(len=lenStr)         :: unit=""                !Unidade
        character(len=lenStr)         :: gridMapping=""         !Grid Mapping
        character(len=lenStr)         :: comment=""             !Comentario
        real                          :: valMax=valMaxD,valMin=valMinD
        real                          :: FillValue=FillValueRef
        real                          :: add_offset=0.
        real                          :: scale_factor=1.
                                                                !Valor máximo e mínimo
        integer                       :: nDimensions=-1         !Número de dimensións
        logical                       :: timeDim=.false.        !Flag da existencia da dimensión temporal
        real,pointer                  :: value2D(:,:)           !Array 2D de valores
        real,pointer                  :: value3D(:,:,:)         !Array 3D de valores

    endType
!----------------------------------------------------------------------------!


!Fin variabeis
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!Interfaces:


 

! Fin interfaces
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Subrutinas:


contains

!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_CAB
!
!  PURPOSE:  O proposito desta subrutina é crear a primeira parte de cabeceira
!            dun ficheiro NetCDF, e defini-las dimensións
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************

   Subroutine NCDF_ESC_CAB(ficheiro)
!    
     

! Variables
 
    
    integer                    :: ncid,status !Variable do netCDF
    integer                    :: xid, yid, zid, timid
    !integer                    :: tmp1
    

    type (T_fileCF)            :: ficheiro



! Corpo

    ! Comprobacións

        !!Comproba que o ficheiro ten extensión nc
        !tmp1=len(trim(ficheiro%name))
        !write (*,*) 'Extensio',(trim(ficheiro%name))(tmp1-3:tmp1)

        !Comproba que o ficheiro ten dimensions impostas

        if(ficheiro%nDimensions<0)                                                               &
        stop 'Cuando se crea un fichero se deben de especificar las dimensiones'
     

        if (len_trim(ficheiro%title)==0)  ficheiro%title="Titulo"  !Titulo
        if (len_trim(ficheiro%convention)==0)  ficheiro%convention="CF-1.0" ! Convencion
        if (len_trim(ficheiro%version)==0)  ficheiro%version="3.5.1" ! Version
        if (len_trim(ficheiro%history)==0) ficheiro%history="historia" ! historia
        if (len_trim(ficheiro%origin)==0)  ficheiro%origin="orixe" ! Orixe
        if (len_trim(ficheiro%institution)==0)  ficheiro%institution="MeteoGalicia" ! Institucion
        if (trim(ficheiro%references)=='')  ficheiro%references="http://www.meteogalicia.es" ! referencias


    !
 
    !Crea un dataset e entra en modo define
    
    status=nf90_create(path=ficheiro%name, cmode=nf90_write,ncid=ncid)
    if(status/=nf90_noerr) call erro(status,101)


    !Define as dimensions

    if(ficheiro%projection%pType==0) then 
     
       status=nf90_def_dim(ncid,'lon',ficheiro%xAxis%size,XID)                            !Lonxitude
       if(status/=nf90_noerr) call erro(status,102)

       status=nf90_def_dim(ncid,'lat',ficheiro%yAxis%size,YID)                            !Latitude
       if(status/=nf90_noerr) call erro(status,103)

       
    else
       status=nf90_def_dim(ncid,'x',ficheiro%xAxis%size,XID)                              !Lonxitude
       if(status/=nf90_noerr) call erro(status,104)

       status=nf90_def_dim(ncid,'y',ficheiro%yAxis%size,YID)                              !Latitude
       if(status/=nf90_noerr) call erro(status,105)

       
    endif

    if (ficheiro%timeDim) then 
        if (ficheiro%nDimensions==4) then
            status=nf90_def_dim(ncid,'lev',ficheiro%zAxis%size,ZID)                             !Latitude
            if(status/=nf90_noerr) call erro(status,106)
        endif    
        status=nf90_def_dim(ncid,'time',nf90_unlimited,TIMID)                               !Tempo
        if(status/=nf90_noerr) call erro(status,107)

    else
        if (ficheiro%nDimensions==3) then
            status=nf90_def_dim(ncid,'lev',ficheiro%zAxis%size,ZID)                             !Latitude
            if(status/=nf90_noerr) call erro(status,108)
        endif  
    endif


    
    

! Poe atributos globais


    status=nf90_put_att(ncid,nf90_global,'Title',ficheiro%title)  !Titulo
    if(status/=nf90_noerr) call erro(status,109)

    status=nf90_put_att(ncid,nf90_global,'Conventions', ficheiro%convention) !Convencion
    if(status/=nf90_noerr) call erro(status,110)
    
    status=nf90_put_att(ncid,nf90_global,'netcdf_version_id',ficheiro%version)!Version 
    if(status/=nf90_noerr) call erro(status,111)
    
    status=nf90_put_att(ncid,nf90_global,'history',ficheiro%history)         !Informacion
    if(status/=nf90_noerr) call erro(status,112)
         
    status=nf90_put_att(ncid,nf90_global,'date',ficheiro%idate)          !Data 
    if(status/=nf90_noerr) call erro(status,113)

    status=nf90_put_att(ncid,nf90_global,'source',ficheiro%origin)                  !Orixe dos datos
    if(status/=nf90_noerr) call erro(status,114)

    status=nf90_put_att(ncid,nf90_global,'institution',ficheiro%institution)                    !Institucion responsavel
    if(status/=nf90_noerr) call erro(status,115)
    
    status=nf90_put_att(ncid,nf90_global,'references',ficheiro%references)       !Referencias
    if(status/=nf90_noerr) call erro(status,116)

    !Poe os eixos

    if (len_trim(ficheiro%xAxis%name)==0)            ficheiro%xAxis%name="lon"
    if (len_trim(ficheiro%xAxis%completeName)==0)    ficheiro%xAxis%completeName="longitude"
    if (len_trim(ficheiro%xAxis%standardName)==0)    ficheiro%xAxis%standardName="longitude"
    if (len_trim(ficheiro%xAxis%unit)==0)            ficheiro%xAxis%unit="degrees_east"
    if (len_trim(ficheiro%xAxis%aType)==0)           ficheiro%xAxis%aType="X"

    ! Define os pontos dos eixos:

    
    
    call NCDF_ESC_EJE_ATT (ncid,xid,ficheiro%xAxis)

    if (len_trim(ficheiro%yAxis%name)==0)             ficheiro%yAxis%name="lat"
    if (len_trim(ficheiro%yAxis%completeName)==0)     ficheiro%yAxis%completeName="latitude"
    if (len_trim(ficheiro%yAxis%standardName)==0)     ficheiro%yAxis%standardName="latitude"
    if (len_trim(ficheiro%yAxis%unit)==0)             ficheiro%yAxis%unit="degrees_north"
    if (len_trim(ficheiro%yAxis%aType)==0)            ficheiro%yAxis%aType="Y"

    
    call NCDF_ESC_EJE_ATT (ncid,yid,ficheiro%yAxis)

    if ((ficheiro%timeDim .and. ficheiro%nDimensions==4) .or.                            &
       (.not.ficheiro%timeDim .and. ficheiro%nDimensions==3)) then

        if (len_trim(ficheiro%zAxis%name)==0)             ficheiro%zAxis%name="lev"
        if (len_trim(ficheiro%zAxis%completeName)==0)     ficheiro%zAxis%completeName="Vertical level"
        if (len_trim(ficheiro%zAxis%standardName)==0)      ficheiro%zAxis%standardName="model_level_number"
        if (len_trim(ficheiro%zAxis%unit)==0)             ficheiro%zAxis%unit=""
        if (len_trim(ficheiro%zAxis%aType)==0)            ficheiro%zAxis%aType="Z"
    
        call NCDF_ESC_EJE_ATT (ncid,zid,ficheiro%zAxis)
    endif

    if(ficheiro%timeDim) then 
        if (len_trim(ficheiro%tAxis%name)==0)            ficheiro%tAxis%name="time"
        if (len_trim(ficheiro%tAxis%completeName)==0)    ficheiro%tAxis%completeName="time"
        if (len_trim(ficheiro%tAxis%standardName)==0)    ficheiro%tAxis%standardName="time"
        if (len_trim(ficheiro%tAxis%unit)==0)         ficheiro%tAxis%unit="seconds since 2000-1-1"
        if (len_trim(ficheiro%tAxis%aType)==0)            ficheiro%tAxis%aType="T"

        call NCDF_ESC_EJE_ATT (ncid,timid,ficheiro%tAxis)
    endif


     

! Sae do modo define

    status = nf90_enddef(ncid)
    if(status/=nf90_noerr) call erro(status,117) 

!Pecha
 
    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,118) 

! Escribe a proxeccion

  if(ficheiro%projection%pType/=0)   call NCDF_ESC_PROJ(ficheiro)  
    

End Subroutine NCDF_ESC_CAB
!----------------------------------------------------------------------------!


!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_EJE_ATT
!
!  PURPOSE:  O proposito desta subrutina é crear os atributos dos eixos
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
    Subroutine NCDF_ESC_EJE_ATT(ncid,dimid,eixo)

!   Variabeis
    type(T_axisCF)          :: eixo 
    integer                 :: ncid, dimid,vareid
    
    integer                 :: status

!   Corpo

    status=nf90_def_var(ncid,eixo%name,nf90_float,dimid,vareid)        !ptos eixo  
    if(status/=nf90_noerr) call erro(status,201)
        
    status=nf90_put_att(ncid,vareid,'long_name',eixo%completeName)      !Nome completo
    if(status/=nf90_noerr) call erro(status,202)

    status=nf90_put_att(ncid,vareid,'standard_name',eixo%standardName)  !Nome estandar
    if(status/=nf90_noerr) call erro(status,203)

    status=nf90_put_att(ncid,vareid,'units',eixo%unit)      !Unidades
    if(status/=nf90_noerr) call erro(status,204)

    status=nf90_put_att(ncid,vareid,'valid_min',eixo%minVal)   !Valor minimo
    if(status/=nf90_noerr) call erro(status,205)

    status=nf90_put_att(ncid,vareid,'valid_max',eixo%maxVal)   !Valor maximo
    if(status/=nf90_noerr) call erro(status,206)

    status=nf90_put_att(ncid,vareid,'step', eixo%step)        !Paso de malla
    if(status/=nf90_noerr) call erro(status,207)

   ! status=nf90_put_att(ncid,vareid,'axis',eixo%tipo)                   !Eixo
   ! if(status/=nf90_noerr) call erro(status,208) 

End Subroutine NCDF_ESC_EJE_ATT
!-------------------------------------------------------------------------!

!----------------------------------------------------------------------------!

!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_EJE
!
!  PURPOSE:  O proposito desta subrutina é escrebe-los eixos
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
    Subroutine NCDF_ESC_EJE (ficheiro)

    !   Variabeis

        type(T_fileCF)                :: ficheiro
        integer                       :: ncid,varid
        integer                       :: status
    integer                       :: ndim, unlimdim


    
    !Abre o dataset
    
    status=nf90_open(ficheiro%name, nf90_write,ncid)
    if(status/=nf90_noerr) call erro(status,301)

    !Olla cantos eixos ten


    
    status=nf90_inquire(ncid,ndim,unlimdim)
    if(status/=nf90_noerr) call erro(status,302)

    !Escrebo os valores

    if (unlimdim/=-1) ndim=ndim-1

    if (ndim>=1) then

        !Eixo X
 status= nf90_inq_varid (ncid, ficheiro%xAxis%name, varid)
        if(status/=nf90_noerr) call erro(status,303)

 status= nf90_put_var(ncid, varid,ficheiro%xAxis%value)
 if(status/=nf90_noerr) call erro(status,304)

    endif


    if (ndim>=2) then
    
        !Eixo Y
 status= nf90_inq_varid (ncid, ficheiro%yAxis%name, varid)
        if(status/=nf90_noerr) call erro(status,305)

 status= nf90_put_var(ncid, varid,ficheiro%yAxis%value)
 if(status/=nf90_noerr) call erro(status,306)

    endif

    if (ndim>=3) then
        !Eixo Z
 status= nf90_inq_varid (ncid, ficheiro%zAxis%name, varid)
        if(status/=nf90_noerr) call erro(status,307)

 status= nf90_put_var(ncid, varid,ficheiro%zAxis%value)
 if(status/=nf90_noerr) call erro(status,308)

    endif





    !Pecha


    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,309) 

!


End Subroutine NCDF_ESC_EJE


!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_TIEMPO
!
!  PURPOSE:  O proposito desta subrutina é escreber o tempo, sexa todo o eixo
!            ou instante determinado
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
    Subroutine NCDF_ESC_TIEMPO (ficheiro,tempo)

    !   Variabeis
 

    type(T_fileCF)                :: ficheiro
    integer                       :: ncid,varid
    integer                       :: status
    real,optional                 :: tempo !Tempo en segundos do instante
    integer                       :: n
    integer                       :: timeId

    
    !Abre o dataset
    
    status=nf90_open(ficheiro%name, nf90_write,ncid)
    if(status/=nf90_noerr) call erro(status,401)

    !Escrebo os valores

    

    status= nf90_inq_varid (ncid, ficheiro%tAxis%name, varid)
    if(status/=nf90_noerr) call erro(status,402)

    if (present(tempo)) then 
        status=nf90_Inquire(ncid, unlimitedDimID=timeId)
        if(status/=nf90_noerr) call erro(status,403)

        status=nf90_Inquire_Dimension(ncid, timeId,len=n)
        if(status/=nf90_noerr) call erro(status,404)


 status= nf90_put_var(ncid, varid,tempo,(/n+1,1/))
 if(status/=nf90_noerr) call erro(status,405)
    else
        status= nf90_put_var(ncid, varid,ficheiro%tAxis%value)
 if(status/=nf90_noerr) call erro(status,406)
    endif

    !Pecha


    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,407) 

!


End Subroutine NCDF_ESC_TIEMPO

!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_VAR_CAB
!
!  PURPOSE:  O proposito desta subrutina é crear a  cabeceira dunha variable
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
Subroutine NCDF_ESC_VAR_CAB (ficheiro,var)

!   Variabeis
    type(T_varCF)             :: var 
    type(T_fileCF)            :: ficheiro

    integer                   :: ncid, varid
    integer                   :: xid,yid,zid,timid
    integer                   :: status

!   Corpo

!   Comprobacion de que a variable ten o numero de dimensions posto

    if (var%nDimensions<0) then
       write(*,'("ERRO: A variable ",a," non ten o numero de dimensions postas")') var%name
       write(*,'("Corrixa o codigo")')
       stop
    endif

!    Abre o dataset e entra no modo define
    
    status=nf90_open(ficheiro%name, nf90_write,ncid)
    if(status/=nf90_noerr) call erro(status,501)

    status=nf90_redef(ncid)
    if(status/=nf90_noerr) call erro(status,502)

    if(ficheiro%projection%pType==0) then  

        status=nf90_inq_varid(ncid,'lon',xid)
        if(status/=nf90_noerr) call erro(status,503)

    status=nf90_inq_varid(ncid,'lat',yid)
        if(status/=nf90_noerr) call erro(status,504)

    else

        status=nf90_inq_varid(ncid,'x',xid)
        if(status/=nf90_noerr) call erro(status,505)

    status=nf90_inq_varid(ncid,'y',yid)
        if(status/=nf90_noerr) call erro(status,506)
        
    endif

    if (var%timeDim) then 
        status=nf90_inq_varid(ncid,'time',timid)
        if(status/=nf90_noerr) call erro(status,507)

        if (var%nDimensions==4) then

            status=nf90_inq_varid(ncid,'lev',zid)
            if(status/=nf90_noerr) call erro(status,508)

            status=nf90_def_var(ncid,var%name,nf90_float,(/xid,yid,zid,timid/),varid) !ptos eixo x
            if(status/=nf90_noerr) call erro(status,509)
        
        elseif (var%nDimensions==3) then
           status=nf90_def_var(ncid,var%name,nf90_float,(/xid,yid,timid/),varid) !ptos eixo x
           if(status/=nf90_noerr) call erro(status,510)
        
        endif

    elseif (var%ndimensions==3) then
        status=nf90_inq_varid(ncid,'lev',zid)
        if(status/=nf90_noerr) call erro(status,511)
        status=nf90_def_var(ncid,var%name,nf90_float,(/xid,yid,zid/),varid) !ptos eixo x
        if(status/=nf90_noerr) call erro(status,512)
    else
        status=nf90_def_var(ncid,var%name,nf90_float,(/xid,yid/),varid) !ptos eixo x
        if(status/=nf90_noerr) call erro(status,513)

    endif


    


! Escrebe os atributos


        call NCDF_ESC_VAR_ATT(ncid,varid,var)

    ! Sae do modo define

    status = nf90_enddef(ncid)
    if(status/=nf90_noerr) call erro(status,514) 


!Pecha


    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,515) 

End Subroutine NCDF_ESC_VAR_CAB


!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_VAR_ATT
!
!  PURPOSE:  O proposito desta subrutina é crear os atributos da variable
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
Subroutine NCDF_ESC_VAR_ATT(ncid,varid,var)

    !   Variabeis
    implicit none

    type(T_varCF)           :: var 
    
    integer                 :: ncid, varid
    integer                 :: status

!   Corpo
    if(len_trim(var%completeName)/=0) then
        status=nf90_put_att(ncid,varid,'long_name',trim(var%completeName))      !Nome completo
        if(status/=nf90_noerr) call erro(status,601)
    endif

    if(len_trim(var%standardName)/=0) then
        status=nf90_put_att(ncid,varid,'standard_name',trim(var%standardName))  !Nome estandar
        if(status/=nf90_noerr) call erro(status,602)
    endif

    if(len_trim(var%unit)/=0) then
        status=nf90_put_att(ncid,varid,'units',trim(var%unit))      !Unidades
        if(status/=nf90_noerr) call erro(status,603)
    endif

    if (var%valMin/=valMinD) then
         status=nf90_put_att(ncid,varid,'valid_min',var%valMin)   !Valor minimo
         if(status/=nf90_noerr) call erro(status,604)
    endif
    
    if (var%valMax/=valMaxD) then
       status=nf90_put_att(ncid,varid,'valid_max',var%valMax)   !Valor maximo
       if(status/=nf90_noerr) call erro(status,605)
    endif

    if(var%gridMapping/="") then
       status=nf90_put_att(ncid,varid,'grid_mapping',trim(var%gridMapping))   
       if(status/=nf90_noerr) call erro(status,606)
    endif

    if(var%comment/="") then
       status=nf90_put_att(ncid,varid,'comment',trim(var%comment))   
       if(status/=nf90_noerr) call erro(status,607)
    endif

    if(var%FillValue/=FillValueRef) then
       status=nf90_put_att(ncid,varid,'_FillValue',var%FillValue)   
       if(status/=nf90_noerr) call erro(status,609)
    endif

    if(var%scale_factor/=1.) then
       status=nf90_put_att(ncid,varid,'scale_factor',var%scale_factor)   
       if(status/=nf90_noerr) call erro(status,610)
    endif

    if(var%add_offset/=0.)then
       status=nf90_put_att(ncid,varid,'add_offset',var%add_offset)   
       if(status/=nf90_noerr) call erro(status,610)
    endif

End Subroutine NCDF_ESC_VAR_ATT

!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_PROJ
!
!  PURPOSE:  O proposito desta subrutina é crear os atributos da proxeccion
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
    Subroutine NCDF_ESC_PROJ (ficheiro)

!   Variabeis
    type(T_fileCF)          :: ficheiro

    integer                 :: ncid, pid
    integer                 :: status

!   Corpo

!    Abre o dataset e entra no modo define
    
    status=nf90_open(ficheiro%name, nf90_write,ncid)
    if(status/=nf90_noerr) call erro(status,701)

    status=nf90_redef(ncid)
    if(status/=nf90_noerr) call erro(status,702)

    


    status=nf90_def_var(ncid,ficheiro%projection%name,nf90_int,pid) !ptos eixo x
    if(status/=nf90_noerr) call erro(status,703)

! Escrebe os atributos


    status=nf90_put_att(ncid,pid,'grid_mapping_name',ficheiro%projection%gmName)      !Nome completo
    if(status/=nf90_noerr) call erro(status,704)

    
    status=nf90_put_att(ncid,pid,'standard_parallel',ficheiro%projection%trueLat)  
    if(status/=nf90_noerr) call erro(status,705)
   

    status=nf90_put_att(ncid,pid,'longitude_of_central_meridian',ficheiro%projection%trueLon)      
    if(status/=nf90_noerr) call erro(status,706)

    status=nf90_put_att(ncid,pid,'longitude_of_projection_origin',ficheiro%projection%ctrLon)  
    if(status/=nf90_noerr) call erro(status,707)

    status=nf90_put_att(ncid,pid,'latitude_of_projection_origin',ficheiro%projection%ctrLat)   
    if(status/=nf90_noerr) call erro(status,708)

    

    ! Sae do modo define

    status = nf90_enddef(ncid)
    if(status/=nf90_noerr) call erro(status,709) 


!Pecha


    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,710) 

End Subroutine NCDF_ESC_PROJ
!-------------------------------------------------------------------------!



!*****************************************************************************
!
!  SUBRUTINA: NCDF_ESC_VAR
!
!  PURPOSE:  O proposito desta subrutina é escrebe-lo valor dunha variabel
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
Subroutine NCDF_ESC_VAR (ficheiro,var,nInstante)

    !   Variabeis

    type(T_fileCF)                :: ficheiro
    type(T_varCF)                 :: var
    integer                       :: ncid,varid
    integer                       :: status
    integer                       :: timeID
    integer,optional              :: nInstante
    integer                       :: n

    integer, dimension(nf90_max_var_dims) :: dimids

    
    !Abre o dataset
    
    status=nf90_open(ficheiro%name, nf90_write,ncid)
    if(status/=nf90_noerr) call erro(status,801)

    !Escrebo os valores


    status= nf90_inq_varid (ncid, var%name, varid)
    if(status/=nf90_noerr) call erro(status,802)

    status=nf90_inquire_variable (ncid,  varid,ndims=var%nDimensions)
    if(status/=nf90_noerr) call erro(status,803)
    status=nf90_inquire_variable (ncid,  varid,dimids=dimids(:var%nDimensions))
    if(status/=nf90_noerr) call erro(status,8031)

    status=nf90_Inquire(ncid,unlimitedDimID=timeId)
    if(status/=nf90_noerr) call erro(status,804)

    do n=1,var%nDimensions
      if (dimids(n)==timeId) var%timeDim=.true.
    enddo


    if (var%timeDim) then 

          if (present(nInstante)) then 
                n=nInstante
          else
                status=nf90_Inquire_Dimension(ncid, timeId,len=n)
                if(status/=nf90_noerr) call erro(status,805)
          endif

                if(var%nDimensions==3) then
                    status= nf90_put_var(ncid, varid,var%value2D,start=(/1,1,n/),&
                                         count=(/ficheiro%xAxis%size,ficheiro%yAxis%size,1/))
                     if(status/=nf90_noerr) call erro(status,806)
                elseif(var%nDimensions==4) then
             
                    status= nf90_put_var(ncid, varid,var%value3D,start=(/1,1,1,n/),&
                                         count=(/ficheiro%xAxis%size,ficheiro%yAxis%size,ficheiro%zAxis%size,1/))
                    if(status/=nf90_noerr) call erro(status,807)
                endif
     else
                if(var%nDimensions==3) then  
                    status= nf90_put_var(ncid, varid,var%value3D,start=(/1,1,1/),&
                                         count=(/ficheiro%xAxis%size,ficheiro%yAxis%size,ficheiro%zAxis%size/))
                    if(status/=nf90_noerr) call erro(status,808)
                elseif(var%nDimensions==2) then  
                    status= nf90_put_var(ncid, varid,var%value2D,start=(/1,1/),count=(/ficheiro%xAxis%size,ficheiro%yAxis%size/))
                    if(status/=nf90_noerr) call erro(status,809)
                endif

     endif
     
    !Pecha


    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,810) 

!


End Subroutine NCDF_ESC_VAR






!**************************************************************************
!**************************************************************************
!**************************************************************************

!      SUBRUTINAS DE LEITURA

!--------------------------------------------------------------------------!

!*****************************************************************************
!
!  SUBRUTINA: NCDF_LE_CAB
!
!  PURPOSE:  O proposito desta subrutina é ler a cabeceira
!            dun ficheiro NetCDF, para determinar todo dimensións
!
!            Nova versión 26/12/2004
!            Autor: Pedro Montero Vilar
!
!****************************************************************************

Subroutine NCDF_LE_CAB(ficheiro)
!    
     
    implicit none
! Variables

    type(T_fileCF)           :: ficheiro
 
    integer                    :: nDims,nVars,nGlobalAtts,unlimDimID
    
    character(len=lenStr)      :: nomeAtt,att
    
    integer                    :: ncid,status !Variable do netCDF
    
    integer                    :: n

    character(len=lenStr)      :: namevar
    integer                    :: xtype,natts
    integer,dimension(4)       :: dimids

    integer                    :: dimid
    !integer                    :: nx,ny,nt

    integer                    :: pid
    character(len=lenStr)      :: nomeTmp
    integer                    :: intTmp


! Corpo

   
 
    !Abre o dataset
    
    status=nf90_open(ficheiro%name, nf90_nowrite,ncid)
    if(status/=nf90_noerr) call erro(status,901)


    !Olla como e o ficheiro
     
    status=nf90_Inquire(ncid,nDims,nVars,nGlobalAtts,unlimDimID)                             
    if(status/=nf90_noerr) call erro(status,902)

    write(*,*)
    write(*,*)
    write(*,*) " Nome do ficheiro: ", ficheiro%name
    write(*,*) " N de Dimensions : ", nDims
    write(*,*) " N de Variabeis  : ", nVars
    write(*,*) " N de Atributos G: ", nGlobalAtts
    write(*,*) " N de Unlimite   : ", unlimDimID
    write(*,*) 
    write(*,*) 

    ficheiro%nDimensions=NDims
    ficheiro%nVariables=NVars-NDims
    if (unlimDimID/=-1) ficheiro%timeDim=.true.

    do n=1,NGlobalAtts

      status=nf90_inq_attname(ncid,nf90_Global,n,nomeAtt)
      if(status/=nf90_noerr) call erro(status,903)
      att=""
      status=nf90_get_att(ncid,nf90_Global,nomeAtt,att)
      
      write(*,'(a20," ----> ",a40)') trim(nomeAtt),trim(att)
    enddo
    write(*,*)
   
    do n=1,NVars

      status=nf90_inquire_Variable(ncid,n,namevar,xtype,ndims,dimids,natts)
      if(status/=nf90_noerr) call erro(status,904)


      if (namevar=="Lambert_Conformal") then 

        ficheiro%projection%name="Lambert_Conformal"
        ficheiro%projection%pType=2
      end if


      
      write(*,'(a20," ------> ",i5)') namevar,ndims
    enddo

    

    !Procura as dimensións
    do n=1,ficheiro%nDimensions
       status=nf90_Inquire_Dimension(ncid,n,ficheiro%nameDim(n))
    enddo

    


        !lon
           status=nf90_inq_dimid(ncid,ficheiro%nameDim(1),dimid)
        if(status/=nf90_noerr) call erro(status,906)
        status=nf90_Inquire_Dimension(ncid,dimid,len=ficheiro%xAxis%size)
        if(status/=nf90_noerr) call erro(status,907)

        !lat
        status=nf90_inq_dimid(ncid,ficheiro%nameDim(2),dimid)
        if(status/=nf90_noerr) call erro(status,908)
        status=nf90_Inquire_Dimension(ncid,dimid,len=ficheiro%yAxis%size)
        if(status/=nf90_noerr) call erro(status,909)
    
    !depth
    if (ficheiro%nDimensions>3) then
        
        status=nf90_inq_dimid(ncid,ficheiro%nameDim(3),dimid)
        if(status/=nf90_noerr) call erro(status,910)
        status=nf90_Inquire_Dimension(ncid,dimid,len=ficheiro%zAxis%size)
        if(status/=nf90_noerr) call erro(status,911)

        status=nf90_inq_dimid(ncid,ficheiro%nameDim(4),dimid)
       if(status/=nf90_noerr) call erro(status,912)

       status=nf90_Inquire_Dimension(ncid,dimid,len=ficheiro%tAxis%size)
       if(status/=nf90_noerr) call erro(status,913)

    else
       status=nf90_inq_dimid(ncid,ficheiro%nameDim(3),dimid)
       if(status/=nf90_noerr) call erro(status,914)

       status=nf90_Inquire_Dimension(ncid,dimid,len=ficheiro%tAxis%size)
       if(status/=nf90_noerr) call erro(status,915)

    endif
       

    

    write(*,*)
    write(*,*) "Dimensions da malla:"
    write(*,'("nx = ",i4,", ny = ",i4,", nt = ",i4)') ficheiro%xAxis%size,ficheiro%yAxis%size,ficheiro%tAxis%size
    if (ficheiro%nDimensions>3) write(*,'("nz = ",i4)') ficheiro%zAxis%size
    write(*,*)



    ! Procura os atributos principais

    ficheiro%title=""
    status=nf90_get_att(ncid,nf90_global,'Title',ficheiro%title)  !Titulo
    if(status/=nf90_noerr) status=nf90_get_att(ncid,nf90_global,'title',ficheiro%title) 
 !   if(status/=nf90_noerr) call erro(status,916)

    ficheiro%convention=""
    status=nf90_get_att(ncid,nf90_global,'Conventions', ficheiro%convention)    !Convencion
    if(status/=nf90_noerr) status=nf90_get_att(ncid,nf90_global,'conventions', ficheiro%convention)
    !if(status/=nf90_noerr) call erro(status,917)
    
    ficheiro%version=""
    status=nf90_get_att(ncid,nf90_global,'netcdf_version_id',ficheiro%version)  !Version 
    !if(status/=nf90_noerr) call erro(status,918)
    
    ficheiro%history=""
    status=nf90_get_att(ncid,nf90_global,'history',ficheiro%history)           !Informacion
   ! if(status/=nf90_noerr) call erro(status,919)
         
    status=nf90_get_att(ncid,nf90_global,'date',ficheiro%idate)                 !Data 
   ! if(status/=nf90_noerr) call erro(status,920)

    ficheiro%origin=""
    status=nf90_get_att(ncid,nf90_global,'source',ficheiro%origin)               !Orixe dos datos
   ! if(status/=nf90_noerr) call erro(status,921)

     ficheiro%institution=""
    status=nf90_get_att(ncid,nf90_global,'institution',ficheiro%institution)    !Institucion responsavel
    !if(status/=nf90_noerr) call erro(status,922)
    
    ficheiro%references=""
    status=nf90_get_att(ncid,nf90_global,'references',ficheiro%references)     !Referencias
    !if(status/=nf90_noerr) call erro(status,923)

! Procura a proxeccion

   
       
    if(ficheiro%projection%pType==2) then
    ! Escrebe os atributos
      status=nf90_inq_varid(ncid,ficheiro%projection%name,pid)
      if(status/=nf90_noerr) call erro(status,924)

      ficheiro%projection%gmName=""
      status=nf90_get_att(ncid,pid,'grid_mapping_name',ficheiro%projection%gmName)      !Nome completo
      if(status/=nf90_noerr) call erro(status,925)

     
      status=nf90_get_att(ncid,pid,'standard_parallel',ficheiro%projection%trueLat(1))  !Nome estandar
      if(status/=nf90_noerr) call erro(status,926)
      ficheiro%projection%trueLat(2)=ficheiro%projection%trueLat(1)
    
      status=nf90_get_att(ncid,pid,'longitude_of_central_meridian',ficheiro%projection%trueLon)      !Unidades
      if(status/=nf90_noerr) call erro(status,927)

      status=nf90_get_att(ncid,pid,'longitude_of_projection_origin',ficheiro%projection%ctrLon)   !Valor minimo
      if(status/=nf90_noerr) ficheiro%projection%ctrLon=ficheiro%projection%trueLon

      status=nf90_get_att(ncid,pid,'latitude_of_projection_origin',ficheiro%projection%ctrLat)   !Valor minimo
      if(status/=nf90_noerr) call erro(status,928)
    endif

 
!le os  eixos
      

        status=nf90_Inquire_Variable(ncid,1,NomeTmp)
        if(status/=nf90_noerr) call erro(status,929)


        if (NomeTmp=="lat".or.NomeTmp=="latitude".or.nomeTmp=="y") then !Por se estan cambiadas de orden
            ficheiro%yAxis%name=nomeTmp
            call NCDF_LE_EJE(ncid,ficheiro%yAxis)
        else
            
            ficheiro%xAxis%name=nomeTmp
            call NCDF_LE_EJE(ncid,ficheiro%xAxis)
        endif
       

        status=nf90_Inquire_Variable(ncid,2,NomeTmp)
        if(status/=nf90_noerr) call erro(status,930)

        if (NomeTmp=="lon".or.NomeTmp=="longitude".or.nomeTmp=="x") then !Por se estan cambiadas de orden
            ficheiro%xAxis%name=nomeTmp
            call NCDF_LE_EJE(ncid,ficheiro%xAxis)
        else
            
            ficheiro%yAxis%name=nomeTmp
            call NCDF_LE_EJE(ncid,ficheiro%yAxis)
        endif
        
        if(ficheiro%nDimensions> 3) then 
          
           status=nf90_Inquire_Variable(ncid,3,NomeTmp,intTmp,ndims)
           if(status/=nf90_noerr) call erro(status,931)
           if(ndims==1) then
              ficheiro%zAxis%name=nomeTmp
              call NCDF_LE_EJE(ncid,ficheiro%zAxis)
           endif

           status=nf90_Inquire_Variable(ncid,4,NomeTmp,intTmp,ndims)
           if(status/=nf90_noerr) call erro(status,932)
           if(ndims==1) then
              ficheiro%tAxis%name=nomeTmp
              call NCDF_LE_EJE(ncid,ficheiro%tAxis)
           endif
        else
         
           status=nf90_Inquire_Variable(ncid,3,NomeTmp,intTmp,ndims)
           if(status/=nf90_noerr) call erro(status,933)
           if(ndims==1) then
              ficheiro%tAxis%name=nomeTmp
              call NCDF_LE_EJE(ncid,ficheiro%tAxis)
           endif

        endif 
 
    

   


!Pecha
 
    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,934) 



End Subroutine NCDF_LE_CAB

!*****************************************************************************
!
!  SUBRUTINA: NCDF_LE_EJE
!
!  PURPOSE:  O proposito desta subrutina é lear os atributos e os valores dos eixos
!
!            Nova versión 07/01/2005
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
Subroutine NCDF_LE_EJE (ncid,eixo)

!   Variabeis
    implicit none

    type(T_axisCF)          :: eixo 
    integer                 :: ncid, vareid
    
    integer                 :: status

   

    integer                  :: nDim

!   Corpo

    status=nf90_inq_varid(ncid,eixo%name,vareid)
    if(status/=nf90_noerr) call erro(status,1001)

    
    eixo%completeName=""
    status=nf90_get_att(ncid,vareid,'long_name',eixo%completeName)      !Nome completo
    !if(status/=nf90_noerr) call erro(status,1002)

    eixo%standardName=""
    status=nf90_get_att(ncid,vareid,'standard_name',eixo%standardName)  !Nome estandar
    !if(status/=nf90_noerr) call erro(status,1003)

    eixo%unit=""
    status=nf90_get_att(ncid,vareid,'units',eixo%unit)      !Unidades
   ! if(status/=nf90_noerr) call erro(status,1004)

    status=nf90_get_att(ncid,vareid,'valid_min',eixo%minVal)   !Valor minimo
    !if(status/=nf90_noerr) call erro(status,1005)

    status=nf90_get_att(ncid,vareid,'valid_max',eixo%maxVal)   !Valor maximo
    !if(status/=nf90_noerr) call erro(status,1006)

    status=nf90_get_att(ncid,vareid,'step', eixo%step)        !Paso de malla
    !if(status/=nf90_noerr) call erro(status,1007)

   ! status=nf90_get_att(ncid,vareid,'axis',eixo%aType)                   !Eixo
   ! if(status/=nf90_noerr) call erro(status,1008) 

    
    

    if (.not. associated(eixo%value)) allocate(eixo%value(eixo%size))

    status=nf90_inquire_variable(ncid,vareid, ndims=nDim)
    if(status/=nf90_noerr) call erro(status,1009)
    if(nDim==1) then 
      status= nf90_get_var(ncid, vareid,eixo%value)
      if(status/=nf90_noerr) call erro(status,1010)
    elseif(nDim==2) then
      if (vareid==1) then !x

        status= nf90_get_var(ncid, vareid,eixo%value,start=(/1,1/),count=(/eixo%size,1/))
        if(status/=nf90_noerr) call erro(status,1011)
      elseif(vareid==2) then !y
        status= nf90_get_var(ncid, vareid,eixo%value,start=(/1,1/),count=(/1,eixo%size/))
        if(status/=nf90_noerr) call erro(status,1012)
      endif 
    endif




End Subroutine NCDF_LE_EJE




!*****************************************************************************
!
!  SUBRUTINA: NCDF_LE_VAR_CAB
!
!  PURPOSE:  O proposito desta subrutina é ler a cabeceira dunha variable,
!            sexa xa con o nome ou sinalando o identificador (opcional)
!
!            Nova versión 10/01/2005
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
Subroutine NCDF_LE_VAR_CAB(ficheiro, var, varidin)

!   Variabeis
    type(T_fileCF)          :: ficheiro 
    type(T_varCF)           :: var
    integer                 :: ncid, varid,timeid

    integer,optional        :: varidin

    integer, dimension(nf90_max_var_dims) :: varDimIds
    integer                 :: status,n

!   Corpo
   

    !Abre o dataset
    
    status=nf90_open(ficheiro%name, nf90_nowrite,ncid)
    if(status/=nf90_noerr) call erro(status,1100)

    !le Atributos
    if (.not.present(varidin)) then
       
        status=nf90_inq_varid(ncid,trim(var%name),varid)
        if(status/=nf90_noerr) call erro(status,1101)
    else
        varid=varidin

       
    endif
    status=nf90_Inquire_Variable(ncid, varid,name=var%name, ndims=var%nDimensions)
    if(status/=nf90_noerr) call erro(status,1102)

   var%completeName=""
   status=nf90_get_att(ncid,varid,'long_name',var%completeName)      !Nome completo
   ! if(status/=nf90_noerr) call erro(status,1103)

    var%standardName=""
    status=nf90_get_att(ncid,varid,'standard_name',var%standardName)  !Nome estandar
   ! if(status/=nf90_noerr) call erro(status,1104)

    var%unit=""
    status=nf90_get_att(ncid,varid,'units',var%unit)      !Unidades
   !if(status/=nf90_noerr) call erro(status,1105)

    status=nf90_get_att(ncid,varid,'valid_min',var%valMin)   !Valor minimo
    !if(status/=nf90_noerr) call erro(status,1106)

    status=nf90_get_att(ncid,varid,'valid_max',var%valMax)   !Valor maximo
    !if(status/=nf90_noerr) call erro(status,1107)

    status=nf90_get_att(ncid,varid,'_FillValue',var%FillValue)

    status=nf90_get_att(ncid,varid,'scale_factor',var%scale_factor)

    status=nf90_get_att(ncid,varid,'add_offset',var%add_offset)

    !Obsolete mais se deixa 

    status=nf90_get_att(ncid,varid,'missing_value',var%FillValue)

    !Le as dimensions
    status=nf90_inquire_variable (ncid,  varid,dimids=varDimIds(:var%nDimensions))
    if(status/=nf90_noerr) call erro(status,1108)

    status=nf90_Inquire(ncid,unlimitedDimID=timeId)
    if(status/=nf90_noerr) call erro(status,1109)

   !
       do n=1,var%nDimensions
           if(varDimIds(n)==timeId) then
              var%timeDim=.true.
           else
              var%timeDim=.false.
           endif
       enddo
    !Pecha
 
    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,1110) 



 
End Subroutine NCDF_LE_VAR_CAB


!*****************************************************************************
!
!  SUBRUTINA: NCDF_LE_VAR
!
!  PURPOSE:  O proposito desta subrutina é ler os valores dunha
!            variabel.
!            É opcional dar o identificador para identifica-la variabel
!            Tamén é opcional dar un instante para ler ou se non, le todos os
!            valores
!
!            Nova versión 10/01/2005
!            Autor: Pedro Montero Vilar
!
!****************************************************************************
Subroutine NCDF_LE_VAR (ficheiro, var, varidin,nTime)

!   Variabeis
    type(T_fileCF)          :: ficheiro 
    type(T_varCF)           :: var
   
    integer                 :: ncid, varid,timeid
    integer,optional        :: nTime
    integer,optional        :: varidin

    integer, dimension(nf90_max_var_dims) :: varDimIds
    integer                 :: status,n

!   Corpo

    !Abre o dataset
    
    status=nf90_open(ficheiro%name, nf90_nowrite,ncid)
    if(status/=nf90_noerr) call erro(status,1200)

    if (.not.present(varidin)) then
        if(len_trim(var%name)==0) stop 'ERRO: Para ler unha variable sen id tenlle que dar un nome'

        status=nf90_inq_varid(ncid,trim(var%name),varid)
        if(status/=nf90_noerr) call erro(status,1201)
    else
        varid=varidin
    endif

    status=nf90_inquire_variable (ncid,  varid,ndims=var%nDimensions)
    if(status/=nf90_noerr) call erro(status,1202)

    status=nf90_Inquire(ncid,unlimitedDimID=timeId)
    if(status/=nf90_noerr) call erro(status,1203)

    !Le as dimensions
    status=nf90_inquire_variable (ncid,  varid,dimids=varDimIds(:var%nDimensions))
    if(status/=nf90_noerr) call erro(status,120301)

    status=nf90_Inquire(ncid,unlimitedDimID=timeId)
    if(status/=nf90_noerr) call erro(status,120302)

   !
       do n=1,var%nDimensions
           if(varDimIds(n)==timeId) var%timeDim=.true.
       enddo


  

    if (.not.var%timeDim) then 

        if (var%nDimensions==2) then 

            if (.not. associated(var%value2D)) allocate(var%value2D(ficheiro%xAxis%size,ficheiro%yAxis%size))
            status= nf90_get_var(ncid, varid,var%value2D)
            if(status/=nf90_noerr) call erro(status,1204)
        elseif(var%nDimensions==3) then
            if (.not. associated(var%value3D)) allocate(var%value3D(ficheiro%xAxis%size,ficheiro%yAxis%size,ficheiro%zAxis%size))
            status= nf90_get_var(ncid, varid,var%value2D)
            if(status/=nf90_noerr) call erro(status,1205)
        endif
    else
        if(.not.present(nTime)) stop 'ERRO: Debe especificar un instante de tempo para ler unha variable'

        if (var%nDimensions==3) then 
            if (.not. associated(var%value2D)) allocate(var%value2D(ficheiro%xAxis%size,ficheiro%yAxis%size))
            status= nf90_get_var(ncid, varid,var%value2D,start=(/1,1,nTime/),count=(/ficheiro%xAxis%size,ficheiro%yAxis%size,1/))
            if(status/=nf90_noerr) call erro(status,1206)
        elseif(var%nDimensions==4) then
            if (.not. associated(var%value3D)) allocate(var%value3D(ficheiro%xAxis%size,ficheiro%yAxis%size,ficheiro%zAxis%size))
            status= nf90_get_var(ncid, varid,var%value3D,start=(/1,1,1,nTime/), &
              count=(/ficheiro%xAxis%size,ficheiro%yAxis%size,ficheiro%zAxis%size,1/))
            if(status/=nf90_noerr) call erro(status,1207)
        endif
    endif


    !Pecha
 
    status =nf90_close(ncid)
    if(status/=nf90_noerr) call erro(status,1208) 


End Subroutine NCDF_LE_VAR

!-------------------------------------------------------------------------!
!-------------------------------------------------------------------------!
!********************************************************************
!                                                                   *
!       NETCDF: Error: Devuelve los mensajes de NetCDF de error     *
!                                                                   *
!********************************************************************
subroutine erro(status,codigo)
!
      integer, intent(in) :: status,codigo
!
      if(status /=nf90_noerr) then
        write(*,*) trim(nf90_strerror(status))
        write(*,*) "codigo:",codigo
        stop 'Stopped'
      endif
end subroutine erro


!***************************************************************
End Module

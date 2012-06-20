!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : BenthicEcology
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2012
! REVISION      : Isabella
! DESCRIPTION   : Module to compute simple benthic ecology processes
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
!DataFile
!   DT                          : real           [3600]         !Time step compute biogeochemical processes
!   PELAGIC_MODEL               : string                        Name of the pelagic model coupled with benthicecology 
!---------------------------------Biochemical processes parameters  --------------------------------------------------   
!Type and Name                    Label (from file)   value         unit                 Description
!BioChemPar%KmaxNit               KMAXNIT             0.1           [1/day]              Maximum nitrification rate
!BioChemPar%KNitIni               KNITINI             0.1           [-]                  Nitrification Inibition factor
!BioChemPar%Kmin                  KMIN                0.1           [1/day]              Mineralisation rate
!BioChemPar%Min_ON_Conv           MIN_ON_CONV         0.06          [-]                  Ratio Oxygen/Nitrogen in mineralization
!BioChemPar%Nit_ON_Conv           NIT_ON_CONV         4.57          [-]                  Ratio Oxygen to Nitrogen in nitrification
!---------------------------------Producer parameters  ----------------------------------------------------------------------------
!Producer%RespirationRate        RESPIRATION_RATE    User defined   [1/day]              Producer Respiration rate       
!Producer%MortalityRate          MORTALITY_RATE      User defined   [1/day]              Producer Mortality rate  
!Producer%Vmax                   VMAX                User defined   [1/day]              Producer Uptake rate
!Producer%alpha                  ALPHA               User defined   [m2/(Watt*day)]      Producer slope of the PI curve    
!Producer%KN                     KN                  User defined   [KgN/m3]             Half saturation Constant for N Uptake
!Producer%KP                     KP                  User defined   [KgP/m3]             Half saturation Constant for P Uptake
!Producer%NCratio                NCRATIO             User Defined   [-]                  Producer N/C ratio
!Producer%PCratio                PCRATIO             User Defined   [-]                  Producer P/C ratio
!---------------------------------Consumer parameters  -----------------------------------------------------------------------------------------
!Consumer%RespirationRate        RESPIRATION_RATE    User defined   [1/day]              Consumer Respiration rate
!Consumer%MortalityRate          MORTALITY_RATE      User defined   [1/day]              Consumer Mortality rate
!Consumer%NCratio                NCRATIO             User Defined   [-]                  Consumer N/C ratio
!Consumer%Vmax                   VMAX                User defined   [1/day]              maximum specific uptake at 10ºC
!Consumer%Ks                     KS                  User defined   [KgN/m3]             Half saturation Constant for prey Uptake (pelagic prey)
!                                                                   [KgN/m2]             Half saturation Constant for prey Uptake (benthic Prey)
!Consumer%Ass_Efic               ASS_EFIC            User defined   [-]                  Assimilation efficiency
!---------------------------------Prey parameters  ----------------------------------------------------------------------------------------------
!Prey%NCRatio                    NCRATIO             User Defined   [-]                           Prey N/C ratio
!Prey%Name                       NAME                User Defined   [-]                           Prey Name

Module ModuleBenthicEcology

    use ModuleGlobalData
    use ModuleEnterData

    implicit none

    private 

    !Subroutines---------------------------------------------------------------

    !Constructor
    
    public  :: ConstructBenthicEcology
    private ::      AllocateInstance
    private ::      ReadData
    private ::          ConstructGlobalVariables
    private ::          ReadOrganicMatterParameters
    private ::          ReadOxygenParameters
    private ::          ReadNitrogenParameters
    private ::          ReadPhosphorusParameters
    private ::          ReadSilicaParameters
    private ::          ConstructProducers
    private ::              AddProducer
    private ::              ConstructProducerParameters
    private ::          ConstructConsumers
    private ::              AddConsumer
    private ::              ConstructConsumerParameters
   ! private ::          ConstructDecomposers
   ! private ::              AddDecomposer
   ! private ::              ConstructDecomposerParameters
    private ::              ConstructGrazing
    private ::                  AddPrey
    private ::                  ConstructPrey
    private ::          PropertyIndexNumber
    private ::          ConstructPropertyList
        
        
    !Selector
    public  :: GetDTBenthicEcology
    public  :: GetBenthicEcologyPropertyList
    public  :: GetBenthicEcologySize
    public  :: GetBenthicEcologyPropIndex
    public  :: UnGetBenthicEcology
    public  :: GetBenthicEcologyRateFlux  
    public  :: UnGetBenthicEcologyRateFlux  
                     
    
    !Modifier
    public  :: ModifyBenthicEcology
    private ::      ComputeBenthicProducers
    private ::      ComputeBenthicConsumers
    !private ::      Decomposers
   
    
    !Destructor
    public  :: KillBenthicEcology                                                     
    private ::      DeAllocateInstance

    !Management
    private ::      Ready
    private ::          LocateObjBenthicEcology 
    
    
    
  private :: T_PropIndex
    type       T_PropIndex
         
        integer                                 :: Ammonia      = null_int        
        integer                                 :: Nitrate      = null_int
        integer                                 :: POC          = null_int                 
        integer                                 :: PON          = null_int 
        integer                                 :: POP          = null_int
        integer                                 :: Oxygen       = null_int       
        integer                                 :: Phyto        = null_int
        integer                                 :: Phosphate    = null_int
        integer                                 :: Bacteria     = null_int
        integer                                 :: PONr         = null_int
        integer                                 :: DONnr        = null_int
        integer                                 :: Diatoms      = null_int
        integer                                 :: PON1         = null_int
        integer                                 :: PON2         = null_int
        integer                                 :: PON3         = null_int
        integer                                 :: PON4         = null_int
        integer                                 :: PON5         = null_int
        integer                                 :: POP1         = null_int
        integer                                 :: POP2         = null_int
        integer                                 :: POP3         = null_int
        integer                                 :: POP4         = null_int
        integer                                 :: POP5         = null_int
        integer                                 :: DissolvedSilica   = null_int
        integer                                 :: BioSilica         = null_int
    end type T_PropIndex
    
    
     private :: T_External
    type       T_External
        real, pointer, dimension(:  )       :: Temperature
        real, pointer, dimension(:,:)       :: MassInKgFromWater
        real, pointer, dimension(:  )       :: Sediment
        real, pointer, dimension(:  )       :: WaterVolume
        real, pointer, dimension(:  )       :: CellArea
        real, pointer, dimension(:  )       :: ShortWaveTop
        real, pointer, dimension(:  )       :: ShortWaveAverage
        real, pointer, dimension(:  )       :: LightExtCoefField
        real, pointer, dimension(:  )       :: Thickness
        real, pointer, dimension(:,:)       :: Mass    
    end type T_External
  
  
  type     T_ComputeOptions
        logical                             :: Nitrogen                      = .false.
        logical                             :: Phosphorus                    = .false.
        logical                             :: Producers                     = .false.     
        logical                             :: Consumers                     = .false. 
        logical                             :: Silica                        = .false.
        logical                             :: Diatoms                       = .false.  ! CONTROLLARE
        logical                             :: Phyto                         = .false.  ! ! CONTROLLARE
        logical                             :: Pompools                      = .false.
  end type     T_ComputeOptions
  
    type     T_OrganicMatter
        real                                        :: NC_Ratio         = null_real 
        real                                        :: PC_Ratio         = null_real 
    end type T_OrganicMatter

    type     T_Nitrogen
        real                                        :: PONDecayRate     = null_real 
        real                                        :: PONDecayTFactor  = null_real
    end type T_Nitrogen
    
    type     T_Phosphorus
        real                                        :: POPDecayRate     = null_real 
        real                                        :: POPDecayTFactor  = null_real
    end type T_Phosphorus
    
    type     T_Oxygen
        real                                        :: Minimum          = null_real 
    end type T_Oxygen

    type     T_Silica
        real                                        :: BioSiDecayRate   = null_real 
    end type T_Silica

    
   private :: T_BenthicEcology
    type       T_BenthicEcology
        integer                                      :: InstanceID
        type (T_Size1D)                              :: Prop
        type (T_Size1D)                              :: Size
        real,    dimension(:,:,:),  pointer          :: Matrix
        real                                         :: DT, DTDay
        character(len=StringLength)                  :: PelagicModel
        integer, dimension(:), pointer               :: PropertyList
        type(T_Size1D       )                        :: Array
        type(T_PropIndex    )                        :: PropIndex
        type(T_Producer     ), pointer               :: FirstProducer
        type(T_Consumer     ), pointer               :: FirstConsumer
        type(T_External     )                        :: ExternalVar
        integer                                      :: ObjEnterData = 0
        type(T_BenthicEcology), pointer              :: Next
        type(T_ComputeOptions)                       :: ComputeOptions
        type(T_Silica        )                       :: Silica
        type(T_Oxygen        )                       :: Oxygen
        type(T_OrganicMatter )                       :: OrganicMatter
        type(T_Nitrogen      )                       :: Nitrogen
        type(T_Phosphorus    )                       :: Phosphorus
    end type  T_BenthicEcology
  
  
  
  
  private :: T_ID
        type       T_ID
            integer                         :: ID
            character(len=StringLength)     :: Name
            character(len=StringLength)     :: Description
        end type   T_ID
        
  
   
  


    private :: T_PoolIndex
        type       T_PoolIndex
            integer                                 :: Carbon      = null_int
            integer                                 :: Nitrogen    = null_int         
            integer                                 :: Phosphorus  = null_int         
            integer                                 :: Silica      = null_int
            integer                                 :: Chlorophyll = null_int        
        end type T_PoolIndex
   
   
   private :: T_Producer
   type     T_Producer
   type(T_ID)                                       :: ID
   type(T_PoolIndex)                                :: PoolIndex
        real                                        :: RespirationRate  = null_real   
        real                                        :: MortalityRate    = null_real 
        real                                        :: Vmax             = null_real  
        real                                        :: alpha            = null_real
        real                                        :: KN               = null_real 
        real                                        :: KP               = null_real 
        real                                        :: KTRANS           = null_real 
        real                                        :: KLIGHT           = null_real
        real                                        :: NCratio          = null_real   
        real                                        :: PCratio          = null_real
   type(T_Producer    ), pointer                    :: Next 
    end type T_Producer
 
 private :: T_Prey
        type   T_Prey
            type(T_ID)                  :: ID
            type(T_Prey), pointer       :: Next
            real                        :: NCratio
            real                        :: PCratio
            logical                     :: Use_Carbon
            logical                     :: Use_Nitrogen
            logical                     :: Use_Phosphorus
            logical                     :: ParticulateWaterFood
        end type T_Prey


    private :: T_Grazing
        type   T_Grazing
            real                        :: Vmax             = FillValueReal         !maximum specific uptake at 10ºC
            real                        :: Ks               = FillValueReal         !half saturation value for uptake
            real                        :: Ass_Efic         = FillValueReal         !assimilation efficiency
            real                        :: NCRatio          = FillValueReal
            type(T_Prey), pointer       :: FirstPrey
        end type T_Grazing
 
   private :: T_Consumer 
   type     T_Consumer
        type(T_ID)                                       :: ID
        type(T_Prey), pointer                            :: FirstPrey 
        type(T_PoolIndex)                                :: PoolIndex
        type(T_Consumer    ), pointer                    :: Next 
        type(T_Grazing     )                             :: Grazing
        real                                             :: RespirationRate  = null_real  
        real                                             :: O2HX             = null_real 
        real                                             :: O2QX             = null_real
        real                                             :: MortalityRate    = null_real
        real                                             :: Td               = null_real ! Time for death of 99% individuals
        real                                             :: BETA20           = null_real
        real                                             :: TemperatureFactor= null_real
        real                                             :: FILT20           = null_real ! l/day/gC
        real                                             :: RESP20          = null_real  ! [1/day]
        real                                             :: NCratio          = null_real
        real                                             :: PCratio          = null_real    
   end type T_Consumer
 

    
    
    
    

!Global Module Variables
    type (T_BenthicEcology), pointer                          :: FirstObjBenthicEcology
    type (T_BenthicEcology), pointer                          :: Me

  !--------------------------------------------------------------------------
    
     contains


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 
   subroutine ConstructBenthicEcology(ObjBenthicEcologyID, FileName, ILB,IUB, STAT)

        !Arguments---------------------------------------------------------------
        integer                                         :: ObjBenthicEcologyID  
        character(len=StringLength)                     :: FileName
        integer, optional, intent(OUT)                  :: STAT
        integer          , intent(IN )                  :: ILB, IUB

        !External----------------------------------------------------------------
        integer                                         :: ready_, STAT_CALL

        !Local-------------------------------------------------------------------
        integer                                         :: STAT_

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        !Assures nullification of the global variable
        if (.not. ModuleIsRegistered(mBenthicEcology_)) then
            nullify (FirstObjBenthicEcology)
            call RegisterModule (mBenthicEcology_) 
        endif
        
        call Ready(ObjBenthicEcologyID, ready_)    

cd0 :   if (ready_ .EQ. OFF_ERR_) then

            call AllocateInstance
            
            Me%Size%ILB = ILB
            Me%Size%IUB = IUB
            
            call ConstructEnterData(Me%ObjEnterData, FileName, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBenthicEcology - ModuleBenthicEcology - ERROR #1'

            call ReadData

            call PropertyIndexNumber
        
            call ConstructPropertyList
            
            call ConstructMatrix

            call KillEnterData(Me%ObjEnterData, STAT = STAT_CALL) 
            if (STAT_CALL .NE. SUCCESS_) stop 'ConstructBenthicEcology - ModuleBenthicEcology - ERROR #2'

            !Returns ID
            ObjBenthicEcologyID          = Me%InstanceID

            STAT_ = SUCCESS_
        else 
            
            stop 'ModuleBenthicEcology - ConstructBenthicEcology - ERROR #1' 

        end if cd0


        if (present(STAT)) STAT = STAT_

        !----------------------------------------------------------------------

    end subroutine ConstructBenthicEcology
   
   
   
    !----------------------------------------------------------------------
    
        subroutine AllocateInstance

        !Arguments-------------------------------------------------------------
                                                    
        !Local-----------------------------------------------------------------
        type (T_BenthicEcology), pointer                         :: NewObjBenthicEcology
        type (T_BenthicEcology), pointer                         :: PreviousObjBenthicEcology


        !Allocates new instance
        allocate (NewObjBenthicEcology)
        nullify  (NewObjBenthicEcology%Next)

        !Insert New Instance into list and makes Current point to it
        if (.not. associated(FirstObjBenthicEcology)) then
            FirstObjBenthicEcology         => NewObjBenthicEcology
            Me                   => NewObjBenthicEcology
        else
            PreviousObjBenthicEcology      => FirstObjBenthicEcology
            Me                   => FirstObjBenthicEcology%Next
            do while (associated(Me))
                PreviousObjBenthicEcology  => Me
                Me               => Me%Next
            enddo
            Me                   => NewObjBenthicEcology
            PreviousObjBenthicEcology%Next => NewObjBenthicEcology
        endif

        Me%InstanceID = RegisterNewInstance (mBenthicEcology_)

    end subroutine AllocateInstance


!____________________________________________________________________________________________________
!____________________________________________________________________


subroutine ReadData

        !Arguments-------------------------------------------------------------
                                                           
        !Local-----------------------------------------------------------------
        
        call ConstructGlobalVariables
        
        call ReadOrganicMatterParameters
        call ReadOxygenParameters
        call ReadNitrogenParameters
        call ReadPhosphorusParameters
        call ReadSilicaParameters
        
        call ConstructProducers

        call ConstructConsumers
        
       ! call ConstructDecomposers

    end subroutine ReadData
    
    
    
       subroutine PropertyIndexNumber

        !Arguments-------------------------------------------------------------


        !Local-----------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
       ! type(T_Decomposer),    pointer             :: Decomposer
        integer                                    :: Index
        !Local-----------------------------------------------------------------
        
        Me%Prop%ILB = 1
        Me%Prop%IUB = 0

        Index               = 0
        


        !Producer index number      
            Producer => Me%FirstProducer
            do while(associated(Producer))
                
                Index                               = Index + 1
                Producer%PoolIndex%Carbon           = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                Producer%PoolIndex%Nitrogen         = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                Index                               = Index + 1
                Producer%PoolIndex%Phosphorus       = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                Producer => Producer%Next
            end do
        
        !Consumer index number
            Consumer => Me%FirstConsumer
            do while(associated(Consumer))
                
                Index                               = Index + 1
                Consumer%PoolIndex%Carbon           = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                Index                               = Index + 1
                Consumer%PoolIndex%Nitrogen         = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1
                
                
                Index                               = Index + 1
                Consumer%PoolIndex%Phosphorus       = Index
                Me%Prop%IUB                         = Me%Prop%IUB + 1

                Consumer => Consumer%Next
            end do
   
        !Decomposer index number
          !  Decomposer => Me%FirstDecomposer
          !  do while(associated(Decomposer))
                
            !    Index                               = Index + 1
             !    Decomposer%PoolIndex%Carbon         = Index
              !   Me%Prop%IUB                         = Me%Prop%IUB + 1


            !    Decomposer => Decomposer%Next
           ! end do
        
        Me%Prop%IUB                     = Me%Prop%IUB + 1
        Me%PropIndex%Oxygen             = Me%Prop%IUB
        
        if(Me%ComputeOptions%Nitrogen)then

        !Nitrogen index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Ammonia        = Me%Prop%IUB
        
        
        if (Me%PelagicModel == LifeModel) then
                !Particulate organic carbon index number
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POC            = Me%Prop%IUB
        endif

        !OrganicMatter index number

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%PON            = Me%Prop%IUB
      
                  if(Me%ComputeOptions%Pompools)then
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON1           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON2           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON3           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON4           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%PON5           = Me%Prop%IUB
                       
            end if   
   end if   
      
      !end if

        if(Me%ComputeOptions%Phosphorus)then
            
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phosphate      = Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%POP            = Me%Prop%IUB
            
            if(Me%ComputeOptions%Pompools)then
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP1           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP2           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP3           = Me%Prop%IUB
            
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP4           = Me%Prop%IUB
                
                Me%Prop%IUB                 = Me%Prop%IUB + 1
                Me%PropIndex%POP5           = Me%Prop%IUB
                       
            end if

        end if

        if(Me%ComputeOptions%Silica)then

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%DissolvedSilica= Me%Prop%IUB

            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%BioSilica      = Me%Prop%IUB

        end if


        if(Me%ComputeOptions%Phyto)then
            Me%Prop%IUB                 = Me%Prop%IUB + 1
            Me%PropIndex%Phyto          = Me%Prop%IUB
        endif
        !----------------------------------------------------------------------

    end subroutine PropertyIndexNumber

    !--------------------------------------------------------------------------
    
    
    subroutine ConstructPropertyList

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        type(T_Consumer),      pointer             :: Consumer
        !type(T_Decomposer),    pointer             :: Decomposer
        integer                                    :: Index
        !Local-----------------------------------------------------------------
        
        allocate(Me%PropertyList(Me%Prop%ILB: Me%Prop%IUB))
        
        Index = 0
        
         Me%PropertyList(Me%PropIndex%Oxygen)                = Oxygen_

        if(Me%ComputeOptions%Nitrogen)then
            Me%PropertyList(Me%PropIndex%Ammonia)           = Ammonia_
            Me%PropertyList(Me%PropIndex%PON)               = PON_
            
            if(Me%ComputeOptions%Pompools)then
                Me%PropertyList(Me%PropIndex%PON1)               = PON1_
                Me%PropertyList(Me%PropIndex%PON2)               = PON2_
                Me%PropertyList(Me%PropIndex%PON3)               = PON3_
                Me%PropertyList(Me%PropIndex%PON4)               = PON4_
                Me%PropertyList(Me%PropIndex%PON5)               = PON5_
            end if
        end if

        if(Me%ComputeOptions%Phosphorus)then
            Me%PropertyList(Me%PropIndex%Phosphate)         = Inorganic_Phosphorus_
            Me%PropertyList(Me%PropIndex%POP)               = POP_
            
            if(Me%ComputeOptions%Pompools)then
                Me%PropertyList(Me%PropIndex%POP1)               = POP1_
                Me%PropertyList(Me%PropIndex%POP2)               = POP2_
                Me%PropertyList(Me%PropIndex%POP3)               = POP3_
                Me%PropertyList(Me%PropIndex%POP4)               = POP4_
                Me%PropertyList(Me%PropIndex%POP5)               = POP5_
            end if
        end if


 if(Me%ComputeOptions%Silica)then

            if(Me%PelagicModel .eq. WaterQualityModel) then
                Me%PropertyList(Me%PropIndex%DissolvedSilica)   = DSilica_
            endif

            if(Me%PelagicModel .eq. LifeModel) then            
                Me%PropertyList(Me%PropIndex%DissolvedSilica)   = Silicate_
            endif

            Me%PropertyList(Me%PropIndex%BioSilica)         = BioSilica_
        
        end if

        if(Me%ComputeOptions%Phyto)then
            Me%PropertyList(Me%PropIndex%Phyto)             = Phytoplankton_
        end if

        !Producer index number      
            Producer => Me%FirstProducer
            do while(associated(Producer))
                

               
                Me%PropertyList(Producer%PoolIndex%Carbon)     = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" carbon")
                Me%PropertyList(Producer%PoolIndex%Nitrogen)   = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" nitrogen")
                Me%PropertyList(Producer%PoolIndex%Phosphorus) = &
                            GetPropertyIDNumber(trim(Producer%ID%name)//" phosphorus")


                Producer => Producer%Next
            end do
        
        !Consumer index number
            Consumer => Me%FirstConsumer
            do while(associated(Consumer))
                

                Me%PropertyList(Consumer%PoolIndex%Carbon)     = &
                            GetPropertyIDNumber(trim(Consumer%ID%name)//" carbon")
                Me%PropertyList(Consumer%PoolIndex%Nitrogen)   = &
                            GetPropertyIDNumber(trim(Consumer%ID%name)//" nitrogen")
                Me%PropertyList(Consumer%PoolIndex%Phosphorus) = &
                            GetPropertyIDNumber(trim(Consumer%ID%name)//" phosphorus")


                Consumer => Consumer%Next
            end do
   
        !Decomposer index number
          !  Decomposer => Me%FirstDecomposer
          !  do while(associated(Decomposer))
                
             !   Me%PropertyList(Decomposer%PoolIndex%Carbon)     = &
                            !GetPropertyIDNumber(trim(Decomposer%ID%name)//" carbon")
               ! Me%PropertyList(Decomposer%PoolIndex%Nitrogen)   = &
                            !GetPropertyIDNumber(trim(Decomposer%ID%name)//" nitrogen")
               ! Me%PropertyList(Decomposer%PoolIndex%Phosphorus) = &
                            !GetPropertyIDNumber(trim(Decomposer%ID%name)//" phosphorus")

               ! Decomposer => Decomposer%Next
           ! end do

        !Nitrogen index number



        !----------------------------------------------------------------------

    end subroutine ConstructPropertyList


    !----------------------------------------------------------------------
        subroutine ConstructMatrix

        !Begin-----------------------------------------------------------------

        allocate(Me%Matrix(Me%Size%ILB:Me%Size%IUB,                           &
                             Me%Prop%ILB:Me%Prop%IUB, &
                             Me%Prop%ILB:Me%Prop%IUB))

    
    end subroutine ConstructMatrix
     !----------------------------------------------------------------------

    subroutine ConstructGlobalVariables

        !Arguments-------------------------------------------------------------
        
        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromFile 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromFile = FromFile)
        
        
         call GetData(Me%DT,                                                             &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'DT',                                               &
                     Default      = 3600.,                                              &
                     ClientModule = 'ModuleBenthicEcology',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #1'

        Me%DTDay = Me%DT / (3600. * 24.)
        
 
        call GetData(Me%PelagicModel,                                                   &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PELAGIC_MODEL',                                    &
                     ClientModule = 'ModuleBenthicEcology',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #9'
        if(iflag==0)then
            write(*,*)'Please define the pelagic model to couple with ModuleBenthicEcology'
            stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #9'
        end if

        if((Me%PelagicModel .ne. WaterQualityModel .and. Me%PelagicModel .ne. LifeModel))then
            write(*,*)'Pelagic model to couple with ModuleBenthicEcology must be one of the following:'
            write(*,*)trim(WaterQualityModel)
            write(*,*)trim(LifeModel)
            stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #10'
        endif
                 
        call GetData(Me%ComputeOptions%Nitrogen,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NITROGEN',                                         &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERROR #11'

        call GetData(Me%ComputeOptions%Phosphorus,                                      &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHOSPHORUS',                                       &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR50'


        call GetData(Me%ComputeOptions%Silica,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'SILICA',                                           &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR60'
        
         call GetData(Me%ComputeOptions%Phyto,                                           &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PHYTO',                                            &
                     Default      = .false.,                                            &
                     ClientModule = 'ModuleBenthicEcology',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ConstructGlobalVariables - ModuleBenthicEcology - ERR90'

    end subroutine ConstructGlobalVariables



!----------------------------------------------------------------------

    subroutine ReadOrganicMatterParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%OrganicMatter%NC_Ratio,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'NC_RATIO',                                         &
                     Default      = 0.18,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOrganicMatterParameters - ModuleBenthos - ERR01'

        call GetData(Me%OrganicMatter%PC_Ratio,                                         &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PC_RATIO',                                         &
                     Default      = 0.024,                                              &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOrganicMatterParameters - ModuleBenthos - ERR10'

    end subroutine ReadOrganicMatterParameters

    !--------------------------------------------------------------------------
    
     

    subroutine ReadOxygenParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Oxygen%Minimum,                                                 &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'MIN_OXYGEN',                                       &
                     Default      = 1e-5,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadOxygenParameters - ModuleBenthos - ERR01'


    end subroutine ReadOxygenParameters
    
    !--------------------------------------------------------------------------
   
    subroutine ReadNitrogenParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Nitrogen%PONDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PON_DECAY_RATE',                                   &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadNitrogenParameters - ModuleBenthos - ERR01'
        
        call GetData(Me%Nitrogen%PONDecayTFactor,                                       &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'PON_DECAY_TFACTOR',                                &
                     Default      = 1.02,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadNitrogenParameters - ModuleBenthos - ERR02'

    end subroutine ReadNitrogenParameters

    !--------------------------------------------------------------------------

    subroutine ReadPhosphorusParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Phosphorus%POPDecayRate,                                        &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POP_DECAY_RATE',                                   &
                     Default      = 0.03,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhosphorusParameters - ModuleBenthos - ERR01'

        call GetData(Me%Phosphorus%POPDecayTFactor,                                     &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'POP_DECAY_TFACTOR',                                &
                     Default      = 1.08,                                               &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadPhosphorusParameters - ModuleBenthos - ERR02'


    end subroutine ReadPhosphorusParameters

    !--------------------------------------------------------------------------

    subroutine ReadSilicaParameters
        
        !Local-----------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL

        !Begin-----------------------------------------------------------------

        call GetData(Me%Silica%BioSiDecayRate,                                          &
                     Me%ObjEnterData, iflag,                                            &
                     SearchType   = FromFile,                                           &
                     keyword      = 'BIOSI_DECAY_RATE',                                 &
                     Default      = 0.1,                                                &
                     ClientModule = 'ModuleBenthos',                                    &
                     STAT         = STAT_CALL)
        if(STAT_CALL .NE. SUCCESS_) stop 'ReadSilicaParameters - ModuleBenthos - ERR01'


    end subroutine ReadSilicaParameters
    
    !--------------------------------------------------------------------------
        
    subroutine ConstructProducers

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_producer>',   &
                                        block_end       = '<end_producer>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    call AddProducer                    (NewProducer)

                    call ConstructProducerParameters    (NewProducer)

                    nullify(NewProducer)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructProducers - ModuleBenthicEcology - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructProducers - ModuleBenthicEcology - ERROR #2'
            else cd1
                    stop       'ConstructProducers - ModuleBenthicEcology - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructProducers


    !--------------------------------------------------------------------------
subroutine AddProducer (ObjProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: ObjProducer
        !Local-----------------------------------------------------------------
        type (T_Producer),      pointer           :: PreviousProducer
        type (T_Producer),      pointer           :: NewProducer
        integer, save                             :: NextProducerID = 1

        !Allocates new Producer
        allocate (NewProducer)
        nullify  (NewProducer%Next)

        !Insert new Producer into list and makes current algae point to it
        if (.not. associated(Me%FirstProducer)) then
            Me%FirstProducer            => NewProducer
            ObjProducer                 => NewProducer
        else
            PreviousProducer            => Me%FirstProducer
            ObjProducer                 => Me%FirstProducer%Next

            do while (associated(ObjProducer))
                PreviousProducer        => ObjProducer
                ObjProducer             => ObjProducer%Next
            enddo
            ObjProducer                 => NewProducer
            PreviousProducer%Next       => NewProducer
        endif

        !Attributes ID
        ObjProducer%ID%ID               = NextProducerID

        NextProducerID                  = NextProducerID + 1


    end subroutine AddProducer
    
    !--------------------------------------------------------------------------
    
        subroutine ConstructProducerParameters (NewProducer)

        !Arguments-------------------------------------------------------------
        type (T_Producer),      pointer           :: NewProducer
        

        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                   :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewProducer%ID%Name,                       &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NAME',                     &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #1'

        call GetData(NewProducer%ID%Description,                &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'DESCRIPTION',              &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #2'

        call GetData(NewProducer%RespirationRate,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'RESPIRATION_RATE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #3'

        call GetData(NewProducer%MortalityRate,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'MORTALITY_RATE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #4'
     
        call GetData(NewProducer%Vmax,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'VMAX',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #5'

        call GetData(NewProducer%KN,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'KN',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #6'
        
                call GetData(NewProducer%KP,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'KP',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #7'
        
        call GetData(NewProducer%alpha,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'ALPHA',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #8'
        
          call GetData(NewProducer%NCratio,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NCRATIO',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #9'
        
            call GetData(NewProducer%PCratio,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'PCRATIO',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructProducerParameters - ModuleBenthicEcology - ERROR #10'

    end subroutine ConstructProducerParameters


!_______________________________________________________________________________________________________
!__________________________________________________________________________________ 
!_______________________________________________________________________________________________________
!________________________________________________________________________________

    subroutine ConstructConsumers

        !Arguments-------------------------------------------------------------
        
        !Local-----------------------------------------------------------------
        type (T_Consumer),      pointer           :: NewConsumer
        integer                                   :: ClientNumber, STAT_CALL
        logical                                   :: BlockFound

        !Begin-----------------------------------------------------------------


do1 :   do
            call ExtractBlockFromBuffer(Me%ObjEnterData,                        &
                                        ClientNumber    = ClientNumber,         &
                                        block_begin     = '<begin_consumer>',   &
                                        block_end       = '<end_consumer>',     &
                                        BlockFound      = BlockFound,           &
                                        STAT            = STAT_CALL)
cd1 :       if(STAT_CALL .EQ. SUCCESS_)then
cd2 :           if (BlockFound) then                                                  
                    
                    call AddConsumer                    (NewConsumer)

                    call ConstructConsumerParameters    (NewConsumer, ClientNumber)

                    nullify(NewConsumer)

                else cd2
                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 

                    if (STAT_CALL .NE. SUCCESS_)                                &
                        stop       'ConstructConsumer - ModuleBenthicEcology - ERROR #1'

                    exit do1
                end if cd2

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then cd1
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBuffer. '
                    stop       'ConstructConsumer - ModuleBenthicEcology - ERROR #2'
            else cd1
                    stop       'ConstructConsumer - ModuleBenthicEcology - ERROR #3'
            end if cd1
        end do do1

    end subroutine ConstructConsumers


    !--------------------------------------------------------------------------

!--------------------------------------------------------------------------


    subroutine AddConsumer (ObjConsumer)

        !Arguments-------------------------------------------------------------
        type (T_Consumer),      pointer           :: ObjConsumer
        !Local-----------------------------------------------------------------
        type (T_Consumer),      pointer           :: PreviousConsumer
        type (T_Consumer),      pointer           :: NewConsumer
        integer, save                             :: NextConsumerID = 1

        !Allocates new Consumer
        allocate (NewConsumer)
        nullify  (NewConsumer%Next)

        !Insert new Consumer into list and makes current ?? point to it
        if (.not. associated(Me%FirstConsumer)) then
            Me%FirstConsumer            => NewConsumer
            ObjConsumer                 => NewConsumer
        else
            PreviousConsumer            => Me%FirstConsumer
            ObjConsumer                 => Me%FirstConsumer%Next

            do while (associated(ObjConsumer))
                PreviousConsumer        => ObjConsumer
                ObjConsumer             => ObjConsumer%Next
            enddo
            ObjConsumer                 => NewConsumer
            PreviousConsumer%Next       => NewConsumer
        endif

        !Attributes ID
        ObjConsumer%ID%ID               = NextConsumerID

        NextConsumerID                  = NextConsumerID + 1


    end subroutine AddConsumer
    
    !--------------------------------------------------------------------------


    subroutine ConstructConsumerParameters (NewConsumer, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Consumer),      pointer           :: NewConsumer
        integer                                   :: ClientNumber
        !External--------------------------------------------------------------
        integer                                   :: iflag, STAT_CALL
        
        !Local-----------------------------------------------------------------
        integer                                   :: FromBlock 

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(NewConsumer%ID%Name,                       &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlock,                  &
                     keyword      = 'NAME',                     &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #1'
        
        

        if(.not. CheckPropertyName(trim(NewConsumer%ID%Name)//" carbon"))&
       stop 'ConstructConsumerParameters - ModuleBenthicEcology - ERROR #2'
 
        call GetData(NewConsumer%O2HX,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'O2HX',                &   ! mgO2/l
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #2'

       call GetData(NewConsumer%O2QX,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'O2QX',                &   ! mgO2/l
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #3'

       call GetData(NewConsumer%Td,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'TD',                &   ! days
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #4'

       call GetData(NewConsumer%BETA20,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'BETA20',                &   ! 
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #5'
      
       call GetData(NewConsumer%FILT20,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'FILT20',                &   
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #6'
        
      
       
       call GetData(NewConsumer%RESP20,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'RESP20',                &   ! 
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #7'
        
             call GetData(NewConsumer%TemperatureFactor,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'TFAC',                &   ! 
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #8' 
        
        
                     call GetData(NewConsumer%NCRatio,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'NCRATIO',                &   ! 
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #9'  


                     call GetData(NewConsumer%PCRatio,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'PCRATIO',                &   ! 
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #10' 
        call ConstructGrazing               (NewConsumer%Grazing , ClientNumber)
        
    
    
    
 end subroutine ConstructConsumerParameters
!--------------------------------------------------------------------------




subroutine ConstructGrazing (Grazing, ClientNumber)

        !Arguments-------------------------------------------------------------
        type (T_Grazing)                                :: Grazing
        integer                                         :: ClientNumber

        !External--------------------------------------------------------------
        integer                                         :: iflag, STAT_CALL
        logical                                         :: BlockInBlockFound
        
        !Local-----------------------------------------------------------------
        integer                                         :: FromBlock 
        integer                                         :: FirstLine, LastLine
        type (T_Prey), pointer                          :: NewPrey

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlock = FromBlock)

        call GetData(Grazing%Vmax,                                  &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'VMAX',             &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #1'


        call GetData(Grazing%Ks,                                    &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'KS',                   &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #2'


        call GetData(Grazing%Ass_Efic,                              &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlock,                      &
                     keyword      = 'ASS_EFIC',                &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #3'
        
      

        do
            call ExtractBlockFromBlock(Me%ObjEnterData, ClientNumber,         &
                                       '<begin_food>', '<end_food>',          &
                                       BlockInBlockFound,                     &
                                       FirstLine = FirstLine,                 &
                                       LastLine  = LastLine,                  &
                                       STAT      = STAT_CALL)
            if      (STAT_CALL .EQ. SUCCESS_) then    
            
                if (BlockInBlockFound) then
            
                    if     (((LastLine + 1) - (FirstLine - 1)) .GE. 1) then

                        call AddPrey        (Grazing, NewPrey)

                        call ConstructPrey  (NewPrey)

                        nullify(NewPrey)

                    else
                        write(*,*)  
                        write(*,*) 'Error counting preys. '
                        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #4'
                    end if
                else

                    call Block_Unlock(Me%ObjEnterData, ClientNumber, STAT = STAT_CALL) 
                    if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #5'
                    
                    exit 

                end if 

            else if (STAT_CALL .EQ. BLOCK_END_ERR_) then
                write(*,*)  
                write(*,*) 'Error calling ExtractBlockFromBlock. '
                if (STAT_CALL .NE. SUCCESS_) stop 'ConstructGrazing - ModuleBenthicEcology - ERROR #6'
            end if

        end do
        
    end subroutine ConstructGrazing
    
    
    ! -----------------------------------------------------------------------------------------
    subroutine AddPrey (Grazing, Prey)

        !Arguments-------------------------------------------------------------
        type (T_Grazing)                          :: Grazing
        type (T_Prey),          pointer           :: Prey

        !Local-----------------------------------------------------------------
        type (T_Prey),          pointer           :: PreviousPrey
        type (T_Prey),          pointer           :: NewPrey
        integer, save                             :: NextPreyID = 1

        !Allocates new Producer
        allocate (NewPrey)
        nullify  (NewPrey%Next)

        !Insert new Prey into list and makes current algae point to it
        if (.not. associated(Grazing%FirstPrey)) then
            NextPreyID              = 1
            Grazing%FirstPrey       => NewPrey
            Prey                    => NewPrey
        else
            PreviousPrey            => Grazing%FirstPrey
            Prey                    => Grazing%FirstPrey%Next

            do while (associated(Prey))
                PreviousPrey        => Prey
                Prey                => Prey%Next
            enddo
            Prey                    => NewPrey
            PreviousPrey%Next       => NewPrey
        endif

        !Attributes ID
        Prey%ID%ID               = NextPreyID

        NextPreyID               = NextPreyID + 1


    end subroutine AddPrey
!--------------------------------------------------------------------------------
    !--------------------------------------------------------------------------


    subroutine ConstructPrey (NewPrey)

        !Arguments-------------------------------------------------------------
        type (T_Prey), pointer                          :: NewPrey

        !External--------------------------------------------------------------
        integer                                         :: STAT_CALL, iflag
        integer                                         :: FromBlockInBlock      

        !Local-----------------------------------------------------------------
        type(T_Producer), pointer                       :: Producer

        !Begin-----------------------------------------------------------------

        call GetExtractType    (FromBlockInBlock = FromBlockInBlock)

        call GetData(NewPrey%ID%Name,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'NAME',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #1'
        
                call GetData(NewPrey%NCRatio,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'NCRATIO',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #2'
                 
                 call GetData(NewPrey%PCRatio,                               &
                     Me%ObjEnterData, iflag,                        &
                     SearchType   = FromBlockInBlock,               &
                     keyword      = 'PCRATIO',                         &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,      &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #3'
        
        call GetData(NewPrey%Use_Carbon,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'CARBON_USE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #4'
        
        call GetData(NewPrey%Use_Nitrogen,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'NITROGEN_USE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #5'
        
                call GetData(NewPrey%Use_Phosphorus,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'PHOSPHORUS_USE',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #6'
        
        call GetData(NewPrey%ParticulateWaterFood,                    &
                     Me%ObjEnterData, iflag,                    &
                     SearchType   = FromBlockInBlock,                  &
                     keyword      = 'PARTICWATERFOOD',               &
                     ClientModule = MohidModules(mBenthicEcology_)%Name,  &
                     STAT         = STAT_CALL)
        if (STAT_CALL .NE. SUCCESS_) stop 'ConstructPrey - ModuleBenthicEcology - ERROR #7'
        
        if((Me%PelagicModel == LifeModel))then
        if(.not. CheckPropertyName(trim(NewPrey%ID%Name)//" nitrogen"))&
        stop 'ConstructPrey - ModuleBenthicEcology - ERROR #7'
        endif        

        Producer => Me%FirstProducer
        do while(associated(Producer))

            !if(NewPrey%ID%Name == trim(Producer%ID%Name))then

             !   NewPrey%Use_Chl = .true.

              !  if(Producer%Use_Silica)NewPrey%Use_Silica = .true.

            !end if

         

            Producer => Producer%Next
        end do

    end subroutine ConstructPrey
  
   
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SELECTOR SE

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine GetDTBenthicEcology(BenthicEcology_ID, DTDay, DTSecond, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcology_ID
        real,    optional, intent(OUT)      :: DTDay
        real,    optional, intent(OUT)      :: DTSecond
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcology_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(DTDay   )) DTDay    = Me%DTDay
            if (present(DTSecond)) DTSecond = Me%DT

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetDTBenthicEcology
    
    !--------------------------------------------------------------------------

    
    subroutine GetBenthicEcologyPropertyList(Life_ID, PropertyList, STAT)

        !Arguments-------------------------------------------------------------
        integer                                                 :: Life_ID
        integer, dimension(:), pointer                          :: PropertyList
        integer, optional, intent(OUT)                          :: STAT

        !External--------------------------------------------------------------
        integer                                                 :: ready_              

        !Local-----------------------------------------------------------------
        integer                                                 :: STAT_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(Life_ID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBenthicEcology_, Me%InstanceID)

            PropertyList => Me%PropertyList

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthicEcologyPropertyList
    
    
        !--------------------------------------------------------------------------

    subroutine GetBenthicEcologySize(BenthicEcology_ID, PropLB, PropUB, STAT)

        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcology_ID
        integer, optional, intent(OUT)      :: PropLB,PropUB
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
       
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcology_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            if (present(PropLB   )) PropLB    = Me%Prop%ILB
            if (present(PropUB   )) PropUB    = Me%Prop%IUB

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthicEcologySize
    
    !--------------------------------------------------------------------------
    
     subroutine GetBenthicEcologyPropIndex (BenthicEcology_ID, PropertyIDNumber, PropertyIndex, STAT)

                                     

        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcology_ID
        integer,           intent(IN )      :: PropertyIDNumber
        integer,           intent(OUT)      :: PropertyIndex
        integer, optional, intent(OUT)      :: STAT

        !External--------------------------------------------------------------
        integer                             :: ready_              

        !Local-----------------------------------------------------------------
        integer                             :: STAT_, CurrentIndex
        logical                             :: found
               
        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcology_ID, ready_)    
        
cd1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then


            found = .false.
            do CurrentIndex = Me%Prop%ILB,Me%Prop%IUB

                if (PropertyIDNumber.eq. Me%PropertyList(CurrentIndex))then
                    PropertyIndex = CurrentIndex
                    found = .true.
                    exit
                end if
            
            end do

            if(.not. found)then
                STAT_ = NOT_FOUND_ERR_
            else
                STAT_ = SUCCESS_
            endif

        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT))STAT = STAT_

    end subroutine GetBenthicEcologyPropIndex


    !--------------------------------------------------------------------------

   
                
                
    integer function SearchPropIndex (PropIDNumber)

        !Arguments-------------------------------------------------------------
        integer,  intent(IN )                         :: PropIDNumber
        !Local-----------------------------------------------------------------
        integer                                       :: CurrentIndex

        !----------------------------------------------------------------------

        SearchPropIndex = UNKNOWN_

        do CurrentIndex = Me%Prop%ILB, Me%Prop%IUB

            if (PropIDNumber == Me%PropertyList(CurrentIndex))then
                SearchPropIndex = CurrentIndex
                exit
            end if
                    
        end do

    end function SearchPropIndex

    !--------------------------------------------------------------------------

        subroutine UnGetBenthicEcology(BenthicEcologyID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BenthicEcologyID
        integer, dimension(:), pointer                  :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcologyID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBenthicEcology_, Me%InstanceID, "UnGetBenthicEcology3D_I")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBenthicEcology
    
     !--------------------------------------------------------------------------
        subroutine GetBenthicEcologyRateFlux(BenthicEcologyID, FirstProp, SecondProp, RateFlux, STAT)


        !Arguments-------------------------------------------------------------
        integer                             :: BenthicEcologyID
        integer,           intent(IN )      :: FirstProp, SecondProp
        real,    dimension(:), pointer      :: RateFlux
        integer, optional, intent(OUT)      :: STAT

        !Local-----------------------------------------------------------------
        integer                             :: STAT_
        integer                             :: ready_              

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcologyID, ready_)    
        
if1 :   if ((ready_ .EQ. IDLE_ERR_     ) .OR.                                 &
            (ready_ .EQ. READ_LOCK_ERR_)) then

            call Read_Lock(mBenthicEcology_, Me%InstanceID)

            select case(FirstProp)
                

                case default

                    RateFlux => Me%Matrix    (:, FirstProp, SecondProp      )

            end select

            STAT_ = SUCCESS_
        else 
            STAT_ = ready_
        end if if1

        if (present(STAT)) STAT = STAT_

    end subroutine GetBenthicEcologyRateFlux

    !--------------------------------------------------------------------------
    
      subroutine UnGetBenthicEcologyRateFlux(BenthicEcologyID, Array, STAT)

        !Arguments-------------------------------------------------------------
        integer                                         :: BenthicEcologyID
        real, dimension(:), pointer                     :: Array
        integer, intent(OUT), optional                  :: STAT

        !Local-----------------------------------------------------------------
        integer                                         :: STAT_, ready_

        !----------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(BenthicEcologyID, ready_)

        if (ready_ .EQ. READ_LOCK_ERR_) then

            nullify(Array)
            call Read_Unlock(mBenthicEcology_, Me%InstanceID, "UnGetBenthicEcologyRateFlux")

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_

    end subroutine UnGetBenthicEcologyRateFlux
    
    
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    subroutine ModifyBenthicEcology(ObjBenthicEcologyID, Temperature, WaterVolume, &
               CellArea, MassInKgFromWater, Sediment, OpenPoints, Mass, STAT)
    
    
            !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthicEcologyID
        real,    dimension(:  ), pointer            :: Temperature
        real,    dimension(: ,: ), pointer          :: MassInKgFromWater
        real,    dimension(:  ), pointer            :: Sediment
        real,    dimension(:  ), pointer            :: WaterVolume
        real,    dimension(:  ), pointer            :: CellArea
        integer, dimension(:  ), pointer, optional  :: OpenPoints
        real,    dimension(:,:), pointer            :: Mass
        integer, intent(OUT), optional              :: STAT

        !Local-----------------------------------------------------------------
        integer                                     :: STAT_, ready_
        integer                                     :: Index
        logical                                     :: Compute
        !Begin-----------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBenthicEcologyID, ready_)

        if (ready_ .EQ. IDLE_ERR_) then

            Me%ExternalVar%Temperature  => Temperature
            if (.not. associated(Me%ExternalVar%Temperature))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR01'
                
             
             Me%ExternalVar%Sediment  => Sediment
            if (.not. associated(Me%ExternalVar%Sediment))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR03'
                
             Me%ExternalVar%MassInKgFromWater  => MassInKgFromWater
            if (.not. associated(Me%ExternalVar%MassInKgFromWater))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR04'
                
            Me%ExternalVar%WaterVolume  => WaterVolume
            if (.not. associated(Me%ExternalVar%WaterVolume))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR05'
                
            Me%ExternalVar%CellArea  => CellArea
            if (.not. associated(Me%ExternalVar%CellArea))       &
                stop 'ModifyBenthicEcology - ModuleBenthicEcology - ERR06'
              
            Me%ExternalVar%Mass         => Mass
            if (.not. associated(Me%ExternalVar%Mass))              &
                stop 'ModifyBenthos - ModuleBenthos - ERR05'


            do Index = Me%Size%ILB, Me%Size%IUB

                if (present(OpenPoints)) then
                    if (OpenPoints(Index) == OpenPoint) then
                        Compute = .true.
                    else
                        Compute = .false.
                    endif
                else
                    Compute = .true.
                endif

                Me%Matrix(Index,:,:) = 0.
                
                
                                if(Compute)then
                    
                    call ComputeBenthicProducers      (Index)
                    
                     call ComputeBenthicConsumers    (Index)
                
                    if(Me%ComputeOptions%Nitrogen  ) call ComputeBenthicNitrogen    (Index)
                    
                    if(Me%ComputeOptions%Phosphorus) call ComputeBenthicPhosphorus  (Index)
                    
                    if(Me%ComputeOptions%Silica    ) call ComputeBenthicSilica      (Index)
                
                endif

            enddo

            STAT_ = SUCCESS_
        else               
            STAT_ = ready_
        end if

        if (present(STAT)) STAT = STAT_
    
    end subroutine ModifyBenthicEcology
    
    !--------------------------------------------------------------------------------------
    
    
    subroutine ComputeBenthicProducers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Producer),      pointer             :: Producer
        
       
        integer :: Producer_N, Producer_C, Producer_P
        integer :: AM,      NA, IP, O2
        integer :: PON, POP, POC
        real    :: AverageRadiation, TemperatureDependence
        real    :: Lightlim, NLim, PLim, NutLim, GrowthRate,UptakeNA, UptakeAM, UptakeP
        real    :: RespirationC, RespirationN, RespirationP
        real    :: MortalityC, MortalityN, MortalityP
        real    :: x1,x2,x3,x4, AmmoniaPreferenceFactor

    !------------------------------------------------------------------------
       
        AM      = Me%PropIndex%Ammonia
        NA      = Me%PropIndex%Nitrate
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        IP      = Me%PropIndex%Phosphate
        O2      = Me%PropIndex%Oxygen
        
        if(Me%PelagicModel == LifeModel) then
        POC     = Me%PropIndex%POC
        endif


        Producer => Me%FirstProducer

d1:     do while(associated(Producer))
        
        AverageRadiation = Me%ExternalVar%ShortWaveAverage(index) !W/m2
        
        ! TemperatureDependence is the function expressing the dependence on the temperature
        ! Dimensionless factor
        TemperatureDependence=exp(0.07*(Me%ExternalVar%Temperature(Index))) 
        
        Producer_N   = Producer%PoolIndex%Nitrogen
        Producer_C   = Producer%PoolIndex%Carbon
        Producer_P   = Producer%PoolIndex%Phosphorus
        
        ! Evans and Parslow model (1985) model 
        Lightlim =  Producer%alpha* AverageRadiation / &
                    sqrt(Producer%Vmax*Producer%Vmax + (Producer%alpha**2) * ( AverageRadiation**2))
        
        ! Nutrients Limitation
        NLim =   (Me%ExternalVar%MassInKgFromWater(AM,Index)/Me%ExternalVar%WaterVolume(Index)                  + &
                 Me%ExternalVar%MassInKgFromWater(NA,Index)/Me%ExternalVar%WaterVolume(Index) )                 / &
                (Producer%KN + Me%ExternalVar%MassInKgFromWater(AM, Index)/Me%ExternalVar%WaterVolume(index)    + &
                 Me%ExternalVar%MassInKgFromWater(NA, Index)/Me%ExternalVar%WaterVolume(index))
        
        Plim = (Me%ExternalVar%MassInKgFromWater(IP,Index)/Me%ExternalVar%WaterVolume(Index))                 /   &
               (Producer%KP + Me%ExternalVar%MassInKgFromWater(IP, Index)/Me%ExternalVar%WaterVolume(Index)    )
        
        NutLim = min(NLim,PLim)
        
        ! ammonia preference factor (AmmoniaPreferenceFactor,dimensionless)
        ! following the same calculations made in module WaterQuality for phytoplankton
        
         x1 = (Me%ExternalVar%MassInKgFromWater(AM, index) * Me%ExternalVar%MassInKgFromWater(NA, index))        / &
                                                             Me%ExternalVar%WaterVolume(Index)

         x2 = (Producer%KN + Me%ExternalVar%MassInKgFromWater(AM, index)/Me%ExternalVar%WaterVolume(Index))                                   &
               * (Producer%KN + Me%ExternalVar%MassInKgFromWater(NA, index)/Me%ExternalVar%WaterVolume(Index)) 

         x3 = Producer%KN * Me%ExternalVar%MassInKgFromWater(AM, index)/Me%ExternalVar%WaterVolume(Index)

         x4 = (Me%ExternalVar%MassInKgFromWater(AM, index)/Me%ExternalVar%WaterVolume(Index) + &
               Me%ExternalVar%MassInKgFromWater(NA, index)/Me%ExternalVar%WaterVolume(Index)) * &
               (Producer%KN + Me%ExternalVar%MassInKgFromWater(NA, index)/Me%ExternalVar%WaterVolume(Index))

            if ((x1 .EQ. 0.0) .AND. (x3 .EQ. 0.0)) then
                AmmoniaPreferenceFactor = 0.0                 
            else 
                AmmoniaPreferenceFactor = (x1 / x2) + (x3 / x4)
            end if 
        
        ! growth rate 1/day
        GrowthRate=Producer%Vmax*LightLim*NutLim
        
        ! uptake of ammonia (Michaelis Menten kinetics)
        ! KgN     = 1/day  * (KgN/KgC) * KgC*day
        UptakeAM = AmmoniaPreferenceFactor * GrowthRate * Producer%NCRatio * &
                  Me%ExternalVar%Mass(Producer_C, Index)* Me%DTDay ! KgN
        
        UptakeNA =(1.- AmmoniaPreferenceFactor)* GrowthRate * Producer%NCRatio * &
                  Me%ExternalVar%Mass(Producer_C, Index)* Me%DTDay ! KgN
        
        UptakeP = GrowthRate * Producer%PCRatio * Me%ExternalVar%Mass(Producer_C, Index)* Me%DTDay ! KgP
        
        
        !KgC = KgC *day * 1/day * [-]                                                                                       
        MortalityC = Me%ExternalVar%Mass(Producer_C, Index) * Me%DTDay * Producer%MortalityRate*TemperatureDependence 
        MortalityN = MortalityC*Producer%NCRatio  ! KgN
        MortalityP = MortalityC*Producer%PCRatio  ! KgP
        
        ! units of mass
        RespirationC = Me%ExternalVar%Mass(Producer_C, Index) * Producer%RespirationRate*TemperatureDependence * Me%DTDay !KgC
        RespirationN = RespirationC*Producer%NCRatio  ! KgN
        RespirationP = RespirationC*Producer%PCRatio  ! KgP
        
                  
        ! mass balance
        
        ! KgC
        Me%ExternalVar%Mass(Producer_C,     Index) = Me%ExternalVar%Mass(Producer_C,     Index) + &  ! KgC
                                                     GrowthRate                                 * &  ! 1/day
                                                     Me%ExternalVar%Mass(Producer_C, Index)     * &  ! KgC
                                                     Me%DTDay                                   - &  ! day
                                                     MortalityC                                 - &  ! KgC
                                                     RespirationC                                    ! KgC
        ! KgN
        Me%ExternalVar%Mass(Producer_N,     Index) = Me%ExternalVar%Mass(Producer_N,     Index) + &
                                                     GrowthRate                                 * &
                                                     Me%ExternalVar%Mass(Producer_C, Index)     * &
                                                     Producer%NCRatio                           - &
                                                     MortalityN                                 - &
                                                     RespirationN  
        ! KgP
        Me%ExternalVar%Mass(Producer_P,     Index) = Me%ExternalVar%Mass(Producer_P,     Index) + &  ! KgP
                                                     GrowthRate                                 * &  ! 
                                                     Me%ExternalVar%Mass(Producer_C, Index)     * &
                                                     Producer%PCRatio                           * &
                                                     Me%DTDay                                   - &
                                                     MortalityP                                 - &
                                                     RespirationP  
           
        if(Me%ComputeOptions%Nitrogen)then
            
            !what passes from The benthic producer to PON
            Me%ExternalVar%Mass(PON,     Index) = Me%ExternalVar%Mass(PON,     Index) + MortalityN
            !what passes from The benthic producer to Ammonia
            Me%ExternalVar%MassInKgFromWater(AM,     Index) = Me%ExternalVar%MassInKgFromWater(AM,     Index) + RespirationN
            !what passes from nitrate to the benthic producer
            Me%ExternalVar%MassInKgFromWater(NA,     Index) = Me%ExternalVar%MassInKgFromWater(NA,     Index) - UptakeNA
            !what passes from ammonia to benthic producer 
            Me%ExternalVar%MassInKgFromWater(AM,     Index) = Me%ExternalVar%MassInKgFromWater(AM,     Index) - UptakeAM
        end if
           
       if(Me%ComputeOptions%Phosphorus)then
            
            !what passes from The benthic producer to POP
            Me%ExternalVar%Mass(POP,     Index) = Me%ExternalVar%Mass(POP,     Index) + MortalityP
            !what passes from The benthic producer to Phosphate
            Me%ExternalVar%MassInKgFromWater(IP,     Index) = Me%ExternalVar%MassInKgFromWater(IP,     Index) + RespirationP
            !what passes from Phosphate to benthic producer 
            Me%ExternalVar%MassInKgFromWater(IP,     Index) = Me%ExternalVar%MassInKgFromWater(IP,     Index) - UptakeP
       end if
       
       Me%ExternalVar%MassInKgFromWater(O2, Index)=Me%ExternalVar%MassInKgFromWater(O2, Index) + &
                                                   (GrowthRate                                 * &
                                                    Me%ExternalVar%Mass(Producer_C, Index)     - &
                                                    RespirationC)*32./12.

        if(Me%PelagicModel == LifeModel) then
        !what passes from The benthic producer to POC (only if ModuleLife is active)
        Me%ExternalVar%Mass(POC,     Index) = Me%ExternalVar%Mass(POC,     Index) + MortalityC
        endif


            Producer => Producer%Next
        end do d1




        
        end subroutine ComputeBenthicProducers
    
    
    
    
    
    
       subroutine ComputeBenthicConsumers (index)
    
    !Arguments---------------------------------------------------------------

        integer, intent(IN) :: index

    !Local-------------------------------------------------------------------
        type(T_Consumer),      pointer             :: Consumer
        type (T_Grazing)                           :: Grazing
        type(T_Prey),          pointer             :: Prey
        integer                                    :: PON
        integer                                    :: POP
        integer                                    :: POC
        integer                                    :: Phyto
        integer                                    :: AM
        integer                                    :: IP
        integer                                    :: O2
        integer                                    :: Consumer_N
        integer                                    :: Consumer_C
        integer                                    :: Consumer_P
        integer                                    :: PreyIndexC
        integer                                    :: PreyIndexN
        integer                                    :: PreyIndexP
        real                                       :: TemperatureDependence
        real                                       :: Mortality
        real                                       :: MortalityRate
        real                                       :: Respiration
        real                                       :: OxygenLimitation
        real                                       :: SedimentLimitation
        real                                       :: Predation
        real                                       :: PredationRate
        real                                       :: RespirationRate
        real                                       :: IngestionC
        real                                       :: IngestionN
        real                                       :: IngestionP
        real                                       :: IngestionC_tot
        real                                       :: IngestionN_tot
        real                                       :: IngestionP_tot
        real                                       :: EgestionC
        real                                       :: EgestionN
        real                                       :: EgestionP
        real                                       :: EgestionC_tot
        real                                       :: EgestionN_tot
        real                                       :: EgestionP_tot
        real                                       :: FiltrationRate
        real, parameter                            :: atss = 81.
        real, parameter                            :: btsss = 24.
        
        
        

        
        AM      = Me%PropIndex%Ammonia
        PON     = Me%PropIndex%PON
        POP     = Me%PropIndex%POP
        Phyto   = Me%PropIndex%Phyto
        IP      = Me%PropIndex%Phosphate
        O2      = Me%PropIndex%Oxygen
        
        if(Me%PelagicModel == LifeModel) then
        POC     = Me%PropIndex%POC
        endif
        

        
        Consumer => Me%FirstConsumer

d1:     do while(associated(Consumer))
        
        
        Consumer_N   = Consumer%PoolIndex%Nitrogen
        Consumer_C   = Consumer%PoolIndex%Carbon
        Consumer_P   = Consumer%PoolIndex%Phosphorus 
        
        TemperatureDependence =Consumer%TemperatureFactor**(Me%ExternalVar%Temperature(Index)-20.)
       
         ! cohesive sediment is  g/l
        SedimentLimitation=(1.-0.01*(atss+btsss*log10(Me%ExternalVar%Sediment(Index))))
        ! O2 enters the equation as concentration in g/l (or Kg/m3)
        OxygenLimitation = 1./(1.+exp(1.1*(Consumer%O2HX-Me%ExternalVar%MassInKgFromWater(O2,Index)/ &
                             Me%ExternalVar%WaterVolume(index))/(Consumer%O2HX-Consumer%O2QX))) ![dimensionless]
        
        MortalityRate=(1.-OxygenLimitation)*(-log(1./100.)/Consumer%Td) ! 1/day
        ! 1/KgC*day 
        PredationRate = Consumer%BETA20*TemperatureDependence
        
        !
        ! 
        ! m3/day/KgC =      m3/day /KgC    *  [-]                                                                                                                
        FiltrationRate=Consumer%FILT20*  TemperatureDependence* &
                                          OxygenLimitation     * &
                                          SedimentLimitation
        
        !  1/day =           l/day     *[-]                                                                     
        RespirationRate=Consumer%RESP20*TemperatureDependence* &
                                         OxygenLimitation 
        
       
      !     KgC   =          1/day      * KgC                                 *day
        Respiration=RespirationRate*Me%ExternalVar%Mass(Consumer_C, Index)* Me%DTDay  
        
        !     KgC   =          1/day      * KgC                                 *day
        Mortality =  MortalityRate*Me%ExternalVar%Mass(Consumer_C, Index) * Me%DTDay  ! KgC
        !     KgC   =     1/KgC*day      * KgC*KgC                                 *day
        Predation = PredationRate*(Me%ExternalVar%Mass(Consumer_C, Index)**2) * Me%DTDay    
       

       IngestionC_tot=0.
       IngestionN_tot=0.
       IngestionP_tot=0.
       EgestionC_tot=0.
       EgestionN_tot=0.
       EgestionP_tot=0.
       
       
       Prey => Consumer%Grazing%FirstPrey
  
  
       d2:   do while(associated(Prey))


          if ((Prey%ID%Name=='phytoplankton').AND.(Me%PelagicModel == WaterQualityModel)) then
               
               PreyIndexC  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)))
          
          else if ((Prey%ID%Name=='particulate organic').AND.(Me%PelagicModel == WaterQualityModel)) then
               
               if(Prey%Use_Nitrogen) then
               PreyIndexN  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" nitrogen"))
               endif
               if(Prey%Use_Phosphorus) then
               PreyIndexP  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" phosphorus"))
               endif
         else     
               
             if(Prey%Use_Carbon)  then
             
                  PreyIndexC  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" carbon"))
             endif
                  
             if(Prey%Use_Nitrogen)  then
                  PreyIndexN  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" nitrogen"))
             endif     
             
             if(Prey%Use_Phosphorus)  then
                  PreyIndexP  = SearchPropIndex(GetPropertyIDNumber(trim(Prey%ID%Name)//" phosphorus"))
             endif
         
         endif
          
           
        
                         
        if(Prey%Use_Nitrogen)  then
            if((Prey%ID%Name=='particulate organic') .and. (Prey%ParticulateWaterFood)) then
                ! the food eaten is from the water. 
                Me%Matrix(Index, PreyIndexN, Consumer_N) =(FiltrationRate/Prey%NCRatio)*Me%ExternalVar%Mass(Consumer_N, Index)* &
                                                          (Me%ExternalVar%MassInKgFromWater(PreyIndexN, Index)                / &
                                                           Me%ExternalVar%WaterVolume(index))* Me%DTDay ! kgN
                                                           
                Me%ExternalVar%MassInKgFromWater(PreyIndexN, Index)= Me%ExternalVar%MassInKgFromWater(PreyIndexN, Index)      - &
                                                                     Me%Matrix(Index, PreyIndexN, Consumer_N)
                else
                ! the food eaten is on the bottom
                Me%Matrix(Index, PreyIndexN, Consumer_N) =(FiltrationRate/Prey%NCRatio)*Me%ExternalVar%Mass(Consumer_N, Index) * &
                                                          (Me%ExternalVar%Mass(PreyIndexN, Index)                              / &
                                                           Me%ExternalVar%WaterVolume(index))* Me%DTDay ! kgN
                
                Me%ExternalVar%Mass(PreyIndexN, Index)   = Me%ExternalVar%Mass(PreyIndexN, Index)- &
                                                           Me%Matrix(Index, PreyIndexN, Consumer_N)  ! kgN
            endif
        
        IngestionN=Consumer%Grazing%Ass_Efic*Me%Matrix(Index, PreyIndexN, Consumer_N)  ! kgN
        EgestionN=(1.-Consumer%Grazing%Ass_Efic)*Me%Matrix(Index, PreyIndexN, Consumer_N) ! ! kgN
        endif
        
        if(Prey%Use_Carbon)  then
        
            if(Prey%ParticulateWaterFood)then
                ! the food eaten is from the water. 
                Me%Matrix(Index, PreyIndexC, Consumer_C) = FiltrationRate*Me%ExternalVar%Mass(Consumer_C, Index)  * &
                                                          (Me%ExternalVar%MassInKgFromWater(PreyIndexC, Index)    / &
                                                           Me%ExternalVar%WaterVolume(index))* Me%DTDay ! kgN
                
                Me%ExternalVar%MassInKgFromWater(PreyIndexC, Index)= Me%ExternalVar%MassInKgFromWater(PreyIndexC, Index)- &
                                                                    Me%Matrix(Index, PreyIndexC, Consumer_C)
                else
                ! the food eaten is on the bottom
                Me%Matrix(Index, PreyIndexC, Consumer_C) =FiltrationRate*Me%ExternalVar%Mass(Consumer_C, Index)         * &
                                                         (Me%ExternalVar%Mass(PreyIndexC, Index)                        / &
                                                         Me%ExternalVar%WaterVolume(index))* Me%DTDay 
                
                Me%ExternalVar%Mass(PreyIndexC, Index)= Me%ExternalVar%Mass(PreyIndexC, Index)                          - &
                                                        Me%Matrix(Index, PreyIndexC, Consumer_C)
            endif
        
        IngestionC=Consumer%Grazing%Ass_Efic*Me%Matrix(Index, PreyIndexC, Consumer_C)
        EgestionC=(1.-Consumer%Grazing%Ass_Efic)*Me%Matrix(Index, PreyIndexC, Consumer_C)
        
        endif
        
        if(Prey%Use_Phosphorus)  then
        
            if((Prey%ID%Name=='particulate organic') .and. (Prey%ParticulateWaterFood)) then
                ! food eaten from the water
                Me%Matrix(Index, PreyIndexP, Consumer_P) =(FiltrationRate/Prey%PCRatio)                         * &
                                                          Me%ExternalVar%Mass(Consumer_P, Index)                * &
                                                          (Me%ExternalVar%MassInKgFromWater(PreyIndexP, Index)  / &
                                                           Me%ExternalVar%WaterVolume(index))* Me%DTDay ! kgP
                
                Me%ExternalVar%MassInKgFromWater(PreyIndexP, Index)= Me%ExternalVar%MassInKgFromWater(PreyIndexP, Index)- &
                                                                     Me%Matrix(Index, PreyIndexP, Consumer_P)
                else
                ! food eaten from the bottom 
                Me%Matrix(Index, PreyIndexP, Consumer_P) =(FiltrationRate/Prey%PCRatio)                         * &
                                                           Me%ExternalVar%Mass(Consumer_P, Index)*(Me%ExternalVar%Mass(PreyIndexP, Index)/ &
                                                           Me%ExternalVar%WaterVolume(index))* Me%DTDay
                
                Me%ExternalVar%Mass(PreyIndexP, Index)= Me%ExternalVar%Mass(PreyIndexP, Index)- Me%Matrix(Index, PreyIndexP, Consumer_P)
            
            endif
            
        IngestionP=Consumer%Grazing%Ass_Efic*Me%Matrix(Index, PreyIndexP, Consumer_P)
        EgestionP=(1.-Consumer%Grazing%Ass_Efic)*Me%Matrix(Index, PreyIndexP, Consumer_P)
        endif
       
       if ((Prey%ID%Name=='phytoplankton').AND.(Me%PelagicModel == WaterQualityModel)) then
           IngestionN=IngestionC*Consumer%NCRatio 
           IngestionP=IngestionC*Consumer%PCRatio
           EgestionN=EgestionC*Consumer%NCRatio
           EgestionP=EgestionC*Consumer%PCRatio
       endif
       
       
       if ((Prey%ID%Name=='particulate organic').AND.(Me%PelagicModel == WaterQualityModel)) then
       
            if(Prey%Use_Nitrogen) then
                 ! particulate organic nitrogen
               IngestionC=IngestionN/Consumer%NCRatio
               EgestionC=EgestionN/Consumer%NCRatio
               IngestionP=IngestionN*(Consumer%PCRatio/Consumer%NCRatio)
               EgestionP=EgestionN*(Consumer%PCRatio/Consumer%NCRatio)
           endif
           
           if(Prey%Use_Phosphorus) then
               ! particulate organic phosphorus
               IngestionC=IngestionP/Prey%PCRatio
               EgestionC=EgestionP/Prey%PCRatio
               
               IngestionN=IngestionP*(Prey%NCRatio/Prey%PCRatio)
               EgestionN=EgestionP*(Prey%NCRatio/Prey%PCRatio)
           endif
       
           IngestionC_tot=IngestionC_tot+IngestionC   ! Ingestion of each prey is added to the total. 
           IngestionN_tot=IngestionN_tot+IngestionN
           IngestionP_tot=IngestionP_tot+IngestionP
           
           
           EgestionC_tot=EgestionC_tot+EgestionC      ! egestion of each prey is added to the total.
           EgestionN_tot=EgestionN_tot+EgestionN
           EgestionP_tot=EgestionP_tot+EgestionP
           
       endif
       
       Prey => Prey%Next
        end do d2
            
        Me%ExternalVar%Mass(Consumer_C, Index) = Me%ExternalVar%Mass(Consumer_C, Index)+IngestionC_tot &
                                                                                       -Mortality &
                                                                                       -Predation &
                                                                                       -Respiration

        Me%ExternalVar%Mass(Consumer_N, Index) = Me%ExternalVar%Mass(Consumer_N, Index)+IngestionN_tot &
                                                                                       -Mortality*Me%OrganicMatter%NC_Ratio &
                                                                                       -Predation*Me%OrganicMatter%NC_Ratio &
                                                                                       -Respiration*Me%OrganicMatter%NC_Ratio
                                                                                       
       Me%ExternalVar%Mass(Consumer_P, Index) = Me%ExternalVar%Mass(Consumer_P, Index)+ IngestionP_tot &
                                                                                       -Mortality*Me%OrganicMatter%PC_Ratio &
                                                                                       -Predation*Me%OrganicMatter%PC_Ratio &
                                                                                       -Respiration*Me%OrganicMatter%PC_Ratio                                                                              
                   
      
      ! Predation is a loss from the system, to avoid this it was added to organic matter
      
      Me%ExternalVar%Mass(PON, Index)=Me%ExternalVar%Mass(PON, Index)+EgestionN_tot + (Mortality +Predation)*Me%OrganicMatter%NC_Ratio ! PON e POP are preys, so the predation was calculated already before
      Me%ExternalVar%Mass(POP, Index)=Me%ExternalVar%Mass(POP, Index)+EgestionP_tot + (Mortality+Predation)*Me%OrganicMatter%PC_Ratio
      if (Me%PelagicModel == LifeModel) then
        Me%ExternalVar%Mass(POC, Index)=Me%ExternalVar%Mass(POC, Index)+EgestionC_tot + Mortality
      endif
      Me%ExternalVar%MassInKgFromWater(O2, Index)=Me%ExternalVar%MassInKgFromWater(O2, Index)-Respiration*32./12.
      Me%ExternalVar%MassInKgFromWater(AM, Index)=Me%ExternalVar%MassInKgFromWater(AM, Index)+Respiration*Consumer%NCRatio
      Me%ExternalVar%MassInKgFromWater(IP, Index)=Me%ExternalVar%MassInKgFromWater(IP, Index)+Respiration*Consumer%PCRatio
      
   
      
      
      Consumer => Consumer%Next
        end do d1

        
        end subroutine ComputeBenthicConsumers
    
    
    
    !-----------------------------------------------------------------------------------
    
        subroutine ComputeBenthicNitrogen(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: AM, PON, O2
        integer                                     :: PON1, PON2, PON3, PON4, PON5
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------

        AM  = Me%PropIndex%Ammonia
        PON = Me%PropIndex%PON
        O2  = Me%PropIndex%Oxygen
        

        
        if(Me%ComputeOptions%Pompools)then
            PON1 = Me%PropIndex%PON1
            PON2 = Me%PropIndex%PON2
            PON3 = Me%PropIndex%PON3
            PON4 = Me%PropIndex%PON4
            PON5 = Me%PropIndex%PON5        
        end if

        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max((Me%ExternalVar%MassInKgFromWater(O2,Index)/Me%ExternalVar%WaterVolume(index))*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Nitrogen%PONDecayRate       *  &
                             Me%Nitrogen%PONDecayTFactor    ** &
                            (Me%ExternalVar%Temperature(Index) - 20.0)

       
        !kgN * day * day-1 (what passes from PON to ammonia)
        Me%Matrix(Index, PON, AM) = Me%ExternalVar%Mass(PON, Index) * Me%DTDay * &
                                    MineralizationRate * OxygenLimitation

        if(Me%ComputeOptions%Pompools)then
            Me%Matrix(Index, PON1, AM) = Me%ExternalVar%Mass(PON1, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON2, AM) = Me%ExternalVar%Mass(PON2, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON3, AM) = Me%ExternalVar%Mass(PON3, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON4, AM) = Me%ExternalVar%Mass(PON4, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
            
            Me%Matrix(Index, PON5, AM) = Me%ExternalVar%Mass(PON5, Index) * Me%DTDay * &
                                         MineralizationRate * OxygenLimitation
        end if


        if(.NOT. Me%ComputeOptions%Pompools)then

            Me%ExternalVar%MassInKgFromWater(AM,  Index) = Me%ExternalVar%MassInKgFromWater(AM , Index) + Me%Matrix(Index, PON, AM)

            Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) - Me%Matrix(Index, PON, AM)

            !what is consumed of oxygen due to mineralization of PON
            OxygenConsumption               = Me%Matrix(Index, PON, AM) * 1. / Me%OrganicMatter%NC_Ratio * &
                                              32. / 12.

        else
         
            Me%ExternalVar%MassInKgFromWater(AM,  Index) = Me%ExternalVar%MassInKgFromWater(AM , Index) + Me%Matrix(Index, PON, AM) + &
                                              Me%Matrix(Index, PON1, AM) + Me%Matrix(Index, PON2, AM)     + &
                                              Me%Matrix(Index, PON3, AM) + Me%Matrix(Index, PON4, AM)     + &
                                              Me%Matrix(Index, PON5, AM)  

            Me%ExternalVar%Mass(PON, Index) = Me%ExternalVar%Mass(PON, Index) - Me%Matrix(Index, PON, AM)
            
            Me%ExternalVar%Mass(PON1, Index) = Me%ExternalVar%Mass(PON1, Index) - Me%Matrix(Index, PON1, AM)
            Me%ExternalVar%Mass(PON2, Index) = Me%ExternalVar%Mass(PON2, Index) - Me%Matrix(Index, PON2, AM)
            Me%ExternalVar%Mass(PON3, Index) = Me%ExternalVar%Mass(PON3, Index) - Me%Matrix(Index, PON3, AM)
            Me%ExternalVar%Mass(PON4, Index) = Me%ExternalVar%Mass(PON4, Index) - Me%Matrix(Index, PON4, AM)
            Me%ExternalVar%Mass(PON5, Index) = Me%ExternalVar%Mass(PON5, Index) - Me%Matrix(Index, PON5, AM)

            !what is consumed of oxygen due to mineralization of PON
            OxygenConsumption               = (Me%Matrix(Index, PON, AM) + Me%Matrix(Index, PON1, AM)  + &
                                              Me%Matrix(Index, PON2, AM) + Me%Matrix(Index, PON3, AM)  + &
                                              Me%Matrix(Index, PON4, AM) + Me%Matrix(Index, PON5, AM)) * &
                                              1. / Me%OrganicMatter%NC_Ratio * 32. / 12.
        
        end if
        
        

        Me%ExternalVar%MassInKgFromWater(O2, Index ) = Me%ExternalVar%MassInKgFromWater(O2, Index ) - OxygenConsumption


    end subroutine ComputeBenthicNitrogen
    
    !--------------------------------------------------------------------------


    subroutine ComputeBenthicPhosphorus(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index
        
        !Local-----------------------------------------------------------------
        integer                                     :: IP, POP, O2
        integer                                     :: POP1, POP2, POP3, POP4, POP5
        real                                        :: MineralizationRate
        real                                        :: OxygenConsumption
        real                                        :: OxygenLimitation

        !Begin-----------------------------------------------------------------
        
        IP  = Me%PropIndex%Phosphate
        POP = Me%PropIndex%POP
        O2  = Me%PropIndex%Oxygen
        

        
        if(Me%ComputeOptions%Pompools)then
            POP1 = Me%PropIndex%POP1
            POP2 = Me%PropIndex%POP2
            POP3 = Me%PropIndex%POP3
            POP4 = Me%PropIndex%POP4
            POP5 = Me%PropIndex%POP5        
        end if
        
        
        !Multiplication by 1000 because oxygen units are given in g/l
        OxygenLimitation = max((Me%ExternalVar%MassInKgFromWater(O2,Index)/Me%ExternalVar%WaterVolume(index))*1000., Me%Oxygen%Minimum)

        !OxygenLimitation = 1 when Oxygen levels are high 
        !OxygenLimitation = 0 when Oxygen levels are low 
        OxygenLimitation = OxygenLimitation / (OxygenLimitation + 0.5)

        !day-1
        MineralizationRate = Me%Phosphorus%POPDecayRate       *  &
                             Me%Phosphorus%POPDecayTFactor    ** &
                            (Me%ExternalVar%Temperature(Index) - 20.0)


        !kgP * day * day-1 (what passes from POP to inorganic phosphorus)
        Me%Matrix(Index, POP, IP) = Me%ExternalVar%Mass(POP, Index) * Me%DTDay * &
                         MineralizationRate * OxygenLimitation
        
        
            !------------------------------------------POM POOLS                 
            if(Me%ComputeOptions%Pompools)then
            
                Me%Matrix(Index, POP1, IP) = Me%ExternalVar%Mass(POP1, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                
                Me%Matrix(Index, POP2, IP) = Me%ExternalVar%Mass(POP2, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                
                Me%Matrix(Index, POP3, IP) = Me%ExternalVar%Mass(POP3, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                
                Me%Matrix(Index, POP4, IP) = Me%ExternalVar%Mass(POP4, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation
                                             
                Me%Matrix(Index, POP5, IP) = Me%ExternalVar%Mass(POP5, Index) * Me%DTDay * &
                                             MineralizationRate * OxygenLimitation                                    
            end if


        if(.NOT. Me%ComputeOptions%Pompools)then
        
            Me%ExternalVar%MassInKgFromWater(IP,  Index) = Me%ExternalVar%MassInKgFromWater(IP , Index) + Me%Matrix(Index, POP, IP)

            Me%ExternalVar%Mass(POP, Index) = Me%ExternalVar%Mass(POP, Index) - Me%Matrix(Index, POP, IP)

            OxygenConsumption               = Me%Matrix(Index, POP, IP) * 1. / Me%OrganicMatter%PC_Ratio * &
                                              32. / 12.
        
        else
        
            Me%ExternalVar%MassInKgFromWater(IP,  Index) = Me%ExternalVar%MassInKgFromWater(IP , Index) + Me%Matrix(Index, POP, IP) +  &
                                              Me%Matrix(Index, POP1, IP) + Me%Matrix(Index, POP2, IP)     +  &
                                              Me%Matrix(Index, POP3, IP) + Me%Matrix(Index, POP4, IP)     +  &
                                              Me%Matrix(Index, POP5, IP)

            Me%ExternalVar%Mass(POP, Index) = Me%ExternalVar%Mass(POP, Index) - Me%Matrix(Index, POP, IP)
            
            Me%ExternalVar%Mass(POP1, Index) = Me%ExternalVar%Mass(POP1, Index) - Me%Matrix(Index, POP1, IP)
            Me%ExternalVar%Mass(POP2, Index) = Me%ExternalVar%Mass(POP2, Index) - Me%Matrix(Index, POP2, IP)
            Me%ExternalVar%Mass(POP3, Index) = Me%ExternalVar%Mass(POP3, Index) - Me%Matrix(Index, POP3, IP)
            Me%ExternalVar%Mass(POP4, Index) = Me%ExternalVar%Mass(POP4, Index) - Me%Matrix(Index, POP4, IP)
            Me%ExternalVar%Mass(POP5, Index) = Me%ExternalVar%Mass(POP5, Index) - Me%Matrix(Index, POP5, IP)
            

            OxygenConsumption               = (Me%Matrix(Index, POP, IP) + Me%Matrix(Index, POP1, IP)  +  &
                                              Me%Matrix(Index, POP2, IP) + Me%Matrix(Index, POP3, IP)  +  &
                                              Me%Matrix(Index, POP4, IP) + Me%Matrix(Index, POP5, IP)) *  &
                                              1. / Me%OrganicMatter%PC_Ratio * 32. / 12.

        end if
        
        !if(.NOT.Me%ComputeOptions%Nitrogen)then ! if N is activated, O2 consumption was already calculated for PON decomposition
        Me%ExternalVar%MassInKgFromWater(O2, Index ) = Me%ExternalVar%MassInKgFromWater(O2, Index ) - OxygenConsumption 
        !endif                 
                         
    end subroutine ComputeBenthicPhosphorus
    
    !--------------------------------------------------------------------------


    subroutine ComputeBenthicSilica(Index)

        !Arguments-------------------------------------------------------------
        integer, intent(IN)                         :: Index

        !Local-----------------------------------------------------------------
        integer                                     :: BioSi, Sil

        !Begin-----------------------------------------------------------------


        
        Sil   = Me%PropIndex%DissolvedSilica
        BioSi = Me%PropIndex%BioSilica

        !kg * day * day-1 (what passes from biogenic silica to inorganic dissolved silica)
        Me%Matrix(Index, BioSi, Sil) = Me%ExternalVar%Mass(BioSi, Index) * Me%DTDay * Me%Silica%BioSiDecayRate

        Me%ExternalVar%Mass(Sil, Index)   = Me%ExternalVar%Mass(Sil,   Index) + Me%Matrix(Index, BioSi, Sil)

        Me%ExternalVar%Mass(BioSi, Index) = Me%ExternalVar%Mass(BioSi, Index) - Me%Matrix(Index, BioSi, Sil)

    
    end subroutine ComputeBenthicSilica
    
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



        subroutine KillBenthicEcology(ObjBenthicEcologyID, STAT)

        !Arguments---------------------------------------------------------------
        integer                             :: ObjBenthicEcologyID              
        integer, optional, intent(OUT)      :: STAT

        !External----------------------------------------------------------------
        integer                             :: ready_              

        !Local-------------------------------------------------------------------
        integer                             :: STAT_, nUsers           

        !------------------------------------------------------------------------

        STAT_ = UNKNOWN_

        call Ready(ObjBenthicEcologyID, ready_)    

cd1 :   if (ready_ .NE. OFF_ERR_) then

            nUsers = DeassociateInstance(mBenthicEcology_,  Me%InstanceID)

            if (nUsers == 0) then

                deallocate(Me%PropertyList)

                !Deallocates Instance
                call DeallocateInstance ()


                ObjBenthicEcologyID = 0
                STAT_      = SUCCESS_

            end if
        else 
            STAT_ = ready_
        end if cd1

        if (present(STAT)) STAT = STAT_

    end subroutine KillBenthicEcology
        

    !------------------------------------------------------------------------
    
    
    subroutine DeallocateInstance ()

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        type (T_BenthicEcology), pointer          :: AuxObjBenthicEcology
        type (T_BenthicEcology), pointer          :: PreviousObjBenthicEcology

        !Updates pointers
        if (Me%InstanceID == FirstObjBenthicEcology%InstanceID) then
            FirstObjBenthicEcology => FirstObjBenthicEcology%Next
        else
            PreviousObjBenthicEcology => FirstObjBenthicEcology
            AuxObjBenthicEcology      => FirstObjBenthicEcology%Next
            do while (AuxObjBenthicEcology%InstanceID /= Me%InstanceID)
                PreviousObjBenthicEcology => AuxObjBenthicEcology
                AuxObjBenthicEcology      => AuxObjBenthicEcology%Next
            enddo

            !Now update linked list
            PreviousObjBenthicEcology%Next => AuxObjBenthicEcology%Next

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

    subroutine Ready (ObjBenthicEcology_ID, ready_) 

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthicEcology_ID
        integer                                     :: ready_

        !----------------------------------------------------------------------

        nullify (Me)

cd1:    if (ObjBenthicEcology_ID > 0) then
            call LocateObjBenthicEcology (ObjBenthicEcology_ID)
            ready_ = VerifyReadLock (mBenthicEcology_, Me%InstanceID)
        else
            ready_ = OFF_ERR_
        end if cd1

        !----------------------------------------------------------------------

    end subroutine Ready

    !--------------------------------------------------------------------------

    subroutine LocateObjBenthicEcology (ObjBenthicEcologyID)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjBenthicEcologyID

        !Local-----------------------------------------------------------------

        Me => FirstObjBenthicEcology
        do while (associated (Me))
            if (Me%InstanceID == ObjBenthicEcologyID) exit
            Me => Me%Next
        enddo

        if (.not. associated(Me)) stop 'ModuleBenthicEcology - LocateObjBenthicEcology - ERR01'

    end subroutine LocateObjBenthicEcology
    
end module ModuleBenthicEcology

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------



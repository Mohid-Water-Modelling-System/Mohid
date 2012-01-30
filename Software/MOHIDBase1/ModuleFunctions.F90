!------------------------------------------------------------------------------
!        IST/MARETEC, Water Modelling Group, Mohid modelling system
!------------------------------------------------------------------------------
!
! TITLE         : Mohid Model
! PROJECT       : Mohid Base 1
! MODULE        : Functions
! URL           : http://www.mohid.com
! AFFILIATION   : IST/MARETEC, Marine Modelling Group
! DATE          : May 2003
! REVISION      : Frank Braunschweig - v4.0
! DESCRIPTION   : Module which contains common functions
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

Module ModuleFunctions

    use ModuleGlobalData
    use ModuleTime
    use ModuleEnterData,        only : GetData
    use ModuleStopWatch,        only : StartWatch, StopWatch
#ifdef _ENABLE_CUDA    
    !JPW. Include C/C++/CUDA binding
    use ModuleCuda
#endif _ENABLE_CUDA    
    !$ use omp_lib !!
    
#if _USE_MPI
    use mpi
#endif    
 
    

    implicit none

    private

    !Parameter
    real,    parameter                              :: StefanBoltzmann          = 5.669e-08     ![W/m2/K4]
    integer, parameter                              :: SimpleHeight_            = 1
    integer, parameter                              :: ComplexHeight_           = 2

    !types---------------------------------------------------------------------
    !griflet
    public  :: T_THOMAS2D
    public  :: T_THOMAS
    public  :: T_D_E_F
    public  :: T_D_E_F2D
    public  :: T_VECGW

    !Functions-----------------------------------------------------------------
    
    !Matrix Operations
    public  :: SetMatrixValue
    public  :: SetMatrixValueAllocatable
    public  :: GetPointer
#ifdef _USE_SEQASSIMILATION
    public  :: InvSingularDiagMatrix2D
    public  :: CholeskyFactorization
    public  :: EOFAnalysis
    private :: NormalRand
    private :: UniformRand01
#endif
    public :: Pad

    !Linear systems solvers
#ifdef _USE_SEQASSIMILATION
    public  :: CholLinSystemSolver
    private :: TriangLinSystemSolver
#endif

    !Linear tridiagonal sytems solvers (Thomas algorithm)
    public  :: THOMAS_2D
    public  :: THOMAS_3D
    public  :: THOMASZ
    public  :: tridag

    !Linear tridiagonal sytems solvers tridiagonal system with cyclic boundaries
    public  :: tridag_cyclic

    !Linear iterative solver
    public  :: CGS2D
    public  :: LISOLVE
    public  :: SIPSOL
    public  :: BICGSTAB2D

    !Advection routines
    public  :: ComputeAdvectionFace
    public  :: ComputeAdvection1D_V2
    public  :: ComputeAdvection1D
    public  :: ComputeAdvection3D

    public  :: SaturatedVaporPressure
    public  :: LongWaveDownward
    public  :: LongWaveUpward
    public  :: LatentHeat
    public  :: SensibleHeat

    public  :: AerationFlux
    public  :: AerationFlux_CO2

    !Diffusion routines
    public  :: ComputeDiffusion3D

    public  :: OrlanskiCelerity2D
    public  :: OxygenSaturation
    public  :: OxygenSaturationHenry
    public  :: OxygenSaturationCeQualW2
    public  :: CO2PartialPressure
    public  :: CO2_K0
    public  :: Density
    public  :: Sigma
    public  :: SigmaWang
    public  :: SigmaUNESCO
    public  :: SigmaLeendertse
    public  :: SigmaUNESCOPressureCorrection
    public  :: SigmaJMD95PressureCorrection
    public  :: SigmaMel96PressureCorrection
    public  :: ConvertTemperature
    public  :: SpecificHeatUNESCO
    !Converts condutivity in salinity 
    public  :: SAL78 
    !Converts pressure in meters
    public  :: depth
    !Compute the sound of speed in sea water
    public  :: SVEL    

    !Coordinates of grid cells
    public  :: RODAXY
    public  :: FromCartesianToGrid
    public  :: FromGridToCartesian

    !Interpolation of a value in Time
    public  :: InterpolateValueInTime
    public  :: InterpolateMatrix2DInTime
    public  :: InterpolateMatrix3DInTime

    public  :: LinearInterpolation
    public  :: InterpolateLinearyMatrix2D
    public  :: InterpolateLinearyMatrix3D
    
    !Extrapolation 
    public  :: ExtraPol2DNearestCell 
    public  :: ExtraPol3DNearestCell 
    public  :: ExtraPol3DNearestCell_8
    
    public  :: FillMatrix2DNearestCell
    public  :: FillMatrix3DNearestCell
    
    !Reading of Time Keywords
    public  :: ReadTimeKeyWords

    !Reading of Property ID
    public  :: ConstructPropertyID

#ifdef _USE_SEQASSIMILATION
    !Checking if Property is from Hydrodynamic or WaterProperties module
    public  :: Check_Hydrodynamic_Property
    public  :: Check_Water_Property
#endif

    !Reading keywords from OnlineString
    public  :: GetDataOnlineString
    
    !Checks if number is odd or even
    public  :: IsOdd
    
    !Interpolation routines
    public  :: InterpolateProfile
    public  :: InterpolateProfileR8
    public  :: QuadraticInterpolProfile
    public  :: PolIntProfile
    public  :: polint

    !Polygon 
    public  :: RelativePosition4VertPolygon
   
    !Secant
    public  :: Secant

    !Light limitation factors for water quality 
    public  :: PhytoLightLimitationFactor           !Function
    
    !T90 decay time calculation methods
    public  :: ComputeT90_Canteras                  !Function
    public  :: ComputeT90_Chapra                    !Function
    

    public  ::  Normcrossprod

    !Coordinates Functions
    public  ::  ComputeGridZone
    private ::  GetLambda0
    private ::  GetEllipsoid
    public  ::  LatLonToUTM
    public  ::  UTMToLatLon
    public  ::  LatLonToLambertSP2
    public  ::  DistanceBetweenTwoGPSPoints

#ifdef _USE_PROJ4  
    public  ::  GeographicToCartesian
    public  ::  CartesianToGeographic
#endif

    !Compute settling velocity    
    public  :: SettlingVelocity

    !Bathymetry smoother
    public  :: SLPMIN
    public  :: SLPMIN2
    
    !OpenMP Chunk Size
    public  :: Chunk_K
    public  :: Chunk_J
    public  :: Chunk_I

    !TimeToString
    public  :: TimeToString

    !ChangeSuffix
    public  :: ChangeSuffix

    !Ordering functions
    private :: SortNumerically_3D

    private ::  QuadraticInterpolation
    
    public  :: maxival
    public  :: minival
    
    public  :: WGS84toGoogleMaps    
    interface  WGS84toGoogleMaps
        module procedure WGS84toGoogleMaps1D
        module procedure WGS84toGoogleMaps2D
    end interface  WGS84toGoogleMaps
    
    !types -------------------------------------------------------------------

    !griflet
    type T_D_E_F2D
        real   , pointer, dimension(: , :)  :: D
        real(8), pointer, dimension(: , :)  :: E
        real   , pointer, dimension(: , :)  :: F
    end type T_D_E_F2D

    type  T_D_E_F
        real   , pointer, dimension(: , : , :)  :: D
        real(8), pointer, dimension(: , : , :)  :: E
        real   , pointer, dimension(: , : , :)  :: F
#ifdef _USE_PAGELOCKED
        !Pointers needed to allocate pagelocked memory for faster CUDA transfers
        type(C_PTR)                             :: DPtr
        type(C_PTR)                             :: EPtr
        type(C_PTR)                             :: FPtr
#endif _USE_PAGELOCKED
    end type T_D_E_F

    !griflet
    type  T_VECGW
        real(8), pointer, dimension(:)          :: G                     !Auxiliar thomas arrays 
        real(8), pointer, dimension(:)          :: W                     !Auxiliar thomas arrays 
    end type   T_VECGW        

    !griflet
    type   T_THOMAS2D
        type(T_D_E_F2D), pointer                :: COEF2
        real, pointer, dimension(: , : )        :: TI
        type(T_VECGW), pointer, dimension(:)    :: VEC
    end type   T_THOMAS2D

    type ::   T_THOMAS
        type(T_D_E_F), pointer                  :: COEF3
        real, pointer, dimension(: , : , :)     :: TI
        type(T_VECGW), pointer, dimension(:)    :: VEC
    end type   T_THOMAS

    !interfaces -------------------------------------------------------------------
    
    !griflet: out with these non-working interfaces    
    interface THOMAS_2D
        module procedure THOMAS_2D_Original
        module procedure THOMAS_2D_NewType
    end interface THOMAS_2D

    interface THOMAS_3D
        module procedure THOMAS_3D_Original
        module procedure THOMAS_3D_NewType
    end interface THOMAS_3D
    
    interface THOMASZ
        module procedure THOMASZ_Original
        module procedure THOMASZ_NewType
    end interface THOMASZ
    
    interface QuadraticInterpolation
        module procedure QuadraticInterpolationR8
        module procedure QuadraticInterpolationR4
    end interface QuadraticInterpolation
    
    public:: interpolate3D
    interface interpolate3D
        module procedure interpolate3D_R4
        module procedure interpolate3D_R8
    end interface interpolate3D

    interface SetMatrixValue
        module procedure SetMatrixValues1D_I4_FromMatrix
        module procedure SetMatrixValues1D_R4_FromMatrix
        module procedure SetMatrixValues1D_I8_FromMatrix
        module procedure SetMatrixValues1D_R8_FromMatrix
        module procedure SetMatrixValues2D_I4_Constant
        module procedure SetMatrixValues2D_R4_Constant
        module procedure SetMatrixValues2D_R8_Constant
        module procedure SetMatrixValues2D_R4_FromMatrix
        module procedure SetMatrixValues2D_R8_FromMatrix
        module procedure SetMatrixValues2D_I4_FromMatrix
        module procedure SetMatrixValues3D_I4_Constant
        module procedure SetMatrixValues3D_R4_Constant
        module procedure SetMatrixValues3D_R8_Constant
        module procedure SetMatrixValues3D_R4_FromMatrix
        module procedure SetMatrixValues3D_R8_FromMatrix
        module procedure SetMatrixValues3D_I4_FromMatrix
    end interface SetMatrixValue

    interface SetMatrixValueAllocatable
        module procedure SetMatrixValues2D_R4_ConstantAllocatable
        module procedure SetMatrixValues2D_R4_FromMatrixAllocatable
        module procedure SetMatrixValues2D_R8_ConstantAllocatable
        module procedure SetMatrixValues2D_R8_FromMatrixAllocatable
        module procedure SetMatrixValues3D_R4_ConstantAllocatable
        module procedure SetMatrixValues3D_R8_ConstantAllocatable
        module procedure SetMatrixValues3D_R4_FromMatrixAllocatable
        module procedure SetMatrixValues3D_R8_FromMatrixAllocatable
    end interface SetMatrixValueAllocatable

    interface GetPointer
        module procedure GetPointer2D_I4
        module procedure GetPointer2D_R4
        module procedure GetPointer2D_R8
        module procedure GetPointer3D_R4
        module procedure GetPointer3D_R8        
    end interface

#ifdef _USE_MPI
    public ::  MPIKind

    interface MPIKind
        module procedure MPIKind0D
        module procedure MPIKind1D
        module procedure MPIKind2D
        module procedure MPIKind3D
    end interface MPIKind
    !include "mpif.f90"
#endif

    !griflet: ersatz functions to minval and maxval. The goal is to lift
    !potential stacksize bottlenecks.
    interface minival
        module procedure minival1D_R4
        module procedure minival2D_R4
        module procedure minival3D_R4
        module procedure minival1D_R8
        module procedure minival2D_R8
        module procedure minival3D_R8
        module procedure minival1D_I4
        module procedure minival2D_I4
        module procedure minival3D_I4
        module procedure minival1D_I8
        module procedure minival2D_I8
        module procedure minival3D_I8
    end interface minival

    interface maxival
        module procedure maxival1D_R4
        module procedure maxival2D_R4
        module procedure maxival3D_R4
        module procedure maxival1D_R8
        module procedure maxival2D_R8
        module procedure maxival3D_R8
        module procedure maxival1D_I4
        module procedure maxival2D_I4
        module procedure maxival3D_I4
        module procedure maxival1D_I8
        module procedure maxival2D_I8
        module procedure maxival3D_I8
    end interface maxival


    !--------------------------------------------------------------------------

    contains

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues1D_R4_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:), pointer                  :: Matrix
        type (T_Size1D), intent(in)                     :: Size
        real(4), dimension(:), pointer                  :: InMatrix
        integer, dimension(:), pointer, optional        :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_I(Size%ILB, Size%IUB)
                
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i) == 1) then
                    Matrix (i) = InMatrix(i)
                endif
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                Matrix (i) = InMatrix(i)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    


    end subroutine SetMatrixValues1D_R4_FromMatrix

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues1D_I4_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        integer(4), dimension(:), pointer                  :: Matrix
        type (T_Size1D)                                 :: Size
        integer(4), dimension(:), pointer                  :: InMatrix
        integer, dimension(:), pointer, optional        :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_I(Size%ILB, Size%IUB) 
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i) == 1) then
                    Matrix (i) = InMatrix(i)
                endif
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                Matrix (i) = InMatrix(i)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues1D_I4_FromMatrix

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues1D_R8_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:), pointer                  :: Matrix
        type (T_Size1D)                                 :: Size
        real(8), dimension(:), pointer                  :: InMatrix
        integer, dimension(:), pointer, optional        :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i
        integer                                         :: CHUNK 

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_I(Size%ILB, Size%IUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i) == 1) then
                    Matrix (i) = InMatrix(i)
                endif
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                Matrix (i) = InMatrix(i)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues1D_R8_FromMatrix

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues1D_I8_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        integer(8), dimension(:), pointer                  :: Matrix
        type (T_Size1D)                                 :: Size
        integer(8), dimension(:), pointer                  :: InMatrix
        integer, dimension(:), pointer, optional        :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i
        integer                                         :: CHUNK 

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_I(Size%ILB, Size%IUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i) == 1) then
                    Matrix (i) = InMatrix(i)
                endif
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do i = Size%ILB, Size%IUB
                Matrix (i) = InMatrix(i)
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues1D_I8_FromMatrix

    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_I4_Constant (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        integer(4), dimension(:, :), pointer            :: Matrix
        type (T_Size2D)                                 :: Size
        integer(4), intent (IN)                         :: Value
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = Value
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = Value
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_I4_Constant
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R4_Constant (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :), pointer               :: Matrix
        type (T_Size2D)                                 :: Size
        real(4), intent (IN)                            :: Value
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = Value
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = Value
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R4_Constant

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R4_ConstantAllocatable (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :), allocatable           :: Matrix
        type (T_Size2D)                                 :: Size
        real(4), intent (IN)                            :: Value
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = Value
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = Value
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R4_ConstantAllocatable
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R8_Constant (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), pointer               :: Matrix
        type (T_Size2D)                                 :: Size
        real(8), intent (IN)                            :: Value
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = Value
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = Value
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R8_Constant
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R8_ConstantAllocatable (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), allocatable           :: Matrix
        type (T_Size2D)                                 :: Size
        real(8), intent (IN)                            :: Value
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = Value
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = Value
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R8_ConstantAllocatable    

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R4_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :), pointer               :: Matrix
        type (T_Size2D)                                 :: Size
        real(4), dimension(:, :), pointer               :: InMatrix
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = InMatrix(i, j)
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = InMatrix(i, j)
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R4_FromMatrix

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R4_FromMatrixAllocatable (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :), allocatable           :: Matrix
        type (T_Size2D)                                 :: Size
        real(4), dimension(:, :), allocatable           :: InMatrix
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = InMatrix(i, j)
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = InMatrix(i, j)
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R4_FromMatrixAllocatable
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_I4_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        integer, dimension(:, :), pointer               :: Matrix
        type (T_Size2D)                                 :: Size
        integer, dimension(:, :), pointer               :: InMatrix
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = InMatrix(i, j)
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = InMatrix(i, j)
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_I4_FromMatrix

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R8_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), pointer               :: Matrix
        type (T_Size2D)                                 :: Size
        real(8), dimension(:, :), pointer               :: InMatrix
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = InMatrix(i, j)
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = InMatrix(i, j)
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R8_FromMatrix
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues2D_R8_FromMatrixAllocatable (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :), allocatable           :: Matrix
        type (T_Size2D)                                 :: Size
        real(8), dimension(:, :), allocatable           :: InMatrix
        integer, dimension(:, :), pointer, optional     :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_J(Size%JLB, Size%JUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j) == 1) then
                    Matrix (i, j) = InMatrix(i, j)
                endif
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j) = InMatrix(i, j)
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues2D_R8_FromMatrixAllocatable
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R8_FromMatrixAllocatable (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :, :), allocatable        :: Matrix
        type (T_Size3D), intent(in)                     :: Size
        real(8), dimension(:, :, :), allocatable        :: InMatrix
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK
        
        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)

        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = InMatrix(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = InMatrix(i, j, k)
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_R8_FromMatrixAllocatable
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_I4_Constant (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        integer(4), dimension(:, :, :), pointer         :: Matrix
        type (T_Size3D)                                 :: Size
        integer(4), intent (IN)                         :: Value
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = Value
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = Value
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_I4_Constant
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R4_Constant (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :, :), pointer            :: Matrix
        type (T_Size3D)                                 :: Size
        real(4), intent (IN)                            :: Value
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_K(Size%KLB, Size%KUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = Value
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = Value
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif

    end subroutine SetMatrixValues3D_R4_Constant
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R4_ConstantAllocatable (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :, :), allocatable        :: Matrix
        type (T_Size3D)                                 :: Size
        real(4), intent (IN)                            :: Value
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------
        
        CHUNK = CHUNK_K(Size%KLB, Size%KUB)
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = Value
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = Value
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif

    end subroutine SetMatrixValues3D_R4_ConstantAllocatable
        
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R8_Constant (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :, :), pointer            :: Matrix
        type (T_Size3D)                                 :: Size
        real(8), intent (IN)                            :: Value
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)   
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = Value
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, chunk)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = Value
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_R8_Constant

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R8_ConstantAllocatable (Matrix, Size, Value, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :, :), allocatable        :: Matrix
        type (T_Size3D)                                 :: Size
        real(8), intent (IN)                            :: Value
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK

        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)   
        
        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = Value
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, chunk)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = Value
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_R8_ConstantAllocatable
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R4_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :, :), pointer            :: Matrix
        type (T_Size3D)                                 :: Size
        real(4), dimension(:, :, :), pointer            :: InMatrix
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK
        
        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)

        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J, K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = InMatrix(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J, K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = InMatrix(i, j, k)
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_R4_FromMatrix

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R4_FromMatrixAllocatable (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(4), dimension(:, :, :), allocatable        :: Matrix
        type (T_Size3D)                                 :: Size
        real(4), dimension(:, :, :), allocatable        :: InMatrix
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK
        
        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)

        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J, K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = InMatrix(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J, K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = InMatrix(i, j, k)
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_R4_FromMatrixAllocatable
    
    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_R8_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:, :, :), pointer            :: Matrix
        type (T_Size3D)                                 :: Size
        real(8), dimension(:, :, :), pointer            :: InMatrix
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK
        
        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)

        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = InMatrix(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = InMatrix(i, j, k)
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif    

    end subroutine SetMatrixValues3D_R8_FromMatrix

    !--------------------------------------------------------------------------

    subroutine SetMatrixValues3D_I4_FromMatrix (Matrix, Size, InMatrix, MapMatrix)

        !Arguments-------------------------------------------------------------
        integer, dimension(:, :, :), pointer            :: Matrix
        type (T_Size3D)                                 :: Size
        integer, dimension(:, :, :), pointer            :: InMatrix
        integer, dimension(:, :, :), pointer, optional  :: MapMatrix

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        integer                                         :: CHUNK
        
        !Begin-----------------------------------------------------------------

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)

        if (present(MapMatrix)) then
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                if (MapMatrix(i, j, k) == 1) then
                    Matrix (i, j, k) = InMatrix(i, j, k)
                endif
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL PRIVATE(I,J,K)
            !$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%JLB, Size%JUB
            do i = Size%ILB, Size%IUB
                Matrix (i, j, k) = InMatrix(i, j, k)
            enddo
            enddo
            enddo
            !$OMP END DO NOWAIT
            !$OMP END PARALLEL
        endif

    end subroutine SetMatrixValues3D_I4_FromMatrix

    !--------------------------------------------------------------------------

    function GetPointer2D_I4(matrix) result(ptr)
    
        !Arguments-------------------------------------------------------------
        integer(4), dimension(:,:), allocatable, intent(in), target    :: matrix
        
        !Local-----------------------------------------------------------------
        integer(4), dimension(:,:), pointer                            :: ptr
        
        ptr => matrix
    
    end function GetPointer2D_I4

    !--------------------------------------------------------------------------

    function GetPointer2D_R4(matrix) result(ptr)
    
        !Arguments-------------------------------------------------------------
        real(4), dimension(:,:), allocatable, intent(in), target    :: matrix
        
        !Local-----------------------------------------------------------------
        real(4), dimension(:,:), pointer                            :: ptr
        
        ptr => matrix
    
    end function GetPointer2D_R4
    
    !--------------------------------------------------------------------------

    function GetPointer2D_R8(matrix) result(ptr)
    
        !Arguments-------------------------------------------------------------
        real(8), dimension(:,:), allocatable, intent(in), target    :: matrix
        
        !Local-----------------------------------------------------------------
        real(8), dimension(:,:), pointer                            :: ptr
        
        ptr => matrix
    
    end function GetPointer2D_R8

    !--------------------------------------------------------------------------

    function GetPointer3D_R4(matrix) result(ptr)
    
        !Arguments-------------------------------------------------------------
        real(4), dimension(:,:,:), allocatable, intent(in), target  :: matrix
        
        !Local-----------------------------------------------------------------
        real(4), dimension(:,:,:), pointer                          :: ptr
        
        ptr => matrix
    
    end function GetPointer3D_R4
    
    !--------------------------------------------------------------------------

    function GetPointer3D_R8(matrix) result(ptr)

        !Arguments-------------------------------------------------------------
        real(8), dimension(:,:,:), allocatable, intent(in), target :: matrix
        
        !Local-----------------------------------------------------------------
        real(8), dimension(:,:,:), pointer             :: ptr
        
        ptr => matrix
    
    end function GetPointer3D_R8
    
    !--------------------------------------------------------------------------
    ! Function Pad
    ! Returns an upperbound value for padding a matrix.
    ! Should be applied to the first dimension of the matrix to get a performance gain.
    ! Note that this uses slightly more memory.
    ! Written by Jonathan van der Wielen, 2012-01-27, jonathan.vanderwielen@hydrologic.com
    !--------------------------------------------------------------------------
    function Pad(LowerBound, UpperBound) result (PaddedUpperBound)
    
        !Arguments-------------------------------------------------------------
        integer, intent(in)     :: LowerBound, UpperBound
        
        !Local-----------------------------------------------------------------
        integer                 :: PaddedUpperBound
        
#ifdef _PAD_MATRICES
        ! Padding a matrix to a multiple of 128 bytes gives a significant performance gain
        ! Padding to 32 values ensures a padding of 128 bytes for integer / single and 256 bytes for double
        ! Also take LowerBound into account when padding (is 0 by default, but who knows what might change in the future)
        PaddedUpperBound = int(ceiling((UpperBound - LowerBound) / 32.0)) * 32 + LowerBound
#else
        PaddedUpperBound = UpperBound
#endif _PAD_MATRICES

    end function Pad

#ifdef _USE_SEQASSIMILATION

    !--------------------------------------------------------------------------

    subroutine InvSingularDiagMatrix2D(Matrix, Dim, InvMatrix) 

        !Arguments---------------------------------------------------------------
        real, dimension (:, :), pointer             :: Matrix
        integer                                     :: Dim
        real, dimension (:, :), pointer             :: InvMatrix

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: i

        !------------------------------------------------------------------------

        !Output matrix is already allocated
        do i = 1, Dim

            InvMatrix(i,i) = 1/Matrix(i,i)

        enddo

    end subroutine InvSingularDiagMatrix2D

    !----------------------------------------------------------------------

    subroutine CholLinSystemSolver(A, B, Dim, X) 

        !Arguments---------------------------------------------------------------
        !real(8), dimension (:, :), pointer          :: A
        real, dimension (:, :), pointer             :: A
        integer                                     :: Dim
        real(8), dimension (:, :), pointer          :: B, X

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        real(8), dimension (:, :), pointer          :: AuxVector
        !real(8), dimension (:, :), pointer          :: AuxDiagMatrix1, AuxDiagMatrix2
        real, dimension (:, :), pointer             :: AuxDiagMatrix1, AuxDiagMatrix2
        logical                                     :: Lower

        !------------------------------------------------------------------------

        !Solves a linear system of equations, A*X = B based on the Cholesky decomposition 
        !of the system matrix A (symmetric positive definite) in a lower triangular 
        !square root matrix, AuxDiagMatrix1.

        !allocate aux variables
        allocate (AuxVector         (1:Dim, 1:1))
        allocate (AuxDiagMatrix1    (1:Dim, 1:Dim))
        allocate (AuxDiagMatrix2    (1:Dim, 1:Dim))

        AuxVector       (:,:)   = FillValueReal
        AuxDiagMatrix1  (:,:)   = 0.
        AuxDiagMatrix2  (:,:)   = 0.

        !Factorize A by Cholesky method
        !(no test is made for symmetric positive definite character)
        call CholeskyFactorization(A, Dim, AuxDiagMatrix1)

        !Solve AuxDiagMatrix1*AuxVector = B
        Lower = .true. !AuxDiagMatrix1 is lower triangular
        call TriangLinSystemSolver(AuxDiagMatrix1, B, Dim, Lower, AuxVector)

        !Solve AuxDiagMatrix1'*x = AuxVector
        AuxDiagMatrix2 = TRANSPOSE(AuxDiagMatrix1)
        Lower = .false. !AuxDiagMatrix2 is upper triangular
        call TriangLinSystemSolver(AuxDiagMatrix2, AuxVector, Dim, Lower, X) 

        !deallocate aux variables
        deallocate (AuxVector)
        deallocate (AuxDiagMatrix1)
        deallocate (AuxDiagMatrix2)

    end subroutine CholLinSystemSolver

    !----------------------------------------------------------------------

    subroutine CholeskyFactorization(A, Dim, L) 

        !Arguments---------------------------------------------------------------
        !real(8), dimension (:, :), pointer          :: A, L
        real, dimension (:, :), pointer             :: A, L
        integer                                     :: Dim

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: j, i, k
        real                                        :: AuxSum

        !------------------------------------------------------------------------

        !Computes the Cholesky factorization of a real symmetric 
        !positive definite matrix A = L*L'
        !(no test is made for symmetric positive definite character)

        !First column of L
        L(1,1) = sqrt(A(1,1))

        do i = 2, Dim

            L(i,1) = A(1,i)/L(1,1)

        enddo

        !Next columns of L
        do j = 2, Dim

            do i = j, Dim
            
                if (i == j) then
                
                    AuxSum = 0.

                    do k = 1, j-1
                    
                        AuxSum = AuxSum + (L(i,k))**2        
 
                    end do

                    L(i,i) = sqrt(A(i,i) - AuxSum)

                else

                    AuxSum = 0.

                    do k = 1, j-1

                        AuxSum = AuxSum + L(j,k)*L(i,k)

                    end do

                    L(i,j) = (A(j,i) - AuxSum)/L(j,j)

                endif

            end do

        end do

        !Reference: math-linux.com, http://www.math-linux.com/spip.php?article43

    end subroutine CholeskyFactorization

    !----------------------------------------------------------------------

    subroutine TriangLinSystemSolver(L, B, Dim, Lower, X) 

        !Arguments---------------------------------------------------------------
        !real(8), dimension (:, :), pointer          :: L
        real, dimension (:, :), pointer             :: L
        real(8), dimension (:, :), pointer          :: B, X
        integer                                     :: Dim
        logical                                     :: Lower

        !External----------------------------------------------------------------

        !Local-------------------------------------------------------------------
        integer                                     :: i, k
        real(8)                                     :: AuxSum

        !------------------------------------------------------------------------

        !Solves a linear system of equations L*X = B where L is a triangular matrix.
 
        if (Lower) then

            do i = 1, Dim

                AuxSum = 0.

                do k = 1, i-1

                    AuxSum = AuxSum + L(i,k)*X(k,1)

                enddo

                X(i,1) = (B(i,1) - AuxSum)/L(i,i)

            enddo

        else !Upper triangular matrix

            do i = Dim, 1, -1

                AuxSum = 0.

                do k = i+1, Dim
        
                    AuxSum = AuxSum + L(i,k)*X(k,1)

                enddo

                X(i,1) = (B(i,1) - AuxSum)/L(i,i)

            enddo

        endif

    end subroutine TriangLinSystemSolver

    !--------------------------------------------------------------------------

    subroutine EOFAnalysis(Covariance, LMatrix, UMatrix, rank, dim)

        !Arguments-------------------------------------------------------------
        real, dimension(:, :),  pointer                 :: Covariance
        real(8), dimension (:,:), pointer               :: LMatrix
        real, dimension (:,:), pointer                  :: UMatrix
        integer                                         :: rank, dim

        !Local-----------------------------------------------------------------
        integer                                         :: idum, ieig, k, iter
        real(8)                                         :: angle
        real                                            :: eps
        integer                                         :: itermax
        real(8), dimension(:), pointer                  :: vec2, vecr
        integer                                         :: j
        real(8)                                         :: tmp, a1, a2, b2, b, t, t2

        integer, parameter                              :: miter = 10000
        real, parameter                                 :: prec = 1e-8

        !----------------------------------------------------------------------

        !EOFAnalysis performs the EOF analysis of Covariance
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output,
        !
        ! Local parameters:
        !
        !

        allocate(vecr(1:dim))
        allocate(vec2(1:dim))

        idum = 1
        itermax = miter - 1

        write(*,*)
        write(*,*) 'Performing EOF analysis...' 

        do ieig = 1, rank
            
            eps = ieig*prec

            !Initialize each EOF as random
            do k = 1, dim !number of covariance matrix lines/columns
                LMatrix(k,ieig) = NormalRand(idum)
            end do

            !Calculate each EOF
            !call IEigenPower(ieig, miter - 1, ieig*prec, angle, iter)

            vecr(:) = FillValueReal
            vec2(:) = FillValueReal

            !1. Make new ieig-th eigen vector ortogonal with the others (are normalized)
            do j = 1,ieig - 1
                tmp = 0
                do k = 1, dim
                    tmp = tmp + LMatrix(k,ieig)*LMatrix(k,j)
                end do
                do k = 1, dim
                    LMatrix(k,ieig) = LMatrix(k,ieig) - tmp*LMatrix(k,j)
                end do
            end do

            b2 = 0.

            do k = 1, dim
                b2 = b2 + LMatrix(k,ieig)**2
            end do

            b = dsqrt(b2)
  
            do k = 1, dim
                LMatrix(k,ieig) = LMatrix(k,ieig)/b
            end do

            vecr(:) = MATMUL(Covariance, LMatrix(:,ieig))

            a1 = 0.
            do k = 1, dim
                a1 = a1 + vecr(k)*LMatrix(k,ieig)
            end do

            do k = 1, dim
                vecr(k) = vecr(k) - a1*LMatrix(k,ieig)
            end do

            iter = 1

            !2. Make vecr ortogonal with the other eigen vectors (are normalized)
do1 :       do 
                do j = 1, ieig - 1 !500
                    tmp = 0
                    do k = 1, dim
                        tmp = tmp + vecr(k)*LMatrix(k,j)

                    end do
                    do k = 1, dim
                        vecr(k) = vecr(k) - tmp*LMatrix(k,j)

                    end do
                end do

                b2 = 0.
        
                do k = 1, dim
                    b2 = b2 + vecr(k)**2
                end do

                b = dsqrt(b2)

                !stoping test
                angle = b/a1
        
                if (angle .lt. eps) exit do1 !goto 600

                iter = iter + 1

                if (iter .gt. itermax) exit do1 !goto 600

                do k = 1, dim
                    vec2(k) = vecr(k)/b
                end do

                vecr = MATMUL(Covariance, vec2)

                a2 = 0
                do k = 1, dim
                    a2 = a2 + vec2(k)*vecr(k)
                end do

                t = 2*b/(dabs(a1 - a2) + dsqrt((a1 - a2)**2 + 4*b2))

                if ((a1.le.a2) .and. ((a1.ne. a2).or.(a1+a2.lt.0))) then
                    t = -t
                endif

                t2 = 1 + t**2

                if (a1**2 .gt. a2**2) then
                    a1 = (a1 + (t**2)*a2 + 2*t*b)/t2 
                    t2 = dsqrt(t2)
                    tmp = t/t2
            
                    do k = 1, dim
                       vecr(k) = (vecr(k) - a2*vec2(k) - b*LMatrix(k,ieig))*tmp
                       LMatrix(k,ieig) = (LMatrix(k,ieig) + t*vec2(k))/t2
                    end do
                else
                    a1 = (a2 + (t**2)*a1 - 2*t*b)/t2
                    t2 = dsqrt(t2)
            
                    do k = 1, dim
                       vecr(k) = (vecr(k) - a2*vec2(k) - b*LMatrix(k,ieig))/t2
                       LMatrix(k,ieig) = (vec2(k) - t*LMatrix(k,ieig))/t2
                    end do
                endif
            !goto 500
            end do do1

            tmp = dsqrt(a1*a1 + b2) !600
      
            do k = 1, dim
                LMatrix(k,ieig) = (a1*LMatrix(k,ieig) + vecr(k))/tmp
            end do

            UMatrix(ieig,ieig) = a1

        end do

        !Deallocate aux vectors      
        deallocate(vecr)
        deallocate(vec2)
        nullify(vecr)
        nullify(vec2)

    end subroutine EOFAnalysis

    !--------------------------------------------------------------------------

    real function UniformRand01(idum)

        !Arguments-------------------------------------------------------------
        integer                                         :: idum

        !Local-----------------------------------------------------------------
        integer                                         :: k
        integer, parameter                              :: IA = 16807
        integer, parameter                              :: IM = 2147483647
        integer, parameter                              :: IQ=127773
        integer, parameter                              :: IR=2836
        real, parameter                                 :: AM=1./IM

        !----------------------------------------------------------------------

        !UniformRand01 returns a random number uniformely distributed in interval [0, 1]
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output, real UniformRand01, random real in the range 0 to 1.
        !
        ! Local parameters:
        !
        !

        k = idum/IQ
        idum = IA*(idum - k*IQ) - IR*k
        if (idum .lt. 0) then 
            idum = idum + IM
        endif

        UniformRand01 = AM*idum

    end function UniformRand01

    !--------------------------------------------------------------------------

    real function NormalRand(idum)

        !Arguments-------------------------------------------------------------
        integer                                         :: idum

        !Local-----------------------------------------------------------------
        integer , save                                  :: iset
        real    , save                                  :: gset
        real                                            :: v1, v2, fac, rsq

        !----------------------------------------------------------------------

        !NormalRand returns a random number following normal distribution N(0,1)
        !
        ! Author: 
        !
        !   Adapted from Fortran 77 code by:
        !   Ibrahim Hoteit,
        !   Scripps Institution of Oceanography,
        !   Nierenberg Hall, Room 410,
        !   9500 Gilman Drive, Dept 0230,
        !   La Jolla, CA 92093-0230, USA,
        !   ihoteit@ucsd.edu
        !       
        ! Parameters:
        !
        !    Input,
        !
        !    Input/output,
        !
        !    Output, real NormalRand, random real with distribution N(0,1)
        !
        ! Local parameters:
        !
        !

        DATA     iset /0/ !initializes iset


        if (iset .eq. 0) then
!10          v1 = 2.0 * UniformRand01(idum) - 1.0
            !following the same as(Fortran 77): if ((rsq .ge. 1.0).or.(rsq .eq. 0)) goto 10
            do
                v1 = 2.0 * UniformRand01(idum) - 1.0
                v2 = 2.0 * UniformRand01(idum) - 1.0
            
                rsq = v1*v1 + v2*v2
                if ((rsq .lt. 1.0) .and. (rsq /= 0)) exit
            end do
            
            fac = sqrt(- 2.0*alog(rsq)/rsq) !alog = neperian logarithm
            iset = 1
            gset = v2*fac
            NormalRand = v1*fac
        else
            iset = 0
            NormalRand = gset
        endif

    end function NormalRand

#endif

    !--------------------------------------------------------------------------
    
    subroutine THOMAS_2D_Original(IJmin, IJmax,                                        &
                         JImin, JImax,                                        &
                         di,    dj,                                           &
                         DCoef_2D, ECoef_2D, FCoef_2D, TiCoef_2D,             &
                         ANSWER,                                              &
                         VECG, VECW)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: di,    dj
        real,    dimension(:,:), pointer            :: DCoef_2D, FCoef_2D, TiCoef_2D
        real(8), dimension(:,:), pointer            :: ECoef_2D
        real,    dimension(:,:), pointer            :: ANSWER
        real(8), dimension(:  ), pointer            :: VECG, VECW

        !Local-----------------------------------------------------------------
        integer :: IJ, JI, II, MM, I, J

        !----------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "THOMAS_2D")

do2 :   do IJ = IJmin, IJmax
            I = IJmin-1 + IJ*dj + di
            J = JImin-1 + IJ*di + dj
            VECW(JImin) = -FCoef_2D (I, J)/ECoef_2D(I, J)
            VECG(JImin) =  TiCoef_2D(I, J)/ECoef_2D(I, J)

do3 :       do JI=JImin+1,JImax+1
                I        = IJ*dj + JI*di
                J        = IJ*di + JI*dj

                VECW(JI) = -FCoef_2D(I,J) / (ECoef_2D(I,J) + DCoef_2D(I,J) * VECW(JI-1))

                VECG(JI) = (TiCoef_2D(I,J) - DCoef_2D(I,J) * VECG(JI-1))/   &
                           (ECoef_2D (I,J) + DCoef_2D(I,J) * VECW(JI-1))
            end do do3

            I = IJ * dj + (JImax+1) * di
            J = IJ * di + (JImax+1) * dj

            ANSWER(I, J) = VECG(JImax+1)

do1 :       do II = JImin+1, JImax+1
                MM = JImax+JImin+1-II
                I  = IJ*dj + MM*di
                J  = IJ*di + MM*dj

                ANSWER(I,J) = VECW(MM) * ANSWER(I+di,J+dj) + VECG(MM)
            end do do1
        end do do2

        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "THOMAS_2D")

    end subroutine THOMAS_2D_Original

    !--------------------------------------------------------------------------
    
    subroutine THOMAS_2D_NewType(IJmin, IJmax,                                &
                         JImin, JImax,                                        &
                         di,    dj,                                           &
                         THOMAS,                                              &
                         ANSWER)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: di,    dj
        !real,    dimension(:,:), pointer            :: DCoef_2D, FCoef_2D, TiCoef_2D
        !real(8), dimension(:,:), pointer            :: ECoef_2D
        type(T_THOMAS2D), pointer                   :: THOMAS
        real,    dimension(:,:), pointer            :: ANSWER

        !Local-----------------------------------------------------------------
        integer :: IJ, JI, II, MM, I, J

        !griflet
        type(T_VECGW), pointer                      :: VEC
        integer                                     :: TID
        !$ integer                                  :: CHUNK !

        !Begin-----------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "THOMAS_2D")
        
        !$ CHUNK = CHUNK_J(IJmin,IJmax) !

        !$OMP PARALLEL PRIVATE(TID,VEC,IJ,I,J,JI,II,MM)
        TID = 1
        !$ TID = 1 + omp_get_thread_num() !
        VEC => THOMAS%VEC(TID)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do2 :   do IJ = IJmin, IJmax
            I = IJmin-1 + IJ*dj + di
            J = JImin-1 + IJ*di + dj
            VEC%W(JImin) = -THOMAS%COEF2%F(I, J)/THOMAS%COEF2%E(I, J)
            VEC%G(JImin) =  THOMAS%TI(I, J)/THOMAS%COEF2%E(I, J)

do3 :       do JI=JImin+1,JImax+1
                I        = IJ*dj + JI*di
                J        = IJ*di + JI*dj

                VEC%W(JI) = -THOMAS%COEF2%F(I,J) / (THOMAS%COEF2%E(I,J) + THOMAS%COEF2%D(I,J) * VEC%W(JI-1))

                VEC%G(JI) = (THOMAS%TI(I,J) - THOMAS%COEF2%D(I,J) * VEC%G(JI-1))/   &
                           (THOMAS%COEF2%E(I,J) + THOMAS%COEF2%D(I,J) * VEC%W(JI-1))
            end do do3

            I = IJ * dj + (JImax+1) * di
            J = IJ * di + (JImax+1) * dj

            ANSWER(I, J) = VEC%G(JImax+1)

do1 :       do II = JImin+1, JImax+1
                MM = JImax+JImin+1-II
                I  = IJ*dj + MM*di
                J  = IJ*di + MM*dj

                ANSWER(I,J) = VEC%W(MM) * ANSWER(I+di,J+dj) + VEC%G(MM)
            end do do1
        end do do2
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "THOMAS_2D")

    end subroutine THOMAS_2D_NewType

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !griflet: New interface
    subroutine THOMAS_3D_Original(IJmin, IJmax,                                        &
                                  JImin, JImax,                                        &
                                  Kmin,  Kmax,                                         &
                                  di,    dj,                                           &
                                  DCoef_3D, ECoef_3D, FCoef_3D, TiCoef_3D,             &
                                  ANSWER,                                              &
                                  VECG, VECW)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: Kmin , Kmax
        integer,                         intent(IN) :: di,    dj
        real,    dimension(:,:,:), pointer          :: DCoef_3D, FCoef_3D, TiCoef_3D
        real(8), dimension(:,:,:), pointer          :: ECoef_3D
        real,    dimension(:,:,:), pointer          :: ANSWER
        real(8), dimension(:  ), pointer            :: VECG, VECW

        !Local-----------------------------------------------------------------
        
        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "THOMAS_3D")

        if (di == 0 .and. dj == 1) then

            call THOMAS_3D_i0_j1(IJmin, IJmax,                                        &
                                 JImin, JImax,                                        &
                                 Kmin,  Kmax,                                         &
                                 DCoef_3D, ECoef_3D, FCoef_3D, TiCoef_3D,             &
                                 ANSWER,                                              &
                                 VECG, VECW)

        else if (di == 1 .and. dj == 0) then

            call THOMAS_3D_i1_j0(IJmin, IJmax,                                        &
                                 JImin, JImax,                                        &
                                 Kmin,  Kmax,                                         &
                                 DCoef_3D, ECoef_3D, FCoef_3D, TiCoef_3D,             &
                                 ANSWER,                                              &
                                 VECG, VECW)

        else
            
            stop 'THOMAS_3D - Module Functions - ERR01'

        endif

        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "THOMAS_3D")

    end subroutine THOMAS_3D_Original

    !--------------------------------------------------------------------------

    subroutine THOMAS_3D_i0_j1 (IJmin, IJmax,                                        &
                                JImin, JImax,                                        &
                                Kmin,  Kmax,                                         &
                                DCoef_3D, ECoef_3D, FCoef_3D, TiCoef_3D,             &
                                ANSWER,                                              &
                                VECG, VECW)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: Kmin , Kmax
        real,    dimension(:,:,:), pointer          :: DCoef_3D, FCoef_3D, TiCoef_3D
        real(8), dimension(:,:,:), pointer          :: ECoef_3D
        real,    dimension(:,:,:), pointer          :: ANSWER
        real(8), dimension(:  ), pointer            :: VECG, VECW

        !Local-----------------------------------------------------------------
        integer                                     :: IJ, JI, II, K

        !Begin-----------------------------------------------------------------

do4:    do K  = Kmin, Kmax
do2 :   do IJ = IJmin, IJmax

            VECW(JImin) =-FCoef_3D (IJ, JImin, K)/ECoef_3D(IJ, JImin, K)
            VECG(JImin) = TiCoef_3D(IJ, JImin, K)/ECoef_3D(IJ, JImin, K)

do3 :       do JI=JImin+1,JImax+1
                VECW(JI) = - FCoef_3D(IJ, JI, K) / (ECoef_3D(IJ, JI, K) +                    &
                             DCoef_3D(IJ, JI, K) * VECW(JI-1))

                VECG(JI) =  (TiCoef_3D(IJ, JI, K) - DCoef_3D(IJ, JI, K) * VECG(JI-1))/      &
                            (ECoef_3D (IJ, JI, K) + DCoef_3D(IJ, JI, K) * VECW(JI-1))
            end do do3

            ANSWER(IJ, (JImax+1), K) = VECG(JImax+1)

do1 :       do II = JImin+1, JImax+1
                ANSWER(IJ, JImax+JImin+1-II, K) = VECW(JImax+JImin+1-II) * ANSWER(IJ, JImax+JImin+1-II+1, K) + &
                                                  VECG(JImax+JImin+1-II)
            end do do1
        end do do2
        end do do4

    end subroutine THOMAS_3D_i0_j1

    !--------------------------------------------------------------------------

    subroutine THOMAS_3D_i1_j0(IJmin, IJmax,                                        &
                               JImin, JImax,                                        &
                               Kmin,  Kmax,                                         &
                               DCoef_3D, ECoef_3D, FCoef_3D, TiCoef_3D,             &
                               ANSWER,                                              &
                               VECG, VECW)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: Kmin , Kmax
        real,    dimension(:,:,:), pointer          :: DCoef_3D, FCoef_3D, TiCoef_3D
        real(8), dimension(:,:,:), pointer          :: ECoef_3D
        real,    dimension(:,:,:), pointer          :: ANSWER
        real(8), dimension(:  ), pointer            :: VECG, VECW

        !Local-----------------------------------------------------------------
        integer                                     :: IJ, JI, II, K

        !Begin-----------------------------------------------------------------

do4:    do K  = Kmin, Kmax
do2 :   do IJ = IJmin, IJmax

            VECW(JImin) =-FCoef_3D (JImin, IJ, K)/ECoef_3D(JImin, IJ, K)
            VECG(JImin) = TiCoef_3D(JImin, IJ, K)/ECoef_3D(JImin, IJ, K)

do3 :       do JI=JImin+1,JImax+1
                VECW(JI) = - FCoef_3D(JI, IJ, K) / (ECoef_3D(JI, IJ, K) +                    &
                             DCoef_3D(JI, IJ, K) * VECW(JI-1))

                VECG(JI) = (TiCoef_3D(JI, IJ, K) - DCoef_3D(JI, IJ, K) * VECG(JI-1))/      &
                            (ECoef_3D (JI, IJ, K) + DCoef_3D(JI, IJ, K) * VECW(JI-1))
            end do do3

            ANSWER(JImax+1, IJ, K) = VECG(JImax+1)

do1 :       do II = JImin+1, JImax+1
                ANSWER(JImax+JImin+1-II, IJ, K) = VECW(JImax+JImin+1-II) * ANSWER(JImax+JImin+1-II+1, IJ, K) + &
                                                  VECG(JImax+JImin+1-II)
            end do do1
        end do do2
        end do do4

    end subroutine THOMAS_3D_i1_j0

    !--------------------------------------------------------------------------
    !griflet: New interface
    subroutine THOMAS_3D_NewType(IJmin, IJmax,                                          &
                                  JImin, JImax,                                         &
                                  Kmin,  Kmax,                                          &
                                  di,    dj,                                            &
                                  Thomas,                                               &
                                  ANSWER                                                &
#ifdef _ENABLE_CUDA
                                 , CudaID                                               &
                                 , SaveResults                                          &
#endif _ENABLE_CUDA 
                                )
        
        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: Kmin , Kmax
        integer,                         intent(IN) :: di,    dj
        type(T_THOMAS), pointer                     :: Thomas
        real,    dimension(:,:,:), pointer          :: ANSWER
#ifdef _ENABLE_CUDA
        ! Solve Thomas on a CUDA device
        integer                                     :: CudaID
        logical                                     :: SaveResults
#endif _ENABLE_CUDA

         !Local-----------------------------------------------------------------
        integer                                     :: Dim
       
        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "THOMAS_3D")
        
        if (di == 0 .and. dj == 1) then
            Dim = 1
            
#ifdef _USE_CUDA
            write(*,*) ' Solving Thomas for Y'
            ! If Y dimension, JImin / JImax = JLB / JUB. Dim = 1 = Y
            call SolveThomas(CudaID, IJmin, IJmax, JImin, JImax, Kmin, Kmax,          &
                Thomas%COEF3%D, Thomas%COEF3%E, Thomas%COEF3%F, Thomas%TI, ANSWER, Dim)
#else
            call THOMAS_3D_i0_j1_NewType(IJmin, IJmax,                                  &
                                 JImin, JImax,                                          &
                                 Kmin,  Kmax,                                           &
                                 Thomas,                                                &
                                 ANSWER)
#endif _USE_CUDA

        else if (di == 1 .and. dj == 0) then
            Dim = 0
#ifdef _USE_CUDA
            write(*,*) ' Solving Thomas for X'
            ! If X dimension, JImin / JImax = ILB / IUB. Dim = 0 = X
            call SolveThomas(CudaID, JImin, JImax, IJmin, IJmax, Kmin, Kmax,          &
                Thomas%COEF3%D, Thomas%COEF3%E, Thomas%COEF3%F, Thomas%TI, ANSWER, Dim)

#else
            call THOMAS_3D_i1_j0_NewType(IJmin, IJmax,                                  &
                                 JImin, JImax,                                          &
                                 Kmin,  Kmax,                                           &
                                 Thomas,                                                &
                                 ANSWER)
#endif _USE_CUDA

        else

            stop 'THOMAS_3D - Module Functions - ERR01'

        endif

        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "THOMAS_3D")

#ifdef _ENABLE_CUDA
        if(SaveResults) then
            ! Correctness test only
            ! call SaveThomas(CudaID, ANSWER, Dim)
        endif
#endif _ENABLE_CUDA

    end subroutine THOMAS_3D_NewType

    !--------------------------------------------------------------------------
    
    !griflet: new Thomas3D subroutine for openmp
    !--------------------------------------------------------------------------

    subroutine THOMAS_3D_i0_j1_NewType (IJmin, IJmax,                                   &
                                JImin, JImax,                                           &
                                Kmin,  Kmax,                                            &
                                Thomas,                                                 &
                                ANSWER)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: Kmin , Kmax
        type(T_THOMAS), pointer                     :: Thomas
        real,    dimension(:,:,:), pointer          :: ANSWER

        !Local-----------------------------------------------------------------
        integer                                     :: IJ, JI, II, K
        
        !griflet
        type(T_VECGW), pointer                      :: VEC
        integer                                     :: TID
        !$ integer                                  :: CHUNK !

        !Begin-----------------------------------------------------------------

        !$ CHUNK = CHUNK_K(Kmin,Kmax) !

        !$OMP PARALLEL PRIVATE(TID,VEC,K,IJ,JI,II)
        TID = 1
        !$ TID = 1 + omp_get_thread_num() !
        VEC => Thomas%VEC(TID)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do4:    do K  = Kmin, Kmax
do2 :   do IJ = IJmin, IJmax

            VEC%W(JImin) =-Thomas%COEF3%F (IJ, JImin, K)/Thomas%COEF3%E(IJ, JImin, K)
            VEC%G(JImin) = Thomas%TI(IJ, JImin, K)/Thomas%COEF3%E(IJ, JImin, K)

do3 :       do JI=JImin+1,JImax+1
                VEC%W(JI) = - Thomas%COEF3%F(IJ, JI, K) / (Thomas%COEF3%E(IJ, JI, K) +                    &
                             Thomas%COEF3%D(IJ, JI, K) * VEC%W(JI-1))

                VEC%G(JI) =  (Thomas%TI(IJ, JI, K) - Thomas%COEF3%D(IJ, JI, K) * VEC%G(JI-1))/      &
                            (Thomas%COEF3%E (IJ, JI, K) + Thomas%COEF3%D(IJ, JI, K) * VEC%W(JI-1))
            end do do3

            ANSWER(IJ, (JImax+1), K) = VEC%G(JImax+1)

do1 :       do II = JImin+1, JImax+1
                ANSWER(IJ, JImax+JImin+1-II, K) = VEC%W(JImax+JImin+1-II) * ANSWER(IJ, JImax+JImin+1-II+1, K) + &
                                                  VEC%G(JImax+JImin+1-II)
            end do do1
        end do do2
        end do do4
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine THOMAS_3D_i0_j1_NewType

    !--------------------------------------------------------------------------

    subroutine THOMAS_3D_i1_j0_NewType(IJmin, IJmax,                                &
                               JImin, JImax,                                        &
                               Kmin,  Kmax,                                         &
                               Thomas,                                              &
                               ANSWER)

        !Arguments-------------------------------------------------------------
        integer,                         intent(IN) :: IJmin, IJmax
        integer,                         intent(IN) :: JImin, JImax
        integer,                         intent(IN) :: Kmin , Kmax
        type(T_THOMAS),pointer                      :: Thomas
        real,    dimension(:,:,:), pointer          :: ANSWER

        !Local-----------------------------------------------------------------
        integer                                     :: IJ, JI, II, K

        !griflet
        type(T_VECGW), pointer                      :: VEC
        integer                                     :: TID
        !$ integer                                  :: CHUNK !

        !Begin-----------------------------------------------------------------

        !$ CHUNK = CHUNK_K(Kmin,Kmax) !

        !$OMP PARALLEL PRIVATE(TID,VEC,K,IJ,JI,II)
        TID = 1
        !$ TID = 1 + omp_get_thread_num() !
        VEC => Thomas%VEC(TID)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do4:    do K  = Kmin, Kmax
do2 :   do IJ = IJmin, IJmax

            VEC%W(JImin) =-Thomas%COEF3%F(JImin, IJ, K)/Thomas%COEF3%E(JImin, IJ, K)
            VEC%G(JImin) = Thomas%TI(JImin, IJ, K)/Thomas%COEF3%E(JImin, IJ, K)

do3 :       do JI=JImin+1,JImax+1
                VEC%W(JI) = - Thomas%COEF3%F(JI, IJ, K) / (Thomas%COEF3%E(JI, IJ, K) +                    &
                             Thomas%COEF3%D(JI, IJ, K) * VEC%W(JI-1))

                VEC%G(JI) = (Thomas%TI(JI, IJ, K) - Thomas%COEF3%D(JI, IJ, K) * VEC%G(JI-1))/      &
                            (Thomas%COEF3%E(JI, IJ, K) + Thomas%COEF3%D(JI, IJ, K) * VEC%W(JI-1))
            end do do3

            ANSWER(JImax+1, IJ, K) = VEC%G(JImax+1)

do1 :       do II = JImin+1, JImax+1
                ANSWER(JImax+JImin+1-II, IJ, K) = VEC%W(JImax+JImin+1-II) * ANSWER(JImax+JImin+1-II+1, IJ, K) + &
                                                  VEC%G(JImax+JImin+1-II)
            end do do1

        end do do2
        end do do4
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL

    end subroutine THOMAS_3D_i1_j0_NewType
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !griflet: New public interface
    subroutine THOMASZ_Original (ILB, IUB,                                        &
                                 JLB, JUB,                                        &
                                 KLB, KUB,                                        &
                                 AW, BW, CW, TIW,                                 &
                                 RES,                                             &
                                 VECG, VECW)

        !Arguments---------------------------------------------------------------
        integer,                         intent(IN) :: ILB, IUB
        integer,                         intent(IN) :: JLB, JUB
        integer,                         intent(IN) :: KLB, KUB
        real,    dimension(:,:,:), pointer          :: AW, CW, TIW
        real(8), dimension(:,:,:), pointer          :: BW
        real,    dimension(:,:,:), pointer          :: RES
        real(8), dimension(:    ), pointer          :: VECG, VECW

        !Local-------------------------------------------------------------------
        integer :: I, J, K
        integer :: II, MM

        !------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "THOMASZ")

do2 :   DO J = JLB, JUB
do1 :   DO I = ILB, IUB
            VECW(KLB) =-CW (I, J, 1) / BW(I, J, 1)
            VECG(KLB) = TIW(I, J, 1) / BW(I, J, 1)

do3 :       DO K  = KLB+1, KUB+1
                VECW(K) = -CW(I, J, K) / (BW(I, J, K) + AW(I, J, K) * VECW(K-1))

                VECG(K) = (TIW(I, J, K) - AW(I, J, K) * VECG(K-1))            &
                         / (BW (I, J, K) + AW(I, J, K) * VECW(K-1))
            END DO do3

            RES(I, J, KUB+1) = VECG(KUB+1)

do4 :       DO II = KLB+1, KUB+1
                MM            = KUB + KLB + 1 - II
                RES(I, J, MM) = VECW(MM) * RES(I, J, MM+1) + VECG(MM)
            END DO do4
        END DO do1
        END DO do2
        
        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "THOMASZ")

    end subroutine THOMASZ_Original

    !--------------------------------------------------------------------------
    !griflet: New public interface
    subroutine THOMASZ_NewType (ILB, IUB,                                         &
                                 JLB, JUB,                                        &
                                 KLB, KUB,                                        &
                                 Thomas,                                          &
                                 RES                                              &
#ifdef _ENABLE_CUDA
                                 , CudaID                                         &
                                 , SaveResults                                    &
#endif _ENABLE_CUDA 
                                )

        !Arguments---------------------------------------------------------------
        integer,                         intent(IN) :: ILB, IUB
        integer,                         intent(IN) :: JLB, JUB
        integer,                         intent(IN) :: KLB, KUB
        real,    dimension(:,:,:), pointer          :: RES
        type(T_THOMAS), pointer                     :: Thomas

#ifdef _ENABLE_CUDA
        ! Solve Thomas on a CUDA device
        integer                                     :: CudaID
        logical                                     :: SaveResults
#endif _ENABLE_CUDA

        !Local-------------------------------------------------------------------
        type(T_VECGW), pointer                      :: VEC
        integer                                     :: TID
        !$ integer                                  :: CHUNK !
        integer :: I, J, K
        integer :: II, MM     

        !------------------------------------------------------------------------

        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "THOMASZ")
        
            
#ifdef _USE_CUDA
            write(*,*) ' Solving Thomas for Z'
            
        ! This method can solve Thomas for any dimension. Dim 0 = X, 1 = Y, 2 = Z
        call SolveThomas(CudaID, ILB, IUB, JLB, JUB, KLB, KUB,                     &
                          Thomas%COEF3%D, Thomas%COEF3%E, Thomas%COEF3%F,           &
                          Thomas%TI, RES, 2)
#else
        !$ CHUNK = CHUNK_J(JLB,JUB) !

        !$OMP PARALLEL PRIVATE(J,I,K,II,MM,TID,VEC)
        TID = 1
        !$ TID = 1 + omp_get_thread_num() !
        VEC => Thomas%VEC(TID)
        !$OMP DO SCHEDULE(DYNAMIC,CHUNK)
do2 :   DO J = JLB, JUB
do1 :   DO I = ILB, IUB
            ! JPW: Changed 1 to KLB for consistency. Results should be the same as long as KLB = 1
            !VEC%W(KLB) =-Thomas%COEF3%F (I, J, KLB) / Thomas%COEF3%E(I, J, KLB)
            !VEC%G(KLB) = Thomas%TI(I, J, KLB) / Thomas%COEF3%E(I, J, KLB)
            ! JPW: Original
            VEC%W(KLB) =-Thomas%COEF3%F (I, J, 1) / Thomas%COEF3%E(I, J, 1)
            VEC%G(KLB) = Thomas%TI(I, J, 1) / Thomas%COEF3%E(I, J, 1)

do3 :       DO K  = KLB+1, KUB+1
                VEC%W(K) = -Thomas%COEF3%F(I, J, K) / (Thomas%COEF3%E(I, J, K) + Thomas%COEF3%D(I, J, K) * VEC%W(K-1))

                VEC%G(K) = (Thomas%TI(I, J, K) - Thomas%COEF3%D(I, J, K) * VEC%G(K-1))            &
                         / (Thomas%COEF3%E (I, J, K) + Thomas%COEF3%D(I, J, K) * VEC%W(K-1))
            END DO do3

            RES(I, J, KUB+1) = VEC%G(KUB+1)

do4 :       DO II = KLB+1, KUB+1
                MM            = KUB + KLB + 1 - II
                RES(I, J, MM) = VEC%W(MM) * RES(I, J, MM+1) + VEC%G(MM)
            END DO do4
        END DO do1
        END DO do2
        !$OMP END DO NOWAIT
        !$OMP END PARALLEL
       
#endif
        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "THOMASZ")

#ifdef _ENABLE_CUDA
        if(SaveResults) then
            ! Correctness test only
            ! call SaveThomas(CudaID, RES, 2);
        endif
#endif
    end subroutine THOMASZ_NewType

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    Subroutine OrlanskiCelerity2D(NewField, OldField, ReferenceField, ComputePoints,     &
                                  IMin, IMax, JMin, JMax, di, dj, i, j, k,               &
                                  EastNorthBoundary, LimitMax,                           &
                                  DT, TrelaxIn, TrelaxOut, FlowVelX,                     &
                                  InternalCelerity, WaveDirectionX, WaveDirectionY,      & 
                                  ConstantCelerity, NormalRadiation, Explicit,           &
                                  WaveCelerityX, WaveCelerityY, NewValue)

        !Arguments------------------------------------------------------------------------
        real   , pointer, dimension(:,:,:)           :: NewField, OldField
        real   , pointer, dimension(:,:,:), optional :: ReferenceField
        integer, pointer, dimension(:,:,:)           :: ComputePoints
        real   , intent(IN)                          :: DT, LimitMax
        integer, intent(IN)                          :: IMin, IMax, JMin, JMax, di, dj, i, j, k
        logical, intent(IN)                          :: EastNorthBoundary
        logical, intent(IN),                optional :: ConstantCelerity, NormalRadiation, Explicit
        real   , intent(IN),                optional :: InternalCelerity, WaveDirectionX, WaveDirectionY
        real   , intent(IN),                optional :: TrelaxIn, TrelaxOut, FlowVelX
        real   , intent(OUT),               optional :: WaveCelerityX, WaveCelerityY, NewValue


        !Local----------------------------------------------------------------------------
        real                                         :: Interior1New, Interior1Old, Interior2New, &
                                                        Interior2Old, Interior3Old, Interior3New, &    
                                                        Boundary3Old, Boundary4Old, Boundary5Old, &
                                                        Boundary32Old, Boundary33Old,             &
                                                        Boundary42Old, Boundary43Old,             &
                                                        Interior6Old, Interior7Old
        real                                         :: Trelax, TrelaxIn_, TrelaxOut_, TimeVariability
        real                                         :: AdjacentPropX, AdjacentPropY, ReferenceValue,  &   
                                                        SpaceSquare, SpaceVariabilityX, SpaceVariabilityY 
        real                                         :: AuxInt, AuxIntY, AuxBound
        real                                         :: WaveCelerityX_, WaveCelerityY_
        logical                                      :: NoSouthWest, NoNorthEast
        logical                                      :: ConstantCelerity_, NormalRadiation_, Explicit_


        !Begin----------------------------------------------------------------------------

        NoSouthWest = .false.
        NoNorthEast = .false.

        if (EastNorthBoundary) then
    
            Interior1New  = NewField (i - di, j - dj, k)

            Interior1Old  = OldField (i - di, j - dj, k) 

            Interior2New  = NewField (i - 2 * di, j - 2 * dj, k) 

            Interior3New  = NewField (i - 3 * di, j - 3 * dj, k) 

            Interior2Old  = OldField (i - 2 * di, j - 2 * dj, k) 

            Interior3Old  = OldField (i - 3 * di, j - 3 * dj, k) 

            Interior6Old  = OldField (i - di + dj, j - dj + di, k) 

            Interior7Old  = OldField (i - di - dj, j - dj - di, k) 

            if ((i - di + 3 * dj) > IMax .or. (j - dj + 3 * di) > JMax) then
            
                NoNorthEast = .true. 

            else

                if ((ComputePoints(i - di +     dj, j - dj +     di, k)   *              &
                     ComputePoints(i - di + 2 * dj, j - dj + 2 * di, k)   *              & 
                     ComputePoints(i - di + 3 * dj, j - dj + 3 * di, k)) == 0 )          &
                     NoNorthEast = .true.

            endif

            if ((i - di - 3 * dj) < IMin .or. (j - dj - 3 * di) < JMin) then
            
                NoSouthWest = .true. 

            else

                if ((ComputePoints(i - di -     dj, j - dj -     di, k)   *              &
                     ComputePoints(i - di - 2 * dj, j - dj - 2 * di, k)   *              & 
                     ComputePoints(i - di - 3 * dj, j - dj - 3 * di, k)) == 0 )          &
                     NoSouthWest = .true.

            endif


        else 

            Interior1New  = NewField (i + di, j + dj, k)   

            Interior1Old  = OldField (i + di, j + dj, k) 
        
            Interior2New  = NewField (i + 2 * di, j + 2 * dj, k) 

            Interior3New  = NewField (i + 3 * di, j + 3 * dj, k) 

            Interior2Old  = OldField (i + 2 * di, j + 2 * dj, k) 

            Interior3Old  = OldField (i + 3 * di, j + 3 * dj, k) 

            Interior6Old  = OldField (i + di + dj, j + dj + di, k) 

            Interior7Old  = OldField (i + di - dj, j + dj - di, k) 

            if ((i + di + 3 * dj) > IMax .or. (j + dj + 3 * di) > JMax) then
            
                NoNorthEast = .true. 

            else

                if ((ComputePoints(i + di +     dj, j + dj +     di, k)   *              &
                     ComputePoints(i + di + 2 * dj, j + dj + 2 * di, k)   *              & 
                     ComputePoints(i + di + 3 * dj, j + dj + 3 * di, k)) == 0 )          &
                     NoNorthEast = .true.

            endif

            if ((i + di - 3 * dj) < IMin .or. (j + dj - 3 * di) < JMin) then
            
                NoSouthWest = .true. 

            else

                if ((ComputePoints(i + di -     dj, j + dj -     di, k)   *                  &
                     ComputePoints(i + di - 2 * dj, j + dj - 2 * di, k)   *                  & 
                     ComputePoints(i + di - 3 * dj, j + dj - 3 * di, k)) == 0 )              &
                     NoSouthWest = .true.

            endif


        endif

        if (.not. NoNorthEast) then

            Boundary3Old  = OldField (i +     dj, j +     di, k) 

            Boundary32Old = OldField (i + 2 * dj, j + 2 * di, k) 

            Boundary33Old = OldField (i + 3 * dj, j + 3 * di, k) 

        endif

        if (.not. NoSouthWest) then

            Boundary4Old  = OldField (i -     dj, j -     di, k) 

            Boundary42Old = OldField (i - 2 * dj, j - 2 * di, k) 

            Boundary43Old = OldField (i - 3 * dj, j - 3 * di, k) 

        endif

        Boundary5Old  = OldField (i     , j     , k) 

        if (present (ConstantCelerity)) then

            ConstantCelerity_ = ConstantCelerity

        else

            ConstantCelerity_ = .false.

        endif 

        if (present (Explicit)) then

            Explicit_ = Explicit

        else

            Explicit_ = .false.

        endif

        if (present (NormalRadiation)) then

            NormalRadiation_ = NormalRadiation

        else

            NormalRadiation_ = .false.

        endif 

        if (.not. ConstantCelerity_) then

            !The numerical scheme of Marchesiello et al. (2000)

            TimeVariability   = Interior1New - Interior1Old

            SpaceVariabilityX = Interior1New - Interior2New

            AuxIntY = Interior6Old - Interior7Old

            if ((AuxIntY * TimeVariability) > 0 ) then

                SpaceVariabilityY =  Interior1Old - Interior7Old
                
                if (NoSouthWest) SpaceVariabilityY = 0.

            else

                SpaceVariabilityY =  Interior6Old - Interior1Old 

                if (NoNorthEast) SpaceVariabilityY = 0.

            endif


        
            SpaceSquare = SpaceVariabilityX**2. + SpaceVariabilityY**2
        

            if (SpaceSquare > 0.) then

                WaveCelerityX_ = - (TimeVariability * SpaceVariabilityX) / SpaceSquare
                WaveCelerityY_ = - (TimeVariability * SpaceVariabilityY) / SpaceSquare

                if (abs(WaveCelerityX_) > LimitMax) WaveCelerityX_ = LimitMax * WaveCelerityX_/ abs(WaveCelerityX_)
                if (abs(WaveCelerityY_) > LimitMax) WaveCelerityY_ = LimitMax * WaveCelerityY_/ abs(WaveCelerityY_)
            
            else

                WaveCelerityX_ = 0.
                WaveCelerityY_ = 0. 

            endif


        else

            WaveCelerityX_ = InternalCelerity * WaveDirectionX
            WaveCelerityY_ = InternalCelerity * WaveDirectionY

        endif


        if (present(ReferenceField)) then

            ReferenceValue = ReferenceField(i, j, k)

        else

            ReferenceValue = 0.

        endif

        if (present(TrelaxOut)) then

            TrelaxOut_ = TrelaxOut

        else

            TrelaxOut_ = 86400 * 300

        endif



        if (present(TrelaxIn)) then

            TrelaxIn_ = TrelaxIn

        else

            TrelaxIn_ = 86400 * 300

        endif

        if(present(FlowVelX)) WaveCelerityX_ = WaveCelerityX_ + FlowVelX


        !Implicit or explicit calculation of the Baroclinic velocity in a boundary point

        if (WaveCelerityX_ > 0) then

            WaveCelerityX_ = 4 * WaveCelerityX_
        
            if (Explicit_) then
                
                !AdjacentPropX    =   Interior1Old

                AdjacentPropX    =   0.0546875 * Interior3Old - 0.2578125 * Interior2Old +    &
                                    0.6015625 * Interior1Old + 0.6015625 * OldField (i, j, k)

            else

               ! AdjacentPropX    =   Interior1New
                AdjacentPropX    =   0.0546875 * Interior3New - 0.2578125 * Interior2New +    &
                                    0.6015625 * Interior1New 

            endif

            Trelax          =   TrelaxOut_


        else

            AdjacentPropX     =   ReferenceValue
            WaveCelerityX_   =   0.
            WaveCelerityY_   =   0.
            Trelax           =   TrelaxIn_


        endif


        AuxInt = WaveCelerityX_ * AdjacentPropX

        if (.not. NormalRadiation_) then

            if      (WaveCelerityY_ >= 0 .and. .not. NoSouthWest) then
                
        
                !AuxInt         = AuxInt - WaveCelerityY_ * (Boundary5Old - Boundary4Old)

                AdjacentPropY  =   0.0546875 * Boundary43Old - 0.2578125 * Boundary42Old +    &
                                   0.6015625 * Boundary4Old  + 0.6015625 * Boundary5Old
                AuxInt         =   AuxInt - 4 * WaveCelerityY_ * (Boundary5Old - AdjacentPropY)

            else if (WaveCelerityY_ <  0 .and. .not. NoNorthEast) then

!                 AuxInt         = AuxInt - WaveCelerityY_ * (Boundary3Old - Boundary5Old)

                AdjacentPropY  =   0.0546875 * Boundary33Old - 0.2578125 * Boundary32Old +    &
                                   0.6015625 * Boundary3Old  + 0.6015625 * Boundary5Old
                AuxInt         =   AuxInt - 4 * WaveCelerityY_ * (AdjacentPropY - Boundary5Old)

            endif

        endif


        OldField (i, j, k)=  NewField (i, j, k)

        AuxBound = 1

        if (Explicit_) then

            AuxInt   = AuxInt  - (WaveCelerityX_ + DT / Trelax) * OldField (i, j, k)

        else
        
            !AuxBound = AuxBound + WaveCelerityX_  + DT / Trelax
            AuxBound = AuxBound + WaveCelerityX_ * (1 - 0.6015625) + DT / Trelax


        endif


        NewField (i, j, k)= (OldField (i, j, k) +  AuxInt + &
                             ReferenceValue * DT / Trelax) /AuxBound 


        if (present(WaveCelerityX)) WaveCelerityX = WaveCelerityX_
        if (present(WaveCelerityY)) WaveCelerityY = WaveCelerityY_

        if (present(NewValue     )) NewValue      = NewField (i, j, k)

    end subroutine OrlanskiCelerity2D

    !--------------------------------------------------------------------------



    !C02 Parcial Pressure in the water
    !
    !With the values of dissolved CO2 concentration (mg/l), temperature (oC), Salinity (ppt), and Pressure (atm)
    !this function calculates the partial pressure of CO2 in the water (uatm)
    !
    ! M Mateus by Nov2009
    
    real function CO2PartialPressure(CO2, Temperature, Salinity, Pressure)   !uatm
    
        !Arguments-------------------------------------------------------------

        real, intent(IN) :: Temperature     !oC
        real, intent(IN) :: Salinity
        real, intent(IN) :: Pressure        !atm
        real, intent(IN) :: CO2             !mg/l

        !Local-----------------------------------------------------------------

        real :: TKelvin, Ff, Fp, K0
        real :: CO2mass

        !----------------------------------------------------------------------

        TKelvin = Temperature + 273.15
        
        CO2mass = (CO2 / 44 / 1025)               !mg/l to mol/kg    water density = 1025 kg m-3
        
        Ff = -1636.75 + 12.0408 * TKelvin - 0.0327957 * (TKelvin**2) + 3.16528 * (1E-5) * (TKelvin**3)
        
        Fp = EXP(((Ff + 2 *(57.7-0.118 * TKelvin)) * Pressure) / (82.05675 * TKelvin))
    
        K0 = CO2_K0(Temperature, Salinity)
    
        CO2PartialPressure = (CO2mass / (Fp * K0)) * 1E6
    
    end function CO2PartialPressure

    !--------------------------------------------------------------------------

    
    real function CO2_K0(Temperature, Salinity)   !mol Kg-1 atm-1
        
        !Arguments-------------------------------------------------------------
        real, intent(IN) :: Temperature     !oC
        real, intent(IN) :: Salinity
        
        !Local-----------------------------------------------------------------
        real :: TKelvin
        
        !----------------------------------------------------------------------

        TKelvin = Temperature + 273.15
        
            CO2_K0 = EXP(-60.2409 + 93.4517 * 100 / TKelvin + 23.3585 * LOG(TKelvin / 100) + Salinity * &
                     (0.023517 - 0.023656 * TKelvin / 100 + 0.0047036 * (TKelvin / 100)**2))
        
       end function CO2_K0



    !Saturation Oxygen concentration
    !
    !With the values of temperature (oC) and Salinity (ppt), this function calculates the 
    !saturation Oxygen concentration.
    !
    !Based in Portela's thesis (1996), following APHA (1992).

    real function OxygenSaturation(Temperature, Salinity)

        !Arguments-------------------------------------------------------------

        real, intent(IN) :: Temperature     !oC
        real, intent(IN) :: Salinity        !ppt

        !Local-----------------------------------------------------------------
        
        real :: TKelvin
        logical :: OutOfRange

        !----------------------------------------------------------------------
        
        OutOfRange = .false.

        if (Temperature < 0 .or. Temperature > 50 .or. Salinity < 0 .or. Salinity > 80) then
            OutOfRange = .true.
        endif
        
        if (OutOfRange) then
            !mgO2/L
            OxygenSaturation = 6.
        
        else

            TKelvin = Temperature + AbsoluteZero

            OxygenSaturation = exp(-139.34411 + 1.575701E05 / TKelvin                                  &
                                              - 6.642308E07 /(TKelvin * TKelvin                    )   &
                                              + 1.243800E10 /(TKelvin * TKelvin * TKelvin          )   &
                                              - 8.621949E11 /(TKelvin * TKelvin * TKelvin * TKelvin)   &
                                   - (Salinity / 1.80655) * (3.1929E-02 - 19.428 / TKelvin             &
                                                             + 3.8673E03/ (TKelvin * TKelvin)))
                                                             
        endif

    end function OxygenSaturation

    !--------------------------------------------------------------------------

    
    !Saturation Oxygen concentration [mol/lw]
    !
    !With the values of temperature (ºC), this function calculates the 
    !saturation Oxygen concentration.
    !
    !Based in Henry's law (Metcalf & Eddy) Pedro Galvão 2002

    real function OxygenSaturationHenry (Temperature)

        !Arguments-------------------------------------------------------------

        real, intent(IN) :: Temperature     !ºC
        !Local-----------------------------------------------------------------

        real :: Henry
        real :: m
        real :: b
        logical :: OutOfRange

        !----------------------------------------------------------------------
        
        OutOfRange = .false.

        if (Temperature < 0 .or. Temperature > 50) then
            OutOfRange = .true.
        endif
        
        if (OutOfRange) then
            !mgO2/L
            OxygenSaturationHenry = 6.
        
        else
        
            !----------------------------------------------------------------------

            if      (Temperature .GE. 0 .AND. Temperature .LE. 10   ) then
                m       = (3.27 - 2.55) /10
                b       = 3.27 - 10 * m 

                Henry   = m *  Temperature + b

            else if (Temperature .GE. 10 .AND. Temperature .LE. 20  ) then
                m       = (4.01 - 3.27) /10
                b       = 4.01-20 * m

                Henry   = m *  Temperature + b

            else if (Temperature .GE. 20 .AND. Temperature .LE. 30  ) then
                m       = (4.75 - 4.01) /10
                b       = 4.75 - 30 * m

                Henry   = m *  Temperature + b

            else if (Temperature .GE. 30 .AND. Temperature .LE. 40  ) then
                m       = (5.35 - 4.75) /10
                b       = 5.35 - 40 * m

                Henry   = m *  Temperature + b

            else if (Temperature .GE. 40 .AND. Temperature .LE. 50  ) then
                m       = (5.88 - 5.35) /10
                b       = 5.88 - 50 * m

                Henry   = m *  Temperature + b

            else

            end if
                
            OxygenSaturationHenry = 0.21 / (Henry * 1. * 10. ** 4.) * 55.6
        
        endif            

    end function OxygenSaturationHenry

    !--------------------------------------------------------------------------

    !Saturation Oxygen concentration
    !
    !With the values of temperature (oC) and Salinity (ppt), this function calculates the 
    !saturation Oxygen concentration.
    !PALT Altitude correction
    !SALTWATER True/False  
    !Based CeQualW2 V3.1 (Mortimer 1981)

    real function OxygenSaturationCeQualW2(Temperature, Salinity, PALT)

        !Arguments-------------------------------------------------------------
        real,    intent(IN)             :: Temperature     !oC
        real,    intent(IN)             :: Salinity        !ppt
        real,    intent(IN)             :: PALT 

        !Local-----------------------------------------------------------------
        real                            :: Aux
        logical :: OutOfRange

        !Begin-----------------------------------------------------------------

        OutOfRange = .false.

        if (Temperature < 0 .or. Temperature > 50 .or. Salinity < 0 .or. Salinity > 80) then
            OutOfRange = .true.
        endif
        
        if (OutOfRange) then
            !mgO2/L
            OxygenSaturationCeQualW2 = 6.
        
        else
            Aux = exp(7.7117-1.31403*(log(Temperature+45.93)))*PALT

            if (Salinity.gt.0.5) then  
          
                Aux = exp(log(aux)-Salinity *                       &
                      (1.7674E-2 - 1.0754E-1/(Temperature+273.15) + &
                      2.1407E3/(Temperature+273.15)**2))                
            endif

            OxygenSaturationCeQualW2 = Aux

        endif
        
    end function OxygenSaturationCeQualW2

    !--------------------------------------------------------------------------

    real function Sigma(Method, PressCorrec, T, S, Depth)

        !Arguments----------------------------------------------------
        integer                             ::  Method
        logical                             ::  PressCorrec
        real                                ::  T, S
        real, optional                      ::  Depth

        select case(Method)

            case(LeendertseState_)
                Sigma = SigmaLeendertse(T, S)

            case(UNESCOState_)
                Sigma = SigmaUNESCO(T, S)

            case(Mel96State_)
                Sigma = SigmaUNESCO(T, S)

            case(JMD95State_)
                Sigma = SigmaUNESCO(T, S)

            case(WangState_)
                Sigma = SigmaWang(T, S)

            case(Linear_)
                Sigma = 1025. -  dble(SigmaDensityReference) + 0.78 * (S - 33.75)

            case default
                Sigma = 1025. -  dble(SigmaDensityReference)

        end select

        if (PressCorrec) then

            select case(Method)

                case(UNESCOState_)
                    Sigma = SigmaUNESCOPressureCorrection(T, S, Depth, Sigma)

                case(Mel96State_)
                    Sigma = SigmaMel96PressureCorrection(T, S, Depth, Sigma)

                case(JMD95State_)
                    Sigma = SigmaJMD95PressureCorrection(T, S, Depth, Sigma)

               case default

            end select

        end if

    end function Sigma

    !--------------------------------------------------------------------------

    real function Density(Method, PressCorrec, T, S, Depth)

        !Arguments----------------------------------------------------
        integer                             ::  Method
        logical                             ::  PressCorrec
        real                                ::  T, S
        real, optional                      ::  Depth

        Density = SigmaDensityReference + Sigma(Method, PressCorrec, T, S, Depth)

    end function Density

    !--------------------------------------------------------------------------

    real function SigmaWang (T, S)
        ! Modified Jan/2006     :Guillaume Riflet
        !  THIS FUNCTION COMPUTES DENSITY-1; FROM WANG, JPO '84, 1191-1199.
        
        !Arguments-------------------------------------------------------------
        real                                                  :: T, S
 
        !Local-----------------------------------------------------------------
        REAL R1,R2,RR3,R4,R5,R6

        DATA R1,R2,RR3,R4,R5,R6/28.152,7.35E-2,4.69E-3,8.02E-1,2.0E-3,35.0/

        SigmaWang = R1 - R2*T-RR3*T*T + (R4-R5*T)*(S-R6)

    end function SigmaWang

    !--------------------------------------------------------------------------

    real function SigmaUNESCO (T, S)

        ! Modified jul/2004     :João Nogueira
        ! Modified Mar/2005     :Guillaume Riflet
        !                       Now it gives a sigma density variable (1000 kg/m3)
        
        !Arguments-------------------------------------------------------------
        real                                              :: T, S

        !Local-----------------------------------------------------------------
        real                                              :: x, T1, S1, T2, T3, T4, T5
        real                                              :: S2, S3
        real(8)                                           :: Aux
        
        !Begin-----------------------------------------------------------------

        T1 = T
        S1 = S

        T2 = T1*T1 
        T3 = T1*T2 
        T4 = T2*T2 
        T5 = T1*T4
        S2 = S1*S1 
        S3 = S1*S2 

!       Calculate density for atmosfere pressure:
!
!
!       RO_ST0 = 999.842594+.06793952*T-.00909529*T**2+1.001685E-4*T**3-1.120083E-6*T**4
!                + 6.536332E-9*T**5  +(.824493-.0040899*T +7.6438E-5*T**2-8.2467E-7*T**3 
!                + 5.3875E-9*T4)*S+ (-.00572466+1.0227E-4*T-1.6546E-6*T2).S**1.5
!                +4.8314E-4.S**2
!
!
        Aux = 999.842594 - dble(SigmaDensityReference)
        x= real(Aux)+6.793952e-02*T1-9.09529e-03*T2+1.001685e-04*T3 
        x=x-1.120083e-06*T4+6.536332e-09*T5 
        x=x+S1*(0.824493-4.0899e-03*T1+7.6438e-05*T2-8.2467e-07*T3) 
        x=x+S1*5.3875e-09*T4 
        x=x+sqrt(S3)*(-5.72466e-03+1.0227e-04*T1-1.6546e-06*T2) 
        x=x+4.8314e-04*S2

        SigmaUNESCO = x

    end function SigmaUNESCO                   

    !--------------------------------------------------------------------------
                 
    real function SigmaLeendertse (T, S)

        !Arguments-------------------------------------------------------------
        real                                        :: T, S

        !Local-----------------------------------------------------------------

        SigmaLeendertse = 1000.0 * ((5890.0 + 38.0 * T - 0.375 * T * T + 3.0 * S)       &
                          / ((1779.5 + 11.25 * T - 0.0745 * T * T) -(3.8 + 0.01 * T) * S &
                          + 0.698 * (5890.0 + 38.0 * T - 0.375 * T * T + 3.0 * S)) - 1.)

    end function SigmaLeendertse                   

    !--------------------------------------------------------------------------

    real function SigmaUNESCOPressureCorrection (T, S, Depth, Sigma_P0)

        ! Created jul/2004     :João Nogueira
        !
        ! The correction is made with the UNESCO International Equation of State

        !Arguments-------------------------------------------------------------
        real                                                  :: T, S
        real                                                  :: Depth, Sigma_P0

        !Local-----------------------------------------------------------------
        real                                              :: T1, S1, T2, T3, T4, T5
        real                                              :: S2, S3
        real                                              :: pressureBars
        real                                              :: RO_ST0, RO_STP, K_ST0, K_STP
        real                                              :: CoefA, CoefB, x1, x2, x3
        !Begin-----------------------------------------------------------------

        T1 = T
        S1 = S

        T2 = T1*T1 
        T3 = T1*T2 
        T4 = T2*T2 
        T5 = T1*T4
        S2 = S1*S1 
        S3 = S1*S2 

        RO_ST0 = Sigma_P0


!!      Include pressure effect (K is the bulk modulus; pressureBars is the pressure in Bars)
!!
!!      RO_STP=RO_ST0/(1-p/K_STP)
!!
!!      K_STP=K_ST0+pressureBars*(CoefA+pressureBars*CoefB)
!!
!!
!!      Calculate pressure in Bars


            pressureBars = 0.1*Depth


!!      Calculate bulk modulus for atmosfere pressure      
!!
!!
!!      K_ST0 = 19652.21+148.4206*T-2.327105*T**2+.01360477*T**3-5.155288E-5*T**4
!!              +(54.6746-.603459*T+.0109987*T**2-6.167E-5*T**3)*S+(.07944+.016483*T
!!              -5.3009E-4*T**2)*S**1.5


            x1 = 19652.21+148.4206*T1-2.327105*T2+.01360477*T3-5.155288E-5*T4
            x1 = x1+(54.6746-.603459*T1+.0109987*T2-6.167E-5*T3)*S1
            x1 = x1+(.07944+.016483*T1-5.3009E-4*T2)*sqrt(S3)
            K_ST0 = x1



!!      Calculate pressure coeficients CoefA and CoefB
!!
!!
!!      CoefA = 3.239908+.00143713*T+1.16092E-4*T**2-5.77905E-7*T**3+(.0022838
!!              -1.0981E-5*T-1.6078E-6*T**2)*S + 1.91075E-4*S**1.5
!!
!!      CoefB = 8.50935E-5-6.12293E-6*T+5.2787E-8.T**2+(-9.9348E-7
!!              +2.0816E-8*T+9.1697E-10*T**2)*S


            x2 = 3.239908+.00143713*T1+1.16092E-4*T2-5.77905E-7*T3
            x2 = x2+(.0022838-1.0981E-5*T1-1.6078E-6*T2)*S1 + 1.91075E-4*sqrt(S3)
            CoefA = x2
 
!
            x3 = 8.50935E-5-6.12293E-6*T1+5.2787E-8*T2
            x3 = x3+(-9.9348E-7+2.0816E-8*T1+9.1697E-10*T2)*S1
            CoefB = x3
 
 
!!      Calculate bulk modulus
 
            K_STP=K_ST0+pressureBars*(CoefA+pressureBars*CoefB)
 
 
!!      Calculate corrected density
 
!            RO_STP=RO_ST0/(1.0-pressureBars/K_STP)
            RO_STP= (SigmaDensityReference*pressureBars + RO_ST0*K_STP) /(K_STP-pressureBars)
 
            SigmaUNESCOPressureCorrection = RO_STP
  

    end function SigmaUNESCOPressureCorrection
    
    !--------------------------------------------------------------------------
 
    real function SigmaJMD95PressureCorrection (T, S, Depth, Sigma_P0)

        ! Created jul/2004     :João Nogueira
        !
        ! The correction is made with the UNESCO International Equation of State
        
        !Arguments-------------------------------------------------------------
        real                                                  :: T, S
        real                                                  :: Depth, Sigma_P0

        !Local-----------------------------------------------------------------
        real                                              :: T1, S1, T2, T3, T4, T5
        real                                              :: S2, S3
        real                                              :: pressureBars
        real                                              :: RO_ST0, RO_STP, K_ST0, K_STP
        real                                              :: CoefA, CoefB, x1, x2, x3
       
        !Begin-----------------------------------------------------------------

        T1 = T
        S1 = S

        T2 = T1*T1 
        T3 = T1*T2 
        T4 = T2*T2 
        T5 = T1*T4
        S2 = S1*S1 
        S3 = S1*S2 

        RO_ST0 = Sigma_P0

!!      Calculate pressure in Bars
        pressureBars = 0.1*Depth

        x1 = 19659.33 + 144.4304*T1 - 1.706103*T2 + .009648704*T3 - 4.190253E-5*T4
        x1 = x1+( 52.84855 - .3101089*T1 + .006283263*T2 - 5.084188E-5*T3)*S1
        x1 = x1+( .3886640 + .009085835*T1 - 4.619924E-4*T2)*sqrt(S3)
        K_ST0 = x1

        x2 = 3.186519 + .02212276*T1 - 2.984642E-4*T2 + 1.9564155E-6*T3
        x2 = x2+(.006704388 - 1.847318E-4*T1 + 2.059331E-7*T2)*S1 + 1.480266E-4*sqrt(S3)
        CoefA = x2

        x3 = 2.102898E-4 - 1.202016E-5*T1 + 1.394680E-7*T2
        x3 = x3+( - 2.040237E-6 + 6.128773E-8*T1 + 6.207323E-10*T2)*S1
        CoefB = x3
 
        K_STP=K_ST0+pressureBars*(CoefA+pressureBars*CoefB)
 
        RO_STP= (SigmaDensityReference*pressureBars + RO_ST0*K_STP) /(K_STP-pressureBars)
 
        SigmaJMD95PressureCorrection = RO_STP
  

    end function SigmaJMD95PressureCorrection
    
    !--------------------------------------------------------------------------

    real function SigmaMel96PressureCorrection(T, S, Depth, Sigma_P0)
    ! Created Mar/2005     :Guillaume Riflet
    !
    !The correction is made in decibars 
    !
    !Mellor 96 formula.
    !
    !T is the POTENTIAL TEMPERATURE!
        
     !Arguments-------------------------------------------------------------
      real                                                  :: T, S
      real                                                  :: Depth, Sigma_P0

     !Local-----------------------------------------------------------------
      real                                              :: p, T2, c, cc, K

        T2 = T*T

!       Calculate pressure in db (decibars)
        p = Depth

        c = 1449.1 + .00821*p + 4.55*T - .045*T2 + 1.34*(S-35.);
        cc = c*c;
        cc = 1./cc;
        K = 1.e4*p*cc*(1.-.20*p*cc); 
        SigmaMel96PressureCorrection = Sigma_P0 + K;

    end function SigmaMel96PressureCorrection
    
    !--------------------------------------------------------------------------

    real function ConvertTemperature(T,S,Depth)
    !
    !Created: Mar/2005
    !Author: Guillaume Riflet
    !
    !Function: converts in-situ temperature into potential temperature
    !Reference: G.Mellor 96 and Bryden73
    !
        !Arguments-------------------------------------------------------------
        real, intent (in)       :: T,S,Depth     !in-situ temperature, salinity and depth
    
        !Local-----------------------------------------------------------------
        real                    :: p

        p = Depth

        ConvertTemperature = T  + p*( -1*(3.6504e-5 + T*(8.3198e-6 - T*(5.4065e-8 + 4.0274e-10*T))) &
            - (1.7439e-6 - 2.9778e-8*T)*(S-35.)                 &
            + p*(4.1057e-11*(S-35.)                             &
            - (8.9309e-9 - T*(3.1628e-10 + 2.1987e-12*T))       & 
            - (-1.6056e-13 + 5.0484e-15*T)*p) );
    
    end function ConvertTemperature
    
    !--------------------------------------------------------------------------

    real function SpecificHeatUNESCO(s, t, Depth)
    !Calculates the specific heat of water
    !Created: Feb'05    Guillaume Riflet
    
    !ref: millero et al,1973,jgr,78,4499-4507
    !      millero et al, unesco report no. 38 1981 pp. 99-188.
    !UNESCO formulas: reference pressure p = 0 at sealevel
        
        !Arguments-------------------------------------------------------------
        real, intent (in)       :: s, t, Depth     !Salinity, temperature and pressure in decibars or depth

        !Local-----------------------------------------------------------------
        real                    :: a, b, c
        real                    :: cp0, cp1, cp2
        real                    :: p, sr

        !Scale pressure in bars from depth
        p = Depth*0.100650603
        sr = sqrt(abs(s))

        !specific heat for p=0
        a = (-1.38385e-3*t + 0.1072763)*t - 7.643575
        b = (5.148e-5*t - 4.07718e-3)*t + 0.1770383
        c = (((2.093236e-5*t - 2.654387e-3)*t + 0.1412855)*t - 3.720283)*t + 4217.4
        cp0 = (b*sr + a)*s + c

        !specific heat for s=0
        a = (((1.7168e-8*t + 2.0357e-6)*t - 3.13885e-4)*t + 1.45747e-2)*t - 0.49592
        b = (((2.2956e-11*t-4.0027e-9)*t + 2.87533e-7)*t - 1.08645e-5)*t + 2.4931e-4
        c = ((6.136e-13*t - 6.5637e-11)*t + 2.6380e-9)*t - 5.422e-8
        cp1 = ((c*p + b)*p + a)*p

        !specific heat for s>0
        a = (((-2.9179e-10*t + 2.5941e-8)*t + 9.802e-7)*t - 1.28315e-4)*t + 4.9247e-3
        b = (3.122e-8*t - 1.517e-6)*t - 1.2331e-4
        a = (a + b*sr)*s
        b = ((1.8448e-11*t-2.3905e-9)*t+1.17054e-7)*t-2.9558e-6
        b = (b+9.971e-8*sr)*s
        c = (3.513e-13*t-1.7682e-11)*t+5.540e-10
        c = (c-1.4300e-12*t*sr)*s
        cp2 = ((c*p + b)*p + a)*p

        !specific heat return
        SpecificHeatUNESCO = cp0 + cp1 + cp2

    end function SpecificHeatUNESCO

!--------------------------------------------------------------------------
!---------------------------------------------------------------

real function SAL78(CND,T,P)

! ***********************************************
! SAL78 converts conductivity to salinity
! the conductivity ratio (CND)=1.0000000 for salinity=35
! PSS-78. temperature=15.0 deg. celcius, and atmospheric
! pressure.
! CND = Condutivity / Condutivity(S=35,T=15,P=0)
! The CTD Sea-Birds of  use for salinity=35
! PSS-78. temperature=15.0 deg. celcius, and atmospheric
! pressure a value of 42.914 mmho/cm (mmho/cm == 1000 microS/cm)
! http://www.seabird.com/application_notes/AN14.htm
! ***********************************************
! references: also located in UNESCO Report NO. 37 1981
! Practical Salinity Scale 1978: E.L. Lewis IEEE Ocean Eng.
! Jan. 1980
! ***********************************************
! units:
!     pressure P decibars
!     temperature T deg. celcius (IPTS-68)
!     conductivity CND ratio (M=0)
!     salinity SAL78 (PSS-78) (M=0)
! ***********************************************
! SAL78 ratio: returns zero for conductivity ratio: <0.0005
! SAL78: returns zero for salinity: <0.02
! ***********************************************
! Practical Salinity Scale 1978 definition with temperature
! correction
! ***********************************************
! convert conductivity to salinity

!Arguments----------------------------------------
real                :: CND,T,P
!Local -------------------------------------------
real                :: R, DT, RT, RT35, A, B, C
!Begin -------------------------------------------

R=CND

DT=T-15.0

RT35=(((1.0031E-9*T-6.9698E-7)*T+1.104259E-4)*T +2.00564E-2)*T+0.6766097

A=-3.107E-3*T+0.4215

B=(4.464E-4*T+3.426E-2)*T+1.0

C=((3.989E-15*P-6.370E-10)*P+2.070E-5)*P

RT=R/(RT35*(1.0+C/(B+A*R)))

RT=sqrt(ABS(RT))

SAL78=((((2.7081*RT-7.0261)*RT+14.0941)*RT+25.3851)*RT -0.1692)*                        &
          RT+0.0080 +(DT/(1.0+0.0162*DT))*(((((-0.0144*RT+ 0.0636)*                     &
          RT-0.0375)*RT-0.0066)*RT-0.0056)*RT+0.0005)
End function

!-------------------------------------------------------------------

real function depth(P,lat)

! ****************************
! depth in meters from pressure in decibars using
! Saunder's and Fofonoff's method.
! deep-sea res., 1976,23,109-111.
! formula refitted for 1980 equation of state
! units:
!     pressure p decibars
!     latitude lat degrees
!     depth depth meters
!     checkvalue: depth=9712.653 M for P=10000 decibars,
!     latitude=30 deg
!     above for standard ocean: T=0 deg. celcius;
!     S=35 (PSS-78)
! ****************************

!Arguments----------------------------------------
real                :: P, lat
!Local -------------------------------------------
real                :: x, gr
!Begin -------------------------------------------

x=sin(lat/57.29578)

x=x*x

! gr=gravity variation with latitude: Anon (1970)
! Bulletin Geodesique

gr=9.780318*(1.0+(5.2788E-3+2.36E-5*x)*x)+1.092E-6*p

depth=(((-1.82E-15*P+2.279E-10)*P-2.2512E-5)*P+9.72659)*P

depth=depth/gr

end function

real function SVEL(S,T,PO)

! ***********************************************
! Sound Speed Seawater Chen and Millero 1977,JASA,62,1129-35
! units:
!     pressure PO decibars
!     temperature T deg celcius (IPTS-68)
!     salinity S (PSS-78)
!     sound speed SVEL meters/second
!     checkvalue:SVEL=1731,995m/s,S=40(PSS-78),T=40degC,P=10000db
! ***********************************************
! scale pressure to bars

!Arguments----------------------------------------
real                :: S,T,PO
!Local -------------------------------------------
real                :: SR,P, D, B1, B0, B4, A0, A1, A2, A3, A4, &
                       C0, C1, C2, C3, C4
!Begin -------------------------------------------

P=PO/10.

SR=sqrt(ABS(S))

! S**2 term

D=1.727E-3-7.9836E-6*P

! S**3/2 term

B1=7.3637E-5+1.7945E-7*T

B0=-1.922E-2-4.42E-5*T

B4=B0+B1*P

! S**1 term

A3=(-3.389E-13*T+6.649E-12)*T+1.100E-10

A2=((7.988E-12*T-1.6002E-10)*T+9.1041E-9)*T-3.9064E-7

A1=(((-2.0122E-10*T+1.0507E-8)*T-6.4885E-8)*T-1.2580E-5)*T+9.4742E-5

A0=(((-3.21E-8*T+2.006E-6)*T+7.164E-5)*T-1.262E-2) *T+1.389

A4=((A3*P+A2)*P+A1)*P+A0

! s**0 term

C3=(-2.3643E-12*T+3.8504E-10)*T-9.7729E-9

C2=(((1.0405E-12*T-2.5335E-10)*T+2.5974E-8)*T-1.7107E-6) *T+3.1260E-5

C1=(((-6.1185E-10*T+1.3621E-7)*T-8.1788E-6)*T+6.8982E-4) *T+0.153563

C0=((((3.1464E-9*T-1.47800E-6)*T+3.3420E-4)*T-5.80852E-2)*T+5.03711)*T+1402.388

C4=((C3*P+C2)*P+C1)*P+C0

! sound speed return

SVEL=C4+(A4+B4*SR+D*S)*S

end function

!----------------------------------------------------------------------------------------------------------- 
   !**********************************************************************
    !                                                                     *
    ! RODAXY - Conversao das coordenadas (x,y) de um ponto num referencial*
    !          cartesiano para outro que se encontra rodado e translacio- *
    !          nado em relacao ao primeiro                                *
    !                                                                     *
    !---------------------------------------------------------------------*
    !                                                                     *
    !   Ultima alteracao:  95-11-24                   Flavio              *
    !                      May99 - Code alignment     Frank               *
    !                                                                     *
    !**********************************************************************      
    subroutine RODAXY (XORIG, YORIG, ALPHA, X_PONTO, Y_PONTO)

        !Arguments-------------------------------------------------------------
        real, intent (in)   :: XORIG, YORIG, ALPHA
        real, intent (inout)  :: X_PONTO, Y_PONTO         

        !Local-----------------------------------------------------------------
        real(8), parameter  :: PI_DBLE       = 3.1415926536 !PI
        real(8)             :: Radianos, A11, A12, A21, A22, X_TEMP, Y_TEMP

        !----------------------------------------------------------------------


        !Coordenadas do ponto no referencial rodado
        if (ALPHA == 0.) then
            X_PONTO = X_PONTO + XORIG
            Y_PONTO = Y_PONTO + YORIG

        else  

            !Calcula os cosenos directores
            Radianos = ALPHA * PI_DBLE / 180.

            A11 = COS(Radianos)
            A12 =-COS(PI_DBLE/2. - Radianos)
            A21 = SIN(Radianos)
            A22 = SIN(PI_DBLE/2. - Radianos)

            X_TEMP = dble(XORIG) + A11 * dble(X_PONTO) + A12 * dble(Y_PONTO)
            Y_TEMP = dble(YORIG) + A21 * dble(X_PONTO) + A22 * dble(Y_PONTO)

            X_PONTO= real(X_TEMP)
            Y_PONTO= real(Y_TEMP)

        endif

    end subroutine RodaXY

    !--------------------------------------------------------------------------

    subroutine FromCartesianToGrid (Xcart, Ycart, Tetha1, Tetha2, Xgrid, Ygrid)

        !Arguments-------------------------------------------------------------
        real, intent (in)     :: Xcart, Ycart, Tetha1, Tetha2
        real, intent (out)    :: Xgrid, Ygrid

        !Local-----------------------------------------------------------------
        real                  :: a, b, c, d, e, f, g

        a = Xcart
        b = Ycart
        c = cos(Tetha1)
        d = cos(Tetha2)
        e = sin(Tetha1)
        f = sin(Tetha2)

        g = (b - a*e/c) / (f - d*e/c)

        !X grid component
        Xgrid = (a - d*g) / c

        !Y grid component
        Ygrid = g

    end subroutine FromCartesianToGrid

    !--------------------------------------------------------------------------

    subroutine FromGridToCartesian (Xgrid, Ygrid, Tetha1, Tetha2, Xcart, Ycart)

        !Arguments-------------------------------------------------------------
        real, intent (in)     :: Xgrid, Ygrid, Tetha1, Tetha2
        real, intent (out)    :: Xcart, Ycart

        !Local-----------------------------------------------------------------

        Xcart = Xgrid * cos(Tetha1) + Ygrid * cos(Tetha2)
        Ycart = Xgrid * sin(Tetha1) + Ygrid * sin(Tetha2)

    end subroutine FromGridToCartesian

    !--------------------------------------------------------------------------


    subroutine InterpolateValueInTime(ActualTime, Time1, Value1, Time2, Value2, ValueOut)

        !Arguments-------------------------------------------------------------
        type(T_Time),      intent(IN)               :: ActualTime
        type(T_Time),      intent(IN)               :: Time1
        real,              intent(IN)               :: Value1
        type(T_Time),      intent(IN)               :: Time2
        real,              intent(IN)               :: Value2
        
        real,              intent(OUT)              :: ValueOUT

        !Local-----------------------------------------------------------------
        real                                        :: DT1, DT2, DTtotal


        DT1      = ActualTime - Time1
        DT2      = Time2      - ActualTime
        DTtotal  = DT1 + DT2

        ValueOUT = (DT1 * Value2 + DT2 * Value1) / DTtotal

    end subroutine InterpolateValueInTime

    !--------------------------------------------------------------------------

    subroutine InterpolateMatrix2DInTime(ActualTime, Size, Time1, Matrix1, &
                                         Time2, Matrix2, MatrixOUT, PointsToFill2D)

        !Arguments-------------------------------------------------------------
        type(T_Time),      intent(IN)                   :: ActualTime
        type(T_Size2D)                                  :: Size
        type(T_Time),      intent(IN)                   :: Time1
        real, dimension(:,:), pointer                   :: Matrix1
        type(T_Time),      intent(IN)                   :: Time2
        real, dimension(:,:), pointer                   :: Matrix2
        real, dimension(:,:), pointer                   :: MatrixOUT
        integer, dimension(:, :), pointer, optional     :: PointsToFill2D

        !Local-----------------------------------------------------------------
        real                                            :: X1, X, X2

        !Begin-----------------------------------------------------------------

        !Time1 = 0
        X1       = 0.
        X        = ActualTime - Time1
        X2       = Time2      - Time1

        call InterpolateLinearyMatrix2D(X, Size, X1, Matrix1,                           &
                                        X2, Matrix2, MatrixOUT, PointsToFill2D)

    end subroutine InterpolateMatrix2DInTime

    !--------------------------------------------------------------------------

    subroutine InterpolateMatrix3DInTime(ActualTime, Size, Time1, Matrix1, &
                                         Time2, Matrix2, MatrixOUT, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        type(T_Time),      intent(IN)                   :: ActualTime
        type(T_Size3D)                                  :: Size
        type(T_Time),      intent(IN)                   :: Time1
        real, dimension(:,:,:), pointer                 :: Matrix1
        type(T_Time),      intent(IN)                   :: Time2
        real, dimension(:,:,:), pointer                 :: Matrix2
        real, dimension(:,:,:), pointer                 :: MatrixOUT
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local-----------------------------------------------------------------
        real                                            :: X1, X, X2
        
        !Begin-----------------------------------------------------------------

        !Time1 = 0
        X1       = 0.
        X        = ActualTime - Time1
        X2       = Time2 - Time1

        call InterpolateLinearyMatrix3D(X, Size, X1, Matrix1,                           &
                                        X2, Matrix2, MatrixOUT, PointsToFill3D)

    end subroutine InterpolateMatrix3DInTime

    !--------------------------------------------------------------------------
    
    real function LinearInterpolation (x1, y1, x2, y2, xc)

        !Arguments-------------------------------------------------------------
        real                                        :: x1,y1,x2,y2,xc
        !Local-----------------------------------------------------------------
        !real                                        :: dist1, dist2
        !Begin-----------------------------------------------------------------

!        if (x2 > x1) Then
!            dist1 = xc - x1
!            dist2 = x2 - xc
!        else
!            dist1 = xc - x2
!            dist2 = x1 - xc
!        end If
!        LinearInterpolation = (dist2 * y1 + dist1 * y2) / (dist1 + dist2)
        
        LinearInterpolation = y1 + (xc - x1) / (x2 - x1) * (y2 - y1)

    end function
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine InterpolateLinearyMatrix2D(X, Size, X1, Matrix1, &
                                         X2, Matrix2, MatrixOUT, PointsToFill2D)

        !Arguments-------------------------------------------------------------
        real                                            :: X1, X2, X
        type(T_Size2D)                                  :: Size
        real, dimension(:,:), pointer                   :: Matrix1
        real, dimension(:,:), pointer                   :: Matrix2
        real, dimension(:,:), pointer                   :: MatrixOUT
        integer, dimension(:, :), pointer, optional     :: PointsToFill2D

        !Local-----------------------------------------------------------------
        integer                                         :: i, j
        real                                            :: DT1, DT2, DTtotal
        integer                                         :: CHUNK 

        !Begin-----------------------------------------------------------------

        DT1      = X - X1
        DT2      = X2 - X
        DTtotal  = DT1 + DT2

        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "InterpolateLinearyMatrix2D")
       
        CHUNK = CHUNK_J(Size%JLB, Size%JUB)

        if(present(PointsToFill2D))then
            
            !!$OMP PARALLEL SHARED(CHUNK, DT1, Matrix2, DT2, Matrix1, DTtotal, PointsToFill2D) PRIVATE(I,J)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%ILB, Size%JUB
            do i = Size%ILB, Size%IUB
            
                if(PointsToFill2D(i,j) == 1)then

                    MatrixOUT(i,j) = (DT1 * Matrix2(i,j) + DT2 * Matrix1(i,j)) / DTtotal
            
                endif

            enddo
            enddo
            !!$OMP END DO NOWAIT
            !!$OMP END PARALLEL

        else
            
            !!$OMP PARALLEL SHARED(CHUNK, DT1, Matrix2, DT2, Matrix1, DTtotal) PRIVATE(I,J)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do j = Size%ILB, Size%JUB
            do i = Size%ILB, Size%IUB
            
                MatrixOUT(i,j) = (DT1 * Matrix2(i,j) + DT2 * Matrix1(i,j)) / DTtotal

            enddo
            enddo
            !!$OMP END DO NOWAIT
            !!$OMP END PARALLEL

        endif

        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "InterpolateLinearyMatrix2D")

    end subroutine InterpolateLinearyMatrix2D

    !--------------------------------------------------------------------------

    subroutine InterpolateLinearyMatrix3D(X, Size, X1, Matrix1, &
                                         X2, Matrix2, MatrixOUT, PointsToFill3D)

        !Arguments-------------------------------------------------------------
        real                                            :: X1, X2, X
        type(T_Size3D)                                  :: Size
        real, dimension(:,:,:), pointer                 :: Matrix1
        real, dimension(:,:,:), pointer                 :: Matrix2
        real, dimension(:,:,:), pointer                 :: MatrixOUT
        integer, dimension(:, :, :), pointer, optional  :: PointsToFill3D

        !Local-----------------------------------------------------------------
        integer                                         :: i, j, k
        real                                            :: DT1, DT2, DTtotal
        integer                                         :: CHUNK 
        
        !Begin-----------------------------------------------------------------


        DT1      = X - X1
        DT2      = X2 - X
        DTtotal  = DT1 + DT2

        if (MonitorPerformance) call StartWatch ("ModuleFunctions", "InterpolateLinearyMatrix3D")

        CHUNK = CHUNK_K(Size%KLB, Size%KUB)

        if(present(PointsToFill3D))then

            !!$OMP PARALLEL SHARED(CHUNK, DT1, Matrix2, DT2, Matrix1, DTtotal, PointsToFill3D) PRIVATE(I,J)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%ILB, Size%JUB
            do i = Size%ILB, Size%IUB

                if(PointsToFill3D(i,j,k) == 1)then

                    MatrixOUT(i,j,k) = (DT1 * Matrix2(i,j,k) + DT2 * Matrix1(i,j,k)) / DTtotal

                endif

            enddo
            enddo
            enddo
            !!$OMP END DO NOWAIT
            !!$OMP END PARALLEL

        else
            
            !!$OMP PARALLEL SHARED(CHUNK, DT1, Matrix2, DT2, Matrix1, DTtotal, PointsToFill3D) PRIVATE(I,J)
            !!$OMP DO SCHEDULE(DYNAMIC, CHUNK)
            do k = Size%KLB, Size%KUB
            do j = Size%ILB, Size%JUB
            do i = Size%ILB, Size%IUB

                MatrixOUT(i,j,k) = (DT1 * Matrix2(i,j,k) + DT2 * Matrix1(i,j,k)) / DTtotal

            enddo
            enddo
            enddo
            !!$OMP END DO NOWAIT
            !!$OMP END PARALLEL

        end if

        if (MonitorPerformance) call StopWatch ("ModuleFunctions", "InterpolateLinearyMatrix3D")

    end subroutine InterpolateLinearyMatrix3D

    !--------------------------------------------------------------------------

    subroutine ExtraPol2DNearestCell (ILB, IUB, JLB, JUB, KUB, ComputePoints3D, OutValues2D)

        !Arguments-------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB, KUB        
        real,         dimension(:,:),   pointer     :: OutValues2D
        integer,      dimension(:,:,:), pointer     :: ComputePoints3D

        !Local-----------------------------------------------------------------

        integer                                     :: dij, Count, i, j, NumberOfCells
        integer                                     :: jj, ii, dijmax, dimax, djmax
        real                                        :: SumValues
        !$ integer                                  :: CHUNK !
        
        !Begin-----------------------------------------------------------------
            
        !$ CHUNK = CHUNK_J(JLB,JUB) !

        NumberOfCells =  Sum(ComputePoints3D(ILB:IUB, JLB:JUB, KUB))

        if (NumberOfCells > 0) then

            !!!$OMP PARALLEL PRIVATE(j,i,dimax,djmax,dijmax,SumValues,Count,dij,jj,ii)
            !!!$OMP DO SCHEDULE(DYNAMIC,CHUNK)        
            do j = JLB, JUB
            do i = ILB, IUB
            
                if (OutValues2D(i, j) < FillValueReal/4. .and. ComputePoints3D(i, j, KUB) == 1) then
            
                    dimax = IUB-ILB + 1
                    djmax = JUB-JLB + 1

                    dijmax = max(dimax, djmax)
                
                    SumValues   = 0
                    Count = 0

                    do dij=1,dijmax

                        do jj=j-dij,j+dij
                        do ii=i-dij,i+dij

                            if (jj < JLB) cycle
                            if (jj > JUB) cycle
                            if (ii < ILB) cycle
                            if (ii > IUB) cycle

                            if (OutValues2D(ii, jj) > FillValueReal/4.) then
                                SumValues   = SumValues   + OutValues2D(ii, jj) 
                                Count = Count + 1
                            endif

                        enddo
                        enddo

                        if (Count > 0) exit

                    enddo

                    if (Count > 0) then

                        OutValues2D(i, j) = SumValues / real(Count)

                    else
                        stop 'ExtraPol2DNearestCell - ModuleFunctions - ERR10'
                    endif

                endif

            enddo
            enddo
            !!!$OMP END DO NOWAIT
            !!!$OMP END PARALLEL

        endif

    end subroutine ExtraPol2DNearestCell

    !-------------------------------------------------------------------------------------

    subroutine ExtraPol3DNearestCell (ILB, IUB, JLB, JUB, KLB, KUB, ComputePoints3D, OutValues3D)

        !Arguments-------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB, KLB, KUB  
        real,         dimension(:,:,:), pointer     :: OutValues3D
        integer,      dimension(:,:,:), pointer     :: ComputePoints3D

        !Local-----------------------------------------------------------------
        integer                                     :: dij, dk, Count, i, j, k, NumberOfCells
        integer                                     :: jj, ii, kk, dijmax, dimax, djmax
        real                                        :: SumValues
        logical                                     :: NoMapping, OkMap
        
        !$ integer                                  :: CHUNK !
        
        !Begin-----------------------------------------------------------------

        !$ CHUNK = CHUNK_J(JLB,JUB) !
        
        !griflet
        !!!$OMP PARALLEL PRIVATE( k,NumberOfCells,NoMapping,j,i,OkMap,&
        !!!$OMP                   dimax,djmax,dijmax,SumValues,Count, &
        !!!$OMP                   dij,jj,ii,dk,kk)
d1:     do k = KLB, KUB

            if (associated(ComputePoints3D)) then

                NumberOfCells =  Sum(ComputePoints3D(ILB:IUB, JLB:JUB, k))
                NoMapping     = .false.
            
            else
            
                NoMapping     = .true. 
                NumberOfCells = 1
            
            endif 

            if (NumberOfCells > 0) then

                !!!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = JLB, JUB
                do i = ILB, IUB
                
                    if (NoMapping) then
                    
                        OkMap = .true.
                        
                                   
                    else
                    
                        if (ComputePoints3D(i, j, k) == 1) then
                            OkMap = .true.
                        else
                            OkMap = .false.
                        endif    
                    
                    endif

                    if (OutValues3D(i, j, k) < FillValueReal/4. .and. OkMap) then

                        dimax = IUB-ILB + 1
                        djmax = JUB-JLB + 1

                        dijmax = max(dimax, djmax)
                    
                        SumValues   = 0
                        Count = 0

                        do dij=1,dijmax

                            do jj=j-dij,j+dij
                            do ii=i-dij,i+dij

                                if (jj < JLB) cycle
                                if (jj > JUB) cycle
                                if (ii < ILB) cycle
                                if (ii > IUB) cycle

                                if (OutValues3D(ii, jj, k) > FillValueReal/4.) then
                                    SumValues   = SumValues   + OutValues3D(ii, jj, k) 
                                    Count = Count + 1
                                endif

                            enddo
                            enddo

                            if (Count > 0) exit

                        enddo

                        if (Count > 0) then

                            OutValues3D(i, j, k) = SumValues / real(Count)

                        else
                            do dk = 1, KUB - KLB + 1
                            
                                do kk = k - dk, k + dk

                                    if (kk < KLB) cycle
                                    if (kk > KUB) cycle
                                    
                                    if (OutValues3D(i, j, kk) > FillValueReal/4.) then
                                        SumValues   = SumValues   + OutValues3D(i, j, kk) 
                                        Count = Count + 1
                                    endif

                                enddo  
                                
                                if (Count >0) exit
                                
                            enddo 
                                
                            if (Count > 0) then

                                OutValues3D(i, j, k) = SumValues / real(Count)
                            
                            else

                                if (NoMapping) then
                                
                                    OutValues3D(i, j, k) = FillValueReal
                                                                       
                                else

                                    stop 'ExtraPol3DNearestCell - ModuleFunctions - ERR10'
                                
                                endif
                                
                            endif
                        endif

                    endif

                enddo
                enddo
                !!!$OMP END DO NOWAIT

            endif

        enddo d1
        !!!$OMP END PARALLEL

    end subroutine ExtraPol3DNearestCell

    !-------------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------------

    subroutine ExtraPol3DNearestCell_8 (ILB, IUB, JLB, JUB, KLB, KUB, ComputePoints3D, OutValues3D)

        !Arguments-------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB, KLB, KUB  
        real(8),      dimension(:,:,:), pointer     :: OutValues3D
        integer,      dimension(:,:,:), pointer     :: ComputePoints3D

        !Local-----------------------------------------------------------------
        integer                                     :: dij, dk, Count, i, j, k, NumberOfCells
        integer                                     :: jj, ii, kk, dijmax, dimax, djmax
        real                                        :: SumValues
        logical                                     :: NoMapping, OkMap
        
        !$ integer                                  :: CHUNK !
        
        !Begin-----------------------------------------------------------------

        !$ CHUNK = CHUNK_J(JLB,JUB) !
        
        !griflet
        !!!$OMP PARALLEL PRIVATE( k,NumberOfCells,NoMapping,j,i,OkMap,&
        !!!$OMP                   dimax,djmax,dijmax,SumValues,Count, &
        !!!$OMP                   dij,jj,ii,dk,kk)
d1:     do k = KLB, KUB

            if (associated(ComputePoints3D)) then

                NumberOfCells =  Sum(ComputePoints3D(ILB:IUB, JLB:JUB, k))
                NoMapping     = .false.            
            else
            
                NoMapping     = .true. 
                NumberOfCells = 1
            
            endif 

            if (NumberOfCells > 0) then

                !!!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                do j = JLB, JUB
                do i = ILB, IUB
                
                    if (NoMapping) then
                    
                        OkMap = .true.
                                                           
                    else
                    
                        if (ComputePoints3D(i, j, k) == 1) then
                            OkMap = .true.
                        else
                            OkMap = .false.
                        endif    
                    
                    endif

                    if (OutValues3D(i, j, k) < FillValueReal/4. .and. OkMap) then

                        dimax = IUB-ILB + 1
                        djmax = JUB-JLB + 1

                        dijmax = max(dimax, djmax)
                    
                        SumValues   = 0
                        Count = 0

                        do dij=1,dijmax

                            do jj=j-dij,j+dij
                            do ii=i-dij,i+dij

                                if (jj < JLB) cycle
                                if (jj > JUB) cycle
                                if (ii < ILB) cycle
                                if (ii > IUB) cycle

                                if (OutValues3D(ii, jj, k) > FillValueReal/4.) then
                                    SumValues   = SumValues   + OutValues3D(ii, jj, k) 
                                    Count = Count + 1
                                endif

                            enddo
                            enddo

                            if (Count > 0) exit

                        enddo

                        if (Count > 0) then

                            OutValues3D(i, j, k) = SumValues / real(Count)

                        else
                            do dk = 1, KUB - KLB + 1
                            
                                do kk = k - dk, k + dk

                                    if (kk < KLB) cycle
                                    if (kk > KUB) cycle
                                    
                                    if (OutValues3D(i, j, kk) > FillValueReal/4.) then
                                        SumValues   = SumValues   + OutValues3D(i, j, kk) 
                                        Count = Count + 1
                                    endif

                                enddo  
                                
                                if (Count >0) exit
                                
                            enddo 
                                
                            if (Count > 0) then

                                OutValues3D(i, j, k) = SumValues / real(Count)
                            
                            else

                                if (NoMapping) then
                                
                                    OutValues3D(i, j, k) = FillValueReal
                                                                       
                                else

                                    stop 'ExtraPol3DNearestCell_8 - ModuleFunctions - ERR10'
                                
                                endif
                                    
                            endif
                        endif

                    endif

                enddo
                enddo
                !!!$OMP END DO NOWAIT

            endif

        enddo d1
        !!!$OMP END PARALLEL

    end subroutine ExtraPol3DNearestCell_8

    !-------------------------------------------------------------------------------------
    !--------------------------------------------------------------------------

    subroutine FillMatrix2DNearestCell (ILB, IUB, JLB, JUB, ComputePoints2D, OutValues2D)

        !Arguments-------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB
        real,         dimension(:,:),   pointer     :: OutValues2D
        integer,      dimension(:,:),   pointer     :: ComputePoints2D

        !Local-----------------------------------------------------------------

        integer                                     :: dij, Count, i, j, NumberOfCells
        integer                                     :: jj, ii, dijmax, dimax, djmax
        real                                        :: SumValues
       
        !Begin-----------------------------------------------------------------

        NumberOfCells =  Sum(ComputePoints2D(ILB:IUB, JLB:JUB))

i1:     if (NumberOfCells > 0) then
        
d1:         do j = JLB, JUB
            do i = ILB, IUB
            
i2:             if (ComputePoints2D(i, j) == 0) then
            
                    dimax = max(IUB-i + 1, i - ILB + 1)
                    djmax = max(JUB-j + 1, j - JLB + 1)

                    dijmax = max(dimax, djmax)
                
                    SumValues   = 0
                    Count = 0

d2:                 do dij=1,dijmax

                        do jj=j-dij,j+dij
                        do ii=i-dij,i+dij

                            if (jj < JLB) cycle
                            if (jj > JUB) cycle
                            if (ii < ILB) cycle
                            if (ii > IUB) cycle

                            if (ComputePoints2D(ii, jj) == 1) then
                                SumValues   = SumValues   + OutValues2D(ii, jj) 
                                Count = Count + 1
                            endif

                        enddo
                        enddo

                        if (Count > 0) exit

                    enddo d2

                    OutValues2D(i, j) = SumValues / real(Count)

                endif i2

            enddo
            enddo d1
        
        else  i1
        
        !Maintain the original values because there are no nearst cell with values

        endif i1

    end subroutine FillMatrix2DNearestCell

    !-------------------------------------------------------------------------------------

    subroutine FillMatrix3DNearestCell (ILB, IUB, JLB, JUB, KLB, KUB, ComputePoints3D, OutValues3D)

        !Arguments-------------------------------------------------------------
        integer                                     :: JLB, JUB, ILB, IUB, KLB, KUB  
        real,         dimension(:,:,:), pointer     :: OutValues3D
        integer,      dimension(:,:,:), pointer     :: ComputePoints3D

        !Local-----------------------------------------------------------------
        integer                                     :: dij, Count, i, j, k, NumberOfCells
        integer                                     :: dijmax, dimax, djmax
        integer                                     :: imin, imax, jmin, jmax
        real                                        :: SumValues
        
        !Begin-----------------------------------------------------------------

d1:     do k = KLB, KUB

            NumberOfCells =  Sum(ComputePoints3D(ILB:IUB, JLB:JUB, k))


i1:         if (NumberOfCells > 0) then

d2:             do j = JLB, JUB
                do i = ILB, IUB
                    
i2:                 if (ComputePoints3D(i, j, k) == 0) then

                        dimax = max(IUB-i + 1, i - ILB + 1)
                        djmax = max(JUB-j + 1, j - JLB + 1)

                        dijmax = max(dimax, djmax)
                    
                        SumValues = 0
                        Count     = 0
                        
d3:                     do dij=1,dijmax
                            
                            jmin =max(j-dij,JLB)
                            jmax =min(j+dij,JUB)

                            imin =max(i-dij,ILB)
                            imax =min(i+dij,IUB)
                            
                            Count = sum(ComputePoints3D(imin:imax, jmin:jmax, k))
                            
                            if (Count == 0) then
                                cycle
                            else
                                SumValues=sum(ComputePoints3D(imin:imax, jmin:jmax, k) * OutValues3D(imin:imax, jmin:jmax, k)) 
                                OutValues3D(i, j, k) = SumValues / real(Count)
                                exit
                            endif


                        enddo d3

                       
                    endif i2
                    
                enddo
                enddo d2
            
            else  i1

                !Maintain the original values in the layer where are no nearst cell with values
            
            endif i1

        enddo d1

    end subroutine FillMatrix3DNearestCell

    !-------------------------------------------------------------------------------------

   
    !-------------------------------------------------------------------------------------    
    subroutine ReadTimeKeyWords(ObjEnterData, ExtractTime, BeginTime, EndTime, DT,       &
                                VariableDT, ClientModule, MaxDT, GmtReference,           &
                                DTPredictionInterval)

        !Arguments-------------------------------------------------------------
        integer                                     :: ObjEnterData
        integer                                     :: ExtractTime
        type (T_Time)                               :: BeginTime, EndTime
        real                                        :: DT
        logical                                     :: VariableDT
        character(*)                                :: ClientModule
        real, optional                              :: MaxDT
        real, optional                              :: GmtReference
        real, optional                              :: DTPredictionInterval

        !Local-----------------------------------------------------------------
        real(8)                                     :: aux
        real                                        :: MinError
        integer                                     :: STAT_CALL, iflag

        !Reads Begin Time
        call GetData(BeginTime, ObjEnterData, iflag, keyword = 'START',                  &
                     SearchType   = ExtractTime,                                         &
                     ClientModule = ClientModule,                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadTimeKeyWords - ModuleFunctions - ERR01'   

        !Reads End Time 
        call GetData(EndTime,   ObjEnterData, iflag, keyword = 'END',                    &
                     SearchType   = ExtractTime,                                         &
                     ClientModule = ClientModule,                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadTimeKeyWords - ModuleFunctions - ERR02'   

        !Reads DT
        call GetData(DT,   ObjEnterData, iflag, keyword = 'DT',                          &
                     SearchType   = ExtractTime,                                         &
                     ClientModule = ClientModule,                                   &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadTimeKeyWords - ModuleFunctions - ERR03'
        if (DT == 0) then
            write (*,*) 'Time Interval can not be zero'
            write (*,*) 'Module :',ClientModule
            stop 'ReadTimeKeyWords - ModuleFunctions - ERR04'
        endif

        !Reads VariableDT
        call GetData(VariableDT,   ObjEnterData, iflag, keyword = 'VARIABLEDT',          &
                     SearchType   = ExtractTime,                                         &
                     ClientModule = ClientModule,                                   &
                     default      = .false.,                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadTimeKeyWords - ModuleFunctions - ERR05'


        !Reads MaxDT
        if (VariableDT) then

            if (present(MaxDT)) then 
                call GetData(MaxDT,   ObjEnterData, iflag, keyword = 'MAXDT',                    &
                             SearchType   = ExtractTime,                                         &
                             ClientModule = ClientModule,                                   &
                             STAT         = STAT_CALL)
                if (STAT_CALL /= SUCCESS_ .or. iflag == 0) stop 'ReadTimeKeyWords - ModuleFunctions - ERR03a'
                if (MaxDT == 0) then
                    write (*,*) 'Time Interval can not be zero'
                    write (*,*) 'Module :',ClientModule
                    stop 'ReadTimeKeyWords - ModuleFunctions - ERR04a'
                endif
            else
                !MAXDT will be read somewhere else...                
            endif

        endif

        !Reads GmtReference
        if (present(GmtReference)) then 
            call GetData(GmtReference,   ObjEnterData, iflag, keyword = 'GMTREFERENCE',             &
                         SearchType   = ExtractTime,              &
                         ClientModule = ClientModule,             &
                         Default      = 0.,                       &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadTimeKeyWords - ModuleFunctions - ERR22'
            if ( (GmtReference .LT. -12.) .OR. (GmtReference .GT. 12.) ) then
                write (*,*) 'GmtReference cannot exceed 12 modulo'
                write (*,*) 'Module :',ClientModule
                stop 'ReadTimeKeyWords - ModuleFunctions - ERR22a'
            endif
        else
            !GmtReference will be read somewhere else...                
        endif
        
        
        if (present(DTPredictionInterval)) then
            call GetData(DTPredictionInterval, ObjEnterData,      &
                         iflag,                                   &
                         keyword      = 'DT_PREDICTION_INTERVAL', &
                         SearchType   = ExtractTime,              &
                         ClientModule = ClientModule,             &
                         Default      = 60.,                      &
                         STAT         = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadTimeKeyWords - ModuleFunctions - ERR22'
            if (DTPredictionInterval <= 0.0) then
                write (*,*) 'DT_PREDICTION_INTERVAL must be greater then 0'
                write (*,*) 'Module :',ClientModule
                stop 'ReadTimeKeyWords - ModuleFunctions - ERR25a'
            endif
        endif
        

        !Verifies Time Variables
        if (EndTime .le. BeginTime) then
            write (*,*) 'End Time is BEFORE model Start Time'
            write (*,*) 'Module :',ClientModule
            stop 'ReadTimeKeyWords - ModuleFunctions - ERR06'
        endif

        if (EndTime .lt. BeginTime+DT) then
            write (*,*) 'End Time is BEFORE model Start Time + DT'
            write (*,*) 'Module :',ClientModule
            stop 'ReadTimeKeyWords - ModuleFunctions - ERR07'
        endif

        !Verifies the DT
        if (.not. VariableDT) then

            !Run period in seconds
            aux  = EndTime - BeginTime
        
            !The run period must be a multiple of the model DT
            !The abs function is used, to avoid rounding erros
            !The old way was removed, to be able to run with Timesteps lower then 1 sec
            !Frank Dec - 2000
            MinError = min (abs(mod (aux, dble(DT))), abs(dble(DT) - mod (aux, dble(DT))))
            if (MinError >= 1.e-5.and.(.not.VariableDT)) then
                write(*,*) 
                write(*,*)' Time step error - Run period must be a multiple of DT'
                write(*,*)' Shorten the run time by ', mod (aux, dble(DT))
                write(*,*)' or increase the run time by ', DT - mod (aux, dble(DT))
                stop 'ReadTimeKeyWords - ModuleFunctions - ERR08'
            endif
        endif

    end subroutine ReadTimeKeyWords

    !--------------------------------------------------------------------------

    subroutine ConstructPropertyID (PropertyID, ObjEnterData, ExtractType)

        !Arguments-------------------------------------------------------------
        type (T_PropertyID)                         :: PropertyID
        integer                                     :: ObjEnterData
        integer                                     :: ExtractType

        !Local-----------------------------------------------------------------
        integer                                     :: flag
        integer                                     :: STAT_CALL

        !Property Name
        call GetData(PropertyID%Name, ObjEnterData, flag,                                &
                     SearchType   = ExtractType,                                         &
                     keyword      = 'NAME',                                              &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyID - ModuleFunctions - ERR01'
        if (flag==0) then
            write (*,*)'Property without name'
            stop 'ConstructPropertyID - ModuleFunctions - ERR02'
        endif

        if (.not. CheckPropertyName (PropertyID%Name, PropertyID%IDnumber)) then
            write (*, *)'The property isnt recognized by the model :', trim(PropertyID%Name)
            stop 'ConstructPropertyID - ModuleFunctions - ERR03'
        endif
 
        !Units
        call GetData(PropertyID%Units, ObjEnterData, flag,                               &
                     SearchType   = ExtractType,                                         &
                     keyword      = 'UNITS',                                             &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyID - ModuleFunctions - ERR04'

        !Description
        call GetData(PropertyID%Description, ObjEnterData, flag,                               &
                     SearchType   = ExtractType,                                         &
                     keyword      = 'DESCRIPTION',                                       &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ConstructPropertyID - ModuleFunctions - ERR05'


    end subroutine ConstructPropertyID

    !----------------------------------------------------------------------

#ifdef _USE_SEQASSIMILATION
    logical function Check_Hydrodynamic_Property(Property)

        !Arguments-------------------------------------------------------------
        integer, intent (IN) :: Property

        !----------------------------------------------------------------------

cd1 :   if ((Property == WaterLevel_            ) .OR.  (Property == VelocityModulus_       ) .OR.          &
            (Property == FlowModulus_           ) .OR.  (Property == VelocityU_             ) .OR.          &
            (Property == VelocityV_             ) .OR.  (Property == VelocityW_             ) .OR.          &
            (Property == WaterFluxX_            ) .OR.  (Property == WaterFluxY_            ) .OR.          &
            (Property == CoriolisX_             ) .OR.  (Property == BaroclinicForceX_      ) .OR.          &
            (Property == HorizontalTransportX_  ) .OR.  (Property == CoriolisY_             ) .OR.          &
            (Property == BaroclinicForceY_      ) .OR.  (Property == HorizontalTransportY_  ) .OR.          &
            (Property == BarotropicVelocityU_   ) .OR.  (Property == BarotropicVelocityV_   ) .OR.          &
            (Property == BaroclinicVelocityU_   ) .OR.  (Property == BaroclinicVelocityV_   ) .OR.          &
            (Property == ObstacleDragCoef_      ) .OR.  (Property == Vorticity_             ) .OR.          &
            (Property == BaroclinicKE_          ) .OR.  (Property == PerturbationPE_        ) .OR.          &
            (Property == KineticEnergy_         ) .OR.  (Property == ShearVelocity_         ) .OR.          &
            (Property == ShearStress_           )) then
            !review to include only possible state variables
            Check_Hydrodynamic_Property = .TRUE.

        else

            Check_Hydrodynamic_Property = .FALSE.

        end if cd1

    end function Check_Hydrodynamic_Property

    !----------------------------------------------------------------------

    logical function Check_Water_Property(Property)

        !Arguments-------------------------------------------------------------
        integer, intent (IN) :: Property

        !----------------------------------------------------------------------

cd1 :   if ((Property == Temperature_           ) .OR.  (Property == Salinity_              ) .OR.          &
            (Property == Phytoplankton_         ) .OR.                                                      &
            (Property == Zooplankton_           ) .OR.  (Property == DOPRefractory_         ) .OR.          &
            (Property == DOPNon_Refractory_     ) .OR.  (Property == DONRefractory_         ) .OR.          &
            (Property == DONNon_Refractory_     ) .OR.  (Property == Inorganic_Phosphorus_  ) .OR.          &
            (Property == POC_                   ) .OR.  (Property == POP_                   ) .OR.          &
            (Property == DOC_                   ) .OR.  (Property == DOP_                   ) .OR.          &            
            (Property == DON_                   ) .OR.  (Property == DOCsl_                 ) .OR.          &  
            (Property == DOPsl_                 ) .OR.  (Property == DONsl_                 ) .OR.          &  
            (Property == Ammonia_               ) .OR.  (Property == Nitrate_               ) .OR.          &  
            (Property == Silicate_              ) .OR.  (Property == BioSilica_             ) .OR.          &  
            (Property == CarbonDioxide_         ) .OR.  (Property == Oxygen_                ) .OR.          &
            (Property == DissolO2PercentSat_    ) .OR.  (Property == Diatom_C_              ) .OR.          &
            (Property == Diatom_N_              ) .OR.  (Property == Diatom_P_              ) .OR.          &
            (Property == Diatom_Si_             ) .OR.  (Property == Diatom_Chl_            ) .OR.          &
            (Property == Mix_Flagellate_C_      ) .OR.  (Property == Mix_Flagellate_N_      ) .OR.          &
            (Property == Mix_Flagellate_P_      ) .OR.  (Property == Mix_Flagellate_Chl_    ) .OR.          &
            (Property == Picoalgae_C_           ) .OR.  (Property == Picoalgae_N_           ) .OR.          &
            (Property == Picoalgae_P_           ) .OR.  (Property == Picoalgae_Chl_         ) .OR.          &  
            (Property == Flagellate_C_          ) .OR.  (Property == Flagellate_N_          ) .OR.          &  
            (Property == Flagellate_P_          ) .OR.  (Property == Flagellate_Chl_        ) .OR.          &
            (Property == Microzooplankton_C_    ) .OR.  (Property == Microzooplankton_N_    ) .OR.          &
            (Property == Microzooplankton_P_    ) .OR.  (Property == Het_Nanoflagellate_C_  ) .OR.          &  
            (Property == Het_Nanoflagellate_N_  ) .OR.  (Property == Het_Nanoflagellate_P_  ) .OR.          &  
            (Property == Mesozooplankton_C_     ) .OR.  (Property == Mesozooplankton_N_     ) .OR.          &  
            (Property == Mesozooplankton_P_     ) .OR.  (Property == Het_Bacteria_C_        ) .OR.          &
            (Property == Het_Bacteria_N_        ) .OR.  (Property == Het_Bacteria_P_        ) .OR.          &
            (Property == Nitrite_               ) .OR.  (Property == BOD_                   ) .OR.          &
            (Property == Cohesive_Sediment_     ) .OR.  (Property == Fecal_Coliforms_       ) .OR.          &
            (Property == E_Coli_                ) .OR.  (Property == T90_                   ) .OR.          &
            (Property == Oil_                   ) .OR.  (Property == Ciliate_               ) .OR.          &
            (Property == Bacteria_              ) .OR.  (Property == ParticulateArsenic_    ) .OR.          &
            (Property == DissolvedArsenic_      ) .OR.  (Property == Larvae_                ) .OR.          &
            (Property == Age_                   ) .OR.  (Property == Fish_                  ) .OR.          &
            (Property == FishFood_              ) .OR.  (Property == MacroAlgae_            ) .OR.          &
            (Property == DriftingMacroAlgae_    ) .OR.  (Property == MicroPhytoBenthos_     ) .OR.          &
            (Property == AdsorbedAmmonia_       ) .OR.  (Property == RefreactaryOrganicN_   ) .OR.          &
            (Property == Ngas_                  ) .OR.  (Property == HeterotrophicN_        ) .OR.          &
            (Property == AnaerobicN_            ) .OR.  (Property == AutotrophicN_          ) .OR.          &
            (Property == AnaerobicC_            ) .OR.  (Property == AutotrophicC_          ) .OR.          &
            (Property == HeterotrophicC_        ) .OR.  (Property == LabileOrganicC_        ) .OR.          &
            (Property == RefreactaryOrganicC_   ) .OR.  (Property == GenericProperty_       ) .OR.          &
            (Property == GrossProd_             ) .OR.  (Property == NetProd_               ) .OR.          &
            (Property == NutrientLim_           ) .OR.  (Property == NLim_                  ) .OR.          &
            (Property == PLim_                  ) .OR.  (Property == LightLim_              ) .OR.          &
            (Property == TemperatureLim_        ) .OR.  (Property == SalinityLim_           ) .OR.          &
            (Property == DiaGrossProd_          ) .OR.  (Property == DiaNutrientLim_        ) .OR.          &
            (Property == DiaNLim_               ) .OR.  (Property == DiaPLim_               ) .OR.          &
            (Property == DiaSiLim_              ) .OR.  (Property == DiaLightLim_           ) .OR.          &
            (Property == DiaTemperatureLim_     ) .OR.  (Property == Diatoms_               ) .OR.          &
            (Property == ParticulateContaminant_) .OR.  (Property == DissolvedContaminant_  ) .OR.          &
            (Property == Sediment               ) .OR.  (Property == DissolvedSodium_       ) .OR.          &
            (Property == DissolvedCalcium_      ) .OR.  (Property == ParticulateSodium_     ) .OR.          &
            (Property == ParticulateCalcium_    ) .OR.  (Property == RPOM_                  ) .OR.          &
            (Property == LPOM_                  ) .OR.  (Property == LDOM_                  ) .OR.          &
            (Property == RDOM_                  ) .OR.  (Property == PSilica_               ) .OR.          &
            (Property == DSilica_               ) .OR.  (Property == ICarbon_               ) .OR.          &
            (Property == pH_                    ) .OR.  (Property == HCO3_                  ) .OR.          &
            (Property == CO3_                   ) .OR.  (Property == Algae_1_               ) .OR.          &
            (Property == Algae_2_               ) .OR.  (Property == Algae_3_               ) .OR.          &
            (Property == Algae_4_               ) .OR.  (Property == Algae_5_               ) .OR.          &
            (Property == Epiphyton_1_           ) .OR.  (Property == Epiphyton_2_           ) .OR.          &
            (Property == Epiphyton_3_           ) .OR.  (Property == Epiphyton_4_           ) .OR.          &
            (Property == Epiphyton_5_           ) .OR.  (Property == Alkalinity_            ) .OR.          &
            (Property == Detritus_              ) .OR.  (Property == ANLim_                 ) .OR.          &
            (Property == APLim_                 ) .OR.  (Property == ASLim_                 ) .OR.          &
            (Property == ALightLim_             ) .OR.  (Property == AOverallLim_           ) .OR.          &
            (Property == ENLim_                 ) .OR.  (Property == EPLim_                 ) .OR.          &
            (Property == ESLim_                 ) .OR.  (Property == ELightLim_             ) .OR.          &
            (Property == EOverallLim_           ) .OR.  (Property == NH4D_                  ) .OR.          &
            (Property == NO3D_                  ) .OR.  (Property == LDOMD_                 ) .OR.          &
            (Property == RDOMD_                 ) .OR.  (Property == LPOMD_                 ) .OR.          &
            (Property == RPOMD_                 ) .OR.  (Property == LRDOMD_                ) .OR.          &
            (Property == LRPOMD_                ) .OR.  (Property == CBODD_                 ) .OR.          &
            (Property == PO4ER_                 ) .OR.  (Property == PO4EG_                 ) .OR.          &
            (Property == PO4AR_                 ) .OR.  (Property == PO4AG_                 ) .OR.          &
            (Property == PO4OM_                 ) .OR.  (Property == PO4BOD_                ) .OR.          &
            (Property == NH4ER_                 ) .OR.  (Property == NH4EG_                 ) .OR.          &
            (Property == NH4AR_                 ) .OR.  (Property == NH4AG_                 ) .OR.          &
            (Property == NH4OM_                 ) .OR.  (Property == NH4BOD_                ) .OR.          &
            (Property == NO3AG_                 ) .OR.  (Property == NO3EG_                 ) .OR.          &
            (Property == DSIAG_                 ) .OR.  (Property == DSIEG_                 ) .OR.          &
            (Property == DSID_                  ) .OR.  (Property == PSIAM_                 ) .OR.          &
            (Property == PSID_                  ) .OR.  (Property == LDOMAP_                ) .OR.          &
            (Property == LDOMEP_                ) .OR.  (Property == LPOMAP_                ) .OR.          &
            (Property == DOAP_                  ) .OR.  (Property == DOEP_                  ) .OR.          &
            (Property == DOAR_                  ) .OR.  (Property == DOER_                  ) .OR.          &
            (Property == DOOM_                  ) .OR.  (Property == DONIT_                 ) .OR.          &
            (Property == ICarbonAP_             ) .OR.  (Property == ICarbonEP_             ) .OR.          &
            (Property == ICarbonBOD_            ) .OR.  (Property == Diameter_              ) .OR.          &
            (Property == Percentage_            ) .OR.  (Property == D35_                   ) .OR.          &
            (Property == D50_                   ) .OR.  (Property == D90_                   ) .OR.          &
            (Property == BedRock_               ) .OR.  (Property == SandTauCritic_         ) .OR.          &
            (Property == Sand_                  ) .OR.  (Property == TransportCapacity_     ) .OR.          &
            (Property == TransportCapacityX_    ) .OR.  (Property == TransportCapacityY_    ) .OR.          &
            (Property == ConsolidationFlux_     ) .OR.  (Property == Porosity_              ) .OR.          &
            (Property == TSS_                   ) .OR.  (Property == COHSED_FINE_           ) .OR.          &
            (Property == COHSED_MEDIUM_         ) .OR.  (Property == COHSED_COARSE_         ) .OR.          &
            (Property == DissolvedMetal_        ) .OR.  (Property == DissolvedCopper_       ) .OR.          &
            (Property == DissolvedCadmium_      ) .OR.  (Property == DissolvedLead_         ) .OR.          &
            (Property == DissolvedMercury_      ) .OR.  (Property == DissolvedZinc_         )   ) then

            Check_Water_Property = .TRUE.
        
        else

            Check_Water_Property = .FALSE.

        end if cd1
           
    end function Check_Water_Property
#endif

    !--------------------------------------------------------------------------

    subroutine GetDataOnlineString (Keyword, CharData, ArrayData, RealData, IntData)

        !Arguments------------------------------------------------
        character (Len=*), intent (IN)              :: Keyword
        character (Len=*), intent (OUT), optional   :: CharData
        real,  pointer,  dimension(:)  , optional   :: ArrayData
        real,              intent (OUT), optional   :: RealData
        integer,           intent (OUT), optional   :: IntData
        !Local----------------------------------------------------
        character (Len=StringLength)                :: AuxString
        logical                                     :: FoundKeyword
        integer                                     :: STAT_CALL
        integer                                     :: n, nmax, s, lk, ls, sk, asl

        !Begin----------------------------------------------------

        lk = len_trim(Keyword     )
        ls = len_trim(OnLineString)

        FoundKeyword = .false.

        do s = 1, ls - lk
        
            if ('&'//Keyword(1:lk)//'=' == OnLineString(s:s+lk+1)) then
                FoundKeyword = .true.
                sk = s + lk 
                exit
            endif

        enddo

        if (.not. FoundKeyword) then
            write(*,*) ' Not found ', trim(Keyword),' in online string'
            stop 'Modulefunctions - GetDataOnlineString - ERR10'
        endif

        s = Scan(OnLineString(sk:ls),'&')

        s = sk + s - 1

        AuxString=OnLineString(sk+2:s-1)

        asl = len_trim(AuxString)

        if (present(CharData)) then
            CharData = AuxString(1:asl)
        endif

        do s = 1, asl
            if (AuxString(s:s) == '_') AuxString(s:s) = ' '
        enddo

        if (present(ArrayData)) then

            nmax = size(ArrayData)

            read(AuxString(1:asl),*,IOSTAT=STAT_CALL) (ArrayData(n),n=1,Nmax)
            if (STAT_CALL /= SUCCESS_) stop 'Modulefunctions - GetDataOnlineString - ERR30'
        endif


        if (present(IntData)) then
            read(AuxString(1:asl),*,IOSTAT=STAT_CALL) IntData
            if (STAT_CALL /= SUCCESS_) stop 'Modulefunctions - GetDataOnlineString - ERR40'
        endif

        if (present(RealData)) then
            read(AuxString(1:asl),*,IOSTAT=STAT_CALL) RealData
            if (STAT_CALL /= SUCCESS_) stop 'Modulefunctions - GetDataOnlineString - ERR50'
        endif

    end subroutine GetDataOnlineString

    !--------------------------------------------------------------------------

    logical function IsOdd(Number)

        !Arguments-------------------------------------------------------------
        integer                                         :: Number

        !Begin-----------------------------------------------------------------

        IsOdd = .false.

        if(Number/2. > int(Number/2.)) IsOdd = .true.

    end function

    !--------------------------------------------------------------------------


    subroutine Secant (x,x0,d,coefa,coefb)

        !Arguments-------------------------------------------------------------
        integer  :: maxiterations=100
        real     :: eps=5.0e-5
        real     :: x1,x2,f1,f2,csi

        !Local-----------------------------------------------------------------
        real     :: x,x0,d,coefa,coefb
        integer  :: i

        !Begin-----------------------------------------------------------------

        x  = x0
        x2 = x0 + 0.01
        f2 = x2 - coefa * tanh(coefb * d / x2)

        do i = 1, maxiterations

            f1 = x - coefa * tanh(coefb * d / x)

            ! ---> convergence test

            if (abs(f1) .lt. eps) return
                x1 = x
            if (abs(f1) .lt. abs(f2)) then
                csi= f1/f2
                x = x1 + csi/(csi - 1.0)*(x2 - x1)
            else
                csi = f2/f1
                x = x1 + (x2 - x1)/(1.0 - csi)
            endif

            f2 = f1
            x2 = x1

        end do

    end subroutine Secant

    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    subroutine tridag(a,b,c,r,x,gam, jmin, jmax)

        !Arguments
        real(8), dimension(:), pointer:: b, gam
        real   , dimension(:), pointer:: a, c, r, x
        integer                       :: jmin, jmax
        !Local 
        integer                       :: j
        real(8)                       :: bet

        !Begin

        if(b(jmin).eq.0.) stop 'tridag: rewrite equations'
        bet=b(jmin)
        x(jmin)=r(jmin)/bet

        do j=jmin+1,jmax
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j)* gam(j)
            if(bet.eq.0.) stop 'tridag failed'
            x(j)=(r(j)-a(j)*x(j-1))/bet
        enddo

        do j=jmax-1,jmin,-1
            x(j)=x(j)-gam(j+1)*x(j+1)
        enddo

    end subroutine tridag

    
    subroutine Tridag_cyclic(a, b, c, alpha, beta, r, x, gam, bb, u, z, imin, imax)
        
        !Arguments-------------------------------------------------------------
        real                          :: alpha, beta
        real(8), dimension(:), pointer:: b, bb, gam
        real,    dimension(:), pointer:: a, c, r, x, u, z
        integer                       :: imin, imax

        !Local----------------------------------------------------------------- 
        integer                       :: i
        real                          :: fact,gamma

        !Begin-----------------------------------------------------------------

        if((imax-imin +1) < 2)stop 'imax too small in cyclic'

        gamma=-b(imin)
        bb(imin)=b(imin)-gamma
        bb(imax)=b(imax)-alpha*beta/gamma
        do  i=imin+1,imax-1
            bb(i)=b(i)
        enddo

        call tridag(a,bb, c, r, x, gam, imin, imax)

        u(imin)=gamma
        u(imax)=alpha

        do i=imin+1,imax-1
            u(i)=0.
        enddo

        call tridag(a,bb, c, u, z, gam, imin, imax)
                   
        fact=(x(imin)+beta*x(imax)/gamma)/(1.+z(imin)+beta*z(imax)/gamma)

        do i=imin,imax
            x(i)=x(i)-fact*z(i)
        enddo

    end subroutine Tridag_cyclic

    !--------------------------------------------------------------------------

    !Computes de value for a depth based in a profile. 
    real function InterpolateProfile (CellDepth, NDEPTHS, Depth, Values)

        !Arguments-------------------------------------------------------------
        real                         :: CellDepth
        integer                      :: NDEPTHS
        real, dimension (:), pointer :: Depth, Values

        !Local-----------------------------------------------------------------
        real(8)                      :: dx
        integer                      :: i

        !Begin-----------------------------------------------------------------

        if (Depth(1) > Depth(NDEPTHS)) then

            if (CellDepth >= Depth(1)) then

                InterpolateProfile = Values(1)

            else if (CellDepth <= Depth(NDEPTHS)) then
        
                InterpolateProfile = Values(NDEPTHS)

            else  

                do i = 1 , NDEPTHS

                    if (CellDepth <= Depth(i) .and. CellDepth >  Depth(i + 1)) then

                        dx = dble(Depth(i) - CellDepth) / dble(Depth(i) - Depth(i + 1))

                        InterpolateProfile = Values(i+1) * dx + Values(i  ) *  (1. - dx) 
                        
                        exit

                    endif


                enddo

            endif     
        else

            if (CellDepth <= Depth(1)) then

                InterpolateProfile = Values(1)

            else if (CellDepth >= Depth(NDEPTHS)) then
        
                InterpolateProfile = Values(NDEPTHS)

            else  

                do i = 1 , NDEPTHS

                    if (CellDepth > Depth(i) .and. CellDepth <=  Depth(i + 1)) then

                        dx = (Depth(i) - CellDepth) / (Depth(i) - Depth(i + 1))

                        InterpolateProfile = Values(i+1) * dx + Values(i  ) *  (1 - dx)

                        exit 

                    endif


                enddo

            endif     

        endif


    end function InterpolateProfile


    !--------------------------------------------------------------------------

    !Computes de value for a depth based in a profile. 
    real(8) function InterpolateProfileR8 (CellDepth, NDEPTHS, Depth, Values, FoundBottom, FoundSurface)

        !Arguments-------------------------------------------------------------
        real(8)                         :: CellDepth
        integer                         :: NDEPTHS
        real(8), dimension (:)          :: Depth, Values
        logical, optional               :: FoundBottom, FoundSurface

        !Local-----------------------------------------------------------------
        real(8)                         :: dx
        integer                         :: i

        !Begin-----------------------------------------------------------------

        if (present(FoundBottom )) FoundBottom  = .false.
        if (present(FoundSurface)) FoundSurface = .false.

        if (CellDepth >= Depth(1)) then

            InterpolateProfileR8 = Values(1)

            if (present(FoundBottom ))  FoundBottom = .true.

        else if (CellDepth <= Depth(NDEPTHS)) then
    
            InterpolateProfileR8 = Values(NDEPTHS)
            if (present(FoundSurface))  FoundSurface = .true.
        else  

            do i = 1 , NDEPTHS

                if (CellDepth <= Depth(i) .and. CellDepth >  Depth(i + 1)) then

                    dx = (Depth(i) - CellDepth) / (Depth(i) - Depth(i + 1))

                    InterpolateProfileR8 = Values(i+1) * dx + Values(i  ) *  (1. - dx) 

                    exit 
                    
                endif


            enddo

        endif     

    end function InterpolateProfileR8


    !--------------------------------------------------------------------------

    !Computes de value for a depth based in a profile. 
    real(8) function QuadraticInterpolProfile (CellDepth, NDEPTHS, Depth, Values, FoundBottom, FoundSurface)

        !Arguments-------------------------------------------------------------
        real(8)                         :: CellDepth
        integer                         :: NDEPTHS
        real(8), dimension (:)          :: Depth, Values
        logical                         :: FoundBottom, FoundSurface

        !Local-----------------------------------------------------------------
        real(8)                         :: Depth1, Depth2
        integer                         :: i

        !Begin-----------------------------------------------------------------

        FoundBottom  = .false.
        FoundSurface = .false.

        if (CellDepth < Depth(NDEPTHS-1)) then

            QuadraticInterpolProfile = QuadraticInterpolationR8(Values(NDEPTHS),   Values(NDEPTHS-1), &
                                                                Values(NDEPTHS-2), Depth(NDEPTHS),    &
                                                                Depth(NDEPTHS-1),  Depth(NDEPTHS-2), CellDepth)
            if  (CellDepth <= Depth(NDEPTHS)) FoundSurface = .true.

        else if (CellDepth >  Depth(2      )) then

            QuadraticInterpolProfile = QuadraticInterpolationR8(Values(3), Values(2), Values(1), &
                                                                Depth(3), Depth(2), Depth(1), CellDepth)
            if  (CellDepth >= Depth(1      )) FoundBottom  = .true.
        else 

            do i = 1 , NDEPTHS

                Depth1 = (Depth(i    ) + Depth(i + 1)) / 2.
                Depth2 = (Depth(i + 1) + Depth(i + 2)) / 2.

                if (CellDepth <= Depth1 .and. CellDepth >  Depth2) then
                    QuadraticInterpolProfile = QuadraticInterpolationR8(Values(i+2), Values(i+1),   &
                                                                        Values(i), Depth(i+2),      &
                                                                        Depth(i+1), Depth(i), CellDepth)
                    exit
                endif


            enddo

        endif
    

    end function QuadraticInterpolProfile


    !--------------------------------------------------------------------------
!--------------------------------------------------------------------------

    !Computes de value for a depth based in a profile. 
    real(8) function PolIntProfile (CellDepth, NDEPTHS, Depth, Values, Npoli, IsEven, Error)

        !Arguments-------------------------------------------------------------
        real(8)                         :: CellDepth
        integer                         :: NDEPTHS
        real(8), dimension (:)          :: Depth, Values
        integer                         :: Npoli
        logical                         :: IsEven !Even = Par
        real(8), optional               :: Error

        !Local-----------------------------------------------------------------
        real(8)                         :: Depth1, Depth2, Error_
        integer                         :: i, di, iupper, ilower

        !Begin-----------------------------------------------------------------

        di = int(NPoli/2.) 

        if (.not. IsEven) then
            
            di = di + 1

        endif


        if (CellDepth < Depth(NDEPTHS - di)) then
            
            iupper = NDEPTHS
            ilower = NDEPTHS - Npoli

        else if (CellDepth >  Depth(1 + di)) then
            iupper = Npoli + 1
            ilower = 1
        else 

            if (IsEven) then

                do i = 1 , NDEPTHS

                    Depth1 = (Depth(i    ) + Depth(i + 1)) / 2.
                    Depth2 = (Depth(i + 1) + Depth(i + 2)) / 2.

                    if (CellDepth <= Depth1 .and. CellDepth >  Depth2) then
                        iupper = i + 1 + di
                        ilower = i + 1 - di
                        exit
                    endif
                enddo

            else
                do i = 1 , NDEPTHS
                    if (CellDepth <= Depth(i) .and.  CellDepth > Depth(i + 1)) then
                        iupper = i + 1 + int(NPoli/2.)
                        ilower = i     - int(NPoli/2.)
                        exit
                    endif
                enddo
            endif

        endif
        
        call PolInt(Depth(ilower:iupper), Values(ilower:iupper),Npoli + 1,CellDepth,PolIntProfile,Error_)

        if (present(Error)) Error = Error_
    

    end function PolIntProfile


    !--------------------------------------------------------------------------

    subroutine RelativePosition4VertPolygon(Xa, Ya, Xb, Yb, Xc, Yc, Xd, Yd, Xe, Ye, Xex, Yey)

        !Arguments---------------------------------------------------------------
        real                :: Xa, Ya, Xb, Yb, Xc, Yc, Xd, Yd, Xe, Ye, Xex, Yey
        
        !Local-------------------------------------------------------------------
        real(8)             :: DXdc, DYac, DXba, DYbd, DXef, DYeg 
        real(8)             :: MinDx, SumAux
        real(8)             :: a1, b1, a2, b2, a3, b3, a4, b4
        real(8)             :: Seg_ac, Seg_dc,Seg_hc, Seg_ic
        real(8)             :: Xf, Yf, Xg, Yg, Xh, Yh, Xi, Yi, TgX, TgY
        real(8)             :: XaR8, YaR8, XbR8, YbR8, XcR8, YcR8, XdR8, YdR8, XeR8, YeR8
        !Begin-------------------------------------------------------------------

        XaR8 = dble(Xa); YaR8 = dble(Ya); XbR8 = dble(Xb); YbR8 = dble(Yb); XcR8 = dble(Xc)
        YcR8 = dble(Yc); XdR8 = dble(Xd); YdR8 = dble(Yd); XeR8 = dble(Xe); YeR8 = dble(Ye)

        !the four segments of the cell
        DXdc = XdR8 - XcR8 
        DXba = XbR8 - XaR8
        DYac = YaR8 - YcR8 
        DYbd = YbR8 - YdR8
        

        SumAux = abs(DXdc + DXba + DYac + DYbd)
        MinDx  = SumAux * 1.e-16

        if (abs(DXdc)==0.) DXdc = MinDx
        if (abs(DXba)==0.) DXba = MinDx
        if (abs(DYac)==0.) DYac = MinDx
        if (abs(DYbd)==0.) DYbd = MinDx

        a1 = (XdR8*YcR8 - XcR8*YdR8) / DXdc
        b1 = (   YdR8 -    YcR8) / DXdc

        a2 = (XbR8*YaR8 - XaR8*YbR8) / DXba
        b2 = (   YbR8 -    YaR8) / DXba

        a3 = (XcR8*YaR8 - XaR8*YcR8) / DYac
        b3 = (   XaR8 -    XcR8) / DYac

        a4 = (XdR8*YbR8 - XbR8*YdR8) / DYbd
        b4 = (   XbR8 -    XdR8) / DYbd

        !intersection points 
        if (b2/=b1) then

            !F point intersection point of X faces
            Xf = (a1 - a2) / (b2 - b1)
            Yf = a1 + b1*Xf

            !H point intersection with segment CA
            DXef = XeR8 - Xf

            if (abs(DXef)==0.) DXef = MinDx
            TgX = (YeR8 - Yf) / DXef
        else
            !H point intersection with segment CA
            TgX = b1
        endif

        Yh  = (YeR8 + (a3 - XeR8) * TgX) / (1. - b3*TgX)
        Xh  = a3 + b3 * Yh


        if (b4/=b3) then

            !G point intersection point of X faces
            Yg = (a3 - a4) / (b4 - b3)
            Xg = a3 + b3*Yg

            !i point intersection with segment CA
            DYeg = YeR8 - Yg

            if (abs(DYeg)==0.) DYeg = MinDx
            TgY = (XeR8 - Xg) / DYeg
        else

            !i point intersection with segment CA
            TgY = b3
        endif

        Xi  = (XeR8 + (a1 - YeR8) * TgY) / (1. - b1*TgY)
        Yi  =  a1 + b1 * Xi


        Seg_ac = sqrt((XaR8-XcR8)*(XaR8-XcR8) + (YaR8-YcR8)*(YaR8-YcR8))
        Seg_dc = sqrt((XdR8-XcR8)*(XdR8-XcR8) + (YdR8-YcR8)*(YdR8-YcR8))
        Seg_hc = sqrt((Xh  -XcR8)*(Xh  -XcR8) + (Yh  -YcR8)*(Yh  -YcR8))
        Seg_ic = sqrt((Xi  -XcR8)*(Xi  -XcR8) + (Yi  -YcR8)*(Yi  -YcR8))

        Xex = Seg_ic / Seg_dc
        Yey = Seg_hc / Seg_ac

    end subroutine RelativePosition4VertPolygon
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------

#ifdef _USE_MPI

    integer function MPIKind0D(Variable0D)
        
        !Arguments-------------------------------------------------------------
        real                                :: Variable0D
        
        !Begin-----------------------------------------------------------------

        if(    kind(Variable0D) == 4)then
            MPIKind0D = MPI_REAL
        elseif(kind(Variable0D) == 8)then
            MPIKind0D = MPI_DOUBLE_PRECISION
        end if

    end function MPIKind0D

    !--------------------------------------------------------------------------


    integer function MPIKind1D(Variable1D)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:), pointer             :: Variable1D
        
        !Begin-----------------------------------------------------------------

        if(    kind(Variable1D) == 4)then
            MPIKind1D = MPI_REAL
        elseif(kind(Variable1D) == 8)then
            MPIKind1D = MPI_DOUBLE_PRECISION
        end if    
    
    end function MPIKind1D


    !--------------------------------------------------------------------------


    integer function MPIKind2D(Variable2D)
        
        !Arguments-------------------------------------------------------------
        real, dimension(:,:), pointer           :: Variable2D
        
        !Begin-----------------------------------------------------------------

        if(    kind(Variable2D) == 4)then
            MPIKind2D = MPI_REAL
        elseif(kind(Variable2D) == 8)then
            MPIKind2D = MPI_DOUBLE_PRECISION
        end if

    end function MPIKind2D


    !--------------------------------------------------------------------------

    integer function MPIKind3D(Variable3D)

        !Arguments-------------------------------------------------------------
        real, dimension(:,:,:), pointer         :: Variable3D
        
        !Begin-----------------------------------------------------------------

        if(    kind(Variable3D) == 4)then
            MPIKind3D = MPI_REAL
        elseif(kind(Variable3D) == 8)then
            MPIKind3D = MPI_DOUBLE_PRECISION
        end if

    end function MPIKind3D

#endif



    !----------------------------------------------------------------------------
    !function WQPhytoLightLimitationFactor computes the Phytoplankton light limitation factor
    !of each water mass.
    !

    function PhytoLightLimitationFactor(Thickness,                        &
                                        TopRadiation,                     &
                                        PExt,                             &
                                        Photoinhibition)      
        real ::  PhytoLightLimitationFactor

        !Arguments---------------------------------------------------------------
        real, intent(IN) :: Thickness                   !Thickness of the water mass
        real, intent(IN) :: TopRadiation                !Radiation at the water mass top
        real, intent(IN) :: PExt
        real, intent(IN) :: Photoinhibition     

        !Local variables---------------------------------------------------------

        real             :: xd, xe, xf

        !------------------------------------------------------------------------

        !Parsons, Takahashi & Hargrave -> Biological Oceanographic Processes

        xd   = exp(1.0) / (PExt * Thickness)

        xe   = -(TopRadiation / Photoinhibition) * exp(-PExt * Thickness)

        xf   = -(TopRadiation / Photoinhibition) 

        PhytoLightLimitationFactor = xd * (exp(xe) - exp(xf))

        !EPA---------------------------------------------------------------------

        !    xd   = exp(1.0) / (PExt * Thickness)
        !    xe   = -(TopRadiation / Photoinhibition) * (PExt * Thickness)
        !    xf   = -(TopRadiation / Photoinhibition)

        !    WQPhytoLightLimitationFactor = xd * exp(-exp(xe) - exp(xf))

        !------------------------------------------------------------------------

        !The light limitation factor can't be negative
cd1 :   if (PhytoLightLimitationFactor .LT. 0.0) then
            write(*, *) 'Light Limitation Factor negative', PhytoLightLimitationFactor
            write(*, *) 'Set to zero'
            write(*, *) 'Function PhytoLightLimitationFactor - WRN01'
            PhytoLightLimitationFactor = 0.
        end if cd1

    !------------------------------------------------------------------------

    end function PhytoLightLimitationFactor


    !----------------------------------------------------------------------------
    !Computes T90 deacy time according to Canteras et al. (1995)
    
    function ComputeT90_Canteras       (Temperature, Salinity, Radiation)      

        real ::  ComputeT90_Canteras

        !Arguments---------------------------------------------------------------
        real, intent(IN)    :: Temperature,  Salinity, Radiation

        !Local variables---------------------------------------------------------
        real                :: Mortality

        !Mortality per day
        Mortality  = 2.533 * (1.04**(Temperature - 20.)) * (1.012**Salinity)    &
                   + (0.113 * Radiation)

        ComputeT90_Canteras  =  (2.303 / Mortality) * 24. * 3600.

    end function ComputeT90_Canteras


    !------------------------------------------------------------------------
    !Fecal decay according to Chapra (1997)
    
    function ComputeT90_Chapra       (Temperature, Salinity, Light)      

        real ::  ComputeT90_Chapra

        !Arguments---------------------------------------------------------------
        real, intent(IN) :: Temperature,  Salinity, Light

        !Local variables---------------------------------------------------------
        real                :: Mortality

        !Mortality per day
        Mortality  = (0.8 + 0.02 * Salinity) * 1.07 ** (Temperature - 20.) + Light

        ComputeT90_Chapra  =  (2.303 / Mortality) * 24. * 3600.


    end function ComputeT90_Chapra

    !------------------------------------------------------------------------
    !Downward Longwave Radiation (Wunderlich et. al. 1968)

    real function LongWaveDownward (CloudCover, AirTemp)

        !Arguments---------------------------------------------------------------
        real, intent(IN)    :: CloudCover, AirTemp
        
        LongWaveDownward = 0.97 * StefanBoltzmann * 0.937e-5 * (1.0 + 0.17 * CloudCover) * (AirTemp + 273.15) ** 6.0

    end function LongWaveDownward

    !------------------------------------------------------------------------
    !Downward Longwave Radiation (Wunderlich et. al. 1968)

    real function LongWaveUpward (WaterTemp)

        !Arguments---------------------------------------------------------------
        real, intent(IN)    :: WaterTemp
        
        LongWaveUpward = -0.97 * StefanBoltzmann * (WaterTemp + 273.15) ** 4.0

    end function LongWaveUpward

    !--------------------------------------------------------------------------
    !Clausius Clapeyron - Water Modeling Review - Central Valley 2000 eq. 2-35
    real function SaturatedVaporPressure(Temperature)

        !Arguments-------------------------------------------------------------
        real                                        :: Temperature

        !Local-----------------------------------------------------------------
        real                                        :: a = 6.108   ![mb]
        real                                        :: b = 17.27
        real                                        :: c = 273.2   ![ºC]

        SaturatedVaporPressure = a * exp(b * Temperature / (Temperature + c))

    end function SaturatedVaporPressure

    !--------------------------------------------------------------------------
    !Water Modeling Review - Central Valley 2000 Table 2.6

    real function LatentHeatOfVaporization (Temperature)

        !Arguments-------------------------------------------------------------
        real                                        :: Temperature

        LatentHeatOfVaporization = 1000 * (2499 - 2.36 * Temperature)

    end function LatentHeatOfVaporization

    !--------------------------------------------------------------------------
    !Water Modeling Review - Central Valley 2000 Table 2.8

    real function WindFunction (WindSpeed)

        !Arguments-------------------------------------------------------------
        real                                        :: WindSpeed

        WindFunction = 4.18e-9 + 0.95e-9 * WindSpeed

    end function WindFunction

    !--------------------------------------------------------------------------
    !Water Modeling Review - Central Valley 2000 eq. 2-45

    real function LatentHeat (WaterDensity, WaterTemp, AirTemp, RelHumidity, WindVel)

        !Arguments-------------------------------------------------------------
        real                                        :: WaterDensity, WaterTemp, AirTemp
        real                                        :: RelHumidity, WindVel

        LatentHeat     = - WaterDensity * LatentHeatOfVaporization(WaterTemp) * &
                           (SaturatedVaporPressure(WaterTemp) -  RelHumidity *  &
                            SaturatedVaporPressure(AirTemp)) * WindFunction(WindVel)

    end function LatentHeat

    !--------------------------------------------------------------------------
    !Water Modeling Review - Central Valley 2000 eq. 2-47

    real function SensibleHeat (WaterDensity, WaterTemp, AirTemp, WindVel)

        !Arguments-------------------------------------------------------------
        real                                        :: WaterDensity, WaterTemp, AirTemp
        real                                        :: WindVel

        !Local-----------------------------------------------------------------
        real, parameter                             :: Cb = 0.61 ![mb]

        SensibleHeat     = WaterDensity * LatentHeatOfVaporization(WaterTemp) * &
                           WindFunction(WindVel) * Cb * (AirTemp - WaterTemp)

    end function SensibleHeat



    !--------------------------------------------------------------------------
    !Several Formulations to calculate Aeration Flux (Oxygen)
    real function AerationFlux (AerationEquation, WindVelocity, WaterTemperature)

        !Arguments-------------------------------------------------------------
        integer                                     :: AerationEquation
        real                                        :: WindVelocity, WaterTemperature
        
        !Local-----------------------------------------------------------------
        real,    parameter                          :: ReaerationCoefficient    = 1.024
        real                                        :: a, bcoef, KL, DMO2

        select case(AerationEquation)


!            case(Broecker_et_al_1978)

!                write(*,*)'Aeration method Broecker et al, 1978, not available.'
!                stop 'AerationFlux - ModuleInterfaceWaterAir - ERR01'

            case(Gelda_et_al_1996)

                !Coeficient from Gelda
                if(WindVelocity <= 3.5) then
                    a     = 0.2
                    bcoef = 1.0
                else
                    a     = 0.057
                    bcoef = 2.0
                endif

                !KL in m/day    
                KL = a * WindVelocity ** bcoef

            case(Banks_Herrera_1977)

                KL = (0.728*sqrt(WindVelocity)-0.317*WindVelocity+0.0372*WindVelocity**2)

            case(Wanninkhof_et_al_1991)
                
                KL = 0.0986*WindVelocity**1.64

            case(Chen_Kanwisher_1963)

                DMO2 = 2.04E-9
                KL   = 86400.0*DMO2/((200.0-60.0*sqrt(min(WindVelocity,11.0)))*1.E-6)


            case(Cole_Buchak_1993)
                
                KL = (0.5+0.05*WindVelocity*WindVelocity)

            case(Banks_1975)
                
                if (WindVelocity <= 5.5) then
                    KL = 0.362 * sqrt(WindVelocity)
                else
                    KL = 0.0277* WindVelocity**2
                endif


            case(Smith_1978)

                KL = 0.64+0.128*WindVelocity**2

            case(Liss_1973)

                if (WindVelocity <= 4.1) then
                    KL = 0.156  * WindVelocity**0.63
                else
                    KL = 0.0269 * WindVelocity**1.9
                endif

            case(Downing_Truesdale_1955)

                KL = 0.0276*WindVelocity**2 

            case(Kanwisher_1963)

                KL = 0.0432*WindVelocity**2
            
            case(Yu_et_al_1977)

                KL = 0.319*WindVelocity

            case(Weiler_1974)

                if (WindVelocity <= 1.6) then                      
                    KL = 0.398                    
                else
                    KL = 0.155 * WindVelocity**2
                endif

            case default

                write(*,*)'Unknown Aeration method'
                stop 'AerationFlux - ModuleFunctions - ERR01'

        end select

        !minimum value of KL in m/day
        if (KL <= 0.6) KL = 0.6

        !temperature correction
        KL = KL * ReaerationCoefficient**(WaterTemperature - 20.0)

        !convert from m/d to m/s
        KL = KL / 86400.             

        AerationFlux = KL
        return

    end function AerationFlux



    !--------------------------------------------------------------------------
    !Several Formulations to calculate Aeration Flux (Carbon Dioxide)
    real function AerationFlux_CO2 (CO2AerationEquation, WindVelocity, WaterVelocity, WaterTemperature, WaterSalinity)
    
        !Arguments-------------------------------------------------------------
            integer                                     :: CO2AerationEquation
            real                                        :: WindVelocity, WaterVelocity, WaterTemperature
            real                                        :: WaterDepth, WaterSalinity
            real                                        :: fTS, Sc0, Sc35, Alfa, ko, k
            
            !Local-----------------------------------------------------------------
    
    
        WaterDepth = 1
        
        ! OConnor_Dobbins_1958
        ko = 1.719 * ((WaterVelocity * 100/WaterDepth)**(1./2.))
        
        
        select case(CO2AerationEquation)
        
            case(Borges_et_al_2004)
            
                Sc0  = 1800.6 - (120.1 * WaterTemperature) + (3.7818 * WaterTemperature**2.)     &
                       - (0.047608 * WaterTemperature**3.)
                      
                Sc35 = 1953.4 - (128.0 * WaterTemperature) + (3.9918 * WaterTemperature**2.)     &
                       - (0.050091 * WaterTemperature**3.)
                 
                fTS  = ((Sc35 - Sc0) * WaterSalinity) / 35.0 + Sc0       
                                
                Alfa = (600/fTS)**(1./2.)
            
                k = ((1. + 2.58 * WindVelocity) + ko) * Alfa
            
            
                        
            case(Carini_et_al_1996)
            
                fTS = 2073.1 - (125.62 * WaterTemperature) + (3.6276 * (WaterTemperature**2.))    &
                     - (0.043219 * (WaterTemperature**3.))
                 
                Alfa = (600/fTS)**(1./2.)
            
                k = ((0.045 + 2.0277 * WindVelocity) + ko) * Alfa 
            
            
            
            case(Raimond_Cole_2001)
            
                 fTS = 2073.1 - (125.62 * WaterTemperature) + (3.6276 * (WaterTemperature**2.))    &
                     - (0.043219 * (WaterTemperature**3.))
                     
                 Alfa = (600/fTS)**(1./2.)
                 
                 k = ((1.91 * Exp(0.35 * WindVelocity)) + ko) * Alfa 
                
            
        
            case default

                write(*,*)'Unknown Aeration method'
                stop 'AerationFlux_CO2 - ModuleFunctions - ERR01'

        end select
            
        AerationFlux_CO2 = k / (100. * 3600.)     !conversion from cm/h to m/s
        return
        
    end function AerationFlux_CO2 



    !------------------------------------------------------------------------

    subroutine ComputeAdvection1D(ilb, iub, dt, du, Prop, Q, V, ComputePoints,          &
                                  Ticoef, Ecoef, DCoef, Fcoef,                          &
                                  Method, TVD_Limitation, TetaExplicit, VolumeRelMax, Upwind2)

        !Arguments---------------------------------------------------
        real(8), dimension(:), intent(IN)  :: Q, V
        real,    dimension(:), intent(IN)  :: du, Prop
        integer, dimension(:), intent(IN)  :: ComputePoints
        real   ,               intent(IN)  :: dt, VolumeRelMax
        integer,               intent(IN)  :: ilb, iub, Method, TVD_Limitation
        real,                  intent(IN)  :: TetaExplicit !(0 - fully explicit, 1 - fully implicit) 
        logical,               intent(IN)  :: Upwind2

        real(8), dimension(:), intent(OUT) :: Ecoef
        real,    dimension(:), intent(OUT) :: Ticoef, DCoef, Fcoef
        !Local-------------------------------------------------------
        real(8), dimension(4)              :: V4
        real,    dimension(4)              :: CFace, Prop4, du4
        real(8)                            :: QFace
        real                               :: Aux, CrLeft, CrRight
        integer                            :: i
        logical                            :: NearBoundary

        !Begin-------------------------------------------------------


d1:     do i = ilb, iub

i1:         if (ComputePoints(i-1) == Compute .and. ComputePoints(i) == Compute) then

                QFace = Q(i)

                NearBoundary = .false.


                if (QFace > 0) then
                    if (ComputePoints(i-2) /= Compute) NearBoundary = .true.
                else
                    if (ComputePoints(i+1) /= Compute) NearBoundary = .true.
                endif

                Prop4(1) = Prop(i-2);Prop4(2) = Prop(i-1);Prop4(3) = Prop(i);Prop4(4) = Prop(i+1);
                du4  (1) = du  (i-2);du4  (2) = du  (i-1);du4  (3) = du  (i);du4  (4) = du  (i+1);
                V4   (1) = V   (i-2);V4   (2) = V   (i-1);V4   (3) = V   (i);V4   (4) = V   (i+1)

                call ComputeAdvectionFace(Prop4, V4, du4, dt, QFace, VolumeRelMax,         &
                                          Method, TVD_Limitation, NearBoundary, Upwind2, CFace)

                CrLeft  = QFace * dt / V4(2)
                CrRight = QFace * dt / V4(3)

                
                Aux           = CFace(1) * Prop(i-2) + CFace(2) * Prop(i-1) +           &
                                CFace(3) * Prop(i  ) + CFace(4) * Prop(i+1)

                TiCoef(i - 1) = TiCoef(i - 1) - (1 - TetaExplicit) * CrLeft  * Aux
                TiCoef(i    ) = TiCoef(i    ) + (1 - TetaExplicit) * CrRight * Aux

                if (TetaExplicit > 0. .and. (abs(CFace(1)) > 0. .or. abs(CFace(4)) > 0.)) then
                
                    stop 'ModuleFunctions - ComputeAdvection1D - ERR01'

                endif

                ECoef(i - 1) = ECoef(i - 1) + TetaExplicit * CrLeft  * CFace(2)
                FCoef(i - 1) = FCoef(i - 1) + TetaExplicit * CrLeft  * CFace(3)

                DCoef(i    ) = DCoef(i    ) - TetaExplicit * CrRight * CFace(2)
                ECoef(i    ) = ECoef(i    ) - TetaExplicit * CrRight * CFace(3)

            endif i1


        enddo d1

    end subroutine ComputeAdvection1D

    subroutine ComputeAdvection1D_V2(ilb, iub, dt, du, Prop, Q, V, ComputePoints,       &
                                     C_FLux, D_flux, E_flux, F_flux,                    &
                                     Method, TVD_Limitation, VolumeRelMax, Upwind2)

        !Arguments---------------------------------------------------
        real(8), dimension(:), intent(IN)  :: Q, V
        real,    dimension(:), intent(IN)  :: du, Prop
        integer, dimension(:), intent(IN)  :: ComputePoints
        real   ,               intent(IN)  :: dt, VolumeRelMax
        integer,               intent(IN)  :: ilb, iub, Method, TVD_Limitation
        logical,               intent(IN)  :: Upwind2

        real,    dimension(:), intent(OUT) :: C_flux, D_flux, E_flux, F_flux
        !Local-------------------------------------------------------
        real(8), dimension(4)              :: V4
        real,    dimension(4)              :: CFace, Prop4, du4
        real(8)                            :: QFace
        integer                            :: i
        logical                            :: NearBoundary

        !Begin-------------------------------------------------------


d1:     do i = ilb, iub

i1:         if (ComputePoints(i-1) == Compute .and. ComputePoints(i) == Compute) then

                QFace = Q(i)

                NearBoundary = .false.


                if (QFace > 0) then
                    if (ComputePoints(i-2) /= Compute) NearBoundary = .true.
                else
                    if (ComputePoints(i+1) /= Compute) NearBoundary = .true.
                endif

                Prop4(1) = Prop(i-2);Prop4(2) = Prop(i-1);Prop4(3) = Prop(i);Prop4(4) = Prop(i+1);
                du4  (1) = du  (i-2);du4  (2) = du  (i-1);du4  (3) = du  (i);du4  (4) = du  (i+1);
                V4   (1) = V   (i-2);V4   (2) = V   (i-1);V4   (3) = V   (i);V4   (4) = V   (i+1)

                call ComputeAdvectionFace(Prop4, V4, du4, dt, QFace, VolumeRelMax,         &
                                          Method, TVD_Limitation, NearBoundary, Upwind2, CFace)

                C_Flux(i    ) = QFace * CFace(1)
                D_Flux(i    ) = QFace * CFace(2)
                E_Flux(i    ) = QFace * CFace(3)
                F_Flux(i    ) = QFace * CFace(4)


            endif i1


        enddo d1
    

    end subroutine ComputeAdvection1D_V2



    subroutine ComputeAdvectionFace(Prop, V, du, dt, QFace, VolumeRelMax,                &
                                    Method, TVD_Limitation, NearBoundary, Upwind2, CFace)

        !Arguments---------------------------------------------------

        real(8), dimension(4), intent(IN)   :: V
        real,    dimension(4), intent(IN)   :: du, Prop
        real   ,               intent(IN)   :: dt, VolumeRelMax
        real(8),               intent(IN)   :: QFace
        integer,               intent(IN)   :: Method, TVD_Limitation
        logical              , intent(IN)   :: NearBoundary, Upwind2
        real,   dimension(4) , intent(OUT)  :: CFace

        !Local-------------------------------------------------------
        real                                :: Cr, Theta, dC, r, a, b, AuxLeft, AuxRight
        real,   dimension(4)                :: Cup1, CupHighOrder
        real(8)                             :: Aux, VolumeRel
                                            
        !Begin-------------------------------------------------------

        call FaceConcUpFirstOrder(Cup1, QFace)

        if (.not. NearBoundary) then
            if (QFace > 0) then
                Cr          = Courant (QFace,V(2),dt)  
                
                Aux         = min(V(1),V(2),V(3))
                VolumeRel   = max(V(1),V(2),V(3)) / Aux 
            else
                Cr  = Courant (QFace,V(3),dt)
                
                Aux         = min(V(2),V(3),V(4))
                VolumeRel   = max(V(2),V(3),V(4)) / Aux 
            endif
        endif

i2 :    if (Method == UpwindOrder1 .or. (NearBoundary.and.Upwind2)) then 
            
            CupHighOrder(1:4) = Cup1(1:4)
            Theta = 0.

            !if (QFace > 0) then
            !    AuxLeft  =  1.
            !    AuxRight =  0.
            !else
            !    AuxLeft  =  0.
            !    AuxRight =  1.
            !endif

        else if (Method == UpwindOrder2 .and. .not. NearBoundary) then i2

            call FaceConcUpSecondOrder(CupHighOrder, QFace)
            
            !1st order upwind is used in this case
            if (VolumeRel > VolumeRelMax) then
                Theta = 0. 
            else
                Theta = 1.
            endif
            
        else if (Method ==  UpwindOrder3 .and. .not. NearBoundary) then i2

            call FaceConcUpThirdOrder(CupHighOrder, QFace, Cr)
            !1st order upwind is used in this case
            if (VolumeRel > VolumeRelMax) then
                Theta = 0. 
            else
                Theta = 1.
            endif


        else if (Method == CentralDif .or. Method == LeapFrog) then i2

            AuxLeft  =  du(3  ) / (du(2) + du(3))
            AuxRight =  du(2  ) / (du(2) + du(3))

            CupHighOrder(1) = 0.
            CupHighOrder(2) = AuxLeft
            CupHighOrder(3) = AuxRight 
            CupHighOrder(4) = 0.

            Theta = 1.

        else if (Method == P2_TVD  .and. .not. NearBoundary) then i2

            CupHighOrder(1:4) = 0.

            if (QFace > 0.) then

                CupHighOrder(3) = 1.

                dC  = (Prop(3)-Prop(2)) / (du(3) + du(2))

                if (abs(dC)< MinValue ) then
                    if (dc>= 0) then
                        dc =   MinValue
                    else
                        dc = - MinValue
                    endif
                endif
             
                r = (Prop(2)-Prop(1))/ (du(2) + du(1))/ dC

            else

                CupHighOrder(2) = 1.

                dC  = (Prop(2)-Prop(3)) / (du(3) + du(2))

                if (abs(dC)< MinValue ) then
                    if (dc>= 0) then
                        dc =   MinValue
                    else
                        dc = - MinValue
                    endif
                endif
         
                r = (Prop(3)-Prop(4))/ (du(3) + du(4)) / dC
            endif

            
i5:         if      (TVD_Limitation == MinMod) then
                !minmod
                Theta = max(0.,min(1.,r))
            else if (TVD_Limitation == VanLeer) then
                !VanLeer
                if (r<0) then
                    Theta = 0.
                else
                    Theta = 2.*r/(1+r)
                endif

            else if (TVD_Limitation == Muscl) then
                !Muscl
                Theta = max(0.,min(2.,2.*r,(1+r)/2.))
            else if (TVD_Limitation == SuperBee) then

                !SuperBee
                Theta = max(0.,min(1.,2.*r),min(r, 2.))

            else if (TVD_Limitation == PDM) then

                a = 0.5 + (1 - 2.*abs(cr))/6.

                b = 0.5 - (1 - 2.*abs(cr))/6.

                Aux = a + b * r

                if (abs(Cr)< MinValue) Cr = MinValue
                !PDM
                Theta = max(0.,min(real(Aux),2./(1.-Cr),2.*r/Cr))

            else i5
                stop "This TVD Limitation option is not valid to compute Advection1D"
            endif i5

            Theta = 0.5 * Theta * (1. - Cr)


            !if (QFace > 0) then
            !    AuxLeft  =  1. - Theta
            !    AuxRight =  Theta
            !else
            !    AuxLeft  =  Theta
            !    AuxRight =  1. - Theta
            !endif

        else i2

            stop "This method is not valid to compute Advection1D"

        endif i2
        
   
        CFace(1) =  (1. - Theta) * Cup1(1) + Theta * CupHighOrder(1)
        CFace(2) =  (1. - Theta) * Cup1(2) + Theta * CupHighOrder(2)
        CFace(3) =  (1. - Theta) * Cup1(3) + Theta * CupHighOrder(3)
        CFace(4) =  (1. - Theta) * Cup1(4) + Theta * CupHighOrder(4)

    

    end subroutine ComputeAdvectionFace



    real function Courant (QFace,V,dt)

        !Arguments---------------------------------------------------
        real(8) :: QFace,V
        real    :: dt

        Courant = QFace * dt / V

    end function Courant

    Subroutine FaceConcUpSecondOrder (Coef, QFace)

        !Arguments---------------------------------------------------
        real(8)                        :: QFace
        real, dimension(4)             :: Coef
        !Local-------------------------------------------------------

        Coef(1:4) = 0.
        
        if(QFace>0) then
        
            !In a grid with constant spatial step
            Coef(1) = -1./8.
            Coef(2) =  6./8.
            Coef(3) =  3./8.


        else if(QFace<0) then

           !In a grid with constant spatial step
            Coef(4) = -1./8.
            Coef(3) =  6./8.
            Coef(2) =  3./8.

        else
            !Important for the land boundaries
            call FaceConcUpFirstOrder (Coef, QFace)
        endif

    end Subroutine FaceConcUpSecondOrder

    Subroutine FaceConcUpThirdOrder (Coef, QFace, Cr)

        !Arguments---------------------------------------------------
        real(8)                         :: QFace
        real, dimension(4)              :: Coef
        real                            :: Cr
        !Local-------------------------------------------------------
        real                            :: a, b, c, d
        !Begin-------------------------------------------------------

        c = (1 - 2.*abs(Cr))/6.

        a = 0.5 + c

        b = 0.5 - c

        d = (1 - abs(Cr))/2.
       
        if(QFace>0) then
        
            
            Coef(1) = - d * b
            Coef(2) = 1 + d * (b - a)
            Coef(3) =   d * a
            Coef(4) = 0.

        else if(QFace<0) then

            Coef(4) = - d * b
            Coef(3) = 1 + d * (b - a)
            Coef(2) =   d * a
            Coef(1) = 0.

        else
            !Important for the land boundaries
            call FaceConcUpFirstOrder (Coef, QFace)

        endif

    end Subroutine FaceConcUpThirdOrder

    Subroutine FaceConcUpFirstOrder (Coef, QFace)

        !Arguments---------------------------------------------------
        real, dimension(4)             :: Coef
        real(8)                        :: QFace


        Coef(1:4) = 0.                
        if (QFace>0.) then
            Coef(2) = 1.
        else
            Coef(3) = 1.
        endif           

    end Subroutine FaceConcUpFirstOrder

    real(4) function QuadraticInterpolationR4(y1, y2, y3, x1In, x2In, x3In, xIn)

        !Arguments---------------------------------------------------
        real(4) :: y1, y2, y3, x1In, x2In, x3In, xIn
        !Local-------------------------------------------------------
        real(4) :: x2, x3, x
        real(4) :: a, b, c

        x2 = x2In - x1In
        x3 = x3In - x1In
        x  =  xIn - x1In

        !y(x) = ax^2+bx+c
        c = y1
        b = (y3-y1) + (y1-y2)*x3*x3/x2/x2
        b = b / (x3 - x3*x3/x2)
        a = (y2 - y1 - b*x2) / x2 / x2
        
        QuadraticInterpolationR4 = a * x * x + b * x + c

    end function QuadraticInterpolationR4


    real(8) function QuadraticInterpolationR8(y1, y2, y3, x1In, x2In, x3In, xIn)

        !Arguments---------------------------------------------------
        real(8) :: y1, y2, y3, x1In, x2In, x3In, xIn
        !Local-------------------------------------------------------
        real(8) :: x2, x3, x
        real(8) :: a, b, c

        x2 = x2In - x1In
        x3 = x3In - x1In
        x  =  xIn - x1In


        !y(x) = ax^2+bx+c
        c = y1
        b = (y3-y1) + (y1-y2)*x3*x3/x2/x2
        b = b / (x3 - x3*x3/x2)
        a = (y2 - y1 - b*x2) / x2 / x2
        
        QuadraticInterpolationR8 = a * x * x + b * x + c

    end function QuadraticInterpolationR8

!
!Given arrays xa and ya, each of length n, and given a value x, this routine returns a
!value y, and an error estimate dy. If P(x) is the polynomial of degree N - 1 such that
!P(xai) = yai, i = 1, . . . , n, then the returned value y = P(x).


    subroutine polint(xa,ya,n,x,y,dy, STAT)
        !Arguments---------------------------------------------------
        integer, intent(IN) ::  n
        real(8), intent(IN) ::  x,xa(n),ya(n)
        real(8), intent(OUT)::  dy,y
        integer, intent(OUT), optional :: STAT
        
        !Largest anticipated value of n.
        integer, PARAMETER  :: NMAX=10
        !Local-------------------------------------------------------
        integer :: i,m,ns
        real    :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
        
        if(present(STAT))STAT = SUCCESS_

        ns=1
        dif=abs(x-xa(1))
d1:     do i=1,n ! Here we find the index ns of the closest table entry,
            dift=abs(x-xa(i))
            if (dift.lt.dif) then
                ns=i
                dif=dift
            endif
            c(i)=ya(i) !and initialize the tableau of c’s and d’s.
            d(i)=ya(i)
        enddo d1
        y=ya(ns) !This is the initial approximation to y.
        ns=ns-1
d3:     do m=1,n-1 !For each column of the tableau,
d2:         do i=1,n-m ! we loop over the current c’s and d’s and update them.
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
                if(abs(den)<1e-12) then
                    if (present(STAT)) then 
                        STAT = UNKNOWN_
                        exit
                    else
                        stop "failure in polint"
                    endif
                    
                endif 
                !This error can occur only if two input xa’s are (to within roundoff) identical.
                den=w/den
                d(i)=hp*den ! Here the c’s and d’s are updated.
                c(i)=ho*den
            enddo d2
            
            if(present(STAT))then
                if (STAT /= SUCCESS_) exit 
            endif
            
            if (2*ns.lt.n-m)then !After each column in the tableau is completed, we decide
                                 !which correction, c or d, we want to add to our accumulating
                                 !value of y, i.e., which path to take through
                                 !the tableau—forking up or down. We do this in such a
                                 !way as to take the most “straight line” route through the
                                 !tableau to its apex, updating ns accordingly to keep track
                                 !of where we are. This route keeps the partial approximations
                                 !centered (insofar as possible) on the target x. T he
                                 !last dy added is thus the error indication.
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            endif
            y=y+dy
        enddo d3

    end subroutine
    !End------------------------------------------------------------ 
    Subroutine ComputeDiffusion1D(ilb, iub, dt, du, Prop, k, v, ComputePoints, &
                                  Ticoef, Ecoef, DCoef, Fcoef, Theta)
        !Arguments---------------------------------------------------
        real(8), dimension(:), intent(IN)  :: v, k
        real,    dimension(:), intent(IN)  :: du, Prop 
        integer, dimension(:), intent(IN)  :: ComputePoints
        real   ,               intent(IN)  :: dt
        integer,               intent(IN)  :: ilb, iub
        real,                  intent(IN)  :: Theta
        real(8), dimension(:), intent(OUT) :: Ecoef
        real,    dimension(:), intent(OUT) :: Ticoef, DCoef, Fcoef
        !Locals---------------------------------------------------
        integer                            :: i
        real(8)                            :: aux, auxR, auxL
        !Begin---------------------------------------------------
        do i = ilb, iub
            if((ComputePoints(i)==Compute).and. &
               (ComputePoints(i-1)==Compute)) then                              
                aux = - k(i) * dt / du(i) !k = K * A !
                auxL = aux / v(i-1)
                auxR = aux / v(i) 
                DCoef(i)    = DCoef(i)   - Theta * auxR
                ECoef(i)    = ECoef(i)   + Theta * auxR
                ECoef(i-1)  = ECoef(i-1) + Theta * auxL
                FCoef(i-1)  = FCoef(i-1) - Theta * auxL
                TiCoef(i)   = TiCoef(i)   + (1 - Theta) * auxR * (Prop(i) - Prop(i-1))
                TiCoef(i-1) = TiCoef(i-1) - (1 - Theta) * auxL * (Prop(i) - Prop(i-1))
            endif
        enddo
    End Subroutine ComputeDiffusion1D
    !End------------------------------------------------------------
    Subroutine ComputeDiffusion3D(ilb, iub, jlb, jub, klb, kub, &
                                  dx, dy, dz,                   &
                                  cpx, cpy, cpz,                &
                                  kx, ky, kz,                   &
                                  dt, p, v,                     &
                                  Ticoef, Ecoef, DCoef, Fcoef,  &
                                  ThetaX, ThetaY, ThetaZ,       &
                                  CalcX, CalcY, CalcZ)
        !Arguments---------------------------------------------------
        integer,                   intent(IN)  :: ilb, iub,jlb, jub, klb, kub
        real, dimension(:,:), intent(IN)       :: dx, dy
        real, dimension(:,:,:), intent(IN)     :: dz    
        integer, dimension(:,:,:), intent(IN)  :: cpx, cpy, cpz    
        real(8), dimension(:,:,:), intent(IN)  :: kx, ky, kz    
        real(8), dimension(:,:,:), intent(IN)  :: v
        real,    dimension(:,:,:), intent(IN)  :: p
        real   ,                   intent(IN)  :: dt
        real,                      intent(IN)  :: ThetaX, ThetaY, ThetaZ
        real(8), dimension(:,:,:), intent(OUT) :: Ecoef
        real,    dimension(:,:,:), intent(OUT) :: Ticoef, DCoef, Fcoef
        logical                                :: CalcX, CalcY, CalcZ
        !Locals---------------------------------------------------
        integer                                 :: i, j, k
        real, pointer, dimension(:)             :: pint
        !Begin---------------------------------------------------
        !
        !----Allocate pint
        !
        allocate(pint(1:max(iub,jub,kub)+1))
        !----X diffusion
        if(CalcX) then 
            do k = klb+1, kub+1
            do i = ilb+1, iub+1
                !interpolate on faces
                do j = jlb+1, jub+1
                    if((cpx(i,j,k)==compute).and.(cpx(i,j-1,k)==compute)) then
                        pint(j) = (p(i,j,k) + p(i,j-1,k)) / 2    
                    endif
                enddo
                call ComputeDiffusion1D(jlb+1, jub+1, dt, dx(i,:), pint,                 &
                                        kx(i,:,k), v(i,:,k), cpx(i,:,k),                 &
                                        TiCoef(i,:,k), ECoef(i,:,k), DCoef(i,:,k),       &
                                        FCoef(i,:,k), ThetaX)
            enddo
            enddo
        endif
        !----Y diffusion
        if(CalcY) then 
            do k = klb+1, kub+1
            do j = jlb+1, jub+1
                !interpolate on faces
                do i = ilb+1, iub+1
                    if((cpy(i,j,k)==compute).and.(cpy(i-1,j,k)==compute)) then
                        pint(j) = (p(i-1,j,k) + p(i,j,k)) / 2    
                    endif
                enddo
                call ComputeDiffusion1D(ilb+1, iub+1, dt, dy(:,j), pint,                 &
                                        ky(:,j,k), v(:,j,k),  cpy(:,j,k),                &
                                        TiCoef(:,j,k), ECoef(:,j,k), DCoef(:,j,k),       &
                                        FCoef(:,j,k), ThetaY)
            enddo
            enddo
        endif
        !----Z diffusion
        if(CalcZ) then 
            do j = jlb+1, jub+1
            do i = ilb+1, iub+1
                !interpolate on faces
                do k = klb+1, kub+1
                    if((cpz(i,j,k)==compute).and.(cpz(i,j,k-1)==compute)) then
                        pint(j) = (p(i,j,k) + p(i,j,k-1)) / 2    
                    endif
                enddo
                call ComputeDiffusion1D(klb+1, kub+1, dt, dz(i,j,:), pint,               &
                                        kz(i,j,:), v(i,j,:), cpz(i,j,:),                 &
                                        TiCoef(i,j,:), ECoef(i,j,:), DCoef(i,j,:),       &
                                        FCoef(i,j,:), ThetaZ)
            enddo !do i
            enddo !do j
        endif
        deallocate(pint)
    End Subroutine ComputeDiffusion3D   
    !End------------------------------------------------------------
    Subroutine ComputeAdvection3D(ilb, iub, jlb, jub, klb, kub,     &
                                  dx, dy, dz,                       &
                                  cpx, cpy, cpz,                    &
                                  fx, fy, fz,                       &
                                  dt, p, v,                         &
                                  ti, e, d, f,                      &
                                  ThetaX, ThetaY, ThetaZ,           &
                                  CalcX, CalcY, CalcZ,              &
                                  methodH, methodV,                 &
                                  TVD_LimitationH, TVD_LimitationV, &
                                  VolumeRelMax, Upwind2H, Upwind2V)
        !Arguments---------------------------------------------------
        integer,                   intent(IN)  :: ilb, iub,jlb, jub, klb, kub
        real, dimension(:,:), intent(IN)       :: dx, dy
        real, dimension(:,:,:), intent(IN)     :: dz    
        integer, dimension(:,:,:), intent(IN)  :: cpx, cpy, cpz    
        real(8), dimension(:,:,:), intent(IN)  :: fx, fy, fz    
        real(8), dimension(:,:,:), intent(IN)  :: v
        real,    dimension(:,:,:), intent(IN)  :: p
        real   ,                   intent(IN)  :: dt
        real,                      intent(IN)  :: ThetaX, ThetaY, ThetaZ
        real(8), dimension(:,:,:), intent(OUT) :: e
        real,    dimension(:,:,:), intent(OUT) :: ti, d, f
                integer                        :: methodH, methodV 
        integer                                :: TVD_LimitationH, TVD_LimitationV
        real                                   :: VolumeRelMax
        logical                                :: Upwind2H, Upwind2V
        logical                                :: CalcX, CalcY, CalcZ
        !Locals---------------------------------------------------
        integer                                 :: i, j, k
        !Begin-----------------------------------------------------------
        !
        !----X advection
        if(CalcX) then 
            do k = klb+1, kub+1
            do i = ilb+1, iub+1
                call ComputeAdvection1D (jlb+1 , jub+1,                                                   &
                                dt, dx(i, :), p(i, :, k), fx(i, :, k), v(i, :, k), &
                                cpx(i,:,k), &
                                ti(i,:,k), e(i,:,k), d(i,:,k), f(i,:,k),                              &
                                MethodH, TVD_LimitationH, ThetaX, VolumeRelMax, Upwind2H)
            enddo
            enddo
        endif
        !----Y advection
        if(CalcY) then 
            do k = klb+1, kub+1
            do j = jlb+1, jub+1
                call ComputeAdvection1D (ilb+1 , iub+1,                                                       &
                                dt, dy(:, j), p(:, j, k), fy(:, j, k), v(:, j, k), &
                                cpy(:, j, k), &
                                ti(:, j, k), e(:, j, k), d(:, j, k), f(:, j, k),                              &
                                MethodH, TVD_LimitationH, ThetaY, VolumeRelMax, Upwind2H)
            enddo
            enddo
        endif
        !----Z advection
        if(CalcZ) then 
            do j = jlb+1, jub+1
            do i = ilb+1, iub+1
                call ComputeAdvection1D (klb +1, kub+1,                                                            &
                                dt, dz(i, j, :), p(i, j, :), fz(i, j, :),v(i, j, :), &
                                cpz(i, j, :), &
                                ti(i, j, :), e(i, j, :), d(i, j, :), f(i, j, :),       &
                                MethodV, TVD_LimitationV, ThetaZ , VolumeRelMax, Upwind2V)
            enddo
            enddo
        endif
    EndSubroutine
    !End------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !    
    ! Interpolate3D
    !
    ! This routine makes interpolation of a value in a 3D grid.
    ! For the moment only a tri-linear scheme is used.
    ! This can be improved, but I cannot see any reason for it.
    ! It reduces to bilinear or simple interpolation if one or more 
    ! interpolation coefficients are equal to zero.
    ! The true improvement would be to avoid unecessary operations for
    ! those cases.
    !
    ! This routine gives us the interpolated values at 26 different 
    ! locations on the faces of a grid. We should be very happy.
    !
    ! Note that dx, dy and dz are distances BETWEEN FACES measured trough cell centers
    ! fi, fj, fk are cell coordinates 
    ! (0,0,0) represents the center
    ! a +/- 1 represents a displacement to one of the faces.
    ! e.g. : (1,1,1) is upper corner
    !        (0,1,0) is center of x positive face of the cell
    !        (-1,-1,-1) is the lower corner
    !                                                            
    ! Author: Hernani Theias (05/2004)                                                       !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(4) Function interpolate3D_R4(p, dx, dy, dz, i, j, k, fi, fj, fk)
    !Arguments-----------------------------------------------------------------------------
        real(4), dimension(:, :, :), pointer    :: p            ! property that we interpolate
        real,    dimension(:, :   ), pointer    :: dx, dy       ! XY distance between faces through centers
        real,    dimension(:, :, :), pointer    :: dz           ! Z distance between faces in the vertical (also through centers)
        integer                                 :: i, j, k      ! cell at which interpolation is made
        integer                                 :: fi, fj, fk   ! cell coordinates (see above) of the interpolation location
        !Local---------------------------------------------------------------------------------
        real(8)                                 :: eta, xsi, zeta !interpolation coefficients
        real(8)                                 :: dz1, dz2     ! z spacing at interpolation location
        !Begin---------------------------------------------------------------------------------
        if((fi==0).and.(fj==0).and.(fk==0))then !save some time if there is no interpolation
            interpolate3D_R4 = p(i, j, k)
            return
        endif 
        ! eta = x / DX
        if((fj == 0).or.(abs(p(i, j + fj, k)/FillValueReal)>0.01) &
                    .or.(abs(p(i + fi, j + fj, k)/FillValueReal)>0.01)) then
            eta = 0 !save some time here
        else
            eta  = (dx(i, j) + dx(i + fi, j))                           &
                   / (                                                  &          
                         dx(i, j) + dx(i + fi, j)                       &
                       + dx(i, j + fj) + dx(i + fi, j + fj)             &
                     )
        endif
        ! xsi = y / DY
        if((fi == 0).or.(abs(p(i + fi, j, k)/FillValueReal)>0.01) &
                    .or.(abs(p(i + fi, j + fj, k)/FillValueReal)>0.01)) then
            xsi = 0 !again...
        else
            xsi  = (dy(i, j) + dy(i, j + fj))                           &
                   / (                                                  &          
                         dy(i, j) + dy(i, j + fj)                       &
                       + dy(i + fi, j) + dy(i + fi, j + fj)             &
                     )
        endif
        ! zeta = z / DZ
        if((fk == 0) .or.(abs(p(i + fi, j, k)/FillValueReal)>0.01)          &
                    .or.(abs(p(i + fi, j, k +fk)/FillValueReal)>0.01)      &
                    .or.(abs(p(i , j + fj, k)/FillValueReal)>0.01)         &
                    .or.(abs(p(i, j + fj, k + fk)/FillValueReal)>0.01)     &
                    .or.(abs(p(i + fi, j + fj, k)/FillValueReal)>0.01)     &
                    .or.(abs(p(i + fi, j + fj, k +fk)/FillValueReal)>0.01) &
                    .or.(abs(p(i, j, k +fk)/FillValueReal)>0.01)) then
           zeta = 0
        else
           dz1 =   (1 - eta) * (1 - xsi) * dz(i     , j     , k     )     &
                 + (    eta) * (1 - xsi) * dz(i     , j + fj, k     )     &
                 + (1 - eta) * (    xsi) * dz(i + fi, j     , k     )     &
                 + (    eta) * (    xsi) * dz(i + fi, j + fj, k     )
           dz2 =   (1 - eta) * (1 - xsi) * dz(i     , j     , k + fk)     &
                 + (    eta) * (1 - xsi) * dz(i     , j + fj, k + fk)     &
                 + (1 - eta) * (    xsi) * dz(i + fi, j     , k + fk)     &
                 + (    eta) * (    xsi) * dz(i + fi, j + fj, k + fk)
           zeta = dz1 / (dz1 + dz2) 
        endif 
        !interpolate (tri-linear)
        interpolate3D_R4 =   (1 - xsi) * (1 - eta) * (1 - zeta) * p(i     , j     , k     )  &
                           + (    xsi) * (1 - eta) * (1 - zeta) * p(i + fi, j     , k     )  &
                           + (1 - xsi) * (    eta) * (1 - zeta) * p(i     , j + fj, k     )  &
                           + (1 - xsi) * (1 - eta) * (    zeta) * p(i     , j     , k + fk)  &
                           + (    xsi) * (    eta) * (1 - zeta) * p(i + fi, j + fj, k     )  &
                           + (    xsi) * (1 - eta) * (    zeta) * p(i + fi, j     , k + fk)  &
                           + (1 - xsi) * (    eta) * (    zeta) * p(i     , j + fj, k + fk)  &
                           + (    xsi) * (    eta) * (    zeta) * p(i + fi, j + fj, k + fk)
    end function interpolate3D_R4
    !End------------------------------------------------------------
    real(8) function interpolate3D_R8(p, dx, dy, dz, i, j, k, fi, fj, fk)
    !Arguments-----------------------------------------------------------------------------
        real(8), dimension(:, :, :), pointer    :: p            ! property that we interpolate
        real,    dimension(:, :), pointer       :: dx, dy       ! XY distance between faces through centers
        real,    dimension(:, :, :), pointer    :: dz           ! Z distance between faces in the vertical (also through centers)
        integer                                 :: i, j, k      ! cell at which interpolation is made
        integer                                 :: fi, fj, fk   ! cell coordinates (see above) of the interpolation location
        !Local---------------------------------------------------------------------------------
        real(8)                                 :: eta, xsi, zeta !interpolation coefficients
        real(8)                                 :: dz1, dz2     ! z spacing at interpolation location
        real, parameter                         :: FillValueReal=-9.9e15
        !Begin---------------------------------------------------------------------------------
        if((fi==0).and.(fj==0).and.(fk==0))then !save some time if there is no interpolation
            interpolate3D_R8 = p(i, j, k)
            return
        endif 
        ! eta = x / DX
        if((fj == 0).or.(abs(p(i, j + fj, k)/FillValueReal)>0.01) &
                    .or.(abs(p(i + fi, j + fj, k)/FillValueReal)>0.01)) then
            eta = 0 !save some time here
        else
            eta  = (dx(i, j) + dx(i + fi, j))                           &
                   / (                                                  &          
                         dx(i, j) + dx(i + fi, j)                       &
                       + dx(i, j + fj) + dx(i + fi, j + fj)             &
                     )
        endif
        ! xsi = y / DY
        if((fi == 0).or.(abs(p(i + fi, j, k)/FillValueReal)>0.01) &
                    .or.(abs(p(i + fi, j + fj, k)/FillValueReal)>0.01)) then
            xsi = 0 !again...
        else
            xsi  = (dy(i, j) + dy(i, j + fj))                           &
                   / (                                                  &          
                         dy(i, j) + dy(i, j + fj)                       &
                       + dy(i + fi, j) + dy(i + fi, j + fj)             &
                     )
        endif
        ! zeta = z / DZ
        if((fk == 0) .or.(abs(p(i + fi, j, k)/FillValueReal)>0.01)          &
                    .or.(abs(p(i + fi, j, k +fk)/FillValueReal)>0.01)      &
                    .or.(abs(p(i , j + fj, k)/FillValueReal)>0.01)         &
                    .or.(abs(p(i, j + fj, k + fk)/FillValueReal)>0.01)     &
                    .or.(abs(p(i + fi, j + fj, k)/FillValueReal)>0.01)     &
                    .or.(abs(p(i + fi, j + fj, k +fk)/FillValueReal)>0.01) &
                    .or.(abs(p(i, j, k +fk)/FillValueReal)>0.01)) then
           zeta = 0
        else
           dz1 =   (1 - eta) * (1 - xsi) * dz(i     , j     , k     )     &
                 + (    eta) * (1 - xsi) * dz(i     , j + fj, k     )     &
                 + (1 - eta) * (    xsi) * dz(i + fi, j     , k     )     &
                 + (    eta) * (    xsi) * dz(i + fi, j + fj, k     )
           dz2 =   (1 - eta) * (1 - xsi) * dz(i     , j     , k + fk)     &
                 + (    eta) * (1 - xsi) * dz(i     , j + fj, k + fk)     &
                 + (1 - eta) * (    xsi) * dz(i + fi, j     , k + fk)     &
                 + (    eta) * (    xsi) * dz(i + fi, j + fj, k + fk)
           zeta = dz1 / (dz1 + dz2) 
        endif 
        !interpolate (tri-linear)
        interpolate3D_R8 =   (1 - xsi) * (1 - eta) * (1 - zeta) * p(i     , j     , k     )  &
                           + (    xsi) * (1 - eta) * (1 - zeta) * p(i + fi, j     , k     )  &
                           + (1 - xsi) * (    eta) * (1 - zeta) * p(i     , j + fj, k     )  &
                           + (1 - xsi) * (1 - eta) * (    zeta) * p(i     , j     , k + fk)  &
                           + (    xsi) * (    eta) * (1 - zeta) * p(i + fi, j + fj, k     )  &
                           + (    xsi) * (1 - eta) * (    zeta) * p(i + fi, j     , k + fk)  &
                           + (1 - xsi) * (    eta) * (    zeta) * p(i     , j + fj, k + fk)  &
                           + (    xsi) * (    eta) * (    zeta) * p(i + fi, j + fj, k + fk)
    end function interpolate3D_R8
    !End------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !    
    ! CGS2D
    !
    ! Pre-Conditionned Conjugate Gradient Solver for symmetric 
    ! matrices (e.g., pressure or pressure correction equation, heat conduction, etc.)
    !
    ! Taken from 'Computational Methods for Fluid Dynamics'
    ! by J.H. Ferziger and M.Peric
    ! Original code by I.Demirdzic - Masinski Fakultet, Sarajevo (1987)
    ! 
    ! Adapted for MOHID : Hernani Theias (05/04)                                                       !
    !                                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine CGS2D(fi, ap, as, an, aw, ae, q,        &
                     ilb, iub, jlb, jub, klb, kub, &
                     dIi, dIj, dIk,                &
                     dJi, dJj, dJk,                &
                     dKi, dKj, dKk,                &
                     maxit, resmax, mapping, normalized)
        !Arguments----------------------------------------------------------
        real, dimension(:,:,:), pointer             :: fi !solution vector
        real, dimension(:,:,:), pointer             :: q !independent term
        real(8), dimension(:,:,:), pointer          :: ap !diagonal coefficents
        real, dimension(:,:,:), pointer             :: an, as !coefficients corresponding to dJ
        real, dimension(:,:,:), pointer             :: ae, aw !coefficients corresponding to dI
        integer                                     :: ilb, iub, jlb, jub, klb, kub !bounds
        integer                                     :: maxit !max number of iterations
        integer                                     :: dIi, dIj, dIk !index change in W-E direction
        integer                                     :: dJi, dJj, dJk !index change in N-S direction
        integer                                     :: dKi, dKj, dKk !index change in B-T direction
        real                                        :: resmax !maximum value allowed for residual
        integer, dimension(:,:,:), pointer          :: mapping !cells mapping
        logical                                     :: normalized !specifies if abs or norm residuals must 
                                                                  !be used in convergence test
        !Local--------------------------------------------------------------
        integer                                     :: IL, IU, JL, JU, KL, KU, II, JJ, KK
        integer                                     :: i, j, k, n
        real                                        :: res0, resn, rsm
        real                                        :: s0, sk, s2, alf, bet
        real, dimension(:,:,:), pointer             :: res
        real(8), dimension(:,:,:), pointer          :: d, zk, pk
        logical                                     :: EndIterate
        !Begin--------------------------------------------------------------
        !
        !----Allocate local matrixes
        !
        allocate(res(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(d(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(zk(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(pk(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        !
        !----Calculate index bounds
        IL = ilb * dIi + jlb * dIj + klb * dIk
        IU = iub * dIi + jub * dIj + kub * dIk
        JL = ilb * dJi + jlb * dJj + klb * dJk
        JU = iub * dJi + jub * dJj + kub * dJk
        KL = ilb * dKi + jlb * dKj + klb * dKk
        KU = iub * dKi + jub * dKj + kub * dKk
        !----Solve for each K (in the non-solving direction)
        do KK = KL, KU
        !--------Calculate intial residual vector 
            res0 = 0.0
            do II = IL, IU
            do JJ = JL, JU
                i = II * dIi + JJ * dJi + KK * dKi
                j = II * dIj + JJ * dJj + KK * dKj
                k = II * dIk + JJ * dJk + KK * dKk
                if(mapping(i, j, k) == OpenPoint) then
                    res(i, j , k) =   q(i, j, k)                                                &
                                    - ap(i, j, k) * fi(i, j, k)                                 &
                                    - aw(i, j, k) * fi(i - dIi, j - dIj, k - dIk)               &
                                    - ae(i, j, k) * fi(i + dIi, j + dIj, k + dIk)               &
                                    - as(i, j, k) * fi(i - dJi, j - dJj, k - dJk)               &
                                    - an(i, j, k) * fi(i + dJi, j + dJj, k + dJk)
                    res0 = res0 + abs(res(i, j, k))
                endif
            enddo
            enddo
        !--------Pre-Conditioning Matrix Diagonal
            do II = IL, IU
            do JJ = JL, JU
                i = II * dIi + JJ * dJi + KK * dKi
                j = II * dIj + JJ * dJj + KK * dKj
                k = II * dIk + JJ * dJk + KK * dKk
                if(mapping(i, j, k) == OpenPoint) then
                    d(i, j, k) = 1.0 / ap(i, j, k)                                             &
                                  - aw(i, j, k) ** 2 * d(i - dIi, j - dIj, k - dIk)             &
                                  - as(i, j, k) ** 2 * d(i - dJi, j - dJj, k - dJk) 
                endif
            enddo
            enddo
        !--------Calculation of Zk, inner iteration loop
            s0 = 1.e20
            n = 0
            EndIterate = .false.
            do while((n < maxit).and.(.not.EndIterate))
                n = n + 1
        !------------Forward substitution
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        zk(i, j, k) =   d(i, j, k)                                              &
                                      * (                                                       &
                                           res(i ,j , k)                                        &
                                         - aw(i, j, k) * zk(i - dIi, j - dIj, k - dIk)          &
                                         - as(i, j, k) * zk(i - dJi, j - dJj, k - dJk)          &
                                        )
                    endif
                enddo
                enddo
                !
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        zk(i, j, k) =   zk(i, j, k) / (d(i, j, k) + 1.e-20)
                    endif
                enddo
                enddo
        !------------Backward substitution
                sk = 0.0
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        zk(i, j, k) =   d(i, j, k)                                              &
                                      * (                                                       &
                                           zk(i ,j , k)                                         &
                                         - ae(i, j, k) * zk(i + dIi, j + dIj, k + dIk)          &
                                         - an(i, j, k) * zk(i + dJi, j + dJj, k + dJk)          &
                                        )
                        sk = sk + res(i, j, k) * zk(i, j ,k)
                    endif
                enddo
                enddo
        !------------Calculate beta and new search vector
                bet = sk / s0
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        pk(i, j, k) = zk(i, j, k) + bet * pk(i, j, k)
                    endif
                enddo
                enddo
        !------------Calculate scalar product (pk. A pk) and alpha
                s2 = 0.0
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        zk(i, j, k) =  ap(i, j, k) * pk(i, j, k)                                &
                                     + aw(i, j, k) * pk(i - dIi, j - dIj, k - dIk)              &
                                     + ae(i, j, k) * pk(i + dIi, j + dIj, k + dIk)              &
                                     + an(i, j, k) * pk(i + dJi, j + dJj, k + dJk)              &
                                     + as(i, j, k) * pk(i - dJi, j - dJj, k - dJk)
                        s2 = s2 + pk(i, j, k) * zk(i, j, k)
                    endif
                enddo
                enddo
        ! 
                alf = sk / (s2 + 1.e-20)
        !
        !------------Calculate new residual and update variable
                resn = 0.0
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        fi(i, j, k) = fi(i, j, k) + alf * pk(i, j, k)
                        res(i, j, k) = res(i, j, k) - alf * zk(i, j, k)
                        resn = resn + abs(res(i, j, k))
                    endif
                enddo
                enddo
        !------------Check convergence
                if(normalized) then
                    rsm = resn / (res0 + 1.e-20)
                else
                    rsm = resn
                endif
                if(rsm.LT.resmax) EndIterate = .true.
            enddo !while - inner loop
            if(n>=maxit) then
                write(*,*) 'Solution not converging.'
                stop
            endif
        enddo ! do kk
        !----De-Allocate local matrixes
        !
        deallocate(res)
        nullify(res)
        deallocate(d)
        nullify(d)
        deallocate(zk)
        nullify(zk)
        deallocate(pk)
        nullify(pk)
        !
        return            
    End Subroutine ! CSG2D
    !End------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                          
    ! LISOLVE
    !
    ! TDMA algortihm applied line-by-line alterning along J and I lines
    !
    ! ----------------------------
    !
    ! author : Hernani Theias
    ! last modified  : 06/2004 
    ! e-mail : hernanitheias@netcabo.pt     
    !                                                                                      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine LISOLVE( fi, ap, as, an, aw, ae, q,          &
                        ilb, iub, jlb, jub, klb, kub,       &
                        dIi, dIj, dIk,                      &
                        dJi, dJj, dJk,                      &
                        dKi, dKj, dKk,                      &
                        maxit, resmax, mapping)
        !Arguments----------------------------------------------------------
        real, dimension(:,:,:), pointer             :: fi !solution vector
        real, dimension(:,:,:), pointer             :: q !independent term
        real(8), dimension(:,:,:), pointer          :: ap !diagonal coefficents
        real, dimension(:,:,:), pointer             :: an, as !coefficients corresponding to dJ
        real, dimension(:,:,:), pointer             :: ae, aw !coefficients corresponding to dI
        integer                                     :: ilb, iub, jlb, jub, klb, kub !bounds
        integer                                     :: maxit !max number of iterations
        integer                                     :: dIi, dIj, dIk !index change in W-E direction
        integer                                     :: dJi, dJj, dJk !index change in N-S direction
        integer                                     :: dKi, dKj, dKk !index change in B-T direction
        real                                        :: resmax !maximum value allowed for residual
        integer, dimension(:,:,:), pointer          :: mapping !cells mapping
        !Local--------------------------------------------------------------
        integer                                     :: IL, IU, JL, JU, KL, KU
        integer                                     :: II, JJ, KK
        integer                                     :: i, j, k, n
        real                                        :: res0, resn, rsm, res
        real(8), dimension(:), pointer              :: vecg, vecw
        logical                                     :: EndIterate
        !Begin--------------------------------------------------------------
        !
        !----Calculate index bounds
        IL = ilb * dIi + jlb * dIj + klb * dIk
        IU = iub * dIi + jub * dIj + kub * dIk
        JL = ilb * dJi + jlb * dJj + klb * dJk
        JU = iub * dJi + jub * dJj + kub * dJk
        KL = ilb * dKi + jlb * dKj + klb * dKk
        KU = iub * dKi + jub * dKj + kub * dKk
        !----Allocate and initialize local matrixes
        allocate(vecg(min(IL-1,JL-1):max(IU+1,JU+1)))
        allocate(vecw(min(IL-1,JL-1):max(IU+1,JU+1)))
        vecg(:) = 0.0
        vecw(:) = 0.0
        !--------Solve for each K (in the non-solving direction)
        do KK = KL, KU
        !--------Inner loop for each IJ plane
            n = 0
            EndIterate = .false.
            res0 = 0.0
            do while((n < maxit).and.(.not.EndIterate))
        !-------------Solve with 'TMDA' along J-lines
                n = n + 1
                do II = IL, IU
                    i = II * dIi + JL * dJi + KK * dKi
                    j = II * dIj + JL * dJj + KK * dKj
                    k = II * dIk + JL * dJk + KK * dKk
                    vecw(JL) = - an(i, j, k) / ap(i, j, k)
                    vecg(JL) = - q(i, j, k) / ap(i, j, k)
                    do JJ = JL + 1, JU
                        i = II * dIi + JJ * dJi + KK * dKi
                        j = II * dIj + JJ * dJj + KK * dKj
                        k = II * dIk + JJ * dJk + KK * dKk
                        vecw(JJ) = - an(i, j, k) / (ap(i, j, k) + as(i, j, k) * vecw(JJ - 1))
                        vecg(JJ) =   (q(i, j, k) - as(i, j, k) * vecg(JJ - 1))              &
                                   / (ap(i, j, k) + as(i, j, k) * vecw(JJ - 1))
                    enddo !JJ
                    i = II * dIi + JU * dJi + KK * dKi
                    j = II * dIj + JU * dJj + KK * dKj
                    k = II * dIk + JU * dJk + KK * dKk
                    fi(i, j, k) = vecg(JU)
                    do JJ = JU - 1, JL + 1, -1
                        i = II * dIi + JJ * dJi + KK * dKi
                        j = II * dIj + JJ * dJj + KK * dKj
                        k = II * dIk + JJ * dJk + KK * dKk
                        fi(i, j, k) = vecw(JJ) * fi(i + dJi, j + dJj, k + dJk) + vecg(JJ)
                    enddo !JJ
                enddo !II
        !------------Solve with 'TMDA' along I-lines
                do JJ = JL, JU
                    i = IL * dIi + JJ * dJi + KK * dKi
                    j = IL * dIj + JJ * dJj + KK * dKj
                    k = IL * dIk + JJ * dJk + KK * dKk
                    vecw(IL) = - an(i, j, k) / ap(i, j, k)
                    vecg(IL) = - q(i, j, k) / ap(i, j, k)
                    do II = IL + 1, IU
                        i = II * dIi + JJ * dJi + KK * dKi
                        j = II * dIj + JJ * dJj + KK * dKj
                        k = II * dIk + JJ * dJk + KK * dKk
                        vecw(II) = - an(i, j, k) / (ap(i, j, k) + as(i, j, k) * vecw(II - 1))
                        vecg(II) =   (q(i, j, k) - as(i, j, k) * vecg(II - 1))              &
                                   / (ap(i, j, k) + as(i, j, k) * vecw(II - 1))
                    enddo !JJ
                    i = IU * dIi + JJ * dJi + KK * dKi
                    j = IU * dIj + JJ * dJj + KK * dKj
                    k = IU * dIk + JJ * dJk + KK * dKk
                    fi(i, j, k) = vecg(IU)
                    do II = IU - 1, IL + 1, -1
                        i = II * dIi + JJ * dJi + KK * dKi
                        j = II * dIj + JJ * dJj + KK * dKj
                        k = II * dIk + JJ * dJk + KK * dKk
                        fi(i, j, k) = vecw(II) * fi(i + dIi, j + dIj, k + dIk) + vecg(II)
                    enddo !JJ
                enddo !II
        !--------Calculate new residual and update variable
                resn = 0.0
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        res  =   q(i, j, k)                                                &
                               - ap(i, j, k) * fi(i, j, k)                                 &
                               - aw(i, j, k) * fi(i - dIi, j - dIj, k - dIk)               &
                               - ae(i, j, k) * fi(i + dIi, j + dIj, k + dIk)               &
                               - as(i, j, k) * fi(i - dJi, j - dJj, k - dJk)               &
                               - an(i, j, k) * fi(i + dJi, j + dJj, k + dJk)
                        resn = resn + abs(res)
                    endif
                enddo
                enddo
        !------------Check convergence
                if(n == 1) res0 = resn
                rsm = resn / (res0 + 1.e-20)
                if(rsm.LT.resmax) EndIterate = .true.
            enddo !while - inner loop
            if(n>=maxit) then
                write(*,*) 'Solution not converging.'
                stop
            endif
        enddo !KK
        !----De-Allocate local matrixes
        !
        deallocate(vecg)
        nullify(vecg)
        deallocate(vecw)
        nullify(vecw)
        !      
    End Subroutine
    !End------------------------------------------------------------

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                      !    
    ! SIP_SOL
    !
    ! This is the ILU solver after Stone.
    !
    ! Taken from 'Computational Methods for Fluid Dynamics'
    ! by J.H. Ferziger and M.Peric, Institut fuer Schiffbau, Hamburg, 1995
    ! Original code by M. Peric
    ! 
    ! Adapted for MOHID : Hernani Theias (05/04)                                                       !
    !                                                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine SIPSOL(  fi, ap, as, an, aw, ae, q,          &
                        ilb, iub, jlb, jub, klb, kub,       &
                        dIi, dIj, dIk,                      &
                        dJi, dJj, dJk,                      &
                        dKi, dKj, dKk,                      &
                        alpha,                              &
                        maxit, resmax, mapping)
        !Arguments----------------------------------------------------------
        real, dimension(:,:,:), pointer             :: fi !solution vector
        real, dimension(:,:,:), pointer             :: q !independent term
        real(8), dimension(:,:,:), pointer          :: ap !diagonal coefficents
        real, dimension(:,:,:), pointer             :: an, as !coefficients corresponding to dJ
        real, dimension(:,:,:), pointer             :: ae, aw !coefficients corresponding to dI
        integer                                     :: ilb, iub, jlb, jub, klb, kub !bounds
        integer                                     :: maxit !max number of iterations
        integer                                     :: dIi, dIj, dIk !index change in W-E direction
        integer                                     :: dJi, dJj, dJk !index change in N-S direction
        integer                                     :: dKi, dKj, dKk !index change in B-T direction
        real                                        :: resmax !maximum value allowed for residual
        real                                        :: alpha
        integer, dimension(:,:,:), pointer          :: mapping
        !Local--------------------------------------------------------------
        integer                                     :: IL, IU, JL, JU, KL, KU, II, JJ, KK
        integer                                     :: i, j, k, n
        real                                        :: res0, resn, rsm
        real, dimension(:,:,:), pointer             :: res
        real(8), dimension(:,:,:), pointer          :: lw, ls, un, ue, lpr
        logical                                     :: EndIterate
        !Begin--------------------------------------------------------------
        !
        !----Allocate local matrixes
        allocate(lw(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(ls(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(un(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(ue(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(lpr(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        allocate(res(ilb-1:iub+1, jlb-1:jub+1, klb-1:kub+1))
        lw(:,:,:) = 0.
        ls(:,:,:) = 0.
        un(:,:,:) = 0.
        ue(:,:,:) = 0.
        lpr(:,:,:) = 1.
        res(:,:,:) = 0.
        !----Calculate index bounds
        IL = ilb * dIi + jlb * dIj + klb * dIk
        IU = iub * dIi + jub * dIj + kub * dIk
        JL = ilb * dJi + jlb * dJj + klb * dJk
        JU = iub * dJi + jub * dJj + kub * dJk
        KL = ilb * dKi + jlb * dKj + klb * dKk
        KU = iub * dKi + jub * dKj + kub * dKk
        !----Solve for each K (in the non-solving direction)
        do KK = KL, KU
        !--------Calculate elements of [L] ans [U] matrices 
            res0 = 0.0
            do II = IL, IU
            do JJ = JL, JU
                i = II * dIi + JJ * dJi + KK * dKi
                j = II * dIj + JJ * dJj + KK * dKj
                k = II * dIk + JJ * dJk + KK * dKk
                lw(i, j, k) = aw(i, j, k) / (1 + alpha * un(i - dJi, j - dJj, k - dJk)) !north
                ls(i, j, k) = as(i, j, k) / (1 + alpha * ue(i - dIi, j - dIj, k - dIk)) !north
                lpr(i, j, k) =   ap(i, j, k)                                           &
                               + alpha * lw(i, j, k) * un(i - dJi, j - dJj, k - dJk)   &
                               + alpha * ls(i, j, k) * ue(i - dIi, j - dIj, k - dIk)   &
                               -         lw(i, j, k) * ue(i - dJi, j - dJj, k - dJk)   &
                               -         ls(i, j, k) * un(i - dIi, j - dIj, k - dIk)
                un(i, j, k) = (   an(i, j, k)                                                 &
                                - alpha * lw(i, j, k) * un(i - dJi, j - dJj, k - dJk)         &
                              ) / lpr(i, j, k)
                ue(i, j, k) = (   ae(i, j, k)                                                 &
                                - alpha * ls(i, j, k) * ue(i - dIi, j - dIj, k - dIk)         &
                              ) / lpr(i, j, k)
            enddo
            enddo
        !--------Inner iteration loop.
            n = 0
            EndIterate = .false.
            do while((n < maxit).and.(.not.EndIterate))
                n = n + 1
        !------------Calculate residual
                resn = 0.
                do II = IL, IU
                do JJ = JL, JU
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        res(i, j , k) =   q(i, j, k)                                                &
                                        - ap(i, j, k) * fi(i, j, k)                                 &
                                        - aw(i, j, k) * fi(i - dIi, j - dIj, k - dIk)               &
                                        - ae(i, j, k) * fi(i + dIi, j + dIj, k + dIk)               &
                                        - as(i, j, k) * fi(i - dJi, j - dJj, k - dJk)               &
                                        - an(i, j, k) * fi(i + dJi, j + dJj, k + dJk)
                        resn = resn + abs(res(i, j, k))
                        res(i, j, k) = (   res(i, j ,k)                                             &
                                         - ls(i, j, k) * res(i - dIi, j - dIj, k - dIk)             &
                                         - lw(i, j, k) * res(i - dJi, j - dJj, k - dJk)             &
                                        ) / lpr(i, j, k)

                    endif
                enddo
                enddo
                if(n==1) res0 = resn
        !------------Calculate increment and correct variable (backwards)
                do II = IU, IL, -1
                do JJ = JU, JL, -1
                    i = II * dIi + JJ * dJi + KK * dKi
                    j = II * dIj + JJ * dJj + KK * dKj
                    k = II * dIk + JJ * dJk + KK * dKk
                    if(mapping(i, j, k) == OpenPoint) then
                        res(i, j, k) =   res(i, j, k)                                                   &
                                       - un(i, j, k) * res(i + dIi, j + dIj, k + dIk)                   &
                                       - ue(i, j, k) * res(i + dJi, j + dJj, k + dJk)
                        fi(i, j, k) = fi(i, j, k) + res(i, j, k)
                    endif
                enddo
                enddo
        !-----------Check convergence
                if(n == 1) res0 = resn
                rsm = resn / (res0 + 1.e-20)
                if(rsm.LT.resmax) EndIterate = .true.
            enddo
            if(n>=maxit) then
                write(*,*) 'Solution not converging in the specified number of iterations.'
                stop
            endif
        enddo !do K
        !----De-Allocate local matrixes
        !
        deallocate(res)
        nullify(res)
        deallocate(ls)
        nullify(ls)
        deallocate(lw)
        nullify(lw)
        deallocate(lpr)
        nullify(lpr)
        deallocate(un)
        nullify(un)
        deallocate(ue)
        nullify(ue)
        !      
    End Subroutine ! SIP_SOLV
    !End------------------------------------------------------------
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   
    ! BICGSTAB2D
    !
    ! ----------------------------
    !
    ! Bi-Conjugate Gradient Stabilized Algorithm.
    ! See :
    ! G. Sleijpen, H. Van der Vorst, 'Hybrid Bi-Conjugate Gradient Methods for CFD Problems' 
    !
    ! author : Hernani Theias
    ! last modified  : 06/2004 
    ! e-mail : hernanitheias@netcabo.pt     
    !
    !                                                                                      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Subroutine BICGSTAB2D(fi, cp, cs, cn, cw, ce, qq,       &
                          ilb, iub, jlb, jub, klb, kub,     &
                          normal,                           &
                          maxit, resmax, mapping, alphaLU, normalized)
        !Arguments----------------------------------------------------------
        real, dimension(:,:,:), pointer             :: fi !solution vector
        real, dimension(:,:,:), pointer             :: qq !independent term
        real(8), dimension(:,:,:), pointer          :: cp !diagonal coefficents
        real, dimension(:,:,:), pointer             :: cn, cs !coefficients corresponding to I - should be X or Y dir
        real, dimension(:,:,:), pointer             :: ce, cw !coefficients corresponding to J - should be Z
        integer                                     :: ilb, iub, jlb, jub, klb, kub !bounds
        integer                                     :: normal ! direction normal to the plan in which !
                                                              !we are solving 0 = X, 1 = Y, 2 = Z
        integer                                     :: maxit !max number of iterations
        real                                        :: resmax !maximum value allowed for residual
        integer, dimension(:,:,:), pointer          :: mapping !cells mapping
        real                                        :: alphaLU !alpha coeffient of Stone's method
        logical                                     :: normalized !specifies if abs or norm residuals must 
                                                                  !be used in convergence test
        !Local--------------------------------------------------------------
        integer                                     :: IU, JU, KL, KU
        integer                                     :: i, j, k, n
        real                                        :: res0, resn, rsm
        real                                        :: rho0, rho1, alpha, omega, gama, tt, ts, beta
        real, dimension(:,:), pointer               :: r, r0, v, p, s, t, sk, pk
        real, dimension(:,:), pointer               :: lpr, lw, ls, ue, un
        real(8), dimension(:,:), pointer            :: ap
        real, dimension(:,:), pointer               :: aw, ae, an, as, q
        real, dimension(:,:), pointer               :: x
        integer, dimension(:,:), pointer            :: mapp
        logical                                     :: EndIterate
        real                                        :: global_residual 
        integer                                     :: global_n !numb. of it.
        !Begin--------------------------------------------------------------
        !
        !----Calculate index bounds
        if(normal == 0) then
           IU = iub - ilb + 2
           JU = kub - klb + 2
           KL = jlb; KU = jub
        elseif(normal == 1) then
           IU = jub - jlb + 2
           JU = kub - klb + 2
           KL = ilb; KU = iub
        elseif(normal == 2) then
           IU = iub - ilb + 2
           JU = jub - jlb + 2
           KL = klb; KU = kub
        endif
        !----Allocate remaining matrices
        allocate(r(1:IU + 1,1:JU + 1))
        allocate(r0(1:IU + 1,1:JU + 1))
        allocate(t(1:IU + 1,1:JU + 1))
        allocate(s(1:IU + 1,1:JU + 1))
        allocate(sk(1:IU + 1,1:JU + 1))
        allocate(p(1:IU + 1,1:JU + 1))
        allocate(pk(1:IU + 1,1:JU + 1))
        allocate(v(1:IU + 1,1:JU + 1))
        allocate(lpr(1:IU + 1,1:JU + 1))
        allocate(lw(1:IU + 1,1:JU + 1))
        allocate(ls(1:IU + 1,1:JU + 1))
        allocate(un(1:IU + 1,1:JU + 1))
        allocate(ue(1:IU + 1,1:JU + 1))
        global_residual = 0.
        global_n = 0
        !----Solve for each K (in the non-solving direction)
        do k = KL, KU
        !--------Initialize values
            v(:, :) = 0.
            r(:, :) = 0.
            r0(:,:) = 0.
            s(:, :) = 0.
            sk(:, :) = 0.
            t(:, :) = 0.
            p(:, :) = 0.
            pk(:, :) = 0.
            lw(:,:) = 0.
            ls(:,:) = 0.
            un(:,:) = 0.
            ue(:,:) = 0.
            lpr(:,:) = 1.
        !--------Address Pointers
        if(normal == 0) then
            ap => cp(:, k, :)
            an => cn(:, k, :)
            as => cs(:, k, :)
            ae => ce(:, k, :)
            aw => cw(:, k, :)
            q => qq(:, k, :)
            x  => fi(:, k, :)
            mapp  => mapping(:, k, :)
        elseif(normal == 1) then
            ap => cp(k, :, :)
            an => cn(k, :, :)
            as => cs(k, :, :)
            ae => ce(k, :, :)
            aw => cw(k, :, :)
            q => qq(k, :, :)
            x  => fi(k, :, :)
            mapp  => mapping(k, :, :)
        elseif(normal == 2) then
            ap => cp(:, :, k)
            an => cn(:, :, k)
            as => cs(:, :, k)
            ae => ce(:, :, k)
            aw => cw(:, :, k)
            x => fi(:, :, k)
            q => qq(:, :, k)
            mapp  => mapping(:, :, k)
        end if
        !--------Calculate intial residual vector
            res0 = 0.0
            do i = 2, IU
            do j = 2, JU
                r(i, j) =   q(i, j)                              &
                          - ap(i, j) * x(i    , j    )           &
                          - aw(i, j) * x(i    , j - 1)           &
                          - ae(i, j) * x(i    , j + 1)           &
                          - as(i, j) * x(i - 1, j    )           &
                          - an(i, j) * x(i + 1, j    )
                r0(i, j) = r(i, j) 
                res0 = res0 + abs(r(i, j))
            enddo
            enddo
        !--------Calculate LU Coefficients (Pre-Conditioning Step)
            do i = 2, IU
            do j = 2, JU
                lw(i, j) = aw(i, j) / (1 + alphaLU * un(i    , j - 1))
                ls(i, j) = as(i, j) / (1 + alphaLU * ue(i - 1, j    ))
                lpr(i, j) =   ap(i, j)                                &
                            + alphaLU * lw(i, j) * un(i    , j - 1)   &
                            + alphaLU * ls(i, j) * ue(i - 1, j    )   &
                            -           lw(i, j) * ue(i    , j - 1)   &
                            -           ls(i, j) * un(i - 1, j    )
                un(i, j)  =   (an(i, j) - alphaLU * lw(i, j) * un(i - 1, j    ))    &
                            / lpr(i, j)
                ue(i, j)  =   (ae(i, j) - alphaLU * ls(i, j) * ue(i    , j - 1))    &
                            / lpr(i, j)
            enddo
            enddo
        !--------Inner iteration loop
            n = 0
            EndIterate = .false.
            do while((n < maxit).and.(.not.EndIterate))
                n = n + 1
                rho1 = - omega * rho0
                rho0 = 0.
        !--------Dot product of residual vectors <r,r0>
                do i = 2, IU
                do j = 2, JU
                    rho0 = rho0 + r(i, j) * r0(i, j)
                enddo
                enddo
                beta = alpha * rho0 / (rho1  + 1.e-37)
        !--------Calculation of p0
                do i = 2, IU
                do j = 2, JU
                    p(i, j) = r(i, j) - beta * (p(i, j) - omega * v(i, j))
                enddo
                enddo
        !--------Solve [M][p] = [p0] with [M] = [L][U]
        !
        !------------Calculate [INV([L])][p0]
                pk = p
                do i = 2, IU
                do j = 2, JU
                    pk(i, j) = (   pk(i, j)                          &
                                - ls(i, j) * pk(i - 1, j    )        &
                                - lw(i, j) * pk(i    , j - 1)        &
                                ) / lpr(i, j)
                enddo
                enddo
        !------------Backward Substitution [U][p] = [INV([L])][p0]
                do i = IU, 2, -1
                do j = JU, 2, -1
                    pk(i, j) =   pk(i, j)                               &
                               - un(i, j) * pk(i + 1, j    )            &
                               - ue(i, j) * pk(i    , j + 1)
                enddo
                enddo
        !--------Calculate [v] = [A][p]
                do i = 2, IU
                do j = 2, JU
                    v(i, j) =   ap(i, j) * pk(i, j)                       &
                              + as(i, j) * pk(i - 1, j    )               &
                              + an(i, j) * pk(i + 1, j    )               &
                              + aw(i, j) * pk(i    , j - 1)               &
                              + ae(i, j) * pk(i    , j + 1)
                    
                enddo
                enddo
        !--------Dot product <v,r0>
                gama = 0.
                do i = 2, IU
                do j = 2, JU
                    gama = gama + v(i, j) * r0(i, j)
                enddo
                enddo
                alpha = rho0 / (gama + 1.e-20)
        !--------Calculation of s0
                do i = 2, IU
                do j = 2, JU
                    s(i, j) = r(i, j) - alpha * v(i, j)
                enddo
                enddo
        !--------Solve [M][sk] = [s0] with [M] = [L][U]
        !
        !------------Calculate [INV([L])][s0]
                sk = s
                do i = 2, IU
                do j = 2, JU
                    sk(i, j) = (   sk(i, j)                               &
                                - ls(i, j) * sk(i - 1, j    )             &
                                - lw(i, j) * sk(i    , j - 1)             &
                                ) / lpr(i, j)
                enddo
                enddo
        !------------Backward Substitution [U][sk] = [INV([L])][s0]
                do i = IU, 2, -1
                do j = JU, 2, -1
                        sk(i, j) =   sk(i, j)                                     &             
                                  - un(i, j) * sk(i + 1, j    )                   &
                                  - ue(i, j) * sk(i    , j + 1)
                enddo
                enddo
        !--------Calculate [t] = [A][s]
                do i = 2, IU
                do j = 2, JU
                        t(i, j) =   ap(i, j) * sk(i, j)                       &
                                  + as(i, j) * sk(i - 1, j    )               &
                                  + an(i, j) * sk(i + 1, j    )               &
                                  + aw(i, j) * sk(i    , j - 1)               &
                                  + ae(i, j) * sk(i    , j + 1)
                enddo
                enddo
        !--------Calculate omega = <t,s> / <t,t>
                ts = 0.
                tt = 0.
                do i = 2, IU
                do j = 2, JU
                        ts = ts + t(i, j) * s(i, j)
                        tt = tt + t(i, j) * t(i, j)
                enddo
                enddo
                omega = ts / (tt + 1.e-20)
        !--------Re-Evaluate solution and calculate residual
                resn = 0.
                do i = 2, IU
                do j = 2, JU
                    x(i, j) = x(i, j) + alpha * pk(i, j) + omega * sk(i, j)
                    r(i, j) = s(i, j) - omega * t(i, j)
                    resn = resn + abs(r(i, j))
                enddo
                enddo
        !--------Check convergence
                if(normalized) then
                    rsm = resn / (res0 + 1.e-20)
                else
                    rsm = resn
                endif
                if(rsm.LT.resmax) EndIterate = .true.
            enddo !while - inner loop
            global_residual = global_residual + resn
            global_n = global_n + n
            if(n>=maxit) then
                write(*,*) 'ERR - NH MODEL - Solution not converging to specified value.' 
                write(*,*) " - Residual : ", resn
                stop
            endif
        enddo !do KK
        write(*, *)" - NH MODEL - Global Residual : ", global_residual 
        write(*, *)" - NH MODEL - Number of iterations : ", global_n - (KU - KL + 1)
        !----De-Allocate Local Matrices
        deallocate(r)
        nullify(r)
        deallocate(r0)
        nullify(r0)
        deallocate(s)
        nullify(s)
        deallocate(v)
        nullify(v)
        deallocate(sk)
        nullify(sk)
        deallocate(t)
        nullify(t)
        deallocate(p)
        nullify(p)
        deallocate(pk)
        nullify(pk)
        deallocate(ls)
        nullify(ls)
        deallocate(lw)
        nullify(lw)
        deallocate(ue)
        nullify(ue)
        deallocate(un)
        nullify(un)
        deallocate(lpr)
        nullify(lpr)
    End Subroutine BICGSTAB2D

    !--------------------------------------------------------------------------

!*************************************************************************
!/*
! * Peter Daly
! * MIT Ocean Acoustics
! * pmd@mit.edu
! * 25-MAY-1998
! * 
! Revisions:
!   Jan. 25, 1999 DHG  Port to Fortran 90
!   Mar. 23, 1999 DHG  To add Lewis Dozier's fix to "rr1" calculation 
! * 
! Description:
! * 
! * These routines convert UTM to Lat/Longitude and vice-versa,
! * using the WGS-84 (GPS standard) or Clarke 1866 Datums.
! * 
! * The formulae for these routines were originally taken from
! * Chapter 10 of "GPS: Theory and Practice," by B. Hofmann-Wellenhof,
! * H. Lictenegger, and J. Collins. (3rd ed) ISBN: 3-211-82591-6,
! * however, several errors were present in the text which
! * made their formulae incorrect.
! *
! * Instead, the formulae for these routines was taken from
! * "Map Projections: A Working Manual," by John P. Snyder
! * (US Geological Survey Professional Paper 1395)
! *
! * Copyright (C) 1998 Massachusetts Institute of Technology
! *               All Rights Reserved
! *
! * RCS ID: $Id: convert_datum.c,v 1.2 1998/06/04 20:50:47 pmd Exp pmd $
! */
!
!*************************************************************************
!
    subroutine ComputeGridZone (longitude, latitude, grid_zone, lambda0)

        !Arguments------------------------------------------------------
        real (kind=8), intent (IN)              ::  longitude, latitude
        integer      , intent (OUT)             ::  grid_zone(2)
        real (kind=8), intent (OUT), optional   :: lambda0

        !Local----------------------------------------------------------
        integer  zone_long, zone_lat

        real (kind=8) M_PI

        m_pi = ACOS (-1.0)

        !  /* Solve for the grid zone, returns the central meridian */

        zone_long = INT ((longitude + 180.0) / 6.0) + 1
        zone_lat = NINT ((latitude + 80.0) / 8.0)
        grid_zone(1) = zone_long
        grid_zone(2) = zone_lat

        !  /* First, let's take care of the polar regions */

        if (present(lambda0)) then

            if ((latitude < -80.0) .OR. (latitude > 84.0)) then
                lambda0 = 0.0 * M_PI / 180.0
                return
            endif

            !  /* Now the special "X" grid */

            if (latitude .GT. 72.0 .AND. &
                longitude .GT. 0.0 .AND. longitude .LT. 42.0) then
                if (longitude .LT. 9.0) then
                    lambda0 = 4.5 * M_PI / 180.0
                elseif (longitude .LT. 21.0) then
                    lambda0 = 15.0 * M_PI / 180.0
                elseif (longitude .LT. 33.0) then
                    lambda0 = 27.0 * M_PI / 180.0
                elseif (longitude .LT. 42.0) then
                    lambda0 = 37.5 * M_PI / 180.0
                endif
                return
            endif

            !  /* Handle the special "V" grid */

            if (latitude .GT. 56.0 .AND. latitude .LT. 64.0 .AND. &
                longitude .GT. 0.0 .AND. longitude .LT. 12.0) then
                if (longitude .LT. 3.0) then
                    lambda0 = 1.5 * M_PI / 180.0
                elseif (longitude .LT. 12.0) then
                    lambda0 = 7.5 * M_PI / 180.0
                endif
                return
            endif

            !  /* The remainder of the grids follow the standard rule */

            lambda0 = (FLOAT (zone_long - 1) * 6.0 + (-180.0) + 3.0) * M_PI / 180.0

        endif

    end subroutine ComputeGridZone

    !--------------------------------------------------------------------------

    subroutine GetLambda0 (grid_zone, lambda0, ierr)


        integer       grid_zone(2)
        real (kind=8) lambda0
        integer       ierr

        integer zone_long
        integer zone_lat
        real (kind=8) latitude, longitude

        real (kind=8) M_PI
 
        m_pi = ACOS (-1.0)

        !/* Given the grid zone, then set the central meridian, lambda0 */

        !/* Check the grid zone format */

        zone_long = grid_zone(1)
        zone_lat = grid_zone(2)
        if ((zone_long .LT. 1) .OR. (zone_long .GT. 61)) then
            write (*,*) 'Invalid grid zone format: ', zone_long, zone_lat
            ierr = -1
            return 
        endif

        longitude = (FLOAT (zone_long - 1) * 6.0) - 180.0
        latitude = (FLOAT (zone_lat) * 8.0) - 80.0

        !/* Take care of special cases */

        if ((latitude .LT. -80.0) .OR. (latitude .GT. 84.0)) then
            lambda0 = 0.0
            ierr = 0
            return 
        endif

        if (latitude .GT. 56.0 .AND. latitude .LT. 64.0 .AND. &
            longitude .GT. 0.0 .AND. longitude .LT. 12.0) then
        if (longitude .LT. 3.0) then
            lambda0 = 1.5 * M_PI / 180.0
        elseif (longitude .LT. 12) then
            lambda0 = 7.5 * M_PI / 180.0
        endif
            ierr = 0
            return
        endif

        if (latitude .GT. 72.0 .AND. &
            longitude .GT. 0.0 .AND. longitude < 42.0) then
            if (longitude .LT. 9.0) then
                lambda0 = 4.5 * M_PI / 180.0
            elseif (longitude .LT. 21.0) then
                lambda0 = 15.0 * M_PI / 180.0
            elseif (longitude .LT. 33.0) then
                lambda0 = 27.0 * M_PI / 180.0
            elseif (longitude .LT. 42.0) then
                lambda0 = 37.5 * M_PI / 180.0
            endif
            ierr = 0
            return
        endif

        !/* Now handle standard cases */

        lambda0 = (FLOAT (zone_long - 1) * 6.0 + (-180.0) + 3.0) * M_PI / 180.0

        !/* All done */

        ierr = 0
        return

    end subroutine GetLambda0

   !--------------------------------------------------------------------------

    subroutine GetEllipsoid (datum, a,b)

        !Arguments-------------------------------------------------------------
        integer, intent(in)                         :: datum
        real(8), intent(out)                        :: a,b


        if (datum == CLARKE_1866_DATUM) then      ! CLARKE_1866_DATUM:
            a = 6378206.4
            b = 6356583.8
        elseif (datum == GRS_80_DATUM) then       ! GRS_80_DATUM:
            a = 6378137
            b = 6356752.3
        elseif (datum == WGS_84_DATUM) then       ! WGS_84_DATUM:
            a = 6378137.0             !/* semimajor axis of ellipsoid (meters) */
            b = 6356752.31425         !/* semiminor axis of ellipsoid (meters) */
        elseif (datum == ED_50_DATUM) then      ! ED_50_DATUM: ellipsoid International 1926 (Hayford)
            a = 6378388.0           
            b = 6356911.94613  
        elseif (datum == SPHERE_DATUM) then       !normal sphere (r=6370997m)
            a = 6370997.0
            b = 6370997.0
        else
            write (*,*) 'Unknown datum: ', datum
            return
        endif

    end subroutine GetEllipsoid


    !*************************************************************************

    subroutine LatLonToUTM (longitude, latitude, utm_x, utm_y, grid_zone, datum)

        real (kind=8) latitude, longitude
        real (kind=8) utm_x, utm_y
        integer       grid_zone(2)
        integer       datum

        real (kind=8)  a, b, f, e, e2, e4, e6
        real (kind=8)  phi, lambda, lambda0, phi0, k0
        real (kind=8)  t, rho, x, y, mm, mm0
        real (kind=8)  aa, aa2, aa3, aa4, aa5, aa6
        real (kind=8)  ep2, nn, tt, cc

        real (kind=8) M_PI
!!!   parameter (M_PI = 3.141592654)

        !---------------------------------------------------------------------------


        m_pi = ACOS (-1.0)

        !/* Converts lat/long to UTM, using the specified datum */

        call GetEllipsoid (datum, a, b)

        !/* Calculate flatness and eccentricity */

        f = 1 - (b / a)
        e2 = 2 * f - f * f
        e = sqrt (e2)
        e4 = e2 * e2
        e6 = e4 * e2

        !/* Convert latitude/longitude to radians */

        phi = latitude * M_PI / 180.0
        lambda = longitude * M_PI / 180.0

        !/* Figure out the UTM zone, as well as lambda0 */

        call ComputeGridZone (longitude, latitude, grid_zone, lambda0)
        phi0 = 0.0

        !/* See if this will use UTM or UPS */

        if (latitude .GT. 84.0) then

            !/* use Universal Polar Stereographic Projection (north polar aspect) */

            k0 = 0.994

            t = sqrt ( ((1 - sin (phi)) / (1 + sin (phi))) * &
                   (((1 + e * sin (phi)) / (1 - e * sin (phi))) ** e) )
            rho = 2.0 * a * k0 * t / sqrt ( ((1.0 + e) ** (1.0 + e)) * ((1.0 - e) ** (1.0 - e)) )
            !!! Not needed (dhg) m = cos (phi) / sqrt (1.0 - e2 * sin (phi) * sin (phi))

            x = rho * sin (lambda - lambda0)
            y = -rho * cos (lambda - lambda0)
            !!! Not needed (dhg) k = rho * a * m

            !/* Apply false easting/northing */

            x = x + 2000000.0
            y = y + 2000000.0

        elseif (latitude .LT. -80.0) then

            !/* use Universal Polar Stereographic Projection (south polar aspect) */

            phi = -phi
            lambda = -lambda
            lambda0 = -lambda0

            k0 = 0.994

            t = sqrt (((1.0 - sin (phi)) / (1.0 + sin (phi))) * &
                 ( ( (1.0 + e * sin (phi)) / (1.0 - e * sin (phi)) ** e) ) )
            rho = 2.0 * a * k0 * t / sqrt ( ((1+e) ** (1+e)) * ((1-e) ** (1-e)) )
            !!! Not needed (dhg) m = cos (phi) / sqrt (1.0 - e2 * sin (phi) * sin (phi))

            x = rho * sin (lambda - lambda0)
            y = -rho * cos (lambda - lambda0)
            !!! Not needed (dhg) k = rho * a * m

            x = -x
            y = -y

            !/* Apply false easting/northing */

            x = x + 2000000.0
            y = y + 2000000.0

        else

            !/* Use UTM */

            !/* set scale on central median (0.9996 for UTM) */

            k0 = 0.9996

            mm = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi - &
                  (3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0) * sin (2.0*phi) + &
                  (15.0*e4/256.0 + 45.0*e6/1024.0) * sin (4.0*phi) - &
                  (35.0*e6/3072.0) * sin (6.0*phi))

            mm0 = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi0 - &
                   (3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0) * sin (2.0*phi0) + &
                   (15.0*e4/256.0 + 45.0*e6/1024.0) * sin (4.0*phi0) - &
                   (35.0*e6/3072.0) * sin (6.0*phi0))

            aa = (lambda - lambda0) * cos(phi)
            aa2 = aa * aa
            aa3 = aa2 * aa
            aa4 = aa2 * aa2
            aa5 = aa4 * aa
            aa6 = aa3 * aa3

            ep2 = e2 / (1.0 - e2)
            nn = a / sqrt (1.0 - e2 * sin (phi) * sin (phi))
            tt = tan (phi) * tan (phi)
            cc = ep2 * cos (phi) * cos (phi)

            !!! Not needed (dhg) k = k0 * (1 + (1+cc)*aa2/2 + (5-4*tt+42*cc+13*cc*cc-28*ep2) * aa4 / 24.0 + &
            !!! Not needed (dhg)         (61-148*tt+16*tt*tt) * aa6 / 720.0)
            x = k0 * nn * (aa + (1-tt+cc) * aa3 / 6 + &
                      (5-18*tt+tt*tt+72*cc-58*ep2) * aa5 / 120.0)
            y = k0 * (mm - mm0 + nn * tan (phi) * &
                     (aa2 / 2 + (5-tt+9*cc+4*cc*cc) * aa4 / 24.0 + &
                 (61 - 58*tt + tt*tt + 600*cc - 330*ep2) * aa6 / 720))

            !/* Apply false easting and northing */

            x = x + 500000.0
            if (y .LT. 0.0) then
               y = y + 10000000.0
            endif
        endif

        !/* Set entries in UTM structure */

        utm_x = x
        utm_y = y

        !/* done */

        return

    end subroutine LatLonToUTM

!*************************************************************************

    subroutine UTMToLatLon (utm_x, utm_y, longitude, latitude, grid_zone, datum)


        real (kind=8) utm_x, utm_y
        real (kind=8) latitude, longitude
        integer       grid_zone(2)
        integer       datum

        integer ierr
        real (kind=8)  a, b, f, e, e2, e4, e6, e8
        real (kind=8)  lambda0, x, y, k0, rho, t, chi, phi, phi1, phit
        real (kind=8)  lambda, phi0, e1, e12, e13, e14
        real (kind=8)  mm, mm0, mu, ep2, cc1, tt1, nn1, rr1
        real (kind=8)  dd, dd2, dd3, dd4, dd5, dd6

        real (kind=8) M_PI
        real (kind=8) LOWER_EPS_LIMIT
        parameter (LOWER_EPS_LIMIT = 1.0e-14)
        real (kind=8) M_PI_2


    !---------------------------------------------------------------------------


        m_pi = ACOS (-1.0)

        M_PI_2 = M_PI * 2.0

        !/* Converts UTM to lat/long, using the specified datum */

        call GetEllipsoid (datum, a, b)

        !/* Calculate flatness and eccentricity */

        f = 1.0 - (b / a)
        e2 = (2.0 * f) - (f * f)
        e = sqrt (e2)
        e4 = e2 * e2
        e6 = e4 * e2
        e8 = e4 * e4

        !/* Given the UTM grid zone, generate a baseline lambda0 */

        call GetLambda0 (grid_zone, lambda0, ierr)
        if (ierr .NE. 0) then
            write (*,*) 'Unable to translate UTM to LL'
            return
        endif

        latitude = (FLOAT (grid_zone(2)) * 8.0) - 80.0

        !/* Take care of the polar regions first. */

        if (latitude .GT. 84.0) then !/* north polar aspect */

            !/* Subtract the false easting/northing */

            x = utm_x - 2000000.0
            y = utm_y - 2000000.0

            !/* Solve for inverse equations */

            k0 = 0.994
            rho = sqrt (x*x + y*y)
            t = rho * sqrt ( ((1+e) ** (1+e)) * ((1-e) ** (1-e)) ) / (2*a*k0)

            !/* Solve for latitude and longitude */

            chi = M_PI_2 - 2 * atan (t)
            phit = chi + (e2/2 + 5*e4/24 + e6/12 + 13*e8/360) * sin(2*chi) + &
                     (7*e4/48 + 29*e6/240 + 811*e8/11520) * sin(4*chi) + &
                     (7*e6/120 + 81*e8/1120) * sin(6*chi) + &
                     (4279*e8/161280) * sin(8*chi)

            do while (ABS (phi-phit) .GT. LOWER_EPS_LIMIT)
                phi = phit
                phit = M_PI_2 - 2 * atan ( t * (((1 - e * sin (phi)) / (1 + e * sin (phi))) ** (e / 2)) )
            enddo

            lambda = lambda0 + atan2 (x, -y)

        elseif (latitude .LT. -80.0) then !/* south polar aspect */

            !/* Subtract the false easting/northing */

            x = -(utm_x - 2000000)
            y = -(utm_y - 2000000)

            !/* Solve for inverse equations */

            k0 = 0.994
            rho = sqrt (x*x + y*y)
            t = rho * sqrt ( ((1+e) ** (1+e)) * ((1-e) ** (1-e)) ) / (2*a*k0)

            !/* Solve for latitude and longitude */

            chi = M_PI_2 - 2 * atan (t)
            phit = chi + (e2/2 + 5*e4/24 + e6/12 + 13*e8/360) * sin (2*chi) + &
               (7*e4/48 + 29*e6/240 + 811*e8/11520) * sin (4*chi) + &
               (7*e6/120 + 81*e8/1120) * sin (6*chi) + &
               (4279*e8/161280) * sin (8*chi)

            do while (ABS (phi-phit) .GT. LOWER_EPS_LIMIT)
                phi = phit;
                phit = M_PI_2 - 2 * atan (t * ( ((1-e*sin(phi)) / (1+e*sin(phi)) ) ** (e/2)))
            enddo

            phi = -phi
            lambda = -(-lambda0 + atan2 (x , -y))

        else

            !/* Now take care of the UTM locations */

            k0 = 0.9996

            !/* Remove false eastings/northings */

            x = utm_x - 500000.0
            y = utm_y

            if (latitude .LT. 0.0) then  !/* southern hemisphere */
                y = y - 10000000.0
            endif

            !/* Calculate the footpoint latitude */

            phi0 = 0.0
            e1 = (1.0 - sqrt (1.0-e2)) / (1.0 + sqrt (1.0-e2))
            e12 = e1 * e1
            e13 = e1 * e12
            e14 = e12 * e12

            mm0 = a * ((1.0-e2/4.0 - 3.0*e4/64.0 - 5.0*e6/256.0) * phi0 - &
               (3.0*e2/8.0 + 3.0*e4/32.0 + 45.0*e6/1024.0) * sin (2.0*phi0) + &
               (15.0*e4/256.0 + 45.0*e6/1024.0) * sin (4.0*phi0) - &
               (35.0*e6/3072.0) * sin (6.0*phi0))
            mm = mm0 + y/k0;
            mu = mm / (a * (1.0-e2/4.0-3.0*e4/64.0-5.0*e6/256.0))

            phi1 = mu + (3.0*e1/2.0 - 27.0*e13/32.0) * sin (2.0*mu) + &
               (21.0*e12/16.0 - 55.0*e14/32.0) * sin (4.0*mu) + &
               (151.0*e13/96.0) * sin (6.0*mu) + &
               (1097.0*e14/512.0) * sin (8.0*mu)

            !/* Now calculate lambda and phi */

            ep2 = e2 / (1.0 - e2)
            cc1 = ep2 * cos (phi1) * cos (phi1)
            tt1 = tan (phi1) * tan (phi1)
            nn1 = a / sqrt (1.0 - e2 * sin (phi1) * sin (phi1))
            !!!DHG Old Code rr1 = a * (1.0 - e2) / ((1.0 - e2 * sin (phi) * sin (phi)) ** 1.5)
            !!!DHG L.Dozier's fix is next
            rr1 = a * (1.0 - e2) / ((1.0 - e2 * sin (phi1) * sin (phi1)) ** 1.5)
            dd = x / (nn1 * k0)

            dd2 = dd * dd
            dd3 = dd * dd2
            dd4 = dd2 * dd2
            dd5 = dd3 * dd2
            dd6 = dd4 * dd2

            phi = phi1 - (nn1 * tan (phi1) / rr1) * &
              (dd2/2.0 - (5.0+3.0*tt1+10.0*cc1-4.0*cc1*cc1-9.0*ep2) * dd4 / 24.0 + &
              (61.0+90.0*tt1+298.0*cc1+45.0*tt1*tt1-252.0*ep2-3.0*cc1*cc1) * dd6 / 720.0)
            lambda = lambda0 + &
                 (dd - (1.0+2.0*tt1+cc1) * dd3 / 6.0 + &
                 (5.0-2.0*cc1+28.0*tt1-3.0*cc1*cc1+8.0*ep2+24.0*tt1*tt1) * dd5 / 120.0) / cos (phi1)
        endif

        !/* Convert phi/lambda to degrees */

        latitude = phi * 180.0 / M_PI
        longitude = lambda * 180.0 / M_PI

        !/* All done */

        return

    end subroutine UTMToLatLon


    !--------------------------------------------------------------------------

    subroutine LatLonToLambertSP2 (lat,lon, lat_ref,lon_ref, sp1, sp2, datum, x,y)

        !could be replaced by library proj4
        !http://www.posc.org/Epicentre.2_2/DataModel/ExamplesofUsage/eu_cs34e.html

        !Arguments-------------------------------------------------------------
        real(8)                                     :: lat,lon, lat_ref,lon_ref
        real(8)                                     :: sp1, sp2  !standard parallels
        integer                                     :: datum
        real(8), intent(out)                        :: x,y

        !Locals----------------------------------------------------------------
        real(8)                                     :: a, b, flat, e, e2, e05
        real(8)                                     :: torad, es1, es2, esr, es
        real(8)                                     :: m1, m2, dt1,dt2, dtr, dt
        real(8)                                     :: t1, t2, tr, t
        real(8)                                     :: aux1, n, F, r, r_ref, theta


        if (sp1 .GE. sp2) then
            write(*,*)'sp1 must be lower than sp2'
            write(*,*)'sp1 = ', sp1
            write(*,*)'sp2 = ', sp2
            stop 'LambertConformalConicProj_ToMetric - ModuleFunctions - ERR01' 
        endif

        call GetEllipsoid (datum, a,b)

        !/* Calculate flatness and eccentricity */

        flat = 1.0 - (b / a)
        e2 = (2.0 * flat) - (flat * flat)
        e = sqrt (e2)
        e05 = e / 2.
 
        torad = Pi / 180.

        lat     = lat * torad
        lon     = lon * torad
        lat_ref = lat_ref * torad
        lon_ref = lon_ref * torad
        sp1     = sp1 * torad
        sp2     = sp2 * torad

        es1 = e * sin(sp1)
        es2 = e * sin(sp2)
        esr = e * sin(lat_ref)
        es  = e * sin(lat) 
        
        m1 = cos(sp1) / sqrt((1.- es1**2.))
        m2 = cos(sp2) / sqrt((1.- es2**2.))
        
        dt1 = (1. - es1) / (1. + es1)
        dt2 = (1. - es2) / (1. + es2) 
        dtr = (1. - esr) / (1. + esr) 
        dt  = (1. - es)  / (1. + es) 

        t1 = tan(0.25*Pi - 0.5*sp1)     / dt1**(e05)
        t2 = tan(0.25*Pi - 0.5*sp2)     / dt2**(e05)
        tr = tan(0.25*Pi - 0.5*lat_ref) / dtr**(e05)    
        t  = tan(0.25*Pi - 0.5*lat)     / dt **(e05)
        
        aux1    = log(t1) - log(t2)
        
        n       = (log(m1) - log(m2)) / aux1

        F       = m1 / (n * (t1**n))

        r       = a * F * (t**n)
        r_ref   = a * F * (tr**n)

        theta   = n * (lon - lon_ref)

        x = r * sin(theta)
        y = r_ref - r * cos(theta)


    end subroutine LatLonToLambertSP2

    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    
    !Calculates the distance in m between two GPS coordinates using the Haversine formulation (aprox)
    !Equations from http://mathforum.org/library/drmath/view/51879.html
    real function DistanceBetweenTwoGPSPoints(lon1, lat1, lon2, lat2)
        
        !Arguments-------------------------------------------------------------
        real, intent(in)                            :: lon1, lat1, lon2, lat2
        
        !Local-----------------------------------------------------------------
        real                                        :: degreesToRadians = Pi / 180.0
        real                                        :: dlon, dlat, a, c

        dlon = lon2 * degreesToRadians - lon1 * degreesToRadians
        dlat = lat2 * degreesToRadians - lat1 * degreesToRadians
            
        a = Sin(dlat / 2.0)**2.0 + Cos(lat1 * degreesToRadians) * Cos(lat2 * degreesToRadians) * Sin(dlon / 2.0) ** 2.0
        c = 2 * Atan2(Sqrt(a), Sqrt(1 - a))
        
        DistanceBetweenTwoGPSPoints = MeanRadius_ * c;

       

    end function        

#ifdef _USE_PROJ4

    subroutine GeographicToCartesian(lat,lon, params, x, y)

        use proj4

        !Arguments-------------------------------------------------------------
        real(8)                                     :: lat,lon
        character(len=20), dimension(:)             :: params 
        real(8), intent(out)                        :: x,y
        
        !Internal--------------------------------------------------------------
        integer                                     :: status
        type(prj90_projection)                      :: proj

        status=prj90_init(proj,params)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'GeographicToCartesian - ModuleFunctions - ERR01'
        endif

        status = prj90_fwd(proj,lon,lat,x,y)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'GeographicToCartesian - ModuleFunctions - ERR02'
        end if

        status = prj90_free(proj)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'GeographicToCartesian - ModuleFunctions - ERR03'
        end if
    
    end subroutine GeographicToCartesian

    !--------------------------------------------------------------------------

    subroutine CartesianToGeographic (x, y, params, lat,lon)

        use proj4

        !Arguments-------------------------------------------------------------
        real(8)                                     :: x,y
        character(len=20), dimension(:)             :: params 
        real(8), intent(out)                        :: lat,lon
 
        !Internal--------------------------------------------------------------
        integer                                     :: status
        type(prj90_projection)                      :: proj

        status=prj90_init(proj,params)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'CartesianToGeographic - ModuleFunctions - ERR01'
        endif

        status = prj90_inv(proj,x,y,lon,lat)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'CartesianToGeographic - ModuleFunctions - ERR02'
        end if

        status = prj90_free(proj)
        if (status.ne.PRJ90_NOERR) then
            write(*,*) prj90_strerrno(status)
            stop 'CartesianToGeographic - ModuleFunctions - ERR03'
        end if
    
    end subroutine CartesianToGeographic

    !--------------------------------------------------------------------------

#endif

    !--------------------------------------------------------------------------
    
    function Normcrossprod(x, y, z)
        
        !Arguments-------------------------------------------------------------
        real(8), dimension(3)                       :: normcrossprod
        real   , dimension(3), intent(in)           :: x, y, z
        real(8)                                     :: t1(3), t2(3), norm

        t1(1) = x(2) - x(1)
        t1(2) = y(2) - y(1)
        t1(3) = z(2) - z(1)
        t2(1) = x(3) - x(1)
        t2(2) = y(3) - y(1)
        t2(3) = z(3) - z(1)
        normcrossprod(1) = t1(2)*t2(3) - t1(3)*t2(2)
        normcrossprod(2) = t1(3)*t2(1) - t1(1)*t2(3)
        normcrossprod(3) = t1(1)*t2(2) - t1(2)*t2(1)
        norm = sqrt(dot_product(normcrossprod,normcrossprod))
        if (norm /= 0._8) normcrossprod = normcrossprod/norm

    end function normcrossprod

    !--------------------------------------------------------------------------

    !Settling velocity computed as function of a hindered settling concentration
    !Used only for cohesive sediments. Units in this formulation in kg/m3
    real function SettlingVelocity (Concentration, CHS, KL, KL1, M, ML, I,J,K)
        
        !Arguments-------------------------------------------------
        real,       intent(IN) :: Concentration     !kg/m3 or g/l 
        real,       intent(IN) :: CHS               !Hindered settling concentration
        real,       intent(IN) :: KL, KL1, M, ML
        integer,    intent(IN) :: i, j, k

        !Local-----------------------------------------------------
        real                    :: Aux
        
        !Begin-----------------------------------------------------


        Aux = KL1 * (Concentration - CHS)

        if (Concentration < CHS .and. Concentration >= 0.) then 

            SettlingVelocity = KL*(Concentration)**M

        elseif(Aux < 1. .and. Aux >= 0.) then

            SettlingVelocity = KL*(CHS)**M*(1.0-Aux)**ML 

        elseif(Aux > 1. .and. Concentration < 10.) then

            SettlingVelocity = 0. !if concentration is to high settling velocity is null

        else

            write(*,*)'Concentration (g/l)          = ', Concentration
            write(*,*)'KL1 * (Concentration - CHS)  = ', Aux
            write(*,*)'Cell(i, j, k)                = ', i, j, k
            stop 'Error computing the settling velocity - SettlingVelocity - ModuleFreeVerticalMovement'

        endif

    end function SettlingVelocity



    subroutine SLPMIN(Hin,imax,jmax,SLMIN, WaterPoints2D, Slope)

        !Arguments------------------------------------------------
        real,    dimension(:,:), pointer             :: Hin
        integer, dimension(:,:), pointer             :: WaterPoints2D
        integer                         , intent(IN) :: imax, jmax
        real                            , intent(IN) :: SLMIN 
        real,    dimension(:,:), pointer             :: Slope

        !local-----------------------------------------------------
        integer                                      :: i, j
        real                                         :: DH, SN, H1, H2, H3, Hmin


!     This subroutine limits the maximum difference in H divided
!     by twice the average of H of adjacent cells.
!     The maximum possible value is unity.

        Hmin = 0.5
    
!    sweep right
D7:     do I=2,imax-1
D6:     do J=2,jmax-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I,J+1).EQ.1.) then
                H1 = max(Hin(I,J)  ,Hmin)
                H2 = max(Hin(I,J+1),Hmin)
                H3 = H1+H2
                Slope(i,j)=ABS(Hin(I,J+1)-Hin(I,J))/H3
                if(Slope(i,j).GE.SLMIN) then
                  DH=0.5*(Slope(i,j)-SLMIN)*H3
                  SN=-1.
                  if (Hin(I,J+1).GT.Hin(I,J)) SN=1.
                  Hin(I,J+1)=Hin(I,J+1)-SN*DH
                  Hin(I,J)=Hin(I,J)+SN*DH
                endif        
            endif
        enddo D6
!    sweep left
 
D5:     do J=jmax-1,2,-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I,J+1).EQ.1.) then
                H1 = max(Hin(I,J)  ,Hmin)
                H2 = max(Hin(I,J+1),Hmin)
                H3 = H1+H2
                Slope(i,j)=ABS(Hin(I,J+1)-Hin(I,J))/H3
                if (Slope(i,j).GE.SLMIN) then
                  DH=0.5*(Slope(i,j)-SLMIN)*H3
                  SN=-1.
                  IF(Hin(I,J+1).GT.Hin(I,J)) SN=1.
                  Hin(I,J+1)=Hin(I,J+1)-SN*DH
                  Hin(I,J)=Hin(I,J)+SN*DH
                endif
            endif
        enddo D5
        enddo D7

!   sweep up     
D3:     do J=2,jmax-1
D1:     do I=2,imax-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I+1,J).EQ.1.) then
                H1 = max(Hin(I  ,J),Hmin)
                H2 = max(Hin(I+1,J),Hmin)
                H3 = H1+H2
                Slope(i,j)=ABS(Hin(I+1,J)-Hin(I,J))/H3
                if (Slope(i,j).GE.SLMIN) then
                    DH=0.5*(Slope(i,j)-SLMIN)*H3
                    SN=-1.
                    IF(Hin(I+1,J).GT.Hin(I,J)) SN=1.
                    Hin(I+1,J)=Hin(I+1,J)-SN*DH
                    Hin(I,J)=Hin(I,J)+SN*DH
                endif
            endif
        enddo D1

!   sweep down
D2:     do I=imax-1,2,-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I+1,J).EQ.1.) then
                H1 = max(Hin(I  ,J),Hmin)
                H2 = max(Hin(I+1,J),Hmin)
                H3 = H1+H2
                Slope(i,j)=ABS(Hin(I+1,J)-Hin(I,J))/H3
                if (Slope(i,j).GE.SLMIN) then
                    DH=0.5*(Slope(i,j)-SLMIN)*H3
                    SN=-1.
                    IF(Hin(I+1,J).GT.Hin(I,J)) SN=1.
                    Hin(I+1,J)=Hin(I+1,J)-SN*DH
                    Hin(I,J)=Hin(I,J)+SN*DH
                endif
            endif
        enddo D2
        enddo D3

    end subroutine SLPMIN

!--------------------------------------------------------------------------------------


    subroutine SLPMIN2(Hin,DZX, DZY, imax,jmax,SLMIN, WaterPoints2D, Slope)

        !Arguments------------------------------------------------
        real,    dimension(:,:), pointer             :: Hin,DZX, DZY
        integer, dimension(:,:), pointer             :: WaterPoints2D
        integer                         , intent(IN) :: imax, jmax
        real                            , intent(IN) :: SLMIN 
        real,    dimension(:,:), pointer             :: Slope

        !local-----------------------------------------------------
        integer                                      :: i, j
        real                                         :: DH, SN, H1, H2, H3, Hmin


!     This subroutine limits the maximum slope of adjacent cells in H .

        Hmin = 0.5
    
!    sweep right
D7:     do I=2,imax-1
D6:     do J=2,jmax-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I,J+1).EQ.1.) then
                H1 = max(Hin(I,J)  ,Hmin)
                H2 = max(Hin(I,J+1),Hmin)
                H3 = DZX(I, J)
                Slope(i,j)=ABS(Hin(I,J+1)-Hin(I,J))/H3
                if(Slope(i,j).GE.SLMIN) then
                  DH=(Slope(i,j)-SLMIN)*H3
                  SN=-1.
                  if (Hin(I,J+1).GT.Hin(I,J)) SN=1.
                  Hin(I,J+1)=Hin(I,J+1)-SN*DH
                  Hin(I,J)=Hin(I,J)+SN*DH
                endif        
            endif
        enddo D6
!    sweep left
 
D5:     do J=jmax-1,2,-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I,J+1).EQ.1.) then
                H1 = max(Hin(I,J)  ,Hmin)
                H2 = max(Hin(I,J+1),Hmin)
                H3 = DZX(I, J)
                Slope(i,j)=ABS(Hin(I,J+1)-Hin(I,J))/H3
                if (Slope(i,j).GE.SLMIN) then
                  DH=(Slope(i,j)-SLMIN)*H3
                  SN=-1.
                  IF(Hin(I,J+1).GT.Hin(I,J)) SN=1.
                  Hin(I,J+1)=Hin(I,J+1)-SN*DH
                  Hin(I,J)=Hin(I,J)+SN*DH
                endif
            endif
        enddo D5
        enddo D7

!   sweep up     
D3:     do J=2,jmax-1
D1:     do I=2,imax-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I+1,J).EQ.1.) then
                H1 = max(Hin(I  ,J),Hmin)
                H2 = max(Hin(I+1,J),Hmin)
                H3 = DZY(I, J)
                Slope(i,j)=ABS(Hin(I+1,J)-Hin(I,J))/H3
                if (Slope(i,j).GE.SLMIN) then
                    DH=(Slope(i,j)-SLMIN)*H3
                    SN=-1.
                    IF(Hin(I+1,J).GT.Hin(I,J)) SN=1.
                    Hin(I+1,J)=Hin(I+1,J)-SN*DH
                    Hin(I,J)=Hin(I,J)+SN*DH
                endif
            endif
        enddo D1

!   sweep down
D2:     do I=imax-1,2,-1
            if (WaterPoints2D(I,J).EQ.1..and.WaterPoints2D(I+1,J).EQ.1.) then
                H1 = max(Hin(I  ,J),Hmin)
                H2 = max(Hin(I+1,J),Hmin)
                H3 = DZY(I, J)
                Slope(i,j)=ABS(Hin(I+1,J)-Hin(I,J))/H3
                if (Slope(i,j).GE.SLMIN) then
                    DH=(Slope(i,j)-SLMIN)*H3
                    SN=-1.
                    IF(Hin(I+1,J).GT.Hin(I,J)) SN=1.
                    Hin(I+1,J)=Hin(I+1,J)-SN*DH
                    Hin(I,J)=Hin(I,J)+SN*DH
                endif
            endif
        enddo D2
        enddo D3

    end subroutine SLPMIN2

!--------------------------------------------------------------------------------------


    real(8) function CalcPotentialEnergy(Size, Density, Volume, Height)

        !Arguments-----------------------------------------------
        real, dimension(:), pointer                                 :: Density, Volume, Height
        integer                                                     :: Size
        
        !locals--------------------------------------------------
        integer                                                     :: m

        CalcPotentialEnergy = 0.
        do m = 1, Size
            CalcPotentialEnergy = CalcPotentialEnergy + Density(m)*Volume(m)*Height(m)
        enddo
        CalcPotentialEnergy = Gravity * CalcPotentialEnergy

    end function CalcPotentialEnergy

    integer function GetNumberWaterpoints3D(Waterpoints3D, Size3D)

        !Arguments---------------------------------------------
        type(T_Size3D), pointer                                     :: Size3D
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D

        !Locals------------------------------------------------
        integer                                                     :: i,j,k,m

        m = 0;
        do i = Size3D%ILB, Size3D%IUB
        do j = Size3D%JLB, Size3D%JUB
        do k = Size3D%KLB, Size3D%KUB

            if(Waterpoints3D(i,j,k) == 1) then
                m = m + 1
            end if
            
        end do 
        end do
        end do

        GetNumberWaterpoints3D = m

    end function

    subroutine CalculateStageAreas( VolumeZ, VolumeHeight, Areas, Size3D, Waterpoints3D)
    
        !Arguments-------------------------------------------------------------------------
        type(T_Size3D), pointer                                     :: Size3D
        real, dimension(:,:,:), pointer                             :: VolumeZ, VolumeHeight    
        real, dimension(:    ), pointer                             :: Areas    
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D        

        !Locals----------------------------------------------------------------------------
        integer                                                     :: i,j,k,m,p,pstart,pstop
        real                                                        :: Area

        m = 0;
        pstart = 1;
        pstop = 1;
        do k = Size3D%KLB, Size3D%KUB

            Area = 0.;

            do j = Size3D%JLB, Size3D%JUB
            do i = Size3D%ILB, Size3D%IUB

                if (Waterpoints3D(i,j,k) == 1) then

                    m = m + 1
                    Area = Area + VolumeZ(i,j,k)/VolumeHeight(i,j,k)

                end if

            enddo
            enddo

            pstop = m
            do p = pstart, pstop
                Areas(p) = Area
            enddo
            pstart = pstop

        enddo

    end subroutine CalculateStageAreas

    subroutine SortNumerically_3D( TargetReference, TargetAcolyte, SourceReference, SourceAcolyte, Size3D, Waterpoints3D)
        
        !Arguments-------------------------------------------------------------------------
        type(T_Size3D), pointer                                     :: Size3D
        real, dimension(:,:,:), pointer                             :: SourceReference, SourceAcolyte
        real, dimension(:    ), pointer                             :: TargetReference, TargetAcolyte
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D

        !Locals----------------------------------------------------------------------------
        integer                                                     :: i,j,k,m,p

        m = 0;
        do i = Size3D%ILB, Size3D%IUB
        do j = Size3D%JLB, Size3D%JUB
        do k = Size3D%KLB, Size3D%KUB

            if(Waterpoints3D(i,j,k) == 1) then
                m = m + 1
                TargetReference(m) = SourceReference(i,j,k)
                TargetAcolyte(m) = SourceAcolyte(i,j,k)
            end if
            
            p = m
            do while (p>1 .AND. TargetReference(p) > TargetReference(p-1)) 
                call swap(TargetReference, p, p-1)
                call swap(TargetAcolyte, p, p-1)
                p=p-1
            end do

        end do 
        end do
        end do

    end subroutine SortNumerically_3D

    subroutine DetermineCentreHeight_3D( TargetHeight, SourceHeight, Size3D, Waterpoints3D, method, VolumeZ, Areas)
        
        !Arguments-------------------------------------------------------------------------
        type(T_Size3D), pointer                                     :: Size3D
        real, dimension(:,:,:), pointer                             :: SourceHeight
        real, dimension(:    ), pointer                             :: TargetHeight
        integer, dimension(:,:,:), pointer                          :: Waterpoints3D
        integer                                                     :: method
        real, dimension(:,:,:), pointer, optional                   :: VolumeZ
        real, dimension(:), pointer, optional                       :: Areas

        !Locals----------------------------------------------------------------------------
        integer                                                     :: m
        integer                                                     :: i,j,k
        real                                                        :: HeightEdge, HeightEdge_old
        real                                                        :: ElementHeight


        select case (method)

            case(SimpleHeight_) 
                m = 0;
                do i = Size3D%ILB, Size3D%IUB
                do j = Size3D%JLB, Size3D%JUB
                do k = Size3D%KLB, Size3D%KUB

                    if(Waterpoints3D(i,j,k) == 1) then
                        m = m + 1
                        TargetHeight(m) = SourceHeight(i,j,k)
                    end if
            
                end do 
                end do
                end do

            case(ComplexHeight_)
                HeightEdge_old  = 0.; 
                HeightEdge      = 0.; 
                m = 0;
                do i = Size3D%ILB, Size3D%IUB
                do j = Size3D%JLB, Size3D%JUB
                do k = Size3D%KLB, Size3D%KUB

                    if(Waterpoints3D(i,j,k) == 1) then
                        m = m + 1
                        ElementHeight = VolumeZ(i,j,k) / Areas(m)
                        if (m .GT. 1) then
                            HeightEdge = HeightEdge_old + ElementHeight
                        end if
                        TargetHeight(m) = HeightEdge + 0.5 * ElementHeight
                        HeightEdge_old = HeightEdge
                    end if
            
                end do 
                end do
                end do

        end select

    end subroutine DetermineCentreHeight_3D
    
    subroutine swap(Vec, m, n)

        !Arguments-------------------------------------------------------
        real, dimension(:)                              :: Vec
        integer                                         :: m, n
        
        !Locals----------------------------------------------------------
        real                                            :: a

        a = Vec(m)
        Vec(m) = Vec(n)
        Vec(n) = a

    end subroutine
    
    !--------------------------------------------------------------------------------------
    !
    !Centralizes way to calculate chunk size. For future optimization
    !
    integer function Chunk_K(KLB, KUB)
    
        !Arguments-------------------------------------------------------
        integer, intent(IN)                         :: KLB, KUB
        
        !Locals----------------------------------------------------------
        
        Chunk_K = max((KUB - KLB) / 3, 1)

    end function Chunk_K
        
    !--------------------------------------------------------------------------------------
    !
    !Centralizes way to calculate chunk size. For future optimization
    !
    integer function Chunk_J(JLB, JUB)
    
        !Arguments-------------------------------------------------------
        integer, intent(IN)                         :: JLB, JUB
        
        !Locals----------------------------------------------------------
        
        Chunk_J = max((JUB - JLB) / 10, 1)

    end function Chunk_J

    !--------------------------------------------------------------------------------------
    !
    !Centralizes way to calculate chunk size. For future optimization
    !
    integer function Chunk_I(ILB, IUB)
    
        !Arguments-------------------------------------------------------
        integer, intent(IN)                         :: ILB, IUB
        
        !Locals----------------------------------------------------------
        
        Chunk_I = max((IUB - ILB) / 10, 1)

    end function Chunk_I

    !--------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------

    character(len=15) function TimeToString(Date)

        !Arguments-------------------------------------------------------------
        type(T_Time)                            :: Date
        real,    dimension(6)                   :: AuxTime
        character(len=4)                        :: CharYear
        character(len=2)                        :: CharMonth
        character(len=2)                        :: CharDay
        character(len=2)                        :: CharHour
        character(len=2)                        :: CharMinute
        character(len=2)                        :: CharSecond

        !Begin-----------------------------------------------------------------

        call ExtractDate(Date, Year     = AuxTime(1), Month  = AuxTime(2), &
                               Day      = AuxTime(3), Hour   = AuxTime(4), &
                               Minute   = AuxTime(5), Second = AuxTime(6))
        
        write(CharYear,  '(i4)')int(AuxTime(1))
        write(CharMonth, '(i2)')int(AuxTime(2))
        write(CharDay,   '(i2)')int(AuxTime(3))
        write(CharHour,  '(i2)')int(AuxTime(4))
        write(CharMinute,'(i2)')int(AuxTime(5))
        write(CharSecond,'(i2)')int(AuxTime(6))

        if(len_trim(trim(adjustl(CharMonth)))   < 2)then 
            CharMonth = "0"//trim(adjustl(CharMonth))
        endif
        
        if(len_trim(trim(adjustl(CharDay)))     < 2)then 
            CharDay = "0"//trim(adjustl(CharDay))
        endif

        if(len_trim(trim(adjustl(CharHour)))    < 2)then 
            CharHour = "0"//trim(adjustl(CharHour))
        endif

        if(len_trim(trim(adjustl(CharMinute)))  < 2)then 
            CharMinute = "0"//trim(adjustl(CharMinute))
        endif

        if(len_trim(trim(adjustl(CharSecond)))  < 2)then 
            CharSecond = "0"//trim(adjustl(CharSecond))
        endif

        !Output Format: YYYYMMDD-hhmmss
        TimeToString = CharYear//CharMonth//CharDay//"-"//&
                       CharHour//CharMinute//CharSecond

    end function

    !--------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine WGS84toGoogleMaps2D(lon, lat, ILB, IUB, JLB, JUB, x, y) 

        !Arguments-------------------------------------------------------------
         real,    dimension(:,:), pointer :: lon, lat
         integer                          :: ILB, IUB, JLB, JUB
         real(8), dimension(:,:), pointer :: x, y
        !Local-----------------------------------------------------------------
        integer                           :: i, j
        !Begin-----------------------------------------------------------------
        
        do j=JLB, JUB+1
        do i=ILB, IUB+1        
        
            x(i,j) = lon(i,j) * 20037508.34 / 180;
            if (abs(lat(i,j))<90) then
                y(i,j) = log(tan((90 + lat(i,j)) * Pi / 360)) / (Pi / 180)
                y(i,j) = y(i,j) * 20037508.34 / 180
            else
                write(*,*) 'Out of range - Lat >= 90 or Lat <=-90'
                write(*,*) 'i=',i,' j=',j,' Lat=',Lat(i,j)
                stop 'WGS84toGoogleMaps2D - ModuleHorizontalGrid - ERR10'
            endif
            
        enddo
        enddo
        
    end subroutine WGS84toGoogleMaps2D
  
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    subroutine WGS84toGoogleMaps1D(lon, lat, Dim,  x, y) 

        !Arguments-------------------------------------------------------------
         real(8), dimension(:  ), pointer :: lon, lat
         integer                          :: Dim
         real(8), dimension(:  ), pointer :: x, y
        !Local-----------------------------------------------------------------
        integer                           :: i
        !Begin-----------------------------------------------------------------
        
        do i=1, Dim
            x(i) = lon(i) * 20037508.34 / 180;
        enddo
        
        do i=1, Dim
            if (abs(lat(i))<90) then
                y(i) = log(tan((90 + lat(i)) * Pi / 360)) / (Pi / 180)
                y(i) = y(i) * 20037508.34 / 180
            else
                write(*,*) 'Out of range - Lat >= 90 or Lat <=-90'
                stop 'WGS84toGoogleMap1D - ModuleHorizontalGrid - ERR10'
            endif
        enddo
        
    end subroutine WGS84toGoogleMaps1D
  
!------------------------------------------------------------------------------  

    !--------------------------------------------------------------------------

    character(Len = Pathlength) function ChangeSuffix(Filename, NewSuffix)

        !Arguments-------------------------------------------------------------
        character(Len = *)                  :: Filename
        character(Len = *)                  :: NewSuffix

        !Local-----------------------------------------------------------------
        Integer                             :: Position
        Logical                             :: Back
        
        !Initializing our variable
        ChangeSuffix = ""       

        !Search the first "." char in the string filename (ex:"d:\aplica\filename.dat")
        !counting backwards.
        Back = .true.
        Position = scan(Filename, ".", Back)

        !Make sure variable Position is below StringLength
        if (Position > Pathlength)                            &
            stop 'ChangeSuffix - Module Functions - ERR01'

        ChangeSuffix = Filename(1:Position-1)//trim(NewSuffix)
        
    end function

    !--------------------------------------------------------------------------
    !griflet: adding minval and maxval ersatz functions in order to lift
    !stacksize and heapsize limitations
    
    function minival1D_R4(array, size1D)
            
        real(4), dimension(:), pointer     :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        real(4)                           :: minival1D_R4
        
        minival1D_R4 = 1e15
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueReal .and. minival1D_R4 > array(i) ) then
                minival1D_R4 = array(i)
            end if
        end do
        
    end function minival1D_R4

    function minival2D_R4(array, size2D)
            
        real(4), dimension(:,:), pointer    :: array
        type(T_size2D)                      :: size2D
        integer                             :: i,j
        real(4)                             :: minival2D_R4
        
        minival2D_R4 = 1e15
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if (array(i,j) .ne. FillValueReal .and. minival2D_R4 > array(i,j)) then
                minival2D_R4 = array(i,j)
            end if
        end do
        end do
        
    end function minival2D_R4

    function minival3D_R4(array, size3D)
            
        real(4), dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        real(4)                            :: minival3D_R4
        
        minival3D_R4 = 1e15
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueReal .and. minival3D_R4 > array(i,j,k) ) then
                minival3D_R4 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function minival3D_R4

    function maxival1D_R4(array, size1D)
            
        real(4), dimension(:), pointer     :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        real(4)                            :: maxival1D_R4
        
        maxival1D_R4 = -1e15
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueReal .and. maxival1D_R4 < array(i) ) then
                maxival1D_R4 = array(i)
            end if
        end do
        
    end function maxival1D_R4

    function maxival2D_R4(array, size2D)
            
        real(4), dimension(:,:), pointer   :: array
        type(T_size2D)                  :: size2D
        integer                         :: i,j
        real(4)                            :: maxival2D_R4
        
        maxival2D_R4 = -1e15
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueReal .and. maxival2D_R4 < array(i,j) ) then
                maxival2D_R4 = array(i,j)
            end if
        end do
        end do
        
    end function maxival2D_R4

    function maxival3D_R4(array, size3D)
            
        real(4), dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        real(4)                            :: maxival3D_R4
        
        maxival3D_R4 = -1e15
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueReal .and. maxival3D_R4 < array(i,j,k) ) then
                maxival3D_R4 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function maxival3D_R4

    function minival1D_R8(array, size1D)
            
        real(8), dimension(:), pointer     :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        real(8)                            :: minival1D_R8
        
        minival1D_R8 = 1e15
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueReal .and. minival1D_R8 > array(i) ) then
                minival1D_R8 = array(i)
            end if
        end do
        
    end function minival1D_R8

    function minival2D_R8(array, size2D)
            
        real(8), dimension(:,:), pointer   :: array
        type(T_size2D)                  :: size2D
        integer                         :: i,j
        real(8)                         :: minival2D_R8
        
        minival2D_R8 = 1e15
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueReal .and. minival2D_R8 > array(i,j) ) then
                minival2D_R8 = array(i,j)
            end if
        end do
        end do
        
    end function minival2D_R8

    function minival3D_R8(array, size3D)
            
        real(8), dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        real(8)                         :: minival3D_R8
        
        minival3D_R8 = 1e15
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueReal .and. minival3D_R8 > array(i,j,k) ) then
                minival3D_R8 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function minival3D_R8

    function maxival1D_R8(array, size1D)
            
        real(8), dimension(:), pointer  :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        real(8)                         :: maxival1D_R8
        
        maxival1D_R8 = -1e15
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueReal .and. maxival1D_R8 < array(i) ) then
                maxival1D_R8 = array(i)
            end if
        end do
        
    end function maxival1D_R8

    function maxival2D_R8(array, size2D)
            
        real(8), dimension(:,:), pointer   :: array
        type(T_size2D)                  :: size2D
        integer                         :: i,j
        real(8)                         :: maxival2D_R8
        
        maxival2D_R8 = -1e15
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueReal .and. maxival2D_R8 < array(i,j) ) then
                maxival2D_R8 = array(i,j)
            end if
        end do
        end do
        
    end function maxival2D_R8

    function maxival3D_R8(array, size3D)
            
        real(8), dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        real(8)                         :: maxival3D_R8
        
        maxival3D_R8 = -1e15
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueReal .and. maxival3D_R8 < array(i,j,k) ) then
                maxival3D_R8 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function maxival3D_R8

    function minival1D_I4(array, size1D)
            
        integer, dimension(:), pointer     :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        integer                            :: minival1D_I4
        
        minival1D_I4 = 1e9
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueInt .and. minival1D_I4 > array(i) ) then
                minival1D_I4 = array(i)
            end if
        end do
        
    end function minival1D_I4

    function minival2D_I4(array, size2D)
            
        integer, dimension(:,:), pointer   :: array
        type(T_size2D)                     :: size2D
        integer                            :: i,j
        integer                            :: minival2D_I4
        
        minival2D_I4 = 1e9
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueInt .and. minival2D_I4 > array(i,j) ) then
                minival2D_I4 = array(i,j)
            end if
        end do
        end do
        
    end function minival2D_I4

    function minival3D_I4(array, size3D)

        integer, dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        integer                            :: minival3D_I4

        minival3D_I4 = 1e9
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueInt .and. minival3D_I4 > array(i,j,k) ) then
                minival3D_I4 = array(i,j,k)
            end if
        end do
        end do
        end do

    end function minival3D_I4

    function maxival1D_I4(array, size1D)

        integer, dimension(:), pointer  :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        integer                         :: maxival1D_I4
        
        maxival1D_I4 = -1e9
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueInt .and. maxival1D_I4 < array(i) ) then
                maxival1D_I4 = array(i)
            end if
        end do
        
    end function maxival1D_I4

    function maxival2D_I4(array, size2D)
            
        integer, dimension(:,:), pointer   :: array
        type(T_size2D)                  :: size2D
        integer                         :: i,j
        integer                            :: maxival2D_I4
        
        maxival2D_I4 = -1e9
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueInt .and. maxival2D_I4 < array(i,j) ) then
                maxival2D_I4 = array(i,j)
            end if
        end do
        end do
        
    end function maxival2D_I4

    function maxival3D_I4(array, size3D)
            
        integer, dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        integer                         :: maxival3D_I4
        
        maxival3D_I4 = -1e9
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueInt .and. maxival3D_I4 < array(i,j,k) ) then
                maxival3D_I4 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function maxival3D_I4

    function minival1D_I8(array, size1D)
            
        integer(8), dimension(:), pointer     :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        integer(8)                            :: minival1D_I8
        
        minival1D_I8 = 1e15
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueInt .and. minival1D_I8 > array(i) ) then
                minival1D_I8 = array(i)
            end if
        end do
        
    end function minival1D_I8

    function minival2D_I8(array, size2D)
            
        integer(8), dimension(:,:), pointer   :: array
        type(T_size2D)                  :: size2D
        integer                         :: i,j
        integer(8)                         :: minival2D_I8
        
        minival2D_I8 = 1e15
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueInt .and. minival2D_I8 > array(i,j) ) then
                minival2D_I8 = array(i,j)
            end if
        end do
        end do
        
    end function minival2D_I8

    function minival3D_I8(array, size3D)
            
        integer(8), dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        integer(8)                         :: minival3D_I8
        
        minival3D_I8 = 1e15
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueInt .and. minival3D_I8 > array(i,j,k) ) then
                minival3D_I8 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function minival3D_I8

    function maxival1D_I8(array, size1D)
            
        integer(8), dimension(:), pointer  :: array
        type(T_size1D)                  :: size1D
        integer                         :: i
        integer(8)                         :: maxival1D_I8
        
        maxival1D_I8 = -1e15
        do i = size1D%ILB,size1D%IUB
            if ( array(i) .ne. FillValueInt .and. maxival1D_I8 < array(i) ) then
                maxival1D_I8 = array(i)
            end if
        end do
        
    end function maxival1D_I8

    function maxival2D_I8(array, size2D)
            
        integer(8), dimension(:,:), pointer   :: array
        type(T_size2D)                  :: size2D
        integer                         :: i,j
        integer(8)                         :: maxival2D_I8
        
        maxival2D_I8 = -1e15
        do j = size2D%JLB,size2D%JUB
        do i = size2D%ILB,size2D%IUB
            if ( array(i,j) .ne. FillValueInt .and. maxival2D_I8 < array(i,j) ) then
                maxival2D_I8 = array(i,j)
            end if
        end do
        end do
        
    end function maxival2D_I8

    function maxival3D_I8(array, size3D)
            
        integer(8), dimension(:,:,:), pointer :: array
        type(T_size3D)                  :: size3D
        integer                         :: i,j,k
        integer(8)                         :: maxival3D_I8
        
        maxival3D_I8 = -1e15
        do k = size3D%KLB,size3D%KUB
        do j = size3D%JLB,size3D%JUB
        do i = size3D%ILB,size3D%IUB
            if ( array(i,j,k) .ne. FillValueInt .and. maxival3D_I8 < array(i,j,k) ) then
                maxival3D_I8 = array(i,j,k)
            end if
        end do
        end do
        end do
        
    end function maxival3D_I8

end module ModuleFunctions

!----------------------------------------------------------------------------------------------------------
!MOHID Water Modelling System.
!Copyright (C) 1985, 1998, 2002, 2005. Instituto Superior Técnico, Technical University of Lisbon. 
!----------------------------------------------------------------------------------------------------------


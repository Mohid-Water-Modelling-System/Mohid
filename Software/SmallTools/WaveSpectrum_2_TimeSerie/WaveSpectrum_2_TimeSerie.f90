!  WaveSpectrum_2_TimeSerie.f90 
!
!  FUNCTIONS:
!  WaveSpectrum_2_TimeSerie - Convert a wind wave spectrum in sea level and horizontal velocity average in depth
!

!****************************************************************************
!
!  PROGRAM: WaveSpectrum_2_TimeSerie
!
!  PURPOSE:  Convert a wind wave spectrum in sea level and horizontal velocity average in depth
!
!****************************************************************************

program WaveSpectrum_2_TimeSerie

    Use ModuleFunctions

    implicit none
    
    ! Variables
    real                             :: Hs,Tp,gamma,hmax, hmin, dh,fmin,fmax,tmax,dt
    real, dimension(:), allocatable  :: ETA, UU, f, s    
    integer                          :: i, imax, nf
    logical                          :: EqualEnergy

    ! Initialize variables
    Hs = 2.0
    Tp = 18.0
    gamma = 3.3
    hmax = 30.0
    hmin = 0
    dh = 0.5
    fmin = 0.03
    fmax = 0.3
    dt = 0.1  ! Initialize dt before using it to compute tmax
    tmax = 900
    nf   = 90

    ! Allocate arrays based on imax
    !allocate(ETA(1:imax))
    !allocate(UU(1:imax))

    ! Body of Console1
    print *, 'From spectrum to time serie'
    print *, 'SeaLevel_Vel.txt - Time vs Hydrodynamics'
    print *, 'Frequency_SpectralDensity.txt - Frequency vs Spectral Density'

    call JONSWAP2(Hs, Tp, gamma, hmax, hmin, dh, fmin, fmax, nf, tmax, dt, ETA, UU, f, s,EqualEnergy = .false.)

    open(unit=12, file='SeaLevel_Vel.txt', status='unknown')
    
    imax = int(tmax/dt)+1
    
    write(12,*) 't,n[m],vel[m/s]'
    do i = 1, imax
        write(12,'(f12.4,A1,f8.2,A1,f8.2)') real(i)*dt, ',', ETA(i), ',', UU(i)
    end do
    close(12)
    
    open(unit=12, file='Frequency_SpectralDensity.txt', status='unknown')    
    
    write(12,*) 'f[s-1],s[m2*s]'
    do i = 1, nf
        write(12,'(f12.4,A1,f12.4)') f(i), ',', s(i)
    end do
    close(12)

    deallocate(ETA, UU,f,s)

    end program WaveSpectrum_2_TimeSerie
    
    

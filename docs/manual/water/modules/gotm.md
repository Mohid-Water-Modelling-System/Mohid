# Module GOTM

## General

Module GOTM is a "wrap-up" module containing [GOTM](http://www.gotm.net) routines which were coupled into MOHID in 2001. This routines consist of a set of turbulence models including a k-ε model and Mellor-Yamada second order turbulent closure model (Mellor and Yamada, 1982). GOTM stands for General Ocean Turbulence Model and consists of a one dimensional water column for most important hydrodynamic and thermodynamic processes related to vertical mixing in natural waters. It is subsequently coupled to a three dimensional model for the computation of vertical mixing.

### 1.Introduction

GOTM consists mainly of one dimensional transport equations for momentum, salt and heat. The idea of turbulent closure modules is to have a simple tool to parametrize physical processes on smaller and larger scales as to isolate and study the physical phenomenon in which we are really interested. This procedure is called turbulent closure and GOTM incorporates a vast number of well tested turbulence models which can be applied, e.g.:

- Empirical models

- Energy models

- Two equation models

- Algebraic Stress models

- K - profile parametrisations

In order to apply one of these models, GOTM uses a modular approach. It consists of:

1\. GOTM - The Main Program

2\. The Modular Interface GOTM

3\. The Mean Flow model

4\. The Module Interface Mean Flow

5\. The Turbulence Model

6\. The Module Air - Sea Interaction

#### 1\. The Main Program

The main program merely calls other subroutines in the right order.

#### 2\. The Modular Interface GOTM

This is the general framework which reads the turbulence namelists to initialize the model and manages the global time - stepping and the clean - up after the run is over. As GOTM is programmed in FORTRAN, it will require the definition of a range of units for different purposes;

• unit=10 is reserved for reading namelists

• units 20-29 are reserved for the airsea module

• units 30-39 are reserved for the meanflow module

• units 40-49 are reserved for the turbulence module

• units 50-59 are reserved for the output module

• units 60-69 are reserved for the extra modules like those dealing with sediments or sea-grass

• units 70- are not reserved and can be used as you wish

##### Initialise the model

The model is initialised by the subroutine init\_gotm( )which reads the turbulence namelists and specialised subroutines.

##### Manage globale time - stepping

The routine subroutine time\_loop( ) contains the heart of the code, It contains the main time-loop inside of which all routines required during the time step are called.

Because of this routine following steps are triggered:

- The model is updated and the output is prepared

- Air - sea interactions (flux, SST) are computed

- Time step is performed on mean flow equations (momentum, temperature)

- Some quantities related to shear and stratification are updated(e.g. shear number, buoyancy frequency)

- Turbulence is updated depending on what turbulence model has been specified by the user

- Results are written in the output files

After this the the subroutine clean\_up( )will be called. This function is just a wrapper for the external routine close\_output(). All open files will be closed after this call.

##### Mean Flow Model

This module contains the definitions of the most important mean flow variables used in geophysical models. In GOTM, these are

• the mean horizontal velocity components,U and V

• the mean potential temperature,θ, (or the mean buoyancy,B)

• the mean salinity,S

Each variable Φ can be decomposed into a mean value and a perturbation:

Φ = <Φ> + Φ'

Due to the one-dimensional character of GOTM, the state-variables listed above are assumed to be horizontally homogeneous, depending only on the vertical z-coordinate. As a consequence, all horizontal gradients have to be taken from observations, or they have to be estimated, parameterised or neglected.

## References

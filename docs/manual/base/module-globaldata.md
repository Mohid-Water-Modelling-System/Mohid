# Module GlobalData

## Overview

Module GlobalData, is responsible for global variables such as [properties names](../reference/properties-names.md) and ID numbers, error types ID, constants and parameters, classes registration numbers and some derived types used frequently in other classes (e.g. type T\_Size, a derived type containing the matrixes bounds).

This is mainly static information needed and used equally by all classes. Some methods are also provided by this class, mainly related to checking properties names spelling and attributing ID numbers. It also handles error and used keywords logging and I/O units.

## Sub-routines

### Internal names

- CheckPropertyName
- GetPropertyName
- GetPropertyIDNumber
- Check\_Particulate\_Property
- TranslateTypeZUV

### I/O management

- UnitsManager

### Errors and logs management

- SetError
- WriteErrorMessage
- LogKeyWord

### Start/End management

- StartupMohid
- ShutdownMohid

### Object oriented management routines

- RegisterNewInstance
- AssociateInstance
- DeassociateInstance
- RegisterModule
- ModuleIsRegistered
- Read\_Lock
- Read\_Unlock
- VerifyReadLock

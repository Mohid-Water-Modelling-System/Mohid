#ifndef IPHREEQC_NO_FORTRAN_MODULE
MODULE IPhreeqc
implicit none

! GetSelectedOutputValue TYPES
INTEGER(KIND=4),PARAMETER :: TT_EMPTY  = 0
INTEGER(KIND=4),PARAMETER :: TT_ERROR  = 1
INTEGER(KIND=4),PARAMETER :: TT_DOUBLE = 3
INTEGER(KIND=4),PARAMETER :: TT_STRING = 4

! ERROR RETURN TYPES
INTEGER(KIND=4),PARAMETER :: IPQ_OK           =  0
INTEGER(KIND=4),PARAMETER :: IPQ_OUTOFMEMORY  = -1
INTEGER(KIND=4),PARAMETER :: IPQ_BADVARTYPE   = -2
INTEGER(KIND=4),PARAMETER :: IPQ_INVALIDARG   = -3
INTEGER(KIND=4),PARAMETER :: IPQ_INVALIDROW   = -4
INTEGER(KIND=4),PARAMETER :: IPQ_INVALIDCOL   = -5
INTEGER(KIND=4),PARAMETER :: IPQ_BADINSTANCE  = -6

!!!SAVE
CONTAINS

INTEGER FUNCTION AccumulateLine(id, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION AccumulateLineF(id, line) &
            BIND(C, NAME='AccumulateLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: line(*)
        END FUNCTION AccumulateLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: line
    AccumulateLine = AccumulateLineF(id, trim(line)//C_NULL_CHAR)
    return
END FUNCTION AccumulateLine

INTEGER FUNCTION AddError(id, error_msg)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION AddErrorF(id, error_msg) &
            BIND(C, NAME='AddErrorF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: error_msg(*)
        END FUNCTION AddErrorF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: error_msg
    AddError = AddErrorF(id, trim(error_msg)//C_NULL_CHAR)
    return
END FUNCTION AddError

INTEGER FUNCTION AddWarning(id, warn_msg)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION AddWarningF(id, warn_msg) &
            BIND(C, NAME='AddWarningF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: warn_msg(*)
        END FUNCTION AddWarningF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: warn_msg
    AddWarning = AddWarningF(id, trim(warn_msg)//C_NULL_CHAR)
    return
END FUNCTION AddWarning

INTEGER FUNCTION ClearAccumulatedLines(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION ClearAccumulatedLinesF(id) &
            BIND(C, NAME='ClearAccumulatedLinesF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION ClearAccumulatedLinesF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    ClearAccumulatedLines = ClearAccumulatedLinesF(id)
    return
END FUNCTION ClearAccumulatedLines

INTEGER FUNCTION CreateIPhreeqc()
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION CreateIPhreeqcF() &
            BIND(C, NAME='CreateIPhreeqcF')
            USE ISO_C_BINDING
            IMPLICIT NONE
        END FUNCTION CreateIPhreeqcF
    END INTERFACE
    CreateIPhreeqc = CreateIPhreeqcF()
    return
END FUNCTION CreateIPhreeqc

INTEGER FUNCTION DestroyIPhreeqc(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION DestroyIPhreeqcF(id) &
            BIND(C, NAME='DestroyIPhreeqcF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION DestroyIPhreeqcF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    DestroyIPhreeqc = DestroyIPhreeqcF(id)
    return
END FUNCTION DestroyIPhreeqc

INTEGER FUNCTION GetComponentCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetComponentCountF(id) &
            BIND(C, NAME='GetComponentCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetComponentCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetComponentCount = GetComponentCountF(id)
    return
END FUNCTION GetComponentCount

SUBROUTINE GetComponent(id, n, comp)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetComponentF(id, n, comp, l) &
            BIND(C, NAME='GetComponentF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: comp(*)
        END SUBROUTINE GetComponentF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: comp
    call GetComponentF(id, n, comp, len(comp))
    return
END SUBROUTINE GetComponent

INTEGER FUNCTION GetCurrentSelectedOutputUserNumber(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetCurrentSelectedOutputUserNumberF(id) &
            BIND(C, NAME='GetCurrentSelectedOutputUserNumberF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetCurrentSelectedOutputUserNumberF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetCurrentSelectedOutputUserNumber = GetCurrentSelectedOutputUserNumberF(id)
    return
END FUNCTION GetCurrentSelectedOutputUserNumber

SUBROUTINE GetDumpFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetDumpFileNameF(id, fname, l) &
            BIND(C, NAME='GetDumpFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: fname(*)
        END SUBROUTINE GetDumpFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: fname
    call GetDumpFileNameF(id, fname, len(fname))
    return
END SUBROUTINE GetDumpFileName

LOGICAL FUNCTION GetDumpFileOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetDumpFileOnF(id) &
            BIND(C, NAME='GetDumpFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetDumpFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetDumpFileOn = (GetDumpFileOnF(id) .ne. 0)
    return
END FUNCTION GetDumpFileOn

INTEGER FUNCTION GetDumpStringLineCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetDumpStringLineCountF(id) &
            BIND(C, NAME='GetDumpStringLineCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetDumpStringLineCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetDumpStringLineCount = GetDumpStringLineCountF(id)
    return
END FUNCTION GetDumpStringLineCount
    
SUBROUTINE GetDumpStringLine(id, n, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetDumpStringLineF(id, n, line, l) &
            BIND(C, NAME='GetDumpStringLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: line(*)
        END SUBROUTINE GetDumpStringLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: line
    call GetDumpStringLineF(id, n, line, len(line))
    return
END SUBROUTINE GetDumpStringLine


LOGICAL FUNCTION GetDumpStringOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetDumpStringOnF(id) &
            BIND(C, NAME='GetDumpStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetDumpStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetDumpStringOn = (GetDumpStringOnF(id) .ne. 0)
    return
END FUNCTION GetDumpStringOn

SUBROUTINE GetErrorFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetErrorFileNameF(id, fname, l) &
            BIND(C, NAME='GetErrorFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: fname(*)
        END SUBROUTINE GetErrorFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: fname
    call GetErrorFileNameF(id, fname, len(fname))
    return
END SUBROUTINE GetErrorFileName

LOGICAL FUNCTION GetErrorFileOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetErrorFileOnF(id) &
            BIND(C, NAME='GetErrorFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetErrorFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetErrorFileOn = (GetErrorFileOnF(id) .ne. 0)
    return
END FUNCTION GetErrorFileOn

INTEGER FUNCTION GetErrorStringLineCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetErrorStringLineCountF(id) &
            BIND(C, NAME='GetErrorStringLineCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetErrorStringLineCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetErrorStringLineCount = GetErrorStringLineCountF(id)
    return
END FUNCTION GetErrorStringLineCount

SUBROUTINE GetErrorStringLine(id, n, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetErrorStringLineF(id, n, line, l) &
            BIND(C, NAME='GetErrorStringLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: line(*)
        END SUBROUTINE GetErrorStringLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: line
    call GetErrorStringLineF(id, n, line, len(line))
    return
END SUBROUTINE GetErrorStringLine

LOGICAL FUNCTION GetErrorStringOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetErrorStringOnF(id) &
            BIND(C, NAME='GetErrorStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetErrorStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetErrorStringOn = (GetErrorStringOnF(id) .ne. 0)
    return
END FUNCTION GetErrorStringOn

SUBROUTINE GetLogFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetLogFileNameF(id, fname, l) &
            BIND(C, NAME='GetLogFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: fname(*)
        END SUBROUTINE GetLogFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: fname
    call GetLogFileNameF(id, fname, len(fname))
    return
END SUBROUTINE GetLogFileName

LOGICAL FUNCTION GetLogFileOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetLogFileOnF(id) &
            BIND(C, NAME='GetLogFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetLogFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetLogFileOn = (GetLogFileOnF(id) .ne. 0)
    return
END FUNCTION GetLogFileOn

INTEGER FUNCTION GetLogStringLineCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetLogStringLineCountF(id) &
            BIND(C, NAME='GetLogStringLineCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetLogStringLineCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetLogStringLineCount = GetLogStringLineCountF(id)
    return
END FUNCTION GetLogStringLineCount

SUBROUTINE GetLogStringLine(id, n, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetLogStringLineF(id, n, line, l) &
            BIND(C, NAME='GetLogStringLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: line(*)
        END SUBROUTINE GetLogStringLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: line
    call GetLogStringLineF(id, n, line, len(line))
    return
END SUBROUTINE GetLogStringLine
   
LOGICAL FUNCTION GetLogStringOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetLogStringOnF(id) &
            BIND(C, NAME='GetLogStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetLogStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetLogStringOn = (GetLogStringOnF(id) .ne. 0)
    return
END FUNCTION GetLogStringOn

INTEGER FUNCTION GetNthSelectedOutputUserNumber(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetNthSelectedOutputUserNumberF(id, n) &
            BIND(C, NAME='GetNthSelectedOutputUserNumberF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n
        END FUNCTION GetNthSelectedOutputUserNumberF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    GetNthSelectedOutputUserNumber = GetNthSelectedOutputUserNumberF(id, n)
    return
END FUNCTION GetNthSelectedOutputUserNumber

SUBROUTINE GetOutputFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetOutputFileNameF(id, fname, l) &
            BIND(C, NAME='GetOutputFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: fname(*)
        END SUBROUTINE GetOutputFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: fname
    call GetOutputFileNameF(id, fname, len(fname))
    return
END SUBROUTINE GetOutputFileName

LOGICAL FUNCTION GetOutputFileOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetOutputFileOnF(id) &
            BIND(C, NAME='GetOutputFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetOutputFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetOutputFileOn = (GetOutputFileOnF(id) .ne. 0)
    return
END FUNCTION GetOutputFileOn

INTEGER FUNCTION GetOutputStringLineCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetOutputStringLineCountF(id) &
            BIND(C, NAME='GetOutputStringLineCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetOutputStringLineCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetOutputStringLineCount = GetOutputStringLineCountF(id)
    return
END FUNCTION GetOutputStringLineCount

SUBROUTINE GetOutputStringLine(id, n, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetOutputStringLineF(id, n, line, l) &
            BIND(C, NAME='GetOutputStringLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: line(*)
        END SUBROUTINE GetOutputStringLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: line
    call GetOutputStringLineF(id, n, line, len(line))
    return
END SUBROUTINE GetOutputStringLine

LOGICAL FUNCTION GetOutputStringOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetOutputStringOnF(id) &
            BIND(C, NAME='GetOutputStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetOutputStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetOutputStringOn = (GetOutputStringOnF(id) .ne. 0)
    return
END FUNCTION GetOutputStringOn

INTEGER FUNCTION GetSelectedOutputColumnCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputColumnCountF(id) &
            BIND(C, NAME='GetSelectedOutputColumnCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetSelectedOutputColumnCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetSelectedOutputColumnCount = GetSelectedOutputColumnCountF(id)
    return
END FUNCTION GetSelectedOutputColumnCount

INTEGER FUNCTION GetSelectedOutputCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputCountF(id) &
            BIND(C, NAME='GetSelectedOutputCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetSelectedOutputCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetSelectedOutputCount = GetSelectedOutputCountF(id)
    return
END FUNCTION GetSelectedOutputCount

SUBROUTINE GetSelectedOutputFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetSelectedOutputFileNameF(id, fname, l) &
            BIND(C, NAME='GetSelectedOutputFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: fname(*)
        END SUBROUTINE GetSelectedOutputFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(inout) :: fname
    call GetSelectedOutputFileNameF(id, fname, len(fname))
    return
END SUBROUTINE GetSelectedOutputFileName

LOGICAL FUNCTION GetSelectedOutputFileOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputFileOnF(id) &
            BIND(C, NAME='GetSelectedOutputFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetSelectedOutputFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetSelectedOutputFileOn = (GetSelectedOutputFileOnF(id) .ne. 0)
    return
END FUNCTION GetSelectedOutputFileOn

INTEGER FUNCTION GetSelectedOutputStringLineCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputStringLineCountF(id) &
            BIND(C, NAME='GetSelectedOutputStringLineCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetSelectedOutputStringLineCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetSelectedOutputStringLineCount = GetSelectedOutputStringLineCountF(id)
    return
END FUNCTION GetSelectedOutputStringLineCount

SUBROUTINE GetSelectedOutputStringLine(id, n, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetSelectedOutputStringLineF(id, n, line, l) &
            BIND(C, NAME='GetSelectedOutputStringLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: line(*)
        END SUBROUTINE GetSelectedOutputStringLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: line
    call GetSelectedOutputStringLineF(id, n, line, len(line))
    return
END SUBROUTINE GetSelectedOutputStringLine

LOGICAL FUNCTION GetSelectedOutputStringOn(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputStringOnF(id) &
            BIND(C, NAME='GetSelectedOutputStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetSelectedOutputStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetSelectedOutputStringOn = (GetSelectedOutputStringOnF(id) .ne. 0)
    return
END FUNCTION GetSelectedOutputStringOn

INTEGER FUNCTION GetSelectedOutputRowCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputRowCountF(id) &
            BIND(C, NAME='GetSelectedOutputRowCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetSelectedOutputRowCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetSelectedOutputRowCount = GetSelectedOutputRowCountF(id)
    return
END FUNCTION GetSelectedOutputRowCount

INTEGER FUNCTION GetSelectedOutputValue(id, row, col, vtype, dvalue, svalue)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, l) &
            BIND(C, NAME='GetSelectedOutputValueF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, row, col, l
            INTEGER(KIND=C_INT), INTENT(out) :: vtype 
            REAL(KIND=C_DOUBLE), INTENT(out) :: dvalue
            CHARACTER(KIND=C_CHAR), INTENT(out) :: svalue(*)
        END FUNCTION GetSelectedOutputValueF
    END INTERFACE
    INTEGER, INTENT(in) :: id, row, col
    INTEGER, INTENT(out) :: vtype
    DOUBLE PRECISION, INTENT(out) :: dvalue
    CHARACTER(len=*), INTENT(out) :: svalue
    INTEGER :: l
    l = len(svalue)
    GetSelectedOutputValue = GetSelectedOutputValueF(id, row, col, vtype, dvalue, svalue, l)
    return
END FUNCTION GetSelectedOutputValue

SUBROUTINE GetVersionString(version)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetVersionStringF(version, l) &
            BIND(C, NAME='GetVersionStringF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: version(*)
        END SUBROUTINE GetVersionStringF
    END INTERFACE
    CHARACTER(len=*), INTENT(inout) :: version
    call GetVersionStringF(version, len(version))
    return
END SUBROUTINE GetVersionString

INTEGER FUNCTION GetWarningStringLineCount(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION GetWarningStringLineCountF(id) &
            BIND(C, NAME='GetWarningStringLineCountF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION GetWarningStringLineCountF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    GetWarningStringLineCount = GetWarningStringLineCountF(id)
    return
END FUNCTION GetWarningStringLineCount


SUBROUTINE GetWarningStringLine(id, n, line)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE GetWarningStringLineF(id, n, line, l) &
            BIND(C, NAME='GetWarningStringLineF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n, l
            CHARACTER(KIND=C_CHAR), INTENT(out) :: line(*)
        END SUBROUTINE GetWarningStringLineF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    CHARACTER(len=*), INTENT(inout) :: line
    call GetWarningStringLineF(id, n, line, len(line))
    return
END SUBROUTINE GetWarningStringLine

INTEGER FUNCTION LoadDatabase(id, filename)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION LoadDatabaseF(id, filename) &
            BIND(C, NAME='LoadDatabaseF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: filename(*)
        END FUNCTION LoadDatabaseF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: filename
    LoadDatabase = LoadDatabaseF(id, trim(filename)//C_NULL_CHAR)
    return
END FUNCTION LoadDatabase

INTEGER FUNCTION LoadDatabaseString(id, input)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION LoadDatabaseStringF(id, input) &
            BIND(C, NAME='LoadDatabaseStringF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: input(*)
        END FUNCTION LoadDatabaseStringF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: input
    LoadDatabaseString = LoadDatabaseStringF(id, trim(input)//C_NULL_CHAR)
    return
END FUNCTION LoadDatabaseString

SUBROUTINE OutputAccumulatedLines(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE OutputAccumulatedLinesF(id) &
            BIND(C, NAME='OutputAccumulatedLinesF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END SUBROUTINE OutputAccumulatedLinesF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    call OutputAccumulatedLinesF(id)
    return
END SUBROUTINE OutputAccumulatedLines

SUBROUTINE OutputErrorString(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE OutputErrorStringF(id) &
            BIND(C, NAME='OutputErrorStringF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END SUBROUTINE OutputErrorStringF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    call OutputErrorStringF(id)
    return
END SUBROUTINE OutputErrorString

SUBROUTINE OutputWarningString(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        SUBROUTINE OutputWarningStringF(id) &
            BIND(C, NAME='OutputWarningStringF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END SUBROUTINE OutputWarningStringF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    call OutputWarningStringF(id)
    return
END SUBROUTINE OutputWarningString

INTEGER FUNCTION RunAccumulated(id)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RunAccumulatedF(id) &
            BIND(C, NAME='RunAccumulatedF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
        END FUNCTION RunAccumulatedF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    RunAccumulated = RunAccumulatedF(id)
    return
END FUNCTION RunAccumulated

INTEGER FUNCTION RunFile(id, filename)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RunFileF(id, filename) &
            BIND(C, NAME='RunFileF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: filename(*)
        END FUNCTION RunFileF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: filename
    RunFile = RunFileF(id, trim(filename)//C_NULL_CHAR)
    return
END FUNCTION RunFile

INTEGER FUNCTION RunString(id, input)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION RunStringF(id, input) &
            BIND(C, NAME='RunStringF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: input(*)
        END FUNCTION RunStringF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: input
    RunString = RunStringF(id, trim(input)//C_NULL_CHAR)
    return
END FUNCTION RunString
#ifdef IPHREEQC_NO_FORTRAN_MODULE
INTEGER FUNCTION SetBasicFortranCallback(id, fcn)
    INTERFACE
        INTEGER FUNCTION SetBasicFortranCallbackF(id, fcn)
            IMPLICIT NONE
            INTEGER, INTENT(in) :: id
            INTERFACE
                DOUBLE PRECISION FUNCTION fcn(x1, x2, str, l)
                    INTEGER, INTENT(in)   :: l
                    DOUBLE PRECISION, INTENT(in) :: x1, x2
                    CHARACTER, INTENT(in) :: str(*)
                END FUNCTION fcn
           END INTERFACE
        END FUNCTION SetBasicFortranCallbackF
    END INTERFACE
    INTEGER, INTENT(in) :: id     
    INTERFACE
        DOUBLE PRECISION FUNCTION fcn(x1, x2, str, l) 
            INTEGER, INTENT(in)          :: l
            DOUBLE PRECISION, INTENT(in) :: x1, x2
            CHARACTER, INTENT(in) :: str(*)
        END FUNCTION fcn
    END INTERFACE
    SetBasicFortranCallback = SetBasicFortranCallbackF(id, fcn)
    return
END FUNCTION SetBasicFortranCallback
#else
INTEGER FUNCTION SetBasicFortranCallback(id, fcn)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetBasicFortranCallbackF(id, fcn) &
            BIND(C, NAME='SetBasicFortranCallbackF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            INTERFACE
                REAL(KIND=C_DOUBLE) FUNCTION fcn(x1, x2, str, l) BIND(C)
                    USE ISO_C_BINDING
                    IMPLICIT none
                    REAL(KIND=C_DOUBLE),    INTENT(in) :: x1, x2
                    CHARACTER(KIND=C_CHAR), INTENT(in) :: str(*)
                    INTEGER(KIND=C_INT),    INTENT(in), value :: l
                END FUNCTION fcn
           END INTERFACE
        END FUNCTION SetBasicFortranCallbackF
    END INTERFACE
    INTEGER, INTENT(in) :: id     
    INTERFACE
        REAL(KIND=C_DOUBLE) FUNCTION fcn(x1, x2, str, l) BIND(C)
            USE ISO_C_BINDING
            IMPLICIT none
            REAL(KIND=C_DOUBLE),    INTENT(in) :: x1, x2
            CHARACTER(KIND=C_CHAR), INTENT(in) :: str(*)
            INTEGER(KIND=C_INT),    INTENT(in), value :: l
        END FUNCTION fcn
    END INTERFACE
    SetBasicFortranCallback = SetBasicFortranCallbackF(id, fcn)
    return
END FUNCTION SetBasicFortranCallback
#endif

INTEGER FUNCTION SetCurrentSelectedOutputUserNumber(id, n)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetCurrentSelectedOutputUserNumberF(id, n) &
            BIND(C, NAME='SetCurrentSelectedOutputUserNumberF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, n
        END FUNCTION SetCurrentSelectedOutputUserNumberF
    END INTERFACE
    INTEGER, INTENT(in) :: id, n
    SetCurrentSelectedOutputUserNumber = SetCurrentSelectedOutputUserNumberF(id, n)
    return
END FUNCTION SetCurrentSelectedOutputUserNumber

INTEGER FUNCTION SetDumpFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetDumpFileNameF(id, fname) &
            BIND(C, NAME='SetDumpFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: fname(*)
        END FUNCTION SetDumpFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: fname
    SetDumpFileName = SetDumpFileNameF(id, trim(fname)//C_NULL_CHAR)
    return
END FUNCTION SetDumpFileName

INTEGER FUNCTION SetDumpFileOn(id, dump_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetDumpFileOnF(id, dump_on) &
            BIND(C, NAME='SetDumpFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, dump_on
        END FUNCTION SetDumpFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: dump_on
    INTEGER :: tf = 0
    tf = 0
    if (dump_on) tf = 1
    SetDumpFileOn = SetDumpFileOnF(id, tf)
    return
END FUNCTION SetDumpFileOn

INTEGER FUNCTION SetDumpStringOn(id, dump_string_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetDumpStringOnF(id, dump_string_on) &
            BIND(C, NAME='SetDumpStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, dump_string_on
        END FUNCTION SetDumpStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: dump_string_on
    INTEGER :: tf = 0
    tf = 0
    if (dump_string_on) tf = 1
    SetDumpStringOn = SetDumpStringOnF(id, tf)
    return
END FUNCTION SetDumpStringOn

INTEGER FUNCTION SetErrorFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetErrorFileNameF(id, fname) &
            BIND(C, NAME='SetErrorFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: fname(*)
        END FUNCTION SetErrorFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: fname
    SetErrorFileName = SetErrorFileNameF(id, trim(fname)//C_NULL_CHAR)
    return
END FUNCTION SetErrorFileName

INTEGER FUNCTION SetErrorFileOn(id, error_file_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetErrorFileOnF(id, error_file_on) &
            BIND(C, NAME='SetErrorFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, error_file_on
        END FUNCTION SetErrorFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: error_file_on
    INTEGER :: tf = 0
    tf = 0
    if (error_file_on) tf = 1
    SetErrorFileOn = SetErrorFileOnF(id, tf)
    return
END FUNCTION SetErrorFileOn

INTEGER FUNCTION SetErrorStringOn(id, error_string_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetErrorStringOnF(id, error_string_on) &
            BIND(C, NAME='SetErrorStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, error_string_on
        END FUNCTION SetErrorStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: error_string_on
    INTEGER :: tf = 0
    tf = 0
    if (error_string_on) tf = 1
    SetErrorStringOn = SetErrorStringOnF(id, tf)
    return
END FUNCTION SetErrorStringOn

INTEGER FUNCTION SetLogFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetLogFileNameF(id, fname) &
            BIND(C, NAME='SetLogFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: fname(*)
        END FUNCTION SetLogFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: fname
    SetLogFileName = SetLogFileNameF(id, trim(fname)//C_NULL_CHAR)
    return
END FUNCTION SetLogFileName

INTEGER FUNCTION SetLogFileOn(id, log_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetLogFileOnF(id, log_on) &
            BIND(C, NAME='SetLogFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, log_on
        END FUNCTION SetLogFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: log_on
    INTEGER :: tf = 0
    tf = 0
    if (log_on) tf = 1
    SetLogFileOn = SetLogFileOnF(id, tf)
    return
END FUNCTION SetLogFileOn

INTEGER FUNCTION SetLogStringOn(id, log_string_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetLogStringOnF(id, log_string_on) &
            BIND(C, NAME='SetLogStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, log_string_on
        END FUNCTION SetLogStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: log_string_on
    INTEGER :: tf = 0
    tf = 0
    if (log_string_on) tf = 1
    SetLogStringOn = SetLogStringOnF(id, tf)
    return
END FUNCTION SetLogStringOn

INTEGER FUNCTION SetOutputFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetOutputFileNameF(id, fname) &
            BIND(C, NAME='SetOutputFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: fname(*)
        END FUNCTION SetOutputFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: fname
    SetOutputFileName = SetOutputFileNameF(id, trim(fname)//C_NULL_CHAR)
    return
END FUNCTION SetOutputFileName

INTEGER FUNCTION SetOutputFileOn(id, output_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetOutputFileOnF(id, output_on) &
            BIND(C, NAME='SetOutputFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, output_on
        END FUNCTION SetOutputFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: output_on
    INTEGER :: tf
    tf = 0
    if (output_on) tf = 1
    SetOutputFileOn = SetOutputFileOnF(id, tf)
    return
END FUNCTION SetOutputFileOn

INTEGER FUNCTION SetOutputStringOn(id, output_string_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetOutputStringOnF(id, output_string_on) &
            BIND(C, NAME='SetOutputStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, output_string_on
        END FUNCTION SetOutputStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: output_string_on
    INTEGER :: tf = 0
    tf = 0
    if (output_string_on) tf = 1
    SetOutputStringOn = SetOutputStringOnF(id, tf)
    return
END FUNCTION SetOutputStringOn

INTEGER FUNCTION SetSelectedOutputFileName(id, fname)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetSelectedOutputFileNameF(id, fname) &
            BIND(C, NAME='SetSelectedOutputFileNameF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id
            CHARACTER(KIND=C_CHAR), INTENT(in) :: fname(*)
        END FUNCTION SetSelectedOutputFileNameF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: fname
    SetSelectedOutputFileName = SetSelectedOutputFileNameF(id, trim(fname)//C_NULL_CHAR)
    return
END FUNCTION SetSelectedOutputFileName

INTEGER FUNCTION SetSelectedOutputFileOn(id, sel_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetSelectedOutputFileOnF(id, sel_on) &
            BIND(C, NAME='SetSelectedOutputFileOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, sel_on
        END FUNCTION SetSelectedOutputFileOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: sel_on
    INTEGER :: tf = 0
    tf = 0
    if (sel_on) tf = 1
    SetSelectedOutputFileOn = SetSelectedOutputFileOnF(id, tf)
    return
END FUNCTION SetSelectedOutputFileOn

INTEGER FUNCTION SetSelectedOutputStringOn(id, selected_output_string_on)
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTERFACE
        INTEGER(KIND=C_INT) FUNCTION SetSelectedOutputStringOnF(id, selected_output_string_on) &
            BIND(C, NAME='SetSelectedOutputStringOnF')
            USE ISO_C_BINDING
            IMPLICIT NONE
            INTEGER(KIND=C_INT), INTENT(in) :: id, selected_output_string_on
        END FUNCTION SetSelectedOutputStringOnF
    END INTERFACE
    INTEGER, INTENT(in) :: id
    LOGICAL, INTENT(in) :: selected_output_string_on
    INTEGER :: tf = 0
    tf = 0
    if (selected_output_string_on) tf = 1
    SetSelectedOutputStringOn = SetSelectedOutputStringOnF(id, tf)
    return
END FUNCTION SetSelectedOutputStringOn

END MODULE
#endif

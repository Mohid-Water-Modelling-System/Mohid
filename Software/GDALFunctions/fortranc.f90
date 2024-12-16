!    Copyright 2011 Davide Cesari <dcesari69 at gmail dot com>
!
!    This file is part of FortranGIS.
!
!    FortranGIS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    FortranGIS is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with FortranGIS.  If not, see
!    <http://www.gnu.org/licenses/>.
!#include "config.h"

!> Utility module for supporting Fortran 2003 C language interface module.
!! This module contains various utilties for simplifying the exchange
!! of character variables between Fortran and C when using the
!! <tt>ISO_C_BINDING</tt> intrinsic module of Fortran 2003.
!!
!! For an example of application of the \a fortranc module, please
!! refer to the following test program, which, among the other
!! operations, decodes the output of a C function returning a
!! <tt>char**</tt> result:
!! \include fortranc_test.F90
!!
!! \ingroup libfortranc
MODULE fortranc
USE,INTRINSIC :: ISO_C_BINDING
#ifdef WITH_VARYING_STRING
USE iso_varying_string
#endif
IMPLICIT NONE


!> Fortran derived type for handling <tt>void**</tt>, <tt>char**</tt>,
!! etc C objects (pointer to pointer or array of pointers).  The array
!! of pointers is assumed to be terminated by a <tt>NULL</tt>
!! pointer. Each pointer of the array typically points to a
!! null-terminated string, although this is not always the
!! case. Methods are provided both for receiving the data structure
!! from C and unpacking it in Fortran as well as for creating it in
!! Fortran and passing it to C.
!!
!! Example of <tt>char**</tt> object created in C and unpacked in Fortran:
!! \code
!! TYPE(c_ptr_ptr) :: envp
!! INTEGER :: i
!! ...
!! envp = c_ptr_ptr_new(interfaced_c_procedure())
!! DO i = 1, c_ptr_ptr_getsize(envp)
!!   PRINT*,i,TRIM(strtofchar(c_ptr_ptr_getptr(envp, i),100))
!! ENDDO
!! CALL delete(envp)
!! \endcode
!!
!! Example of <tt>char**</tt> object created in Fortran and passed to C:
!! \code
!! TYPE(c_ptr_ptr) :: envp
!! ...
!! envp = c_ptr_ptr_new((/'APPLE=3  ','PEAR=2  ','ORANGE=20'/))
!! CALL interfaced_c_procedure(c_ptr_ptr_getobject(envp))
!! CALL delete(envp)
!! ...
!! \endcode
TYPE c_ptr_ptr
  PRIVATE
  TYPE(c_ptr),POINTER :: elem(:) => NULL()
  CHARACTER(len=1),POINTER :: buffer(:) => NULL()
END TYPE c_ptr_ptr

!> Equivalent of the strlen C function.
!!
!! \param string null-terminated C-style string to test
INTERFACE strlen
  MODULE PROCEDURE strlen_char, strlen_chararr, strlen_intarr, &
   strlen_ptr
#ifdef WITH_VARYING_STRING
  MODULE PROCEDURE strlen_var_str
#endif
END INTERFACE

!> Convert a null-terminated C string into a Fortran <tt>CHARACTER</tt>
!! variable of the proper length.  The input can be provided as a
!! Fortran <tt>CHARACTER</tt> scalar of any length, as a Fortran array
!! of <tt>CHARACTER</tt> of length one, as an array of 1-byte integers or as
!! a C pointer to char (<tt>char*</tt>).
!!
!! It is typically used for:
!!
!!  - converting a string created/modified by a C function and passed
!!    as a <tt>char *</tt> argument, interfaced as
!!    <tt>CHARACTER(kind=c_char,len=*) :: fchar</tt> for its
!!    subsequent use in Fortran
!!
!!  - (more frequently) converting a string returned by a C function
!!    declared as <tt>char*</tt>, interfaced as <tt>TYPE(c_ptr)</tt>
!!    for its subsequent use in Fortran
!!
!!  - converting a string contained in an C-interoperable derived
!!    type, declared in C as <tt>char*</tt>, interfaced as
!!    <tt>TYPE(c_ptr)</tt> for its subsequent use in Fortran
!!
!! \param string null-terminated C-style string to convert
INTERFACE strtofchar
  MODULE PROCEDURE strtofchar_char, strtofchar_chararr, strtofchar_intarr, &
   strtofchar_ptr_2
END INTERFACE

!> Constructor for a \a c_ptr_ptr object.
!! An object of this type can be constructed either from a pointer
!! returned by a C procedure, (either as an argument, interfaced as
!! <tt>TYPE(c_ptr),VALUE</tt> or as the result of a function,
!! interfaced as <tt>TYPE(c_ptr)</tt>) or from Fortran array of
!! character variables (<tt>char **</tt> objects only).
INTERFACE c_ptr_ptr_new
  MODULE PROCEDURE c_ptr_ptr_new_from_c, c_ptr_ptr_new_from_fchar
END INTERFACE c_ptr_ptr_new

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE strtofchararr_assign
END INTERFACE ASSIGNMENT(=)

PRIVATE
PUBLIC strlen, strtofchar, fchartostr, fchartrimtostr, ASSIGNMENT(=)
PUBLIC c_ptr_ptr, c_ptr_ptr_new, c_ptr_ptr_getsize, c_ptr_ptr_getptr, c_ptr_ptr_getobject

CONTAINS


PURE FUNCTION strlen_char(string) RESULT(strlen)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strlen_char
#endif
CHARACTER(kind=c_char,len=*),INTENT(in) :: string
INTEGER :: strlen

INTEGER :: i

DO i = 1, LEN(string)
  IF (string(i:i) == CHAR(0)) EXIT
ENDDO
strlen = i - 1

END FUNCTION strlen_char


PURE FUNCTION strlen_chararr(string) RESULT(strlen)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strlen_chararr
#endif
CHARACTER(kind=c_char,len=1),INTENT(in) :: string(:)
INTEGER :: strlen

INTEGER :: i

DO i = 1, SIZE(string)
  IF (string(i) == CHAR(0)) EXIT
ENDDO
strlen = i - 1

END FUNCTION strlen_chararr


PURE FUNCTION strlen_intarr(string) RESULT(strlen)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strlen_intarr
#endif
INTEGER(kind=c_signed_char),INTENT(in) :: string(:)
INTEGER :: strlen

INTEGER :: i

DO i = 1, SIZE(string)
  IF (string(i) == 0) EXIT
ENDDO
strlen = i - 1

END FUNCTION strlen_intarr


FUNCTION strlen_ptr(string) RESULT(strlen)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strlen_ptr
#endif
TYPE(c_ptr),INTENT(in) :: string
INTEGER :: strlen

INTEGER(kind=c_signed_char),POINTER :: pstring(:)
INTEGER :: i

IF (C_ASSOCIATED(string)) THEN ! conflicts with PURE
! null C pointer does not produce unassociated Fortran pointer with Intel
  CALL C_F_POINTER(string, pstring, (/HUGE(i)/))
! IF (ASSOCIATED(pstring)) THEN
  DO i = 1, SIZE(pstring)
    IF (pstring(i) == 0) EXIT
  ENDDO
  strlen = i - 1
ELSE
  strlen = 0
ENDIF

END FUNCTION strlen_ptr


#ifdef WITH_VARYING_STRING
PURE FUNCTION strlen_var_str(string) RESULT(strlen)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strlen_var_str
#endif
TYPE(varying_string),INTENT(in) :: string
INTEGER :: strlen

strlen = len(string)

END FUNCTION strlen_var_str
#endif


FUNCTION strtofchar_char(string) RESULT(fchar)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strtofchar_char
#endif
CHARACTER(kind=c_char,len=*),INTENT(in) :: string
CHARACTER(len=strlen(string)) :: fchar

fchar(:) = string(1:LEN(fchar))

END FUNCTION strtofchar_char


FUNCTION strtofchar_chararr(string) RESULT(fchar)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strtofchar_chararr
#endif
CHARACTER(kind=c_char,len=1),INTENT(in) :: string(:)
CHARACTER(len=strlen(string)) :: fchar

INTEGER :: i

DO i = 1, LEN(fchar)
  fchar(i:i) = string(i)
ENDDO

END FUNCTION strtofchar_chararr


FUNCTION strtofchar_intarr(string) RESULT(fchar)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strtofchar_intarr
#endif
INTEGER(kind=c_signed_char),INTENT(in) :: string(:)
CHARACTER(len=strlen(string)) :: fchar

fchar(:) = TRANSFER(string(1:LEN(fchar)), fchar)

END FUNCTION strtofchar_intarr


! this unfortunately works only with gfortran where c_f_pointer is
! "erroneously" declared as PURE thus strlen_ptr can be PURE as well

!FUNCTION strtofchar_ptr(string) RESULT(fchar)
!TYPE(c_ptr),INTENT(in) :: string
!CHARACTER(len=strlen(string)) :: fchar
!
!CHARACTER(len=strlen(string)),POINTER :: pfchar
!
!IF (C_ASSOCIATED(string)) THEN
!  CALL c_f_pointer(string, pfchar)
!  fchar(:) = pfchar(:)
!!ELSE
!! silently return an empty string probably useless because
!! strlen is zero in this case (to be tested)
!!  fchar = ''
!ENDIF
!
!END FUNCTION strtofchar_ptr


FUNCTION strtofchar_ptr_2(string, fixlen) RESULT(fchar)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strtofchar_ptr_2
#endif
TYPE(c_ptr),INTENT(in) :: string
INTEGER,INTENT(in) :: fixlen
CHARACTER(len=fixlen) :: fchar

CHARACTER(len=fixlen),POINTER :: pfchar
INTEGER :: safelen

safelen = MIN(strlen(string), fixlen)

fchar = ''
IF (C_ASSOCIATED(string)) THEN
  CALL c_f_pointer(string, pfchar)
  fchar(1:safelen) = pfchar(1:safelen)
ENDIF

END FUNCTION strtofchar_ptr_2


!> Convert a Fortran \a CHARACTER variable into a null-terminated C
!! string. The result is still of type \a CHARACTER but it is
!! interoperable with a C null-terminated string argument <tt>const
!! char*</tt> interfaced as <tt>CHARACTER(kind=c_char) :: cstr</tt>.
FUNCTION fchartostr(fchar) RESULT(string)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: fchartostr
#endif
CHARACTER(len=*),INTENT(in) :: fchar !< Fortran \a CHARACTER variable to convert
CHARACTER(kind=c_char,len=LEN(fchar)+1) :: string

string = fchar//CHAR(0)

END FUNCTION fchartostr


!> Trim trailing blanks and convert a Fortran \a CHARACTER variable
!! into a null-terminated C string. The result is still of type \a
!! CHARACTER but it is interoperable with a C null-terminated string
!! argument <tt>const char*</tt> interfaced as
!! <tt>CHARACTER(kind=c_char) :: cstr</tt>.
FUNCTION fchartrimtostr(fchar) RESULT(string)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: fchartrimtostr
#endif
CHARACTER(len=*),INTENT(in) :: fchar !< Fortran \a CHARACTER variable to convert
CHARACTER(kind=c_char,len=LEN_TRIM(fchar)+1) :: string

string = TRIM(fchar)//CHAR(0)

END FUNCTION fchartrimtostr


SUBROUTINE strtofchararr_assign(fchar, string)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: strtofchararr_assign
#endif
CHARACTER(kind=c_char,len=1),ALLOCATABLE,INTENT(out) :: fchar(:)
TYPE(c_ptr),INTENT(in) :: string

CHARACTER(kind=c_char),POINTER :: pstring(:)
INTEGER :: l

l = strlen(string)
CALL C_F_POINTER(string, pstring, (/l/))
ALLOCATE(fchar(l))
fchar(:) = pstring(:)

END SUBROUTINE strtofchararr_assign


!> Constructor for a \a c_ptr_ptr object.
!! The argument, a generic C pointer, must be a C array of pointers
!! (<tt>char** c_ptr_ptr_c</tt> or <tt>char* c_ptr_ptr_c[n]</tt>),
!! typically the result of a C function. The resulting object can be
!! queried by means of the \a c_ptr_ptr_getsize and \a
!! c_ptr_ptr_getptr methods, but it should not be modified by Fortran.
FUNCTION c_ptr_ptr_new_from_c(c_ptr_ptr_c) RESULT(this)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: c_ptr_ptr_new_from_c
#endif
TYPE(c_ptr),VALUE :: c_ptr_ptr_c !< pointer returned by a C procedure
TYPE(c_ptr_ptr) :: this

INTEGER :: i
TYPE(c_ptr),POINTER :: charp(:)

IF (C_ASSOCIATED(c_ptr_ptr_c)) THEN
  ! HUGE() here is ugly, but we must set a finite size
  CALL C_F_POINTER(c_ptr_ptr_c, charp, (/HUGE(1)/))
  DO i = 1, SIZE(charp)
    IF (.NOT.C_ASSOCIATED(charp(i))) THEN
      CALL C_F_POINTER(c_ptr_ptr_c, this%elem, (/i/))
      RETURN
    ENDIF
  ENDDO
ENDIF
END FUNCTION c_ptr_ptr_new_from_c


!> Constructor for a \a c_ptr_ptr object.
!! The argument is an array of Fortran character variables which will
!! be trimmed and stored in the resulting object. The object can be
!! passed to a C procedure as a <tt>char **</tt> argument after
!! applying the \a c_ptr_prt_getptr method, but it should not be
!! modified by the C procedure.
FUNCTION c_ptr_ptr_new_from_fchar(fchar) RESULT(this)
CHARACTER(len=*) :: fchar(:) !< array of characters that will compose the object
TYPE(c_ptr_ptr) :: this

INTEGER :: i, j, totlen

totlen = 0
DO i = 1, SIZE(fchar)
  totlen = totlen + LEN_TRIM(fchar(i)) + 1
ENDDO
ALLOCATE(this%buffer(totlen), this%elem(SIZE(fchar) + 1))
totlen = 1
DO i = 1, SIZE(fchar)
  this%elem(i) = C_LOC(this%buffer(totlen))
  DO j = 1, LEN_TRIM(fchar(i))
    this%buffer(totlen) = fchar(i)(j:j)
    totlen = totlen + 1
  ENDDO
  this%buffer(totlen) = CHAR(0)
  totlen = totlen + 1
ENDDO
this%elem(i) = C_NULL_PTR

END FUNCTION c_ptr_ptr_new_from_fchar


!> Return the number of valid pointers in the array pointer \a this.
!! If the object has not been initialized or has been initialized with
!! errors, zero is returned.
FUNCTION c_ptr_ptr_getsize(this)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: c_ptr_ptr_getsize
#endif
TYPE(c_ptr_ptr),INTENT(in) :: this
INTEGER :: c_ptr_ptr_getsize

IF (ASSOCIATED(this%elem)) THEN
  c_ptr_ptr_getsize = SIZE(this%elem) - 1
ELSE
  c_ptr_ptr_getsize = 0
ENDIF

END FUNCTION c_ptr_ptr_getsize

!> Return the n-th pointer in the array pointer \a this.
!! Ths method is useful if the object \a this has been created from C.
!! If the object has not been initialized, or \a n is out of bounds, a
!! NULL pointer is returned, this condition can be checked by means of
!! the <tt>C_ASSOCIATED()</tt> function. If \a this is an array of
!! pointers to C null-terminated strings, the string can be returned
!! as a Fortran \a CHARACTER variable of the proper length by using
!! the \a strtofchar function.
FUNCTION c_ptr_ptr_getptr(this, n)
#ifdef DLL_EXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: c_ptr_ptr_getptr
#endif
TYPE(c_ptr_ptr),INTENT(in) :: this !< object to query
INTEGER,INTENT(in) :: n !< the number of pointer to get (starting from 1)
TYPE(c_ptr) :: c_ptr_ptr_getptr

c_ptr_ptr_getptr = C_NULL_PTR
IF (ASSOCIATED(this%elem)) THEN
  IF (n > 0 .AND. n <= SIZE(this%elem)) THEN
    c_ptr_ptr_getptr = this%elem(n)
  ENDIF
ENDIF

END FUNCTION c_ptr_ptr_getptr


!> Return the C pointer to the first pointer in the array pointer \a this.
!! This method is useful if the object \a this has been created from
!! Fortran and it has to be passed to a C procedure.
FUNCTION c_ptr_ptr_getobject(this)
TYPE(c_ptr_ptr),INTENT(in) :: this !< object to query
TYPE(c_ptr) :: c_ptr_ptr_getobject

c_ptr_ptr_getobject = C_NULL_PTR
IF (ASSOCIATED(this%elem)) THEN
  c_ptr_ptr_getobject = C_LOC(this%elem(1))
ENDIF

END FUNCTION c_ptr_ptr_getobject

END MODULE fortranc

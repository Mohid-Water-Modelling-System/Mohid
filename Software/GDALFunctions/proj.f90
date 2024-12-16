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

!> Fortran 2003 interface to the proj.4 https://github.com/OSGeo/proj.4 library.
!! The following functions are directly interfaced to their
!! corresponding C version, so they are undocumented here, please
!! refer to the original proj C API documentation, e.g. at the address
!! https://github.com/OSGeo/proj.4/wiki/ProjAPI , for their use:
!!  - pj_init_plus() -> FUNCTION pj_init_plus()
!!  - pj_transform()() -> FUNCTION pj_transform()()
!!  - pj_datum_transform() -> FUNCTION pj_datum_transform()
!!  - pj_fwd() -> FUNCTION pj_fwd()
!!  - pj_inv() -> FUNCTION pj_inv()
!!  - pj_geocentric_to_geodetic() -> FUNCTION pj_geocentric_to_geodetic()
!!  - pj_geodetic_to_geocentric() -> FUNCTION pj_geodetic_to_geocentric()
!!  - pj_compare_datums() -> FUNCTION pj_compare_datums()
!!  - pj_is_latlong() -> FUNCTION pj_is_latlong()
!!  - pj_is_geocent() -> FUNCTION pj_is_geocent()
!!  - pj_get_def() -> FUNCTION pj_get_def()
!!  - pj_latlong_from_proj() -> FUNCTION pj_latlong_from_proj()
!!  - pj_free() -> SUBROUTINE pj_free()
!!
!! Notice that, if relevant, the result of functions returning an
!! integer has to be interpreted as 0=false, nonzero=true or 0=ok,
!! nonzero=error.
!!
!! Some of these functions have also a more Fortran-friendly interface
!! explicitely documented here, with an \a _f appended to the name.
!!
!! For an example of application of the \a proj module, please refer
!! to the following test program, which performs a forward and
!! backward transformation:
!! \include proj_test.F90
!!
!! \ingroup libfortrangis
MODULE proj
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE

!> Object describing a cartographic projection.
!! Its components are private so they should not be manipulated
!! directly.
TYPE,BIND(C) :: pj_object
  PRIVATE
  TYPE(c_ptr) :: ptr = C_NULL_PTR
END TYPE pj_object

!> Object representing a null cartographic projection.
!! It should be used for comparison of function results
!! with the \a == operator for error checking.
TYPE(pj_object),PARAMETER :: pj_object_null=pj_object(C_NULL_PTR)

!> Object describing a coordinate pair.
!! It can indicate either projected coordinate (u=x, v=y) or
!! spherical coordinates (u=lon, v=lat).
TYPE,BIND(C) :: pjuv_object
  REAL(kind=c_double) :: u=HUGE(1.0_c_double)
  REAL(kind=c_double) :: v=HUGE(1.0_c_double)
END TYPE pjuv_object

REAL(kind=c_double),PARAMETER :: pj_deg_to_rad=.0174532925199432958D0 !< equivalent to the C api symbol DEG_TO_RAD
REAL(kind=c_double),PARAMETER :: pj_rad_to_deg=57.29577951308232D0 !< equivalent to the C api symbol RAD_TO_DEG

!> Test whether an opaque object is valid.
!! Please use the interface name pj_associated, but for the documentation
!! see the specific function pj_associated_object.
INTERFACE pj_associated
  MODULE PROCEDURE pj_associated_object
END INTERFACE pj_associated

!> Initialize a projection from a string.
!! It returns an object of pj_object type.
INTERFACE
  FUNCTION pj_init_plus(name) BIND(C,name='pj_init_plus')
  IMPORT
  CHARACTER(kind=c_char) :: name(*) !< Projection string, must be terminated by //CHAR(0)
  TYPE(pj_object) :: pj_init_plus
  END FUNCTION pj_init_plus
END INTERFACE

INTERFACE
  FUNCTION pj_transform(src, dst, point_count, point_offset, x, y, z) &
   BIND(C,name='pj_transform')
  IMPORT
  TYPE(pj_object),VALUE :: src
  TYPE(pj_object),VALUE :: dst
  INTEGER(kind=c_long),VALUE :: point_count
  INTEGER(kind=c_int),VALUE :: point_offset
  TYPE(C_PTR), VALUE :: x !< Must be C-pointer to an array real(kind=c_double) :: x(:)
  TYPE(C_PTR), VALUE :: y !< Must be C-pointer to an array real(kind=c_double) :: y(:)
  TYPE(C_PTR), VALUE :: z !< Must be C-pointer to an array real(kind=c_double) :: z(:)
  INTEGER(kind=c_int) :: pj_transform
  END FUNCTION pj_transform
END INTERFACE

INTERFACE
  FUNCTION pj_datum_transform(src, dst, point_count, point_offset, x, y, z) &
   BIND(C,name='pj_datum_transform')
  IMPORT
  TYPE(pj_object),VALUE :: src
  TYPE(pj_object),VALUE :: dst
  INTEGER(kind=c_long),VALUE :: point_count
  INTEGER(kind=c_int),VALUE :: point_offset
  REAL(kind=c_double) :: x(*)
  REAL(kind=c_double) :: y(*)
  REAL(kind=c_double) :: z(*)
  INTEGER(kind=c_int) :: pj_datum_transform
  END FUNCTION pj_datum_transform
END INTERFACE

INTERFACE
  FUNCTION pj_fwd(val, proj) BIND(C,name='pj_fwd')
  IMPORT
  TYPE(pjuv_object),VALUE :: val
  TYPE(pj_object),VALUE :: proj
  TYPE(pjuv_object) :: pj_fwd
  END FUNCTION pj_fwd
END INTERFACE

INTERFACE
  FUNCTION pj_inv(val, proj) BIND(C,name='pj_inv')
  IMPORT
  TYPE(pjuv_object),VALUE :: val
  TYPE(pj_object),VALUE :: proj
  TYPE(pjuv_object) :: pj_inv
  END FUNCTION pj_inv
END INTERFACE

INTERFACE
  FUNCTION pj_geocentric_to_geodetic(a, es, point_count, point_offset, x, y, z) &
   BIND(C,name='pj_geocentric_to_geodetic')
  IMPORT
  REAL(kind=c_double),VALUE :: a
  REAL(kind=c_double),VALUE :: es
  INTEGER(kind=c_long),VALUE :: point_count
  INTEGER(kind=c_int),VALUE :: point_offset
  REAL(kind=c_double) :: x(*)
  REAL(kind=c_double) :: y(*)
  REAL(kind=c_double) :: z(*)
  INTEGER(kind=c_int) :: pj_geocentric_to_geodetic
  END FUNCTION pj_geocentric_to_geodetic
END INTERFACE

INTERFACE
  FUNCTION pj_geodetic_to_geocentric(a, es, point_count, point_offset, x, y, z) &
   BIND(C,name='pj_geodetic_to_geocentric')
  IMPORT
  REAL(kind=c_double),VALUE :: a
  REAL(kind=c_double),VALUE :: es
  INTEGER(kind=c_long),VALUE :: point_count
  INTEGER(kind=c_int),VALUE :: point_offset
  REAL(kind=c_double) :: x(*)
  REAL(kind=c_double) :: y(*)
  REAL(kind=c_double) :: z(*)
  INTEGER(kind=c_int) :: pj_geodetic_to_geocentric
  END FUNCTION pj_geodetic_to_geocentric
END INTERFACE

INTERFACE
  FUNCTION pj_compare_datums(srcdefn, dstdefn) BIND(C,name='pj_compare_datums')
  IMPORT
  TYPE(pj_object),VALUE :: srcdefn
  TYPE(pj_object),VALUE :: dstdefn
  INTEGER(kind=c_int) :: pj_compare_datums
  END FUNCTION pj_compare_datums
END INTERFACE

INTERFACE
  FUNCTION pj_is_latlong(proj) BIND(C,name='pj_is_latlong')
  IMPORT
  TYPE(pj_object),VALUE :: proj
  INTEGER(kind=c_int) :: pj_is_latlong
  END FUNCTION pj_is_latlong
END INTERFACE

INTERFACE
  FUNCTION pj_is_geocent(proj) BIND(C,name='pj_is_geocent')
  IMPORT
  TYPE(pj_object),VALUE :: proj
  INTEGER(kind=c_int) :: pj_is_geocent
  END FUNCTION pj_is_geocent
END INTERFACE

INTERFACE
  FUNCTION pj_get_def(proj, options) BIND(C,name='pj_get_def')
  IMPORT
  TYPE(pj_object),VALUE :: proj
  INTEGER(kind=c_int) :: options
  TYPE(c_ptr) :: pj_get_def
  END FUNCTION pj_get_def
END INTERFACE
  

INTERFACE
  FUNCTION pj_latlong_from_proj(proj) BIND(C,name='pj_latlong_from_proj')
  IMPORT
  TYPE(pj_object),VALUE :: proj
  TYPE(pj_object) :: pj_latlong_from_proj
  END FUNCTION pj_latlong_from_proj
END INTERFACE

!INTERFACE
!  FUNCTION pj_apply_gridshift(c, i, point_count, point_offset, x, y, z) &
!   BIND(C,name='pj_apply_gridshift')
!  IMPORT
!  CHARACTER(kind=c_char) :: c(*)
!  INTEGER(kind=c_int),VALUE :: i
!  INTEGER(kind=c_long),VALUE :: point_count
!  INTEGER(kind=c_int),VALUE :: point_offset
!  REAL(kind=c_double) :: x(*)
!  REAL(kind=c_double) :: y(*)
!  REAL(kind=c_double) :: z(*)
!  INTEGER(kind=c_int) :: pj_apply_gridshift
!  END FUNCTION pj_apply_gridshift
!END INTERFACE
!
!INTERFACE
!  SUBROUTINE pj_deallocate_grids() BIND(C,name='pj_deallocate_grids')
!  END SUBROUTINE pj_deallocate_grids
!END INTERFACE

INTERFACE
  SUBROUTINE pj_free(proj) BIND(C,name='pj_free')
  IMPORT
  TYPE(pj_object),VALUE :: proj
  END SUBROUTINE pj_free
END INTERFACE

!void pj_pr_list(projPJ);
!void pj_set_finder( const char *(*)(const char *) );
!void pj_set_searchpath ( int count, const char **path );
!projPJ pj_init(int, char **);
!char *pj_get_def(projPJ, int);
!void *pj_malloc(size_t);
!void pj_dalloc(void *);
!char *pj_strerrno(int);
!int *pj_get_errno_ref(void);
!const char *pj_get_release(void);

CONTAINS

!> Test whether the result of a pj_init_plus is a valid projection.
!! Returns a logical value. If the second argument is provided, it
!! checks whether they point to the same projection.
FUNCTION pj_associated_object(arg1, arg2) RESULT(associated_)
TYPE(pj_object),INTENT(in) :: arg1 !< projecton object to test
TYPE(pj_object),INTENT(in),OPTIONAL :: arg2 !< optional second argument for equality test instead of validity
LOGICAL :: associated_
IF(PRESENT(arg2)) THEN
  associated_ = C_ASSOCIATED(arg1%ptr, arg2%ptr)
ELSE
  associated_ = C_ASSOCIATED(arg1%ptr)
ENDIF
END FUNCTION pj_associated_object

! Fortran specific version of some functions

!> Fortran version of \a pj_transform proj API function.
!! This is the Fortran version of \a pj_transform function, the array
!! arguments are assumed-shape Fortran arrays of equal length so no
!! array size nor offset need to be passed, see the original C API
!! documentation for the use of the function.
FUNCTION pj_transform_f(src, dst, x, y, z)
TYPE(pj_object),VALUE :: src !< source coordinate system
TYPE(pj_object),VALUE :: dst !< destination coordinate system
REAL(kind=c_double), TARGET           :: x(:) !< array of x coordinates
REAL(kind=c_double), TARGET           :: y(:) !< array of y coordinates
REAL(kind=c_double), TARGET, OPTIONAL :: z(:) !< optional array of z coordinates
INTEGER(kind=c_int) :: pj_transform_f

REAL(kind=c_double),POINTER :: px, py, pz

! a fortran pointer is required to avoid compilation errors with some
! versions of gfortran due to x,y,z being C-incompatible assumed-shape
! arrays
px => x(1)
py => y(1)

IF (PRESENT(z)) THEN
  pz => z(1)
  pj_transform_f = pj_transform(src, dst, &
   INT(MIN(SIZE(x),SIZE(y),SIZE(z)), kind=c_long), 1_c_int, C_LOC(px), C_LOC(py), C_LOC(pz))
ELSE
  pj_transform_f = pj_transform(src, dst, &
   INT(MIN(SIZE(x),SIZE(y)), kind=c_long), 1_c_int, C_LOC(px), C_LOC(py), C_NULL_PTR)
ENDIF

END FUNCTION pj_transform_f


!> Fortran version of \a pj_datum_transform proj API function.
!! This is the Fortran version of \a pj_datum_transform function, the array
!! arguments are assumed-shape Fortran arrays of equal length so no
!! array size nor offset need to be passed, see the original C API
!! documentation for the use of the function.
FUNCTION pj_datum_transform_f(src, dst, x, y, z)
TYPE(pj_object),VALUE :: src !< source coordinate system
TYPE(pj_object),VALUE :: dst !< destination coordinate system
REAL(kind=c_double) :: x(:) !< array of x coordinates
REAL(kind=c_double) :: y(:) !< array of y coordinates
REAL(kind=c_double),OPTIONAL :: z(:) !< optional array of z coordinates
INTEGER(kind=c_int) :: pj_datum_transform_f

REAL(kind=c_double),POINTER :: dummyz(:)

IF (PRESENT(z)) THEN
  pj_datum_transform_f = pj_datum_transform(src, dst, &
   INT(MIN(SIZE(x),SIZE(y),SIZE(z)), kind=c_long), 1_c_int, x, y, z)
ELSE
  NULLIFY(dummyz)
  pj_datum_transform_f = pj_datum_transform(src, dst, &
   INT(MIN(SIZE(x),SIZE(y)), kind=c_long), 1_c_int, x, y, dummyz)
ENDIF

END FUNCTION pj_datum_transform_f


!> Fortran version of \a pj_geocentric_to_geodetic proj API function.
!! This is the Fortran version of \a pj_geocentric_to_geodetic function, the array
!! arguments are assumed-shape Fortran arrays of equal length so no
!! array size nor offset need to be passed, see the original C API
!! documentation for the use of the function.
FUNCTION pj_geocentric_to_geodetic_f(a, es, x, y, z)
REAL(kind=c_double),VALUE :: a !< Earth semi-major axis
REAL(kind=c_double),VALUE :: es !< Earth flattening
REAL(kind=c_double) :: x(:) !< array of x coordinates
REAL(kind=c_double) :: y(:) !< array of y coordinates
REAL(kind=c_double) :: z(:) !< array of z coordinates
INTEGER(kind=c_int) :: pj_geocentric_to_geodetic_f

pj_geocentric_to_geodetic_f = &
 pj_geocentric_to_geodetic(a, es, &
 INT(MIN(SIZE(x),SIZE(y),SIZE(x)), kind=c_long), 1_c_int, x, y, z)

END FUNCTION pj_geocentric_to_geodetic_f


!> Fortran version of \a pj_geodetic_to_geocentric proj API function.
!! This is the Fortran version of \a pj_geodetic_to_geocentric function, the array
!! arguments are assumed-shape Fortran arrays of equal length so no
!! array size nor offset need to be passed, see the original C API
!! documentation for the use of the function.
FUNCTION pj_geodetic_to_geocentric_f(a, es, x, y, z)
REAL(kind=c_double),VALUE :: a !< Earth semi-major axis
REAL(kind=c_double),VALUE :: es !< Earth flattening
REAL(kind=c_double) :: x(:) !< array of x coordinates
REAL(kind=c_double) :: y(:) !< array of y coordinates
REAL(kind=c_double) :: z(:) !< array of z coordinates
INTEGER(kind=c_int) :: pj_geodetic_to_geocentric_f

pj_geodetic_to_geocentric_f = &
 pj_geodetic_to_geocentric(a, es, &
 INT(MIN(SIZE(x),SIZE(y),SIZE(x)), kind=c_long), 1_c_int, x, y, z)

END FUNCTION pj_geodetic_to_geocentric_f


END MODULE proj

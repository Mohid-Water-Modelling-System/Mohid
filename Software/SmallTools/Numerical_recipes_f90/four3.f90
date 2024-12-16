	SUBROUTINE four3(data,isign)
	USE nrtype
	USE nr, ONLY : fourrow_3d
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:,:,:), INTENT(INOUT) :: data
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: dat2,dat3
	call fourrow_3d(data,isign)
	allocate(dat2(size(data,2),size(data,3),size(data,1)))
	dat2=reshape(data,shape=shape(dat2),order=(/3,1,2/))
	call fourrow_3d(dat2,isign)
	allocate(dat3(size(data,3),size(data,1),size(data,2)))
	dat3=reshape(dat2,shape=shape(dat3),order=(/3,1,2/))
	deallocate(dat2)
	call fourrow_3d(dat3,isign)
	data=reshape(dat3,shape=shape(data),order=(/3,1,2/))
	deallocate(dat3)
	END SUBROUTINE four3

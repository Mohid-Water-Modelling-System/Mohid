# Fortran

**Fortran95** is the language of choice to program in MOHID. Although a primitive form of language it proves to be the most efficient and performant in the area of scientific computing. Still widely used in numerical computation. Here's a [crash course](fortran-crash-course.md). Here's a free [manual](http://www.mod-ambiental.mohid.com/Exemplos/138fortran90.pdf). Check out the [DoLoops performance in Fortran](doloops-performance-in-fortran.md).

## Examples

### Hello World program

```
program helloworld

implicit none

write (*,*) "Hello world!"

end program helloworld
```

### Hello world in a module

```
module modulehelloworld

implicit none

private

public :: showhelloworld

contains

subroutine showhelloworld()

        write(*,*) "Hello world!"

end subroutine showhelloworld

end module modulehelloworld
```

```
program helloworld

use modulehelloworld, only: showhelloworld

implicit none

call showhelloworld

end program helloworld
```

## Subroutines & functions

### Subroutine

```
program subfunc

implicit none

real :: a,b,c,root1,root2
logical :: realroots
write (*,*) "Hello world!"

write (*,10)
read (*,*) a,b,c

call solvit(a,b,c,root1,root2,realroots)
write (*,20) root1,root2

if (realroots) then
        write(*,*) 'Sorry, there are no real roots'
end if

10 format('Enter 3 coefficients')
20 format('The roots are ', 2f12.6)
end program subfunc
```

```
subroutine solvit(a,b,c,root1,root2,realroots)

!Arguments
real ::a,b,c,root1,root2
logical :: realroots

!Locals
real :: test

test = b**2-4*a*c

if(test>=0.0) then
        root1 = (-b + sqrt(test))/(2.0*a)
        root2 = (-b - sqrt(test))/(2.0*a)
        realroots = .true.
else
        realroots = .false.
end if

return
end subroutine solvit
```

### Function

```
program subfunc

implicit none

real :: a,b,c,root1,root2,bigroot
logical :: realroots
write (*,*) "Hello world!"

write (*,10)
read (*,*) a,b,c

!call solvit(a,b,c,root1,root2,realroots)
!write (*,20) root1,root2
write(*,20) bigroot(a,b,c)

if (realroots) then
        write(*,*) 'Sorry, there are no real roots'
end if

10 format('Enter 3 coefficients')
!20 format('The roots are ', 2f12.6)
20 format('The biggest root is ', f12.6)
end program subfunc
```

```
function bigroot(a,b,c)

!Arguments
real :: a,b,c,bigroot

!Local
real :: root1, root2, test

if(test>=0.0) then
        root1 = (-b + sqrt(test))/(2.0*a)
        root2 = (-b - sqrt(test))/(2.0*a)
        if (root1 .gt. root2) then
                bigroot = root1
        else
                bigroot = root2
        end if
else
        bigroot = -9.0e35
end if

return
end function bigroot
```

### Interface procedure

## Type variables & arrays

### Variables

The variable declaration is as such

```
real            :: x, y
real, parameter :: z
integer         :: m, n
real(8)         :: w
character(len=30) :: frase
```

### Arrays

```
real, dimension(3) :: vector
real, dimension(3,2) :: matriz3x2
integer, dimension(:)   :: vector-n
integer, dimension(:,:) :: matriz-mxn

!alocação do espaço em memória
allocate(vector-n(21))
allocate(matriz-mxn(4,5))
...
!Limpar a casa no final do programa
deallocate(vector-n)
deallocate(matriz-mxn)
```

### Types

```
public :: T_Pais
type T_Pais
    character(len=128) :: nome
    character(len=128) :: capital
    integer            :: pop
end type T_Pais
```

## Flow control

### If

```
if (.not. associated(pais)) then
    allocate(pais)
else
    write(*,*) 'The country ', pais%nome, 'is already allocated.'
end if
```

### Do

Loop that increments the counter.

```
do i=1,10
   write(*,*) 'Couting ',i
end do
```

Loop that increments the counter by two.

```
do i=1,21,2
   write(*,*) 'Counting pair ',i-1
end do
```

Loop that decrements the counter.

```
do i=20,1,-1
   write(*,*) 'Reverse counting of even ',i*2-1
end do
```

### Case

## Write I/O

Here's sample that writes formatted output to the console.

```
program WriteIOTest

implicit none

character(len=30)               :: myformat
character(len=128)              :: ModuleName
integer                         :: iter
real                            :: DT
integer                         :: i, j, k

ModuleName = "Hello world"
iter = 25
DT = 13.2
i = 113
j = 77
k =43

write (*, 211) trim(ModuleName), iter, DT, i, j, k
pause

211 format(a, i8, f10.3, i5, i5, i5)

end program
```

## Samples

### Average program

```
program average

  ! Read in some numbers and take the average
  ! As written, if there are no data points, an average of zero is returned
  ! While this may not be desired behavior, it keeps this example simple

  implicit none

  real, dimension(:), allocatable :: points
  integer                         :: number_of_points
  real                            :: average_points=0., positive_average=0., negative_average=0.

  write (*,*) "Input number of points to average:"
  read  (*,*) number_of_points

  allocate (points(number_of_points))

  write (*,*) "Enter the points to average:"
  read  (*,*) points

  ! Take the average by summing points and dividing by number_of_points
  if (number_of_points > 0) average_points = sum(points) / number_of_points

  ! Now form average over positive and negative points only
  if (count(points > 0.) > 0) then
     positive_average = sum(points, points > 0.) / count(points > 0.)
  end if
  if (count(points < 0.) > 0) then
     negative_average = sum(points, points < 0.) / count(points < 0.)
  end if

  deallocate (points)

  ! Print result to terminal
  write (*,'(a,g12.4)') 'Average = ', average_points
  write (*,'(a,g12.4)') 'Average of positive points = ', positive_average
  write (*,'(a,g12.4)') 'Average of negative points = ', negative_average

end program average
```

### Countries program (simple)

```
    program Console8

    use modulePais

    implicit none

    type(T_Pais), pointer :: portugal, espanha, france

    integer :: i

    call lerPais(portugal)
    call lerPais(espanha)
    call lerPais(france)

    write(*,*) ''

    call escrevePais(portugal)
    call escrevePais(espanha)
    call escrevePais(france)

    !call ordenaPaisesPorGolos(portugal, espanha, france)
    !Já não dá porque fizemos deallocate dentro da sub escrevePais

    pause

    end program Console8
```

```
 module modulePais

 implicit none

 !Todas as variáveis, os tipos, as funções e as subrotinas
 !ficam privadas por defeito.
 private

 !Listagem das funções e subrotinas públicas e privadas
 public :: lerPais
 public :: escrevePais

 !Declaração de tipos, variáveis e funções/subrotinas publicas e privadas
 !Um novo tipo de variável estruturada
 !public :: StringLength
 integer, parameter  :: StringLength = 128

 public :: T_Pais
 type T_Pais
    character(len=StringLength) :: Nome
    character(len=StringLength) :: Capital
    integer            :: Pop
 end type T_Pais

 !Daqui para a frente só subrotinas e funções
 contains

 subroutine lerPais(pais)

    type(T_Pais), pointer :: pais

    allocate(pais)

    !subroutine lerPais
    write(*,*) 'Por favor insira o nome do País'
    read(*,*) pais%Nome

    write(*,*) 'Por favor insira a capital'
    read(*,*) pais%Capital

    write(*,*) 'Por favor insira a população'
    read(*,*) pais%pop

 end subroutine lerPais

 subroutine escrevePais(um_pais_qq)

    type(T_Pais), pointer :: um_pais_qq

    write(*,*) 'O pais de nome ', trim(um_pais_qq%Nome), ','
    write(*,*) 'cuja capital eh ', trim(um_pais_qq%Capital), ','
    write(*,*) 'tem uma população de ', um_pais_qq%Pop, ' habitantes.'

    deallocate(um_pais_qq)

 end subroutine escrevePais

 end module modulePais
```

### Countries program (linked lists)

```
        program progtypes

        use modulePais

        implicit none

        type(T_Pais), pointer :: portugal

        !1a aula de tipos
        call lerPais(portugal)
        call escrevePais(portugal)
        write(*,*) ''

        !2a aula de tipos
        call novoPais
        write(*,*) ''

        call novoPais
        write(*,*) ''

        call listarPaises

        pause

        call removeAtlas

        end program progtypes
```

```
        module modulePais

        implicit none

        private

        !Declaração do alcance das subrotinas
        !Rotinas que mexem com o tipo "T_Pais"
        public :: lerPais
        public :: escrevePais

        !Rotinas que mexem com o tipo "T_Atlas"
        public :: novoPais
        public :: listarPaises
        private :: novoAtlas
        public :: removeAtlas

        !Declaração de tipos e variáveis globais do módulo
        public :: StringLength
        integer, parameter  :: StringLength = 128

        public :: T_Pais
        type T_Pais
                character(len=StringLength) :: nome
                character(len=StringLength) :: capital
                integer                     :: pop
                real                        :: pib
        end type T_Pais

        public :: T_Atlas
        type T_Atlas
                type(T_Pais), pointer   :: pais
                type(T_Atlas), pointer   :: prim
                type(T_Atlas), pointer   :: prox
        end type T_Atlas

        type(T_Atlas), pointer          :: Me

        contains

        subroutine lerPais(pais)

                type(T_Pais), pointer :: pais

                allocate(pais)

                write(*,*) 'Insira o nome do pais'
                read(*,*) pais%nome

                write(*,*) 'Insira a sua capital'
                read(*,*) pais%capital

                write(*,*) 'Insira a população'
                read(*,*) pais%pop

                write(*,*) 'Insira o PIB em milhões de euros'
                read(*,*) pais%pib

        end subroutine lerPais

        subroutine escrevePais(pais)

                type(T_Pais), pointer :: pais

                write(*,*) 'O pais ', trim(pais%nome),','
                write(*,*) 'cuja capital eh ', trim(pais%capital),','
                write(*,*) 'tem uma populacao de', pais%pop,'habitantes'
                write(*,*) 'e um pib de ', pais%pib, 'milhoes de euros.'

        end subroutine escrevePais

        subroutine novoPais(neoPais)

                !Variável pais
                type(T_Pais), pointer, optional   :: neoPais
                type(T_Pais), pointer		  :: novPais
                type(T_Atlas), pointer            :: novAtlas

                !Um novo país é criado
                if (.not. present(neoPais)) then
                        call lerPais(novPais)
		else
			novPais => neoPais
                end if

                !Se o Atlas ainda não existe, então é criado
                if (.not. associated(Me)) then
                        call novoAtlas(Me, novPais, Me)
                !Senão o Atlas é incrementado
                else
                        call novoAtlas(novAtlas, novPais, Me%prim)
                        Me%prox => novAtlas
                        Me => Me%prox
                end if

        end subroutine novoPais

        subroutine novoAtlas(atlas, pais, prim)

                type(T_Atlas), pointer :: atlas, prim
                type(T_Pais), pointer :: pais

                allocate(atlas)
                atlas%pais => pais
                atlas%prim => prim
                nullify(atlas%prox)

        end subroutine novoAtlas

        subroutine removeAtlas

                type(T_Atlas), pointer :: oldAtlas

                !Rewind the Me global variable
                Me => Me%prim

                !Loop through the Atlas items
                do while(associated(Me))
                        oldAtlas => Me
                        Me => Me%prox
                        deallocate(oldAtlas%pais)
                        deallocate(oldAtlas)
                end do

        end subroutine removeAtlas

        subroutine listarPaises

                type(T_Atlas), pointer :: atlas

                atlas => Me%prim

                do while(associated(atlas))
                        write(*,*) atlas%pais%nome
                        atlas => atlas%prox
                end do

        end subroutine listarPaises

        end module modulePais
```

## MOHID templates

### Module

### Program

## Other languages

- [C/C++](../../archive/linux/cpp.md)
- [Perl](../../archive/linux/perl.md)
- [Python](../../archive/linux/python.md)
- [Ruby](../../archive/linux/ruby.md)

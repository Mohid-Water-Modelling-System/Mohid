# Ping

!------------------------------------------------------------------------------

! IST/MARETEC, Water Modelling Group, Mohid modelling system

!------------------------------------------------------------------------------

!

! TITLE  : ping

! PROJECT  : Example of the actor model using Fortran and MPI

! URL  : [http://www.mohid.com](http://www.mohid.com)

! AFFILIATION  : IST/MARETEC, Marine Modelling Group

! DATE  : Nov 2014

! REVISION  : Rciardo Miranda - v1.0

!

!------------------------------------------------------------------------------

!

!This program is free software; you can redistribute it and/or

!modify it under the terms of the GNU General Public License

!version 2, as published by the Free Software Foundation.

!

!This program is distributed in the hope that it will be useful,

!but WITHOUT ANY WARRANTY; without even the implied warranty of

!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the

!GNU General Public License for more details.

!

!You should have received a copy of the GNU General Public License

!along with this program; if not, write to the Free Software

!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

!

!------------------------------------------------------------------------------

!

! Compile:

! mpif90 -o /myDirectory/projects/lang/fortran/MPI/pingPong/ping /myDirectory/projects/lang/fortran/MPI/pingPong/moduleMPImanagement.f90 /myDirectory/projects/lang/fortran/MPI/pingPong/ping.f90

!

! Run:

! /opt/mpich/bin/mpiexec -n 1 /myDirectory/projects/lang/fortran/MPI/pingPong/pong : -n /myDirectory/projects/lang/fortran/MPI/pingPong/ping

program ping

```
   use mpi
   use moduleMPImanagemnt
```

```
   implicit none
```

```
   integer, parameter :: UNDEFINED   =-99
   integer, parameter :: SUCCESS     = 0
   integer, parameter :: NULL_INT    =-999999
```

```
   type T_pingPong
       integer :: myMPI_id      = NULL_INT
       logical :: hasBall       =.FALSE.
       logical :: gameON        =.TRUE.
       integer :: counting      = NULL_INT
       integer :: otherCounting = NULL_INT
```

```
       integer :: otherMPI_id   = NULL_INT
   end type T_pingPong
```

```
   call main()
```

```
   contains
```

```
   !---------------------------------------------------------------------------
```

```
   subroutine main()
       type(T_pingPong), pointer :: pingPong
       integer :: STAT_CALL
```

```
       STAT_CALL = UNDEFINED
```

```
       pingPong => constructPingPong()
```

```
       STAT_CALL = startGame(pingPong)
```

if1 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "subroutine main, error calling startGame, ERR01"
       end if if1
```

```
       call loop(pingPong)
```

```
       call EXIT(SUCCESS)
```

```
   end subroutine main
```

```
   !---------------------------------------------------------------------------
```

```
   function constructPingPong()
       type(T_pingPong), pointer   :: constructPingPong
       type(T_pingPong), pointer   :: newPingPong
```

```
       allocate(newPingPong)
       newPingPong%myMPI_id = startMPI(newPingPong)
```

```
       constructPingPong => newPingPong
   end function constructPingPong
```

```
   !---------------------------------------------------------------------------
```

```
   integer function startGame(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
       integer                     :: counting
       integer,dimension(8)        :: values
       integer                     :: sec, rand, aux, I
```

```
       STAT_CALL = UNDEFINED
```

```
       call DATE_AND_TIME(VALUES=values)
       sec = values(7)
       call SRAND(sec)
```

```
       aux = rand(0)
```

do1 : do I = -10, aux

```
           counting  = MAX(5, rand(0) / 100000)
       end do do1
```

```
       STAT_CALL = setCounting(pingPong, counting = counting)
```

if2 : if (STAT\_CALL .NE. SUCCESS) then

```
               print*, "STAT_CALL = ", STAT_CALL
               stop "function startGame, error calling setCounting, ERR01"
           end if if2
```

```
       STAT_CALL = sendBall1(pingPong)
```

if3 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function startGame, error calling sendBall1, ERR02"
       end if if3
```

```
       startGame = SUCCESS
   end function startGame
```

```
   !---------------------------------------------------------------------------
```

```
   integer function startMPI(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL, myMPI_id, numprocs
```

```
       STAT_CALL = UNDEFINED
       myMPI_id  = NULL_INT
       numprocs  = NULL_INT
```

```
       call MPI_INIT(IERROR = STAT_CALL)
```

if1 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function startMPI, error calling MPI_INIT, ERR01"
       end if if1
```

```
       call MPI_COMM_RANK(MPI_COMM_WORLD,                             &
                          myMPI_id,                                   &
                          IERROR = STAT_CALL)
```

if2 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function startMPI, error calling MPI_COMM_RANK, ERR02.1"
       end if if2
```

```
       STAT_CALL = setMyMPI_id(pingPong, myMPI_id)
```

if21 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function startMPI, error calling setMyMPI_id, ERR02.2"
       end if if21
```

```
       call MPI_COMM_SIZE(MPI_COMM_WORLD,                               &
                          numprocs,                                     &
                          IERROR = STAT_CALL)
```

if3 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function startMPI, error calling MPI_COMM_SIZE, ERR03"
       end if if3
```

```
       print *, 'Process ', myMPI_id, ' of ', numprocs, ' is alive'
```

if5 : if (getMyMPI\_id(pingPong) .EQ. 0) then

```
           STAT_CALL = setOtherMPI_id(pingPong, otherMPI_id = 1)
```

if51 : if (STAT\_CALL .NE. SUCCESS) then

```
               print*, "STAT_CALL = ", STAT_CALL
               stop "function startMPI, error calling setOtherMPI_id, ERR05"
           end if if51
       else if (getMyMPI_id(pingPong) .EQ. 1) then if5
           STAT_CALL = setOtherMPI_id(pingPong, otherMPI_id = 0)
```

if52 : if (STAT\_CALL .NE. SUCCESS) then

```
               print*, "STAT_CALL = ", STAT_CALL
               stop "function startMPI, error calling setOtherMPI_id, ERR06"
           end if if52
       else if5
           stop "function startMPI, ERR07"
       end if if5
```

```
       startMPI = myMPI_id
   end function startMPI
```

```
   !---------------------------------------------------------------------------
```

```
   integer function getCounting(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       getCounting = pingPong%counting
   end function getCounting
   !---------------------------------------------------------------------------
```

```
   integer function getOtherCounting(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       getOtherCounting = pingPong%otherCounting
   end function getOtherCounting
```

```
   !---------------------------------------------------------------------------
```

```
   integer function getMyMPI_id(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       getMyMPI_id = pingPong%myMPI_id
   end function getMyMPI_id
```

```
   !---------------------------------------------------------------------------
```

```
   integer function getOtherMPI_id(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       getOtherMPI_id = pingPong%otherMPI_id
   end function getOtherMPI_id
```

```
   !---------------------------------------------------------------------------
```

```
   integer function setOtherCounting(pingPong, otherCounting)
       type(T_pingPong), pointer   :: pingPong
       integer, intent(IN)         :: otherCounting
```

```
       pingPong%otherCounting = otherCounting
```

```
       setOtherCounting = SUCCESS
   end function setOtherCounting
   !---------------------------------------------------------------------------
```

```
   integer function setCounting(pingPong, counting)
       type(T_pingPong), pointer   :: pingPong
       integer, intent(IN)         :: counting
```

```
       pingPong%counting = counting
```

```
       setCounting = SUCCESS
   end function setCounting
```

```
   !---------------------------------------------------------------------------
```

```
   integer function decreaseCounting(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       pingPong%counting = getCounting(pingPong)-1
```

```
       decreaseCounting = SUCCESS
   end function decreaseCounting
```

```
   !---------------------------------------------------------------------------
```

```
   integer function setMyMPI_id(pingPong, myMPI_id)
       type(T_pingPong), pointer   :: pingPong
       integer, intent(IN)         :: myMPI_id
```

```
       pingPong%myMPI_id = myMPI_id
```

```
       setMyMPI_id = SUCCESS
   end function setMyMPI_id
```

```
   !---------------------------------------------------------------------------
```

```
   integer function setOtherMPI_id(pingPong, otherMPI_id)
       type(T_pingPong), pointer   :: pingPong
       integer, intent(IN)         :: otherMPI_id
```

```
       pingPong%otherMPI_id = otherMPI_id
```

```
       setOtherMPI_id = SUCCESS
   end function setOtherMPI_id
```

```
   !---------------------------------------------------------------------------
```

```
   integer function setHasBall(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       pingPong%hasBall =.TRUE.
```

```
       setHasBall = SUCCESS
   end function setHasBall
```

```
   !---------------------------------------------------------------------------
```

```
   integer function setHasNoBall(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       pingPong%hasBall =.FALSE.
```

```
       setHasNoBall = SUCCESS
   end function setHasNoBall
```

```
   !---------------------------------------------------------------------------
```

```
   integer function setGameOFF(pingPong)
       type(T_pingPong), pointer   :: pingPong
```

```
       pingPong%gameON =.FALSE.
```

```
       setGameOFF = SUCCESS
   end function setGameOFF
```

```
   !---------------------------------------------------------------------------
```

```
   recursive subroutine loop(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
       integer                     :: STATUS(MPI_STATUS_SIZE)
```

```
       STAT_CALL = UNDEFINED
```

```
       call MPI_PROBE(MPI_ANY_SOURCE,                                          &
                      MPI_ANY_TAG,                                             &
                      MPI_COMM_WORLD,                                          &
                      STATUS,                                                  &
                      STAT_CALL)
```

if25 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "subroutine loop, error calling MPI_PROBE, ERR01"
       end if if25
```

if41 : if (STATUS(MPI\_SOURCE) .EQ. getOtherMPI\_id(pingPong)) then if42 : if (STATUS(MPI\_TAG) .EQ. getMsgPlayBall3Tag()) then

```
               STAT_CALL = receiveBall(pingPong, STATUS(MPI_TAG))
```

if4 : if (STAT\_CALL .NE. SUCCESS) then

```
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "subroutine loop, error calling receiveBall, ERR02"
               end if if4
```

```
           else if (STATUS(MPI_TAG) .EQ. getMsgEndGameTag()) then if42
               STAT_CALL = killPingPong(pingPong)
```

if2 : if (STAT\_CALL .NE. SUCCESS) then

```
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "subroutine main, error calling killPingPong, ERR02"
               end if if2
           end if if42
       end if if41
```

```
       if (pingPong%gameON) call loop(pingPong)
   end subroutine loop
```

```
   !---------------------------------------------------------------------------
```

```
   integer function receiveBall(pingPong, msgTag)
       type(T_pingPong), pointer   :: pingPong
       integer, intent(IN)         :: msgTag
       integer                     :: STATUS(MPI_STATUS_SIZE)
       integer                     :: counting
       integer                     :: STAT_CALL
       integer                     :: i1, i2
```

```
       STAT_CALL  = UNDEFINED
       counting   = NULL_INT
```

```
       call MPI_RECV(counting,                                                 &
                     1,                                                        &
                     MPI_INTEGER,                                              &
                     getOtherMPI_id(pingPong),                                 &
                     getMsgPlayBall3Tag(),                                     &
                     MPI_COMM_WORLD,                                           &
                     STATUS,                                                   &
                     STAT_CALL)
```

if5 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function receiveBall, error calling MPI_RECV, ERR01b"
       end if if5
```

```
       STAT_CALL = setOtherCounting(pingPong, otherCounting = counting)
```

if6 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function receiveBall, error calling setOtherCounting, ERR04"
       end if if6
```

```
       print*, "ping, receiveBall getMyMPI_id =", getMyMPI_id(pingPong), &
               "getCounting =", getCounting(pingPong),                   &
               "getOtherCounting =", getOtherCounting(pingPong)
```

```
       call doSomething()
```

if71 : if (getCounting(pingPong) .GT. 0) then if70 : if ( mod(getCounting(pingPong), 2) .EQ. 0) then

```
               STAT_CALL = sendBall1(pingPong)
```

if7 : if (STAT\_CALL .NE. SUCCESS) then

```
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "function receiveBall, error calling sendBall1, ERR01"
               end if if7
```

```
           else if ((mod(getCounting(pingPong), 3) .EQ. 0) .AND.         &
                    (mod(getCounting(pingPong), 2) .NE. 0)) then if70
               STAT_CALL = sendBall2(pingPong)
```

if75 : if (STAT\_CALL .NE. SUCCESS) then

```
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "function receiveBall, error calling sendBall2, ERR01"
               end if if75
```

```
           else if70
               STAT_CALL = sendBall3(pingPong)
```

if76 : if (STAT\_CALL .NE. SUCCESS) then

```
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "function receiveBall, error calling sendBall3, ERR01"
               end if if76
           end if if70
       else if71
           ! Program terminates because there are no more plays
           STAT_CALL = sendKillGame(pingPong)
```

if3 : if (STAT\_CALL .NE. SUCCESS) then

```
               print*, "STAT_CALL = ", STAT_CALL
               stop "subroutine main, error calling killPingPong, ERR02"
           end if if3
```

```
           STAT_CALL = killPingPong(pingPong)
```

if2 : if (STAT\_CALL .NE. SUCCESS) then

```
               print*, "STAT_CALL = ", STAT_CALL
               stop "subroutine main, error calling killPingPong, ERR03"
           end if if2
       end if if71
```

```
       receiveBall = SUCCESS
   end function receiveBall
```

```
   !---------------------------------------------------------------------------
```

```
   subroutine doSomething
       integer                     :: i,j,k, res1, res2, res3
```

```
       call sleep(10)
```

do1 : do i=1,1000 do2 : do j=1,1000 do3 : do k=1,1000

```
           res1=i*j+j/(k**2)*(-1.0)
           res2=i*j+j/(k**2)*(-21.0)
           res3=MAX(res1, res2)
       end do do3
       end do do2
       end do do1
```

```
   end subroutine doSomething
```

```
   !---------------------------------------------------------------------------
```

```
   integer function sendBall1(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
```

```
       STAT_CALL  = UNDEFINED
```

```
       STAT_CALL = decreaseCounting(pingPong)
```

if1 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendBall1, error calling decreaseCounting, ERR01"
       end if if1
```

```
       call MPI_SEND(getCounting(pingPong),                                    &
                     1,                                                        &
                     MPI_INTEGER,                                              &
                     getOtherMPI_id(pingPong),                                 &
                     getMsgPlayBall1Tag(),                                     &
                     MPI_COMM_WORLD,                                           &
                     STAT_CALL)
```

if16 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendBall1, error calling MPI_SEND, ERR012"
       end if if16
```

```
       sendBall1 = SUCCESS
   end function sendBall1
```

```
   !---------------------------------------------------------------------------
```

```
   integer function sendBall2(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
```

```
       STAT_CALL  = UNDEFINED
```

```
       STAT_CALL = decreaseCounting(pingPong)
```

if1 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendBall2, error calling decreaseCounting, ERR01"
       end if if1
```

```
       call MPI_SEND(getCounting(pingPong),                                    &
                     1,                                                        &
                     MPI_INTEGER,                                              &
                     getOtherMPI_id(pingPong),                                 &
                     getMsgPlayBall2Tag(),                                     &
                     MPI_COMM_WORLD,                                           &
                     STAT_CALL)
```

if16 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendBall2, error calling MPI_SEND, ERR012"
       end if if16
```

```
       sendBall2 = SUCCESS
   end function sendBall2
```

```
   !---------------------------------------------------------------------------
```

```
   integer function sendBall3(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
```

```
       STAT_CALL  = UNDEFINED
```

```
       STAT_CALL = decreaseCounting(pingPong)
```

if1 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendBall3, error calling decreaseCounting, ERR01"
       end if if1
```

```
       call MPI_SEND(getCounting(pingPong),                                    &
                     1,                                                        &
                     MPI_INTEGER,                                              &
                     getOtherMPI_id(pingPong),                                 &
                     getMsgPlayBall3Tag(),                                     &
                     MPI_COMM_WORLD,                                           &
                     STAT_CALL)
```

if16 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendBall3, error calling MPI_SEND, ERR012"
       end if if16
```

```
       sendBall3 = SUCCESS
   end function sendBall3
```

```
   !---------------------------------------------------------------------------
```

```
   integer function sendKillGame(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
```

```
       STAT_CALL  = UNDEFINED
```

```
       call MPI_SEND(NULL_INT,                                                 &
                     1,                                                        &
                     MPI_INTEGER,                                              &
                     getOtherMPI_id(pingPong),                                 &
                     getMsgEndGameTag(),                                       &
                     MPI_COMM_WORLD,                                           &
                     STAT_CALL)
```

if16 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function sendKillGame, error calling MPI_SEND, ERR01"
       end if if16
```

```
       sendKillGame = SUCCESS
   end function sendKillGame
```

```
   !---------------------------------------------------------------------------
```

```
   integer function killPingPong(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
```

```
       STAT_CALL = UNDEFINED
```

```
       STAT_CALL = stopMPI()
```

if2 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function killPingPong, error calling stopMPI, ERR01"
       end if if2
```

```
       print*, "Ping Pong Terminated, ping"
```

```
       STAT_CALL = setGameOFF(pingPong)
```

if3 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function killPingPong, error calling setGameOFF, ERR01"
       end if if3
```

```
       killPingPong = SUCCESS
   end function killPingPong
   !---------------------------------------------------------------------------
```

```
   integer function stopMPI()
       integer                     :: STAT_CALL
```

```
       STAT_CALL = UNDEFINED
```

```
       call MPI_BARRIER(MPI_COMM_WORLD,                               &
                        IERROR = STAT_CALL)
```

if6 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function stopMPI, error calling MPI_BARRIER, ERR01"
       end if if6
```

```
       call MPI_FINALIZE(STAT_CALL)
```

if2 : if (STAT\_CALL .NE. SUCCESS) then

```
           print*, "STAT_CALL = ", STAT_CALL
           stop "function stopMPI, error calling MPI_FINALIZE, ERR02"
       end if if2
```

```
       stopMPI =  SUCCESS
   end function stopMPI
```

end program ping

# Fortran implementation of the Actor Model using MPI

A very simple ping-pong program is shown here. There are 2 programs, *[ping](../../reference/ping.md)* and *[pong](../../reference/pong.md)*, that shoot messages among them, starting each with a random number of plays. The first program reaching 0 will send a message for the other program to stop and stops itself. File *[moduleMPImanagement](../modulempimanagement.md)* is necessary to compile both programs.

If you examine the code there is no single MPI\_BARRIER or any other explicit synchronization point. The main subroutine is:

```
   subroutine main()
       type(T_pingPong), pointer :: pingPong
       integer :: STAT_CALL
       STAT_CALL = UNDEFINED
       pingPong => constructPingPong()
       STAT_CALL = startGame(pingPong)
       if (STAT_CALL .NE. SUCCESS) then
           print*, "STAT_CALL = ", STAT_CALL
           stop "subroutine main, error calling startGame, ERR01"
       end if
       call loop(pingPong) !This loop reads MPI messages
       call EXIT(SUCCESS)
   end subroutine main
```

There is a call to an infinite loop. This routine reads the MPI message queue identifying messages to itself (MPI\_PROBE).

```
   recursive subroutine loop(pingPong)
       type(T_pingPong), pointer   :: pingPong
       integer                     :: STAT_CALL
       integer                     :: STATUS(MPI_STATUS_SIZE)
       STAT_CALL = UNDEFINED
       call MPI_PROBE(MPI_ANY_SOURCE,                                          &
                      MPI_ANY_TAG,                                             &
                      MPI_COMM_WORLD,                                          &
                      STATUS,                                                  &
                      STAT_CALL)
       if (STAT_CALL .NE. SUCCESS) then
           print*, "STAT_CALL = ", STAT_CALL
           stop "subroutine loop, error calling MPI_PROBE, ERR01"
       end if
       if (STATUS(MPI_SOURCE) .EQ. getOtherMPI_id(pingPong)) then
           if (STATUS(MPI_TAG) .EQ. getMsgPlayBall3Tag()) then
               STAT_CALL = receiveBall(pingPong, STATUS(MPI_TAG))
               if (STAT_CALL .NE. SUCCESS) then
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "subroutine loop, error calling receiveBall, ERR02"
               end if
           else if (STATUS(MPI_TAG) .EQ. getMsgEndGameTag()) then
               STAT_CALL = killPingPong(pingPong)
              if (STAT_CALL .NE. SUCCESS) then
                   print*, "STAT_CALL = ", STAT_CALL
                   stop "subroutine main, error calling killPingPong, ERR02"
               end if
           end if
       end if
       if (pingPong%gameON) call loop(pingPong)
   end subroutine loop
```

In the recursive subroutine loop 2 messages are processed: getMsgPlayBall3Tag and getMsgEndGameTag. As promised by the actor model:

- internal state is encapsulated;
- messages are processed one at a time; and
- communication is asynchronous.

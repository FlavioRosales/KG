module mpi_lib
  use mpi
  implicit none
  !
  ! return an error value
  !
  integer :: ierr
  !
  ! rank of the calling process in the group of comm
  !
  integer :: rank
  !
  ! number of processes in the group of comm
  !
  integer :: nproc
  !
  integer :: master = 0
  !



contains
  !============================================================================!
  subroutine MPI_CREATE
    implicit none

    !
    ! the number of processes in the group associated with the communicator MPI_COMM_WORLD.
    !
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    !
    ! the rank of the calling process in the group associated with the communicator MPI_COMM_WORLD.
    !
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    !
    ! converts process rank to proces grid coordinates.
    !

  end subroutine MPI_CREATE
  !============================================================================!
 

end module mpi_lib

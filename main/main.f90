program main
  use global_numbers
  implicit none

  integer :: level, i, j, k

  !
  ! Initialize MPI interfaz
  !
  call MPI_INIT(ierr)
  !
  ! Read and create initial data and parameters
  !
  call read_data
  call create_data
  !
  ! evolution using berger-oliger algorithm by rk3 method
  !
  call evolve
  !
  ! Close MPI interfaz
  !
  call MPI_FINALIZE(ierr)

end program main

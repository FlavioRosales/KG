!=============================================================================!
! Autor: Iván Álvarez-Rios
! Date: September 10, 2021
!
! Description:
! Evenly spaced concentric mesh refinement in each spatial direction.
!=============================================================================!
module mesh_refinement_lib
  use mpi_lib
  implicit none

  integer, private :: count, status(MPI_STATUS_SIZE), request


  !============================================================================!
  type, public :: mesh_refinement
    integer :: n_refinements, Nx, Ny, Nz, g_pts
    integer :: iL, iR, jL, jR, kL, kR
    integer, allocatable, dimension(:,:) :: index_bounds
    real(kind=8) :: xmin, xmax, dx
    real(kind=8) :: dy
    real(kind=8) :: dz
    real(kind=8), allocatable, dimension(:,:,:) :: x, y, z
    real(kind=8), allocatable, dimension(:,:,:) :: r, theta, phi
  contains
    procedure :: mesh_refinement_memory
    procedure :: mesh_refinement_boundary
    procedure :: mesh_refinement_create
    procedure :: mesh_refinement_info
  end type mesh_refinement
  !============================================================================!
  interface mesh_refinement
    module procedure :: mesh_refinement_constructor
  end interface mesh_refinement
  !============================================================================!
contains
  !============================================================================!
  function mesh_refinement_constructor(xmin, xmax, Nx, Ny, Nz, g_pts) result(this)
    implicit none
    integer, intent(in) :: Nx, Ny, Nz, g_pts
    real(kind=8), intent(in) :: xmin, xmax
    type(mesh_refinement) :: this
    call MPI_CREATE

    if(mod(Nx, nproc).ne.0) stop "Error Nx must be divisible by nproc."

    this%g_pts = g_pts
    this%Nx = Nx / nproc; this%Ny = Ny ; this%Nz = Nz
    this%iL = 0; this%iR = this%Nx
    this%jL = 0; this%jR = this%Ny
    this%kL = 0; this%kR = this%Nz
    this%xmin = xmin
    this%xmax = xmax


    call this%mesh_refinement_memory
    call this%mesh_refinement_boundary(xmin, xmax)
    call this%mesh_refinement_create

    this%index_bounds(1,1) = this%iL - this%g_pts !index x left
    this%index_bounds(1,2) = this%iR + this%g_pts !index x right
    this%index_bounds(2,1) = this%jL - this%g_pts !index y left
    this%index_bounds(2,2) = this%jR + this%g_pts !index y right
    this%index_bounds(3,1) = this%kL !index z left
    this%index_bounds(3,2) = this%kR !index z right

  end function mesh_refinement_constructor
  !============================================================================!
  subroutine mesh_refinement_destructor(this)
    implicit none
    type(mesh_refinement), intent(in out) :: this


  end subroutine mesh_refinement_destructor
  !============================================================================!
  subroutine mesh_refinement_memory(this)
    implicit none
    class(mesh_refinement), intent(in out) :: this
    allocate(this%x(&
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts, &
    this%kL:this%kR))
    allocate(this%y(&
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts, &
    this%kL:this%kR))
    allocate(this%z(&
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts, &
    this%kL:this%kR))
    allocate(this%r(this%Nx,this%Ny,this%Nz))
    allocate(this%theta(this%Nx,this%Ny,this%Nz))
    allocate(this%phi(this%Nx,this%Ny,this%Nz))


    allocate(this%index_bounds(3,2)) 

  end subroutine mesh_refinement_memory
  !============================================================================!
  subroutine mesh_refinement_boundary(this, xmin, xmax)
    implicit none
    class(mesh_refinement), intent(in out) :: this
    real(kind=8), intent(in) :: xmin, xmax;     integer :: i


    real(kind=8) :: dx


    dx = (xmax - xmin) / dble(nproc)


    this%xmin  = xmin + dble(rank + 0) * dx; this%xmax = xmin + dble(rank + 1) * dx


    this%dx = (this%xmax - this%xmin) / dble(this%Nx)
    this%dy =  acos(-1.0d0)/dble(this%Ny-1)
    this%dz = 2.0d0*acos(-1.0d0)/dble(this%Nz -1)


  end subroutine mesh_refinement_boundary
  !============================================================================!
  subroutine mesh_refinement_create(this)
    implicit none
    class(mesh_refinement), intent(in out) :: this
    integer :: i, j, k

    do i=this%iL - this%g_pts, this%iR + this%g_pts
      do j=this%jL - this%g_pts,  this%jR + this%g_pts 
        do k=this%kL, this%kR
          this%x(i,j,k) =  this%xmin + dble(i+this%g_pts)*this%dx
          this%y(i,j,k) = -this%dy/2.0d0 + dble(j) * this%dy
          this%z(i,j,k) =  dble(k) * this%dz
        end do
      end do
    end do

    do i=1, this%Nx
      do j=1, this%Ny 
        do k=1, this%Nz
          this%r(i,j,k) = this%xmin + dble(i-1)*this%dx 
          this%theta(i,j,k) = dble(j-1)*this%dy  
          this%phi(i,j,k) = dble(k-1)*this%dz  
        end do 
      end do 
    end do 

  end subroutine mesh_refinement_create
  !============================================================================!
  subroutine mesh_refinement_info(this)
    implicit none
    class(mesh_refinement), intent(in out) :: this
    character(len=100) :: FMT
    integer :: i

    FMT = "(A, I2, A, F6.2, A, F6.2, A, F6.2, A, F6.2, A, F6.2, A, F6.2, A)"


      if(rank.eq.0) then
        print*, '|======================================================================================|'
        print*, '|  rank  |   r_min  |   r_max  |   theta_min  |   theta_max  |   phi_min  |   phi_max  |'
        print*, '|======================================================================================|'
      end if

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

      write(*, FMT) ' | ', rank, '     | ', &
      this%x(this%iL-this%g_pts,this%jL - this%g_pts ,this%kL),  &
      '  | ',this%x(this%iR+this%g_pts,this%jR + this%g_pts ,this%kR), '  | ', &
      this%y(this%iL-this%g_pts,this%jL - this%g_pts ,this%kL),  &
      '  | ',this%y(this%iR+this%g_pts,this%jR + this%g_pts ,this%kR), '  | ', &
      this%z(this%iL-this%g_pts,this%jL ,this%kL), '  | ',this%z(this%iR+this%g_pts,this%jR ,this%kR), '  |'

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)


    if(rank.eq.0) print*,  '|======================================================================================|'
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine mesh_refinement_info
  !============================================================================!
end module mesh_refinement_lib

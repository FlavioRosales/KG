!=============================================================================!
! Autor: Iván Álvarez-Rios
! Date: September 10, 2021
!
! Description:
! runge-kutta method 4.
!=============================================================================!
module lines_method_lib
  use geometry_lib
  implicit none

  integer, private :: count, status(MPI_STATUS_SIZE), request

  !============================================================================!
  type, extends(geometry) :: lines_method

    integer :: nvars, Nt, n_det
    real(kind=8) :: tmin, tmax, CFL
    real(kind=8) :: m !masa del bosón
    real(kind=8) :: t, dt
    real(kind=8), allocatable, dimension(:) :: r_det, flux
    real(kind=8), allocatable, dimension(:,:,:,:) :: rhs_u
    real(kind=8), allocatable, dimension(:,:,:,:) :: u, u_p
    character(len=20) :: boundary_type_rmax

  contains

    procedure :: lines_method_info
    procedure :: rk3
    procedure :: diagnostic
    procedure, private :: lines_method_memory

    procedure, private :: send_recv_xL, send_recv_xR

    procedure, private :: boundary_xL, boundary_xR
    procedure, private :: boundary_yL, boundary_yR
    procedure, private :: boundary_zL, boundary_zR

    procedure :: check_point, read_from_check_point

  end type lines_method
  !============================================================================!
  interface lines_method
    module procedure :: lines_method_constructor
  end interface lines_method
  !============================================================================!
contains
  !============================================================================!
  function lines_method_constructor(&
     nvars, g_pts, xmin, xmax, Nx, Ny, Nz, tmin, tmax, CFL, &
     boundary_type_rmax, metric, m, n_det, first_det, space_det) result(this)
    implicit none
    integer, intent(in) ::  nvars, g_pts, Nx, Ny, Nz
    integer, intent(in) :: n_det
    real(kind=8), intent(in) :: first_det, space_det
    real(kind=8), intent(in) :: xmin, xmax, tmin, tmax, CFL, m
    character(len=*) :: boundary_type_rmax, metric
    type(lines_method) :: this
    integer i

    call MPI_CREATE
    !
    ! constructor for mesh refinement
    !
    if(mod(Nx, nproc).ne.0) stop "Error Nx must be divisible by nproc."

    this%nvars =  nvars
    this%g_pts = g_pts
    this%Nx = Nx / nproc; this%Ny = Ny ; this%Nz = Nz
    this%iL = 0; this%iR = this%Nx
    this%jL = 0; this%jR = this%Ny
    this%kL = 0; this%kR = this%Nz
    this%m = m

    print*, ' boson mass m_b = ', this%m

    this%xmin = xmin
    this%xmax = xmax

    this%boundary_type_rmax =  boundary_type_rmax

    call this%mesh_refinement_memory()
    call this%mesh_refinement_boundary(xmin, xmax)
    call this%mesh_refinement_create()
    !
    call this%lines_method_memory()
    this%tmin = tmin; this%tmax = tmax; this%CFL = CFL
    !this%dt = CFL * min(this%dx, this%dy, this%dz)
    this%dt = CFL * min(this%dx, (xmin+this%dx)*this%dy, &
      (xmin+this%dx)*sin(this%dy/2.0d0)*this%dz)

    !Max aristas 

    print*, 'dr_max', this%dx 
    print*, 'dtheta_max', this%x(this%iR + this%g_pts,1,1)*this%dy
    print*, 'dphi_max', this%x(this%iR + this%g_pts,1,1)*this%dz


    this%Nt = int((this%tmax - this%tmin) / this%dt)
    if(mod(this%Nt,2).ne.0) this%Nt = this%Nt + 1
    this%t = this%tmin
    call this%mesh_refinement_info()
    call this%lines_method_info()

    ! goemetry
    this%metric = metric
    call this%geometry_memory
    call this%geometry_load_metric
    call this%geometry_load_inverse

    !radius detectors
    this%n_det = n_det

    do i=1, this%n_det
      this%r_det(i) = first_det + dble(i-1)*space_det
    end do


  end function lines_method_constructor
  !============================================================================!
  subroutine lines_method_memory(this)
    implicit none
    class(lines_method), intent(in out) :: this

    allocate(this%u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%u_p(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%rhs_u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%r_det(this%n_det))
    allocate(this%flux(this%n_det))


  end subroutine lines_method_memory
  !============================================================================!
  subroutine lines_method_info(this)
    implicit none
    class(lines_method), intent(in out) :: this
    character(len=100) :: FMT

    FMT = "(A, F6.1, A, F9.1, A, F6.3, A, ES8.2, A, I10, A)"

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if(rank.eq.master) then
      print*, '|====================================================================|'
      print*, '|                    Temporary this information                    |'
      print*, '|====================================================================|'
      print*, '|    tmin    |     tmax    |   CFL   |     dt    |         Nt        |'
      print*, '|====================================================================|'
      write(*,FMT) ' |  ', this%tmin, '    |  ', this%tmax, '  | ', this%CFL, '  | ', this%dt, '  |  ', this%Nt, '       |'
      print*, '|====================================================================|'
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)


  end subroutine lines_method_info
  !============================================================================!
  subroutine rk3(this, rhs)
    implicit none
    class(lines_method), intent(in out) :: this
    !external :: rhs
    interface
      subroutine rhs(refinement)
        implicit none
        integer, intent(in) :: refinement
      end subroutine
    end interface
    integer :: rk

    this%u_p = this%u

    do rk=1, 3

      call rhs(1)

      if(rk.eq.1) then
        this%u = this%u_p + this%rhs_u * this%dt
      else if(rk.eq.2) then
        this%u = 0.750d0 * this%u_p &
        + 0.250d0 * (this%u + this%rhs_u * this%dt)
      else
        this%u = this%u_p/ 3.0d0  &
        + 2.0d0 * (this%u + this%rhs_u * this%dt) / 3.0d0
      end if

      call boundary

    end do

    this%t = this%t + this%dt

  contains

    subroutine boundary

      call this%boundary_zL()
      call this%boundary_zR()
      call this%boundary_yL()
      call this%boundary_yR()

      if(rank /=0) then
        call this%send_recv_xL()
      else
        call this%boundary_xL()
      end if

      if(rank /=nproc-1) then
        call this%send_recv_xR()
      else
        call this%boundary_xR()
      end if


    end subroutine boundary

  end subroutine rk3
  !============================================================================!
  subroutine send_recv_xL(this)
    implicit none
    class(lines_method), intent(in out) :: this


    count = size(this%u(:,this%iL:this%iL+this%g_pts,:,:))

    call MPI_SENDRECV(&
    this%u(:,this%iL:this%iL+this%g_pts,:,:), &
    count, MPI_DOUBLE_PRECISION, rank -1, 1, &  ! send
    this%u(:,this%iL-this%g_pts:this%iL,:,:),&
    count, MPI_DOUBLE_PRECISION, rank -1, 0, &  ! recv
    MPI_COMM_WORLD, status, ierr)

  end subroutine send_recv_xL
  !============================================================================!
  subroutine send_recv_xR(this)
    implicit none
    class(lines_method), intent(in out) :: this


    count = size(this%u(:,this%iR:this%iR+this%g_pts,:,:))

    call MPI_SENDRECV(&
    this%u(:,this%iR-this%g_pts:this%iR,:,:), &
    count, MPI_DOUBLE_PRECISION, rank +1, 0, &  ! send
    this%u(:,this%iR:this%iR+this%g_pts,:,:), &
    count, MPI_DOUBLE_PRECISION, rank +1, 1, &  ! recv
    MPI_COMM_WORLD, status, ierr)

  end subroutine send_recv_xR
  !============================================================================!

  subroutine boundary_xL(this)
    implicit none
    class(lines_method), intent(in out) :: this
    integer :: i

    ! boundary xL == r_min, radius of excision
    ! third order extrapolations

    ! this%u(:,this%iL-this%g_pts,:,:) = this%u(:,this%iL-this%g_pts -1,:,:)

      this%u(:,this%iL-this%g_pts,:,:) = 3.0d0*this%u(:,this%iL-this%g_pts+1,:,:) &
            -3.0d0*this%u(:,this%iL-this%g_pts+2,:,:) + this%u(:,this%iL-this%g_pts+3,:,:)

  end subroutine boundary_xL
  !============================================================================!
  subroutine boundary_xR(this)
    implicit none
    class(lines_method), intent(in out) :: this
    integer i

    if(this%boundary_type_rmax.eq.'outflow') then
      this%u(:,this%iR+this%g_pts,:,:) = this%u(:,this%iR+this%g_pts-1,:,:)
    else if(this%boundary_type_rmax.eq.'inflow') then
      this%u(:,this%iR+this%g_pts,:,:) = this%u_p(:,this%iR+this%g_pts,:,:)
    else if(this%boundary_type_rmax.eq.'wind') then
      do i=this%jL , this%jR 
      if(this%y(this%iR,i,this%kR).lt.2.0d0*atan(1.0d0)) then
      this%u(:,this%iR+this%g_pts,i,:) = this%u(:,this%iR+this%g_pts-1,i,:)
      else
      this%u(:,this%iR+this%g_pts,i,:) = this%u_p(:,this%iR+this%g_pts,i,:)
      end if
      end do
    else if(this%boundary_type_rmax.eq.'wind_y') then
      do i=this%kL, this%kR
      if(this%z(this%iR,this%jR,i).lt.acos(-1.0d0)) then
      this%u(:,this%iR+this%g_pts,:,i) = this%u(:,this%iR+this%g_pts-1,:,i)
      else
      this%u(:,this%iR+this%g_pts,:,i) = this%u_p(:,this%iR+this%g_pts,:,i)
      end if
      end do
    else
      print*, 'This boundary condition type for r_max is not implemented yet'
      print*, 'try with: outflow or inflow'
      stop
    end if

  end subroutine boundary_xR
  !============================================================================!
  subroutine boundary_yL(this)
    implicit none
    class(lines_method), intent(in out) :: this
      integer :: i,k,dis,j, g

      dis = (this%Nz+1)/2

      do k = this%kL, this%kR
        do g=this%g_pts, 0, -1
        if(k.ge.abs(dis)) dis = -abs(dis)
          this%u(1,:,this%jL-g,k) = this%u(1,:,this%jL+g+1,k+dis)
          this%u(2,:,this%jL-g,k) = this%u(2,:,this%jL+g+1,k+dis)
          this%u(3,:,this%jL-g,k) = -this%u(3,:,this%jL+g+1,k+dis)
          this%u(4,:,this%jL-g,k) = this%u(4,:,this%jL+g+1,k+dis)
          this%u(5,:,this%jL-g,k) = this%u(5,:,this%jL+g+1,k+dis)
        end do
      end do




  end subroutine boundary_yL
  !============================================================================!
  subroutine boundary_yR(this)
    implicit none
    class(lines_method), intent(in out) :: this
      integer :: i,k,dis,j,g

      dis = (this%Nz + 1)/2

      do k = this%kL, this%kR
        do g= this%g_pts, 0, -1
        if(k.ge.abs(dis)) dis = -abs(dis)
          this%u(1,:,this%jR + g,k) = this%u(1,:,this%jR-g-1,k+dis)
          this%u(2,:,this%jR + g,k) = this%u(2,:,this%jR-g-1,k+dis)
          this%u(3,:,this%jR + g,k) = -this%u(3,:,this%jR-g-1,k+dis)
          this%u(4,:,this%jR + g,k) = this%u(4,:,this%jR-g-1,k+dis)
          this%u(5,:,this%jR + g,k) = this%u(5,:,this%jR-g-1,k+dis)
        end do
      end do


  end subroutine boundary_yR
  !============================================================================!
  subroutine boundary_zL(this)
    implicit none
    class(lines_method), intent(in out) :: this
      integer :: i

        this%u(:,:,:,this%kL) = this%u(:,:,:,this%kR-1)

  end subroutine boundary_zL
  !============================================================================!
  subroutine boundary_zR(this)
    implicit none
    class(lines_method), intent(in out) :: this
      integer :: i

      this%u(:,:,:,this%kR) = this%u(:,:,:,this%kL+1)


  end subroutine boundary_zR
  !============================================================================!
  subroutine check_point(this)
    use strings_lib
    implicit none
    class(lines_method), intent(in out) :: this
    type(string) :: name
    integer :: i,j,k,refinement,ios

    name = 'check_point_' + str(rank+1) + '.data'

    call system('rm -rf ' // trim(name%string_data))

    open(unit=100, file=name%string_data, iostat=ios, status="new", action="write")
    if ( ios /= 0 ) stop "Error opening file check_point."

      write(unit=100, fmt=*, iostat=ios) &
      this%nvars, this%Nx, this%Ny, this%Nz, this%g_pts
      if ( ios /= 0 ) stop "Write error in file unit check_point."
      write(unit=100, fmt=*, iostat=ios) this%t
      if ( ios /= 0 ) stop "Write error in file unit check_point."

        do i = this%iL - this%g_pts, this%iR + this%g_pts
          do j = this%jL , this%jR 
            do k = this%kL, this%kR
              write(unit=100, fmt=*, iostat=ios) &
              this%x(i,j,k), this%y(i,j,k), this%z(i,j,k), this%u(:,i,j,k)
              if ( ios /= 0 ) stop "Write error in file unit check_point."
            end do
          end do
        end do

    close(unit=100, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file unit check_point."

  end subroutine check_point
  !============================================================================!
  subroutine read_from_check_point(this)
    use strings_lib
    implicit none
    class(lines_method), intent(in out) :: this
    type(string) :: name
    integer :: i,j,k,refinement,ios

    name = 'check_point_' + str(rank+1) + '.data'

    call system('rm -rf ' // trim(name%string_data))

    open(unit=100, file=name%string_data, iostat=ios, status="old", action="read")
    if ( ios /= 0 ) stop "Error opening file check_point."

      read(unit=100, fmt=*, iostat=ios) &
      this%nvars, this%Nx, this%Ny, this%Nz, this%g_pts
      if ( ios /= 0 ) stop "Write error in file unit check_point."
      read(unit=100, fmt=*, iostat=ios) this%t
      if ( ios /= 0 ) stop "Write error in file unit check_point."

        do i = this%iL - this%g_pts, this%iR + this%g_pts
          do j = this%jL , this%jR 
            do k = this%kL, this%kR
              read(unit=100, fmt=*, iostat=ios) &
              this%x(i,j,k), this%y(i,j,k), this%z(i,j,k), this%u(:,i,j,k)
              if ( ios /= 0 ) stop "Write error in file unit check_point."
            end do
          end do
        end do

    close(unit=100, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file unit check_point."

  end subroutine read_from_check_point

  subroutine diagnostic(this)
    use integral
    implicit none
    class(lines_method), intent(in out) :: this
    real(kind=8), dimension( & 
        this%jL - this%g_pts:this%jR + this%g_pts, &
        this%kL:this%kR) :: flux_shell
    integer :: i_det
    integer :: l, i

    do l=1, this%n_det
    if((this%x(this%iL,1,1).le.this%r_det(l)).and.(this%r_det(l).ge.this%x(this%iR,1,1))) then 
        do i=this%iL, this%iR
            if((this%x(i,1,1).le.this%r_det(l)).and.(this%x(i+1,1,1).gt.this%r_det(l))) then
                i_det = i
                print*,rank
            end if
        end do
    end if
    flux_shell = this%u(2,i_det,:,:) * this%alpha(i_det,:,:) * this%sqrtg(i_det,:,:) 
    this%flux(l) = trapezium_2D(flux_shell,this%dy,this%dz)
    end do
  end subroutine diagnostic


end module lines_method_lib

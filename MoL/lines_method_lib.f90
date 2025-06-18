!=============================================================================!
! Autor: Iván Álvarez-Rios
! Date: September 10, 2021
!
! Description:
! runge-kutta method 4.
!=============================================================================!
module lines_method_lib
  use geometry_lib
  use FFT_sph
  implicit none

  integer, parameter :: l_max = 10
  integer, parameter :: Nk = 100

  integer, private :: count, status(MPI_STATUS_SIZE), request

  !============================================================================!
  type, extends(geometry) :: lines_method

    integer :: nvars, Nt, n_det
    real(kind=8) :: tmin, tmax, CFL
    real(kind=8) :: m !masa del bosón
    real(kind=8) :: t, dt
    real(kind=8) :: r_max
    complex(kind=8), allocatable, dimension(:,:,:,:) :: rhs_u
    complex(kind=8), allocatable, dimension(:,:,:,:) :: u, u_p, u_stage, u_stage_2
    complex(kind=8), allocatable, dimension(:,:,:,:) :: k1_u, k2_u, k3_u, k4_u
    character(len=20) :: boundary_type_rmax

  contains

    procedure :: lines_method_info
    procedure :: rk3, rk4, pirk2, pirk3
    procedure :: rhs_phi, rhs_psi1, rhs_psi2, rhs_psi3, rhs_pi
    procedure :: lines_method_memory
    procedure :: send_recv_xL, send_recv_xR 
    procedure :: boundary_xL, boundary_xR
    procedure :: boundary_yL, boundary_yR
    procedure :: boundary_zL, boundary_zR

    procedure :: check_point, read_from_check_point

  end type lines_method

  !============================================================================!
contains
  !============================================================================!
  function lines_method_constructor(&
     nvars, g_pts, xmin, xmax, Nx, Ny, Nz, tmin, tmax, CFL, &
     boundary_type_rmax, metric, a, M_bh, m, save1d,r_max) result(this)
    implicit none
    integer, intent(in) ::  nvars, g_pts, Nx, Ny, Nz, save1d
    real(kind=8), intent(in) :: a, M_bh
    real(kind=8), intent(in) :: xmin, xmax, tmin, tmax, CFL, m
    real(kind=8) :: dx_aux, dxx, r_max 
    integer :: dis, g, k, k_target 
    character(len=*) :: boundary_type_rmax, metric
    type(lines_method) :: this, r_min
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
    this%a = a
    this%M_bh = M_bh 


    this%n_det = 5


    this%xmin = xmin
    this%xmax = xmax

    this%boundary_type_rmax =  boundary_type_rmax

    print*, 'Boundary type r_max = ', this%boundary_type_rmax 

    call this%mesh_refinement_memory()
    call this%mesh_refinement_boundary(xmin, xmax)
    call this%mesh_refinement_create()
    !
    call this%lines_method_memory()
    this%tmin = tmin; this%tmax = tmax; this%CFL = CFL

      dxx = this%x(this%iR+this%g_pts,1,1)
      call MPI_Bcast(dxx, 1, MPI_DOUBLE, nproc-1, MPI_COMM_WORLD, ierr)
      this%r_max = dxx


    

    this%index_bounds(1,1) = this%iL - this%g_pts !index x left
    this%index_bounds(1,2) = this%iR + this%g_pts !index x right
    this%index_bounds(2,1) = this%jL - this%g_pts !index y left
    this%index_bounds(2,2) = this%jR + this%g_pts !index y right
    this%index_bounds(3,1) = this%kL !index z left
    this%index_bounds(3,2) = this%kR !index z left

    dxx = this%x(this%iL-this%g_pts+1,1,1) - this%x(this%iL-this%g_pts,1,1)
    call MPI_Bcast(dxx, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
    this%dt = CFL * min(dxx, exp(xmin)*this%dy, &
      exp(xmin)*dsin(this%dy/2.0d0)*this%dz) 

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
    call this%validate_inverse_metric

  end function lines_method_constructor
  !============================================================================!
  subroutine lines_method_memory(this)
    implicit none
    class(lines_method), intent(in out) :: this

    allocate(this%u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%u_p(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))


    allocate(this%u_stage(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))


    allocate(this%u_stage_2(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))


    allocate(this%k1_u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))


    allocate(this%k2_u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%k3_u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))


    allocate(this%k4_u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%rhs_u(this%nvars, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts : this%jR + this%g_pts , &
    this%kL:this%kR))

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


      call this%boundary_yL()
      call this%boundary_yR()
      call this%boundary_zL()
      call this%boundary_zR()

    end subroutine boundary

  end subroutine rk3
  !============================================================================!
  subroutine pirk2(this)
    implicit none
    class(lines_method), intent(inout) :: this

    this%u_p = this%u
  
    ! Etapa 1
    call this%rhs_phi(this%rhs_u(1,:,:,:), this%u)
    call this%rhs_psi1(this%rhs_u(2,:,:,:), this%u)
    call this%rhs_psi2(this%rhs_u(3,:,:,:), this%u)
    call this%rhs_psi3(this%rhs_u(4,:,:,:), this%u)
  
    this%u_stage(1,:,:,:) = this%u(1,:,:,:) + this%dt * this%rhs_u(1,:,:,:)
    this%u_stage(2,:,:,:) = this%u(2,:,:,:) + this%dt * this%rhs_u(2,:,:,:)
    this%u_stage(3,:,:,:) = this%u(3,:,:,:) + this%dt * this%rhs_u(3,:,:,:)
    this%u_stage(4,:,:,:) = this%u(4,:,:,:) + this%dt * this%rhs_u(4,:,:,:)
  
    call this%rhs_pi(this%rhs_u(5,:,:,:), this%u_stage)
    this%u_stage(5,:,:,:) = this%u(5,:,:,:) + this%dt * this%rhs_u(5,:,:,:)
  
    call boundary()
  
    ! Etapa 2
    call this%rhs_phi(this%rhs_u(1,:,:,:), this%u_stage)
    call this%rhs_psi1(this%rhs_u(2,:,:,:), this%u_stage)
    call this%rhs_psi2(this%rhs_u(3,:,:,:), this%u_stage)
    call this%rhs_psi3(this%rhs_u(4,:,:,:), this%u_stage)
    call this%rhs_pi(this%rhs_u(5,:,:,:), this%u_stage)
  
    this%u(1,:,:,:) = 0.5d0 * (this%u_p(1,:,:,:) + this%u_stage(1,:,:,:) + this%dt * this%rhs_u(1,:,:,:))
    this%u(2,:,:,:) = 0.5d0 * (this%u_p(2,:,:,:) + this%u_stage(2,:,:,:) + this%dt * this%rhs_u(2,:,:,:))
    this%u(3,:,:,:) = 0.5d0 * (this%u_p(3,:,:,:) + this%u_stage(3,:,:,:) + this%dt * this%rhs_u(3,:,:,:))
    this%u(4,:,:,:) = 0.5d0 * (this%u_p(4,:,:,:) + this%u_stage(4,:,:,:) + this%dt * this%rhs_u(4,:,:,:))
    this%u(5,:,:,:) = 0.5d0 * (this%u_p(5,:,:,:) + this%u_stage(5,:,:,:) + this%dt * this%rhs_u(5,:,:,:))
  
    call boundary()
  
    this%t = this%t + this%dt
  
  
  contains
  
    subroutine boundary
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
  
      call this%boundary_yL()
      call this%boundary_yR()
      call this%boundary_zL()
      call this%boundary_zR()
    end subroutine boundary
  
  end subroutine pirk2
  !============================================================================!
  subroutine pirk3(this)
    implicit none
    class(lines_method), intent(inout) :: this
    this%u_p = this%u

    ! Etapa 1
    call this%rhs_phi(this%rhs_u(1,:,:,:), this%u)
    call this%rhs_psi1(this%rhs_u(2,:,:,:), this%u)
    call this%rhs_psi2(this%rhs_u(3,:,:,:), this%u)
    call this%rhs_psi3(this%rhs_u(4,:,:,:), this%u)

    this%u_stage(1,:,:,:) = this%u(1,:,:,:) + this%dt * this%rhs_u(1,:,:,:)
    this%u_stage(2,:,:,:) = this%u(2,:,:,:) + this%dt * this%rhs_u(2,:,:,:)
    this%u_stage(3,:,:,:) = this%u(3,:,:,:) + this%dt * this%rhs_u(3,:,:,:)
    this%u_stage(4,:,:,:) = this%u(4,:,:,:) + this%dt * this%rhs_u(4,:,:,:)

    call this%rhs_pi(this%rhs_u(5,:,:,:), this%u_stage)
    this%u_stage(5,:,:,:) = this%u(5,:,:,:) + this%dt * this%rhs_u(5,:,:,:)

    call boundary()

    ! Etapa 2
    call this%rhs_phi(this%rhs_u(1,:,:,:), this%u_stage)
    call this%rhs_psi1(this%rhs_u(2,:,:,:), this%u_stage)
    call this%rhs_psi2(this%rhs_u(3,:,:,:), this%u_stage)
    call this%rhs_psi3(this%rhs_u(4,:,:,:), this%u_stage)

    this%u_stage_2(1,:,:,:) = 0.75d0 * this%u_p(1,:,:,:) + 0.25d0 * (this%u_stage(1,:,:,:) + this%dt * this%rhs_u(1,:,:,:))
    this%u_stage_2(2,:,:,:) = 0.75d0 * this%u_p(2,:,:,:) + 0.25d0 * (this%u_stage(2,:,:,:) + this%dt * this%rhs_u(2,:,:,:))
    this%u_stage_2(3,:,:,:) = 0.75d0 * this%u_p(3,:,:,:) + 0.25d0 * (this%u_stage(3,:,:,:) + this%dt * this%rhs_u(3,:,:,:))
    this%u_stage_2(4,:,:,:) = 0.75d0 * this%u_p(4,:,:,:) + 0.25d0 * (this%u_stage(4,:,:,:) + this%dt * this%rhs_u(4,:,:,:))

    call this%rhs_pi(this%rhs_u(5,:,:,:), this%u_stage_2)
    this%u_stage_2(5,:,:,:) = 0.75d0 * this%u_p(5,:,:,:) + 0.25d0 * (this%u_stage(5,:,:,:) + this%dt * this%rhs_u(5,:,:,:))

    call boundary()

    ! Etapa 3
    call this%rhs_phi(this%rhs_u(1,:,:,:), this%u_stage_2)
    call this%rhs_psi1(this%rhs_u(2,:,:,:), this%u_stage_2)
    call this%rhs_psi2(this%rhs_u(3,:,:,:), this%u_stage_2)
    call this%rhs_psi3(this%rhs_u(4,:,:,:), this%u_stage_2)
    call this%rhs_pi(this%rhs_u(5,:,:,:), this%u_stage_2)

    this%u(1,:,:,:) = (1.0d0/3.0d0) * this%u_p(1,:,:,:) + (2.0d0/3.0d0) * (this%u_stage_2(1,:,:,:) + this%dt * this%rhs_u(1,:,:,:))
    this%u(2,:,:,:) = (1.0d0/3.0d0) * this%u_p(2,:,:,:) + (2.0d0/3.0d0) * (this%u_stage_2(2,:,:,:) + this%dt * this%rhs_u(2,:,:,:))
    this%u(3,:,:,:) = (1.0d0/3.0d0) * this%u_p(3,:,:,:) + (2.0d0/3.0d0) * (this%u_stage_2(3,:,:,:) + this%dt * this%rhs_u(3,:,:,:))
    this%u(4,:,:,:) = (1.0d0/3.0d0) * this%u_p(4,:,:,:) + (2.0d0/3.0d0) * (this%u_stage_2(4,:,:,:) + this%dt * this%rhs_u(4,:,:,:))
    this%u(5,:,:,:) = (1.0d0/3.0d0) * this%u_p(5,:,:,:) + (2.0d0/3.0d0) * (this%u_stage_2(5,:,:,:) + this%dt * this%rhs_u(5,:,:,:))

    call boundary()
    this%t = this%t + this%dt

    contains

    subroutine boundary
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
  
      call this%boundary_yL()
      call this%boundary_yR()
      call this%boundary_zL()
      call this%boundary_zR()
    end subroutine boundary

  end subroutine pirk3

  !============================================================================!
  subroutine rk4(this,rhs)
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
    do rk=1, 4

      call rhs(1)

      if(rk.eq.1) then
        this%k1_u = this%rhs_u * this%dt
        this%u = this%u_p + 0.50d0 * this%k1_u
      else if(rk.eq.2) then
        this%k2_u = this%rhs_u * this%dt
        this%u = this%u_p + 0.50d0 * this%k2_u
      else if(rk.eq.3) then
        this%k3_u = this%rhs_u * this%dt
        this%u = this%u_p + this%k3_u
      else
        this%k4_u = this%rhs_u * this%dt
        this%u = this%u_p + (this%k1_u + 2.0d0 * this%k2_u + 2.0d0 * this%k3_u + this%k4_u) / 6.0d0
      end if
  
      call boundary
    end do
  
    this%t = this%t + this%dt

    contains

    subroutine boundary


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
      call this%boundary_yL()
      call this%boundary_yR()
      call this%boundary_zL()
      call this%boundary_zR()
      

    end subroutine boundary



  end subroutine rk4 
  subroutine send_recv_xL(this)
    implicit none
    class(lines_method), intent(in out) :: this


    count = size(this%u(:,this%iL:this%iL+this%g_pts,:,:))

    call MPI_SENDRECV(&
    this%u(:,this%iL:this%iL+this%g_pts,:,:), &
    count, MPI_DOUBLE_COMPLEX, rank -1, 1, &  ! send
    this%u(:,this%iL-this%g_pts:this%iL,:,:),&
    count, MPI_DOUBLE_COMPLEX, rank -1, 0, &  ! recv
    MPI_COMM_WORLD, status, ierr)

  end subroutine send_recv_xL
  !============================================================================!
  subroutine send_recv_xR(this)
    implicit none
    class(lines_method), intent(in out) :: this


    count = size(this%u(:,this%iR:this%iR+this%g_pts,:,:))

    call MPI_SENDRECV(&
    this%u(:,this%iR-this%g_pts:this%iR,:,:), &
    count, MPI_DOUBLE_COMPLEX, rank +1, 0, &  ! send
    this%u(:,this%iR:this%iR+this%g_pts,:,:), &
    count, MPI_DOUBLE_COMPLEX, rank +1, 1, &  ! recv
    MPI_COMM_WORLD, status, ierr)

  end subroutine send_recv_xR
  !============================================================================!
  subroutine boundary_xL(this)
    implicit none
    class(lines_method), intent(in out) :: this
    integer :: i
    real(kind=8) :: x0,x1,x2,x3
    real(kind=8) :: L1,L2,L3
    ! boundary xL == r_min, radius of excision
    ! third order extrapolations

    this%u(:,this%iL-this%g_pts,:,:) = 3.0d0*this%u(:,this%iL-this%g_pts+1,:,:) &
          -3.0d0*this%u(:,this%iL-this%g_pts+2,:,:) + this%u(:,this%iL-this%g_pts+3,:,:)


  end subroutine boundary_xL
  !============================================================================!
  subroutine boundary_xR(this)
    implicit none
    class(lines_method), intent(in out) :: this
    if(this%boundary_type_rmax.eq.'outflow') then
      this%u(:,this%iR+this%g_pts,:,:) = this%u(:,this%iR+this%g_pts-1,:,:)
    else if(this%boundary_type_rmax.eq.'inflow') then
      this%u(:,this%iR+this%g_pts,:,:) = this%u_p(:,this%iR+this%g_pts,:,:)
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
    integer :: k, g, dis, k_target

    dis = (this%Nz) / 2
    do k = this%kL , this%kR 
      k_target = mod(k + dis, this%Nz - 1) 
      do g = this%g_pts, 0, -1
        ! if(rank.eq.master) then 
        !   print*, this%jL-g, this%jL+g+1, k, k_target
        !   print*, this%y(1,this%jL-g,k), this%y(1,this%jL+g+1,k_target), &
        !           this%z(1,this%jL-g,k) - this%z(1,this%jL+g+1,k_target)
        ! end if 
        this%u(:,:,this%jL-g,k) = this%u(:,:,this%jL+g+1,k_target) 
        !this%u(3,:,this%jL-g,k) = -this%u(3,:,this%jL+g+1,k_target) 
        !this%u(4,:,this%jL-g,k) = -this%u(4,:,this%jL+g+1,k_target) 
      end do
    end do


  end subroutine boundary_yL
  !============================================================================!
  subroutine boundary_yR(this)
    implicit none
    class(lines_method), intent(in out) :: this
    integer :: k, g, dis, k_target

    dis = (this%Nz) / 2
    do k = this%kL , this%kR
      k_target = mod(k + dis, this%Nz - 1)  
      do g = this%g_pts, 0, -1
        ! if(rank.eq.master) then 
        !   print*, this%jR +g, this%jR -g-1, k, k_target
        !   print*, this%y(1,this%jR +g,k),this%y(1,this%jR -g-1,k_target), & 
        !           this%z(1,this%jR +g,k) - this%z(1,this%jR -g-1,k_target)
        ! end if 
        this%u(:,:,this%jR+g,k) = this%u(:,:,this%jR-g-1,k_target)
        !this%u(3,:,this%jR+g,k) = -this%u(3,:,this%jR-g-1,k_target) 
        !this%u(4,:,this%jR+g,k) = -this%u(4,:,this%jR-g-1,k_target) 
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
          do j = this%jL - this%g_pts , this%jR + this%g_pts 
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
          do j = this%jL - this%g_pts , this%jR + this%g_pts 
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
  !============================================================================!

    ! rhs para \phi
  subroutine rhs_phi(this, rhs, vars)
    complex(kind=8), intent(out) :: rhs(:,:,:)
    complex(kind=8), intent(in)  :: vars(:,:,:,:)
    class(lines_method), intent(in out) :: this 
    rhs = this%alpha * vars(5,:,:,:) + this%beta(1,:,:,:) * vars(2,:,:,:)
  end subroutine rhs_phi

  ! rhs para \psi_1 := \partial_1 \phi
  subroutine rhs_psi1(this, rhs, vars)
    use finite_differences, &
      dx1 => first_derivative_x_2, &
      ad_1 => advec_x
    implicit none
    class(lines_method), intent(in out) :: this 
    complex(kind=8), intent(out) :: rhs(:,:,:)
    complex(kind=8), intent(in)  :: vars(:,:,:,:)
    rhs = dx1(this%alpha * vars(5,:,:,:), this%dx, this%index_bounds) &
          + this%beta(1,:,:,:) * ad_1(vars(2,:,:,:), this%beta(1,:,:,:), this%dx, this%index_bounds) &
         + this%dx_beta * vars(2,:,:,:)
  end subroutine rhs_psi1

  ! rhs para \psi_2 := \partial_2 \phi
  subroutine rhs_psi2(this, rhs, vars)
    use finite_differences, &
    dx2 => first_derivative_y_2
    implicit none
    class(lines_method), intent(in out) :: this 
    complex(kind=8), intent(out) :: rhs(:,:,:)
    complex(kind=8), intent(in)  :: vars(:,:,:,:)
      rhs = dx2( this%alpha * vars(5,:,:,:) + this%beta(1,:,:,:) * vars(2,:,:,:), this%dy, this%index_bounds, .false.)
  end subroutine rhs_psi2

  ! rhs para \psi_3 := \partial_3 \phi
  subroutine rhs_psi3(this, rhs, vars)
    use finite_differences, &
    dx3 => first_derivative_z_2
    implicit none
    class(lines_method), intent(in out) :: this 
    complex(kind=8), intent(out) :: rhs(:,:,:)
    complex(kind=8), intent(in)  :: vars(:,:,:,:)
      rhs = dx3(this%alpha * vars(5,:,:,:) + this%beta(1,:,:,:) * vars(2,:,:,:), this%dz, this%index_bounds)

  end subroutine rhs_psi3

  ! rhs para \pi
  subroutine rhs_pi(this, rhs, vars)
    use finite_differences, &
    dx1 => first_derivative_x_2, &
    dx2 => first_derivative_y_2, &
    dx3 => first_derivative_z_2, &
    ad_1 => advec_x
    implicit none
    class(lines_method), intent(in out) :: this 
    complex(kind=8), intent(out) :: rhs(:,:,:)
    complex(kind=8), intent(in)  :: vars(:,:,:,:)

    if(this%metric =='Kerr_Schild') then
        rhs = dx1(this%alpha*this%sqrtg*(this%gamma_up(1,1,:,:,:)*vars(2,:,:,:) &
        + this%gamma_up(1,3,:,:,:)*vars(4,:,:,:)), this%dx, this%index_bounds)/this%sqrtg &
            + dx2(this%alpha*this%sqrtg*this%gamma_up(2,2,:,:,:)*vars(3,:,:,:) &
                  , this%dy, this%index_bounds, .false.)/this%sqrtg &
            + dx3(this%alpha*this%sqrtg*(this%gamma_up(1,3,:,:,:)*vars(2,:,:,:) & 
                  + this%gamma_up(3,3,:,:,:)*vars(4,:,:,:)), this%dz, this%index_bounds)/this%sqrtg &
            + vars(5,:,:,:) * this%dx_beta &
            + this%beta(1,:,:,:) * ad_1(vars(5,:,:,:)*this%sqrtg, this%beta(1,:,:,:), this%dx, this%index_bounds)/this%sqrtg &
            - this%alpha * this%m**2 * vars(1,:,:,:)
    else if(this%metric =='Schw_EF') then
        rhs = dx1(this%alpha*this%sqrtg*this%gamma_up(1,1,:,:,:)*vars(2,:,:,:),this%dx, this%index_bounds)/this%sqrtg &
            + dx2(this%alpha*this%sqrtg*this%gamma_up(2,2,:,:,:)*vars(3,:,:,:), this%dy, this%index_bounds, .false.)/this%sqrtg &
            + dx3(this%alpha*this%gamma_up(3,3,:,:,:)*vars(4,:,:,:), this%dz, this%index_bounds) &
            + vars(5,:,:,:) * this%dx_beta &
            + this%beta(1,:,:,:) * ad_1(vars(5,:,:,:)*this%sqrtg, this%beta(1,:,:,:), this%dx, this%index_bounds)/this%sqrtg &
            - this%alpha * this%m**2 * vars(1,:,:,:)
    else 
      print*, "Metric: ", this%metric, "not implemented yet"
      print*, 'try with: Kerr_Schild or Schw_EF'
      stop
      
    end if



  end subroutine rhs_pi
  !============================================================================!

  


end module lines_method_lib

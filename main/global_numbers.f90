module global_numbers
  use strings_lib
  use diagnostic_lib
  use hdf5_lib
  implicit none

  !
  ! for this
  !
  integer :: iL, iR, jL, jR, kL, kR
  integer :: n_refinements, nvars
  integer :: Nx, Ny, Nz, g_pts
  real(kind=8) :: xmin, xmax
  real(kind=8) :: ymin, ymax
  real(kind=8) :: zmin, zmax
  real(kind=8) :: tmin, tmax, CFL
  character(len=20) :: reconstructor
  character(len=20) :: boundary_type_rmax
  !
  ! save data
  !
  integer :: save_0d, save_1d, save_2d, save_3d
  integer :: i_save, j_save, k_save
  real(kind=8) :: x_save, y_save, z_save
  real(kind=8) :: save_
  logical :: proc_x_save, proc_y_save, proc_z_save
  !
  ! cpu time
  !
  real(kind=8) :: t_initial, t_final, t_cpu
  !
  ! methods and system
  !
  type(scalar_field) :: MoL
  !
  ! checkpoint
  !
  logical :: read_from_checkpoint
  !
  ! intial conditions
  !
  character(len=20) :: initial_conditions_type
  integer :: l_max_intial
  !
  !parameters for space-time
  !
  real(kind=8) :: a, M_bh
  character(len=20) :: metric
  !
  ! boson mass
  !
  real(kind=8) :: m
  !
  ! Cos modulado gaussiana
  !
  real(kind=8) :: mu, sigma, k_f
  !
  ! Detectors
  !
  integer :: n_det
  real(kind=8) :: first_det
  real(kind=8) :: space_det
contains

  subroutine read_data
    use ODE, only: linspace
    implicit none
    character(len=40) :: input_file
    type(string) :: filename
    type(string) :: filename_split(2)
    type(string) :: output_folder

    !
    ! create the list of parameters
    !
    namelist/input/ nvars, Nx, Ny, Nz, g_pts, &
    xmin, xmax, tmax, CFL, &
    reconstructor, &
    save_0d, save_1d, save_2d, save_3d, &
    save_, &
    x_save, y_save, z_save, &
    boundary_type_rmax, &
    l_max_intial, &
    initial_conditions_type, &
    metric, & 
    a, M_bh, &
    m, &
    mu, sigma, k_f, &
    n_det, first_det, space_det, &
    read_from_checkpoint
    !
    ! read the name of the parameter file
    !
    call get_command_argument(1, input_file)
    input_file = trim(input_file)
    !
    ! now the list of parameters is read
    !
    open(100, file = input_file, action = 'read')
      read(100, nml = input)
    close(100)
    !
    ! eliminate the spaces
    !
    filename = remove(input_file, '')
    !
    ! remove the extension ".something"
    !
    filename_split = filename%split('.')
    output_folder = filename_split(1) + '/'
    !
    ! create the folder where the output data is saved.
    !
    call system('mkdir -p '//trim(output_folder%string_data))
    call CHDIR(output_folder%string_data)

  end subroutine read_data

  subroutine create_data
    implicit none
    integer :: i,l
    !
    ! lines method is created
    !
    if(metric=='MKS') then 
      xmin = log(xmin)
      xmax = log(xmax)
    end if
      ! xmin = log(xmin)
      ! xmax = log(xmax)

    MoL = scalar_field_constructor(nvars = nvars, g_pts = g_pts, &
    xmin = xmin, xmax = xmax, Nx = Nx, &
    Ny = Ny, &
    Nz = Nz, &
    tmin = tmin, tmax = tmax, CFL = CFL, &
    boundary_type_rmax=boundary_type_rmax, &
    metric = metric, a=a, M_bh = M_bh, & 
    m=m, & 
    save=save_, r_max = xmax, n_det=n_det)
    !
    ! for save data
    !
    iL = MoL%iL; iR = MoL%iR
    jL = MoL%jL; jR = MoL%jR
    kL = MoL%kL; kR = MoL%kR

    proc_x_save = .false.
    proc_y_save = .false.
    proc_z_save = .false.

    do l=1, MoL%n_det
        
        do i=1, MoL%Nx-1
          if(MoL%r(i,1,1)<=MoL%r_det(l) .and. MoL%r_det(l)<=MoL%r(i+1,1,1)) then
            MoL%index_det(l) = i
            print*, 'rank = ',rank, 'index = ', MoL%index_det(l), 'radius = ', MoL%r_det(l)
            proc_x_save = .true.
            exit
          end if
        end do
        

    end do


    if(MoL%metric == 'MKS') then 
    call create_hdf5_file('phi.h5')
      call create_hdf5_group('phi.h5', '/refinement_1' )

      call write_hdf5('phi.h5', '/refinement_1', 'Xcoord', exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('phi.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('phi.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('rho.h5')
      call create_hdf5_group('rho.h5', '/refinement_1' )

      call write_hdf5('rho.h5', '/refinement_1', 'Xcoord', exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('rho.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('rho.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      call create_hdf5_file('psix.h5')
      call create_hdf5_group('psix.h5', '/refinement_1' )

      call write_hdf5('psix.h5', '/refinement_1', 'Xcoord', exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psix.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psix.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('psiy.h5')
      call create_hdf5_group('psiy.h5', '/refinement_1' )

      call write_hdf5('psiy.h5', '/refinement_1', 'Xcoord', exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiy.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiy.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      call create_hdf5_file('psiz.h5')
      call create_hdf5_group('psiz.h5', '/refinement_1' )

      call write_hdf5('psiz.h5', '/refinement_1', 'Xcoord', exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiz.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiz.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('pi.h5')
      call create_hdf5_group('pi.h5', '/refinement_1' )

      call write_hdf5('pi.h5', '/refinement_1', 'Xcoord', exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('pi.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('pi.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('jr.h5')
      call create_hdf5_group('jr.h5', '/refinement_1' )

      call write_hdf5('jr.h5', '/refinement_1', 'Xcoord',exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('jr.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('jr.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('jtheta.h5')
      call create_hdf5_group('jtheta.h5', '/refinement_1' )

      call write_hdf5('jtheta.h5', '/refinement_1', 'Xcoord',exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('jtheta.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('jtheta.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('jphi.h5')
      call create_hdf5_group('jphi.h5', '/refinement_1' )

      call write_hdf5('jphi.h5', '/refinement_1', 'Xcoord',exp(MoL%x(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1)) )
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('jphi.h5', '/refinement_1', 'Ycoord', MoL%y(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('jphi.h5', '/refinement_1', 'Zcoord', MoL%z(MoL%index_det(1): &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      
    else 
      

      call create_hdf5_file('phi.h5')
      call create_hdf5_group('phi.h5', '/refinement_1' )

      call write_hdf5('phi.h5', '/refinement_1', 'Xcoord', MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('phi.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('phi.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      call create_hdf5_file('rho.h5')
      call create_hdf5_group('rho.h5', '/refinement_1' )

      call write_hdf5('rho.h5', '/refinement_1', 'Xcoord', MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('rho.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('rho.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('psix.h5')
      call create_hdf5_group('psix.h5', '/refinement_1' )

      call write_hdf5('psix.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1)) 
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psix.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psix.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('psiy.h5')
      call create_hdf5_group('psiy.h5', '/refinement_1' )

      call write_hdf5('psiy.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1)) 
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiy.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiy.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      call create_hdf5_file('psiz.h5')
      call create_hdf5_group('psiz.h5', '/refinement_1' )

      call write_hdf5('psiz.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1)) 
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiz.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('psiz.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('pi.h5')
      call create_hdf5_group('pi.h5', '/refinement_1' )

      call write_hdf5('pi.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1)) 
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('pi.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('pi.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call create_hdf5_file('rePhi.h5')
      call create_hdf5_group('rePhi.h5', '/refinement_1' )

      call write_hdf5('rePhi.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1)) 
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('rePhi.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('rePhi.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      call create_hdf5_file('imPhi.h5')
      call create_hdf5_group('imPhi.h5', '/refinement_1' )

      call write_hdf5('imPhi.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      iR,jL:jR,kL:kR-1)) 
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('imPhi.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('imPhi.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR,jL:jR,kL:kR-1))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


      ! call create_hdf5_file('jr.h5')
      ! call create_hdf5_group('jr.h5', '/refinement_1' )

      ! call write_hdf5('jr.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1)) 
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! call write_hdf5('jr.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1))
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! call write_hdf5('jr.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1))
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! call create_hdf5_file('jtheta.h5')
      ! call create_hdf5_group('jtheta.h5', '/refinement_1' )

      ! call write_hdf5('jtheta.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1)) 
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! call write_hdf5('jtheta.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1))
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! call write_hdf5('jtheta.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1))
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      ! call create_hdf5_file('jphi.h5')
      ! call create_hdf5_group('jphi.h5', '/refinement_1' )

      ! call write_hdf5('jphi.h5', '/refinement_1', 'Xcoord',MoL%x(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1)) 
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! call write_hdf5('jphi.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1))
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      ! call write_hdf5('jphi.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      ! iR,jL:jR,kL:kR-1))
      ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    end if

  end subroutine create_data

end module global_numbers

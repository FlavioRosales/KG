module global_numbers
  use strings_lib
  use lines_method_lib
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
  logical :: proc_x_save, proc_y_save, proc_z_save
  !
  ! cpu time
  !
  real(kind=8) :: t_initial, t_final, t_cpu
  !
  ! methods and system
  !
  type(lines_method) :: MoL
  !
  ! checkpoint
  !
  logical :: read_from_checkpoint
  !
  ! intial conditions
  !
  character(len=20) :: initial_conditions_type
  !
  !parameters for space-time
  !
  real(kind=8) :: a, M_bh
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
    x_save, y_save, z_save, &
    boundary_type_rmax, &
    initial_conditions_type, & 
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
    real(kind=8), dimension(n_det) :: r_det
    integer :: i,l
    !
    ! lines method is created
    !
    MoL = lines_method(nvars = nvars, g_pts = g_pts, &
    xmin = xmin, xmax = xmax, Nx = Nx, &
    Ny = Ny, &
    Nz = Nz, &
    tmin = tmin, tmax = tmax, CFL = CFL, &
    boundary_type_rmax=boundary_type_rmax, &
    metric = 'Kerr_Schild', & 
    !metric = 'Minkowski', & 
    m=m, & 
    n_det = n_det, first_det = first_det, space_det = space_det)
    !
    ! for save data
    !
    iL = MoL%iL; iR = MoL%iR
    jL = MoL%jL; jR = MoL%jR
    kL = MoL%kL; kR = MoL%kR

    proc_x_save = .false.
    proc_y_save = .false.
    proc_z_save = .false.

    print*, 'number of detectors', MoL%n_det
    !print*,'location of detectors', MoL%r_det(:)


    if(rank.eq.0) then 
     do i=MoL%iL-MoL%g_pts, MoL%iR
           do l=1, n_det
           if(MoL%x(i,1,1)<=MoL%r_det(l) .and. MoL%r_det(l)<=MoL%x(i+1,1,1)) then
            ! i_save = i
             proc_x_save = .true.
             print*, rank
           end if
            end do
     end do
    else 
      do i=MoL%iL, MoL%iR
        do l=1, n_det
        if(MoL%x(i,1,1)<=MoL%r_det(l) .and. MoL%r_det(l)<=MoL%x(i+1,1,1)) then
         ! i_save = i
          proc_x_save = .true.
          print*, rank
        end if
      end do
  end do
    end if
    !
    ! hdf5 files
    !    !final     :: mesh_refinement_destructor

    !
    ! create files
    !
    call create_hdf5_file('phi.h5')
      call create_hdf5_group('phi.h5', '/refinement_1' )

      call write_hdf5('phi.h5', '/refinement_1', 'Xcoord', MoL%x(iL-g_pts: &
      iR+g_pts,jL:jR,kL:kR))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('phi.h5', '/refinement_1', 'Ycoord', MoL%y(iL-g_pts: &
      iR+g_pts,jL:jR,kL:kR))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call write_hdf5('phi.h5', '/refinement_1', 'Zcoord', MoL%z(iL-g_pts: &
      iR+g_pts,jL:jR,kL:kR))
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)


    !
    ! Detectors
    ! 

  end subroutine create_data

end module global_numbers

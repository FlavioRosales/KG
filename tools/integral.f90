module integral
  implicit none

contains
  !=============================================================================!
  function trapezium(f3d, dx, dy, dz) result(f0d)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(:,:,:), intent(in) :: f3d(:,:,:)
    real(kind=8) :: dx, dy, dz
    !
    ! output
    !
    real(kind=8) :: f0d
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: f2d(size(f3d,1), size(f3d,2)), f1d(size(f3d,1))

    f2d = 0.0d0
    do i=1, size(f3d,3)-1
      f2d = f2d + f3d(:,:,i) + f3d(:,:,i+1)
    end do
    f2d = 0.50d0 * f2d * dz

    f1d = 0.0d0
    do i=1, size(f3d,2)-1
      f1d = f1d + f2d(:,i) + f2d(:,i+1)
    end do
    f1d = 0.50d0 * f1d * dy

    f0d = 0.0d0
    do i=1, size(f3d,1)-1
      f0d = f0d + f1d(i) + f1d(i+1)
    end do
    f0d = 0.50d0 * f0d * dx

  end function trapezium
  !=============================================================================!
  function trapezium_2D(f, dx, dy) result(f0d)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(:,:), intent(in) :: f
    real(kind=8), intent(in) :: dx, dy
    !
    ! output
    !
    real(kind=8) :: f0d
    !
    ! internal
    !
    integer :: i,j,m,n
    
    m = size(f, 1)
    n = size(f, 2)

    ! Inicializar la integral
    f0d = 0.0d0

    do i = 1, m
      do j = 1, n
          if ((i == 1 .or. i == m) .and. (j == 1 .or. j == n)) then
              f0d = f0d + f(i,j) * 0.25* dx * dy
          else if ((i == 1 .or. i == m) .or. (j == 1 .or. j == n)) then
              f0d = f0d + f(i,j) * 0.5* dx * dy
          else
              f0d = f0d + f(i,j)* dx * dy
          end if
      end do
  end do

  end function trapezium_2D
  !=============================================================================!
  function simpson_one_third(f3d, dx, dy, dz) result(f0d)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f3d(:,:,:), dx, dy, dz
    !
    ! output
    !
    real(kind=8) :: f0d
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: f2d(size(f3d,1),size(f3d,2)), f1d(size(f3d,1))

    f2d = 0.0d0
    do i=2, size(f3d,3)-1, 2
      f2d = f2d + f3d(:,:,i-1) + 4.0d0 * f3d(:,:,i) + f3d(:,:,i+1)
    end do
    f2d = f2d * dz / 3.0d0

    f1d = 0.0d0
    do i=2, size(f3d,2)-1, 2
      f1d = f1d + f2d(:,i-1) + 4.0d0 * f2d(:,i) + f2d(:,i+1)
    end do
    f1d = f1d * dy / 3.0d0

    f0d = 0.0d0
    do i=2, size(f3d,1)-1, 2
      f0d = f0d + f1d(i-1) + 4.0d0 * f1d(i) + f1d(i+1)
    end do
    f0d = f0d * dx / 3.0d0

  end function simpson_one_third
  !=============================================================================!
  function simpson_three_eighths(f3d, dx, dy, dz) result(f0d)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f3d(:,:,:), dx, dy, dz
    !
    ! output
    !
    real(kind=8) :: f0d
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: f2d(size(f3d,1),size(f3d,2)), f1d(size(f3d,1))

    f2d = 0.0d0
    do i=3, size(f3d,3)-1, 3
      f2d = f2d + f3d(:,:,i-2) + 3.0d0 * f3d(:,:,i-1) + 3.0d0 * f3d(:,:,i) + f3d(:,:,i+1)
    end do
    f2d = f2d * dz * 3.0d0 / 8.0d0

    f1d = 0.0d0
    do i=3, size(f3d,2)-1, 3
      f1d = f1d + f2d(:,i-2) + 3.0d0 * f2d(:,i-1) + 3.0d0 * f2d(:,i) + f2d(:,i+1)
    end do
    f1d = f1d * dy * 3.0d0 / 8.0d0

    f0d = 0.0d0
    do i=3, size(f3d,1)-1, 3
      f0d = f0d + f1d(i-2) + 3.0d0 * f1d(i-1) + 3.0d0 * f1d(i) + f1d(i+1)
    end do
    f0d = f0d * dx * 3.0d0 / 8.0d0

  end function simpson_three_eighths
  !=============================================================================!
  function integrate(f3d, dx, dy, dz, g_pts) result(f0d)
    use mpi_lib
    implicit none
    !
    ! input
    !
    integer, intent(in) :: g_pts
    real(kind=8), dimension(:,:,:), intent(in) :: f3d(:,:,:)
    real(kind=8) :: dx, dy, dz
    !
    ! output
    !
    real(kind=8) :: f0d
    !
    ! internal
    !
    real(kind=8) :: int
    integer :: iL, iR, jL, jR, kL, kR

    iL = 1 + g_pts; iR = size(f3d,1) - g_pts
    jL = 1 + g_pts; jR = size(f3d,2) - g_pts
    kL = 1 + g_pts; kR = size(f3d,3) - g_pts

    f0d = 0.0d0
    int = trapezium(f3d = f3d(iL:iR,jL:jR,kL:kR),dx =  dx, dy = dy, dz = dz)

    call MPI_ALLREDUCE(int, f0d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  end function integrate
  !=============================================================================!

end module integral

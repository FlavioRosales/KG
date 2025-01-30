module finite_differences
  implicit none

contains

  !=============================================================================!
  ! first derivatives X
  !=============================================================================!
  function first_derivative_x_2(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,1) - 1
    res = 0.50d0 / h

    df(0,:,:) = (-3.0d0 * f(0,:,:) + 4.0d0 * f(1,:,:) - f(2,:,:)) * res
    do i=1, N-1
      df(i,:,:) = (-f(i-1,:,:) + f(i+1,:,:)) * res
    end do
    df(N,:,:) = (+3.0d0 * f(N,:,:) - 4.0d0 * f(N-1,:,:) + f(N-2,:,:)) * res

  end function first_derivative_x_2
  !=============================================================================!
  function first_derivative_x_4(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,1) - 1
    res = 1.0d0 / (12.0d0 * h)

    df(0,:,:) = (-25.0d0 * f(0,:,:) + 48.0d0 * f(1,:,:) - 36.0d0 * f(2,:,:) + 16.0d0 * f(3,:,:) - 3.0d0 * f(4,:,:)) * res
    df(1,:,:) = (-03.0d0 * f(0,:,:) - 10.0d0 * f(1,:,:) + 18.0d0 * f(2,:,:) - 06.0d0 * f(3,:,:) + 1.0d0 * f(4,:,:)) * res
    do i=2, N-2
      df(i,:,:) = (f(i-2,:,:) - 8.0d0 * f(i-1,:,:) + 8.0d0 * f(i+1,:,:) - f(i+2,:,:)) * res
    end do
    df(N-1,:,:) = (+03.0d0 * f(N,:,:) + 10.0d0 * f(N-1,:,:) - 18.0d0 * f(N-2,:,:) + 06.0d0 * f(N-3,:,:) - 1.0d0 * f(N-4,:,:)) * res
    df(N-0,:,:) = (+25.0d0 * f(N,:,:) - 48.0d0 * f(N-1,:,:) + 36.0d0 * f(N-2,:,:) - 16.0d0 * f(N-3,:,:) + 3.0d0 * f(N-4,:,:)) * res

  end function first_derivative_x_4
  !=============================================================================!
  function first_derivative_x_6(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,1) - 1
    res = 1.0d0 / (60.0d0 * h)

    df(0,:,:) = (-147.0d0 * f(0,:,:) + 360.0d0 * f(1,:,:) - 450.0d0 * f(2,:,:) &
    + 400.0d0 * f(3,:,:) - 225.0d0 * f(4,:,:) + 72.0d0 * f(5,:,:) - 10.0d0 * f(6,:,:)) * res
    df(1,:,:) = (-010.0d0 * f(0,:,:) - 077.0d0 * f(1,:,:) + 150.0d0 * f(2,:,:) &
    - 100.0d0 * f(3,:,:) + 050.0d0 * f(4,:,:) - 15.0d0 * f(5,:,:) + 02.0d0 * f(6,:,:)) * res
    df(2,:,:) = (+002.0d0 * f(0,:,:) - 024.0d0 * f(1,:,:) - 035.0d0 * f(2,:,:) &
    + 080.0d0 * f(3,:,:) - 030.0d0 * f(4,:,:) + 08.0d0 * f(5,:,:) - 01.0d0 * f(6,:,:)) * res
    do i=3, N-3
      df(i,:,:) = res * &
      (f(i-3,:,:) + 9.0d0 * f(i-2,:,:) - 45.0d0 * f(i-1,:,:) + 45.0d0 * f(i+1,:,:) - 9.0d0 * f(i+2,:,:) + f(i+3,:,:))
    end do
    df(N-2,:,:) = (-002.0d0 * f(N,:,:) + 024.0d0 * f(N-1,:,:) + 035.0d0 * f(N-2,:,:) &
    - 080.0d0 * f(N-3,:,:) + 030.0d0 * f(N-4,:,:) - 08.0d0 * f(N-5,:,:) + 01.0d0 * f(N-6,:,:)) * res
    df(N-1,:,:) = (+010.0d0 * f(N,:,:) + 077.0d0 * f(N-1,:,:) - 150.0d0 * f(N-2,:,:) &
    + 100.0d0 * f(N-3,:,:) - 050.0d0 * f(N-4,:,:) + 15.0d0 * f(N-5,:,:) - 02.0d0 * f(N-6,:,:)) * res
    df(N-0,:,:) = (+147.0d0 * f(N,:,:) - 360.0d0 * f(N-1,:,:) + 450.0d0 * f(N-2,:,:) &
    - 400.0d0 * f(N-3,:,:) + 225.0d0 * f(N-4,:,:) - 72.0d0 * f(N-5,:,:) + 10.0d0 * f(N-6,:,:)) * res

  end function first_derivative_x_6
  !=============================================================================!
  ! first derivatives Y
  !=============================================================================!
  function first_derivative_y_2(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,2) - 1
    res = 0.50d0 / h

    df(:,0,:) = (-3.0d0 * f(:,0,:) + 4.0d0 * f(:,1,:) - f(:,2,:)) * res
    do i=1, N-1
      df(:,i,:) = (-f(:,i-1,:) + f(:,i+1,:)) * res
    end do
    df(:,N,:) = (+3.0d0 * f(:,N,:) - 4.0d0 * f(:,N-1,:) + f(:,N-2,:)) * res

  end function first_derivative_y_2
  !=============================================================================!
  function first_derivative_y_4(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,2) - 1
    res = 1.0d0 / (12.0d0 * h)

    df(:,0,:) = (-25.0d0 * f(:,0,:) + 48.0d0 * f(:,1,:) - 36.0d0 * f(:,2,:) + 16.0d0 * f(:,3,:) - 3.0d0 * f(:,4,:)) * res
    df(:,1,:) = (-03.0d0 * f(:,0,:) - 10.0d0 * f(:,1,:) + 18.0d0 * f(:,2,:) - 06.0d0 * f(:,3,:) + 1.0d0 * f(:,4,:)) * res
    do i=2, N-2
      df(:,i,:) = (f(:,i-2,:) - 8.0d0 * f(:,i-1,:) + 8.0d0 * f(:,i+1,:) - f(:,i+2,:)) * res
    end do
    df(:,N-1,:) = (+03.0d0 * f(:,N,:) + 10.0d0 * f(:,N-1,:) - 18.0d0 * f(:,N-2,:) + 06.0d0 * f(:,N-3,:) - 1.0d0 * f(:,N-4,:)) * res
    df(:,N-0,:) = (+25.0d0 * f(:,N,:) - 48.0d0 * f(:,N-1,:) + 36.0d0 * f(:,N-2,:) - 16.0d0 * f(:,N-3,:) + 3.0d0 * f(:,N-4,:)) * res

  end function first_derivative_y_4
  !=============================================================================!
  function first_derivative_y_6(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,2) - 1
    res = 1.0d0 / (60.0d0 * h)

    df(:,0,:) = (-147.0d0 * f(:,0,:) + 360.0d0 * f(:,1,:) - 450.0d0 * f(:,2,:) &
    + 400.0d0 * f(:,3,:) - 225.0d0 * f(:,4,:) + 72.0d0 * f(:,5,:) - 10.0d0 * f(:,6,:)) * res
    df(:,1,:) = (-010.0d0 * f(:,0,:) - 077.0d0 * f(:,1,:) + 150.0d0 * f(:,2,:) &
    - 100.0d0 * f(:,3,:) + 050.0d0 * f(:,4,:) - 15.0d0 * f(:,5,:) + 02.0d0 * f(:,6,:)) * res
    df(:,2,:) = (+002.0d0 * f(:,0,:) - 024.0d0 * f(:,1,:) - 035.0d0 * f(:,2,:) &
    + 080.0d0 * f(:,3,:) - 030.0d0 * f(:,4,:) + 08.0d0 * f(:,5,:) - 01.0d0 * f(:,6,:)) * res
    do i=3, N-3
      df(:,i,:) = res * &
      (f(:,i-3,:) + 9.0d0 * f(:,i-2,:) - 45.0d0 * f(:,i-1,:) + 45.0d0 * f(:,i+1,:) - 9.0d0 * f(:,i+2,:) + f(:,i+3,:))
    end do
    df(:,N-2,:) = (-002.0d0 * f(:,N,:) + 024.0d0 * f(:,N-1,:) + 035.0d0 * f(:,N-2,:) &
    - 080.0d0 * f(:,N-3,:) + 030.0d0 * f(:,N-4,:) - 08.0d0 * f(:,N-5,:) + 01.0d0 * f(:,N-6,:)) * res
    df(:,N-1,:) = (+010.0d0 * f(:,N,:) + 077.0d0 * f(:,N-1,:) - 150.0d0 * f(:,N-2,:) &
    + 100.0d0 * f(:,N-3,:) - 050.0d0 * f(:,N-4,:) + 15.0d0 * f(:,N-5,:) - 02.0d0 * f(:,N-6,:)) * res
    df(:,N-0,:) = (+147.0d0 * f(:,N,:) - 360.0d0 * f(:,N-1,:) + 450.0d0 * f(:,N-2,:) &
    - 400.0d0 * f(:,N-3,:) + 225.0d0 * f(:,N-4,:) - 72.0d0 * f(:,N-5,:) + 10.0d0 * f(:,N-6,:)) * res

  end function first_derivative_y_6
  !=============================================================================!
  ! first derivatives Z
  !=============================================================================!
  function first_derivative_z_2(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,3) - 1
    res = 0.50d0 / h

    df(:,:,0) = (-3.0d0 * f(:,:,0) + 4.0d0 * f(:,:,1) - f(:,:,2)) * res
    do i=1, N-1
      df(:,:,i) = (-f(:,:,i-1) + f(:,:,i+1)) * res
    end do
    df(:,:,N) = (+3.0d0 * f(:,:,N) - 4.0d0 * f(:,:,N-1) + f(:,:,N-2)) * res

  end function first_derivative_z_2
  !=============================================================================!
  function first_derivative_z_4(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,3) - 1
    res = 1.0d0 / (12.0d0 * h)

    df(:,:,0) = (-25.0d0 * f(:,:,0) + 48.0d0 * f(:,:,1) - 36.0d0 * f(:,:,2) + 16.0d0 * f(:,:,3) - 3.0d0 * f(:,:,4)) * res
    df(:,:,1) = (-03.0d0 * f(:,:,0) - 10.0d0 * f(:,:,1) + 18.0d0 * f(:,:,2) - 06.0d0 * f(:,:,3) + 1.0d0 * f(:,:,4)) * res
    do i=2, N-2
      df(:,:,i) = (f(:,:,i-2) - 8.0d0 * f(:,:,i-1) + 8.0d0 * f(:,:,i+1) - f(:,:,i+2)) * res
    end do
    df(:,:,N-1) = (+03.0d0 * f(:,:,N) + 10.0d0 * f(:,:,N-1) - 18.0d0 * f(:,:,N-2) + 06.0d0 * f(:,:,N-3) - 1.0d0 * f(:,:,N-4)) * res
    df(:,:,N-0) = (+25.0d0 * f(:,:,N) - 48.0d0 * f(:,:,N-1) + 36.0d0 * f(:,:,N-2) - 16.0d0 * f(:,:,N-3) + 3.0d0 * f(:,:,N-4)) * res

  end function first_derivative_z_4
  !=============================================================================!
  function first_derivative_z_6(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,3) - 1
    res = 1.0d0 / (60.0d0 * h)

    df(:,:,0) = (-147.0d0 * f(:,:,0) + 360.0d0 * f(:,:,1) - 450.0d0 * f(:,:,2) &
    + 400.0d0 * f(:,:,3) - 225.0d0 * f(:,:,4) + 72.0d0 * f(:,:,5) - 10.0d0 * f(:,:,6)) * res
    df(:,:,1) = (-010.0d0 * f(:,:,0) - 077.0d0 * f(:,:,1) + 150.0d0 * f(:,:,2) &
    - 100.0d0 * f(:,:,3) + 050.0d0 * f(:,:,4) - 15.0d0 * f(:,:,5) + 02.0d0 * f(:,:,6)) * res
    df(:,:,2) = (+002.0d0 * f(:,:,0) - 024.0d0 * f(:,:,1) - 035.0d0 * f(:,:,2) &
    + 080.0d0 * f(:,:,3) - 030.0d0 * f(:,:,4) + 08.0d0 * f(:,:,5) - 01.0d0 * f(:,:,6)) * res
    do i=3, N-3
      df(:,:,i) = res * &
      (f(:,:,i-3) + 9.0d0 * f(:,:,i-2) - 45.0d0 * f(:,:,i-1) + 45.0d0 * f(:,:,i+1) - 9.0d0 * f(:,:,i+2) + f(:,:,i+3))
    end do
    df(:,:,N-2) = (-002.0d0 * f(:,:,N) + 024.0d0 * f(:,:,N-1) + 035.0d0 * f(:,:,N-2) &
    - 080.0d0 * f(:,:,N-3) + 030.0d0 * f(:,:,N-4) - 08.0d0 * f(:,:,N-5) + 01.0d0 * f(:,:,N-6)) * res
    df(:,:,N-1) = (+010.0d0 * f(:,:,N) + 077.0d0 * f(:,:,N-1) - 150.0d0 * f(:,:,N-2) &
    + 100.0d0 * f(:,:,N-3) - 050.0d0 * f(:,:,N-4) + 15.0d0 * f(:,:,N-5) - 02.0d0 * f(:,:,N-6)) * res
    df(:,:,N-0) = (+147.0d0 * f(:,:,N) - 360.0d0 * f(:,:,N-1) + 450.0d0 * f(:,:,N-2) &
    - 400.0d0 * f(:,:,N-3) + 225.0d0 * f(:,:,N-4) - 72.0d0 * f(:,:,N-5) + 10.0d0 * f(:,:,N-6)) * res

  end function first_derivative_z_6
  !=============================================================================!



  !=============================================================================!
  ! second derivatives XX
  !=============================================================================!
  function second_derivative_xx_2(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,1) - 1
    res = 1.0d0 / h**2

    df(0,:,:) = (f(0,:,:) - 2.0d0 * f(1,:,:) + f(2,:,:)) * res
    do i=1, N-1
      df(i,:,:) = (f(i-1,:,:) - 2.0d0 * f(i,:,:) + f(i+1,:,:)) * res
    end do
    df(N,:,:) = (f(N,:,:) - 2.0d0 * f(N-1,:,:) + f(N-2,:,:)) * res

  end function second_derivative_xx_2
  !=============================================================================!
  function second_derivative_xx_4(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,1) - 1
    res = 1.0d0 / (12.0d0 * h**2)

    df(0,:,:) = (35.0d0 * f(0,:,:) - 104.0d0 * f(1,:,:) + 114.0d0 * f(2,:,:) - &
    56.0d0 * f(3,:,:) + 11.0d0 * f(4,:,:)) * res
    df(1,:,:) = (11.0d0 * f(0,:,:) - 020.0d0 * f(1,:,:) + 006.0d0 * f(2,:,:) + &
    04.0d0 * f(3,:,:) - 01.0d0 * f(4,:,:)) * res
    do i=2, N-2
      df(i,:,:) = (- f(i-2,:,:) + 16.0d0 * f(i-1,:,:) - 30.0d0 * f(i,:,:) + 16.0d0 * f(i+1,:,:) - f(i+2,:,:)) * res
    end do
    df(N-1,:,:) = (11.0d0 * f(N,:,:) - 020.0d0 * f(N-1,:,:) + 006.0d0 * f(N-2,:,:) &
    + 04.0d0 * f(N-3,:,:) - 01.0d0 * f(N-4,:,:)) * res
    df(N-0,:,:) = (35.0d0 * f(N,:,:) - 104.0d0 * f(N-1,:,:) + 114.0d0 * f(N-2,:,:) &
    - 56.0d0 * f(N-3,:,:) + 11.0d0 * f(N-4,:,:)) * res

  end function second_derivative_xx_4
  !=============================================================================!
  function second_derivative_xx_6(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,1) - 1
    res = 1.0d0 / (180.0d0 * h**2)

    df(0,:,:) = (+812.0d0 * f(0,:,:) - 3132.0d0 * f(1,:,:) + 5265.0d0 * f(2,:,:) &
    - 5080.0d0 * f(3,:,:) + 2970.0d0 * f(4,:,:) - 972.0d0 * f(5,:,:) + 137.0d0 * f(6,:,:)) * res
    df(1,:,:) = (+137.0d0 * f(0,:,:) - 0147.0d0 * f(1,:,:) - 0255.0d0 * f(2,:,:) &
    + 0470.0d0 * f(3,:,:) - 0285.0d0 * f(4,:,:) + 093.0d0 * f(5,:,:) - 013.0d0 * f(6,:,:)) * res
    df(2,:,:) = (-013.0d0 * f(0,:,:) + 0228.0d0 * f(1,:,:) - 0420.0d0 * f(2,:,:) &
    + 0200.0d0 * f(3,:,:) + 0015.0d0 * f(4,:,:) - 012.0d0 * f(5,:,:) + 002.0d0 * f(6,:,:)) * res
    do i=3, N-3
      df(i,:,:) = (2.0d0 * f(i-3,:,:) - 27.0d0 * f(i-2,:,:) + 270.0d0 * f(i-1,:,:) &
      - 490.0d0 * f(i,:,:) + 270.0d0 * f(i+1,:,:) - 27.0d0 * f(i+2,:,:) + 2.0d0 * f(i+3,:,:)) * res
    end do
    df(N-2,:,:) = (-013.0d0 * f(N,:,:) + 0228.0d0 * f(N-1,:,:) - 0420.0d0 * f(N-2,:,:) &
    + 0200.0d0 * f(N-3,:,:) + 0015.0d0 * f(N-4,:,:) - 012.0d0 * f(N-5,:,:) + 002.0d0 * f(N-6,:,:)) * res
    df(N-1,:,:) = (+137.0d0 * f(N,:,:) - 0147.0d0 * f(N-1,:,:) - 0255.0d0 * f(N-2,:,:) &
    + 0470.0d0 * f(N-3,:,:) - 0285.0d0 * f(N-4,:,:) + 093.0d0 * f(N-5,:,:) - 013.0d0 * f(N-6,:,:)) * res
    df(N-0,:,:) = (+812.0d0 * f(N,:,:) - 3132.0d0 * f(N-1,:,:) + 5265.0d0 * f(N-2,:,:) &
    - 5080.0d0 * f(N-3,:,:) + 2970.0d0 * f(N-4,:,:) - 972.0d0 * f(N-5,:,:) + 137.0d0 * f(N-6,:,:)) * res

  end function second_derivative_xx_6
  !=============================================================================!
  ! second derivatives YY
  !=============================================================================!
  function second_derivative_yy_2(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,2) - 1
    res = 1.0d0 / h**2

    df(:,0,:) = (f(:,0,:) - 2.0d0 * f(:,1,:) + f(:,2,:)) * res
    do i=1, N-1
      df(:,i,:) = (f(:,i-1,:) - 2.0d0 * f(:,i,:) + f(:,i+1,:)) * res
    end do
    df(:,N,:) = (f(:,N,:) - 2.0d0 * f(:,N-1,:) + f(:,N-2,:)) * res

  end function second_derivative_yy_2
  !=============================================================================!
  function second_derivative_yy_4(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,2) - 1
    res = 1.0d0 / (12.0d0 * h**2)

    df(:,0,:) = (35.0d0 * f(:,0,:) - 104.0d0 * f(:,1,:) + 114.0d0 * f(:,2,:) - &
    56.0d0 * f(:,3,:) + 11.0d0 * f(:,4,:)) * res
    df(:,1,:) = (11.0d0 * f(:,0,:) - 020.0d0 * f(:,1,:) + 006.0d0 * f(:,2,:) + &
    04.0d0 * f(:,3,:) - 01.0d0 * f(:,4,:)) * res
    do i=2, N-2
      df(:,i,:) = (- f(:,i-2,:) + 16.0d0 * f(:,i-1,:) - 30.0d0 * f(:,i,:) + 16.0d0 * f(:,i+1,:) - f(:,i+2,:)) * res
    end do
    df(:,N-1,:) = (11.0d0 * f(:,N,:) - 020.0d0 * f(:,N-1,:) + 006.0d0 * f(:,N-2,:) &
    + 04.0d0 * f(:,N-3,:) - 01.0d0 * f(:,N-4,:)) * res
    df(:,N-0,:) = (35.0d0 * f(:,N,:) - 104.0d0 * f(:,N-1,:) + 114.0d0 * f(:,N-2,:) &
    - 56.0d0 * f(:,N-3,:) + 11.0d0 * f(:,N-4,:)) * res

  end function second_derivative_yy_4
  !=============================================================================!
  function second_derivative_yy_6(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,2) - 1
    res = 1.0d0 / (180.0d0 * h**2)

    df(:,0,:) = (+812.0d0 * f(:,0,:) - 3132.0d0 * f(:,1,:) + 5265.0d0 * f(:,2,:) &
    - 5080.0d0 * f(:,3,:) + 2970.0d0 * f(:,4,:) - 972.0d0 * f(:,5,:) + 137.0d0 * f(:,6,:)) * res
    df(:,1,:) = (+137.0d0 * f(:,0,:) - 0147.0d0 * f(:,1,:) - 0255.0d0 * f(:,2,:) &
    + 0470.0d0 * f(:,3,:) - 0285.0d0 * f(:,4,:) + 093.0d0 * f(:,5,:) - 013.0d0 * f(:,6,:)) * res
    df(:,2,:) = (-013.0d0 * f(:,0,:) + 0228.0d0 * f(:,1,:) - 0420.0d0 * f(:,2,:) &
    + 0200.0d0 * f(:,3,:) + 0015.0d0 * f(:,4,:) - 012.0d0 * f(:,5,:) + 002.0d0 * f(:,6,:)) * res
    do i=3, N-3
      df(:,i,:) = (2.0d0 * f(:,i-3,:) - 27.0d0 * f(:,i-2,:) + 270.0d0 * f(:,i-1,:) &
      - 490.0d0 * f(:,i,:) + 270.0d0 * f(:,i+1,:) - 27.0d0 * f(:,i+2,:) + 2.0d0 * f(:,i+3,:)) * res
    end do
    df(:,N-2,:) = (-013.0d0 * f(:,N,:) + 0228.0d0 * f(:,N-1,:) - 0420.0d0 * f(:,N-2,:) &
    + 0200.0d0 * f(:,N-3,:) + 0015.0d0 * f(:,N-4,:) - 012.0d0 * f(:,N-5,:) + 002.0d0 * f(:,N-6,:)) * res
    df(:,N-1,:) = (+137.0d0 * f(:,N,:) - 0147.0d0 * f(:,N-1,:) - 0255.0d0 * f(:,N-2,:) &
    + 0470.0d0 * f(:,N-3,:) - 0285.0d0 * f(:,N-4,:) + 093.0d0 * f(:,N-5,:) - 013.0d0 * f(:,N-6,:)) * res
    df(:,N-0,:) = (+812.0d0 * f(:,N,:) - 3132.0d0 * f(:,N-1,:) + 5265.0d0 * f(:,N-2,:) &
    - 5080.0d0 * f(:,N-3,:) + 2970.0d0 * f(:,N-4,:) - 972.0d0 * f(:,N-5,:) + 137.0d0 * f(:,N-6,:)) * res

  end function second_derivative_yy_6
  !=============================================================================!
  ! second derivatives ZZ
  !=============================================================================!
  function second_derivative_zz_2(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,3) - 1
    res = 1.0d0 / h**2

    df(:,:,0) = (f(:,:,0) - 2.0d0 * f(:,:,1) + f(:,:,2)) * res
    do i=1, N-1
      df(:,:,i) = (f(:,:,i-1) - 2.0d0 * f(:,:,i) + f(:,:,i+1)) * res
    end do
    df(:,:,N) = (f(:,:,N) - 2.0d0 * f(:,:,N-1) + f(:,:,N-2)) * res

  end function second_derivative_zz_2
  !=============================================================================!
  function second_derivative_zz_4(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,3) - 1
    res = 1.0d0 / (12.0d0 * h**2)

    df(:,0,:) = (35.0d0 * f(:,0,:) - 104.0d0 * f(:,1,:) + 114.0d0 * f(:,2,:) - &
    56.0d0 * f(:,3,:) + 11.0d0 * f(:,4,:)) * res
    df(:,1,:) = (11.0d0 * f(:,0,:) - 020.0d0 * f(:,1,:) + 006.0d0 * f(:,2,:) + &
    04.0d0 * f(:,:,3) - 01.0d0 * f(:,:,4)) * res
    do i=2, N-2
      df(:,:,i) = (- f(:,:,i-2) + 16.0d0 * f(:,:,i-1) - 30.0d0 * f(:,:,i) + 16.0d0 * f(:,:,i+1) - f(:,:,i+2)) * res
    end do
    df(:,:,N-1) = (11.0d0 * f(:,:,N) - 020.0d0 * f(:,:,N-1) + 006.0d0 * f(:,:,N-2) &
    + 04.0d0 * f(:,:,N-3) - 01.0d0 * f(:,:,N-4)) * res
    df(:,:,N-0) = (35.0d0 * f(:,:,N) - 104.0d0 * f(:,:,N-1) + 114.0d0 * f(:,:,N-2) &
    - 56.0d0 * f(:,:,N-3) + 11.0d0 * f(:,:,N-4)) * res

  end function second_derivative_zz_4
  !=============================================================================!
  function second_derivative_zz_6(f, h) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), h
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
    !
    ! internal
    !
    integer :: i, N
    real(kind=8) :: res

    N = size(f,3) - 1
    res = 1.0d0 / (180.0d0 * h**2)

    df(:,:,0) = (+812.0d0 * f(:,:,0) - 3132.0d0 * f(:,:,1) + 5265.0d0 * f(:,:,2) &
    - 5080.0d0 * f(:,:,3) + 2970.0d0 * f(:,:,4) - 972.0d0 * f(:,:,5) + 137.0d0 * f(:,:,6)) * res
    df(:,:,1) = (+137.0d0 * f(:,:,0) - 0147.0d0 * f(:,:,1) - 0255.0d0 * f(:,:,2) &
    + 0470.0d0 * f(:,:,3) - 0285.0d0 * f(:,:,4) + 093.0d0 * f(:,:,5) - 013.0d0 * f(:,:,6)) * res
    df(:,:,2) = (-013.0d0 * f(:,:,0) + 0228.0d0 * f(:,:,1) - 0420.0d0 * f(:,:,2) &
    + 0200.0d0 * f(:,:,3) + 0015.0d0 * f(:,:,4) - 012.0d0 * f(:,:,5) + 002.0d0 * f(:,:,6)) * res
    do i=3, N-3
      df(:,:,i) = (2.0d0 * f(:,:,i-3) - 27.0d0 * f(:,:,i-2) + 270.0d0 * f(:,:,i-1) &
      - 490.0d0 * f(:,:,i) + 270.0d0 * f(:,:,i+1) - 27.0d0 * f(:,:,i+2) + 2.0d0 * f(:,:,i+3)) * res
    end do
    df(:,:,N-2) = (-013.0d0 * f(:,:,N) + 0228.0d0 * f(:,:,N-1) - 0420.0d0 * f(:,:,N-2) &
    + 0200.0d0 * f(:,:,N-3) + 0015.0d0 * f(:,:,N-4) - 012.0d0 * f(:,:,N-5) + 002.0d0 * f(:,:,N-6)) * res
    df(:,:,N-1) = (+137.0d0 * f(:,:,N) - 0147.0d0 * f(:,:,N-1) - 0255.0d0 * f(:,:,N-2) &
    + 0470.0d0 * f(:,:,N-3) - 0285.0d0 * f(:,:,N-4) + 093.0d0 * f(:,:,N-5) - 013.0d0 * f(:,:,N-6)) * res
    df(:,:,N-0) = (+812.0d0 * f(:,:,N) - 3132.0d0 * f(:,:,N-1) + 5265.0d0 * f(:,:,N-2) &
    - 5080.0d0 * f(:,:,N-3) + 2970.0d0 * f(:,:,N-4) - 972.0d0 * f(:,:,N-5) + 137.0d0 * f(:,:,N-6)) * res

  end function second_derivative_zz_6
  !=============================================================================!


  !=============================================================================!
  ! Laplacian
  !=============================================================================!
  function Laplacian_2(f, dx, dy, dz) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), dx, dy, dz
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df

    df = &
    second_derivative_xx_2(f, dx) + &
    second_derivative_yy_2(f, dy) + &
    second_derivative_zz_2(f, dz)

  end function Laplacian_2
  !=============================================================================!
  function Laplacian_4(f, dx, dy, dz) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), dx, dy, dz
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df

    df = &
    second_derivative_xx_4(f, dx) + &
    second_derivative_yy_4(f, dy) + &
    second_derivative_zz_4(f, dz)

  end function Laplacian_4
  !=============================================================================!
  function Laplacian_6(f, dx, dy, dz) result(df)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: f(0:,0:,0:), dx, dy, dz
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df, dxxf, dyyf, dzzf

    dxxf = second_derivative_xx_6(f = f, h = dx)
    dyyf = second_derivative_yy_6(f = f, h = dy)
    dzzf = second_derivative_zz_6(f = f, h = dz)

    df = dxxf + dyyf + dzzf


  end function Laplacian_6
  !=============================================================================!

end module finite_differences

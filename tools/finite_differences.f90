module finite_differences
  implicit none

  interface first_derivative_x_2
    module procedure first_derivative_x_2_uniform, first_derivative_x_2_nonuniform
  end interface

  interface advec_x
    module procedure advec_x_uniform, advec_x_nonuniform
  end interface

contains

  !=============================================================================!
  ! first derivatives X
  !=============================================================================!
  function first_derivative_x_2_uniform(f, h, bounds) result(df)
    implicit none
    !
    ! input
    !
    integer, intent(in) :: bounds(3,2)
    complex(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))
    real(kind=8) :: h
    !
    ! output
    !
    complex(kind=8), dimension(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) :: df
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: res


    res = 0.50d0 / h

    df(bounds(1,1),:,:) = (-3.0d0 * f(bounds(1,1),:,:) + 4.0d0 * f(bounds(1,1) + 1,:,:) - f(bounds(1,1)+2,:,:)) * res
    do i=bounds(1,1)+1, bounds(1,2)-1
      df(i,:,:) = (-f(i-1,:,:) + f(i+1,:,:)) * res
    end do
    df(bounds(1,2),:,:) = (+3.0d0 * f(bounds(1,2),:,:) - 4.0d0 * f(bounds(1,2)-1,:,:) + f(bounds(1,2)-2,:,:)) * res

  end function first_derivative_x_2_uniform


  !=============================================================================!
  ! function first_derivative_x_2_nonuniform(f, x, bounds) result(df)
  !   implicit none
  !   !
  !   ! input
  !   !
  !   integer, intent(in) :: bounds(3,2)
  !   real(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))
  !   real(kind=8), intent(in) :: x(bounds(1,1):bounds(1,2)) 
  !   !
  !   ! output
  !   !
  !   real(kind=8), dimension(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) :: df
  !   !
  !   ! internal
  !   !
  !   integer :: i
  !   real(kind=8) :: h_m, h_p, h1, h2

  !   h1 = x(bounds(1,1)+1) - x(bounds(1,1))
  !   h2 = x(bounds(1,1)+2) - x(bounds(1,1))

  !   df(bounds(1,1),:,:) = (f(bounds(1,1)+2,:,:)*h1**2 + f(bounds(1,1),:,:)*(h2**2-h1**2) &
  !                         - f(bounds(1,1)+1,:,:)*h2**2)/(h1*h2*(h2-h1))

  !   do i=bounds(1,1)+1, bounds(1,2)-1
  !     h_p = x(i+1) - x(i)
  !     h_m = x(i) - x(i-1)
  !     df(i,:,:) = (-h_p**2*f(i-1,:,:) + (h_p**2 - h_m**2)*f(i,:,:) + h_m**2*f(i+1,:,:))/(h_m*h_p*(h_m+h_p))
  !   end do

  !   h1 = x(bounds(1,2)) - x(bounds(1,2)-1)
  !   h2 = x(bounds(1,2)) - x(bounds(1,2)-2)
  !   df(bounds(1,2),:,:) = (f(bounds(1,2)-2,:,:)*h1**2 + f(bounds(1,2),:,:)*(h2**2-h1**2) &
  !                         - f(bounds(1,2)-1,:,:)*h2**2)/(h1*h2*(h2-h1))

  ! end function first_derivative_x_2_nonuniform
  function first_derivative_x_2_nonuniform(f, x, bounds) result(df)
    implicit none
    !
    ! input
    !
    integer, intent(in) :: bounds(3,2)
    complex(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2), bounds(2,1):bounds(2,2), bounds(3,1):bounds(3,2))
    real(kind=8), intent(in) :: x(bounds(1,1):bounds(1,2))
    !
    ! output
    !
    complex(kind=8), dimension(bounds(1,1):bounds(1,2), bounds(2,1):bounds(2,2), bounds(3,1):bounds(3,2)) :: df
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: h0, h1, denom
  
    ! Forward difference (i = i_min)
    h0 = x(bounds(1,1)+1) - x(bounds(1,1))
    h1 = x(bounds(1,1)+2) - x(bounds(1,1)+1)
    denom = h0 * (h0 + h1)
    df(bounds(1,1),:,:) = (- (2.0d0*h0 + h1)*f(bounds(1,1),:,:) + (h0 + h1)*f(bounds(1,1)+1,:,:) - h0*f(bounds(1,1)+2,:,:)) / denom
  
    ! Central difference (interior points)
    do i = bounds(1,1)+1, bounds(1,2)-1
      h0 = x(i)   - x(i-1)
      h1 = x(i+1) - x(i)
      denom = h0 * (h0 + h1)
      df(i,:,:) = (-h1*f(i-1,:,:) + (h1 - h0)*f(i,:,:) + h0*f(i+1,:,:)) / (0.5d0*(h0 + h1)*h0*h1)
    end do
  
    ! Backward difference (i = i_max)
    h0 = x(bounds(1,2))   - x(bounds(1,2)-1)
    h1 = x(bounds(1,2)-1) - x(bounds(1,2)-2)
    denom = h0 * (h0 + h1)
    df(bounds(1,2),:,:) = ( (2.0d0*h0 + h1)*f(bounds(1,2),:,:) - (h0 + h1)*f(bounds(1,2)-1,:,:) + h0*f(bounds(1,2)-2,:,:)) / denom
  
  end function first_derivative_x_2_nonuniform
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
  function first_derivative_y_2(f, h, bounds,inpar) result(df)
    implicit none
    !
    ! input
    !
    integer, intent(in) :: bounds(3,2)
    complex(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))
    real(kind=8), intent(in) :: h
    logical, intent(in) :: inpar
    !
    ! output
    !
    complex(kind=8), dimension(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) :: df
    !
    ! internal
    !
    integer :: i
    integer :: k, dis, k_target
    real(kind=8) :: res


    res = 0.50d0 / h
    dis = (bounds(3,2)) / 2

   df(:,bounds(2,1),:) = (-3.0d0 * f(:,bounds(2,1),:) + 4.0d0 * f(:,bounds(2,1)+1,:) - f(:,bounds(2,1)+2,:)) * res
     do k = bounds(3,1) , bounds(3,2)
      k_target = mod(k + dis, bounds(3,2) - 1)  
        if(inpar) then 
          df(:,bounds(2,1),k) = (f(:,bounds(2,1)+2,k_target) + f(:,bounds(2,1)+1,k))*res
        else 
          df(:,bounds(2,1),k) = (-f(:,bounds(2,1)+1,k) + f(:,bounds(2,1)+2,k_target))*res
        end if  
     end do
    
    do i=bounds(2,1)+1, bounds(2,2)-1
      df(:,i,:) = (-f(:,i-1,:) + f(:,i+1,:)) * res
    end do
   df(:,bounds(2,2),:) = (+3.0d0 * f(:,bounds(2,2),:) - 4.0d0 * f(:,bounds(2,2)-1,:) + f(:,bounds(2,2)-2,:)) * res
    ! do k = bounds(3,1) , bounds(3,2)
    !   k_target = mod(k + dis, bounds(3,2) - 1) 
    !     if(inpar) then  
    !       df(:,bounds(2,2),k) = (-f(:,bounds(2,2)-1,k) - f(:,bounds(2,2)-2,k_target))*res
    !     else 
    !       df(:,bounds(2,2),k) = (-f(:,bounds(2,2)-2,k_target) + f(:,bounds(2,2)-1,k))*res
    !     end if
    ! end do
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
  function first_derivative_z_2(f, h, bounds) result(df)
    implicit none
    !
    ! input
    !
    integer, intent(in) :: bounds(3,2)
    complex(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))
    real(kind=8), intent(in) :: h
    !
    ! output
    !
    complex(kind=8), dimension(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) :: df
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: res

    res = 0.50d0 / h

   ! df(:,:,bounds(3,1)) = (-3.0d0 * f(:,:,bounds(3,1)) + 4.0d0 * f(:,:,bounds(3,1)+1) - f(:,:,bounds(3,1)+2)) * res
    df(:,:,bounds(3,1)) = (-f(:,:,bounds(3,2)-2)+f(:,:,bounds(3,1)+1))*res
    do i=bounds(3,1)+1, bounds(3,2)-1
      df(:,:,i) = (-f(:,:,i-1) + f(:,:,i+1)) * res
    end do
    df(:,:,bounds(3,2)) = (-f(:,:,bounds(3,1)+2)+f(:,:,bounds(3,2)-1))*res
    !df(:,:,bounds(3,2)) = (+3.0d0 * f(:,:,bounds(3,2)) - 4.0d0 * f(:,:,bounds(3,2)-1) + f(:,:,bounds(3,2)-2)) * res

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
  function advec_x_uniform(f, beta, dx, bounds) result(df)
    implicit none
    ! Input
    integer, intent(in) :: bounds(3,2)
    complex(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2))
    real(kind=8), intent(in) :: beta(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)), dx
    ! Output
    complex(kind=8), dimension(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) :: df
    ! Internal
    real(kind=8) :: res
    integer :: i, j, k

    res = 0.5d0 / dx

    df(bounds(1,1),:,:) = (-3.0d0 * f(bounds(1,1),:,:) + 4.0d0 * f(bounds(1,1)+1,:,:) - f(bounds(1,1)+2,:,:)) * res
    df(bounds(1,1)+1,:,:) = (-3.0d0 * f(bounds(1,1)+1,:,:) + 4.0d0 * f(bounds(1,1)+2,:,:) - f(bounds(1,1)+3,:,:)) * res

    do k = bounds(3,1), bounds(3,2)
      do j = bounds(2,1), bounds(2,2)
        do i = bounds(1,1)+2, bounds(1,2)-2
          if (beta(i,j,k) > 0.0d0) then
            df(i,j,k) = (-3.0d0 * f(i,j,k) + 4.0d0 * f(i+1,j,k) - f(i+2,j,k)) * res
          else if (beta(i,j,k) < 0.0d0) then
            df(i,j,k) = (3.0d0 * f(i,j,k) - 4.0d0 * f(i-1,j,k) + f(i-2,j,k)) * res 
          else
            df(i,j,k) = (f(i+1,j,k) - f(i-1,j,k)) * res
          end if
        end do
      end do
    end do

    df(bounds(1,2)-1,:,:) = (f(bounds(1,2),:,:)-f(bounds(1,2)-2,:,:))*res
    df(bounds(1,2),:,:) = (3.0d0 * f(bounds(1,2),:,:) - 4.0d0 * f(bounds(1,2)-1,:,:) + f(bounds(1,2)-2,:,:)) * res

  end function advec_x_uniform
 !=============================================================================!
  ! function advec_x_nonuniform(f, beta, x, bounds) result(df)
  !   implicit none
  !   ! Input
  !   integer, intent(in) :: bounds(3,2)
  !   real(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)), &
  !                               beta(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)), &
  !                               x(bounds(1,1):bounds(1,2))

  !   ! Output
  !   real(kind=8), dimension(bounds(1,1):bounds(1,2),bounds(2,1):bounds(2,2),bounds(3,1):bounds(3,2)) :: df
  !   ! Internal
  !   integer :: i, j, k
  !   real(kind=8) :: h_m, h_p, h1, h2

  !   h1 = x(bounds(1,1)+1) - x(bounds(1,1))
  !   h2 = x(bounds(1,1)+2) - x(bounds(1,1))

  !   df(bounds(1,1),:,:) = (f(bounds(1,1)+2,:,:)*h1**2 + f(bounds(1,1),:,:)*(h2**2-h1**2) &
  !                         - f(bounds(1,1)+1,:,:)*h2**2)/(h1*h2*(h2-h1))

  !   h1 = x(bounds(1,1)+2) - x(bounds(1,1)+1)
  !   h2 = x(bounds(1,1)+3) - x(bounds(1,1)+1)

  !   df(bounds(1,1)+1,:,:) = (f(bounds(1,1)+3,:,:)*h1**2 + f(bounds(1,1)+1,:,:)*(h2**2-h1**2) &
  !                           - f(bounds(1,1)+2,:,:)*h2**2)/(h1*h2*(h2-h1))

  !   do k = bounds(3,1), bounds(3,2)
  !     do j = bounds(2,1), bounds(2,2)
  !       do i = bounds(1,1)+2, bounds(1,2)-2
  !         if (beta(i,j,k) > 0.0d0) then
  !           h1 = x(i+1) - x(i)
  !           h2 = x(i+2) - x(i)
  !           df(i,j,k) = (f(i+2,j,k)*h1**2 + f(i,j,k)*(h2**2-h1**2) - f(i+1,j,k)*h2**2)/(h1*h2*(h2-h1))
  !         else if (beta(i,j,k) < 0.0d0) then
  !           h1 = x(i) - x(i-1)
  !           h2 = x(i) - x(i-2)
  !           df(i,j,k) = (f(i-2,j,k)*h1**2 + f(i,j,k)*(h2**2-h1**2) - f(i-1,j,k)*h2**2)/(h1*h2*(h2-h1))        
  !         else
  !           h_p = x(i+1) - x(i)
  !           h_m = x(i) - x(i-1)
  !           df(i,j,k) = (-h_p**2*f(i-1,j,k) + (h_p**2 - h_m**2)*f(i,j,k) + h_m**2*f(i+1,j,k))/(h_m*h_p*(h_m+h_p))
  !         end if
  !       end do
  !     end do
  !   end do

  !   h_p = x(bounds(1,2)) - x(bounds(1,2)-1)
  !   h_m = x(bounds(1,2)-1) - x(bounds(1,2)-2)
  !   df(bounds(1,2)-1,:,:) = (-h_p**2*f(bounds(1,2)-2,:,:) + (h_p**2 - h_m**2)*f(bounds(1,2)-1,:,:) & 
  !                         + h_m**2*f(bounds(1,2),:,:))/(h_m*h_p*(h_m+h_p))
  !   h1 = x(bounds(1,2)) - x(bounds(1,2)-1)
  !   h2 = x(bounds(1,2)) - x(bounds(1,2)-2)
  !   df(bounds(1,2),:,:) = (f(bounds(1,2)-2,:,:)*h1**2 + f(bounds(1,2),:,:)*(h2**2-h1**2) &
  !                         - f(bounds(1,2)-1,:,:)*h2**2)/(h1*h2*(h2-h1))
  ! end function advec_x_nonuniform

  function advec_x_nonuniform(f, beta, x, bounds) result(df)
    implicit none
    ! Input
    integer, intent(in) :: bounds(3,2)
    complex(kind=8), intent(in) :: f(bounds(1,1):bounds(1,2), bounds(2,1):bounds(2,2), bounds(3,1):bounds(3,2))
    real(kind=8), intent(in) :: beta(bounds(1,1):bounds(1,2), bounds(2,1):bounds(2,2), bounds(3,1):bounds(3,2))
    real(kind=8), intent(in) :: x(bounds(1,1):bounds(1,2))
    ! Output
    complex(kind=8), dimension(bounds(1,1):bounds(1,2), bounds(2,1):bounds(2,2), bounds(3,1):bounds(3,2)) :: df
    ! Internal
    integer :: i, j, k
    real(kind=8) :: h0, h1, denom
  
    ! Frontera izquierda (i = i_min)
    h0 = x(bounds(1,1)+1) - x(bounds(1,1))
    h1 = x(bounds(1,1)+2) - x(bounds(1,1)+1)
    denom = h0 * (h0 + h1)
    df(bounds(1,1),:,:) = (- (2.0d0*h0 + h1)*f(bounds(1,1),:,:) &
                          + (h0 + h1)*f(bounds(1,1)+1,:,:) - h0*f(bounds(1,1)+2,:,:)) / denom
  
    ! Segundo punto (usamos también forward, ya que el esquema necesita dos antes del centro)
    h0 = x(bounds(1,1)+2) - x(bounds(1,1)+1)
    h1 = x(bounds(1,1)+3) - x(bounds(1,1)+2)
    denom = h0 * (h0 + h1)
    df(bounds(1,1)+1,:,:) = (- (2.0d0*h0 + h1)*f(bounds(1,1)+1,:,:) &
                         + (h0 + h1)*f(bounds(1,1)+2,:,:) - h0*f(bounds(1,1)+3,:,:)) / denom
  
    ! Interior
    do k = bounds(3,1), bounds(3,2)
      do j = bounds(2,1), bounds(2,2)
        do i = bounds(1,1)+2, bounds(1,2)-2
          if (beta(i,j,k) > 0.0d0) then
            h0 = x(i+1) - x(i)
            h1 = x(i+2) - x(i+1)
            denom = h0 * (h0 + h1)
            df(i,j,k) = (- (2.0d0*h0 + h1)*f(i,j,k) + (h0 + h1)*f(i+1,j,k) - h0*f(i+2,j,k)) / denom
          else if (beta(i,j,k) < 0.0d0) then
            h0 = x(i)   - x(i-1)
            h1 = x(i-1) - x(i-2)
            denom = h0 * (h0 + h1)
            df(i,j,k) = ((2.0d0*h0 + h1)*f(i,j,k) - (h0 + h1)*f(i-1,j,k) + h0*f(i-2,j,k)) / denom
          else
            h0 = x(i)   - x(i-1)
            h1 = x(i+1) - x(i)
            df(i,j,k) = (-h1*f(i-1,j,k) + (h1 - h0)*f(i,j,k) + h0*f(i+1,j,k)) / (0.5d0*(h0 + h1)*h0*h1)
          end if
        end do
      end do
    end do
  
    ! Penúltimo punto (central simple entre f_N y f_{N-2})
    h0 = x(bounds(1,2)-1) - x(bounds(1,2)-2)
    h1 = x(bounds(1,2)) - x(bounds(1,2)-1)
    df(bounds(1,2)-1,:,:) = (f(bounds(1,2),:,:) - f(bounds(1,2)-2,:,:)) / (h0 + h1)
  
    ! Último punto (i = i_max)
    h0 = x(bounds(1,2)) - x(bounds(1,2)-1)
    h1 = x(bounds(1,2)-1) - x(bounds(1,2)-2)
    denom = h0 * (h0 + h1)
    df(bounds(1,2),:,:) = ((2.0d0*h0 + h1)*f(bounds(1,2),:,:) - (h0 + h1)*f(bounds(1,2)-1,:,:) + h0*f(bounds(1,2)-2,:,:)) / denom
  
  end function advec_x_nonuniform  
 !=============================================================================!
function advec_x_4(f, beta,dx) result(df)
  implicit none
  ! Input
  real(kind=8), intent(in) :: f(0:,0:,0:), beta(0:,0:,0:), dx
  ! Output
  real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
  ! Internal
  real(kind=8) :: res
  integer :: i, j, k, N

  N = size(f,1)-1

  res = 1.0d0/dx

  df(0,:,:) = (-25.0d0 * f(0,:,:) + 48.0d0 * f(1,:,:) - 36.0d0 * f(2,:,:) + 16.0d0 * f(3,:,:) - 3.0d0 * f(4,:,:)) * res / 12.0d0
  df(1,:,:) = (-3.0d0 * f(0,:,:) - 10.0d0 * f(1,:,:) + 18.0d0 * f(2,:,:) - 6.0d0 * f(3,:,:) + f(4,:,:)) * res / 12.0d0
  df(2,:,:) = (f(0,:,:) - 8.0d0 * f(1,:,:) + 8.0d0 * f(3,:,:) - f(4,:,:)) * res / 12.0d0
  do j=0, size(f,2)-1
  do k=0, size(f,3)-1
  do i=3, N-3
    if(beta(i,j,k) < 0.0d0) then
      df(i,j,k) = (3.0d0 * f(i+1,j,k) + 10.0d0 * f(i,j,k) - 18.0d0 * f(i-1,j,k) + 6.0d0 * f(i-2,j,k) - f(i-3,j,k)) * res / 12.0d0
    else if(beta(i,j,k) > 0.0d0) then
      df(i,j,k) = (-3.0d0 * f(i-1,j,k) - 10.0d0 * f(i,j,k) + 18.0d0 * f(i+1,j,k) - 6.0d0 * f(i+2,j,k) + f(i+3,j,k)) * res / 12.0d0
    else
      df(i,j,k) = (f(i-2,j,k) - 8.0d0 * f(i-1,j,k) + 8.0d0 * f(i+1,j,k) - f(i+2,j,k)) * res / 12.0d0
    end if
  end do
  end do
  end do
  df(N-2,:,:) = (f(N-4,:,:) - 8.0d0 * f(N-3,:,:) + 8.0d0 * f(N-1,:,:) - f(N,:,:)) * res / 12.0d0
  df(N-1,:,:) = (3.0d0 * f(N,:,:) + 10.0d0 * f(N-1,:,:) - 18.0d0 * f(N-2,:,:) + 6.0d0 * f(N-3,:,:) - f(N-4,:,:)) * res / 12.0d0
  df(N  ,:,:) = (25.0d0 * f(N,:,:) - 48.0d0 * f(N-1,:,:) + 36.0d0 * f(N-2,:,:) &
              - 16.0d0 * f(N-3,:,:) + 3.0d0 * f(N-4,:,:)) * res / 12.0d0

end function advec_x_4

!=============================================================================!


function advec_y(f, beta, dy) result(df)
  implicit none
  ! Input
  real(kind=8), intent(in) :: f(0:,0:,0:), beta(0:,0:,0:), dy
  ! Output
  real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
  ! Internal
  real(kind=8) :: res
  integer :: i, j, k, N

  res = 0.5d0 / dy
  N = size(f,2) - 1

  df(:,0,:) = (-3.0d0 * f(:,0,:) + 4.0d0 * f(:,1,:) - f(:,2,:)) * res
  df(:,N,:) = (3.0d0 * f(:,N,:) - 4.0d0 * f(:,N-1,:) + f(:,N-2,:)) * res
  df(:,1,:) = (f(:,2,:)-f(:,0,:))*res
  df(:,N-1,:) = (f(:,N,:)-f(:,N-2,:))*res
  do k = 0, size(f,3)-1
    do i = 0, size(f,1)-1
      do j = 2, N-2
        if (beta(i,j,k) > 0.0d0) then
          df(i,j,k) = (-3.0d0 * f(i,j,k) + 4.0d0 * f(i,j+1,k) - f(i,j+2,k)) * res
        else if (beta(i,j,k) < 0.0d0) then
          df(i,j,k) =  (3.0d0 * f(i,j,k) - 4.0d0 * f(i,j-1,k) + f(i,j-2,k)) * res
        else
          df(i,j,k) = (f(i,j+1,k) - f(i,j-1,k)) * res
        end if
      end do
    end do
  end do

end function advec_y
!=============================================================================!

function advec_z(f, beta, dz) result(df)
  implicit none
  ! Input
  real(kind=8), intent(in) :: f(0:,0:,0:), beta(0:,0:,0:), dz
  ! Output
  real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: df
  ! Internal
  real(kind=8) :: res
  integer :: i, j, k, N

  res = 0.5d0 / dz
  N = size(f,3) - 1

  df(:,:,0) = (-3.0d0 * f(:,:,0) + 4.0d0 * f(:,:,1) - f(:,:,2)) * res
  df(:,:,N) = (3.0d0 * f(:,:,N) - 4.0d0 * f(:,:,N-1) + f(:,:,N-2)) * res
  df(:,:,1) = (f(:,:,2)-f(:,:,0))*res
  df(:,:,N-1) = (f(:,:,N)-f(:,:,N-2))*res
  
  do j = 0, size(f,2)-1
    do i = 0, size(f,1)-1
      do k = 2, N-2
        if (beta(i,j,k) > 0.0d0) then
          df(i,j,k) = (-3.0d0 * f(i,j,k) + 4.0d0 * f(i,j,k+1) - f(i,j,k+2)) * res
        else if (beta(i,j,k) < 0.0d0) then
          df(i,j,k) = (3.0d0 * f(i,j,k) - 4.0d0 * f(i,j,k-1) + f(i,j,k-2)) * res
        else
          df(i,j,k) = (f(i,j,k+1) - f(i,j,k-1)) * res
        end if
      end do
    end do
  end do

end function advec_z
!=============================================================================!


end module finite_differences

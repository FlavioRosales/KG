module ODE
  implicit none

contains
  !=============================================================================!
  function linspace(min, max, npoints) result(this)
    implicit none
    integer, intent(in) :: npoints
    real(kind=8), intent(in) :: min, max
    real(kind=8), dimension(0:npoints) :: this
    integer :: i
    real(kind=8) :: h

    h = (max - min) / dble(npoints)

    do i=0, npoints
      this(i) = min + dble(i) * h
    end do

  end function linspace
  !=============================================================================!
  function solve(initial_conditions, this, rhs) result(solution)
    implicit none
    real(kind=8), dimension(:) :: initial_conditions, this
    interface
      function rhs(t, x) result(dotx)
        real(kind=8), intent(in) :: t, x(:)
        real(kind=8), dimension(size(x)) :: dotx
      end function rhs
    end interface
    real(kind=8), dimension(size(initial_conditions), size(this)) :: solution
    integer :: j
    real(kind=8) :: h
    real(kind=8), dimension(size(initial_conditions)) :: k1, k2, k3, k4

    solution(:,1) = initial_conditions

    do j=2, size(this)

      h = this(j) - this(j-1)

      k1 = rhs(this(j-1), solution(:,j-1)) * h
      k2 = rhs(this(j-1) + 0.50d0 * h, solution(:,j-1) + 0.50d0 * k1) * h
      k3 = rhs(this(j-1) + 0.50d0 * h, solution(:,j-1) + 0.50d0 * k2) * h
      k4 = rhs(this(j-1) + h, solution(:,j-1) + k3) * h

      solution(:,j) = solution(:,j-1) + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0

    end do

  end function solve
  !=============================================================================!
  real(kind=8) function linear_interpolate(f, this, point)
    implicit none
    real(kind=8), intent(in) :: f(:), this(:), point
    integer i

    do i=1, size(f)-1

      if(this(i).le.point .and. point.le.this(i+1)) then
        linear_interpolate = f(i) + (point - this(i)) * (f(i+1) - f(i)) / (this(i+1) - this(i))
        return
      end if

    end do

  end function linear_interpolate
  !=============================================================================!

end module ODE

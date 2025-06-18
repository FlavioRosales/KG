module integral
  implicit none

  interface trapezium_2D
    module procedure trapezium_2D_real, trapezium_2D_complex
  end interface trapezium_2D

  interface trapezium_1D
    module procedure trapezium_1D_real, trapezium_1D_complex
  end interface trapezium_1D

  interface trapezium
      module procedure trapezium_real, trapezium_complex
  end interface
  
  interface integrate
    module procedure integrate_real, integrate_complex
  end interface


  interface simpson_combined_1D
    module procedure simpson_combined_1D_real, simpson_combined_1D_complex
  end interface simpson_combined_1D

  interface simpson_2D
  module procedure simpson_2D_real, simpson_2D_complex
end interface simpson_2D

contains

  !=============================================================================!

function filon_local_corrected(f, k, dx) result(integral)
  implicit none
    ! Inputs
    complex(kind=8), intent(in) :: f(:)
    real(kind=8), intent(in) :: k, dx
    ! Output
    complex(kind=8) :: integral
    ! Internals
    integer :: i_min, i_max, N, i
    real(kind=8) :: x0, coskx0, sinkx0
    real(kind=8) :: alpha, beta, gamma
    complex(kind=8) :: sum_cos, sum_sin

    ! Domain
    i_min = lbound(f,1)
    i_max = ubound(f,1)
    N = i_max - i_min

    if (mod(N,2) /= 0) then
        print *, "Error: Número de subintervalos debe ser par en Filón local."
        stop
    end if

    x0 = 0.0d0
    coskx0 = 1.0d0
    sinkx0 = 0.0d0

    ! Coeficientes clásicos de Filón
    alpha = 2.0d0 * sin(k*dx) / (k*dx)
    beta  = (2.0d0 / (k*dx)) * sin(k*dx/2.0d0)**2
    gamma = (2.0d0 / (k*dx)) * (sin(k*dx) - k*dx*cos(k*dx)) / (k*dx)

    ! Inicializar las sumas
    sum_cos = (0.0d0, 0.0d0)
    sum_sin = (0.0d0, 0.0d0)

    do i = 0, N, 2
        sum_cos = sum_cos + f(i_min+i)
    end do
    do i = 1, N-1, 2
        sum_sin = sum_sin + f(i_min+i)
    end do

    integral = dx * (alpha*sum_cos + beta*(f(i_min) + f(i_max)) + gamma*sum_sin) / 2.0d0

  end function filon_local_corrected


  function filon_or_trapezium_1D_complex(f, x, omega) result(integral)
    implicit none
    complex(kind=8), intent(in) :: f(:)
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(in) :: omega
    complex(kind=8) :: integral
    real(kind=8) :: h
    real(kind=8), parameter :: omega_h_threshold = 1.0d0 ! ← umbral sobre omega*h
    integer :: i, N

    N = size(f)

    if (N < 3) then
        print *, "Error: Se requieren al menos 3 puntos."
        stop
    end if

    h = x(2) - x(1)

    ! Evaluar producto frecuencia × paso
    if (abs(omega * h) < omega_h_threshold) then
        ! Si la frecuencia "percibida" es baja, usar trapecio
        integral = trapezium_1D_complex(f, h)
    else
        ! Si la frecuencia "percibida" es alta, usar Filón
        if (mod(N-1,2) /= 0) then
            print *, "Error: Filón requiere número impar de puntos."
            stop
        end if
        integral = filon_local_corrected(f, omega, h)
    end if

end function filon_or_trapezium_1D_complex

  !=============================================================================!

  function filon_local_exp(f0, f1, f2, theta, h, x0) result(I)
    implicit none
    complex(kind=8), parameter :: i_unreal = (0.0d0,1.0d0)
    complex(kind=8), intent(in) :: f0, f1, f2
    real(kind=8), intent(in) :: theta, h, x0
    complex(kind=8) :: I
    complex(kind=8) :: E0, E1, E2
    complex(kind=8) :: A, B, C
    real(kind=8) :: sin_theta, cos_theta

    ! Fases evaluadas en los tres puntos
    E0 = exp(i_unreal * (0.0d0) * theta) ! exp(i omega x0)
    E1 = exp(i_unreal * (1.0d0) * theta) ! exp(i omega (x0+h))
    E2 = exp(i_unreal * (2.0d0) * theta) ! exp(i omega (x0+2h))

    ! Coeficientes
    sin_theta = dsin(theta)
    cos_theta = dcos(theta)

    A = (2.0d0 * (cos_theta - 1.0d0)) / (theta**2)
    B = (4.0d0 * sin_theta) / (theta)
    C = (2.0d0 * (1.0d0 - cos_theta)) / (theta**2)

    ! Evaluar la contribución local
    I = (h/2.0d0) * (A * f0 * E0 + B * f1 * E1 + C * f2 * E2)

  end function filon_local_exp
  !=============================================================================!
  function filon_exp_1D_complex(f, x, omega) result(integral)
    implicit none
    complex(kind=8), intent(in) :: f(:)
    real(kind=8), intent(in) :: x(:)
    real(kind=8), intent(in) :: omega
    complex(kind=8) :: integral
    integer :: i, N
    real(kind=8) :: h, theta
    complex(kind=8) :: sum

    N = size(f)
    if (N < 3 .or. mod(N-1,2) /= 0) then
        print *, "Error: Filon requiere un número impar de puntos >= 3."
        stop
    end if

    ! Malla uniforme
    h = x(2) - x(1)
    theta = omega * h

    sum = (0.0d0, 0.0d0)

    do i = 1, N-2, 2
        sum = sum + filon_local_exp(f(i), f(i+1), f(i+2), theta, h, x(i))
    end do

    integral = sum

  end function filon_exp_1D_complex
  !=============================================================================!

  function trapezium_1D_real(f, dx) result(integral)
  implicit none
  real(kind=8), intent(in) :: dx
  real(kind=8), intent(in), dimension(:) :: f
  real(kind=8) ::  integral
  real(kind=8) :: sum 
  integer :: i
  
  integral = 0.0d0
  sum = 0.0d0

  ! Sumar los valores de la función en los puntos interiores
  do i = lbound(f,1)+1, ubound(f,1)-1
    sum = sum + f(i)
  end do

  ! Calcular la integral usando la regla del trapecio
  integral = dx * (0.5d0 * f(lbound(f,1)) + sum + 0.5 * f(ubound(f,1)))

  end function trapezium_1D_real
  !=============================================================================!
  function trapezium_1D_complex(f, dx) result(integral)
    implicit none
    real(kind=8), intent(in) :: dx
    complex(kind=8), intent(in), dimension(:) :: f
    complex(kind=8) ::  integral
    complex(kind=8) :: sum 
    integer :: i
    
    integral = 0.0d0
    sum = 0.0d0
  
    ! Sumar los valores de la función en los puntos interiores

    integral = (0.5d0 * (f(1) + f(size(f))))  ! extremos con peso 1/2

    do i = 2, size(f) - 1
        integral = integral + f(i)
    end do

    integral = dx * integral

    end function trapezium_1D_complex
  !=============================================================================!
  function trapezium_real(f3d, dx, dy, dz) result(f0d)
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

  end function trapezium_real
  !=============================================================================!
  function trapezium_complex(f3d, dx, dy, dz) result(f0d)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(:,:,:), intent(in) :: f3d(:,:,:)
    real(kind=8) :: dx, dy, dz
    !
    ! output
    !
    complex(kind=8) :: f0d
    !
    ! internal
    !
    integer :: i
    complex(kind=8) :: f2d(size(f3d,1), size(f3d,2)), f1d(size(f3d,1))

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

  end function trapezium_complex
  !=============================================================================!
  function trapezium_2D_real(f, dx, dy) result(f0d)
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

  end function trapezium_2D_real
  !=============================================================================!
  function trapezium_2D_complex(f, dx, dy) result(integral)
    implicit none
    complex(kind=8), intent(in) :: f(:,:)
    real(kind=8), intent(in) :: dx, dy
    complex(kind=8) :: integral
    integer :: i, j
    integer :: i_min, i_max, j_min, j_max

    i_min = lbound(f,1)
    i_max = ubound(f,1)
    j_min = lbound(f,2)
    j_max = ubound(f,2)

    integral = (0.0d0, 0.0d0)

    do i = i_min, i_max
        do j = j_min, j_max
            select case (merge(1,0,i==i_min .or. i==i_max) + 2*merge(1,0,j==j_min .or. j==j_max))
            case (0)
                ! Punto interior
                integral = integral + f(i,j)
            case (1,2)
                ! Borde (no esquina)
                integral = integral + 0.5d0 * f(i,j)
            case (3)
                ! Esquina
                integral = integral + 0.25d0 * f(i,j)
            end select
        end do
    end do

    integral = dx * dy * integral
end function trapezium_2D_complex

  !=============================================================================!


  function simpson_2D_complex(f, dx, dy) result(integral)
    implicit none
    complex(kind=8), intent(in) :: f(:,:)
    real(kind=8), intent(in) :: dx, dy
    complex(kind=8) :: integral
    integer :: i, j
    integer :: i_min, i_max, j_min, j_max, Nx, Ny
    complex(kind=8) :: sum

    i_min = lbound(f,1)
    i_max = ubound(f,1)
    j_min = lbound(f,2)
    j_max = ubound(f,2)

    Nx = i_max - i_min + 1
    Ny = j_max - j_min + 1

    ! Checar que el número de puntos sea impar en ambas direcciones
    if (mod(Nx,2) == 0 .or. mod(Ny,2) == 0) then
        print *, "Error: Simpson 1/3 requiere número impar de puntos en ambas direcciones."
        print *, "Se recibió Nx =", Nx, ", Ny =", Ny
        stop
    end if

    sum = (0.0d0, 0.0d0)
    integral = (0.0d0, 0.0d0)

    ! Aplicar regla de Simpson 1/3 compuesta en 2D
    do i = i_min, i_max
        do j = j_min, j_max
            if ((i == i_min .or. i == i_max) .and. (j == j_min .or. j == j_max)) then
                ! Esquinas
                sum = sum + f(i,j)
            else if ((i == i_min .or. i == i_max) .or. (j == j_min .or. j == j_max)) then
                ! Bordes
                if (mod(i+j,2) == 0) then
                    sum = sum + 2.0d0 * f(i,j)
                else
                    sum = sum + 4.0d0 * f(i,j)
                end if
            else
                ! Interior
                if (mod(i+j,2) == 0) then
                    sum = sum + 4.0d0 * f(i,j)
                else
                    sum = sum + 8.0d0 * f(i,j)
                end if
            end if
        end do
    end do

    integral = (dx * dy / 9.0d0) * sum

  end function simpson_2D_complex

  function simpson_2D_real(f, dx, dy) result(integral)
    implicit none
    real(kind=8), intent(in) :: f(:,:)
    real(kind=8), intent(in) :: dx, dy
    real(kind=8) :: integral
    integer :: i, j
    integer :: i_min, i_max, j_min, j_max, Nx, Ny
    real(kind=8) :: sum

    i_min = lbound(f,1)
    i_max = ubound(f,1)
    j_min = lbound(f,2)
    j_max = ubound(f,2)

    Nx = i_max - i_min + 1
    Ny = j_max - j_min + 1

    ! Checar que el número de puntos sea impar en ambas direcciones
    if (mod(Nx,2) == 0 .or. mod(Ny,2) == 0) then
        print *, "Error: Simpson 1/3 requiere número impar de puntos en ambas direcciones."
        print *, "Se recibió Nx =", Nx, ", Ny =", Ny
        stop
    end if

    sum = (0.0d0, 0.0d0)
    integral = (0.0d0, 0.0d0)

    ! Aplicar regla de Simpson 1/3 compuesta en 2D
    do i = i_min, i_max
        do j = j_min, j_max
            if ((i == i_min .or. i == i_max) .and. (j == j_min .or. j == j_max)) then
                ! Esquinas
                sum = sum + f(i,j)
            else if ((i == i_min .or. i == i_max) .or. (j == j_min .or. j == j_max)) then
                ! Bordes
                if (mod(i+j,2) == 0) then
                    sum = sum + 2.0d0 * f(i,j)
                else
                    sum = sum + 4.0d0 * f(i,j)
                end if
            else
                ! Interior
                if (mod(i+j,2) == 0) then
                    sum = sum + 4.0d0 * f(i,j)
                else
                    sum = sum + 8.0d0 * f(i,j)
                end if
            end if
        end do
    end do

    integral = (dx * dy / 9.0d0) * sum

  end function simpson_2D_real


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

  function simpson_combined_1D_complex(f, dx) result(integral)
    implicit none
    complex(kind=8), intent(in), dimension(:) :: f
    real(kind=8), intent(in) :: dx
    complex(kind=8) :: integral
    integer :: i, i_min, i_max, N
    complex(kind=8) :: sum

    i_min = lbound(f,1)
    i_max = ubound(f,1)
    N = i_max - i_min + 1
    sum = (0.0d0, 0.0d0)
    integral = (0.0d0, 0.0d0)

    ! Simpson 1/3 requiere un número impar de puntos
    if (mod(N, 2) == 0) then
        print *, "Error: Simpson 1/3 requiere número impar de puntos. Se recibió N =", N
        stop
    end if

    ! Aplicar la regla compuesta de Simpson 1/3
    sum = f(i_min) + f(i_max)
    do i = i_min + 1, i_max - 1, 2
        sum = sum + 4.0d0 * f(i)
    end do
    do i = i_min + 2, i_max - 2, 2
        sum = sum + 2.0d0 * f(i)
    end do

    integral = (dx / 3.0d0) * sum
end function simpson_combined_1D_complex





    function simpson_combined_1D_real(f, dx) result(integral)
      implicit none
      real(kind=8), intent(in) :: dx
      real(kind=8), intent(in) :: f(:)
      real(kind=8) :: integral
      integer :: N, i
      integer :: N1, N2

      N = size(f)
      integral = 0.0d0

      ! Manejo general: combinar Simpson 1/3 y 3/8 si es necesario
      if (mod(N-1,2) == 0) then
          ! Caso perfecto para Simpson 1/3: número de subintervalos par
          integral = f(1) + f(N)
          do i = 2, N-1, 2
              integral = integral + 4.0d0 * f(i)
          end do
          do i = 3, N-2, 2
              integral = integral + 2.0d0 * f(i)
          end do
          integral = integral * dx / 3.0d0

      else if (mod(N-1,3) == 0) then
          ! Caso perfecto para Simpson 3/8: número de subintervalos múltiplo de 3
          integral = f(1) + f(N)
          do i = 2, N-1
              if (mod(i-1,3) == 0) then
                  integral = integral + 2.0d0 * f(i)
              else
                  integral = integral + 3.0d0 * f(i)
              end if
          end do
          integral = integral * 3.0d0 * dx / 8.0d0

      else
          ! Mezcla: usar Simpson 1/3 en la mayor parte, 3/8 al final
          ! Longitud ideal para 1/3
          N1 = N - 3  ! hasta N-3
          N2 = 4      ! últimos 4 puntos para 3/8

          ! Simpson 1/3
          integral = f(1) + f(N1)
          do i = 2, N1-1, 2
              integral = integral + 4.0d0 * f(i)
          end do
          do i = 3, N1-2, 2
              integral = integral + 2.0d0 * f(i)
          end do
          integral = integral * dx / 3.0d0

          ! Simpson 3/8 en los últimos 4 puntos
          integral = integral + (3.0d0 * dx / 8.0d0) * &
              (f(N-3) + 3.0d0*f(N-2) + 3.0d0*f(N-1) + f(N))
      end if
  end function simpson_combined_1D_real


  !=============================================================================!
  function integrate_real(f3d, dx, dy, dz, g_pts) result(f0d)
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
    kL = 1 ; kR = size(f3d,3) 

    f0d = 0.0d0
    int = trapezium(f3d = f3d(iL:iR,jL:jR,kL:kR),dx =  dx, dy = dy, dz = dz)

    call MPI_ALLREDUCE(int, f0d, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  end function integrate_real

  function integrate_complex(f3d, dx, dy, dz, g_pts) result(f0d)
    use mpi_lib
    implicit none
    !
    ! input
    !
    integer, intent(in) :: g_pts
    complex(kind=8), dimension(:,:,:), intent(in) :: f3d(:,:,:)
    real(kind=8) :: dx, dy, dz
    !
    ! output
    !
    complex(kind=8) :: f0d
    !
    ! internal
    !
    complex(kind=8) :: int
    integer :: iL, iR, jL, jR, kL, kR

    iL = 1 + g_pts; iR = size(f3d,1) - g_pts
    jL = 1 + g_pts; jR = size(f3d,2) - g_pts
    kL = 1 ; kR = size(f3d,3) 

    f0d = 0.0d0
    int = trapezium(f3d = f3d(iL:iR,jL:jR,kL:kR),dx =  dx, dy = dy, dz = dz)

    call MPI_ALLREDUCE(int, f0d, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

  end function integrate_complex
  !=============================================================================!

end module integral

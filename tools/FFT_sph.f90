module fft
    !use mpi_lib
    implicit none 

    real(kind=8), parameter  :: pi = acos(-1.0d0)
    complex(kind=8), parameter :: i_unreal = (0.0d0,1.0d0)
    

    contains

    function fft_spherical(fun, r, theta, phi, l_max) result(ft)
        implicit none 
        real(kind=8), intent(in), dimension(:) :: r, theta, phi
        complex(kind=8), intent(in), dimension(:,:,:) :: fun
        integer, intent(in) :: l_max
        complex(kind=8), allocatable, dimension(:,:,:) :: ft
        ! spherical harmonic coeficients
        complex(kind=8), allocatable, dimension(:,:,:) :: alm
        ! spherical harmonic
        complex(kind=8) :: ylm

        integer :: Nr, Ntheta, Nphi 

        integer :: l, m, i, j, k

        Nr = size(r) 
        Ntheta = size(theta) 
        Nphi = size(phi)

        allocate(alm(Nr,0:l_max,-l_max:l_max),ft(Nr,0:l_max,-l_max:l_max))

        ! Harmonic coeficients
        alm =  0.0d0
        do l = 0, l_max
            do m = -l, l
                do i = 1, Nr 
                    do j = 1, Ntheta 
                        do k = 1, Nphi
                            ! Calcular armónico esférico Y_l^m(theta, phi)
                            ylm = compute_spherical_harmonic(l, m, theta(j), phi(k))
                            alm(i, l, m) = alm(i, l, m) + fun(i, j, k) * conjg(ylm) * sin(theta(j))
                        end do
                    end do
                    alm(i, l, m) = alm(i, l, m) * (2.0d0 * pi / Nphi) * (pi / Ntheta)
                end do
            end do
        end do

        ! Radial transform

        do l=0, l_max
            do m=-l, l
                ft(:, l, m) = fft_1D(alm(:,l,m))
            end do
        end do

    end function fft_spherical


    function DFT(f) result(DFT_f)
        implicit none 
        !
        ! input
        !
        complex(kind=8), dimension(:), intent(in) :: f 
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: DFT_f
        !
        ! others
        !
        integer :: j,k, N 

        N = size(f)

        DFT_f = 0.0d0
        do k=0, N-1
            do j=0, N-1
                DFT_f(k+1) = DFT_f(k+1) +  &
                            f(j+1)*exp(2.0d0*pi*i_unreal*dble(j*k)/dble(N))
            end do
        end do

        DFT_f = DFT_f / dble(N) ! Normalización

    end function DFT

    recursive function fft_1D(f) result(fft_f)
        implicit none 
        !
        ! input
        !
        complex(kind=8), dimension(:), intent(in) :: f 
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: fft_f
        complex(kind=8), allocatable, dimension(:) ::  fft_f_even,fft_f_odd
        !
        ! others
        !
        integer :: k, N 

        N = size(f)



        if(mod(N,2).eq.0) then 
            allocate(fft_f_even(N/2),fft_f_odd(N/2))
            fft_f_even = fft_1d(f = f(1::2))
            fft_f_odd  = fft_1d(f = f(2::2))

            do k = 0, N/2 -1
                fft_f(k+1) = fft_f_even(k+1) + &
                exp(2.0d0*pi*i_unreal*dble(k)/dble(N))*fft_f_odd(k+1)

                fft_f(k+1+N/2) = fft_f_even(k+1) - &
                exp(2.0d0*pi*i_unreal*dble(k)/dble(N))*fft_f_odd(k+1)
            end do
        else
            fft_f = DFT(f)
        end if

    end function fft_1D



    function compute_spherical_harmonic(l,m,theta,phi) result(ylm)
        implicit none 
        real(kind=8), intent(in) :: theta, phi
        integer, intent(in) :: l, m
        complex(kind=8) :: ylm

        ylm = dsqrt((dble(factorial(2*l+1))*dble(factorial(l-m)))/&
                    (4.0d0*pi*dble(factorial(l+m)))) &
              *exp(i_unreal*dble(m)*phi)

        ylm = ylm * legendre_associated(l,m,dcos(theta))


    end function compute_spherical_harmonic


    function legendre_associated(l,m,x) result(p_lm)
            ! Calcula el polinomio de Legendre asociado P_l^m(x)
            ! l: grado (entero, entrada)
            ! m: orden (entero, entrada)
            ! x: argumento del polinomio en [-1, 1] (real(kind=8), entrada)
            ! p_lm: valor del polinomio asociado (real(kind=8), salida)
        
        implicit none
        
        ! Entradas
        integer, intent(in) :: l, m
        real(kind=8), intent(in) :: x

        ! Variables locales
        integer :: i
        real(kind=8) :: p0, p1, temp, factor

        real(kind=8) :: p_lm

        ! Inicialización
        p0 = 1.0d0  ! P_0^0(x)
        p1 = x      ! P_1^0(x)

        ! Casos base para m=0
        if (m == 0) then
            if (l == 0) then
                p_lm = 1.0d0   ! P_0^0(x) = 1
                return
            elseif (l == 1) then
                p_lm = x       ! P_1^0(x) = x
                return
            else
                ! Calcular P_l^0(x) usando recurrencia
                do i = 2, l
                    temp = p1
                    p1 = ((2.0d0 * dble(i) - 1.0d0) * x * p1 - (dble(i) - 1.0d0) * p0) / dble(i)
                    p0 = temp
                end do
                p_lm = p1
                return
            end if
        end if


        ! Casos para m > 0
        if (m > 0) then
            ! Calcular el término inicial P_m^m(x)
            p_lm = 1.0d0
            factor = sqrt(1.0d0 - x**2)
            do i = 1, m
                p_lm = -p_lm * (2.0d0 * dble(i) - 1.0d0) * factor
            end do

            ! Si l = m, retornar directamente
            if (l == m) return

            ! Calcular P_{m+1}^m(x) usando recurrencia
            p0 = p_lm
            p1 = x * (2.0d0 * m + 1.0d0) * p_lm
            if (l == m + 1) then
                p_lm = p1
                return
            end if

            ! Continuar la recurrencia para l > m + 1
            do i = m + 2, l
                temp = p1
                p1 = ((2.0d0 * dble(i) - 1.0d0) * x * p1 - (dble(i) + m - 1.0d0) * p0) / (dble(i) - m)
                p0 = temp
            end do
            p_lm = p1
        end if
    end function legendre_associated

    recursive function factorial(n) result(fact)
        implicit none
        integer, intent(in) :: n
        integer :: fact

        if (n < 0) then
            print *, "Error: El factorial no está definido para números negativos."
            stop
        elseif (n == 0) then
            fact = 1  
        else
            fact = n * factorial(n - 1)  
        end if
    end function factorial

    function shift_1D(f) result(shift_f)
        !
        ! input
        !
        complex(kind=8), dimension(:), intent(in) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: shift_f
        !
        ! others
        !

        integer :: N; N = size(f)/2

        shift_f(N+1:2*N) = f(1:N)
        shift_f(1:N) = f(N+1:2*N)
    end function shift_1D



end module
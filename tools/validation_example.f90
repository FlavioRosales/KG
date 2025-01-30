program validation_example
    use fft

    implicit none
    integer, parameter :: n_r = 128, n_theta = 64, n_phi = 64, l_max = 20
    real(kind=8), parameter :: r_min = -4.0d0*pi, r_max = 4.0d0*pi
    real(kind=8), parameter :: a = 0.1d0, k0 = 2.0d0
    real(kind=8), dimension(n_r) :: r, frequencies
    real(kind=8), dimension(n_theta) :: theta
    real(kind=8), dimension(n_phi) :: phi
    real(kind=8) :: dx
    complex(kind=8), dimension(n_r, n_theta, n_phi) :: f
    complex(kind=8), dimension(n_r, 0:l_max, -l_max:l_max) :: ft
    complex(kind=8), dimension(n_r) :: analytic_fft
    integer :: i

    ! Generar las mallas radiales y angulares
    ! Malla radial
    dx =  (r_max - r_min) / dble(n_r)
    do i = 1, n_r
        r(i) = r_min + dble(i-1) * dx
    end do

    ! Mallas angulares
    do i = 1, n_theta
        theta(i) = (dble(i-1)) * pi / dble(n_theta)
    end do
    do i = 1, n_phi
        phi(i) = dble(i-1) * 2.0d0 * pi / dble(n_phi)
    end do

    do i=1, n_r 
        frequencies(i) =2.0d0 * pi * dble(i-1 - n_r/2) / (r_max - r_min)
    end do

    frequencies = dreal(shift_1D(dcmplx(frequencies)))

    !Definir la función f(r, theta, phi)
    call define_function(f, r, theta, phi, a, k0)

    !Transformada de Fourier en coordenadas esféricas
    ft = fft_spherical(dcmplx(f),r,theta,phi,l_max)

    !Calcular la transformada analítica
    do i = 1, n_r
        analytic_fft(i) = sqrt(pi / a) * exp(-((frequencies(i) - k0)**2) / (4.0d0 * a))
    end do

    ! Guardar resultados
    open(unit=10, file="validation_case.dat", status="replace")
    do i = 1, n_r
        write(10, *) frequencies(i), abs(ft(i, 2, 1)), abs(analytic_fft(i))
    end do
    close(10)
end program validation_example

subroutine define_function(f, r, theta, phi, a, k0)
    use fft

    implicit none
    integer, parameter :: n_r = 128, n_theta = 64, n_phi = 64
    real(kind=8), intent(in) :: r(n_r), theta(n_theta), phi(n_phi)
    real(kind=8), intent(in) :: a, k0
    complex(kind=8), intent(out) :: f(n_r, n_theta, n_phi)
    integer :: i, j, k
    complex(kind=8) :: ylm

    ! Calcular f(r, theta, phi)
    do i = 1, n_r
        do j = 1, n_theta
            do k = 1, n_phi
                ylm =  compute_spherical_harmonic(2, 1, theta(j), phi(k))
                f(i, j, k) = exp(-a * r(i)**2) * cos(2.0d0 * pi * k0 * r(i)) * ylm
            end do
        end do
    end do
end subroutine define_function


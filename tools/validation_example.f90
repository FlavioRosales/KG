program validation_example
    use fft

    implicit none
    integer, parameter :: n_r = 128, n_theta = 64, n_phi = 64, l_max = 20
    real(kind=8), parameter :: r_min = 0.0d0, r_max = 4.0d0*pi
    real(kind=8), parameter :: a = 0.1d0
    real(kind=8), dimension(n_r) :: frequencies
    real(kind=8), dimension(n_r, n_theta, n_phi) :: r
    real(kind=8), dimension(n_r, n_theta, n_phi) :: theta
    real(kind=8), dimension(n_r, n_theta, n_phi) :: phi
    real(kind=8) :: dx
    complex(kind=8), dimension(n_r, n_theta, n_phi) :: f
    complex(kind=8), dimension(n_r, 0:l_max, -l_max:l_max) :: ft
    complex(kind=8), dimension(n_r) :: analytic_fft
    integer :: i, l, m

    ! Generate radial and angular grids
    dx = (r_max - r_min) / dble(n_r)
    do i = 1, n_r
        r(i,:,:) = r_min + dble(i-1) * dx
    end do

    do i = 1, n_theta
        theta(:,i,:) = (dble(i-1)) * pi / dble(n_theta)
    end do
    do i = 1, n_phi
        phi(:,:,i) = dble(i-1) * 2.0d0 * pi / dble(n_phi)
    end do

    do i = 1, n_r 
        frequencies(i) = 2.0d0 * pi * dble(i-1 - n_r/2) / (r_max - r_min)
    end do

    ! Apply shift
    frequencies = dreal(shift_1D(dcmplx(frequencies)))

    call define_function(f, r(:,1,1), theta, phi, a)

    ! Compute the spherical Fourier transform    2.0d0*pi*dble(i-1 - Nx/2)/(dble(Nx)*dx)

    ft = fft_spherical(dcmplx(f), r, theta, phi, l_max)

    ! Compute the analytical Fourier transform
    do i = 1, n_r
        analytic_fft(i) = sqrt(pi / a) * exp(-(frequencies(i)**2) / (4.0d0 * a))
    end do

    ! Save results to file
    open(unit=10, file="validation_case.dat", status="replace")
    do i = 1, n_r
        write(10, *) frequencies(i), abs(ft(i, 2, 2)), abs(analytic_fft(i))
    end do
    close(10)

end program validation_example

subroutine define_function(f, r, theta, phi, a)
    use fft

    implicit none
    integer, parameter :: n_r = 128, n_theta = 64, n_phi = 64
    real(kind=8), intent(in) :: r(n_r), theta(n_r, n_theta, n_phi), phi(n_r, n_theta, n_phi)
    real(kind=8), intent(in) :: a
    complex(kind=8), intent(out) :: f(n_r, n_theta, n_phi)
    integer :: i
    complex(kind=8), dimension(n_theta, n_phi) :: ylm

    ! Compute f(r, theta, phi)
    ylm = compute_spherical_harmonic(2, 2, theta(1,:,:), phi(1,:,:))
    do i = 1, n_r
        f(i,:,:) = exp(-a * r(i)**2) * ylm
    end do

end subroutine define_function

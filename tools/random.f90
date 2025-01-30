program random
    use fft_lib
    implicit none 

    integer, parameter ::Nx=256,Ny=256, Nz=256
    real(kind=8) :: xmin, xmax, dx
    real(kind=8) :: ymin, ymax, dy
    real(kind=8) :: zmin, zmax, dz
    real(kind=8) :: p_0, random_value
    real(kind=8), dimension(Nx,Ny,Nz) :: x, y, z, omegax,omegay, omegaz
    complex(kind=8), dimension(Nx,Ny,Nz) :: psi_p, psi_x
    integer i, j, k

    
    xmin = -4.0d0*pi
    xmax = 4.0d0*pi
    dx = (xmax-xmin)/dble(Nx)

    ymin = -4.0d0*pi
    ymax = 4.0d0*pi
    dy = (ymax-ymin)/dble(Ny)

    zmin = -4.0d0*pi
    zmax = 4.0d0*pi
    dz = (zmax-zmin)/dble(Nz)

    do i = 1, Nx
    do j = 1, Ny
    do k = 1, Nz
        x(i,j,k) = xmin + dble(i-1)*dx
        y(i,j,k) = ymin + dble(j-1)*dy
        z(i,j,k) = zmin + dble(k-1)*dz

        omegax(i,j,k) = 2.0d0*pi*dble(i-1 - Nx/2)/(dble(Nx)*dx)
        omegay(i,j,k) = 2.0d0*pi*dble(j-1 - Ny/2)/(dble(Ny)*dy)
        omegaz(i,j,k) = 2.0d0*pi*dble(k-1 - Nz/2)/(dble(Nz)*dz)
    end do
    end do
    end do

    p_0 = 1.0d0

    psi_p = exp(-0.5d0/p_0 * (x**2 + y**2 + z**2))
    do i = 1, Nx 
        do j=1, Ny 
            do k=1, Nz
                call random_number(random_value)
                random_value = random_value * 2.0 * pi  ! Escala el valor al rango [0, 2pi]
                psi_p(i,j,k) = psi_p(i,j,k)*exp(i_unreal*random_value)
            end do 
        end do
     end do
     psi_x = ifft(psi_p)




     open(1,file='random_cart.dat')

     do i=1, Nx 
        do j=1, Ny 
            do k=1, Nz 
                write(1,*) omegax(i,j,k), omegay(i,j,k), omegaz(i,j,k), abs(psi_x(i,j,k))
            end do 
         !   write(1,*)
        end do
        !write(1,*)
     end do

     close(1)

     print*, 'x_max', maxval(omegax)
     print*, 'y_max', maxval(omegay)
     print*, 'z_max', maxval(omegaz)

     print*, 'x_min', minval(omegax)
     print*, 'y_min', minval(omegay)
     print*, 'z_min', minval(omegaz)
end program
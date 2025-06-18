module random
    use fft_lib
    use mpi_lib
    implicit none

    integer, parameter, private ::Nxx=256,Nyy=256, Nzz=256
    real(kind=8), private :: xmin, xmax, dx
    real(kind=8), private :: ymin, ymax, dy
    real(kind=8), private :: zmin, zmax, dz
    real(kind=8), private :: p_0
    real(kind=8), private, dimension(Nxx,Nyy,Nzz) :: x, y, z, omegax,omegay, omegaz
    complex(kind=8), dimension(Nxx,Nyy,Nzz), private :: psi_p, psi_x
    real(kind=8), dimension(Nxx,Nyy,Nzz), private :: random_value

    contains

    subroutine function_cartesian(rmax,m_boson,p_0) 
        implicit none 

        real(kind=8) :: rmax,m_boson, p_0

        integer :: i,j,k

        rmax = rmax + 1.0d0

        xmin = -rmax 
        xmax = rmax  
        dx = (xmax-xmin)/dble(Nxx)
    
        ymin = -rmax
        ymax = rmax
        dy = (ymax-ymin)/dble(Nyy)
    
        zmin = -rmax
        zmax = rmax
        dz = (zmax-zmin)/dble(Nzz)
    
        do i = 1, Nxx
        do j = 1, Nyy
        do k = 1, Nzz
            x(i,j,k) = xmin + dble(i-1)*dx
            y(i,j,k) = ymin + dble(j-1)*dy
            z(i,j,k) = zmin + dble(k-1)*dz
    
            omegax(i,j,k) = 2.0d0*pii*dble(i-1 - Nxx/2)/(dble(Nxx)*dx)
            omegay(i,j,k) = 2.0d0*pii*dble(j-1 - Nyy/2)/(dble(Nyy)*dy)
            omegaz(i,j,k) = 2.0d0*pii*dble(k-1 - Nzz/2)/(dble(Nzz)*dz)
        end do
        end do
        end do

        omegax = dreal(shift(dcmplx(omegax)))
        omegay = dreal(shift(dcmplx(omegay)))
        omegaz = dreal(shift(dcmplx(omegaz)))
    
        if (rank == 0) then
            call random_seed()
            call random_number(random_value)
        end if
    
        call MPI_BCAST(random_value, size(random_value), MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    
    
    
        psi_p = exp(-0.5d0/p_0**2 * (omegax**2 + omegay**2 + omegaz**2))*exp(i_img*2.0d0*pii*random_value)

        psi_x = ifft(psi_p)


    end subroutine

    subroutine spherical_to_cartesian(r, theta, phi, x_c, y_c, z_c)
        implicit none
        real(kind=8), intent(in)  :: r, theta, phi
        real(kind=8), intent(out) :: x_c, y_c, z_c
    
        x_c = r * dsin(theta) * dcos(phi)
        y_c = r * dsin(theta) * dsin(phi)
        z_c = r * dcos(theta)
    end subroutine spherical_to_cartesian


    function trilinear_interpolation(x_aux, y_aux, z_aux) result(psi_interp)
        implicit none
        real(kind=8), intent(in) :: x_aux, y_aux, z_aux
        complex(kind=8) :: psi_interp
        integer :: i, j, k
        real(kind=8) :: xd, yd, zd
        complex(kind=8) :: c000, c100, c010, c001, c110, c101, c011, c111

        do i = 1, Nxx-1
            if (x(i,1,1) <= x_aux .and. x(i+1,1,1) > x_aux) exit
        end do
        do j = 1, Nyy-1
            if (y(1,j,1) <= y_aux .and. y(1,j+1,1) > y_aux) exit
        end do
        do k = 1, Nzz-1
            if (z(1,1,k) <= z_aux .and. z(1,1,k+1) > z_aux) exit
        end do

        ! Calcular los factores de interpolación (fracciones dentro de la celda)
        xd = (x_aux - x(i,1,1)) / (x(i+1,1,1) - x(i,1,1))
        yd = (y_aux - y(1,j,1)) / (y(1,j+1,1) - y(1,j,1))
        zd = (z_aux - z(1,1,k)) / (z(1,1,k+1) - z(1,1,k))

        ! Verificar que xd, yd, zd estén en [0,1] para evitar errores
        xd = max(0.0d0, min(1.0d0, xd))
        yd = max(0.0d0, min(1.0d0, yd))
        zd = max(0.0d0, min(1.0d0, zd))

        ! Obtener valores en los vértices de la celda
        c000 = psi_x(i, j, k)
        c100 = psi_x(i+1, j, k)
        c010 = psi_x(i, j+1, k)
        c001 = psi_x(i, j, k+1)
        c110 = psi_x(i+1, j+1, k)
        c101 = psi_x(i+1, j, k+1)
        c011 = psi_x(i, j+1, k+1)
        c111 = psi_x(i+1, j+1, k+1)
        
        ! Aplicar la interpolación trilineal
        psi_interp = c000*(1-xd)*(1-yd)*(1-zd) + c100*xd*(1-yd)*(1-zd) + &
                    c010*(1-xd)*yd*(1-zd) + c001*(1-xd)*(1-yd)*zd + &
                    c110*xd*yd*(1-zd) + c101*xd*(1-yd)*zd + &
                    c011*(1-xd)*yd*zd + c111*xd*yd*zd


    end function trilinear_interpolation


    subroutine convert_to_spherical(this,p_0)
        use lines_method_lib
        implicit none
        class(lines_method), intent(in out) :: this
        integer :: i, j, k
        real(kind=8) :: x_c, y_c, z_c, r, theta, phi, p_0

        call function_cartesian(this%r_max,this%m,p_0)

        do i = this%iL - this%g_pts, this%iR + this%g_pts
            do j = this%jL - this%g_pts, this%jR + this%g_pts
                do k = this%kL, this%kR
                    r = this%x(i,j,k)
                    theta = this%y(i,j,k)
                    phi = this%z(i,j,k)

                    ! Convertir a coordenadas cartesianas
                    call spherical_to_cartesian(r, theta, phi, x_c, y_c, z_c)

                    ! Interpolar el campo psi_x en la nueva malla esférica
                    this%u(1,i, j, k) = trilinear_interpolation(x_c, y_c, z_c)

                end do
            end do
        end do
    end subroutine convert_to_spherical





end module
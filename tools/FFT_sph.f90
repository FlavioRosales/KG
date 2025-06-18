module FFT_sph
    use iso_c_binding
    use mpi_lib
    use integral
    implicit none 

    real(kind=8), parameter  :: pi = acos(-1.0d0)
    complex(kind=8), parameter :: i_unreal = (0.0d0,1.0d0)
    

    contains

    subroutine generate_full_bessel_basis(j_basis, k_values, r, R_max, Nk, l_max)
        use integral
        implicit none
        integer, intent(in) :: l_max, Nk
        real(kind=8), intent(in) :: r(:), R_max
        real(kind=8), intent(out) :: k_values(0:l_max, Nk)
        real(kind=8), intent(out) :: j_basis(0:l_max, Nk, size(r))
    
        interface
            function gsl_sf_bessel_zero_Jnu(nu, n) bind(C, name="gsl_sf_bessel_zero_Jnu")
                use iso_c_binding
                real(c_double), value :: nu
                integer(c_int), value :: n
                real(c_double) :: gsl_sf_bessel_zero_Jnu
            end function gsl_sf_bessel_zero_Jnu
    
            function gsl_sf_bessel_jl(l, x) bind(C, name="gsl_sf_bessel_jl")
                use iso_c_binding
                integer(c_int), value :: l
                real(c_double), value :: x
                real(c_double) :: gsl_sf_bessel_jl
            end function gsl_sf_bessel_jl
        end interface
    
        integer :: l, n, i, m
        real(kind=8) :: z_n, kr, k_val, dr, int_local, int_global
        real(kind=8), allocatable :: integrando(:)
        character(len=10) :: lstr
    
        dr = abs(r(2) - r(1))
        allocate(integrando(size(r)))
    
        do l = 0, l_max
            do n = 1, Nk                
                z_n = gsl_sf_bessel_zero_Jnu(dble(l) + 0.5d0, n)
                k_val = z_n / R_max
                k_values(l,n) = k_val
    
                do i = 1, size(r)
                    kr = k_val * r(i)
                    j_basis(l,n,i) = gsl_sf_bessel_jl(l, kr)
                end do
    
                do i = 1, size(r)
                    integrando(i) = j_basis(l,n,i)**2 * r(i)**2
                end do
    
                int_local = trapezium_1D(integrando, dr)
                call MPI_ALLREDUCE(int_local, int_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                int_global = sqrt(int_global)

    
                if (int_global > 0.d0) then
                    do i = 1, size(r)
                        j_basis(l,n,i) = j_basis(l,n,i) / int_global
                    end do
                else
                    print *, "[Warning] Normalization failed: int_global=0 at l=", l, " n=", n
                end if
            end do

            write(lstr, '(I0)') l
            open(10, file='bessel_basis_l'//trim(adjustl(lstr))//'.dat', status='replace')

            do n=1, Nk
                do m= 1, Nk 
                    do i = 1, size(r)
                        integrando(i) = j_basis(l,n,i)*j_basis(l,m,i) * r(i)**2
                    end do
                    int_local = trapezium_1D(integrando, dr)
                    call MPI_ALLREDUCE(int_local, int_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                    write(10, *) n, m, int_global
                end do 
                write(10, *) 
            end do 

            close(unit=10)

        end do
    
        deallocate(integrando)
        
    end subroutine generate_full_bessel_basis
    
    function ISphFT_discrete_radial(ft, r, theta, phi, Ylm, j_basis, l_max, Nk) result(ift)
        use integral
        implicit none
        integer, intent(in) :: l_max, Nk
        complex(kind=8), intent(in) :: ft(Nk, 0:l_max, -l_max:l_max)
        real(kind=8), intent(in) :: r(:), theta(:), phi(:)
        complex(kind=8), intent(in) :: Ylm(0:l_max, -l_max:l_max, size(theta), size(phi))
        real(kind=8), intent(in) :: j_basis(0:l_max, Nk, size(r))
        complex(kind=8) :: ift(size(r), size(theta), size(phi))
        complex(kind=8) :: flm_r(size(r), 0:l_max, -l_max:l_max)
        integer :: l, m, i, j, k, n
    
        ift = (0.0d0, 0.0d0)
        flm_r = (0.0d0, 0.0d0)
    
        do l = 0, l_max
                do i = lbound(r,1), ubound(r,1)
                    flm_r(i,l,:) = (0.0d0, 0.0d0)
                    do n = 1, Nk
                        flm_r(i,l,:) = flm_r(i,l,:) + ft(n,l,:) * dcmplx(j_basis(l,n,i))
                    end do
                end do
        end do
    
        do i = lbound(r,1), ubound(r,1)
            do j = lbound(theta,1), ubound(theta,1)
                do k = lbound(phi,1), ubound(phi,1)
                    ift(i,j,k) = (0.0d0, 0.0d0)
                    do l = 0, l_max
                        do m = -l, l
                            ift(i,j,k) = ift(i,j,k) + flm_r(i,l,m) * Ylm(l,m,j,k)
                        end do
                    end do
                end do
            end do
        end do
    end function ISphFT_discrete_radial


    
      


    ! function compute_mps_general(f, r, theta, phi, Ylm, l_max, Nk, index) result(Pk)
    !     use integral
    !     implicit none
      
    !     ! Entradas
    !     integer, intent(in) :: l_max, Nk
    !     integer, intent(in) :: index(3,2)
    !     real(kind=8), intent(in) :: f(index(1,1):index(1,2),index(2,1):index(2,2),index(3,1):index(3,2))
    !     real(kind=8), intent(in) :: r(index(1,1):index(1,2)), theta(index(2,1):index(2,2)), phi(index(3,1):index(3,2))
    !     complex(kind=8), intent(in) :: Ylm(0:l_max, -l_max:l_max, index(2,1):index(2,2), index(3,1):index(3,2))
      
    !     ! Salida
    !     real(kind=8) :: Pk(0:Nk)
      
    !     ! Variables locales
    !     integer :: l, m, n, i, j, k
    !     integer :: g_pts
    !     real(kind=8) :: dr, dtheta, dphi
    !     real(kind=8) :: r_i, sin_theta_j
    !     real(kind=8) :: J_l(0:l_max,0:Nk,index(1,1):index(1,2))
    !     complex(kind=8) :: a_nlm
    !     real(kind=8) :: k_values(0:Nk)
    !     real(kind=8) :: f_vol(index(1,1):index(1,2),index(2,1):index(2,2),index(3,1):index(3,2))
    !     real(kind=8) :: delta(index(1,1):index(1,2),index(2,1):index(2,2),index(3,1):index(3,2))
    !     real(kind=8) :: integrand(index(1,1):index(1,2),index(2,1):index(2,2),index(3,1):index(3,2))
    !     real(kind=8) :: vol_elem(index(1,1):index(1,2),index(2,1):index(2,2),index(3,1):index(3,2))

      
      
    !     dr = abs(r(index(1,1)+1) - r(index(1,1)))
    !     dtheta = abs(theta(index(2,1)+1) - theta(index(2,1)))
    !     dphi = abs(phi(index(3,1)+1) - phi(index(3,1)))
    !     g_pts = abs(index(1,1))

    !     do i=0, Nk 
    !         k_values(i) = 0.01d0 + dble(i)*(pi/dr - 0.01d0)/Nk
    !     end do 

    !     J_l = generate_bessel_basis(r, k_values, l_max, Nk, index)
    !     Pk = 0.0d0

    !     do i=index(1,1), index(1,2)
    !         do j=index(2,1), index(2,2)
    !             do k=index(3,1), index(3,2)
    !                 f_vol(i,j,k) = f(i,j,k)*r(i)**2*dsin(theta(j))
    !                 vol_elem(i,j,k) = r(i)**2 * dsin(theta(j))
    !             end do 
    !         end do 
    !     end do

    !     delta = (f - integrate(f_vol,dr,dtheta,dphi,g_pts))/integrate(vol_elem,dr,dtheta,dphi,g_pts)

    !     do n = 0, Nk
    !         do l = 0, l_max
    !         do m = -l, l
    !             ! Preparar integrando complejo
    !             do i=index(1,1), index(1,2)
    !                 do j=index(2,1), index(2,2)
    !                     do k=index(3,1), index(3,2)
    !                         integrand(i,j,k) = delta(i,j,k)*J_l(l,n,i)*dconjg(Ylm(l,m,j,k))*vol_elem(i,j,k)
    !                     end do 
    !                 end do 
    !             end do
    !             a_nlm = integrate(integrand, dr, dtheta, dphi, g_pts)
    !             Pk(n) = Pk(n) + abs(a_nlm)**2
    !         end do
    !         end do
    !     end do
  

    ! end function compute_mps_general

      
    function SphFT_discrete_radial(f, r, theta, Ylm, l_max,Nk, j_basis, dr, dtheta, dphi) result(ft)
        use integral
        implicit none
        integer, intent(in) :: l_max, Nk
        complex(kind=8), intent(in) :: f(:,:,:)
        real(kind=8), intent(in) :: r(:), theta(:), dr, dtheta, dphi
        complex(kind=8), intent(in) :: Ylm(0:l_max, -l_max:l_max, &
                                           lbound(theta,1):ubound(theta,1), lbound(f,3):ubound(f,3))
        real(kind=8), intent(in) :: j_basis(0:l_max, 0:Nk, lbound(r,1):ubound(r,1))  ! [l, n, r]
        complex(kind=8) :: ft(0:Nk, 0:l_max, -l_max:l_max)
        complex(kind=8) :: int
        complex(kind=8) :: flm(lbound(r,1):ubound(r,1), 0:l_max, -l_max:l_max)
        integer :: i, l, m, n

        ft = (0.0d0, 0.0d0)

        do i = lbound(r,1), ubound(r,1)
            flm(i,:,:) = project_sh(f(i,:,:), Ylm, l_max, theta, dtheta, dphi)
        end do

        do l = 0, l_max
            do m = -l, l
                do n = 0, Nk
                    int = trapezium_1D(flm(:,l,m) * dcmplx(j_basis(l,n,:) * r(:)**2), dr)

                    call MPI_ALLREDUCE(int, ft(n,l,m), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
                
                end do
            end do
        end do
        
    end function SphFT_discrete_radial

    
    function project_sh(f,Ylm,l_max,theta,dtheta,dphi) result(flm)
        use integral
        implicit none 
        !input
        integer, intent(in) :: l_max !maximum of l in spherical harmonics
        complex(kind=8), intent(in) :: f(:,:) !function slice (theta,phi)
        real(kind=8), intent(in) :: theta(:)
        complex(kind=8), intent(in) :: Ylm(0:l_max,-l_max:l_max,&
                                           size(theta), &
                                           size(f,2))
        real(kind=8), intent(in) :: dtheta, dphi
        !output
        complex(kind=8) :: flm(0:l_max,-l_max:l_max)
        !internal
        complex(kind=8) :: fun(size(theta),size(f,2))
        integer :: l,m,j,k
        integer :: Nth, Nph

        Nth = size(theta)
        Nph = size(f,2)

        do l=0, l_max
            do m=-l, l 
                do j=1, Nth
                    do k=1,Nph
                        fun(j,k) = f(j,k)*dconjg(Ylm(l,m,j,k))*dcmplx(dsin(theta(j)))
                    end do
                end do 
                flm(l,m) = trapezium_2D(fun(:,:),dtheta,dphi)
            end do 
        end do

    end function project_sh


    subroutine sph_harmonics(Ylm, theta, phi, l_max)
        use iso_c_binding
        implicit none
        real(kind=8), intent(in) :: theta(:), phi(:)
        integer, intent(in) :: l_max
        complex(kind=8), intent(out) :: Ylm(0:l_max,-l_max:l_max,size(theta),size(phi))


        interface
            function gsl_sf_legendre_Plm(l, m, x) bind(C, name="gsl_sf_legendre_Plm")
                use iso_c_binding
                integer(c_int), value :: l, m
                real(c_double), value :: x
                real(c_double) :: gsl_sf_legendre_Plm
            end function
        end interface
        
        integer :: l, m, i, j
        real(kind=8) :: x, Plm, norm, fact
        complex(kind=8) :: phase
        
        Ylm = (0.0d0, 0.0d0)
        
        
        do l = 0, l_max
            fact = sqrt((2.0d0*l + 1.0d0)/(4.0d0*pi))
            do m = 0, l
                ! Normalization factor
                norm = fact
                if (m > 0) then
                    do i = 1, 2*m
                        norm = norm / sqrt(dble(l-m+i))
                    end do
                end if
                
                do i = 1, size(theta)
                    x = cos(theta(i))
                    ! Handle poles explicitly
                    if (abs(x) > 0.9999d0) then
                        Plm = merge(fact, 0.0d0, m==0)
                    else
                        Plm = gsl_sf_legendre_Plm(l, m, x)
                    end if
                    
                    do j = 1, size(phi)
                        phase = cmplx(cos(m*phi(j)), sin(m*phi(j)), kind=8)
                        Ylm(l,m,i,j) = norm * Plm * phase
                        if (m > 0) then
                            Ylm(l,-m,i,j) = (-1.0d0)**m * conjg(Ylm(l,m,i,j))
                        end if
                    end do
                end do
            end do
        end do
    end subroutine sph_harmonics

    subroutine test_orthonormality(Ylm, theta, phi, l_max)
        use iso_c_binding
        use integral
        implicit none
        integer, intent(in) :: l_max
        real(kind=8), intent(in) :: theta(:), phi(:)
        complex(kind=8), intent(in) :: Ylm(0:l_max, -l_max:l_max, size(theta), size(phi))
        
        real(kind=8) :: dtheta, dphi, integral_val, tolerance
        integer :: l1, m1, l2, m2, i, j, n_errors
        complex(kind=8), allocatable :: integrand(:,:)
        
        ! Configuración
        tolerance = 1.0d-5
        n_errors = 0
        
        ! Calculamos los diferenciales (asumiendo malla uniforme)
        dtheta = theta(2) - theta(1)
        dphi = phi(2) - phi(1)
        
        allocate(integrand(size(theta), size(phi)))
        
        print *, "============================================"
        print *, "  Spherical Harmonics Orthonormality Test   "
        print *, "  Maximum l: ", l_max
        print *, "  Tolerance: ", tolerance
        print *, "  Grid size: ", size(theta), "θ ×", size(phi), "φ"
        print *, "============================================"
        
        do l1 = 0, l_max
            do m1 = -l1, l1
                do l2 = 0, min(l1+1, l_max)  ! Solo probamos l2 cercanos a l1 para eficiencia
                    do m2 = -l2, l2
                        ! Calcular el integrando: Y*_l1m1 · Y_l2m2 · sinθ
                        do i = 1, size(theta)
                            do j = 1, size(phi)
                                integrand(i,j) = dconjg(Ylm(l1,m1,i,j))*Ylm(l2,m2,i,j)*dsin(theta(i))
                            end do
                        end do
                        
                        ! Integración 2D con regla del trapecio
                        integral_val = abs(trapezium_2D(integrand,dtheta,dphi))
                        
                        ! Evaluación
                        if (l1 == l2 .and. m1 == m2) then
                            ! Caso diagonal: debería ser ≈1
                            if (abs(integral_val - 1.0d0) > tolerance) then
                                print *, "FAIL: (l,m)=(", l1, ",", m1, ")", integral_val
                                n_errors = n_errors + 1
                            end if
                        else
                            ! Caso no diagonal: debería ser ≈0
                            if (abs(integral_val) > tolerance) then
                                print *, "FAIL: (", l1, ",", m1, ")|(", l2, ",", m2, ")", integral_val
                                n_errors = n_errors + 1
                            end if
                        end if
                    end do
                end do
            end do
        end do
        
        deallocate(integrand)
        
        print *, "============================================"
        print *, "  Test completed"
        print *, "  Total errors: ", n_errors
        print *, "============================================"
    end subroutine


    subroutine sph_to_domain(f, r, theta, phi, f_d, r_d, theta_d, phi_d, index_bounds)
        implicit none
      
        real(kind=8), intent(in) :: r(:,:,:), theta(:,:,:), phi(:,:,:)
        complex(kind=8), intent(in) :: f(size(r,1), size(r,2), size(r,3))
        integer, intent(in) :: index_bounds(3,2)
        real(kind=8), dimension(index_bounds(1,1):index_bounds(1,2), &
                                index_bounds(2,1):index_bounds(2,2), &
                                index_bounds(3,1):index_bounds(3,2)), intent(in) :: r_d, theta_d, phi_d
      
        complex(kind=8), dimension(index_bounds(1,1):index_bounds(1,2), &
                                   index_bounds(2,1):index_bounds(2,2), &
                                   index_bounds(3,1):index_bounds(3,2)) :: f_d
      
        integer :: i, j, k
        integer :: i_r, i_th, i_ph, i_ph_p1
        real(kind=8) :: rr, tt, pp
        real(kind=8) :: w_r, w_th, w_ph
        complex(kind=8) :: f_interp
      
        integer :: Nr, Nth, Nph
        real(kind=8), allocatable :: r1d(:), th1d(:), ph1d(:)
        real(kind=8) :: rmin, dr, thetamin, dth, phimin, dph, twopi
      
        Nr = size(r,1)
        Nth = size(r,2)
        Nph = size(r,3)
      
        allocate(r1d(Nr), th1d(Nth), ph1d(Nph))
        r1d  = r(:,1,1)
        th1d = theta(1,:,1)
        ph1d = phi(1,1,:)
      
        rmin     = r1d(1)
        dr       = r1d(2) - r1d(1)
        thetamin = th1d(1)
        dth      = th1d(2) - th1d(1)
        phimin   = ph1d(1)
        dph      = ph1d(2) - ph1d(1)
      
        twopi = 2.0d0 * 4.0d0 * atan(1.0d0)  ! 2π
      
        do k = index_bounds(3,1), index_bounds(3,2)
          do j = index_bounds(2,1), index_bounds(2,2)
            do i = index_bounds(1,1), index_bounds(1,2)
      
              rr = r_d(i,j,k)
              tt = theta_d(i,j,k)
              pp = modulo(phi_d(i,j,k), twopi)
      
              ! Estimación inicial
              i_r  = int((rr - rmin)     / dr)  + 1
              i_th = int((tt - thetamin) / dth) + 1
              i_ph = int((pp - phimin)   / dph) + 1
      
              ! Clamping
              if (i_r  < 1)     i_r  = 1
              if (i_r  > Nr-1)  i_r  = Nr-1
              if (i_th < 1)     i_th = 1
              if (i_th > Nth-1) i_th = Nth-1
              if (i_ph < 1)     i_ph = 1
              if (i_ph > Nph)   i_ph = Nph
      
              do while (i_r < Nr-1 .and. rr >= r1d(i_r+1))
                i_r = i_r + 1
              end do
              do while (i_r > 1 .and. rr < r1d(i_r))
                i_r = i_r - 1
              end do
      
              do while (i_th < Nth-1 .and. tt >= th1d(i_th+1))
                i_th = i_th + 1
              end do
              do while (i_th > 1 .and. tt < th1d(i_th))
                i_th = i_th - 1
              end do
      
              do while (i_ph < Nph .and. pp >= ph1d(i_ph+1))
                i_ph = i_ph + 1
              end do
              do while (i_ph > 1 .and. pp < ph1d(i_ph))
                i_ph = i_ph - 1
              end do
      
              ! Pesos
              w_r  = (rr - r1d(i_r))     / dr
              w_th = (tt - th1d(i_th))   / dth
              w_ph = (pp - ph1d(i_ph))   / dph
      
              ! Índice periódico φ+1
              if (i_ph == Nph) then
                i_ph_p1 = 1
              else
                i_ph_p1 = i_ph + 1
              end if
      
              ! Interpolación trilineal con φ periódico
              f_interp = (1.0d0 - w_r)*(1.0d0 - w_th)*(1.0d0 - w_ph)*f(i_r    , i_th    , i_ph    ) + &
                         (1.0d0 - w_r)*(1.0d0 - w_th)*(      w_ph)*f(i_r    , i_th    , i_ph_p1 ) + &
                         (1.0d0 - w_r)*(      w_th)*(1.0d0 - w_ph)*f(i_r    , i_th+1  , i_ph    ) + &
                         (1.0d0 - w_r)*(      w_th)*(      w_ph)*f(i_r    , i_th+1  , i_ph_p1 ) + &
                         (      w_r)*(1.0d0 - w_th)*(1.0d0 - w_ph)*f(i_r+1  , i_th    , i_ph    ) + &
                         (      w_r)*(1.0d0 - w_th)*(      w_ph)*f(i_r+1  , i_th    , i_ph_p1 ) + &
                         (      w_r)*(      w_th)*(1.0d0 - w_ph)*f(i_r+1  , i_th+1  , i_ph    ) + &
                         (      w_r)*(      w_th)*(      w_ph)*f(i_r+1  , i_th+1  , i_ph_p1 )
      
              f_d(i,j,k) = f_interp
      
            end do
          end do
        end do
      
        deallocate(r1d, th1d, ph1d)
      
      end subroutine sph_to_domain

      subroutine domain_to_sph(f_d, r_d, theta_d, phi_d, f, r, theta, phi, index_bounds)
        implicit none
      
        ! === Entrada ===
        integer, intent(in) :: index_bounds(3,2)
        real(kind=8), intent(in), dimension(index_bounds(1,1):index_bounds(1,2), &
                                            index_bounds(2,1):index_bounds(2,2), &
                                            index_bounds(3,1):index_bounds(3,2)) :: r_d, theta_d, phi_d
        complex(kind=8), intent(in), dimension(index_bounds(1,1):index_bounds(1,2), &
                                               index_bounds(2,1):index_bounds(2,2), &
                                               index_bounds(3,1):index_bounds(3,2)) :: f_d
        real(kind=8), intent(in) :: r(:,:,:), theta(:,:,:), phi(:,:,:)
      
        ! === Salida ===
        complex(kind=8) :: f(:, :, :)
      
        ! === Variables internas ===
        integer :: i, j, k
        integer :: i_r, i_th, i_ph, i_ph_p1
        real(kind=8) :: rr, tt, pp
        real(kind=8) :: w_r, w_th, w_ph
        complex(kind=8) :: f_interp
        complex(kind=8) :: f000, f001, f010, f011, f100, f101, f110, f111
        complex(kind=8) :: f00, f01, f10, f11, f0, f1
      
        real(kind=8) :: rmin, dr, thetamin, dth, phimin, dph, twopi
      
        ! === Cálculo de pasos a partir del dominio físico ===
        rmin     = r_d(index_bounds(1,1),1,1)
        dr       = r_d(index_bounds(1,1)+1,1,1) - r_d(index_bounds(1,1),1,1)
      
        thetamin = theta_d(1,index_bounds(2,1),1)
        dth      = abs(theta_d(1,index_bounds(2,1)+1,1) - theta_d(1,index_bounds(2,1),1))
      
        phimin   = phi_d(1,1,index_bounds(3,1))
        dph      = abs(phi_d(1,1,index_bounds(3,1)+1) - phi_d(1,1,index_bounds(3,1)))
      
        ! === Bucle sobre el dominio espectral (con polos y sin zonas fantasma) ===
        do k = 1, size(r,3)
          do j = 1, size(r,2)
            do i = 1, size(r,1)
      
              rr = r(i,j,k)
              tt = theta(i,j,k)
              pp = phi(i,j,k)
      
              ! --- Estimación inicial de índices ---
              i_r  = int((rr - rmin)     / dr)  + index_bounds(1,1)
              i_th = int((tt - thetamin) / dth) + index_bounds(2,1)
              i_ph = int((pp - phimin)   / dph) + index_bounds(3,1)
      
              ! --- Clamping de índices ---
              if (i_r  < index_bounds(1,1))     i_r = index_bounds(1,1)
              if (i_r  > index_bounds(1,2)-1)   i_r = index_bounds(1,2)-1
              if (i_th < index_bounds(2,1))     i_th = index_bounds(2,1)
              if (i_th > index_bounds(2,2)-1)   i_th = index_bounds(2,2)-1
              if (i_ph < index_bounds(3,1))     i_ph = index_bounds(3,1)
              if (i_ph > index_bounds(3,2)-1)   i_ph = index_bounds(3,2)-1
      
              ! --- Corrección robusta ---
              do while (i_r < index_bounds(1,2)-1 .and. rr >= r_d(i_r+1,1,1))
                i_r = i_r + 1
              end do
              do while (i_r > index_bounds(1,1) .and. rr < r_d(i_r,1,1))
                i_r = i_r - 1
              end do
      
              do while (i_th < index_bounds(2,2)-1 .and. tt >= theta_d(1,i_th+1,1))
                i_th = i_th + 1
              end do
              do while (i_th > index_bounds(2,1) .and. tt < theta_d(1,i_th,1))
                i_th = i_th - 1
              end do
      
              do while (i_ph < index_bounds(3,2)-1 .and. pp >= phi_d(1,1,i_ph+1))
                i_ph = i_ph + 1
              end do
              do while (i_ph > index_bounds(3,1) .and. pp < phi_d(1,1,i_ph))
                i_ph = i_ph - 1
              end do
      
              ! --- Manejo de polos ---
              if (tt <= thetamin) then
                f(i,j,k) = f_d(i_r, index_bounds(2,1), i_ph)
                cycle
              else if (tt >= thetamin + dth * (index_bounds(2,2) - index_bounds(2,1))) then
                f(i,j,k) = f_d(i_r, index_bounds(2,2), i_ph)
                cycle
              end if
      
              ! --- Índice siguiente en phi (sin periodicidad) ---
              i_ph_p1 = i_ph + 1
              if (i_ph_p1 > index_bounds(3,2)) cycle  ! evita overflow
      
              ! === Interpolación trilineal explícita ===
              f000 = f_d(i_r    , i_th    , i_ph    )
              f001 = f_d(i_r    , i_th    , i_ph_p1 )
              f010 = f_d(i_r    , i_th+1  , i_ph    )
              f011 = f_d(i_r    , i_th+1  , i_ph_p1 )
              f100 = f_d(i_r+1  , i_th    , i_ph    )
              f101 = f_d(i_r+1  , i_th    , i_ph_p1 )
              f110 = f_d(i_r+1  , i_th+1  , i_ph    )
              f111 = f_d(i_r+1  , i_th+1  , i_ph_p1 )
      
              ! Pesos
              w_r  = (rr - r_d(i_r,1,1))       / dr
              w_th = (tt - theta_d(1,i_th,1))  / dth
              w_ph = (pp - phi_d(1,1,i_ph))    / dph
      
              ! Interpolación
              f00 = (1.0d0 - w_ph)*f000 + w_ph*f001
              f01 = (1.0d0 - w_ph)*f010 + w_ph*f011
              f10 = (1.0d0 - w_ph)*f100 + w_ph*f101
              f11 = (1.0d0 - w_ph)*f110 + w_ph*f111
      
              f0 = (1.0d0 - w_th)*f00 + w_th*f01
              f1 = (1.0d0 - w_th)*f10 + w_th*f11
      
              f_interp = (1.0d0 - w_r)*f0 + w_r*f1
      
              f(i,j,k) = f_interp
      
            end do
          end do
        end do
      
      end subroutine domain_to_sph
      
    
end module
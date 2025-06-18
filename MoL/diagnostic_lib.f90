module diagnostic_lib
    use lines_method_lib
    use FFT_sph


    integer, private :: count, status(MPI_STATUS_SIZE), request

    type, extends(lines_method) :: scalar_field
        real(kind=8) :: Mass_scalar
        complex(kind=8), allocatable, dimension(:,:,:,:) :: j_i
        real(kind=8), allocatable, dimension(:) :: r_det
        real(kind=8), allocatable, dimension(:,:) :: momentum
        real(kind=8), allocatable, dimension(:,:,:) :: rho
        complex(kind=8), allocatable, dimension(:,:,:) :: alm
        complex(kind=8), allocatable, dimension(:,:,:,:) :: Ylm
        real(kind=8), allocatable, dimension(:,:,:) :: j_basis
        real(kind=8), allocatable, dimension(:,:) :: k_values
        complex(kind=8), allocatable, dimension(:,:,:,:) :: u_inter
        integer, allocatable, dimension(:) :: index_det
        integer :: save

    contains

        procedure, private :: scalar_field_memory
        !procedure :: calculate_mass_power_spectrum
        procedure :: currents_ADM
        procedure :: currents_detectors
        procedure :: proyector_detectors
        procedure ::  density
        procedure :: diagnostic

    end type scalar_field

    interface sf_constructor
        procedure :: scalar_field_constructor 
    end interface

    contains 

    function scalar_field_constructor( &
        nvars, g_pts, xmin, xmax, Nx, Ny, Nz, tmin, tmax, CFL, &
        boundary_type_rmax, metric, a, M_bh, m,r_max, save, n_det) result(this)
       implicit none
       integer, intent(in) ::  nvars, g_pts, Nx, Ny, Nz
       integer, intent(in) ::  n_det
       real(kind=8), intent(in) :: a, M_bh
       real(kind=8), intent(in) :: save
       real(kind=8), intent(in) :: xmin, xmax, tmin, tmax, CFL, m
       real(kind=8) :: dx_aux, dxx, r_max, r_min 
       integer :: dis, g, k, k_target 
       character(len=*) :: boundary_type_rmax, metric
       type(scalar_field) :: this
       integer i, n
   
       call MPI_CREATE
       !
       ! constructor for mesh refinement
       !
       if(mod(Nx, nproc).ne.0) stop "Error Nx must be divisible by nproc."
   
       this%nvars =  nvars
       this%g_pts = g_pts
       this%Nx = Nx / nproc; this%Ny = Ny ; this%Nz = Nz
       this%iL = 0; this%iR = this%Nx
       this%jL = 0; this%jR = this%Ny
       this%kL = 0; this%kR = this%Nz
       this%m = m
       this%a = a
       this%M_bh = M_bh 
   
   
       this%n_det = n_det
   
   
       this%xmin = xmin
       this%xmax = xmax
   
       this%boundary_type_rmax =  boundary_type_rmax
   
       print*, 'Boundary type r_max = ', this%boundary_type_rmax 
   
       call this%mesh_refinement_memory()
       call this%mesh_refinement_boundary(xmin, xmax)
       call this%mesh_refinement_create()
       call this%lines_method_memory()
       call this%scalar_field_memory()

       this%tmin = tmin; this%tmax = tmax; this%CFL = CFL
   
         dxx = this%x(this%iR+this%g_pts,1,1)
         call MPI_Bcast(dxx, 1, MPI_DOUBLE, nproc-1, MPI_COMM_WORLD, ierr)
         this%r_max = dxx
   
       this%index_bounds(1,1) = this%iL - this%g_pts !index x left
       this%index_bounds(1,2) = this%iR + this%g_pts !index x right
       this%index_bounds(2,1) = this%jL - this%g_pts !index y left
       this%index_bounds(2,2) = this%jR + this%g_pts !index y right
       this%index_bounds(3,1) = this%kL !index z left
       this%index_bounds(3,2) = this%kR !index z left
   
       dxx = this%x(this%iL-this%g_pts+1,1,1) - this%x(this%iL-this%g_pts,1,1)
       call MPI_Bcast(dxx, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
       this%dt = CFL * min(dxx, exp(xmin)*this%dy, &
         exp(xmin)*dsin(this%dy/2.0d0)*this%dz) 
   
       this%Nt = int((this%tmax - this%tmin) / this%dt)
       if(mod(this%Nt,2).ne.0) this%Nt = this%Nt + 1
       this%t = this%tmin
       call this%mesh_refinement_info()
       call this%lines_method_info()
   
       this%save = this%tmin + ceiling(save/this%dt)

       if (rank.eq.master) print*, 'Saving every: ', this%save, 'iterations', this%save*this%dt, 'time' 
   
   
       ! goemetry
       this%metric = metric
       call this%geometry_memory
       call this%geometry_load_metric
       call this%geometry_load_inverse
       call this%validate_inverse_metric
   
       this%r_det(1) = (1.0d0 + dsqrt(1.0d0 - this%a**2)) !detector horizon
       
       if(metric == 'MKS') then ! detector inifinity
         this%r_det(n_det) = exp(xmax)
       else  
         this%r_det(n_det) = xmax*0.99d0
       end if
   
       do n=2, n_det - 1
        this%r_det(n) = this%r_det(1) + (this%r_det(n_det)-this%r_det(1))/dble(n_det) * dble(n - 1)
       end do 
   
       


       r_min = this%x(this%iL-this%g_pts,1,1)
       call MPI_Bcast(r_min, 1, MPI_DOUBLE, master, MPI_COMM_WORLD, ierr)
   
   
   
       call sph_harmonics(this%ylm,this%theta(1,:,1),this%phi(1,1,:),l_max)

       !call test_orthonormality(this%ylm, this%theta(1,:,1),this%phi(1,1,:), l_max)
   
       call generate_full_bessel_basis(this%j_basis, this%k_values, &
                           this%r(:,1,1), this%r_max, Nk, l_max) 

    end function scalar_field_constructor

    subroutine scalar_field_memory(this)
        implicit none 
        class(scalar_field), intent(in out) :: this

        allocate(this%j_i(3, &
        this%iL - this%g_pts:this%iR + this%g_pts, &
        this%jL - this%g_pts : this%jR + this%g_pts , &
        this%kL:this%kR))

        allocate(this%rho(&
        this%iL - this%g_pts:this%iR + this%g_pts, &
        this%jL - this%g_pts : this%jR + this%g_pts , &
        this%kL:this%kR))

        allocate(this%Ylm(0:l_max, -l_max:l_max, &
        this%Ny, this%Nz))
    
        allocate(this%alm(this%n_det,0:l_max, -l_max:l_max))
    
        allocate(this%momentum(3,this%n_det))
    
        allocate(this%r_det(this%n_det))
        allocate(this%index_det(this%n_det))
    
        allocate(this%j_basis(0:l_max, Nk, this%Nx))
    
        allocate(this%k_values(0:l_max, Nk))

        allocate(this%u_inter(5,this%Nx,this%Ny,this%Nz))

    end subroutine scalar_field_memory

    subroutine diagnostic(this,index)
        use integral
        use save_data_lib
        implicit none
        class(scalar_field), intent(in out) :: this
        integer, intent(in) :: index
        integer :: l, i
        integer, parameter :: Nbins = int(0.1d0*(Nk + 1) * (l_max + 1)**2)
        real(kind=8), dimension(0:Nk) :: Pk, k_bins


        if(mod(index,this%save) == 0) then 
            call this%currents_ADM()
            call this%currents_detectors()
            call this%proyector_detectors()
            !call this%calculate_mass_power_spectrum(MPS,k_bins)
           ! Pk = compute_mps_general(this%m*abs(this%u(1,:,:,:))**2, this%x(:,1,1), this%y(1,:,1), this%z(1,1,:), & 
            !                         this%Ylm, l_max, Nk, this%index_bounds)
           ! do i=0, Nk 
           !   k_bins(i) = 0.01d0 + dble(i)*(pi/this%dx - 0.01d0)/Nk
           ! end do 
          !  if (rank.eq.master) then
          !    call save_k(k_bins, Pk, 'MPS')
          !  end if                                 
        end if
    
      !  call this%density()

    
        this%Mass_scalar = integrate(this%rho*this%sqrtg,this%dx,this%dy,this%dz, this%g_pts)
    
    
    end subroutine diagnostic
    !============================================================================!
    subroutine currents_ADM(this)
        implicit none 
        class(scalar_field), intent(in out) :: this 

        this%j_i(1,:,:,:) = -this%u(5,:,:,:)*(this%gamma_up(1,1,:,:,:)*this%u(2,:,:,:) &
                                            + this%gamma_up(1,3,:,:,:)*this%u(4,:,:,:))
        this%j_i(2,:,:,:) = -this%u(5,:,:,:)*this%gamma_up(2,2,:,:,:)*this%u(3,:,:,:)
        this%j_i(3,:,:,:) = -this%u(5,:,:,:)*(this%gamma_up(1,3,:,:,:)*this%u(2,:,:,:) &
                                            + this%gamma_up(3,3,:,:,:)*this%u(4,:,:,:))

    end subroutine
    !============================================================================!
    subroutine proyector_detectors(this)
        implicit none
        class(scalar_field), intent(in out) :: this 

        real(kind=8), dimension(this%iL-this%g_pts:this%iR+this%g_pts) :: r_aux
        integer :: l


        call domain_to_sph(this%u(1,:,:,:), this%x, this%y, this%z,&
                           this%u_inter(1,:,:,:),this%r,this%theta,this%phi,this%index_bounds)


          do l=1, this%n_det  
              if((this%r(1,1,1).le.this%r_det(l)).and.(this%r(this%Nx,1,1).gt.this%r_det(l))) then 
                this%alm(l,:,:) = project_sh(this%u_inter(1,this%index_det(l),:,:) &
                                            ,this%ylm,l_max,this%theta(1,:,1),this%dy,this%dz)
              end if
          end do
      
    end subroutine proyector_detectors

    subroutine currents_detectors(this)
        use integral
        implicit none 
        class(scalar_field), intent(in out) :: this 
        complex(kind=8), dimension(3,this%Nx,this%Ny,this%Nz) :: integrand
        integer :: l

        call domain_to_sph(this%j_i(1,:,:,:)*this%sqrtg(:,:,:), this%x, this%y, this%z,&
                           integrand(1,:,:,:),this%r,this%theta,this%phi,this%index_bounds)
        call domain_to_sph(this%j_i(2,:,:,:)*this%sqrtg(:,:,:), this%x, this%y, this%z,&
                           integrand(2,:,:,:),this%r,this%theta,this%phi,this%index_bounds)        
        call domain_to_sph(this%j_i(3,:,:,:)*this%sqrtg(:,:,:), this%x, this%y, this%z,&
                           integrand(3,:,:,:),this%r,this%theta,this%phi,this%index_bounds)

        do l=1, this%n_det
            if((this%r(1,1,1).le.this%r_det(l)).and.(this%r(this%Nx,1,1).gt.this%r_det(l))) then 
            this%momentum(1,l) = trapezium_2D(integrand(1,this%index_det(l),:,:),this%dy,this%dz)
            this%momentum(2,l) = trapezium_2D(integrand(2,this%index_det(l),:,:),this%dy,this%dz)
            this%momentum(3,l) = trapezium_2D(integrand(3,this%index_det(l),:,:),this%dy,this%dz)
            end if
        end do

    end subroutine

    !============================================================================!
    subroutine density(this)
    implicit none 
    class(scalar_field), intent(in out) :: this 

    this%rho = 0.5d0*((abs(this%u(5,:,:,:))) &
              + this%gamma_up(1,1,:,:,:)*abs(this%u(2,:,:,:))**2 &
              + 2.0d0*this%gamma_up(1,3,:,:,:)*abs(this%u(2,:,:,:))*this%u(4,:,:,:) &
              + this%gamma_up(2,2,:,:,:)*abs(this%u(3,:,:,:))**2 & 
              + this%gamma_up(3,3,:,:,:)*abs(this%u(4,:,:,:))**2 &
              + this%m**2 * abs(this%u(1,:,:,:)) & 
              )

    end subroutine density
    !============================================================================!
      
end module diagnostic_lib
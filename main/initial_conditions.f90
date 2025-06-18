subroutine initial_conditions
  use global_numbers
  use FFT_sph
  use random
  use mpi_lib
  use integral
  use finite_differences, &
  dx1 => first_derivative_x_2, &
  dx2 => first_derivative_y_2, &
  dx3 => first_derivative_z_2

  implicit none
  complex(kind=8), dimension(MoL%iL - MoL%g_pts:MoL%iR + MoL%g_pts, &
                          MoL%jL - MoL%g_pts :  MoL%jR + MoL%g_pts , & 
                          MoL%kL:MoL%kR) :: partial_t_phi
  real(kind=8) :: random_real, random_img
  real(kind=8) :: Amp
  complex(kind=8), dimension(0:l_max_intial,-l_max_intial:l_max_intial) :: coeff
  complex(kind=8) :: f_spectrum(Nk,0:l_max,-l_max:l_max)
  real(kind=8) :: rand_phase(0:Nk,0:l_max,-l_max:l_max), p0
  integer :: i,n,l

  if(initial_conditions_type.eq.'spherical_shell') then 
    if(MoL%metric == 'MKS') then 
      MoL%u(1,:,:,:) = 1.0d0/dexp(MoL%x) * dexp(-(dexp(MoL%x) - mu)**2/sigma**2)
      MoL%u(1,:,:,:) = Amp*MoL%u(1,:,:,:)
      MoL%u(2,:,:,:) = dx1(MoL%u(1,:,:,:),MoL%dx, MoL%index_bounds)
      MoL%u(3,:,:,:) = 0.0d0
      MoL%u(4,:,:,:) = 0.0d0
      partial_t_phi = MoL%u(2,:,:,:)
      MoL%u(5,:,:,:) = MoL%sqrtg / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) ))

      call MoL%density()

      Amp = 1.0d0/sqrt(integrate(MoL%rho*MoL%alpha*MoL%sqrtg,MoL%dx,MoL%dy,MoL%dz, MoL%g_pts))

      MoL%u(:,:,:,:) = Amp*MoL%u(:,:,:,:)

    else 

      MoL%u(1,:,:,:) = 0.0d0
    ! Solo el proceso maestro inicializa el arreglo con números aleatorios
      if (rank == master) then
        do l=0, l_max_intial
            do n=-l, l
                call random_number(random_real)
                call random_number(random_img)
                coeff(l,n) = complex(random_real, random_img)
            end do
        end do
      end if

    ! Sincronizar el arreglo en todos los procesos
    call MPI_Bcast(coeff, size(coeff), MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)

    MoL%u(1,:,:,:) = 0.0d0

      do l=0, l_max_intial
        do n=-l, l
          do i=MoL%iL-MoL%g_pts, MoL%iR + MoL%g_pts
            MoL%u(1,i,:,:) = 1.0d0/MoL%x(i,:,:) * dexp(-(MoL%x(i,:,:) - mu)**2/sigma**2)! * &
                              !(coeff(l,n)*MoL%ylm(l,n,:,:)) + MoL%u(1,i,:,:)
          end do
        end do
      end do

    !  MoL%u(2,:,:,:) = dx1(MoL%u(1,:,:,:),MoL%dx,MoL%index_bounds)
      
      MoL%u(3,:,:,:) = dx2(MoL%u(1,:,:,:),MoL%dy,MoL%index_bounds,.false.)

      MoL%u(4,:,:,:) = dx3(MoL%u(1,:,:,:),MoL%dz,MoL%index_bounds)

      partial_t_phi = -2.0d0 * (MoL%x-mu)/(MoL%x*sigma**2) *  dexp(-(MoL%x-mu)**2/sigma**2)
!      partial_t_phi = MoL%u(2,:,:,:)
      
      MoL%u(5,:,:,:) = 1.0d0 / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) ))
     
     ! call MoL%boundary_zL()
     ! call MoL%boundary_zR()
     !call MoL%boundary_yL()
     !call MoL%boundary_yR()

      call MoL%density()

      Amp = 1.0d0/sqrt(integrate(MoL%rho*MoL%alpha*MoL%sqrtg,MoL%dx,MoL%dy,MoL%dz, MoL%g_pts))

      MoL%u(:,:,:,:) = Amp*MoL%u(:,:,:,:)


    end if

  else if(initial_conditions_type.eq.'random') then
        
      p0 = 1.0d0

      call convert_to_spherical(MoL,p0)


    
      ! call sph_to_domain(ISphFT_discrete_radial(f_spectrum, &
      !                   MoL%r(:,1,1), MoL%theta(1,:,1), MoL%phi(1,1,:), & 
      !                   MoL%ylm, MoL%j_basis, l_max, Nk), &
      !                   MoL%r, MoL%theta, MoL%phi, &
      !                   MoL%u(1,:,:,:), &
      !                   MoL%x, MoL%y, MoL%z, MoL%index_bounds)
      


      if(rank /=0) then
        call MoL%send_recv_xL()
      else
        call MoL%boundary_xL()
      end if
  
      if(rank /=nproc-1) then
        call MoL%send_recv_xR()
      else
        call MoL%boundary_xR()
      end if

      ! call MoL%boundary_zL()
      ! call MoL%boundary_zR()
      ! call MoL%boundary_yL()
      ! call MoL%boundary_yR()


      MoL%u(2,:,:,:) = dx1(MoL%u(1,:,:,:),MoL%dx,MoL%index_bounds)
        
      MoL%u(3,:,:,:) = dx2(MoL%u(1,:,:,:),MoL%dy,MoL%index_bounds,.false.)

      MoL%u(4,:,:,:) = dx3(MoL%u(1,:,:,:),MoL%dz,MoL%index_bounds)

      partial_t_phi = 0.0d0
      
      MoL%u(5,:,:,:) = 1.0d0 / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) ))


      call MoL%density()

      Amp = 1.0d0/sqrt(integrate(MoL%m*abs(MoL%u(1,:,:,:))**2 *MoL%sqrtg,MoL%dx,MoL%dy,MoL%dz, MoL%g_pts))

      MoL%u(:,:,:,:) = Amp*MoL%u(:,:,:,:)

  else if(initial_conditions_type.eq.'random_shell') then
    
     p0 = 1.0d0

     call convert_to_spherical(MoL,p0)


    ! MoL%u(1,:,:,:) = ISphFT_discrete_radial(f_spectrum, &
    !                                         MoL%x(:,1,1), MoL%y(1,:,1), MoL%z(1,1,:), & 
    !                                         MoL%ylm, MoL%j_basis, l_max, Nk)

      MoL%u(1,:,:,:) = MoL%u(1,:,:,:) * 1.0d0/MoL%x(:,:,:) * dexp(-(MoL%x(:,:,:) - mu)**2/sigma**2)

      MoL%u(2,:,:,:) = dx1(MoL%u(1,:,:,:),MoL%dx,MoL%index_bounds)
        
      MoL%u(3,:,:,:) = dx2(MoL%u(1,:,:,:),MoL%dy,MoL%index_bounds,.false.)

      MoL%u(4,:,:,:) = dx3(MoL%u(1,:,:,:),MoL%dz,MoL%index_bounds)

      partial_t_phi = MoL%u(2,:,:,:)
      
      MoL%u(5,:,:,:) = 1.0d0 / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) ))
      
      if(rank /=0) then
        call MoL%send_recv_xL()
      else
        call MoL%boundary_xL()
      end if
  
      if(rank /=nproc-1) then
        call MoL%send_recv_xR()
      else
        call MoL%boundary_xR()
      end if
      call MoL%boundary_zL()
      call MoL%boundary_zR()
      call MoL%boundary_yL()
      call MoL%boundary_yR()

      call MoL%density()

      Amp = 1.0d0/sqrt(integrate(MoL%m*abs(MoL%u(1,:,:,:))**2 *MoL%sqrtg,MoL%dx,MoL%dy,MoL%dz, MoL%g_pts))

      MoL%u(:,:,:,:) = Amp*MoL%u(:,:,:,:)


    else if(initial_conditions_type.eq.'MPS_test') then 

      do i=MoL%iL - MoL%g_pts, MoL%iR + MoL%g_pts
        MoL%u(1,i,:,:) = MoL%j_basis(5,20,i)*MoL%ylm(5,0,:,:) 
      end do

      ! f_spectrum = (0.0d0,0.0d0)

      !  f_spectrum(3,5,0) = (1.0d0,0.0d0)

      !  MoL%u(1,:,:,:) = ISphFT_discrete_radial(f_spectrum, &
      !                  MoL%x(:,1,1), MoL%y(1,:,1), MoL%z(1,1,:), &
      !                  MoL%ylm, MoL%j_basis, l_max, Nk, &
      !                  MoL%Nx, MoL%Ny, MoL%Nz, MoL%g_pts)

      MoL%u(2,:,:,:) = dx1(MoL%u(1,:,:,:),MoL%dx,MoL%index_bounds)
        
      MoL%u(3,:,:,:) = dx2(MoL%u(1,:,:,:),MoL%dy,MoL%index_bounds,.false.)

      MoL%u(4,:,:,:) = dx3(MoL%u(1,:,:,:),MoL%dz,MoL%index_bounds)

      partial_t_phi = MoL%u(2,:,:,:)
      
      MoL%u(5,:,:,:) = 1.0d0 / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) ))

      call MoL%boundary_zL()
      call MoL%boundary_zR()
      call MoL%boundary_yL()
      call MoL%boundary_yR()


  else 
    

    print*, 'Intial condition is not implemented yet'
    stop
    
  end if


end subroutine initial_conditions


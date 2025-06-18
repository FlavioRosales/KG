subroutine evolve
  use global_numbers
  implicit none

  integer :: l

  t_initial = MPI_Wtime()

  if(read_from_checkpoint) then
    if(rank.eq.master) print*, 'Reading the data from the check point...'
    call MoL%read_from_check_point
    if(rank.eq.master) print*, 'Data read successfully.'
  else
    if(rank.eq.master) print*, 'Start to create the initial conditions...'
    call initial_conditions(mu,sigma, k_f)
    call MoL%diagnostic(0)
    if(rank.eq.master) print*, 'End to create the initial conditions.'
  end if

  call save_gnuplot(0)

  do l=1, MoL%Nt

    !call MoL%rk3(rhs = rhs)
    call MoL%pirk2()
    call MoL%diagnostic(l)
    call save_gnuplot(l)

  end do

contains

  !============================================================================!

  ! ! rhs para \phi
  !   subroutine rhs_phi(rhs, u)
  !     real(kind=8), intent(out) :: rhs(:,:,:)
  !     real(kind=8), intent(in)  :: u(:,:,:,:)
  !     rhs = MoL%alpha * u(5,:,:,:) + MoL%beta(1,:,:,:) * u(2,:,:,:)
  !   end subroutine rhs_phi

  !   ! rhs para \psi_1 := \partial_1 \phi
  !   subroutine rhs_psi1(rhs, u)
  !     use finite_differences, &
  !       dx1 => first_derivative_x_2, &
  !       ad_1 => advec_x
  !     implicit none
  !     real(kind=8), intent(out) :: rhs(:,:,:)
  !     real(kind=8), intent(in)  :: u(:,:,:,:)
  !     rhs = dx1(MoL%alpha * u(5,:,:,:), MoL%dx, MoL%index_bounds) &
  !         + MoL%beta(1,:,:,:) * ad_1(u(2,:,:,:), MoL%beta(1,:,:,:), MoL%dx, MoL%index_bounds) &
  !         + MoL%dx_beta * u(2,:,:,:)
  !   end subroutine rhs_psi1

  !   ! rhs para \psi_2 := \partial_2 \phi
  !   subroutine rhs_psi2(rhs, u)
  !     use finite_differences, &
  !     dx2 => first_derivative_y_2
  !     implicit none
  !     real(kind=8), intent(out) :: rhs(:,:,:)
  !     real(kind=8), intent(in)  :: u(:,:,:,:)
  !     rhs = dx2( MoL%alpha * u(5,:,:,:) + MoL%beta(1,:,:,:) * u(2,:,:,:), MoL%dy, MoL%index_bounds, .false.)
  !   end subroutine rhs_psi2

  !   ! rhs para \psi_3 := \partial_3 \phi
  !   subroutine rhs_psi3(rhs, u)
  !     use finite_differences, &
  !     dx3 => first_derivative_z_2
  !     implicit none
  !     real(kind=8), intent(out) :: rhs(:,:,:)
  !     real(kind=8), intent(in)  :: u(:,:,:,:)
  !     rhs = dx3(MoL%alpha * u(5,:,:,:) + MoL%beta(1,:,:,:) * u(2,:,:,:), MoL%dz, MoL%index_bounds)
  !   end subroutine rhs_psi3

  !   ! rhs para \pi
  !   subroutine rhs_pi(rhs, u)
  !     use finite_differences, &
  !     dx1 => first_derivative_x_2, &
  !     dx2 => first_derivative_y_2, &
  !     dx3 => first_derivative_z_2, &
  !     ad_1 => advec_x
  !     implicit none
  !     real(kind=8), intent(out) :: rhs(:,:,:)
  !     real(kind=8), intent(in)  :: u(:,:,:,:)

  !     rhs = dx1(MoL%alpha*MoL%sqrtg*(MoL%gamma_up(1,1,:,:,:)*u(2,:,:,:) &
  !               + MoL%gamma_up(1,3,:,:,:)*u(4,:,:,:)), MoL%dx, MoL%index_bounds)/MoL%sqrtg &
  !         + dx2(MoL%alpha*MoL%sqrtg*MoL%gamma_up(2,2,:,:,:)*u(3,:,:,:) &
  !               , MoL%dy, MoL%index_bounds, .false.)/MoL%sqrtg &
  !         + dx3(MoL%alpha*MoL%sqrtg*(MoL%gamma_up(1,3,:,:,:)*u(2,:,:,:) & 
  !               + MoL%gamma_up(3,3,:,:,:)*u(4,:,:,:)), MoL%dz, MoL%index_bounds)/MoL%sqrtg &
  !         + u(5,:,:,:) * MoL%dx_beta &
  !         + MoL%beta(1,:,:,:) * ad_1(u(5,:,:,:)*MoL%sqrtg, MoL%beta(1,:,:,:), MoL%dx, MoL%index_bounds)/MoL%sqrtg &
  !         - MoL%alpha * MoL%m**2 * u(1,:,:,:)
  !   end subroutine rhs_pi

  ! subroutine rhs(refinement)
  !   use finite_differences, &
  !   dx1 => first_derivative_x_2, &
  !   dx2 => first_derivative_y_2, &
  !   dx3 => first_derivative_z_2, &
  !   ad_1 => advec_x, &
  !   ad_2 => advec_y, &
  !   ad_3 => advec_z
  !   implicit none
  !   integer, intent(in) :: refinement

  !   !
  !   ! r.h.s. for \phi
  !   !
  !   MoL%rhs_u(1,:,:,:) = &
  ! !  + MoL%alpha/MoL%sqrtg * MoL%u(5,:,:,:) + MoL%beta(1,:,:,:)*MoL%u(2,:,:,:)
  !   + MoL%alpha * MoL%u(5,:,:,:) + MoL%beta(1,:,:,:)*MoL%u(2,:,:,:)
  !   !
  !   ! r.h.s. for \psi_1 := \partial_1 \phi
  !   !
  !   MoL%rhs_u(2,:,:,:) = &
  !   + dx1(MoL%alpha * MoL%u(5,:,:,:),MoL%x(:,1,1),MoL%index_bounds) &
  !   + MoL%beta(1,:,:,:)*ad_1(MoL%u(2,:,:,:),MoL%beta(1,:,:,:),MoL%x(:,1,1),MoL%index_bounds) & 
  !   + MoL%dx_beta*MoL%u(2,:,:,:)
  !   ! r.h.s. for \psi_2 := \partial_2 \phi
  !   !
  !   MoL%rhs_u(3,:,:,:) = dx2(MoL%rhs_u(1,:,:,:),MoL%dy,MoL%index_bounds,.false.)
  !   !
  !   ! r.h.s. for \psi_3 := \partial_3 \phi
  !   !
  !   MoL%rhs_u(4,:,:,:) =  dx3(MoL%rhs_u(1,:,:,:),MoL%dz,MoL%index_bounds)
  !   ! 
  !   ! r.h.s. for \pi := partial_i (\pi beta^i + \alpha \sqrt(gamma) \gamma_ij \psi_j) + \alpha \sqrt(gamma)m^2 \phi
  !   !
  !   MoL%rhs_u(5,:,:,:) = &
  !   +  dx1(MoL%alpha*MoL%sqrtg*(MoL%gamma_up(1,1,:,:,:)*MoL%u(2,:,:,:) + MoL%gamma_up(1,3,:,:,:)*MoL%u(4,:,:,:)), &
  !       MoL%x(:,1,1),MoL%index_bounds)/MoL%sqrtg & 
  !   +  dx2(MoL%alpha*MoL%sqrtg*MoL%gamma_up(2,2,:,:,:)*MoL%u(3,:,:,:),MoL%dy,MoL%index_bounds,.false.)/MoL%sqrtg & 
  !   +  dx3(MoL%alpha*MoL%sqrtg*(MoL%gamma_up(1,3,:,:,:)*MoL%u(2,:,:,:) + MoL%gamma_up(3,3,:,:,:)*MoL%u(4,:,:,:))/MoL%sqrtg, &
  !       MoL%dz,MoL%index_bounds) & 
  !   + MoL%u(5,:,:,:)*MoL%dx_beta  &
  !   + MoL%beta(1,:,:,:)*ad_1(MoL%u(5,:,:,:)*MoL%sqrtg,MoL%beta(1,:,:,:),MoL%x(:,1,1), MoL%index_bounds)/MoL%sqrtg &
  !    - MoL%alpha*MoL%m**2 * MoL%u(1,:,:,:)
    

  ! end subroutine rhs
  !============================================================================!

end subroutine evolve

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
    if(rank.eq.master) print*, 'End to create the initial conditions.'
  end if

  call save_gnuplot(0)

  do l=1, MoL%Nt

    call MoL%rk3(rhs = rhs)
    call MoL%diagnostic()
    call save_gnuplot(l)

  end do

contains

  !============================================================================!
  subroutine rhs(refinement)
    use finite_differences, &
    dx1 => first_derivative_x_2, &
    dx2 => first_derivative_y_2, &
    dx3 => first_derivative_z_2, &
    ad_1 => advec_x, &
    ad_2 => advec_y, &
    ad_3 => advec_z
    implicit none
    integer, intent(in) :: refinement


    !
    ! r.h.s. for \phi
    !
    MoL%rhs_u(1,:,:,:) = &
    MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%beta(2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%beta(3,:,:,:) * MoL%u(4,:,:,:) + &
    MoL%alpha * MoL%u(5,:,:,:) / MoL%sqrtg
    !
    ! r.h.s. for \psi_1 := \partial_1 \phi
    !
    !MoL%rhs_u(2,:,:,:) = dx1(MoL%rhs_u(1,:,:,:), MoL%dx)
    MoL%rhs_u(2,:,:,:) = dx1(MoL%alpha/MoL%sqrtg * MoL%u(5,:,:,:), MoL%dx) & 
                       + dx1(MoL%beta(1,:,:,:), MoL%dx)*MoL%u(2,:,:,:) & 
                       + ad_1(MoL%u(2,:,:,:),MoL%beta(1,:,:,:),MoL%dx)*MoL%beta(1,:,:,:)
    !
    ! r.h.s. for \psi_2 := \partial_2 \phi
    !
    !MoL%rhs_u(3,:,:,:) = dx2(MoL%rhs_u(1,:,:,:), MoL%dy)
    MoL%rhs_u(3,:,:,:) = dx2(MoL%alpha/MoL%sqrtg * MoL%u(5,:,:,:), MoL%dy) & 
    + dx2(MoL%beta(2,:,:,:), MoL%dy)*MoL%u(3,:,:,:) & 
    + ad_2(MoL%u(3,:,:,:),MoL%beta(2,:,:,:),MoL%dy)*MoL%beta(2,:,:,:)
    !
    ! r.h.s. for \psi_3 := \partial_3 \phi
    !
    !MoL%rhs_u(4,:,:,:) = dx3(MoL%rhs_u(1,:,:,:), MoL%dz)
    MoL%rhs_u(4,:,:,:) = dx3(MoL%alpha/MoL%sqrtg * MoL%u(5,:,:,:), MoL%dz) & 
    + dx3(MoL%beta(3,:,:,:), MoL%dz)*MoL%u(4,:,:,:) & 
    + ad_3(MoL%u(4,:,:,:),MoL%beta(3,:,:,:),MoL%dz)*MoL%beta(3,:,:,:)
    !
    ! r.h.s. for \pi := partial_i (\pi beta^i + \alpha \sqrt(gamma) \gamma_ij \psi_j) + \alpha \sqrt(gamma)m^2 \phi
    !

    MoL%rhs_u(5,:,:,:) = &
    MoL%u(5,:,:,:)*dx1(MoL%beta(1,:,:,:), MoL%dx) + MoL%beta(1,:,:,:)*ad_1(MoL%u(5,:,:,:),MoL%beta(1,:,:,:),MoL%dx) + &
    MoL%u(5,:,:,:)*dx2(MoL%beta(2,:,:,:), MoL%dy) + MoL%beta(2,:,:,:)*ad_2(MoL%u(5,:,:,:),MoL%beta(2,:,:,:),MoL%dy) + &
    MoL%u(5,:,:,:)*dx3(MoL%beta(3,:,:,:), MoL%dz) + MoL%beta(3,:,:,:)*ad_3(MoL%u(5,:,:,:),MoL%beta(3,:,:,:),MoL%dz) + &
    dx1(MoL%alpha * MoL%sqrtg * (&
    MoL%gamma_up(1,1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%gamma_up(1,2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%gamma_up(1,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dx) + &
    dx2(MoL%alpha * MoL%sqrtg * (&
    MoL%gamma_up(2,1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%gamma_up(2,2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%gamma_up(2,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dy) + &
    dx3(MoL%alpha * MoL%sqrtg * (&
    MoL%gamma_up(3,1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%gamma_up(3,2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%gamma_up(3,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dz) - &
    MoL%alpha * MoL%sqrtg * MoL%m**2 * MoL%u(1,:,:,:)

    ! dx1(MoL%beta(1,:,:,:) * MoL%u(5,:,:,:) + MoL%alpha * MoL%sqrtg * (&
    ! MoL%gamma_up(1,1,:,:,:) * MoL%u(2,:,:,:) + &
    ! MoL%gamma_up(1,2,:,:,:) * MoL%u(3,:,:,:) + &
    ! MoL%gamma_up(1,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dx) + &
    ! dx2(MoL%beta(2,:,:,:) * MoL%u(5,:,:,:) + MoL%alpha * MoL%sqrtg * (&
    ! MoL%gamma_up(2,1,:,:,:) * MoL%u(2,:,:,:) + &
    ! MoL%gamma_up(2,2,:,:,:) * MoL%u(3,:,:,:) + &
    ! MoL%gamma_up(2,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dy) + &
    ! dx3(MoL%beta(3,:,:,:) * MoL%u(5,:,:,:) + MoL%alpha * MoL%sqrtg * (&
    ! MoL%gamma_up(3,1,:,:,:) * MoL%u(2,:,:,:) + &
    ! MoL%gamma_up(3,2,:,:,:) * MoL%u(3,:,:,:) + &
    ! MoL%gamma_up(3,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dz) - &
    ! MoL%alpha * MoL%sqrtg * MoL%m**2 * MoL%u(1,:,:,:)
    

  end subroutine rhs
  !============================================================================!

end subroutine evolve

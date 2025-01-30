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
    dx1 => first_derivative_x_4, &
    dx2 => first_derivative_y_4, &
    dx3 => first_derivative_z_4
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
    MoL%rhs_u(2,:,:,:) = dx1(MoL%rhs_u(1,:,:,:), MoL%dx)
    !
    ! r.h.s. for \psi_2 := \partial_2 \phi
    !
    MoL%rhs_u(3,:,:,:) = dx2(MoL%rhs_u(1,:,:,:), MoL%dy)
    !
    ! r.h.s. for \psi_3 := \partial_3 \phi
    !
    MoL%rhs_u(4,:,:,:) = dx3(MoL%rhs_u(1,:,:,:), MoL%dz)
    !
    ! r.h.s. for \pi := partial_i (\pi beta^i + \alpha \sqrt(gamma) \gamma_ij \psi_j) + \alpha \sqrt(gamma)m^2 \phi
    !

    MoL%rhs_u(5,:,:,:) = &
    dx1(MoL%beta(1,:,:,:) * MoL%u(5,:,:,:) + MoL%alpha * MoL%sqrtg * (&
    MoL%gamma_up(1,1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%gamma_up(1,2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%gamma_up(1,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dx) + &
    dx2(MoL%beta(2,:,:,:) * MoL%u(5,:,:,:) + MoL%alpha * MoL%sqrtg * (&
    MoL%gamma_up(2,1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%gamma_up(2,2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%gamma_up(2,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dy) + &
    dx3(MoL%beta(3,:,:,:) * MoL%u(5,:,:,:) + MoL%alpha * MoL%sqrtg * (&
    MoL%gamma_up(3,1,:,:,:) * MoL%u(2,:,:,:) + &
    MoL%gamma_up(3,2,:,:,:) * MoL%u(3,:,:,:) + &
    MoL%gamma_up(3,3,:,:,:) * MoL%u(4,:,:,:) ), MoL%dz) + &
    MoL%alpha * MoL%sqrtg * MoL%m**2 * MoL%u(1,:,:,:)**2
    

  end subroutine rhs
  !============================================================================!

end subroutine evolve

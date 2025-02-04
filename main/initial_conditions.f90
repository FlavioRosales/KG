subroutine initial_conditions
  use global_numbers
  use integral
  implicit none
  
  real(kind=8), dimension(MoL%iL - MoL%g_pts:MoL%iR + MoL%g_pts, &
                          MoL%jL - MoL%g_pts:MoL%jR + MoL%g_pts, & 
                          MoL%kL:MoL%kR) :: partial_t_phi
  real(kind=8) :: Amp

  ! MoL%u(1,:,:,:) = exp(-(MoL%x**2 + 10.0d0**2 -2.0d0*10.0d0*MoL%x*dcos(MoL%y))/5.0d0)
  ! MoL%u(2,:,:,:) = MoL%U(1,:,:,:)*(-2.0d0*MoL%x + 2.0d0*10.0d0*dcos(MoL%y))/5.0d0
  ! MoL%u(3,:,:,:) = MoL%u(1,:,:,:)*(-2.0d0*10.0d0*MoL%x*dsin(MoL%y))/5.0d0
  ! MoL%u(4,:,:,:) = 0.0d0
  ! MoL%u(5,:,:,:) = MoL%sqrtg/MoL%alpha * ( MoL%beta(1,:,:,:)*MoL%u(2,:,:,:))

  !MoL%u(1,:,:,:) = exp(-(MoL%x**2 + 10.0d0**2 -2.0d0*10.0d0*MoL%x*dsin(MoL%y)*dsin(MoL%z))/5.0d0)
  !MoL%u(2,:,:,:) = MoL%u(1,:,:,:)*(-2.0d0*MoL%x + 2.0d0*10.0d0*dsin(MoL%y)*dsin(MoL%z))/5.0d0
  !MoL%u(3,:,:,:) = MoL%u(1,:,:,:)*(2.0d0*10.0d0*MoL%x*dcos(MoL%y)*dsin(MoL%z))/5.0d0
  !MoL%u(4,:,:,:) = MoL%u(1,:,:,:)*(2.0d0*10.0d0*MoL%x*dsin(MoL%y)*dcos(MoL%z))/5.0d0
  !MoL%u(5,:,:,:) = MoL%sqrtg/MoL%alpha * ( MoL%beta(1,:,:,:)*MoL%u(2,:,:,:))

  !  print*, 'v_max = ', pi/max(this%dx,this%x(this%iR + this%g_pts)*this%dy,this%x(this%iR + this%g_pts)*this%dz)

  if(initial_conditions_type.eq.'spherical_shell') then 
    if(MoL%metric == 'MKS') then 
      MoL%u(1,:,:,:) = dexp(-0.5d0*(exp(MoL%x) - mu)**2 / sigma**2)/exp(MoL%x) * dcos(k_f*exp(MoL%x))
      Amp = 0.25d0/integrate(MoL%u(1,:,:,:)*MoL%alpha*MoL%sqrtg,MoL%dx,MoL%dy,MoL%dz, MoL%g_pts)
      MoL%u(1,:,:,:) = Amp*MoL%u(1,:,:,:)
      MoL%u(2,:,:,:) =Amp*( -MoL%u(1,:,:,:)/exp(MoL%x) - MoL%u(1,:,:,:) * (exp(MoL%x) - mu)/sigma**2 &
                      -dexp(-0.5d0*(exp(MoL%x) - mu)**2 / sigma**2)/exp(MoL%x) * k_f * dsin(k_f*exp(MoL%x)))*exp(MoL%x)
      MoL%u(3,:,:,:) = 0.0d0
      MoL%u(4,:,:,:) = 0.0d0
      partial_t_phi = &!-MoL%beta(1,:,:,:)/MoL%alpha * MoL%u(1,:,:,:)
              !-(exp(MoL%x) - mu)/sigma**2 * MoL%u(1,:,:,:)  &
              !-  dexp(-0.5d0*(exp(MoL%x) - mu)**2 / sigma**2)/exp(MoL%x) * k_f * dsin(k_f*exp(MoL%x))
            -Amp*MoL%u(1,:,:,:) * (exp(MoL%x) - mu)/sigma**2 
            !MoL%beta(1,:,:,:) * MoL%u(2,:,:,:)
      MoL%u(5,:,:,:) = MoL%sqrtg / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) + &
                      MoL%beta(2,:,:,:) * MoL%u(3,:,:,:) +  &
                      MoL%beta(3,:,:,:) * MoL%u(4,:,:,:)))
    else 
      MoL%u(1,:,:,:) = dexp(-0.5d0*(MoL%x - mu)**2 / sigma**2)/MoL%x * dcos(k_f*MoL%x)
      Amp = 0.25d0/integrate(MoL%u(1,:,:,:)*MoL%alpha*MoL%sqrtg,MoL%dx,MoL%dy,MoL%dz, MoL%g_pts)
      MoL%u(2,:,:,:) = (-MoL%u(1,:,:,:)/MoL%x - MoL%u(1,:,:,:) * (MoL%x - mu)/sigma**2 &
                      -dexp(-0.5d0*(MoL%x - mu)**2 / sigma**2)/MoL%x * k_f * dsin(k_f*MoL%x))*Amp
      MoL%u(3,:,:,:) = 0.0d0
      MoL%u(4,:,:,:) = 0.0d0
      partial_t_phi = &!-MoL%beta(1,:,:,:)/MoL%alpha * MoL%u(1,:,:,:)
              !-(MoL%x - mu)/sigma**2 * MoL%u(1,:,:,:)  &
              !-  dexp(-0.5d0*(MoL%x - mu)**2 / sigma**2)/MoL%x * k_f * dsin(k_f*MoL%x)
            -Amp*MoL%u(1,:,:,:) * (MoL%x - mu)/sigma**2 
            !MoL%beta(1,:,:,:) * MoL%u(2,:,:,:)
      MoL%u(5,:,:,:) = MoL%sqrtg / MoL%alpha * (partial_t_phi - (&
                      MoL%beta(1,:,:,:) * MoL%u(2,:,:,:) + &
                      MoL%beta(2,:,:,:) * MoL%u(3,:,:,:) +  &
                      MoL%beta(3,:,:,:) * MoL%u(4,:,:,:)))
    end if
  else 
    print*, 'Intial condition is not implemented yet'
    stop
  end if

end subroutine initial_conditions
subroutine save_gnuplot(index)
  use save_data_lib
  use global_numbers
  implicit none
  integer, intent(in) :: index
  integer :: i
  real(kind=8), dimension(MoL%iL-MoL%g_pts:MoL%iR+MoL%g_pts) :: r_aux
  character(len=1) :: num_det
  character(len=20) :: name_detector

  r_aux = MoL%x(:,1,1)
  if(MoL%metric == 'MKS') r_aux = dexp(r_aux)


  if((mod(index,save_0d).eq.0).and.(rank.eq.master)) then 
 !   call save_t(t = MoL%t, f=MoL%phi_max(), name = 'phi_max.t')
  !  call save_t(t = MoL%t, f=MoL%phi_min(), name = 'phi_min.t')
    call save_t(t = MoL%t, f=MoL%Mass_scalar, name = 'Mass_scalar.t')

  end if

  if((mod(index,save_0d).eq.0).and.proc_x_save) then
    do i=1, MoL%n_det
      name_detector ='phi_flux'
      if(rank.eq.0) then 
        if((r_aux(MoL%iL - MoL%g_pts).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I1)') i
          call save_t(t = Mol%t, f = MoL%flux(i), name = trim(name_detector) // trim(num_det))
        end if
      else
        if((r_aux(MoL%iL).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I1)') i
          call save_t(t = Mol%t, f = MoL%flux(i), name = trim(name_detector) // trim(num_det))
        end if
      end if

      name_detector ='energy_current'
      if(rank.eq.0) then 
        if((r_aux(MoL%iL - MoL%g_pts).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I1)') i
          call save_t(t = Mol%t, f = MoL%current(i), name = trim(name_detector) // trim(num_det))
        end if
      else
        if((r_aux(MoL%iL).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I1)') i
          call save_t(t = Mol%t, f = MoL%current(i), name = trim(name_detector) // trim(num_det))
        end if
      end if
    end do
  end if

  if(mod(index,save_1d).eq.0 .and. proc_y_save .and. proc_z_save) then
    if(rank.eq.master) then
    print*, '|====================================================================|'
    print*, '|    index     |          time          |        cpu time            |'
    print*, '|====================================================================|'
    end if
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%phi(iL:iR,j_save,k_save), name = 'phi')
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%W(iL:iR,j_save,k_save), name = 'W')
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%sqrtg(iL:iR,j_save,k_save), name = 'e')
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%p(iL:iR,j_save,k_save), name = 'p')
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%vx(iL:iR,j_save,k_save), name = 'vx')
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%vy(iL:iR,j_save,k_save), name = 'vy')
    ! call save_x(axis = MoL%x(iL:iR,j_save,k_save), f = gas%vz(iL:iR,j_save,k_save), name = 'vz')

  end if


  if(mod(index,save_1d).eq.0 .and. proc_x_save .and. proc_z_save) then

    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%phi(i_save,jL:jR,k_save), name = 'phi')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%p(i_save,jL:jR,k_save), name = 'p')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%vx(i_save,jL:jR,k_save), name = 'vx')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%vy(i_save,jL:jR,k_save), name = 'vy')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%vz(i_save,jL:jR,k_save), name = 'vz')
  end if


  if(mod(index,save_1d).eq.0 .and. proc_x_save .and. proc_y_save) then

    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR), f = gas%phi(i_save,j_save,kL:kR), name = 'phi')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR), f = gas%p(i_save,j_save,kL:kR), name = 'p')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR), f = gas%vx(i_save,j_save,kL:kR), name = 'vx')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR), f = gas%vy(i_save,j_save,kL:kR), name = 'vy')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR), f = gas%vz(i_save,j_save,kL:kR), name = 'vz')

  end if

 if(mod(index,save_2d).eq.0) then
   ! call save_xz(&
   ! x = MoL%x(iL:iR,jL:jR,1)*dsin(MoL%y(iL:iR,jL:jR,1))*cos(MoL%z(iL:iR,jL:jR,1)), &
   ! z = MoL%x(iL:iR,jL:jR,1)*dcos(MoL%y(iL:iR,jL:jR,1)), &
   !  f = gas%phi(iL:iR,jL:jR,1), name = 'phi')
 end if



  if(mod(index, save_3d).eq.0) then

      call write_hdf5('phi.h5', &
      '/refinement_1' , 'phi_' // char_integer(index/save_3d), MoL%u(1,iL-g_pts:iR+g_pts,jL:jR,kL:kR))

      call write_hdf5('phi2.h5', &
      '/refinement_1' , 'phi2_' // char_integer(index/save_3d), MoL%u(1,iL-g_pts:iR+g_pts,jL:jR,kL:kR)**2)

      ! call write_hdf5('psi_x.h5', &
      ! '/refinement_1' , 'psi_x_' // char_integer(index/save_3d), MoL%u(2,iL-g_pts:iR+g_pts,jL:jR,kL:kR))

      ! call write_hdf5('psi_y.h5', &
      ! '/refinement_1' , 'psi_y_' // char_integer(index/save_3d), MoL%u(3,iL-g_pts:iR+g_pts,jL:jR,kL:kR))

      ! call write_hdf5('psi_z.h5', &
      ! '/refinement_1' , 'psi_z_' // char_integer(index/save_3d), MoL%u(4,iL-g_pts:iR+g_pts,jL:jR,kL:kR))

      ! call write_hdf5('pi.h5', &
      ! '/refinement_1' , 'pi_' // char_integer(index/save_3d), MoL%u(5,iL-g_pts:iR+g_pts,jL:jR,kL:kR))
  end if

  t_final = MPI_Wtime() - t_initial


  if(rank.eq.master) write(*,*) '|', index, MoL%t, t_final, '   |'

end subroutine save_gnuplot

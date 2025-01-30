subroutine save_gnuplot(index)
  use save_data_lib
  use global_numbers
  implicit none
  integer, intent(in) :: index
  integer :: i
  character(len=1) :: num_det
  character(len=10) :: name_detector

  name_detector ='phi_flux'

  if((mod(index,save_0d).eq.0).and.proc_x_save) then
     do i=1, n_det
       write(num_det, '(I1)') i
       call save_t(t = Mol%t, f = MoL%flux(i), name = trim(name_detector) // trim(num_det))
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
  end if

  t_final = MPI_Wtime() - t_initial


  if(rank.eq.master) write(*,*) '|', index, MoL%t, t_final, '   |'

end subroutine save_gnuplot

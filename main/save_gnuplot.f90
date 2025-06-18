subroutine save_gnuplot(index)
  use save_data_lib
  use global_numbers
  implicit none
  integer, intent(in) :: index
  integer :: i
  real(kind=8), dimension(MoL%iL-MoL%g_pts:MoL%iR+MoL%g_pts) :: r_aux
  character(len=2) :: num_det
  character(len=20) :: name_detector

  r_aux = MoL%x(:,1,1)
  if(MoL%metric == 'MKS') r_aux = dexp(r_aux)


  if((mod(index,save_0d).eq.0).and.(rank.eq.master)) then 
    call save_t(t = MoL%t, f=MoL%Mass_scalar, name = 'Mass_scalar.t')
  end if


  if((mod(index,save_0d).eq.0).and.proc_x_save) then
    do i=1, MoL%n_det

      name_detector ='jr_'
      if(rank.eq.0) then 
        if((r_aux(MoL%iL - MoL%g_pts).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_t(t = Mol%t, f = MoL%momentum(1,i), name = trim(name_detector) // trim(num_det))
        end if
      else
        if((r_aux(MoL%iL).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_t(t = Mol%t, f = MoL%momentum(1,i), name = trim(name_detector) // trim(num_det))
        end if
      end if

      name_detector ='jtheta_'
      if(rank.eq.0) then 
        if((r_aux(MoL%iL - MoL%g_pts).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_t(t = Mol%t, f = MoL%momentum(2,i), name = trim(name_detector) // trim(num_det))
        end if
      else
        if((r_aux(MoL%iL).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_t(t = Mol%t, f = MoL%momentum(2,i), name = trim(name_detector) // trim(num_det))
        end if
      end if

      name_detector ='jphi_'
      if(rank.eq.0) then 
        if((r_aux(MoL%iL - MoL%g_pts).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_t(t = Mol%t, f = MoL%momentum(3,i), name = trim(name_detector) // trim(num_det))
        end if
      else
        if((r_aux(MoL%iL).le.MoL%r_det(i)).and.(r_aux(MoL%iR).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_t(t = Mol%t, f = MoL%momentum(3,i), name = trim(name_detector) // trim(num_det))
        end if
      end if

    end do
  end if

  if((mod(index,MoL%save).eq.0).and.proc_x_save) then
    do i=1, MoL%n_det
        if((MoL%r(1,1,1).le.MoL%r_det(i)).and.(MoL%r(MoL%Nx,1,1).gt.MoL%r_det(i))) then 
          write(num_det, '(I2.2)') i
          call save_a(f = MoL%alm(i,:,:),name = 'Harmonic_coeficients'//trim(num_det),l_max=Ubound(MoL%alm,2) )
        end if
    end do
  end if


  if(mod(index,save_1d).eq.0 .and. proc_y_save .and. proc_z_save) then
    if(rank.eq.master) then
    print*, '|====================================================================|'
    print*, '|    index     |          time          |        cpu time            |'
    print*, '|====================================================================|'
    end if
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%phi(MoL%index_det(1):iR,j_save,k_save), name = 'phi')
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%W(MoL%index_det(1):iR,j_save,k_save), name = 'W')
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%sqrtg(MoL%index_det(1):iR,j_save,k_save), name = 'e')
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%p(MoL%index_det(1):iR,j_save,k_save), name = 'p')
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%vx(MoL%index_det(1):iR,j_save,k_save), name = 'vx')
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%vy(MoL%index_det(1):iR,j_save,k_save), name = 'vy')
    ! call save_x(axis = MoL%x(MoL%index_det(1):iR,j_save,k_save), f = gas%vz(MoL%index_det(1):iR,j_save,k_save), name = 'vz')

  end if


  if(mod(index,save_1d).eq.0 .and. proc_x_save .and. proc_z_save) then

    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%phi(i_save,jL:jR,k_save), name = 'phi')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%p(i_save,jL:jR,k_save), name = 'p')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%vx(i_save,jL:jR,k_save), name = 'vx')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%vy(i_save,jL:jR,k_save), name = 'vy')
    ! call save_y(axis = MoL%y(i_save,jL:jR,k_save), f = gas%vz(i_save,jL:jR,k_save), name = 'vz')
  end if


  if(mod(index,save_1d).eq.0 .and. proc_x_save .and. proc_y_save) then

    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR-1), f = gas%phi(i_save,j_save,kL:kR-1), name = 'phi')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR-1), f = gas%p(i_save,j_save,kL:kR-1), name = 'p')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR-1), f = gas%vx(i_save,j_save,kL:kR-1), name = 'vx')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR-1), f = gas%vy(i_save,j_save,kL:kR-1), name = 'vy')
    ! call save_z(axis = MoL%z(i_save,j_save,kL:kR-1), f = gas%vz(i_save,j_save,kL:kR-1), name = 'vz')

  end if

 if(mod(index,save_2d).eq.0) then
   ! call save_xz(&
   ! x = MoL%x(MoL%index_det(1):iR,jL:jR,1)*dsin(MoL%y(MoL%index_det(1):iR,jL:jR,1))*cos(MoL%z(MoL%index_det(1):iR,jL:jR,1)), &
   ! z = MoL%x(MoL%index_det(1):iR,jL:jR,1)*dcos(MoL%y(MoL%index_det(1):iR,jL:jR,1)), &
   !  f = gas%phi(MoL%index_det(1):iR,jL:jR,1), name = 'phi')
 end if



  if(mod(index, MoL%save).eq.0) then
    call write_hdf5('phi.h5', &
    '/refinement_1' , 'phi_' // char_integer(index/MoL%save), abs(MoL%u(1,iL-g_pts:iR,jL:jR,kL:kR-1)))

    call write_hdf5('rho.h5', &
    '/refinement_1' , 'rho_' // char_integer(index/MoL%save), MoL%rho(iL-g_pts:iR,jL:jR,kL:kR-1))

    call write_hdf5('psix.h5', &
    '/refinement_1' , 'psix_' // char_integer(index/MoL%save), abs(MoL%u(2,iL-g_pts:iR,jL:jR,kL:kR-1)))

    call write_hdf5('psiy.h5', &
    '/refinement_1' , 'psiy_' // char_integer(index/MoL%save), abs(MoL%u(3,iL-g_pts:iR,jL:jR,kL:kR-1)))

    call write_hdf5('psiz.h5', &
    '/refinement_1' , 'psiz_' // char_integer(index/MoL%save), abs(MoL%u(4,iL-g_pts:iR,jL:jR,kL:kR-1)))

    call write_hdf5('pi.h5', &
    '/refinement_1' , 'pi_' // char_integer(index/MoL%save), abs(MoL%u(5,iL-g_pts:iR,jL:jR,kL:kR-1)))

    call write_hdf5('rePhi.h5', &
    '/refinement_1' , 'rePhi_' // char_integer(index/MoL%save), dreal(MoL%u(1,iL-g_pts:iR,jL:jR,kL:kR-1)))

    call write_hdf5('imPhi.h5', &
    '/refinement_1' , 'imPhi_' // char_integer(index/MoL%save), dimag(MoL%u(1,iL-g_pts:iR,jL:jR,kL:kR-1)))

    ! call write_hdf5('jr.h5', &
    ! '/refinement_1' , 'jr_' // char_integer(index/save_3d), MoL%j_i(1,iL-g_pts:iR,jL:jR,kL:kR-1))

    ! call write_hdf5('jtheta.h5', &
    ! '/refinement_1' , 'jtheta_' // char_integer(index/save_3d), MoL%j_i(2,iL-g_pts:iR,jL:jR,kL:kR-1))

    ! call write_hdf5('jphi.h5', &
    ! '/refinement_1' , 'jphi_' // char_integer(index/save_3d), MoL%j_i(3,iL-g_pts:iR,jL:jR,kL:kR-1))
    if(rank.eq.master) write(*,*) '|', index, MoL%t, t_final, '   |'

  end if

  t_final = MPI_Wtime() - t_initial



end subroutine save_gnuplot

module save_data_lib
  use save_data
  use strings_lib
  use mpi_lib
  implicit none

contains

  subroutine save_t(t, f, name)
    implicit none
    real(kind=8), intent(in) :: t, f
    character(len=*), intent(in) :: name
    call save0d(t, f, name + '.t.asc')

  end subroutine save_t

  subroutine save_k(axis, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: axis, f
    character(len=*), intent(in) :: name
      call save1d_k(axis, f, name + '.k.asc')
  end subroutine save_k

  subroutine save_x(axis, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: axis, f
    character(len=*), intent(in) :: name
      call save1d(axis, f, name + '_'  + str(rank+1) + '.x.asc')
  end subroutine save_x

  subroutine save_y(axis, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: axis, f
    character(len=*), intent(in) :: name

    call save1d(axis, f, name + '_' + str(rank+1) + '.y.asc')

  end subroutine save_y

  subroutine save_a(f, name, l_max)
    implicit none
    integer, intent(in) :: l_max
    complex(kind=8),intent(in) ::  f(0:l_max,-l_max:l_max)
    character(len=*), intent(in) :: name

    call savealm( f, name + '.lm.asc',l_max)

  end subroutine save_a

  subroutine save_z(axis, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: axis, f
    character(len=*), intent(in) :: name

    call save1d(axis, f, name + '_'  + str(rank+1) + '.z.asc')

  end subroutine save_z

  subroutine save_xy(x, y, f, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, y, f
    character(len=*), intent(in) :: name

    call save2d(x, y, f, name + '_' + str(rank+1) + '.xy.asc')

  end subroutine save_xy

  subroutine save_xz(x, z, f, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, z, f
    character(len=*), intent(in) :: name

    call save2d(x, z, f, name + '_' + str(rank+1) + '.xz.asc')

  end subroutine save_xz

  subroutine save_yz(y, z, f, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: y, z, f
    character(len=*), intent(in) :: name

    call save2d(y, z, f, name + '_' + str(rank+1) + '.xz.asc')

  end subroutine save_yz

  subroutine save_vect_xy(x, y, fx, fy, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, y, fx, fy
    character(len=*), intent(in) :: name

    call save2d_vect(x, y, fx, fy, name + '_' + str(rank+1) + '.vect_xy.asc')

  end subroutine save_vect_xy

  subroutine save_vect_xz(x, z, fx, fz, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, z, fx, fz
    character(len=*), intent(in) :: name

    call save2d_vect(x, z, fx, fz, name + '_' + str(rank+1) + '.vect_xz.asc')

  end subroutine save_vect_xz

  subroutine save_vect_yz(y, z, fy, fz, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: y, z, fy, fz
    character(len=*), intent(in) :: name

    call save2d_vect(y, z, fy, fz, name + '_' + str(rank+1) + '.vect_yz.asc')

  end subroutine save_vect_yz

end module save_data_lib

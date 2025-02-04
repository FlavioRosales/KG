module save_data
  use strings_lib
  implicit none

  interface save
    module procedure :: save0d, save1d, save2d, save2d_vect
  end interface save

contains
  !=============================================================================!
  subroutine save0d(t, f, name)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: t, f
    type(string) :: name
    !
    ! internal
    !
    integer :: ios, iounit
    logical :: exist


    exist = .false.
    iounit = 100
    inquire(file = name%string_data, exist = exist)



    if(exist) then

      open(unit=iounit, file=name%string_data, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "

        write(iounit, *) t, f

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    else
      open(unit=iounit, file=name%string_data, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file "

        write(iounit, *) t, f

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    end if

  end subroutine save0d
  !=============================================================================!
  subroutine save1d(x, f, name)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(:), intent(in) :: x, f
    type(string) :: name
    !
    ! internal
    !
    integer :: i, j, ios, iounit
    logical :: exist

    

    iounit = 200
    inquire(file = name%string_data, exist = exist)

    if(exist) then

      open(unit=iounit, file=name%string_data, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) then
        stop "Error opening file "
      end if

      do i=1, size(x)
          write(iounit, *) x(i), f(i)
      end do
        write(iounit, *)
        write(iounit, *)

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    else
      open(unit=iounit, file=name%string_data, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) then
        print*, 'falla archivo ', name%string_data
        stop "Error opening file "
      end if
        do i=1, size(x)
            write(iounit, *) x(i), f(i)
        end do
          write(iounit, *)
          write(iounit, *)

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    end if

  end subroutine save1d
  !=============================================================================!
  subroutine save2d(x, y, f, name)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(:,:), intent(in) :: x, y, f
    type(string) :: name
    !
    ! internal
    !
    integer :: i, j, k, ios, iounit

    logical :: exist
    iounit = 300
    inquire(file = name%string_data, exist = exist)

    if(exist) then

      open(unit=iounit, file=name%string_data, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "

          do i=1, size(x,1)
            do j=1, size(x, 2)
              write(iounit, *) x(i,j), y(i,j), f(i,j)
            end do
            write(iounit, *)
          end do
          write(iounit, *)
          write(iounit, *)

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    else

      open(unit=iounit, file=name%string_data, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file "

        do i=1, size(x,1)
          do j=1, size(x, 2)
            write(iounit, *) x(i,j), y(i,j), f(i,j)
          end do
          write(iounit, *)
        end do
        write(iounit, *)
        write(iounit, *)

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    end if

  end subroutine save2d
  !=============================================================================!
  subroutine save2d_vect(x, y, fx, fy, name)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(:,:), intent(in) :: x, y, fx, fy
    type(string) :: name
    !
    ! internal
    !
    integer :: i, j, k, ios, iounit

    logical :: exist
    iounit = 300
    inquire(file = name%string_data, exist = exist)

    if(exist) then

      open(unit=iounit, file=name%string_data, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "

          do i=1, size(x,1)
            do j=1, size(x, 2)
              write(iounit, *) x(i,j), y(i,j), fx(i,j), fy(i,j)
            end do
            write(iounit, *)
          end do
          write(iounit, *)
          write(iounit, *)

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    else

      open(unit=iounit, file=name%string_data, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file "

        do i=1, size(x,1)
          do j=1, size(x, 2)
            write(iounit, *) x(i,j), y(i,j), fx(i,j), fy(i,j)
          end do
          write(iounit, *)
        end do
        write(iounit, *)
        write(iounit, *)

      close(unit=iounit, iostat=ios)
      if ( ios /= 0 ) stop "Error closing file unit "

    end if

  end subroutine save2d_vect
  !=============================================================================!
end module save_data

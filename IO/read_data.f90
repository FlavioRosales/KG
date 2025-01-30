module read_data_lib
  implicit none

  type, public :: read_data
    character(len=:), allocatable :: filename
    integer :: lines_number, columns_number
    real(kind=8), allocatable, dimension(:,:) :: data

  contains

    procedure :: read

  end type read_data

  interface read_data
    module procedure :: read_data_constructor
  end interface read_data

contains

  function read_data_constructor(filename, columns_number) result(this)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: columns_number
    type(read_data) :: this
    this%columns_number = columns_number

    allocate(character(len = len(filename)) :: this%filename)
    this%filename = filename
    call this%read

  end function read_data_constructor

  subroutine read(this)
    implicit none
    class(read_data), intent(in out) :: this
    integer :: i

    open(100, file = this%filename, action = 'read')
      read(100,*) this%lines_number
      allocate(this%data(this%columns_number, 0:this%lines_number))

      do i=0, this%lines_number
        read(100,*) this%data(:,i)
      end do

    close(100)

  end subroutine read

end module read_data_lib

module strings_lib
  implicit none

  !
  ! create the string variable
  !
  type, public :: string
    character(len=:), allocatable :: string_data
  contains
    procedure :: length
    procedure :: split
    !final :: string_destructor
  end type string
  !
  ! assign a string of characters to a string or vice versa
  !
  interface assignment(=)
    module procedure :: string_eq_character, character_eq_string
  end interface
  !
  ! concatenate
  !
  interface operator(+)
    module procedure :: concat_string_string, concat_string_char, concat_char_string, concat_char_char
  end interface
  !
  ! converts an integer to a character string
  !
  interface str
    module procedure :: str_integer, str_double
  end interface str
  !
  ! remove character
  !
  interface remove
    module procedure :: remove_string, remove_char
  end interface remove

contains
  !============================================================================!
  ! create a destructor
  !============================================================================!
  ! subroutine string_destructor(this)
  !   implicit none
  !   type(string), intent(in out) :: this  !
  !   deallocate(this%string_data)
  !
  ! end subroutine string_destructor
  !============================================================================!
  ! number of characters in the string
  !============================================================================!
  integer function length(this)
    implicit none
    class(string), intent(in out) :: this
    length = len(this%string_data)

  end function length
  !============================================================================!
  ! character to string assignment method
  !============================================================================!
  subroutine string_eq_character(str, char)
    implicit none
    type(string), intent(in out) :: str
    character(len=*), intent(in) :: char

    str%string_data = char

  end subroutine string_eq_character
  !============================================================================!
  ! string to character assignment method
  !============================================================================!
  subroutine character_eq_string(char, str)
    implicit none
    character(len=*), intent(in out) :: char
    type(string), intent(in) :: str

    char = str%string_data

  end subroutine character_eq_string
  !============================================================================!
  ! concatenate string with string
  !============================================================================!
  function concat_string_string(str_1, str_2) result(str)
    implicit none
    type(string), intent(in)  :: str_1, str_2
    type(string):: str

    str%string_data = str_1%string_data // str_2%string_data

  end function concat_string_string
  !============================================================================!
  ! concatenate string with char
  !============================================================================!
  function concat_string_char(str_1, char_2) result(str)
    implicit none
    type(string), intent(in)  :: str_1
    character(len=*), intent(in) :: char_2
    type(string):: str

    str%string_data = str_1%string_data // char_2

  end function concat_string_char
  !============================================================================!
  ! concatenate string with char
  !============================================================================!
  function concat_char_string(char_1, str_2) result(str)
    implicit none
    character(len=*), intent(in) :: char_1
    type(string), intent(in)  :: str_2
    type(string):: str

    str%string_data =  char_1 // str_2%string_data

  end function concat_char_string
  !============================================================================!
  ! concatenate char with char
  !============================================================================!
  function concat_char_char(char_1, char_2) result(str)
    implicit none
    character(len=*), intent(in) :: char_1, char_2
    type(string):: str

    str%string_data =  char_1 // char_2

  end function concat_char_char
  !============================================================================!
  ! convert integer to string
  !============================================================================!
  function str_integer(int) result(result)
    implicit none
    integer, intent(in) :: int
    type(string) :: result
    character(len=:), allocatable :: result_char
    integer :: n

    n = count(int)

    if(n.eq.0) n = 1

    allocate(character(len=n) :: result_char)

    write(result_char, '(i' // char(n) // ')') int

    result = result_char

    deallocate(result_char)

  end function str_integer
  !============================================================================!
  function char_integer(int) result(result_char)
    implicit none
    integer, intent(in) :: int
    character(len=:), allocatable :: result_char
    integer :: n

    n = count(int)

    if(n.eq.0) n = 1

    allocate(character(len=n) :: result_char)

    write(result_char, '(i' // char(n) // ')') int

    result_char = trim(result_char)

  end function char_integer
  !============================================================================!
  function str_double(double) result(result)
    implicit none
    real(kind=8), intent(in) :: double
    type(string) :: result
    character(len=13) :: result_char

    write(result_char, '(F6.3)') double

    result = result_char

  end function str_double
  !============================================================================!
  ! convert the integers 0-9 to character
  !============================================================================!
  character(len=1) function char(int)
    implicit none
    integer, intent(in) :: int

    write(char, '(i1.1)') int

  end function char
  !============================================================================!
  ! count the number of digits in a integer
  !============================================================================!
  integer function count(int)
    implicit none
    integer, intent(in) :: int
    integer :: aux

    aux = int

    count = 0

    do while(aux.ne.0)
      aux = aux / 10
      count = count + 1
    end do

  end function count
  !============================================================================!
  subroutine print(str)
    implicit none
    type(string), intent(in) :: str

    print*, str%string_data

  end subroutine print
  !============================================================================!
  function remove_char(input, char) result(output)
    implicit none
    character(len=*), intent(in) :: input
    character(len=*), intent(in) :: char
    type(string) :: output
    integer :: i
    character(len=1) :: input_array(len(input))

    output = ''
    do i=1, len(input)
      if(input(i:i).ne.char) output = output + input(i:i)
    end do

  end function remove_char
  !============================================================================!
  function remove_string(input, char) result(output)
    implicit none
    character(len=1), intent(in) :: char
    type(string) :: input, output
    integer :: i
    character(len=1) :: input_array(len(input%string_data))

    output = ''
    do i=1, len(input%string_data)
      if(input%string_data(i:i).ne.char) output = output + input%string_data(i:i)
    end do

  end function remove_string
  !============================================================================!
  function count_char(input, char) result(output)
    implicit none
    character(len=*) :: input
    character(len=1) :: char
    integer :: output
    integer :: i

    output = 0
    do i=1, len(input)
      if(input(i:i).eq.char) output = output + 1
    end do

  end function count_char
  !============================================================================!
  function split(this, char) result(array)
    class(string), intent(in out) :: this    
    character(len=1) :: char
    type(string), allocatable, dimension(:) :: array
    integer :: i, j, n

    n = count_char(this%string_data, char)

    if(n.eq.0) stop "Error: there is no such character."

    allocate(array(n+1))

    do i=1, n
      array(i) = ''
    end do

    i = 1
    do j=1, this%length()
      if(this%string_data(j:j).ne.char) then
        array(i) = array(i) + this%string_data(j:j)
      else
        i = i + 1
      end if
    end do

  end function split
  !============================================================================!

end module strings_lib

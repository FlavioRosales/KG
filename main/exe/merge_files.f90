program merge_files
    
implicit none

character(len=10) :: file_name, file_extension
integer :: np !Numero de procesadores

!real(kind=8), allocatable, dimension(:) :: axis1, axis2, fun

!character(len=50) :: name
!integer :: ios

!integer points_counter2D
!integer points_counter2D_subgroup
!integer num_points, num_points_sub

!integer i,j

np = 4
file_extension = '.asc'

!--------------------------1D merge------------------------------------
file_name = 'Bx'
call merge1D_group(file_name, file_extension,np)

file_name = 'By'
call merge1D_group(file_name, file_extension,np)

file_name = 'Bz'
call merge1D_group(file_name, file_extension,np)

file_name = 'vx'
call merge1D_group(file_name, file_extension,np)

file_name = 'vy'
call merge1D_group(file_name, file_extension,np)

file_name = 'vz'
call merge1D_group(file_name, file_extension,np)

file_name = 'rho'
call merge1D_group(file_name, file_extension,np)

file_name = 'W'
call merge1D_group(file_name, file_extension,np)

file_name = 'e'
call merge1D_group(file_name, file_extension,np)

file_name = 'p'
call merge1D_group(file_name, file_extension,np)
!--------------------------1D merge------------------------------------

!file_name = 'rho'
!file_plane = 'xy'

!num_points = points_counter2D(file_name,file_extension, file_plane)
!num_points_sub = points_counter2D_subgroup(file_name,file_extension,file_plane)

!allocate(axis1(1:num_points), axis2(1:num_points), fun(1:num_points))

!open(300,file='rho.xy.asc')

!do i=1, np
!    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', i, '.', trim(file_plane), trim(file_extension)
!    open(unit=i, file=name, iostat=ios, status="old", action="read")
!        if(ios/=0) stop 'Error opening file'
!end do

!do while(ios.eq.0)
!    do i=1, np
!        do j=1, num_points_sub
!            read(i,*,iostat=ios) axis1(j), axis2(j) ,fun(j)
!            if(ios/=0) exit
!            write(300,*) axis1(j), axis2(j) ,fun(j)
!        end do
!    end do
!    write(300,*)
!    write(300,*)
!end do

!do i=1, np
!    close(unit=i, iostat=ios)
!    if(ios/=0) stop 'Error closing file'
!end do

!close(300)

end program

integer function points_counter(file_name, file_extension, file_axis)

    character(len=10), intent(in) :: file_name, file_extension, file_axis
    character(len=50)  name
    character(len=50) line
  
    integer ios

    points_counter = -1

    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', 1, '.', trim(file_axis), trim(file_extension)

    open(unit=1, file=name, iostat=ios, status="old", action="read")
    if(ios/=0) stop 'Error opening file'
    read(1, '(A)', iostat=ios) line
    do while(len_trim(line)/=0)
        points_counter = points_counter+1
        read(1, '(A)', iostat=ios) line
    end do
    close(unit=1, iostat=ios)
    if(ios/=0) stop 'Error closing file'

end function points_counter

integer function points_counter2D_subgroup(file_name, file_extension, file_plane)

    character(len=10), intent(in) :: file_name, file_extension, file_plane
    character(len=50)  name
    character(len=50) line
  
    integer ios

    points_counter2D_subgroup = 0

    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', 1, '.', trim(file_plane), trim(file_extension)


    open(unit=1, file=name, iostat=ios, status="old", action="read")
    if(ios/=0) stop 'Error opening file'
    read(1, '(A)', iostat=ios) line
    do while(len_trim(line)/=0)
        points_counter2D_subgroup = points_counter2D_subgroup+1
        read(1, '(A)', iostat=ios) line
    end do
    close(unit=1, iostat=ios)
    if(ios/=0) stop 'Error closing file'

end function points_counter2D_subgroup

integer function points_counter2D(file_name, file_extension, file_plane)

    character(len=10), intent(in) :: file_name, file_extension, file_plane
    character(len=50)  name
    character(len=50) line
  
    integer ios

    integer num_spaces, tot_blank_lines

    points_counter2D = 0
    num_spaces = 0
    tot_blank_lines = 0

    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', 1, '.', trim(file_plane), trim(file_extension)
    open(unit=1, file=name, iostat=ios, status="old", action="read")
    if(ios/=0) stop 'Error opening file'
    read(1, '(A)', iostat=ios) line
    if(len_trim(line)==0) num_spaces = num_spaces + 1
    do while(num_spaces/=3)
        points_counter2D = points_counter2D+1
        read(1, '(A)', iostat=ios) line
        if(len_trim(line)==0) then
            num_spaces = num_spaces + 1
            tot_blank_lines = tot_blank_lines + 1
        end if
        if(len_trim(line)/=0) num_spaces = 0
    end do
    close(unit=1, iostat=ios)
    if(ios/=0) stop 'Error closing file'

    points_counter2D = points_counter2D - tot_blank_lines - 2

end function points_counter2D

integer function files_counter(file_name, file_extension, file_axis, np)

    implicit none

    integer, intent(in) :: np
    character(len=10), intent(in) :: file_name, file_extension, file_axis
    character(len=50)  name
    logical :: non_duplicate 
    real(kind=8), allocatable ,dimension(:) :: first_value_axis


    real(kind=8) :: fun
    integer ios
    integer i,k

    files_counter = 0

    allocate(first_value_axis(1:np))

    do i=1, np 
        write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', i, '.', trim(file_axis), trim(file_extension)
        open(unit=i, file=name, iostat=ios, status="old", action="read")
            if(ios/=0) stop 'Error opening file'
        non_duplicate = .false.
        read(i,*) first_value_axis(i), fun
        if(i.gt.1) then
         do k=i-1, 1, -1
            if((first_value_axis(k).eq.first_value_axis(i))) then
                non_duplicate = .false.
            end if
         end do
        end if
        
        if(non_duplicate) then
            files_counter = files_counter + 1
        end if

        close(unit=i, iostat=ios)
            if(ios/=0) stop 'Error closing file'
    end do

end function files_counter

subroutine index(file_name, file_extension, file_axis, np, index_files, num) 

    implicit none

    integer, intent(in) :: np,num
    character(len=10), intent(in) :: file_name, file_extension, file_axis
    integer, dimension(1:num), intent(inout) :: index_files
    character(len=50)  name
    logical :: non_duplicate 
    real(kind=8), allocatable ,dimension(:) :: first_value_axis

   real(kind=8) :: fun
   integer ios
    integer i,k,j

    j=0

    allocate(first_value_axis(1:np))
    do i=1, np 
        write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', i, '.', trim(file_axis), trim(file_extension)
        open(unit=i, file=name, iostat=ios, status="old", action="read")
            if(ios/=0) stop 'Error opening file'
        non_duplicate = .false.
        read(i,*) first_value_axis(i), fun
        
        if(i.gt.1) then
         do k=i-1, 1, -1
            if((first_value_axis(k).eq.first_value_axis(i))) then
                non_duplicate = .false.
            end if
         end do
        end if
        
        if(non_duplicate) then
            j = j + 1
            index_files(j)  = i 
        end if

        close(unit=i, iostat=ios)
            if(ios/=0) stop 'Error closing file'
    end do

end subroutine index

subroutine merge1D_file(file_name, file_extension, file_axis ,np)

    implicit none

    character(len=10), intent(in) :: file_name, file_extension, file_axis
    integer, intent(in) :: np
    integer  num, files_counter
    integer  num_points, points_counter
    real(kind=8), allocatable, dimension(:) :: axis, fun
    integer, allocatable ,dimension(:) :: index_files
    character(len=50)  name
    character(len=200) :: comando
    integer :: estado

    integer ios
    integer i,j

    num = files_counter(file_name, file_extension, file_axis,np)
    allocate(index_files(1:num))
    call index(file_name, file_extension, file_axis, np, index_files, num)
    num_points = points_counter(file_name,file_extension, file_axis)
    allocate(axis(0:num_points), fun(0:num_points))


    write(name,'(a,a,a,a)') trim(file_name),'.',trim(file_axis), trim(file_extension)

    open(200,file=name)

    do i=1, num
        write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', index_files(i), '.', trim(file_axis), trim(file_extension)
        open(unit=i, file=name, iostat=ios, status="old", action="read")
            if(ios/=0) stop 'Error opening file'
    end do

    do while(ios.eq.0)
        do i=1, num
            do j=0, num_points
                read(i,*,iostat=ios) axis(j), fun(j)
                if(ios/=0) exit
                write(200,*) axis(j), fun(j)
            end do
        end do
        write(200,*)
        write(200,*)
    end do

    do i=1, num
        close(unit=i, iostat=ios)
        if(ios/=0) stop 'Error closing file'
    end do

    do i=1, np
        write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', i, '.', trim(file_axis), trim(file_extension)
        write(comando, '(a,a)') 'rm ', trim(name)
        ! Utiliza la función system para ejecutar el comando
        estado = system(comando)

        ! Verifica si el archivo fue eliminado exitosamente
        if (estado == 0) then
        else
            print *, 'Error al intentar eliminar el archivo.'
        end if
    end do

    close(200)

end subroutine merge1D_file

subroutine merge1D_group(file_name,file_extension,np)

    implicit none
    character(len=10), intent(in) :: file_name, file_extension
    integer, intent(in) :: np
    character(len=50) :: name
    character(len=10) :: file_axis
    logical :: exist

    file_axis='x'
    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', 1, '.', trim(file_axis), trim(file_extension)
    inquire(file=name,exist=exist)
    if(exist) call merge1D_file(file_name,file_extension, file_axis, np)

    file_axis='y'
    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', 1, '.', trim(file_axis), trim(file_extension)
    inquire(file=name,exist=exist)
    if(exist) call merge1D_file(file_name,file_extension, file_axis, np)

    file_axis='z'
    write(name, '(a,a,i0,a,a,a)') trim(file_name), '_', 1, '.', trim(file_axis), trim(file_extension)
    inquire(file=name,exist=exist)
    if(exist) call merge1D_file(file_name,file_extension, file_axis, np)

end subroutine merge1D_group

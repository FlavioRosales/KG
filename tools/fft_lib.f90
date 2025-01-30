module fft_lib
    implicit none 

    real(kind=8), parameter :: pi = acos(-1.0d0)
    complex(kind=8), parameter :: i_unreal = (0.0d0,1.0d0)


    interface save_data
        module procedure :: save_data_1D, save_data_2D 
    end interface

    interface fft
        module procedure :: fft_1D, fft_2D, fft_3D
    end interface

    interface ifft
        module procedure :: ifft_2D, ifft_1D, ifft_3D
    end interface

    interface shift
        module procedure :: shift_1D, shift_2D, shift_3D
    end interface

    contains

    function linspace(min, max, Npoints) result(domain)
        implicit none 
        !
        ! input 
        !
        integer, intent(in) :: Npoints
        real(kind=8), intent(in) :: min, max
        !
        ! output 
        !
        real(kind=8), dimension(Npoints) :: domain 
        !
        ! others 
        !
        integer i 
        real(kind=8) :: h 

        h = (max - min)/dble(Npoints-1)

        do i=1, Npoints 
            domain(i) = min + h*dble(i-1)
        end do
    end function linspace



    subroutine save_data_1D(t,p,filename)
        implicit none 
        ! 
        ! input 
        !
        real(kind=8), dimension(:), intent(in) :: t,p
        character(len=*), intent(in) :: filename 
        ! 
        ! others 
        !
        integer i 
        logical :: exist 

        inquire(file = filename, exist = exist)

        if(exist) then 
        open(123, file = filename, status='old', action='write', position='append')

        else
            open(123, file = filename, status='new', action='write')
        end if  

        do i=1, size(t)
            write(123,*) t(i), p(i)
        end do
    end subroutine save_data_1D

    subroutine save_data_2D(x,y,f,filename)
        implicit none 
        ! 
        ! input 
        !
        real(kind=8), dimension(:,:), intent(in) :: x,y
        complex(kind=8), dimension(:,:), intent(in) :: f
        character(len=*), intent(in) :: filename 
        ! 
        ! others 
        !
        integer i,j        
        
        logical :: exist 

        inquire(file = filename, exist = exist)

        if(exist) then 
        open(123, file = filename, status='old', action='write', position='append')

        else
            open(123, file = filename, status='new', action='write')
        end if  

        do i=1, size(f,1)
            do j=1, size(f,2)
                write(123,*) x(i,j), y(i,j), dreal(f(i,j)), dimag(f(i,j)), abs(f(i,j))
            end do
            write(123,*)
        end do
        write(123,*)
        write(123,*)

    end subroutine save_data_2D

    function DFT(f) result(DFT_f)
        implicit none 
        !
        ! input
        !
        complex(kind=8), dimension(:), intent(in) :: f 
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: DFT_f
        !
        ! others
        !
        integer :: j,k, N 

        N = size(f)

        DFT_f = 0.0d0
        do k=0, N-1
            do j=0, N-1
                DFT_f(k+1) = DFT_f(k+1) +  &
                            f(j+1)*exp(2.0d0*pi*i_unreal*dble(j*k)/dble(N))
            end do
        end do

    end function DFT

    recursive function fft_1D(f) result(fft_f)
    implicit none 
    !
    ! input
    !
    complex(kind=8), dimension(:), intent(in) :: f 
    !
    ! output
    !
    complex(kind=8), dimension(size(f)) :: fft_f
    complex(kind=8), allocatable, dimension(:) ::  fft_f_even,fft_f_odd
    !
    ! others
    !
    integer :: k, N 

    N = size(f)



    if(mod(N,2).eq.0) then 
        allocate(fft_f_even(N/2),fft_f_odd(N/2))
        fft_f_even = fft_1d(f = f(1::2))
        fft_f_odd  = fft_1d(f = f(2::2))

        do k = 0, N/2 -1
            fft_f(k+1) = fft_f_even(k+1) + &
            exp(2.0d0*pi*i_unreal*dble(k)/dble(N))*fft_f_odd(k+1)

            fft_f(k+1+N/2) = fft_f_even(k+1) - &
            exp(2.0d0*pi*i_unreal*dble(k)/dble(N))*fft_f_odd(k+1)
        end do
    else
        fft_f = DFT(f)
    end if


    end function fft_1D

    function fft_2D(f) result(fft_f)
        implicit none 
        !
        ! input 
        !
        complex(kind=8), dimension(:,:), intent(in) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f,1),size(f,2)) :: fft_f
        integer ::  i,j

        fft_f = f

        do j=1, size(f,2)
            fft_f(:,j) = fft_1D(f=fft_f(:,j))
        end do

        do i=1, size(f,1)
            fft_f(i,:) = fft_1D(f=fft_f(i,:))
        end do
    end function fft_2D

    function fft_3D(f) result(fft_f)
        implicit none 
        !
        ! input
        !
        complex(kind=8), intent(in), dimension(:,:,:) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f,1),size(f,2),size(f,3)) :: fft_f
        complex(kind=8), dimension(size(f,1),size(f,2),size(f,3)) :: data_2D
        !
        ! others
        !
        integer :: i,j, k

        fft_f = f 

        do k=1, size(f,3)
            data_2D(:,:,k) = fft_2D(f=f(:,:,k))
        end do
        do i=1, size(f,1)
            do j=1, size(f,2)
                fft_f(i,j,:) = fft_1D(f=data_2D(i,j,:))
            end do
        end do
        
    end function fft_3D

    function ifft_2D(f) result(fft_f)
        implicit none 
        !
        ! input 
        !
        complex(kind=8), dimension(:,:), intent(in) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f,1),size(f,2)) :: fft_f

        fft_f = dconjg(fft_2D(dconjg(f)))/dble(size(f,1)*size(f,2))
    end function ifft_2D

    function ifft_3D(f) result(fft_f)
        implicit none 
        !
        ! input 
        !
        complex(kind=8), dimension(:,:,:), intent(in) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f,1),size(f,2),size(f,3)) :: fft_f

        fft_f = dconjg(fft_3D(dconjg(f)))/dble(size(f,1)*size(f,2)*size(f,3))
    end function ifft_3D

    function ifft_1D(f) result(fft_f)
        implicit none 
        !
        ! input 
        !
        complex(kind=8), dimension(:), intent(in) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: fft_f

        fft_f = dconjg(fft_1D(dconjg(f)))/dble(size(f))
    end function ifft_1D

    function iDFT(f) result(iDFT_f)
        implicit none 
        !
        ! input
        !
        complex(kind=8), dimension(:), intent(in) :: f 
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: iDFT_f
        !
        ! others
        !
        integer :: N 

        N = size(f)

        iDFT_f = dconjg(DFT(dconjg(f)))/dble(N)

    end function iDFT

    function shift_1D(f) result(shift_f)
        !
        ! input
        !
        complex(kind=8), dimension(:), intent(in) :: f
        !
        ! output
        !
        complex(kind=8), dimension(size(f)) :: shift_f
        !
        ! others
        !

        integer :: N; N = size(f)/2

        shift_f(N+1:2*N) = f(1:N)
        shift_f(1:N) = f(N+1:2*N)
    end function shift_1D

    function shift_2D(f) result(shift_f)
        implicit none
        complex(kind=8), dimension(:,:), intent(in) :: f
        complex(kind=8), dimension(size(f,1),size(f,2)) :: shift_f
    
        integer :: Nx,Ny
        Nx = size(f,1)/2
        Ny = size(f,2)/2
    
        shift_f(1:Nx,1:Ny) = f(Nx+1:2*Nx,Ny+1:2*Ny)
        shift_f(1:Nx,Ny+1:2*Ny) = f(Nx+1:2*Nx,1:Ny)
        shift_f(Nx+1:2*Nx,1:Ny) = f(1:Nx,Ny+1:2*Ny)
        shift_f(Nx+1:2*Nx,Ny+1:2*Ny) = f(1:Nx,1:Ny)
    end function shift_2D

    function shift_3D(f) result(shift_f)
        implicit none
        complex(kind=8), dimension(:,:,:), intent(in) :: f
        complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: shift_f
    
        integer :: Nx, Ny, Nz
    
        Nx = size(f,1)/2
        Ny = size(f,2)/2
        Nz = size(f,3)/2
    
        shift_f(1:Nx, 1:Ny, 1:Nz) = f(Nx+1:2*Nx, Ny+1:2*Ny, Nz+1:2*Nz)
        shift_f(1:Nx, 1:Ny, Nz+1:2*Nz) = f(Nx+1:2*Nx, Ny+1:2*Ny, 1:Nz)
        shift_f(1:Nx, Ny+1:2*Ny, 1:Nz) = f(Nx+1:2*Nx, 1:Ny, Nz+1:2*Nz)
        shift_f(1:Nx, Ny+1:2*Ny, Nz+1:2*Nz) = f(Nx+1:2*Nx, 1:Ny, 1:Nz)
        shift_f(Nx+1:2*Nx, 1:Ny, 1:Nz) = f(1:Nx, Ny+1:2*Ny, Nz+1:2*Nz)
        shift_f(Nx+1:2*Nx, 1:Ny, Nz+1:2*Nz) = f(1:Nx, Ny+1:2*Ny, 1:Nz)
        shift_f(Nx+1:2*Nx, Ny+1:2*Ny, 1:Nz) = f(1:Nx, 1:Ny, Nz+1:2*Nz)
        shift_f(Nx+1:2*Nx, Ny+1:2*Ny, Nz+1:2*Nz) = f(1:Nx, 1:Ny, 1:Nz)
    
    end function shift_3D

end module fft_lib
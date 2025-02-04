module geometry_lib
  use mesh_refinement_lib
  implicit none

  type, extends(mesh_refinement)  :: geometry

    character(len=20) :: metric
    real(kind=8) :: a, M_bh
    real(kind=8), allocatable, dimension(:,:,:) :: alpha, sqrtg
    real(kind=8), allocatable, dimension(:,:,:,:) :: beta
    real(kind=8), allocatable, dimension(:,:,:,:,:) :: g_dn, g_up, gamma_up
    real(kind=8), allocatable, dimension(:,:,:,:,:,:) :: Gamma

  contains

    procedure :: D
    procedure :: geometry_memory
    procedure :: geometry_load_metric
    procedure :: geometry_load_inverse
    procedure :: validate_inverse_metric
    procedure :: geometry_load_Christoffel_symbols
    procedure :: boundary_theta_dn, boundary_theta_up

  end type geometry

contains

  subroutine geometry_memory(this)
    implicit none
    class(geometry), intent(in out) :: this

    allocate(this%alpha( &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%sqrtg(&
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%beta(3, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%g_dn(0:3, 0:3, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%gamma_up(3, 3, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%g_up(0:3, 0:3, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

    allocate(this%Gamma(3, 3, 3, &
    this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR))

  end subroutine geometry_memory

  subroutine geometry_load_metric(this)
    implicit none
    class(geometry), intent(in out) :: this
    real(kind=8), dimension(this%iL - this%g_pts:this%iR + this%g_pts, &
    this%jL - this%g_pts:this%jR + this%g_pts , &
    this%kL:this%kR) :: r_ks, theta_ks, d_theta
    real(kind=8) :: h_mks


    if(this%metric.eq.'Minkowski') then

      this%g_dn = 0.0d0
      this%g_dn(0,0,:,:,:) = -1.0d0
      this%g_dn(1,1,:,:,:) = 1.0d0
      this%g_dn(2,2,:,:,:) = this%x**2
      this%g_dn(3,3,:,:,:) = this%x**2 * dsin(this%y)**2

    else if(this%metric.eq.'Kerr_Schild') then

      if(rank.eq.0) then 
        print*, 'metric = ', this%metric, 'a = ',this%a, 'M = ', this%M_bh
      end if


      this%g_dn = 0.0d0

      this%g_dn(0,0,:,:,:) = -1.0d0 + 2.0d0*this%M_bh*this%x/(this%x**2 &
      +this%a**2*dcos(this%y)**2)

      this%g_dn(0,1,:,:,:) = (2.0d0*this%M_bh*this%x/(this%x**2 &
                  + this%a**2 * dcos(this%y)**2))

      this%g_dn(0,2,:,:,:) = 0.0d0

      this%g_dn(0,3,:,:,:) = - 2.0d0*this%a*this%M_bh*this%x*dsin(this%y)**2 &
                  /(this%x**2 + this%a**2 * dcos(this%y)**2)   

      this%g_dn(1,0,:,:,:) = this%g_dn(0,2,:,:,:)

      this%g_dn(1,1,:,:,:) = ( 1.0d0 + 2.0d0*this%M_bh*this%x &
                  / (this%x**2 + this%a**2 * dcos(this%y)**2))
              
      this%g_dn(1,2,:,:,:) = 0.0d0

      this%g_dn(1,3,:,:,:) = (- this%a * (1.0d0 + 2.0d0*this%M_bh*this%x &
                  /(this%x**2 + this%a**2 * dcos(this%y)**2)) &
                  * dsin(this%y)**2)

      this%g_dn(2,0,:,:,:) = 0.0d0 

      this%g_dn(2,1,:,:,:) = 0.0d0 

      this%g_dn(2,2,:,:,:) =( this%x**2 + this%a**2 * dcos(this%y)**2 ) 

      this%g_dn(2,3,:,:,:) = 0.0d0 

      this%g_dn(3,0,:,:,:) = this%g_dn(0,3,:,:,:)

      this%g_dn(3,1,:,:,:) = this%g_dn(1,3,:,:,:)

      this%g_dn(3,2,:,:,:) = 0.0d0

      this%g_dn(3,3,:,:,:) = ((this%a**2 + this%x**2)**2 &
                      - this%a**2 * (-2.0d0*this%M_bh*this%x+this%a**2 + this%x**2) &
                      *dsin(this%y)**2)*dsin(this%y)**2 &
                      /(this%x**2 + this%a**2 * dcos(this%y)**2)
          
    else if(this%metric.eq.'MKS') then
      h_mks = 0.0d0

      if(rank.eq.0) then 
        print*, 'metric = ', this%metric, 'a = ',this%a, 'M = ', this%M_bh
      end if

      r_ks = dexp(this%x)
      theta_ks = this%y + 0.5d0*h_mks*dsin(2.0d0*this%y)
      d_theta = 1.0d0 + h_mks*dcos(2.0d0*this%y)

      this%g_dn(0,0,:,:,:) = -1.0d0 + 2.0d0*this%M_bh*r_ks/(r_ks**2 &
            +this%a**2*dcos(theta_ks)**2)

      this%g_dn(0,1,:,:,:) = (2.0d0*this%M_bh*r_ks/(r_ks**2 &
                  + this%a**2 * dcos(theta_ks)**2))*r_ks

      this%g_dn(0,2,:,:,:) = 0.0d0

      this%g_dn(0,3,:,:,:) = - 2.0d0*this%a*this%M_bh*r_ks*dsin(theta_ks)**2 &
                  /(r_ks**2 + this%a**2 * dcos(theta_ks)**2)   

      this%g_dn(1,0,:,:,:) = this%g_dn(0,1,:,:,:)

      this%g_dn(1,1,:,:,:) = ( 1.0d0 + 2.0d0*this%M_bh*r_ks &
                  / (r_ks**2 + this%a**2 * dcos(theta_ks)**2))*r_ks**2
              
      this%g_dn(1,2,:,:,:) = 0.0d0
    
      this%g_dn(1,3,:,:,:) = (- this%a * (1.0d0 + 2.0d0*this%M_bh*r_ks &
                  /(r_ks**2 + this%a**2 * dcos(theta_ks)**2)) &
                  * dsin(theta_ks)**2)*r_ks

      this%g_dn(2,0,:,:,:) = 0.0d0 

      this%g_dn(2,1,:,:,:) = 0.0d0 

      this%g_dn(2,2,:,:,:) =( r_ks**2 + this%a**2 * dcos(theta_ks)**2 ) * d_theta**2

      this%g_dn(2,3,:,:,:) = 0.0d0 

      this%g_dn(3,0,:,:,:) = this%g_dn(0,3,:,:,:)

      this%g_dn(3,1,:,:,:) = this%g_dn(1,3,:,:,:)

      this%g_dn(3,2,:,:,:) = 0.0d0

      this%g_dn(3,3,:,:,:) = ((this%a**2 + r_ks**2)**2 &
                      - this%a**2 * (-2.0d0*this%M_bh*r_ks+this%a**2 + r_ks**2) &
                      *dsin(theta_ks)**2)*dsin(theta_ks)**2 &
                      /(r_ks**2 + this%a**2 * dcos(theta_ks)**2)
    end if

    call this%boundary_theta_dn()

    ! Spatial metric determinant
    this%sqrtg = dsqrt(this%g_dn(2,2,:,:,:)*(this%g_dn(1,1,:,:,:)*this%g_dn(3,3,:,:,:) - this%g_dn(1,3,:,:,:)**2))

  end subroutine geometry_load_metric

  subroutine geometry_load_inverse(this)
    implicit none
    class(geometry), intent(in out) :: this 
    real(kind=8), dimension(&
    this%iL-this%g_pts:this%iR+this%g_pts, &
    this%jL-this%g_pts:this%jR+this%g_pts, &
    this%kL:this%kR) :: idetg 
    integer :: i,j
  
    idetg = 1.0d0/(((this%g_dn(3,3,:,:,:)*this%g_dn(1,1,:,:,:) &
            - this%g_dn(1,3,:,:,:)**2)*this%g_dn(0,0,:,:,:) &
            - this%g_dn(3,3,:,:,:)*this%g_dn(0,1,:,:,:)**2 &
            + 2.0d0*this%g_dn(1,3,:,:,:)*this%g_dn(0,3,:,:,:)*this%g_dn(0,1,:,:,:)&
            - this%g_dn(1,1,:,:,:)*this%g_dn(0,3,:,:,:)**2)*this%g_dn(2,2,:,:,:))
  
    this%g_up(0,0,:,:,:) = this%g_dn(2,2,:,:,:)*(this%g_dn(1,1,:,:,:)*this%g_dn(3,3,:,:,:) &
                          - this%g_dn(1,3,:,:,:)**2)*idetg
  
    this%g_up(0,1,:,:,:) = this%g_dn(2,2,:,:,:)*(this%g_dn(1,3,:,:,:)*this%g_dn(0,3,:,:,:) &
                          - this%g_dn(0,1,:,:,:)*this%g_dn(3,3,:,:,:))*idetg
  
    this%g_up(0,2,:,:,:) = 0.0d0 
  
    this%g_up(0,3,:,:,:) = this%g_dn(2,2,:,:,:)*(this%g_dn(1,3,:,:,:)*this%g_dn(0,1,:,:,:) &
                          - this%g_dn(1,1,:,:,:)*this%g_dn(0,3,:,:,:))*idetg
  
    this%g_up(1,0,:,:,:) = this%g_up(0,1,:,:,:)
  
    this%g_up(1,1,:,:,:) = this%g_dn(2,2,:,:,:)*(this%g_dn(3,3,:,:,:)*this%g_dn(0,0,:,:,:) &
                          - this%g_dn(0,3,:,:,:)**2)*idetg
  
    this%g_up(1,2,:,:,:) = 0.0d0
  
    this%g_up(1,3,:,:,:) = this%g_dn(2,2,:,:,:)*(-this%g_dn(1,3,:,:,:)*this%g_dn(0,0,:,:,:) &
                          + this%g_dn(0,3,:,:,:)*this%g_dn(0,1,:,:,:))*idetg
  
    this%g_up(2,0,:,:,:) = this%g_up(0,2,:,:,:)
  
    this%g_up(2,1,:,:,:) = this%g_up(1,2,:,:,:)
  
    this%g_up(2,2,:,:,:) = ((this%g_dn(3,3,:,:,:)*this%g_dn(1,1,:,:,:) - this%g_dn(1,3,:,:,:)**2)*this%g_dn(0,0,:,:,:) &
                          - this%g_dn(3,3,:,:,:)*this%g_dn(0,1,:,:,:)**2 &
                          + 2.0d0*this%g_dn(1,3,:,:,:)*this%g_dn(0,3,:,:,:)*this%g_dn(0,1,:,:,:) &
                          - this%g_dn(1,1,:,:,:)*this%g_dn(0,3,:,:,:)**2)*idetg
  
    this%g_up(2,3,:,:,:) = 0.0d0
    
    this%g_up(3,0,:,:,:) = this%g_up(0,3,:,:,:)
  
    this%g_up(3,1,:,:,:) = this%g_up(1,3,:,:,:)
  
    this%g_up(3,2,:,:,:) = this%g_up(2,3,:,:,:)
  
    this%g_up(3,3,:,:,:) = this%g_dn(2,2,:,:,:)*(this%g_dn(1,1,:,:,:)*this%g_dn(0,0,:,:,:) &
                          - this%g_dn(0,1,:,:,:)**2)*idetg

    call this%boundary_theta_up()

    !metric inverse
    this%alpha = dsqrt(-1.0d0/this%g_up(0,0,:,:,:))

    this%beta(1,:,:,:) = this%alpha**2 * this%g_up(0,1,:,:,:)
    this%beta(2,:,:,:) = this%alpha**2 * this%g_up(0,2,:,:,:)
    this%beta(3,:,:,:) = this%alpha**2 * this%g_up(0,3,:,:,:)

    do i=1, 3
      do j=1, 3
        this%gamma_up(i,j,:,:,:) = this%g_up(i,j,:,:,:) + this%beta(i,:,:,:)*this%beta(j,:,:,:)/this%alpha**2
      end do
    end do


  end subroutine geometry_load_inverse

  subroutine geometry_load_Christoffel_symbols(this)
    implicit none
    class(geometry), intent(in out) :: this

    integer :: i, j, k, s

    this%Gamma = 0.0d0

    do i=1, 3
      do j=1, 3
        do k=1, 3
          do s=1, 3

            this%Gamma(i,j,k,:,:,:) = this%Gamma(i,j,k,:,:,:) + &
            0.50d0 * this%g_up(k,s,:,:,:) * (&
            this%D(this%g_dn(j,s,:,:,:), i) + this%D(this%g_dn(s,i,:,:,:), j) - this%D(this%g_dn(i,j,:,:,:), s))

          end do
        end do
      end do
    end do

  end subroutine geometry_load_Christoffel_symbols

  subroutine boundary_theta_dn(this)
    implicit none 
    class(geometry), intent(in out) :: this 
    integer :: k,dis, g

    dis = (this%Nz+ 1)/2

    do k = this%kL, this%kR
      do g=this%g_pts, 0, -1
        if(k.lt.dis) then 
          this%g_dn(1,1,:,this%jL-g,k) = this%g_dn(1,1,:,this%jL+g+1,k+dis)
          this%g_dn(1,2,:,this%jL-g,k) = -this%g_dn(1,2,:,this%jL+g+1,k+dis)
          this%g_dn(1,3,:,this%jL-g,k) = -this%g_dn(1,3,:,this%jL+g+1,k+dis)
          this%g_dn(2,2,:,this%jL-g,k) = this%g_dn(2,2,:,this%jL+g+1,k+dis)
          this%g_dn(2,3,:,this%jL-g,k) = this%g_dn(2,3,:,this%jL+g+1,k+dis)
          this%g_dn(3,3,:,this%jL-g,k) = this%g_dn(3,3,:,this%jL+g+1,k+dis)
        else 
          this%g_dn(1,1,:,this%jL-g,k) = this%g_dn(1,1,:,this%jL+g+1,k-dis)
          this%g_dn(1,2,:,this%jL-g,k) = -this%g_dn(1,2,:,this%jL+g+1,k-dis)
          this%g_dn(1,3,:,this%jL-g,k) = -this%g_dn(1,3,:,this%jL+g+1,k-dis)
          this%g_dn(2,2,:,this%jL-g,k) = this%g_dn(2,2,:,this%jL+g+1,k-dis)
          this%g_dn(2,3,:,this%jL-g,k) = this%g_dn(2,3,:,this%jL+g+1,k-dis)
          this%g_dn(3,3,:,this%jL-g,k) = this%g_dn(3,3,:,this%jL+g+1,k-dis)
        end if
      end do
    end do

    do k = this%kL, this%kR
      do g= this%g_pts, 0, -1
        if(k.lt.dis) then 
          this%g_dn(1,1,:,this%jR + g,k) = this%g_dn(1,1,:,this%jR-g-1,k+dis)
          this%g_dn(1,2,:,this%jR + g,k) = -this%g_dn(1,2,:,this%jR-g-1,k+dis)
          this%g_dn(1,3,:,this%jR + g,k) = -this%g_dn(1,3,:,this%jR-g-1,k+dis)
          this%g_dn(2,2,:,this%jR + g,k) = this%g_dn(2,2,:,this%jR-g-1,k+dis)
          this%g_dn(2,3,:,this%jR + g,k) = this%g_dn(2,3,:,this%jR-g-1,k+dis)
          this%g_dn(3,3,:,this%jR + g,k) = this%g_dn(3,3,:,this%jR-g-1,k+dis)
        else 
          this%g_dn(1,1,:,this%jR + g,k) = this%g_dn(1,1,:,this%jR-g-1,k-dis)
          this%g_dn(1,2,:,this%jR + g,k) = -this%g_dn(1,2,:,this%jR-g-1,k-dis)
          this%g_dn(1,3,:,this%jR + g,k) = -this%g_dn(1,3,:,this%jR-g-1,k-dis)
          this%g_dn(2,2,:,this%jR + g,k) = this%g_dn(2,2,:,this%jR-g-1,k-dis)
          this%g_dn(2,3,:,this%jR + g,k) = this%g_dn(2,3,:,this%jR-g-1,k-dis)
          this%g_dn(3,3,:,this%jR + g,k) = this%g_dn(3,3,:,this%jR-g-1,k-dis)
        end if
      end do
    end do

  end subroutine boundary_theta_dn

  subroutine boundary_theta_up(this)
    implicit none 
    class(geometry), intent(in out) :: this 
    integer :: k,dis, g

    dis = (this%Nz+ 1)/2

    do k = this%kL, this%kR
      do g=this%g_pts, 0, -1
        if(k.lt.dis) then 
          this%g_up(1,1,:,this%jL-g,k) = this%g_up(1,1,:,this%jL+g+1,k+dis)
          this%g_up(1,2,:,this%jL-g,k) = -this%g_up(1,2,:,this%jL+g+1,k+dis)
          this%g_up(1,3,:,this%jL-g,k) = -this%g_up(1,3,:,this%jL+g+1,k+dis)
          this%g_up(2,2,:,this%jL-g,k) = this%g_up(2,2,:,this%jL+g+1,k+dis)
          this%g_up(2,3,:,this%jL-g,k) = this%g_up(2,3,:,this%jL+g+1,k+dis)
          this%g_up(3,3,:,this%jL-g,k) = this%g_up(3,3,:,this%jL+g+1,k+dis)
        else  
          this%g_up(1,1,:,this%jL-g,k) = this%g_up(1,1,:,this%jL+g+1,k-dis)
          this%g_up(1,2,:,this%jL-g,k) = -this%g_up(1,2,:,this%jL+g+1,k-dis)
          this%g_up(1,3,:,this%jL-g,k) = -this%g_up(1,3,:,this%jL+g+1,k-dis)
          this%g_up(2,2,:,this%jL-g,k) = this%g_up(2,2,:,this%jL+g+1,k-dis)
          this%g_up(2,3,:,this%jL-g,k) = this%g_up(2,3,:,this%jL+g+1,k-dis)
          this%g_up(3,3,:,this%jL-g,k) = this%g_up(3,3,:,this%jL+g+1,k-dis)
        end if
      end do
    end do

    do k = this%kL, this%kR
      do g= this%g_pts, 0, -1
        if(k.lt.dis) then 
          this%g_up(1,1,:,this%jR + g,k) = this%g_up(1,1,:,this%jR-g-1,k+dis)
          this%g_up(1,2,:,this%jR + g,k) = -this%g_up(1,2,:,this%jR-g-1,k+dis)
          this%g_up(1,3,:,this%jR + g,k) = -this%g_up(1,3,:,this%jR-g-1,k+dis)
          this%g_up(2,2,:,this%jR + g,k) = this%g_up(2,2,:,this%jR-g-1,k+dis)
          this%g_up(2,3,:,this%jR + g,k) = this%g_up(2,3,:,this%jR-g-1,k+dis)
          this%g_up(3,3,:,this%jR + g,k) = this%g_up(3,3,:,this%jR-g-1,k+dis)
        else 
          this%g_up(1,1,:,this%jR + g,k) = this%g_up(1,1,:,this%jR-g-1,k-dis)
          this%g_up(1,2,:,this%jR + g,k) = -this%g_up(1,2,:,this%jR-g-1,k-dis)
          this%g_up(1,3,:,this%jR + g,k) = -this%g_up(1,3,:,this%jR-g-1,k-dis)
          this%g_up(2,2,:,this%jR + g,k) = this%g_up(2,2,:,this%jR-g-1,k-dis)
          this%g_up(2,3,:,this%jR + g,k) = this%g_up(2,3,:,this%jR-g-1,k-dis)
          this%g_up(3,3,:,this%jR + g,k) = this%g_up(3,3,:,this%jR-g-1,k-dis)
        end if
      end do
    end do

  end subroutine boundary_theta_up

  function D(this, f, i) result(df)
    use finite_differences, only: &
    dx => first_derivative_x_2, &
    dy => first_derivative_y_2, &
    dz => first_derivative_z_2
    implicit none
    class(geometry), intent(in out) :: this
    integer, intent(in) :: i
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1),size(f,2),size(f,3)) :: df

    if(i.eq.1) df = dx(f, this%dx)
    if(i.eq.2) df = dy(f, this%dy)
    if(i.eq.3) df = dz(f, this%dz)

  end function D

  !============================================================================!

subroutine validate_inverse_metric(this)
  implicit none
  class(geometry), intent(in) :: this
  integer :: mu, lambda, nu, ix, iy, iz, errors
  real(kind=8) :: delta, tolerance

  errors = 0
  tolerance = 1.0d-10
  ix = this%iL; iy = this%jL; iz = this%kL  


  do ix = this%iL - this%g_pts, this%iR + this%g_pts
    do iy = this%jL - this%g_pts, this%jR + this%g_pts
      do iz = this%kL, this%kR 
  do mu = 0, 3
    do lambda = 0, 3
      delta = 0.0d0
      do nu = 0, 3
        delta = delta + this%g_dn(mu, nu, ix, iy, iz) * this%g_up(nu, lambda, ix, iy, iz)
      end do
      if (mu == lambda) then
        if (abs(delta - 1.0d0) > tolerance) then
          print *, "Error en Delta(", mu, ",", lambda, ") = ", delta
          errors = errors + 1
        end if
      else
        if (abs(delta) > tolerance) then
          print *, "Error en Delta(", mu, ",", lambda, ") = ", delta
          errors = errors + 1
        end if
      end if
    end do
  end do

      end do
    end do
  end do

   if (errors == 0) then
   print *, "La métrica inversa pasó la validación en el procesador ", rank
  else
   print *, "Se encontraron ", errors, " errores en la métrica inversa."
  end if

end subroutine validate_inverse_metric


end module geometry_lib

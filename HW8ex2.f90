    module polar_coords
   
    private
    public :: polar

    type polar

        real :: z, theta ! THETA is angle in degrees     

    contains

    procedure,private, pass(x) :: turn_complex, turn_polar
    procedure, private :: product_polar, division_polar
    generic,public :: assignment(=) => turn_complex, turn_polar
    generic, public :: operator(*) =>product_polar
    generic, public :: operator(/) => division_polar
    end type polar
    
    contains
    
    
    
    
    type(polar) function product_polar(p1, p2)
    class (polar), intent(in) :: p1, p2
    product_polar%z=p1%z * p2%z
    product_polar%theta=p1%theta + p2%theta
    do while (product_polar%theta > 180)
    product_polar%theta=product_polar%theta-360
    end do
    do while (product_polar%theta < -180)
    product_polar%theta=product_polar%theta+360
    end do
    end function product_polar
    
    
    type(polar) function division_polar(p1, p2)
    class (polar), intent(in) :: p1, p2
    division_polar%z=p1%z / p2%z
    division_polar%theta=p1%theta - p2%theta
    do while (division_polar%theta > 180)
    division_polar%theta=division_polar%theta-360
    end do
    do while (division_polar%theta < -180)
    division_polar%theta=division_polar%theta+360
    end do
    end function division_polar
    
    
    
    
    subroutine turn_complex(y, x)
    class(polar),intent(in):: x
    complex, intent(out) :: y
    real*8 :: pi, a, b
    pi=DACOS(-1.D0) ! defining pi 
    a=x%z * cos(x%theta*pi/180) ! REAL part
    b=x%z * sin(x%theta*pi/180) ! IMAGINARY part
    y=cmplx(a, b)
    end subroutine turn_complex
    
    subroutine turn_polar(x, number)
    complex,intent(in):: number
    class(polar), intent (out) :: x
    real*8 :: pi
    pi=DACOS(-1.D0) ! defining pi 
    x%theta = atan2(aimag(number), real(number)) * 180/pi 
    ! MULTIPLY by 180/pi so that theta is in degrees
    x%z = sqrt(aimag(number)**2 +real(number)**2)    
    end subroutine turn_polar
    
    end module polar_coords
    
    program HW8ex2
    use polar_coords

    implicit none
    
    type(polar) :: p1, p2, p3
    complex :: c1=cmplx(-1.0, 2.0), c2=cmplx(2.0, -0.5)
    
    ! test type conversion
    p1=c1
    c2=p1
    print*, p1
    print*, c2/c1
    
    ! test multiply and divide
    p2=c2
    p3=p1*p2
    print*, p3
    p3=p2 / p1
    print*, p3


    end program HW8ex2


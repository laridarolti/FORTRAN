    module coords
    implicit none

    private
    public :: point, pointT

    type point

        real :: x,y,z      

    contains

    procedure,private :: pointplus, pointminus, pointseparation, dot_prod, cross_prod
    procedure,private, pass(this) :: multiplication1, multiplication2 !, Tplus
    procedure, private :: division1, division2, multiplication3, multiplication4
    procedure,private, pass(this) :: absvec_real32, absvec_real64
    generic,public :: operator(+) => pointplus
    generic,public :: operator(-) => pointminus
    generic,public :: operator(.distance.) => pointseparation
    generic,public :: assignment(=) => absvec_real32, absvec_real64
    generic, public :: operator(.dot.) => dot_prod
    generic, public :: operator(.cross.) => cross_prod
    generic, public :: operator(*) => multiplication3, multiplication4, multiplication1, multiplication2
    generic, public :: operator(/) => division1, division2


    end type point

    ! EXTENDING TYPE POINT AND ADDING TEMPERATURE
    type, extends(point) :: pointT
        real :: T
    contains
    procedure, public :: Tplus, Tminus
    end type pointT
        
    contains

    type(point) function pointplus(this,b)   ! for +
    class(point),intent(in):: this,b
    pointplus%x = this%x + b%x
    pointplus%y = this%y + b%y
    pointplus%z = this%z + b%z
    end function pointplus
    
    real function Tplus(this, b)  !  for PLUS
    class(pointT),intent(in):: b, this
    TPlus = b%T + this%T
    end function Tplus

    real function Tminus(this, b)  !  for MINUS
    class(pointT),intent(in):: b, this
    Tminus = b%T - this%T
    end function Tminus

    type(point) function pointminus(this,b)  !  for -
    class(point),intent(in):: this,b
    pointminus%x = this%x - b%x
    pointminus%y = this%y - b%y
    pointminus%z = this%z - b%z
    end function pointminus

    type(point) function cross_prod(this,b)  !  for .cross.
    class(point),intent(in):: this,b
    cross_prod%x = this%y * b%z - this%z * b%y
    cross_prod%y = this%z * b%x - this%x * b%z
    cross_prod%z = this%x * b%y - this%y * b%x
    end function cross_prod

    type(point) function division1(this,b)  !  for /32
    class(point),intent(in):: this
    real*4, intent(in) :: b
    division1%x = this%x/b
    division1%y = this%y/b
    division1%z = this%z/b
    end function division1

    type(point) function division2(this,b)  !  for /64
    class(point),intent(in):: this
    real(kind=8), intent(in) :: b
    division2%x = this%x/b
    division2%y = this%y/b
    division2%z = this%z/b
    end function division2

    type(point) function multiplication1(b, this)  !  for *, b is 32, b*p
    real(kind=4), intent(in) :: b
    class(point),intent(in):: this
    multiplication1%x = this%x * b
    multiplication1%y = this%y * b
    multiplication1%z = this%z * b
    end function multiplication1
    
    type(point) function multiplication2(b, this)  !  for *, b is 64, b*p
    real(kind=8), intent(in) :: b
    class(point),intent(in):: this
    multiplication2%x = this%x * b
    multiplication2%y = this%y * b
    multiplication2%z = this%z * b
    end function multiplication2

    type(point) function multiplication3(this,b)  !  for *, b is 32, p*b
    real(kind=4), intent(in) :: b
    class(point),intent(in):: this
    multiplication3%x = this%x * b
    multiplication3%y = this%y * b
    multiplication3%z = this%z * b
    end function multiplication3

    type(point) function multiplication4(this,b)  !  for *, b is 64, p*b
    real*8, intent(in) :: b
    class(point),intent(in):: this
    multiplication4%x = this%x * b
    multiplication4%y = this%y * b
    multiplication4%z = this%z * b
    end function  multiplication4


    real function pointseparation(this,b) ! for .distance.
    class(point),intent(in):: this,b
    pointseparation = sqrt(  &
        (this%x-b%x)**2+(this%y-b%y)**2+(this%z-b%z)**2)
    end function pointseparation

    real function dot_prod(this,b) ! for .dot.
    class(point),intent(in):: this,b
    dot_prod = (this%x)*b%x+(this%y)*(b%y)+(this%z)*(b%z)
    end function dot_prod

    subroutine absvec_real32(a, this)  ! for = (distance
    real(kind=4), intent(out):: a     !  from origin)
    class(point),intent(in):: this
    a = sqrt(this%x**2+this%y**2+this%z**2)
    end subroutine absvec_real32

    subroutine absvec_real64(a, this)  ! for = (distance
    real(kind=8),intent(out):: a     !  from origin)
    class(point),intent(in):: this
    a = sqrt(this%x**2+this%y**2+this%z**2)
    end subroutine absvec_real64

    end module coords

    program HW8
    use iso_fortran_env
    use coords
    implicit none
    

    type(point) :: p1, p2, p3
    type(pointT) :: pt1, pt2
    real(real32) :: r4
    real(real64) :: r8
    p1%x=1.2;  p1%y=0.;  p1%z=3.1
    p2%x=0.;  p2%y=1.2;  p2%z=1.7
    !p3=p1+p2; print*, p3
    !p3=p1-p2; print*, p3
    !
    !pt1%T = 3; pt1%x=1.2;  pt1%y=0.;  pt1%z=3.1
    !pt2%T = 7; pt2%x=1.2;  pt2%y=0.;  pt2%z=3.1
    !
    !d=p1; print*, d
    !d2=p1.distance.p2; print*, d2
    !
    !do1 = p1; print*, do1
    !do2 = p1; print*, do2
    !
    !r1 = p1 .cross. p2; print*, r1
    !r3 = p1 .dot. p2; print*, r3
    !
    !r2 = p1*c4; print*, r2
    !r2 = c4*p1; print*, r2
    !
    !r2 = p1*c8; print*, r2
    !r2 = c8*p1; print*, r2
    !
    !print *, pt1%Tplus(pt2)
    
    r4=p1; print*, r4
    r8=p1; print*, r8
    r4=p1.dot.p2; print*, r4
    p3=p1.cross.p2; print*, p3
    p3=p1*r4; print*, p3
    p3=p1*r8; print*, p3
    p3=r4*p1; print*, p3
    p3=r8*p1; print*, p3
    p3=p1/r4; print*, p3
    p3=p1/r8; print*, p3
    r4=pt1%Tplus(pt2); print*, r4
    r4=pt1%Tminus(pt2); print*, r4
    

    end program HW8


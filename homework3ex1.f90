module calculation
    
    integer :: n, i
        real :: r, x, dx
        doubleprecision :: SS, S
        real, allocatable :: a(:), y(:), d2ydx2(:), z(:), d2zdx2(:)

        
        
    contains
    
    
    doubleprecision function mean(a, n)
        implicit none
        integer, intent(in)::n
        real, intent(in)::a(n)
        real :: S=0
         do i=1, n
            S = S + a(i)
          end do
          mean = S/n
    end function mean
    
    doubleprecision function stand_dev(a, n)
        implicit none
        integer, intent(in)::n
        real, intent(in)::a(n)
        real :: S=0, SS=0
         do i=1, n
             S = S + a(i)
             SS = SS + a(i)**2
          end do
          stand_dev = sqrt(SS/n - (S/n)**2)
    end function stand_dev
    
    subroutine derivative2(a, np, h, adoubleprime)
    integer, intent(in) :: np
    real, intent(in) :: a(np), h
    real, intent(out):: adoubleprime (np)
    integer :: i
    
    do i=2, np-1
        adoubleprime (i) = (a(i+1)-2*a(i)+a(i-1))/h**2
    end do
     adoubleprime (1) = 0
     adoubleprime (np) = 0
     end subroutine derivative2
    
    end module calculation
    
    
    program hw3exercise1
    
    use calculation

    implicit none

 
        print*, "Please type how many numbers you want"
        read*, n
        if (n>0) then
        else
         n=0
        do while (n<1)
             read*, n
        end do
        end if
        allocate (a(n))
        print*, n
        print*, "Please type in a series of", n, " numbers"
        
        do i=1, n
             read*, r
             a(i) = r
        end do
        
print*, "The mean is", mean(a, n)
print*, "The standard deviation is", stand_dev(a, n)


write (*, '(a, $)') 'Input numbert of grid points:'
    read*, n
    allocate (y(n), d2ydx2(n), z(n), d2zdx2(n)) ! grid arrays
    
    dx=10.0/(n-1) ! grid spacing
    
    do i=1, n
        x=(i-1)*dx
         y(i)=sin(x)
         z(i)=x**2
    end do
    
    call derivative2 (y, n, dx, d2ydx2)
    call derivative2 (z, n, dx, d2zdx2)
    
    do i=1, n
        x=(i-1)*dx
        print*, d2ydx2(i), -sin(x), -sin(x)-d2ydx2(i), d2zdx2(i), 2, 2-d2zdx2(i)
    end do
    
    deallocate(y, d2ydx2, z, d2zdx2)

    end program hw3exercise1
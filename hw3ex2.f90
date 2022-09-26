module calculation
    
    integer :: n, i
        real :: r
        doubleprecision :: SS, S
        real, allocatable :: a(:)

        
        
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
    
    end module calculation
    
    
    
    program hw3ex2
    use calculation
    implicit none
    integer ::  b
    
    
    open (1, file='file.dat', status='old')
    do
        read (1, *, iostat=i) b
        if (i<0) exit
        if (i/=0) stop
        n=n+1
    end do
    print*, 'found', n, 'values'
    allocate (a(n))
    rewind (1)
    do i=1, n
        read (1, *) a(i)
    end do
    
print*, "The mean is", mean(a, n)
print*, "The standard deviation is", stand_dev(a, n)

    end program hw3ex2


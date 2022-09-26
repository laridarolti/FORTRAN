 program hw2exercise2
    implicit none
        integer :: n, i
        real :: r
        doubleprecision :: SS, S
        real, allocatable :: a(:)
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

    end program hw2exercise2


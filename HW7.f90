    module Poisson_solver
    contains

    recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
    implicit none
    real resV
    real,intent(inout):: u_f(:,:)  ! arguments
    real,intent(in)   :: rhs(:,:),h
    integer         :: nx,ny,nxc,nyc, i  ! local variables
    real,allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
    real            :: alpha=0.7, res_rms

    nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
    if( nx-1/=2*((nx-1)/2) .or. ny-1/=2*((ny-1)/2) ) &
        stop 'ERROR:not a power of 2'
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size

    if (min(nx,ny)>5) then  ! not the coarsest level
        allocate(res_f(nx,ny),corr_f(nx,ny), &
            corr_c(nxc,nyc),res_c(nxc,nyc))
        !---------- take 2 iterations on the fine grid--------------
        res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
        res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
        !---------- restrict the residue to the coarse grid --------
        call residue_2DPoisson(u_f,rhs,h,res_f)
        call restrict(res_f,res_c)
        !---------- solve for the coarse grid correction -----------
        corr_c = 0.
        res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2) ! *RECURSIVE CALL*
        !---- prolongate (interpolate) the correction to the fine grid
        call prolongate(corr_c,corr_f)
        !---------- correct the fine-grid solution -----------------
        u_f = u_f - corr_f
        !---------- two more smoothing iterations on the fine grid---
        res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
        res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
        deallocate(res_f,corr_f,res_c,corr_c)
    else
        !----- coarsest level (ny=5): iterate to get 'exact' solution
        do i = 1,100
            res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
        end do
    end if
    resV = res_rms   ! returns the rms. residue
    end function Vcycle_2DPoisson


    subroutine residue_2DPoisson(u_f, rhs, h, res_f)
    implicit none
    real, intent(in) :: u_f(:,:), rhs(:, :)
    real, intent(in):: h
    real, intent(out):: res_f(:, :)
    integer :: nx, ny, i, j
    nx=size(u_f, 1)
    ny=size(u_f, 2)
    res_f(:, :)=0
    do concurrent (i=2:nx-1, j=2:ny-1)
        res_f(i, j) = (u_f(i, j+1) + u_f(i, j-1) + u_f(i+1, j) + u_f(i-1, j)-4*u_f(i, j))/h**2 - rhs(i, j)
    end do
    end subroutine residue_2DPoisson


    subroutine restrict (res_f, res_c)
    implicit none
    real, intent(in) :: res_f(:, :)
    integer :: i, j, nx, ny, nxc, nyc, i_c, j_c ! indices for restricted grid
    real, allocatable, intent(out) :: res_c (:, :)
    nx=size(res_f, 1)
    ny=size(res_f, 2)
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size
    allocate (res_c(nxc, nyc))
    do concurrent (i=1:nx, j=1:ny)
        if (mod(i, 2)==1 .and. mod(j, 2)==1) then
            i_c = (i+1)/2
            j_c = (j+1)/2
            res_c(i_c, j_c) = res_f(i, j)
        end if
    end do
    end subroutine restrict


    subroutine prolongate (coarse, fine)
    implicit none
    real, intent(in) :: coarse(:, :)
    real, allocatable, intent(out) :: fine(:, :)
    integer :: nx, ny, i, j, nxc, nyc, i_c, j_c
    nxc=size(coarse, 1)
    nyc=size(coarse, 2)
    nx= 1+ (nxc-1)*2
    ny= 1+ (nyc-1)*2
    allocate(fine(nx, ny))
    
    fine(1:nx:2,1:ny:2)=coarse(:,:)
    do i=2,nx-1,2
       fine(i,1:ny:2)=0.5*(fine(i-1,1:ny:2)+fine(i+1,1:ny:2))
    end do
    do j=2,ny-1,2
       fine(:,j)=0.5*(fine(:,j-1)+fine(:,j+1))
    end do
    end subroutine prolongate


    function iteration_2DPoisson ( u_f, rhs, h, alpha) result (res_rms)
    implicit none
    real, intent(in) :: rhs(:,:)
    real, intent(inout) :: u_f(:, :)
    real, intent(in):: h, alpha
    real :: res_rms
    integer :: nx, ny, i, j
    real, allocatable :: res_f(:, :)
    nx=size(u_f, 1)
    ny=size(u_f, 2)
    allocate (res_f(nx, ny))
    res_f(:, :)=0
    !u_fnew(:, :)=0
    ! ITERATION
    do concurrent(i=2:nx-1, j=2:ny-1)
        res_f(i, j) = (u_f(i, j+1) + u_f(i, j-1) + u_f(i+1, j) + u_f(i-1, j)-4*u_f(i, j))/(h**2) - rhs(i, j)
        !u_fnew(i, j)=u_f(i, j)+ alpha*res_f(i, j)*(h**2)/4
    end do
    do concurrent(i=2:nx-1, j=2:ny-1)
        u_f(i, j)=u_f(i, j)+ alpha*res_f(i, j)*(h**2)/4
    end do
    ! CALCULATING RMS RESIDUE
    res_rms=0
    do i=1, nx
        do j=1, ny
            res_rms = res_rms + (res_f(i, j))**2
        end do
    end do
    res_rms = sqrt(res_rms/(nx*ny))
    end function iteration_2DPoisson
    end module Poisson_solver

    module advection_diffusion
    integer, private :: i, nx, ny
    real :: deltax, deltay
    contains


    subroutine derivative2(T, nx, ny, deltax, deltay, del2T)
    integer, intent(in) :: nx, ny
    real, intent(in) :: T(nx, ny), deltax, deltay
    real, intent(out):: del2T (nx, ny)
    integer :: i
    do concurrent (i=1:nx, j=1:ny)
        if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
            del2T (i, j) = 0
        else
            del2T (i, j) = (T(i-1, j)-2*T(i, j)+T(i+1, j))/deltax**2 + (T(i, j-1)-2*T(i, j)+T(i, j+1))/deltay**2
        end if
    end do
    end subroutine derivative2


    subroutine x_and_y_velocities(S, nx, ny, deltax, deltay, vx, vy)
    integer, intent(in) :: nx, ny
    real, intent(in) :: S(nx, ny), deltax, deltay
    real, intent(out):: vx (nx, ny), vy (nx, ny)
    integer:: i
    do concurrent (i=1:nx, j=1:ny)
        if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
            vx (i, j) = 0
            vy (i, j) = 0
        else
            vx (i, j) = (S(i, j+1)-S(i, j-1))/(2*deltay)
            vy (i, j) = (S(i-1, j)-S(i+1, j))/(2*deltax)
        end if
    end do
    end subroutine x_and_y_velocities


    subroutine vgradientT( nx, ny, deltax, deltay, vx, vy, T, vgradT)
    integer, intent(in) :: nx, ny
    real, intent(in) :: vx(nx, ny), vy(nx, ny), T(nx, ny), deltax, deltay
    real, intent(out):: vgradT (nx, ny)
    integer:: i
    real :: vx_dTdx(nx, ny), vy_dTdy(nx, ny)

    do concurrent (i=1:nx, j=1:ny)
        if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
            vx_dTdx(i, j)=0
            vy_dTdy(i, j)=0
        else
            ! calculating vx_dTdx
            if(vx(i, j)>0) then
                vx_dTdx(i, j) = vx(i, j) * (T(i, j)-T(i-1, j))/deltax
            else if (vx(i, j)<0) then
                vx_dTdx(i, j) = vx(i, j) * (T(i+1, j)-T(i, j))/deltax
            else
                vx_dTdx(i, j) = 0
            end if
            ! calculating vy_dTdy
            if(vy(i, j)>0) then
                vy_dTdy(i, j)=vy(i, j)*(T(i, j)-T(i, j-1))/deltay
            else if (vy(i, j)<0) then
                vy_dTdy(i, j)=vy(i, j)*(T(i, j+1)-T(i, j))/deltay
            else
                vy_dTdy(i, j)=0
            end if
        end if
        ! calculating vgradT
        vgradT(i, j) = vx_dTdx(i, j) + vy_dTdy(i, j)
    end do
    end subroutine vgradientT

    end module advection_diffusion







    program HW6
    use Poisson_solver
    use advection_diffusion

    implicit none

    !declaring variables
    integer :: nx=257, ny=65, i, j, niter=0
    logical :: multigrid
    real :: h, res_rms, f_rms=0., r, total_time=0.1, a_adv=0.4, a_dif=0.15, err=1.e-3, pi, dx, width
    real:: k=1., dt_dif, time, dt, dy, vx_max, vy_max, dt_adv, alpha=0.7
    double precision :: Ra=1.e6, Pr=0.01
    real, allocatable :: RHS(:,:), T(:, :), S(:, :), W(:, :), vx(:, :), vy(:, :), del2T(:, :), vgradT(:, :)
    real, allocatable :: dTdx(:, :), del2W(:, :), vgradW(:, :)
    character(len=15) :: Tinit='cosine'
    character(len=100) :: fname='lowPrandt_parameters_EMPTY.txt'


    ! reading control parameters using namelist and command argument
    if(command_argument_count()>0) &
        call GET_COMMAND_ARGUMENT(1, fname)
    print*, fname
    namelist / inputs / Pr, nx, ny, total_time, Ra, err, a_dif, a_adv, Tinit
    open (1, file=fname, status='old')
    read (1, inputs)
    close(1)
    write (*, inputs)
    
    
    ! initialising variables
    h=1./(ny-1); dx=h; dy=h; width=dx*nx
    dt_dif=a_dif*(h**2)/max(Pr, 1.0) ! k is diffusivity
    pi=DACOS(-1.D0) ! defining pi

    allocate (T(nx, ny), S(nx, ny), W(nx, ny), vx(nx, ny), vy(nx, ny), del2T(nx, ny), vgradT(nx, ny), RHS(nx, ny))
    allocate (dTdx(nx, ny), del2W(nx, ny), vgradW(nx, ny))

    ! INITIALISE T
    if (Tinit=='cosine') then
        do concurrent (i=1:nx, j=1:ny)
            T(i, j)=0.5*(1+cos(3*pi*(i-1)*dx/width))
        end do
    else
        call RANDOM_NUMBER(T)
    end if
    !print*, T(:, :)

    ! INITIALISE S, W, v
    S(:, :)=0;  W(:, :)=0;  vx(:, :)=0; vy(:, :)=0; time=0.
    

    ! TIME STEPPING
    do while(time<=total_time)

        ! CALCULATING S
        f_rms=0
        do i=1, nx
            do j=1, ny
                f_rms = f_rms + (W(i, j))**2 ! CALCULATE f_rms
                !print*, 'f_rms=', f_rms
            end do
        end do
        f_rms=sqrt(f_rms/(nx*ny))
        res_rms= iteration_2DPoisson ( S, W, h, alpha); !print*, "res_rms", res_rms
        S (:, ny) = 0; S(:, 1) = 0; S(1, :) = 0;  S(nx, :) = 0
        r =res_rms/f_rms; !print*, "r", r

        do while ( r>=err)
            res_rms= Vcycle_2DPoisson(S, W,h)
            r =res_rms/f_rms; !print*, 'r for S', r
            S (:, ny) = 0; S(:, 1) = 0; S(1, :) = 0;  S(nx, :) = 0
        end do


        ! CALCULATING V(vx vy)
        call x_and_y_velocities(S, nx, ny, dx, dy, vx, vy)

        ! CALCULATE dt_adv -> dt
        vx_max=maxval((abs(vx)))
        vy_max=maxval((abs(vy)))
        dt_adv=a_adv*min((h/vx_max), (h/vy_max))
        dt= min(dt_dif, dt_adv)

        ! CALCULATE del2(T) and vgradT, same for W
        call derivative2(T, nx, ny, dx, dy, del2T)
        call derivative2(W, nx, ny, dx, dy, del2W)

        ! CALCULATING dT/dx for updating W
        dTdx(:, :)=0
        do j=2, ny-1
            do i=2, nx-1
                    dTdx(i, j) = (T(i+1, j)-T(i-1, j))/(2*dx)
            end do
        end do

        call vgradientT( nx, ny, dx, dy, vx, vy, T, vgradT)
        call vgradientT( nx, ny, dx, dy, vx, vy, W, vgradW)

            ! UPDATE T field, W, and time
            T(:, :)=T(:, :)+ dt*(del2T(:, :) - vgradT(:, :))
            T (:, ny) = 0
            T(:, 1) = 1
            T(1, :)=T(2, :)
            T(nx, :)=T(nx-1, :)
            W(:, :)=W(:, :)+ dt*(Pr*del2W(:, :) - vgradW(:, :)-Ra*Pr*dTdx(:, :))
            time=time+dt
            niter=niter+1
            if (mod(niter, 100)==0) then
                print*, time
            end if
        end do
        print*,"Iterations finished", niter

        ! WRITING INTO FILE
        ! PRINTING T & S & W
        open(2, file='test1_T.txt', status='replace')
        do j=1, ny
            write(2, '(257(1pe16.7, " "))') T(:, j)
        end do
        close(2)
        open(3, file='test1_W.txt', status='replace')
        do j=1, ny
            write(3, '(257(1pe20.9, " "))') W(:, j)
        end do
        close(3)
        open(4, file='test1_S.txt', status='replace')
        do j=1, ny
            write(4, '(257(1pe20.9, " "))') S(:, j)
        end do
        close(4)

    end program HW6
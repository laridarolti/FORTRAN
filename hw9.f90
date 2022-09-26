    module everything
    implicit none

    private
    public :: temperature

    type temperature

        real :: h, dt, total_time=0.1, err=1.e-3, a_dif=0.15, a_adv=0.4, width, dt_dif
        integer :: nx=257, ny=65
        double precision :: Ra=1.e6, Pr=0.01
        real, allocatable :: T(:, :), W(:, :), S(:, :), vx(:, :), vy(:, :)
        character(len=30) :: Tinit='cosine'

    contains
    procedure,public, pass(this) :: initialise, velocity_solve, timestep, writing
    procedure,private, nopass :: derivative2, vgradientT
    procedure,private, nopass :: Vcycle_2DPoisson, residue_2DPoisson, restrict
    procedure,private, nopass :: prolongate, iteration_2DPoisson, x_and_y_velocities


    end type temperature
    contains

    ! UNCHANGED SUBROUTINES AND FUNCTIONS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
    implicit none
    real resV
    real,intent(inout):: u_f(:,:)  ! arguments
    real,intent(in)   :: rhs(:,:),h
    integer         :: nx,ny,nxc,nyc, i ! local variables
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
    integer :: nx, ny, nxc, nyc, i_c, j_c, i, j ! indices for restricted grid
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
    integer :: nx, ny, nxc, nyc, i_c, j_c, i, j
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

    ! ITERATION
    do concurrent(i=2:nx-1, j=2:ny-1)
        res_f(i, j) = (u_f(i, j+1) + u_f(i, j-1) + u_f(i+1, j) + u_f(i-1, j)-4*u_f(i, j))/(h**2) - rhs(i, j)
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


    subroutine derivative2(T, nx, ny, h, del2T)
    integer, intent(in) :: nx, ny
    integer :: i, j
    real, intent(in) :: T(nx, ny), h
    real, intent(out):: del2T (nx, ny)
    do concurrent (i=1:nx, j=1:ny)
        if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
            del2T (i, j) = 0
        else
            del2T (i, j) = (T(i-1, j)-2*T(i, j)+T(i+1, j))/h**2 + (T(i, j-1)-2*T(i, j)+T(i, j+1))/h**2
        end if
    end do
    end subroutine derivative2


    subroutine x_and_y_velocities(S, nx, ny, h, vx, vy)
    integer, intent(in) :: nx, ny
    integer :: i, j
    real, intent(in) :: S(nx, ny), h
    real, intent(out):: vx (nx, ny), vy (nx, ny)
    do concurrent (i=1:nx, j=1:ny)
        if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
            vx (i, j) = 0
            vy (i, j) = 0
        else
            vx (i, j) = (S(i, j+1)-S(i, j-1))/(2*h)
            vy (i, j) = (S(i-1, j)-S(i+1, j))/(2*h)
        end if
    end do
    end subroutine x_and_y_velocities


    subroutine vgradientT(nx, ny, h, vx, vy, T, vgradT)
    integer, intent(in) :: nx, ny
    integer :: i, j
    real, intent(in) :: vx(nx, ny), vy(nx, ny), T(nx, ny), h
    real, intent(out):: vgradT (nx, ny)
    real :: vx_dTdx(nx, ny), vy_dTdy(nx, ny)
    do concurrent (i=1:nx, j=1:ny)
        if (i==1 .or. j==1 .or. i==nx .or. j==ny) then
            vx_dTdx(i, j)=0
            vy_dTdy(i, j)=0
        else
            ! calculating vx_dTdx
            if(vx(i, j)>0) then
                vx_dTdx(i, j) = vx(i, j) * (T(i, j)-T(i-1, j))/h
            else if (vx(i, j)<0) then
                vx_dTdx(i, j) = vx(i, j) * (T(i+1, j)-T(i, j))/h
            else
                vx_dTdx(i, j) = 0
            end if
            ! calculating vy_dTdy
            if(vy(i, j)>0) then
                vy_dTdy(i, j)=vy(i, j)*(T(i, j)-T(i, j-1))/h
            else if (vy(i, j)<0) then
                vy_dTdy(i, j)=vy(i, j)*(T(i, j+1)-T(i, j))/h
            else
                vy_dTdy(i, j)=0
            end if
        end if
        ! calculating vgradT
        vgradT(i, j) = vx_dTdx(i, j) + vy_dTdy(i, j)
    end do
    end subroutine vgradientT


    ! NEW SUBROUTINES AND FUNCTIONS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! INITIALISING SUBROUTINE
    subroutine initialise (this)
    implicit none
    class(temperature), intent(inout) :: this
    !real, intent(out) :: dt_dif, time, put this into type temperature????
    real :: pi, Pr, total_time, Ra, err, a_dif, a_adv, time
    integer :: nx, ny, i, j
    character(len=30) :: Tinit, fname='convpars2.txt'

    ! reading control parameters using namelist and command argument
    if(command_argument_count()>0) &
        call GET_COMMAND_ARGUMENT(1, fname)
    print*, fname
    namelist / inputs / Pr, nx, ny, total_time, Ra, err, a_dif, a_adv, Tinit
    open (1, file=fname, status='old')
    read (1, inputs)
    close(1)
    write (*, inputs)
    this%Pr=Pr; this%nx=nx; this%ny=ny; this%total_time=total_time; this%Ra=Ra;
    this%err=err; this%a_dif=a_dif; this%a_adv=a_adv; this%Tinit=Tinit
    this%h=1./(this%ny-1); this%width=this%h*(this%nx-1)
    this%dt_dif=this%a_dif*(this%h**2)/max(this%Pr, 1.0) ! k is diffusivity
    pi=DACOS(-1.D0) ! defining pi
    allocate (this%T(this%nx,this%ny), this%S(this%nx, this%ny))
    allocate(this%W(this%nx, this%ny), this%vx(this%nx, this%ny), this%vy(this%nx, this%ny))

    ! INITIALISE T
    if (this%Tinit=='cosine') then
        do concurrent (i=1:this%nx, j=1:this%ny)
            this%T(i, j)=0.5*(1+cos(3*pi*(i-1)*this%h/this%width))
        end do
    else
        call RANDOM_NUMBER(this%T)
    end if
    ! INITIALISE S, W, V
    this%S(:, :)=0;  this%W(:, :)=0;  this%vx(:, :)=0; this%vy(:, :)=0
    end subroutine initialise


    ! SUBROUTINE TO SOLVE FOR VELOCITY
    subroutine velocity_solve(this)
    implicit none
    class(temperature),intent(inout) :: this
    real :: alpha=0.7
    real :: f_rms=0., r, res_rms
    integer :: i, j
    ! CALCULATING S
    do i=1, this%nx
        do j=1, this%ny
            f_rms = f_rms + (this%W(i, j))**2 ! CALCULATE f_rms
        end do
    end do
    f_rms=sqrt(f_rms/(this%nx*this%ny))
    res_rms= iteration_2DPoisson (this%S, this%W, this%h, alpha);
    this%S (:, this%ny) = 0; this%S(:, 1) = 0; this%S(1, :) = 0;  this%S(this%nx, :) = 0
    if (f_rms/=0.) then
        r =res_rms/f_rms;
        do while ( r>=this%err)
            res_rms= Vcycle_2DPoisson(this%S, this%W, this%h)
            r =res_rms/f_rms; !print*, 'r for S', r
            this%S (:, this%ny) = 0; this%S(:, 1) = 0; this%S(1, :) = 0;  this%S(this%nx, :) = 0
        end do
    end if

    ! CALCULATING V(vx vy)
    call x_and_y_velocities(this%S, this%nx, this%ny, this%h, this%vx, this%vy)
    end subroutine velocity_solve

    ! SUBROUTINE TO SOLVE FOR TIME-STEPPING
    subroutine timestep(this)
    implicit none
    class(temperature),intent(inout) :: this
    real:: k=1., vx_max, vy_max, dt_adv, alpha=0.7
    real :: del2W(this%nx, this%ny), del2T(this%nx, this%ny), dTdx(this%nx,this%ny)
    real :: vgradT(this%nx,this%ny), vgradW(this%nx,this%ny)
    integer :: i, j


    ! CALCULATE dt_adv -> dt
    vx_max=maxval((abs(this%vx)))
    vy_max=maxval((abs(this%vy)))
    if (vx_max/=0. .and. vy_max/=0.) then
        dt_adv=this%a_adv*min((this%h/vx_max), (this%h/vy_max))
        this%dt= min(this%dt_dif, dt_adv)
    else
        this%dt=this%a_dif*(this%h**2)/max(this%Pr, 1.0) ! dt = dt_dif
    end if

    ! CALCULATE del2(T) and vgradT, same for W
    call this%derivative2(this%T, this%nx, this%ny, this%h, del2T)
    call this%derivative2(this%W, this%nx, this%ny, this%h, del2W)

    ! CALCULATING dT/dx for updating W
    dTdx(:, :)=0
    do j=2, this%ny-1
        do i=2, this%nx-1
            dTdx(i, j) = (this%T(i+1, j)-this%T(i-1, j))/(2*this%h)
        end do
    end do

    call vgradientT( this%nx, this%ny, this%h, this%vx, this%vy, this%T, vgradT)
    call vgradientT( this%nx, this%ny, this%h, this%vx, this%vy, this%W, vgradW)

    ! UPDATE T field, W, and time
    this%T(:, :)=this%T(:, :)+ this%dt*(del2T(:, :) - vgradT(:, :))
    this%T (:, this%ny) = 0
    this%T(:, 1) = 1
    this%T(1, :)=this%T(2, :)
    this%T(this%nx, :)=this%T(this%nx-1, :)
    this%W(:, :)=this%W(:, :)+ this%dt*(this%Pr*del2W(:, :) - vgradW(:, :)-this%Ra*this%Pr*dTdx(:, :))
    end subroutine timestep


    ! WRITING SUBROUTINE
    subroutine writing (this)
    implicit none
    class(temperature), intent(inout) :: this
    integer :: i, j
    ! WRITING INTO FILE
    ! PRINTING T & S & W
    open(2, file='test1_T.txt', status='replace')
    do j=1, this%ny
        write(2, '(257(1pe16.7, " "))') this%T(:, j)
    end do
    close(2)
    open(3, file='test1_W.txt', status='replace')
    do j=1, this%ny
        write(3, '(257(1pe20.9, " "))') this%W(:, j)
    end do
    close(3)
    open(4, file='test1_S.txt', status='replace')
    do j=1, this%ny
        write(4, '(257(1pe20.9, " "))') this%S(:, j)
    end do
    close(4)

    end subroutine writing

    end module everything


    program hw9
    use iso_fortran_env
    use everything

    implicit none

    !declaring variables
    integer :: i, j, niter=0
    real :: time=0.
    real:: alpha=0.7
    real, allocatable :: RHS(:,:), vgradT(:, :), vgradW(:, :), dTdx(:, :), del2T(:, :), del2W(:, :)
    type(temperature) :: this



    ! INITIALISING VARIABLES
    call this%initialise()
    allocate (del2T(this%nx, this%ny), vgradT(this%nx,this%ny), RHS(this%nx, this%ny))
    allocate (dTdx(this%nx, this%ny), del2W(this%nx, this%ny), vgradW(this%nx, this%ny))

    ! TIME STEPPING
    do while(time<=this%total_time)
        call this%velocity_solve()
        call this% timestep()
        time=time+this%dt
        niter=niter+1
        if (mod(niter, 100)==0) then
            print*, time
        end if
    end do
    print*,"Iterations finished", niter

    ! WRITING INTO FILE
    ! PRINTING T & S & W
    call this%writing()

    end program hw9


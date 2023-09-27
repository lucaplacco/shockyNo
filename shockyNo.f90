!#####################################################################################
! PARAMETERS AND CONSTANTS
!#####################################################################################

module modshocky

     implicit none
     save
     integer, parameter :: singtype = selected_real_kind(6,37)    ! Single precision
     integer, parameter :: doubtype = selected_real_kind(15,307)  ! Double precision
     integer, parameter :: mykind   = doubtype

     integer :: nx,ng,itf,fsize,itmax,it = 0, nz = 1                                 ! Simulation discretization
     integer :: corder,dorder,methodT,iweno,iflux,itout = 1, vidout = 10            ! Discretization settings
     integer :: dataID = 10, err = 0, icase, bcase, xcase, isibm, ibm_ghs, idbg      ! Simulation settings
     integer :: fep_id

     character(7) :: iter                                     ! Internal timetick
     !character(4) :: chf

     real(mykind) :: ma,re,pr,gm,muinf                        ! Non-dimensional groups
     real(mykind) :: u0,rho0,p0,t0,l0,tme0                    ! Reference values
     real(mykind) :: uinf,pinf,rhoinf,rhouinf,rhoeinf,ubrgs   ! Inflow conditions
     real(mykind) :: uA,pA,rhoA,rhoUA,rhoEA                   ! Shocktube - init left side
     real(mykind) :: uB,pB,rhoB,rhoUB,rhoEB                   ! Shocktube - init right side

     real(mykind) :: l_tot,t_tot,t=0,dt,dx,dz,lz,xdf          ! Time and space domain details
     real(mykind) :: xwall,uw,ampw,frqw,dxr,dxl,dsigma,dalpha ! Interface and IBM distances
     real(mykind) :: pi=acos(-1.)                             ! Math constants

     real(mykind), dimension(3) :: wallpast
 
     real(mykind), dimension(:)  , allocatable :: x,z            ! Eulerian space x and z (PLOT3D)
     real(mykind), dimension(:,:), allocatable :: w,rhs_w,fc,fv  ! State vector and RHS, fluxes

     real(mykind), allocatable, dimension(:,:,:) :: qslice       ! PLOT3D solution array

end module modshocky


!#####################################################################################
! MAIN PROGRAM - v1.2 STABLE
!#####################################################################################

program shockyNo

    use modshocky
    implicit none

    write(*,*)
    write(*,*) 'ShockyNo v1.2'
    write(*,*)
    write(*,*) 'Reading input for the simulation ...'
    call readinput()
    write(*,*) 'Starting up ...'
    call startup()
    write(*,*) 'Initializing the flow domain ...'
    call init()

    ! DEBUG/NORUN routine
    if(idbg.ne.0) then
        call debug()
        if(idbg.ne.1) stop
    endif

    write(*,*)
    write(*,*) 'Simulation start ...'
    write(*,*)

    do while(it.lt.itmax)

        if(isnan(w(1,nx/2))) then
            write(*,*)
            write(*,*) 'Failed, stopping the simulation ...'
            write(*,*)
            exit
        endif

        !call  compute_rkt(w)       ! Shocktube integrator
        call  compute_dbg(w)        ! Walltube integrator

        it=it+1                 ! After RK we update the iteration (iterate from it to it+1)
        t=t+dt                  ! Time advances after integration

        if (mod(it,itout).eq.0) then
            ! Output video
            write(*,*) it,t,w(1,nx/2),w(2,nx/2),w(3,nx/2),xwall
        endif
        if (mod(it,vidout).eq.0) then
            call output_plot3d(w)
            write(*,*) 'Saving 2D fields...'
        endif

    enddo

    ! Save output at the final timestep
    write(*,*) 'Generating closing outputs ...'
    call output_dump(w)

    ! Compute and save solution file of the shocktube
    ! write(*,*) 'Generating analytical solutions from inputs ...'
    !call compute_ext(t)

    ! Shutting down and dealloacating
    write(*,*)
    write(*,*) 'Closing off ...'
    call shutdown()
    write(*,*) 'Bye'
    write(*,*)


end program shockyNo


!#####################################################################################
! SUBROUTINES
!#####################################################################################
! SUBROUTINE INIT
! Initialize the space and time domain / IBM module
!#####################################################################################

subroutine init()

    use modshocky
    implicit none

    integer :: i

    ! Flow reference quantities
    l0   = 1.0_mykind
    t0   = 1.0_mykind
    p0   = 1.0_mykind
    rho0 = 1.0_mykind          ! We fix T0,p0,rho0 to set R0 to 1
    u0   = sqrt(p0/rho0)
    tme0 = l0/u0

    ! Grid and domain init
    l_tot = l_tot/l0
    dx = l_tot/float(nx)
    ng = max(corder,dorder)/2
    dt = dt/tme0
    t_tot = dt*itmax
    ! Plot3D init
    lz = 2.0_mykind         ! Height of the 2D field
    dz = lz/float(nz)

    ! Grid allocation
    do i = -ng,nx+ng
        ! Allocate fluid nodes (Eulerian)
        x(i) = dx*i - (l_tot*0.5_mykind)
    enddo
    ! PLOT3D space allocation
    do i = 0,nz
        z(i) = dz*i
    enddo

    if(xcase.eq.0)  then

        !Basic inflow
        w(1,:) = 1.0_mykind
        w(2,:) = 0.0_mykind
        w(3,:) = 1.0_mykind/(gm-1)

        uinf = sqrt(gm)*Ma
        rhouinf = rhoinf*uinf
        rhoeinf = (0.5_mykind*(uinf**2) + ((pinf/rhoinf)/(gm-1)))*rhoinf

        do i = 0,ng
            w(1,-i) = rhoinf
            w(2,-i) = rhouinf
            w(3,-i) = rhoeinf
        enddo

        open(unit=11, file='init.txt')
        do i = -ng,nx+ng
        write(11,*) i,x(i),w(1,i),w(2,i),w(3,i)
        enddo
        close(11)

    elseif((xcase.eq.1).or.(xcase.eq.4))  then

        ! Shocktube inflow
        ua    = ua/u0
        ub    = ub/u0
        pa    = pa/p0
        pb    = pb/p0
        rhoa  = rhoa/rho0
        rhob  = rhob/rho0

        rhoua = rhoa*ua
        rhoea = (0.5_mykind*(ua**2) + ((pa/rhoa)/(gm-1)))*rhoa
        rhoub = rhob*ub
        rhoeb = (0.5_mykind*(ub**2) + ((pb/rhob)/(gm-1)))*rhob

        open(unit=11,file='init.txt')
        do i = -ng,nx+ng

             ! Shu-Osher Shocktube
             !rhob = 1.0_mykind + 0.2_mykind*sin(5.0_mykind*(dx*(float(i)))-0.5_mykind*l_tot)
             !rhoub = rhoa*ub
             !rhoeb = (0.5_mykind*(ub**2) + ((pb/rhob)/(gm-1)))*rhob

             if(i.le.(xdf*nx)) then
                w(1,i) = rhoa
                w(2,i) = rhoua
                w(3,i) = rhoea
             else
                w(1,i) = rhob
                w(2,i) = rhoub
                w(3,i) = rhoeb
             endif

             if(xcase.eq.4) then
                if(x(i).le.xwall) then
                    w(1,i) = rhoa
                    w(2,i) = rhoua
                    w(3,i) = rhoea
                else
                    w(1,i) = rhob
                    w(2,i) = rhoub
                    w(3,i) = rhoeb
                endif
             endif

             write(11,*) i,x(i),w(1,i),w(2,i),w(3,i)

        enddo
        close(11)

    elseif(xcase.eq.2)  then

        ! Periodic signal
        uinf    = 0.05_mykind
        pinf    = 1.00_mykind

        open(unit=11,file='init.txt')
        do i = -ng,nx+ng

            rhoinf  = 1.0_mykind+0.2_mykind*sin(pi*dx*(float(i)))
            rhouinf = uinf*rhoinf
            rhoeinf = (0.5_mykind*(uinf**2) + ((pinf/rhoinf)/(gm-1)))*rhoinf
            w(1,i) = rhoinf
            w(2,i) = rhouinf
            w(3,i) = rhoeinf
            write(11,*) i,x(i),w(1,i),w(2,i),w(3,i)

        enddo
        close(11)

    else

        ! Burgers advection/Scalar discontinuity
        open(unit=11,file='init.txt')
        do i = -ng,nx+ng

            w(1,i) = 1.0_mykind+0.223_mykind*sin(pi*dx*(float(i)))    ! Phi
            w(2,i) = 0.0_mykind
            w(3,i) = 0.0_mykind

            ! if(i.le.(xdf*nx)) then
            !     w(1,i) = rhoa
            ! else
            !     w(1,i) = rhob
            ! endif
            ! w(2,i) = 0.0_mykind
            ! w(3,i) = 0.0_mykind

            write(11,*) i,x(i),w(1,i),w(2,i),w(3,i)
            !write(11,*) i,dx*(float(i)),w(1,i),w(2,i),w(3,i)

        enddo
        close(11)

    endif

    ! Init IBM module
    if(isibm.eq.1) then
        wallpast = 100.0_mykind
        wallpast(1) = xwall
        ampw = uw
    endif

    ! Plot3D init output
    call output_plot3d(w)

end subroutine init


!#####################################################################################
! SUBROUTINE COMPUTE_EUL
! Compute convective fluxes as output (WENO RECONSTRUCTION), state vector as input
!#####################################################################################

subroutine compute_eul(wr,fcr)

    use modshocky
    implicit none

    integer :: i,j,k,l

    real(mykind) :: alfa,eps
    real(mykind) :: uavg,ravg,aavg,eavg,havg,b1,b2,ai           ! Roe/average interface states

    real(mykind), dimension(0:1)       :: dw
    real(mykind), dimension(3,0:1)     :: aw,bw,ww     ! WENO coefficients  
    real(mykind), dimension(-ng:nx+ng) :: u,p,d,a      ! Internal flow quantities

    real(mykind), dimension(3,3) :: Q,Q1                    ! Eigenvector matrix (charateristic dir.)
    real(mykind), dimension(1:3) :: sigma                   ! Eigenvalue diag. matrix
    real(mykind), dimension(1:3) :: sigmax                  ! Maximum diag. matrix

    real(mykind), dimension(3,-ng:nx+ng)   :: wr                  ! Flow state
    real(mykind), dimension(3,-ng:nx+ng)   :: fcc,fcc_up,fcc_dw   ! Central convective fluxes 

    real(mykind), dimension(3)             :: fcc_pj,wc_pj                      ! Projected quantities
    real(mykind), dimension(3)             :: fcw_up,fcw_dw,gcw_up,gcw_dw       ! WENO reconstruction fluxes
    real(mykind), dimension(3,0:1)         :: sup,sdw                           ! Stencil interpolation
    real(mykind), dimension(3,0:fsize)     :: fcr,fcr_up,fcr_dw                 ! Interface convective fluxes

    fcr = 0.0_mykind
    fcr_up = 0.0_mykind
    fcr_dw = 0.0_mykind

    fcc = 0.0_mykind
    fcc_up = 0.0_mykind
    fcc_dw = 0.0_mykind
    fcc_pj = 0.0_mykind
    wc_pj = 0.0_mykind

    Q = 0.0_mykind
    Q1 = 0.0_mykind

    alfa = 0.0_mykind
    sigma = 0.0_mykind

    aw = 0.0_mykind
    bw = 0.0_mykind
    ww = 0.0_mykind

    uavg = 0.0_mykind
    aavg = 0.0_mykind
    ravg = 0.0_mykind
    havg = 0.0_mykind
    eavg = 0.0_mykind

    dw  = [(1.0/3.0),(2.0/3.0)]       ! WENO3 linear coefficients
    eps = 1.0E-12_mykind

    u(:) = wr(2,:)/wr(1,:)
    d(:) = (u(:)**2)*0.5_mykind
    p(:) = (gm-1)*(wr(3,:)-wr(1,:)*d(:))
    a(:) = sqrt((gm*p(:))/wr(1,:))


    ! Flux splitting at node position
    if(iweno.ne.0) then

        if(xcase.eq.3) then

            fcc(1,:) = wr(1,:)*ubrgs
            if(ubrgs.ge.0) then
                fcc_up(1,:) = fcc(1,:)
            else
                fcc_dw(1,:) = fcc(1,:)
            endif

        else

            ! Compute central convective fluxes at nodes for WENO stencil interpolation
            fcc(1,:) = wr(2,:)
            fcc(2,:) = (wr(2,:) * u(:)) + p(:)
            fcc(3,:) = (wr(3,:) + p(:)) * u(:)

            if(iflux.eq.0) then
                ! Component-wise flux splitting Lax-Friedrichs for central fluxes
                do i=-ng,nx+ng
                    sigma(1) = abs(u(i))
                    sigma(2) = abs(u(i) + a(i))
                    sigma(3) = abs(u(i) - a(i))
                    alfa = maxval(sigma)
                    fcc_up(:,i) = 0.5_mykind*(fcc(:,i) + alfa*wr(:,i))    ! Positive flux
                    fcc_dw(:,i) = 0.5_mykind*(fcc(:,i) - alfa*wr(:,i))    ! Negative flux
                enddo
            endif

        endif

    endif


    ! Convective fluxes computation
    do i = 0,fsize

        if(xcase.ne.3) then 

            ! Complete convective fluxes computing at the cell interface
            fcr(1,i) = (wr(2,i) + wr(2,i-1))*0.5_mykind
            fcr(2,i) = (((wr(2,i) * u(i)) + p(i)) + ((wr(2,i-1) * u(i-1)) + p(i-1)))*0.5_mykind
            fcr(3,i) = (((wr(3,i) + p(i)) * u(i)) + ((wr(3,i-1) + p(i-1)) * u(i-1)))*0.5_mykind

            if(iweno.ne.0) then     ! Fluxes are corrected using WENO

                if((i.eq.0).or.(i.eq.fsize)) cycle      ! First and last flux are excluded
                !if((i.eq.itf).or.(i.eq.(itf+1))) cycle  ! Trans-interface fluxes are excluded (IBM)

                j = i - 1                               ! Sync flux index (i) with node index (j)

                if(iflux.ne.0) then                     ! Compute eigenvectors for characteristic-wise flux splitting

                    ! Average state at cell interface
                    uavg = (u(j)+u(j+1))*0.5_mykind
                    aavg = (a(j)+a(j+1))*0.5_mykind
                    ravg = (wr(1,j)+wr(1,j+1))*0.5_mykind
                    eavg = 0.5_mykind*(uavg**2)
                    havg = eavg + (aavg**2)/(gm-1)

                    ! Compute eigenvectors matrix
                    ai = 1/aavg
                    b2 = (gm-1)/(aavg**2)
                    b1 = b2*eavg

                    Q(1,1) = 1.0_mykind
                    Q(1,2) = 1.0_mykind
                    Q(1,3) = 1.0_mykind
                    Q(2,1) = uavg-aavg
                    Q(2,2) = uavg
                    Q(2,3) = uavg+aavg
                    Q(3,1) = havg-aavg*uavg
                    Q(3,2) = eavg
                    Q(3,3) = havg+aavg*uavg

                    Q1(1,1) =  0.5_mykind*(b1 + uavg*ai)
                    Q1(1,2) = -0.5_mykind*(b2*uavg + ai)
                    Q1(1,3) =  0.5_mykind*b2
                    Q1(2,1) =  1.0_mykind - b1
                    Q1(2,2) =  uavg*b2
                    Q1(2,3) = -b2
                    Q1(3,1) =  0.5_mykind*(b1 - uavg*ai)
                    Q1(3,2) = -0.5_mykind*(b2*uavg - ai)
                    Q1(3,3) =  0.5_mykind*b2

                    ! Local maximum eigenvalues over the WENO stencil
                    sigmax(:) = -100.0_mykind
                    do k = j-1,j+2
                        sigma(1) = (u(k) - a(k))
                        sigma(2) = (u(k))
                        sigma(3) = (u(k) + a(k))
                        do l = 1,3
                            sigmax(l) = max(abs(sigma(l)),sigmax(l))
                        enddo
                    enddo
 
                    ! Characteristic-wise flux splitting
                    do k = j-1,j+2

                        ! Project fluxes (only the ones involved in the stencil) along the characteristics
                        fcc_pj(1) = fcc(1,k)*Q1(1,1) + fcc(2,k)*Q1(1,2) + fcc(3,k)*Q1(1,3)
                        fcc_pj(2) = fcc(1,k)*Q1(2,1) + fcc(2,k)*Q1(2,2) + fcc(3,k)*Q1(2,3)
                        fcc_pj(3) = fcc(1,k)*Q1(3,1) + fcc(2,k)*Q1(3,2) + fcc(3,k)*Q1(3,3)
                        ! Project local flow state
                        wc_pj(1) = wr(1,k)*Q1(1,1) + wr(2,k)*Q1(1,2) + wr(3,k)*Q1(1,3)
                        wc_pj(2) = wr(1,k)*Q1(2,1) + wr(2,k)*Q1(2,2) + wr(3,k)*Q1(2,3)
                        wc_pj(3) = wr(1,k)*Q1(3,1) + wr(2,k)*Q1(3,2) + wr(3,k)*Q1(3,3)

                        ! Split the projected fluxes
                        fcc_up(:,k) = 0.5_mykind*(fcc_pj(:) + sigmax(:)*wc_pj(:))    ! Positive flux
                        fcc_dw(:,k) = fcc_pj(:) - fcc_up(:,k)                        ! Negative flux

                    enddo

                endif

                ! START OF THE WENO ALGORITHM

                ! POSITIVE flux treatment
                ! Stencil computation
                sup(:,0) = (-0.5_mykind)*fcc_up(:,j-1) + (3.0_mykind/2.0_mykind)*fcc_up(:,j)
                sup(:,1) = (0.5_mykind)*fcc_up(:,j) + (0.5_mykind)*fcc_up(:,j+1)
                ! Smoothness coefficients computation
                bw(:,0) = (fcc_up(:,j-1) - fcc_up(:,j))**2
                bw(:,1) = (fcc_up(:,j) - fcc_up(:,j+1))**2
                ! Interpolation weights computation
                aw(:,0) = dw(0)/((bw(:,0) + eps)**2)
                aw(:,1) = dw(1)/((bw(:,1) + eps)**2)
                ww(:,0) = aw(:,0)/(aw(:,0) + aw(:,1))
                ww(:,1) = aw(:,1)/(aw(:,0) + aw(:,1))
                ! Flux convex interpolation
                fcw_up(:) = sup(:,0)*ww(:,0) + sup(:,1)*ww(:,1)

                ! NEGATIVE flux treatment
                ! Stencil computation
                sdw(:,0) = (-0.5_mykind)*fcc_dw(:,j+2) + (3.0_mykind/2.0_mykind)*fcc_dw(:,j+1)
                sdw(:,1) = (0.5_mykind)*fcc_dw(:,j) + (0.5_mykind)*fcc_dw(:,j+1)
                ! Smoothness coefficients computation
                bw(:,0) = (fcc_dw(:,j+1) - fcc_dw(:,j+2))**2
                bw(:,1) = (fcc_dw(:,j) - fcc_dw(:,j+1))**2
                ! Interpolation weights computation
                aw(:,0) = dw(0)/((bw(:,0) + eps)**2)
                aw(:,1) = dw(1)/((bw(:,1) + eps)**2)
                ww(:,0) = aw(:,0)/(aw(:,0)+aw(:,1))
                ww(:,1) = aw(:,1)/(aw(:,0)+aw(:,1))
                ! Flux convex interpolation
                fcw_dw(:) = sdw(:,0)*ww(:,0) + sdw(:,1)*ww(:,1)

                ! Project back to physical space
                if(iflux.ne.0) then
                    gcw_up = fcw_up
                    gcw_dw = fcw_dw
                    fcw_up(1) = gcw_up(1)*Q(1,1) + gcw_up(2)*Q(1,2) + gcw_up(3)*Q(1,3)
                    fcw_up(2) = gcw_up(1)*Q(2,1) + gcw_up(2)*Q(2,2) + gcw_up(3)*Q(2,3)
                    fcw_up(3) = gcw_up(1)*Q(3,1) + gcw_up(2)*Q(3,2) + gcw_up(3)*Q(3,3)
                    fcw_dw(1) = gcw_dw(1)*Q(1,1) + gcw_dw(2)*Q(1,2) + gcw_dw(3)*Q(1,3)
                    fcw_dw(2) = gcw_dw(1)*Q(2,1) + gcw_dw(2)*Q(2,2) + gcw_dw(3)*Q(2,3)
                    fcw_dw(3) = gcw_dw(1)*Q(3,1) + gcw_dw(2)*Q(3,2) + gcw_dw(3)*Q(3,3)
                endif

                ! WENO flux-splitting reconstruction
                fcr(1,i) = fcw_up(1)+fcw_dw(1)
                fcr(2,i) = fcw_up(2)+fcw_dw(2)
                fcr(3,i) = fcw_up(3)+fcw_dw(3)

            endif


        else ! Burgers equation case


            ! Burgers equation convective fluxes
            fcr(1,i) = (wr(1,i)+wr(1,i-1))*ubrgs*0.5_mykind
            fcr(2,i) = 0.0_mykind
            fcr(3,i) = 0.0_mykind

            if(iweno.ne.0) then     ! Fluxes are corrected using WENO

                if((i.eq.0).or.(i.eq.fsize)) cycle      ! First and last flux are excluded
                j = i - 1                               ! Sync flux index (i) with WENO node index (j)
                if(ubrgs.ge.0) then

                    ! Burgers equation - WENO reconstruction
                    ! POSITIVE flux treatment
                    ! Stencil computation
                    sup(:,0) = (-0.5_mykind)*fcc_up(:,j-1) + (3.0_mykind/2.0_mykind)*fcc_up(:,j)
                    sup(:,1) = (0.5_mykind)*fcc_up(:,j) + (0.5_mykind)*fcc_up(:,j+1)
                    ! Smoothness coefficients computation
                    bw(:,0) = (fcc_up(:,j-1) - fcc_up(:,j))**2
                    bw(:,1) = (fcc_up(:,j) - fcc_up(:,j+1))**2
                    ! Interpolation weights computation
                    aw(:,0) = dw(0)/((bw(:,0) + epsilon(pi))**2)
                    aw(:,1) = dw(1)/((bw(:,1) + epsilon(pi))**2)
                    ww(:,0) = aw(:,0)/(aw(:,0) + aw(:,1))
                    ww(:,1) = aw(:,1)/(aw(:,0) + aw(:,1))
                    ! Flux convex interpolation
                    fcw_up(:) = sup(:,0)*ww(:,0) + sup(:,1)*ww(:,1)

                    ! NEGATIVE flux treatment
                    fcw_dw(:) = 0.0_mykind

                else

                    ! NEGATIVE flux treatment
                    ! Stencil computation
                    sdw(:,0) = (-0.5_mykind)*fcc_dw(:,j+2) + (3.0_mykind/2.0_mykind)*fcc_dw(:,j+1)
                    sdw(:,1) = (0.5_mykind)*fcc_dw(:,j) + (0.5_mykind)*fcc_dw(:,j+1)
                    ! Smoothness coefficients computation
                    bw(:,0) = (fcc_dw(:,j+1)-fcc_dw(:,j+2))**2
                    bw(:,1) = (fcc_dw(:,j)-fcc_dw(:,j+1))**2
                    ! Interpolation weights computation
                    aw(:,0) = dw(0)/((bw(:,0) + epsilon(pi))**2)
                    aw(:,1) = dw(1)/((bw(:,1) + epsilon(pi))**2)
                    ww(:,0) = aw(:,0)/(aw(:,0) + aw(:,1))
                    ww(:,1) = aw(:,1)/(aw(:,0) + aw(:,1))
                    ! Flux convex interpolation
                    fcw_dw(:) = sdw(:,0)*ww(:,0) + sdw(:,1)*ww(:,1)

                    ! POSITIVE flux treatment
                    fcw_up(:) = 0.0_mykind

                endif

                fcr(1,i) = fcw_up(1)+fcw_dw(1)
                fcr(2,i) = 0.0_mykind
                fcr(3,i) = 0.0_mykind

            endif


        endif

    enddo


end subroutine compute_eul


!#####################################################################################
! SUBROUTINE COMPUTE_VIS
! Compute viscous fluxes as output, state vector as input
!#####################################################################################

subroutine compute_vis(wr,fvr)

    use modshocky
    implicit none

    integer :: i

    real(mykind) :: flux_mu,flux_tt                                         ! Internal fluxes
    real(mykind), dimension(-ng:nx+ng)   :: u,p,q,h,mu,lmbd                 ! Internal flow quantities

    real(mykind), dimension(3,-ng:nx+ng)   :: wr
    real(mykind), dimension(3,0:fsize)     :: fvr

    fvr = 0.0_mykind

    if(icase.eq.1) then

        u(:) = wr(2,:)/wr(1,:)
        q(:) = (u(:)**2)*0.5_mykind
        p(:) = (gm-1)*(wr(3,:)-wr(1,:)*q(:))
        h(:) = p(:)/wr(1,:)                         ! Temperature

        ! Power-Law
        mu(:) = h(:)**(0.75_mykind)

        ! Thermal conductivity
        lmbd(:) = (gm/(gm-1))*(1/Pr)*(mu(:)/wr(1,:))

        do i = 0,fsize

            ! Viscous fluxes computation ( [1/3] constant omitted!)
            fvr(1,i) = 0.0_mykind
            fvr(2,i) = ((mu(i)+mu(i-1))*0.5_mykind) * ((u(i)-u(i-1))/dx)
            flux_mu  = ((mu(i)+mu(i-1))*0.5_mykind) * ((u(i)+u(i-1))*0.5_mykind) * ((u(i)-u(i-1))/dx)
            flux_tt  = (((lmbd(i)+lmbd(i-1))*0.5_mykind)*gm)/(Pr*(gm-1)) * ((h(i)-h(i-1))/dx)
            fvr(3,i) = flux_mu + flux_tt

        enddo

    endif

end subroutine compute_vis


!#####################################################################################
! SUBROUTINE COMPUTE_IBM
! Current fluxes as input and corrected IBM forces as output - version 2
!#####################################################################################

subroutine compute_ibm(wr,fcr,ibm_wr)

    use modshocky
    implicit none

    integer :: i

    real(mykind) :: pw,tw,rw,ew,xr,xl,xw
    real(mykind), dimension(3) :: fw,ff,f_itf_r,f_itf_l

    real(mykind), dimension(-ng:nx+ng)    :: u,h,p,d
    real(mykind), dimension(3,-ng:nx+ng)  :: wr,ibm_wr
    real(mykind), dimension(3,0:fsize)    :: fcr

    u(:) = wr(2,:)/wr(1,:)
    d(:) = (u(:)**2)*0.5_mykind
    p(:) = (gm-1)*(wr(3,:)-wr(1,:)*d(:))
    h(:) = p(:)/wr(1,:)

    ibm_wr = 0.0_mykind
    xw = xwall

    if(ibm_ghs.eq.1) then

        ! Right-halfside interface
        ! [Left wall (itf)]
        pw = p(itf)
        tw = h(itf)
        rw = pw/tw
        ew = (0.5_mykind*(uw**2) + ((pw/rw)/(gm-1)))*rw

        fw(1) = (rw*uw)
        fw(2) = (rw*uw*uw + pw)
        fw(3) = ((ew + pw)*uw)

        ff(1) = (wr(1,itf)*u(itf))
        ff(2) = (wr(1,itf)*((u(itf)**2)) + p(itf))
        ff(3) = ((wr(3,itf) + p(itf))*u(itf))

        ! Interpolate
        f_itf_l(1) = (fw(1)*(0.5_mykind*dx) + ff(1)*dsigma)/dxl
        f_itf_l(2) = (fw(2)*(0.5_mykind*dx) + ff(2)*dsigma)/dxl
        f_itf_l(3) = (fw(3)*(0.5_mykind*dx) + ff(3)*dsigma)/dxl

        ! [Right wall (itf+1)]
        pw = p(itf+1)
        tw = h(itf+1)
        rw = pw/tw
        ew = (0.5_mykind*(uw**2) + ((pw/rw)/(gm-1)))*rw
        xr = 0.5_mykind*(x(itf) + x(itf+1))

        fw(1) = (rw*uw)
        fw(2) = (rw*uw*uw + pw)
        fw(3) = ((ew + pw)*uw)

        ff(1) = (wr(1,itf+1)*u(itf+1))
        ff(2) = (wr(1,itf+1)*((u(itf+1)**2)) + p(itf+1))
        ff(3) = ((wr(3,itf+1) + p(itf+1))*u(itf+1))

        ! Extrapolate
        f_itf_r(1) = fw(1) + ((xr - xw)/(x(itf+1) - xw))*(ff(1) - fw(1))
        f_itf_r(2) = fw(2) + ((xr - xw)/(x(itf+1) - xw))*(ff(2) - fw(2))
        f_itf_r(3) = fw(3) + ((xr - xw)/(x(itf+1) - xw))*(ff(3) - fw(3))

        ! Compute forcing-correcting terms
        ibm_wr(:,itf) = (fcr(:,itf+1)/dx) - (f_itf_l/dx)
        ibm_wr(:,itf+1) = -(fcr(:,itf+1)/dx) + (f_itf_r/dx)

    else

        ! Left-halfside interface
        ! [Right wall (itf+1)]
        pw = p(itf+1)
        tw = h(itf+1)
        rw = pw/tw
        ew = (0.5_mykind*(uw**2) + ((pw/rw)/(gm-1)))*rw

        fw(1) = (rw*uw)
        fw(2) = (rw*uw*uw + pw)
        fw(3) = ((ew + pw)*uw)

        ff(1) = (wr(1,itf+1)*u(itf+1))
        ff(2) = (wr(1,itf+1)*((u(itf+1)**2)) + p(itf+1))
        ff(3) = ((wr(3,itf+1) + p(itf+1))*u(itf+1))

        ! Interpolate
        f_itf_r(1) = (fw(1)*(0.5_mykind*dx) + ff(1)*dsigma)/dxr
        f_itf_r(2) = (fw(2)*(0.5_mykind*dx) + ff(2)*dsigma)/dxr
        f_itf_r(3) = (fw(3)*(0.5_mykind*dx) + ff(3)*dsigma)/dxr

        ! [Left wall (itf)]
        pw = p(itf)
        tw = h(itf)
        rw = pw/tw
        ew = (0.5_mykind*(uw**2) + ((pw/rw)/(gm-1)))*rw
        xl = 0.5_mykind*(x(itf) + x(itf+1))

        fw(1) = (rw*uw)
        fw(2) = (rw*uw*uw + pw)
        fw(3) = ((ew + pw)*uw)

        ff(1) = (wr(1,itf)*u(itf))
        ff(2) = (wr(1,itf)*((u(itf)**2)) + p(itf))
        ff(3) = ((wr(3,itf) + p(itf))*u(itf))

        ! Extrapolate
        f_itf_l(1) = ff(1) + ((xl - x(itf))/(xw - x(itf)))*(fw(1) - ff(1))
        f_itf_l(2) = ff(2) + ((xl - x(itf))/(xw - x(itf)))*(fw(2) - ff(2))
        f_itf_l(3) = ff(3) + ((xl - x(itf))/(xw - x(itf)))*(fw(3) - ff(3))

        ! Compute forcing-correcting terms
        ibm_wr(:,itf) = (fcr(:,itf+1)/dx) - (f_itf_l/dx)
        ibm_wr(:,itf+1) = -(fcr(:,itf+1)/dx) + (f_itf_r/dx)

    endif


end subroutine compute_ibm


!#####################################################################################
! SUBROUTINE COMPUTE_ITF
! Find closest flow nodes to the IBM interface and update its position
!#####################################################################################

subroutine compute_itf(int_dt,wr)

    use modshocky
    implicit none

    integer :: i
    real(mykind) :: pw,rw,tw,ew
    real(mykind) :: clst,deltaitf,int_dt,int_dx,dw

    real(mykind), dimension(3) :: ww,xp,wp
    real(mykind), dimension(3,-ng:nx+ng)  :: wr

    ! Update interface position
    int_dx = abs(xwall - wallpast(1))     ! Wall movement during this STS and the previous

    ! Evaluate flow-emerging-points - ITF/ITF+1 are still related to the previous TS position
    if((int_dx.gt.wallpast(3)).and.(uw.gt.0)) then

        ! ITF+1 is a flow emerging point (in the flow opposite direction)
        pw = (gm-1)*(wr(3,itf)-wr(1,itf)*(((wr(2,itf)/wr(1,itf))**2)*0.5_mykind))
        tw = pw/wr(1,itf)
        rw = pw/tw
        ew = (0.5_mykind*(uw**2) + ((pw/rw)/(gm-1)))
        dw = abs(int_dx - wallpast(3))
        ww(1) = rw
        ww(2) = rw*uw
        ww(3) = rw*ew

        ! Correct the flow state at the emerging point
        !write(*,*) 'ITF jump forward'
        wr(:,itf+1) = (ww(:)*dx + wr(:,itf)*dw)/(dx + dw)              ! Interpolate from wall
        ! wr(:,itf+1) = wr(:,itf-1) + ((x(itf+1) - x(itf-1))/(x(itf) - x(itf-1)))*(wr(:,itf-1) - wr(:,itf)) ! Extrap.

    elseif((int_dx.gt.wallpast(2)).and.(uw.lt.0)) then

        ! ITF is a flow emerging point (in the flow opposite direction)
        pw = (gm-1)*(wr(3,itf+1)-wr(1,itf+1)*(((wr(2,itf+1)/wr(1,itf+1))**2)*0.5_mykind))
        tw = pw/wr(1,itf+1)
        rw = pw/tw
        ew = (0.5_mykind*(uw**2) + ((pw/rw)/(gm-1)))
        dw = abs(int_dx - wallpast(2))
        ww(1) = rw
        ww(2) = rw*uw
        ww(3) = rw*ew

        ! Correct the flow state at the emerging point
        !write(*,*) 'ITF jump backward'
        wr(:,itf) = (ww(:)*dx + wr(:,itf+1)*dw)/(dx + dw)              ! Interpolate from wall
        !wr(:,itf+1) = wr(:,itf-1) + ((x(itf+1) - x(itf-1))/(x(itf) - x(itf-1)))*(wr(:,itf-1) - wr(:,itf)) ! Extrap.

    endif


    ! Find IBM interface
    clst = l_tot
    if(xwall.gt.x(nx).or.(xwall.lt.(x(0)))) then
        write(*,*) 'IBM interface is out of boundaries!'
        stop
    endif
    do i = -ng,nx-ng
        ! Find interface relative position
        deltaitf = abs(x(i) - xwall)
        if(deltaitf.lt.clst) then
            itf = i
            clst = deltaitf
        endif
    enddo

    ! Update interface position (Eulerian)
    if(xwall.lt.x(itf)) then
        ! Ind itf corresponds to the flow index right to the interface; switch itf index
        itf = itf - 1
    endif
    dxr = abs(xwall - x(itf+1))  ! Right delta to closest node
    dxl = abs(xwall - x(itf))    ! Left delta to closest node

    if(dxl.gt.dxr) then
        ! Interface node is on the right-halfside
        ibm_ghs = 1
        dsigma = abs(dxl - 0.5_mykind*dx)
        dalpha = abs(dxr + 0.5_mykind*dx)
    else
        ! Interface node is on the left-halfside
        ibm_ghs = 0
        dsigma = abs(dxr - 0.5_mykind*dx)
        dalpha = abs(dxl + 0.5_mykind*dx)
    endif

    wallpast(1) = xwall      ! Current (previous, at the beginning of the subroutine) interface position
    wallpast(2) = dxl        ! Current (previous, at the beginning of the subroutine) relative distances
    wallpast(3) = dxr

end subroutine compute_itf


!#####################################################################################
! SUBROUTINE COMPUTE_RHS
! Current fluxes as input and rhs_f as output
!#####################################################################################

subroutine compute_rhs(fcr,fvr,ibm_wr,rhs_wr)

    use modshocky
    implicit none

    integer :: i
    real(mykind), dimension(3,0:fsize)    :: fcr,fvr
    real(mykind), dimension(3,-ng:nx+ng)  :: rhs_wr,ibm_wr

    rhs_wr = 0.0_mykind

    if(icase.eq.0) then
        do i=0,nx
            rhs_wr(:,i) = -(fcr(:,i+1)-fcr(:,i))/dx + ibm_wr(:,i)
            write(10,*) rhs_wr(:,i)
        enddo
    else
        if(xcase.eq.0) then
            do i=0,nx
                rhs_wr(:,i) = -(fcr(:,i+1)-fcr(:,i))/dx + (ma/re)*sqrt(gm)*(fvr(:,i+1)-fvr(:,i))/dx
            enddo
        else
            do i=0,nx
                if(i.eq.itf) then
                    fcr(:,i) = fcr(:,i+1)
                    fvr(:,i) = fvr(:,i+1)
                endif
                rhs_wr(:,i) = -(fcr(:,i+1)-fcr(:,i))/dx + (muinf)*(fvr(:,i+1)-fvr(:,i))/dx
            enddo
        endif
    endif

end subroutine compute_rhs


!#####################################################################################
! SUBROUTINE COMPUTE_BCN
! Apply boundary conditions
!#####################################################################################

subroutine compute_bcn(wr)

    use modshocky
    implicit none

    integer :: i
    real(mykind), dimension(3,-ng:nx+ng) :: wr

    if(bcase.eq.0) then

        ! Supersonic inflow and outflow (transmissive BCs)
        do i = 1,ng
            wr(:,-i)   = wr(:,-i+1)
            wr(:,nx+i) = wr(:,nx+i-1)
        enddo

    elseif(bcase.eq.1) then

        ! Periodic boundary conditions
        wr(:,0) = wr(:,nx)
        do i=1,ng
            wr(:,nx+i) = wr(:,i)
            wr(:,-i)   = wr(:,nx-i)
        enddo

    elseif(bcase.eq.2) then

        ! Open - Right wall
        wr(1,nx) = wr(1,nx-1)
        wr(2,nx) = 0.0_mykind
        wr(3,nx) = ((((gm-1)*wr(3,nx-1))/wr(1,nx-1))/(gm-1))*wr(1,nx-1)

        do i = 1,ng
            wr(:,-i)   = wr(:,-i+1)
            wr(:,nx+i) = wr(:,nx+i-1)
        enddo

    elseif(bcase.eq.3) then

        ! Left wall - Open
        wr(1,0) = wr(1,1)
        wr(2,0) = 0.0_mykind
        wr(3,0) = ((((gm-1)*wr(3,1))/wr(1,1))/(gm-1))*wr(1,1)

        do i = 1,ng
            wr(:,-i)   = wr(:,-i+1)
            wr(:,nx+i) = wr(:,nx+i-1)
        enddo

    elseif(bcase.eq.4) then

        ! Sealed
        wr(1,0) = wr(1,1)
        wr(2,0) = 0.0_mykind
        wr(3,0) = ((((gm-1)*wr(3,1))/wr(1,1))/(gm-1))*wr(1,1)
        wr(1,nx) = wr(1,nx-1)
        wr(2,nx) = 0.0_mykind
        wr(3,nx) = ((((gm-1)*wr(3,nx-1))/wr(1,nx-1))/(gm-1))*wr(1,nx-1)

        do i = 1,ng
            wr(:,-i)   = wr(:,-i+1)
            wr(:,nx+i) = wr(:,nx+i-1)
        enddo

    endif

end subroutine compute_bcn


!#####################################################################################
! SUBROUTINE COMPUTE_RKT
! Advance in time vector w computing RHS
!#####################################################################################

subroutine compute_rkt(wr)

    use modshocky
    implicit none

    real(mykind), dimension(3,-ng:nx+ng) :: wr,wr_k,rhs_wr,rhs_wrk,ibm_wr
    real(mykind), dimension(4) :: alpha,beta

    real(mykind) :: int_dt,err_w
    integer :: i,ik,n_step

    wr_k = 0.0_mykind
    rhs_wr = 0.0_mykind
    rhs_wrk = 0.0_mykind
    ibm_wr = 0.0_mykind

    n_step = methodT
    ! Runge-Kutta coefficients
    if(n_step.eq.1) then

        alpha = [0.,0.,0.,0.]
        beta  = [1.0,0.,0.,0.]

    else if(n_step.eq.3) then

        alpha = [0.,(1.0/3.0),(2.0/3.0),0.]
        beta  = [(1.0/4.0),0.,(3.0/4.0),0.]

    else if(n_step.eq.4) then

        alpha = [0.,(1.0/2.0),(1.0/2.0),1.0]
        beta  = [(1.0/6.0),(2.0/6.0),(2.0/6.0),(1.0/6.0)]

    else

        write(*,*) 'Method not supported!'
        stop

    endif


    do ik = 1,n_step

        int_dt = dt*beta(ik)

        ! Update substep k field (alpha dependance)
        do i=0,nx
            wr_k(:,i) = wr(:,i) + alpha(ik)*int_dt*rhs_wrk(:,i)
        enddo
        call compute_bcn(wr)

        ! Interface track and emerging points correction
        call compute_itf(dt,wr)

        ! Ghost nodes corrections using BCs
        call compute_bcn(wr_k)

        ! Compute convective and viscous fluxes
        call compute_eul(wr_k,fc)
        call compute_vis(wr_k,fv)
        ! Compute IBM forcing term
        ! if(isibm.eq.1) call compute_ibm(wr,fc,ibm_wr)
        ! open(unit=1,file='testrkt.txt')
        ! do i = 0,nx
        !     write(1,*) i,ibm_wr(3,i)
        ! enddo
        ! close(1)

        ! Compute RHS
        call compute_rhs(fc,fv,ibm_wr,rhs_wr)

        ! Solution (time advancement) of NS equation at the current timestep
        do i=0,nx
             wr(:,i)=wr(:,i)+(rhs_wr(:,i)*dt)
        enddo
        ! Ghost nodes corrections using BCs
        call compute_bcn(wr)

        !     do i=0,nx
        !         wr_s(:,i)=wr_s(:,i)+beta(ik)*rhs_wrk(:,i)
        !         !wr_s(:,i)= rhs_wrk(:,i)
        !     enddo

        ! Update wall position after integrating
        xwall = wallpast(1) + uw*dt

    enddo

    ! do ik=1,n_step

    !     int_dt = dt*beta(ik) !dt*(alpha(ik)+beta(ik))

    !     ! Compute w_k for the substep k (alpha dependant) over the whole space domain at the current timestep
    !     do i=0,nx
    !         wr_k(:,i)=wr(:,i)+alpha(ik)*int_dt*rhs_wrk(:,i)
    !         !wr_k(:,i)=wr(:,i)
    !     enddo

    !     ! Interface track and emerging points correction
    !     call compute_itf(int_dt,wr_k)
    !     ! Ghost nodes corrections using BCs
    !     call compute_bcn(wr)

    !     ! Compute convective and viscous fluxes at the intermediate substep
    !     call compute_eul(wr_k,fc)
    !     call compute_vis(wr_k,fv)
    !     ! Compute IBM forcing term
    !     if(isibm.eq.1) call compute_ibm(wr_k,fc,ibm_wrk)
    !     ! Compute RHS of the k intermediate predictor substep
    !     call compute_rhs(fc,fv,ibm_wrk,rhs_wrk)
    !     ! Compute at the rk_step substep w_s value (beta dependant) at the current timestep
    !     do i=0,nx
    !         wr_s(:,i)=wr_s(:,i)+beta(ik)*rhs_wrk(:,i)
    !         !wr_s(:,i)= rhs_wrk(:,i)
    !     enddo

    !     ! Update wall position
    !     xwall = wallpast(1) + uw*int_dt      ! Internal timestep 'int_dt' from RK integration
    !     !write(*,*) 'Moves - END RK - Substep',ik

    ! enddo

    ! ! Solution (time advancement) of NS equation at the current timestep in X (corrector substep for RK)
    ! do i=0,nx
    !       wr(:,i)=wr(:,i)+(wr_s(:,i)*dt)
    ! enddo

    ! Ghost nodes corrections using BCs
    call compute_bcn(wr)


end subroutine compute_rkt


!#####################################################################################
! SUBROUTINE COMPUTE_DBG
! Advance in time vector w using Euler integration (debugging)
!#####################################################################################

subroutine compute_dbg(wr)

    use modshocky
    implicit none

    real(mykind), dimension(3,-ng:nx+ng) :: wr,rhs_wr,ibm_wr
    integer :: i

    ibm_wr = 0.0_mykind
    rhs_wr = 0.0_mykind

    ! Wall velocity sinusoidal modulation
    if(frqw.gt.(0.0_mykind)) then
        uw = ampw*(cos(2*pi*frqw*t))
        !open(unit=1, file='uwall.txt',access = 'append')
        !write(1,*) t,uw,xwall
        !close(1)
    endif

    ! Interface track and emerging points correction
    call compute_itf(dt,wr)
    ! Ghost nodes corrections using BCs
    call compute_bcn(wr)

    ! Compute convective and viscous fluxes
    call compute_eul(wr,fc)
    call compute_vis(wr,fv)
    ! Compute IBM forcing term
    if(isibm.eq.1) call compute_ibm(wr,fc,ibm_wr)
    ! Compute RHS
    call compute_rhs(fc,fv,ibm_wr,rhs_wr)
    ! Solution (time advancement) of NS equation at the current timestep
    do i=0,nx
          wr(:,i)=wr(:,i)+(rhs_wr(:,i)*dt)
    enddo
    ! Ghost nodes corrections using BCs
    call compute_bcn(wr)

    ! Update wall position after integrating
    xwall = wallpast(1) + uw*dt

end subroutine compute_dbg

!#####################################################################################
! SUBROUTINE COMPUTE_EXT
! Save analytical solution for Sod Shocktube problem
!#####################################################################################

subroutine compute_ext(tr)

    use modshocky
    implicit none

    integer :: i

    real(mykind) :: p1,p2,p3,p4,p5
    real(mykind) :: t1,t2,t3,t4,t5
    real(mykind) :: r1,r2,r3,r4,r5
    real(mykind) :: a1,a2,a3,a4,a5
    real(mykind) :: u1,u2,u3,u4,u5
    real(mykind) :: uh,ut,uc,us
    real(mykind) :: xh,xt,xc,xs

    real(mykind) :: p2p1,t2t1

    real(mykind) :: diff,tr,xx,xd

    p1 = pB
    p5 = pA
    t1 = p1/rhoB
    t5 = p5/rhoA
    a1 = sqrt((gm*p1)/rhoB)
    a5 = sqrt((gm*p5)/rhoA)

    u1 = 0.0_mykind
    u5 = 0.0_mykind

    diff = 1.0_mykind
    xd = 0.0_mykind
    p2 = 0.30273_mykind

    p2p1 = p2/p1
    us = a1*sqrt(((gm+1)/(2.0_mykind*gm)*(p2p1-1))+1)
    xs = xd + us*tr

    t2t1 = p2p1*((p2p1 + (gm+1)/(gm-1))/(1 + ((gm+1)/(gm-1))*p2p1))
    t2 = t2t1*t1
    u2 = ((a1**2)/(gm*us))*(p2p1 - 1)
    uc = u2
    u3 = u2
    xc = xd + uc*tr

    p3 = p2
    a3 = a5 - ((gm-1)/2.0_mykind)*uc
    ut = uc - a3
    xt = xd + ut*tr

    uh = -a5
    xh = xd + uh*tr

    ! Plot
    open(unit=13,file='exact.txt')
    do i = 0,nx
        xx = dx*(float(i))
        xx = xx - 5.0_mykind
        if(xx.le.xh) then
            r5 = (gm*p5)/(a5**2)
            t5 = p5/r5
            write(13,*) xx,p5,t5,r5,u5
        else if((xx.gt.xh).and.(xx.le.xt)) then
            ! Fan
            u4 = (2.0_mykind/(gm+1))*(((xx-xd)/tr) + a5)
            a4 = u4 - ((xx-xd)/tr)
            p4 = p5*(a4/a5)**((2.0_mykind*gm)/(gm-1))
            r4 = (gm*p4)/(a4**2)
            t4 = p4/r4
            write(13,*) xx,p4,t4,r4,u4
        else if((xx.gt.xt).and.(xx.le.xc)) then
            r3 = (gm*p3)/(a3**2)
            t3 = p3/r3
            write(13,*) xx,p3,t3,r3,u3
        else if((xx.gt.xc).and.(xx.le.xs)) then
            r2 = p2/t2
            write(13,*) xx,p2,t2,r2,u2
        else
            r1 = (gm*p1)/(a1**2)
            t1 = p1/r1
            write(13,*) xx,p1,t1,r1,u1
        endif
    enddo
    close(13)

    ! Compute p2
    !do while (diff.gt.1e-3)

        !us = a1*sqrt(((gm+1)/(2.0_mykind*gm)*((p2/p1)-1))+1)
        !p1p5 = (p2/p1)*(1-((a1/a5)*((gm-1)/(2*gm)))*(((p2/p1)-1)/(us/a1)))**((-2.0_mykind*gm)/(gm-1))
        !diff = abs(p1p5 - (p1/p5))
        !p2 = p2*1.001

    !enddo


end subroutine compute_ext


!#####################################################################################
! SUBROUTINE OUTPUT_DUMP
! Save outputs from the current timestep
!#####################################################################################

subroutine output_dump(wr)

    use modshocky
    implicit none

    integer :: i,dim1,dim2,qi,qk
    real(mykind) :: wr_err,wr_real,u,q,p,h
    real(mykind), dimension(3,-ng:nx+ng) :: wr

    ! ONE-DIMENSIONAL OUTPUT
    open(unit=1,file='output.txt')
    open(unit=2,file='outflow.txt')

    wr_real = 0.0_mykind
    wr_err = 0.0_mykind

    do i=-ng,nx+ng

        ! Flow quantities
        u = wr(2,i)/wr(1,i)
        q = (u**2)*0.5_mykind
        p = (gm-1)*(wr(3,i)-wr(1,i)*q)
        h = p/wr(1,i)

        write(2,*) i,x(i),u,q,p,h

        if(xcase.eq.3) then
            !wr_real = 1.0_mykind+0.2_mykind*sin(pi*(x(i)-ubrgs*t))
            if(x(i).le.(xdf)) then
                wr_real = rhoa
            else
                wr_real = rhob
            endif
            wr_err  = abs(wr(1,i)-wr_real)
            write(1,*) i,x(i),wr(1,i),wr(2,i),wr(3,i),wr_real,wr_err
        elseif(xcase.eq.2) then
            wr_real = 1.0_mykind+0.2_mykind*sin(pi*(x(i)-uinf*t))
            wr_err  = abs(wr(1,i)-wr_real)
            write(1,*) i,x(i),wr(1,i),wr(2,i),wr(3,i),wr_real,wr_err
        elseif(xcase.eq.1) then
            write(1,*) i,x(i),wr(1,i),wr(2,i),wr(3,i)
        else
            write(1,*) i,x(i),wr(1,i),wr(2,i),wr(3,i)
        endif

    enddo

    close(1)
    close(2)

    ! PLOT3D OUTPUT
    ! Generation of the .xyz file
    dim1 = size(wr,2)
    dim2 = nz+1
    open(unit=3, file='OUT/grid.xyz',form="unformatted")
    write(3) dim1,dim2
    write(3) ((x(qi),qi=-ng,nx+ng),qk=0,nz),((z(qk),qi=-ng,nx+ng),qk=0,nz)
    close(3)

end subroutine output_dump


!#####################################################################################
! SUBROUTINE OUTPUT_PLOT3D 
! Save outputs from the current timestep to PLOT3D slice
!#####################################################################################

subroutine output_plot3d(wr)

    use modshocky
    implicit none

    integer :: i,dim1,dim2,qm,qi,qk
    real(mykind) :: u,q
    real(mykind), dimension(3,-ng:nx+ng) :: wr
    real(mykind), dimension(-ng:nx+ng) :: p,h

    write(iter,'(I7.7)') it

    do i = -ng,nx+ng
        u = wr(2,i)/wr(1,i)
        q = (u**2)*0.5_mykind
        p(i) = (gm-1)*(wr(3,i)-wr(1,i)*q)
        h(i) = p(i)/wr(1,i)
    enddo

    ! PLOT3D OUTPUT
    ! Solution file allocation
    do qk=0,nz
        do qi=-ng,nx+ng
            qslice(1,qi,qk) = wr(1,qi)
            qslice(2,qi,qk) = wr(2,qi)
            qslice(3,qi,qk) = wr(3,qi)
            qslice(4,qi,qk) = p(qi)
            qslice(5,qi,qk) = h(qi)
        enddo
    enddo

    ! Plot3D output
    dim1 = size(wr,2)
    dim2 = nz+1
    open(unit=3, file='OUT/field_'//iter//'.q',form="unformatted")
    write(3) dim1,dim2
    write(3) 2.0,0.0,1000000.0,35.0
    write(3) (((qslice(qm,qi,qk),qi=-ng,nx+1),qk=0,nz),qm=1,5)
    close(3)

end subroutine output_plot3d


!#####################################################################################
! SUBROUTINE READINPUT
! Read input.txt contents and assign simulation parameters
!#####################################################################################

subroutine readinput()

    use modshocky
    implicit none

    ! Read settings file
    open(unit = dataID, file='input.txt',iostat = err)
    if(err .ne. 0) stop ' Error opening data file! '
    ! Skip header
    read(dataID,*)                !1
    read(dataID,*)                !2
    read(dataID,*)                !3
    read(dataID,*)                !4
    read(dataID,*)                !5
    ! Read solver main setting
    read(dataID,*) icase          !6
    read(dataID,*) xcase          !7
    read(dataID,*)                !8
    read(dataID,*)                !9
    ! Read simulation parameters - inflow
    read(dataID,*) Ma             !10
    read(dataID,*) Re             !11
    read(dataID,*) rhoinf         !12
    read(dataID,*) pinf           !13
    read(dataID,*) muinf          !14
    read(dataID,*) Pr             !15
    read(dataID,*) Gm             !16
    read(dataID,*)                !17
    read(dataID,*)                !18
    ! Read simulation parameters - shocktube
    read(dataID,*) rhoA           !19
    read(dataID,*) pA             !20
    read(dataID,*) uA             !21
    read(dataID,*) rhoB           !22
    read(dataID,*) pB             !23
    read(dataID,*) uB             !24
    read(dataID,*) xdf            !25
    ! Read discretization settings
    read(dataID,*)                !26
    read(dataID,*)                !27
    read(dataID,*) l_tot          !28
    read(dataID,*) nx             !29
    read(dataID,*) itmax          !30
    read(dataID,*) dt             !31
    read(dataID,*) methodT        !32
    read(dataID,*) corder         !33
    read(dataID,*) dorder         !34
    read(dataID,*) iweno          !35
    read(dataID,*) iflux          !36
    read(dataID,*)                !37
    ! Read BCs settings
    read(dataID,*)                !38
    read(dataID,*) bcase          !39
    read(dataID,*) isibm          !40
    read(dataID,*) xwall          !41
    read(dataID,*) uw             !42
    read(dataID,*) frqw           !43
    read(dataID,*)                !44
    read(dataID,*)                !45
    ! Read Burgers settings
    read(dataID,*) ubrgs          !46
    read(dataID,*)                !47
    read(dataID,*)                !48
    read(dataID,*) idbg           !49
    close(dataID)


end subroutine readinput


!#################################################################################################################
! SUBROUTINE DEBUG
! Show all variables saved
!#################################################################################################################

subroutine debug()

    use modshocky
    implicit none

     write(*,*)
     write(*,*) '#########################################'
     write(*,*)
     write(*,*) 'Selected simulation parameters'
     if(xcase.eq.0) then
          write(*,*) 'Inlow:'
          write(*,*) 'Ma    = ',ma
          write(*,*) 'Re    = ',re
          write(*,*) 'Rho   = ',rhoinf
          write(*,*) 'p     = ',pinf
     elseif(xcase.eq.1) then
          write(*,*) 'Shocktube:'
          write(*,*) 'RhoA  = ',rhoA
          write(*,*) 'pA    = ',pA
          write(*,*) 'UA    = ',uA
          write(*,*) 'RhoB  = ',rhoB
          write(*,*) 'pB    = ',pB
          write(*,*) 'UB    = ',uB
          write(*,*) 'xDF   = ',xdf
     elseif(xcase.eq.2) then
          write(*,*) 'Periodic signal test'
     elseif(xcase.eq.3) then
          write(*,*) 'Burgers equation debug test'
          write(*,*) 'CC    = ',ubrgs
     else
          write(*,*) 'Walltube - IBM:'
          write(*,*) 'RhoA  = ',rhoA
          write(*,*) 'pA    = ',pA
          write(*,*) 'UA    = ',uA
          write(*,*) 'RhoB  = ',rhoB
          write(*,*) 'pB    = ',pB
          write(*,*) 'UB    = ',uB
          write(*,*) 'xWall   = ',xwall
          write(*,*) 'uWall   = ',uw
          if(isibm.eq.0) then
            write(*,*) 'IBM is off!'
          endif
     endif
     write(*,*) 'Flow constants:'
     write(*,*) 'Pr    = ',pr
     write(*,*) 'Gm    = ',gm
     write(*,*) 'muInf = ',muinf
     write(*,*) 'Space and time discretization:'
     write(*,*) 'Lx    = ',l_tot
     write(*,*) 'Nx    = ',nx
     write(*,*) 'Ng    = ',ng
     write(*,*) 'Dt    = ',dt
     write(*,*)
     write(*,*) 'Simulation settings:'
     if(icase.eq.0) then
          write(*,*) 'Solv# = ',icase, '(EULER EQUATIONS)'
     else
          write(*,*) 'Solv# = ',icase, '(NAVIER-STOKES)'
     endif
     write(*,*) 'Itmax = ',itmax
     if(bcase.eq.0) then
          write(*,*) 'BCnd# = ',bcase, '(SUPERSONIC)'
     elseif(bcase.eq.1) then
          write(*,*) 'BCnd# = ',bcase, '(PERIODICAL)'
     elseif(bcase.eq.2) then
          write(*,*) 'BCnd# = ',bcase, '(SEALED OUTFLOW)'
     endif
     write(*,*) 'Conv# = ',corder
     write(*,*) 'Diff# = ',dorder
     write(*,*) 'Rktt# = ',methodT
     write(*,*) 'WENO# = ',iweno
     if(iflux.eq.0) then
          write(*,*) 'Flux# = ',iflux, '(COMPONENT)'
     else
          write(*,*) 'Flux# = ',iflux, '(CHARACTERISTIC)'
     endif
     if(isibm.eq.1) then
        write(*,*) 'Immersed-Boundary is ON'
        write(*,*)
        write(*,*) 'Interface debug...'
        write(*,*)
        write(*,*) 'xwall = ',xwall
        write(*,*) 'uwall = ',uw
        write(*,*) 'Frequency of oscilation is ',frqw
        write(*,*) 'Closest left node:',x(itf),dxl
        write(*,*) 'Closest right node:',x(itf+1),dxr
     endif
     write(*,*)
     write(*,*)
     write(*,*) '#########################################'
 
end subroutine debug


!#################################################################################################################
! SUBROUTINE STARTUP
! Allocate variables
!#################################################################################################################

subroutine startup()

    use modshocky
    implicit none

    integer :: ngi

    ngi = max(corder,dorder)/2

    ! Allocate flow state
    allocate(w(3,-ngi:nx+ngi),rhs_w(3,-ngi:nx+ngi))
    w = 0.0_mykind
    rhs_w = 0.0_mykind

    fsize = size(w,2)-2

    ! Allocate fluxes
    allocate(fc(3,0:fsize),fv(3,0:fsize))
    fc = 0.0_mykind
    fv = 0.0_mykind

    ! Allocate domain
    allocate(x(-ngi:nx+ngi))
    x = 0.0_mykind

    ! Allocate PLOT3D domain and solution file
    allocate(z(0:nz))
    allocate(qslice(5,-ngi:nx+ngi,-ng:nz+ngi))
    z = 0.0_mykind
    qslice = 0.0_mykind

end subroutine startup


!#################################################################################################################
! SUBROUTINE SHUTDOWN
! Deallocate variables
!#################################################################################################################

subroutine shutdown()

     use modshocky
     implicit none

     deallocate(w)
     deallocate(rhs_w)
     deallocate(fc)
     deallocate(fv)
     deallocate(x)
     deallocate(z)
     deallocate(qslice)

end subroutine shutdown


!#################################################################################################################
!#################################################################################################################

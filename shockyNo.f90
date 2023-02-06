!#####################################################################################
! PARAMETERS AND CONSTANTS
!#####################################################################################

module modshocky

     implicit none
     save
     integer, parameter :: singtype = selected_real_kind(6,37)    ! Single precision
     integer, parameter :: doubtype = selected_real_kind(15,307)  ! Double precision
     integer, parameter :: mykind   = doubtype

     integer :: nx,ng,fsize,itmax,it = 0                                ! Simulation discretization
     integer :: corder,dorder,methodT,iweno,iflux,itout = 100           ! Discretization settings
     integer :: dataID = 10, err = 0, icase, bcase, xcase, idbg         ! Simulation settings

     character(7) :: iter                                     ! Internal timetick

     real(mykind) :: ma,re,pr,gm,muinf                        ! Non-dimensional groups
     real(mykind) :: u0,rho0,p0,t0,l0,tme0                    ! Reference values
     real(mykind) :: uinf,pinf,rhoinf,rhouinf,rhoeinf,ubrgs   ! Inflow conditions
     real(mykind) :: uA,pA,rhoA,rhoUA,rhoEA                   ! Shocktube - init left side
     real(mykind) :: uB,pB,rhoB,rhoUB,rhoEB                   ! Shocktube - init right side

     real(mykind) :: l_tot,t_tot,t=0,dt,dx,xdf                ! Time and space domain details
     real(mykind) :: pi=acos(-1.)                             ! Math constants

     real(mykind), dimension(:,:), allocatable :: w,rhs_w,fc,fv  ! State vector and RHS, fluxes

end module modshocky


!#####################################################################################
! MAIN PROGRAM - v1.1 STABLE 
!#####################################################################################

program shockyNo

    use modshocky
    implicit none

    write(*,*)
    write(*,*) 'ShockyNo v1.1'
    write(*,*)
    write(*,*) 'Reading input for the simulation ...'
    call readinput()
    write(*,*) 'Starting up ...'
    call startup()
    write(*,*) 'Initializing the flow domain ...'
    call init()

    ! DEBUG ROUTINE
    if(idbg.eq.1) then
        call debug()
    endif

    write(*,*)
    write(*,*) 'Simulation start ...'
    write(*,*)

    ! Runge-Kutta order 3/4 integration method
    write(*,*) '... Running Runge-Kutta'
    do while(it.le.itmax)

        if (mod(it,itout).eq.0) then
            ! Output video 
            write(*,*)  it,t,w(1,nx/2),w(2,nx/2),w(3,nx/2)
        endif

        if(isnan(w(1,nx/2))) then
            write(*,*)
            write(*,*) 'Failed, stopping the simulation ...'
            write(*,*)
            exit
        endif

        call  compute_rkt(w)

        it=it+1
        t=t+dt

    enddo
    ! Save output at the final timestep
    write(*,*) 'Generating outputs ...'
    call output_dump(w)

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
! Initialize the space and time domain
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
    dx=l_tot/float(nx)
    ng = max(corder,dorder)/2
    dt = dt/tme0
    t_tot = dt*itmax

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
        write(11,*) i,dx*i,w(1,i),w(2,i),w(3,i)
        enddo
        close(11)

    elseif(xcase.eq.1)  then

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

             ! Shu-Osher shocktube
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
             write(11,*) i,dx*i-(l_tot*0.5_mykind),w(1,i),w(2,i),w(3,i)

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
            write(11,*) i,dx*(float(i)),w(1,i),w(2,i),w(3,i)

        enddo
        close(11)

    else

        ! Burgers advection/Scalar discontinuity
        open(unit=11,file='init.txt')
        do i = -ng,nx+ng

            !w(1,i) = 1.0_mykind+0.2_mykind*sin(pi*dx*(float(i)))    ! Phi
            !w(2,i) = 0.0_mykind
            !w(3,i) = 0.0_mykind
            !write(11,*) i,dx*(float(i)),w(1,i),w(2,i),w(3,i)

            if(i.le.(xdf*nx)) then
                w(1,i) = rhoa
            else
                w(1,i) = rhob
            endif
            w(2,i) = 0.0_mykind
            w(3,i) = 0.0_mykind

            write(11,*) i,dx*i,w(1,i),w(2,i),w(3,i)

        enddo
        close(11)

    endif

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
! SUBROUTINE COMPUTE_RHS
! Current fluxes as input and rhs_f as output
!#####################################################################################

subroutine compute_rhs(fcr,fvr,rhs_wr)

    use modshocky
    implicit none

    integer :: i
    real(mykind), dimension(3,0:fsize)    :: fcr,fvr
    real(mykind), dimension(3,-ng:nx+ng)  :: rhs_wr

    rhs_wr = 0.0_mykind

    if(icase.eq.0) then
        do i=0,nx
            rhs_wr(:,i) = -(fcr(:,i+1)-fcr(:,i))/dx
        enddo
    else
        if(xcase.eq.0) then
            do i=0,nx
                rhs_wr(:,i) = -(fcr(:,i+1)-fcr(:,i))/dx + (ma/re)*sqrt(gm)*(fvr(:,i+1)-fvr(:,i))/dx
            enddo
        else
            do i=0,nx
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

    if(bcase.eq.1) then
        ! Periodic boundary conditions
        wr(:,0) = wr(:,nx)
        do i=1,ng
            wr(:,nx+i) = wr(:,i)
            wr(:,-i)   = wr(:,nx-i)
        enddo
    else

        ! Supersonic inflow and outflow (transmissive BCs)
        do i = 1,ng
            wr(:,-i)  = wr(:,-i+1)
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

    real(mykind), dimension(3,-ng:nx+ng) :: wr,wr_k,rhs_wrk,wr_s
    real(mykind), dimension(4) :: alpha,beta

    integer :: i,ik,n_step

    rhs_wrk = 0.0_mykind
    wr_s = 0.0_mykind
    wr_k = 0.0_mykind

    n_step = methodT
    ! Runge-Kutta coefficients
    if(n_step.eq.3) THEN

        alpha = [0.,(1.0/3.0),(2.0/3.0),0.]
        beta  = [(1.0/4.0),0.,(3.0/4.0),0.]

    else if(n_step.eq.4) THEN

        alpha = [0.,(1.0/2.0),(1.0/2.0),1.0]
        beta  = [(1.0/6.0),(2.0/6.0),(2.0/6.0),(1.0/6.0)]

    else
        write(*,*) 'Method not supported!'
        stop
    endif

    do ik=1,n_step

        ! Compute w_k for the substep k (alpha dependant) over the whole space domain at the current timestep
        do i=0,nx
            wr_k(:,i)=wr(:,i)+alpha(ik)*dt*rhs_wrk(:,i)
        enddo

        ! Ghost nodes corrections using BCs
        call compute_bcn(wr_k)

        ! Compute RHS of the k intermediate predictor substep 
        call compute_eul(wr_k,fc)
        call compute_vis(wr_k,fv)
        call compute_rhs(fc,fv,rhs_wrk)
        ! Compute at the rk_step substep w_s value (beta dependant) at the current timestep
        do i=0,nx
            wr_s(:,i)=wr_s(:,i)+beta(ik)*rhs_wrk(:,i)
        enddo

    enddo

    ! Solution (time advancement) of NS equation at the current timestep in X (corrector substep for RK)
    do i=0,nx
          wr(:,i)=wr(:,i)+(wr_s(:,i)*dt)
    enddo

    ! Ghost nodes corrections using BCs
    call compute_bcn(wr)

end subroutine compute_rkt


!#####################################################################################
! SUBROUTINE OUTPUT_DUMP
! Save outputs from the current timestep
!#####################################################################################

subroutine output_dump(wr)

    use modshocky
    implicit none

    integer :: i
    real(mykind) :: xi,wr_err,wr_real,u,q,p,h
    real(mykind), dimension(3,-ng:nx+ng) :: wr

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

        write(2,*) i,xi-5.0_mykind,u,q,p,h

        xi=dx*(float(i))
        if(xcase.eq.3) then
            !wr_real = 1.0_mykind+0.2_mykind*sin(pi*(xi-ubrgs*t))
            xi = xi-ubrgs*t
            if(xi.le.(xdf)) then
                wr_real = rhoa
            else
                wr_real = rhob
            endif
            xi=dx*(float(i))
            wr_err  = abs(wr(1,i)-wr_real)
            write(1,*) i,xi,wr(1,i),wr(2,i),wr(3,i),wr_real,wr_err
        elseif(xcase.eq.2) then
            wr_real = 1.0_mykind+0.2_mykind*sin(pi*(xi-uinf*t))
            wr_err  = abs(wr(1,i)-wr_real)
            write(1,*) i,xi,wr(1,i),wr(2,i),wr(3,i),wr_real,wr_err
        elseif(xcase.eq.1) then
            write(1,*) i,xi-(l_tot*0.5_mykind),wr(1,i),wr(2,i),wr(3,i)
        else
            write(1,*) i,xi,wr(1,i),wr(2,i),wr(3,i)
        endif

    enddo

    close(1)
    close(2)

end subroutine output_dump


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
    read(dataID,*)                !40
    read(dataID,*)                !41
    ! Read Burgers settings
    read(dataID,*) ubrgs          !42
    read(dataID,*)                !43
    read(dataID,*)                !44
    read(dataID,*) idbg           !45
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
     else
          write(*,*) 'Burgers equation debug test'
          write(*,*) 'CC    = ',ubrgs
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
     else
          write(*,*) 'BCnd# = ',bcase, '(PERIODICAL)'
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

end subroutine shutdown


!#################################################################################################################
!#################################################################################################################

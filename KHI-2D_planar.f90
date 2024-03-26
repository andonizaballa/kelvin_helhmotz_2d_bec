!*************************************************************************
!*** Gross-Pitaevskyy equation in 2D: ground state and dynamics
!*** Vortex creation by barrier removal in a counterflow
!***
!*** Michele Modugno & Chiara Fort, 2023 - 2024
!*************************************************************************

!*************************************************************************
! Our aim is now to be able to modify this code with the adimensional variables
! we defined in the theory part. 
!*************************************************************************
program main
  use,intrinsic :: iso_c_binding
  implicit none
  !!integer, parameter :: n1=128,n2=64
  integer, parameter :: n1=256,n2=128
   integer, dimension(1:2) :: dim=[n1,n2]
  ! physical constants
  complex(8), parameter :: ci=(0.d0,1.d0)
  real(8), parameter :: pi=3.141592653589793238d0
  real(8), parameter :: hbar=1.054572d-34            ! J*s
  real(8), parameter :: amu=1.66053873d-27           ! Kg (atomic mass unit)
  real(8), parameter :: a_bohr=0.5291772083d-10      ! m
  real(8), parameter :: m=12.d0*amu                  ! Kg
  real(8), parameter :: l0=1.d-6                     ! to express coordinates in micron
  real(8), parameter :: e0=hbar**2.d0/(m*l0**2.d0)   ! energy scale (ho correspondence)
  real(8), parameter :: ascat = 1010.0d0*a_bohr/l0 
  real(8), parameter :: nat = 1.d3 
  real(8), parameter :: u = 4.d0*pi*nat*ascat
  !grid
  real(8), parameter, dimension(2) :: xmax=[40.d0,20.d0]
  real(8), parameter, dimension(2) :: xmin=-xmax
  ! waveguide
  real(8), parameter :: A=2*pi*4100*hbar/e0         ! barrier amplitude
  real(8), parameter :: sigma=1.0                   ! barrier width
  real(8), parameter :: B=2.d0*pi*5000*hbar/e0      ! external barrier for the waveguide
  real(8), parameter :: w=1.d0                      ! waveguide stiffness
  real(8), parameter :: r2=xmax(2)-2.d0             ! transverse width of the waveguide
  ! vortices
  real(8), parameter :: Nv=10                       ! number of vortices to be created
  real(8), parameter :: k0=Nv/2.d0*pi/xmax(1)       ! phase factor correspending to a velocity v0=hbar*k0/m
  ! ground state
  real(8), parameter :: target_tol=1.d-8
  ! evolution
  real(8), parameter :: t_evol = 100.0d0            ! ms total evolution time
  real(8), parameter :: dt=0.05d-3                  ! ms (0.05d-3)
  integer, parameter :: shots=nint(t_evol)          ! saves every ms
  integer, parameter :: kmax=t_evol/dt
  real(8), parameter :: dwt=dt*1.d-3*e0/hbar
  real(8), parameter :: t0 = 25.d0                  ! switch-off time of teh barrier
  real(8) :: t,factor
  !
  include 'fftw3.f03'
  integer(8) :: plan_f,plan_b
  integer :: iret,nthreads=1 !6
  integer, parameter :: forth=-1,back=-forth

  ! ----
  complex(8), dimension(n1,n2) :: psi
  complex(8), dimension(n1,n2) :: in,out
  real(8), dimension(n1,n2) :: uext,uWAVEGUIDE,uBARRIER
  real(8), dimension(n1,n2) :: psq
  real(8), dimension(n1) :: x1,p1
  real(8), dimension(n2) :: x2,p2
  real(8), dimension(2) :: dx
  real(8) :: mu,etot
  integer :: i1,i2,iter,icount=0
  character(3) :: file
  !

  call dfftw_init_threads(iret)
  call fftw_plan_with_nthreads(nthreads)

  write(*,*) nthreads, iret

  call dfftw_plan_dft_2d(plan_f,n1,n2,in,out,-1,FFTW_MEASURE)
  call dfftw_plan_dft_2d(plan_b,n1,n2,in,out,1,FFTW_MEASURE)

  call system('mkdir data')

  ! Create the initial state with the potential 
  call init_stat

  ! write the potential in a .dat file
  call wrt_pot

  ! trial wave function
  forall(i1=1:n1,i2=1:n2) psi(i1,i2) = exp(-1/4.d0**2*(x1(i1)**2/10.d0**2 + x2(i2)**2/5.d0**2))

  ! normalize the wave function
  call normalize(psi(:,:))
!!$  call wrt_dat('in-')

  ! find the ground state & write data
  ! call find_the_ground
  call grad_conjg

  ! add a phase for creating the counterflow
  call adding_phase

  call wrt_dat('000')


  ! real time dynamics

  open(UNIT=10,FILE='data/ener.dat',STATUS='unknown')
  open(UNIT=16,FILE='data/xmean.dat',STATUS='unknown')
  open(UNIT=41,FILE='data/ft_x.dat',STATUS='unknown')
  open(UNIT=42,FILE='data/ft_y.dat',STATUS='unknown')

  t=0.d0
  do iter = 0,kmax
     t = dt*dble(iter)

     if(t.lt.t0) then
        factor = t/t0
     else
        factor = 1.d0
     end if

     uext(:,:) = uWAVEGUIDE(:,:) + (1.d0-factor)*uBARRIER(:,:)

     ! evolution
     call one_step_evolution

     if(mod(iter,kmax/1000).eq.0.) then
        call ener(etot)
        write(*,'(2x,f8.4,2x,f10.4)') t,etot
     end if

     if(mod(iter,kmax/shots).eq.0.) then
        write(file,'(i3.3)') icount
        call wrt_dat(file)
        icount = icount+1
     end if

  end do

  close(10)
  close(16)
  close(41)
  close(42)

  call dfftw_destroy_plan(plan_f)
  call dfftw_destroy_plan(plan_b)
  call dfftw_cleanup_threads()

contains

  subroutine init_stat
    implicit none
    real(8), dimension(2) :: dp
    
   ! Initialize the grid
    dx(:)=(xmax(:)-xmin(:))/dim(:)
    dp(:)=2.d0*pi/(xmax(:)-xmin(:))

   ! Define the values for x and p in the two directions
    forall(i1=1:n1) x1(i1) = xmin(1) + (i1-1)*dx(1)
    forall(i1=1:n1/2) p1(i1) = dp(1)*(i1-1)
    forall(i1=n1/2+1:n1) p1(i1) = dp(1)*(i1-1-n1)

    forall(i2=1:n2) x2(i2) = xmin(2) + (i2-1)*dx(2)
    forall(i2=1:n2/2) p2(i2) = dp(2)*(i2-1)
    forall(i2=n2/2+1:n2) p2(i2) = dp(2)*(i2-1-n2)

    forall(i1=1:n1,i2=1:n2) psq(i1,i2) = p1(i1)**2 + p2(i2)**2

    ! waveguide
    forall(i1=1:n1,i2=1:n2) uWAVEGUIDE(i1,i2) = B*(Tanh((x2(i2)-r2)/w)+1) + B*(Tanh((-r2-x2(i2))/w)+1)

    ! barrier
    forall(i1=1:n1,i2=1:n2) uBARRIER(i1,i2) = A*EXP(-(x2(i2))**2/(2.d0*sigma**2))

    ! total potential
    uext(:,:) = uWAVEGUIDE(:,:) + (1.d0 - factor)*uBARRIER(:,:)

  end subroutine init_stat


  subroutine grad_conjg
    implicit none
    complex(8), dimension(:,:), allocatable :: dpsi,gnr,hnr
    real(8) :: ener,norm,tmp,tol_gp
    ! Numerical Recipes
    real(8) :: gg,dgg,dgg_fr
    integer :: niter
    logical :: fin

    allocate(dpsi(n1,n2),gnr(n1,n2),hnr(n1,n2))

    call gradient(psi,dpsi,u,ener,mu)

    gnr=-dpsi
    hnr=-dpsi
    dpsi=-dpsi
    ! debut de l'iteration
    niter=0
    fin=.false.
    do while(.not.fin)

       niter=niter+1

       call linmin(psi,dpsi)
       call gradient(psi,dpsi,u,ener,mu)

       ! renormalization of the wave function (does not change the energy)
       norm=sqrt(product(dx)*sum(abs(psi)**2))
       psi=psi/norm

       gg=product(dx)*sum(abs(gnr)**2)
       dgg_fr=product(dx)*sum(abs(dpsi)**2)  ! Fletcher-Reeves
       dgg=product(dx)*sum(real(conjg(dpsi)*(dpsi+gnr))) ! Polak-Ribiere
       fin=(((sqrt(dgg_fr).lt.target_tol*abs(mu)/norm).or.(gg.eq.0.d0)).and.niter.gt.5)
       tol_gp = sqrt(dgg_fr)/abs(mu)*norm

       if(niter.eq.2) then
          tmp = tol_gp
       end if

       ! keeps track of the path to the minimum
       write(*,'(i5,4(2x,g16.8))') niter,ener,mu,tol_gp

!!$       if(mod(niter,50).eq.0.) then
!!$          call wrt_dat('xxx')
!!$       end if

       if(.not.fin)then
          gnr=-dpsi
          hnr=gnr+(dgg/gg)*hnr
          dpsi=hnr
       else
          tol_gp = tmp
       endif

    enddo

    deallocate(dpsi,gnr,hnr)

    write(*,*)

  end subroutine grad_conjg


  subroutine gradient(psi,dpsi,u,ener,mu)
    ! calcule le gradient : dpsi = h_0 psi/||psi||^2 + u * |psi|^2 psi/||psi||^4  -mu * psi/||psi||^2
    implicit none
    complex(8), dimension(n1,n2) :: psi,dpsi
    real(8) :: u,ener,mu
    real(8) :: norme2,fac2,fac4

    ! norme au carre de psi
    norme2=product(dx)*sum(abs(psi)**2)
    fac2=1.d0/norme2
    fac4=fac2**2

    ! kinetic energy
    in = psi(:,:)
    call fft_transform(forth)
    in = fac2*0.5d0*psq*out
    call fft_transform(back)
    dpsi = out

    ! potential + mean-field
    dpsi = dpsi + (fac2*uext + fac4*u*abs(psi)**2)*psi

    ! chemical potential
    mu = product(dx)*sum(real(conjg(psi)*dpsi))
    dpsi = dpsi-(mu*fac2)*psi

    ! total energy
    ener = mu - fac4*0.5d0*u*product(dx)*sum(abs(psi)**4)

  end subroutine gradient


  subroutine linmin(psi0,dpsi)
    ! bouge psi selon la ligne psi0 + lambda * dpsi de facon a minimiser l'energie
    ! retourne le resultat dans psi0, ne modifie pas dpsi (contrairement a Numerical Recipes)
    implicit none
    integer ::i
    integer :: segno,i_min,i_max
    complex(8), dimension(:,:), allocatable :: grad
    complex(8), dimension(n1,n2) :: psi0,dpsi
    real(8) :: u0,ener,mu
    integer, parameter :: nparam=11
    real(8), dimension(nparam) :: param
    integer, parameter :: degree=6
    real(8), dimension(degree-1) :: q,r
    real(8), dimension(degree+1) :: p
    complex(8), dimension(degree) :: roots
    real(8) :: root
    real(8), external :: ederiv
    logical :: found=.false.
    ! Numerical Recipes¡
    real(8) :: lambda_min=0.d0
    real(8), parameter :: tol=1.d-7
    logical, parameter :: polish=.true.

    ! calcul des trois premiers parametres necessaires a la connaissance de E
    ! sur la droite  psi0 + lambda * dpsi
    param(1) = product(dx)*sum(abs(psi0)**2)
    param(2)  =product(dx)*sum(conjg(psi0)*dpsi + conjg(dpsi)*psi0)
    param(3) = product(dx)*sum(abs(dpsi)**2)
    ! calcul des trois parametres suivants
    ! calcul du gradient de l'energie au point psi0 pour le cas sans interaction
    u0=0.d0
    allocate(grad(n1,n2))
    call gradient(psi0,grad,u0,ener,mu)
    ! compenser le terme de potentiel chimique et la normalisation
    ! si bien que grad = h_0 psi_0
    grad = param(1)*grad + mu*psi0
    param(4) = product(dx)*sum(conjg(psi0)*grad)
    param(5) = product(dx)*sum(conjg(dpsi)*grad + conjg(grad)*dpsi)
    ! calcul de l'energie au point dpsi pour le cas sans interaction
    call gradient(dpsi,grad,u0,ener,mu)   ! grad est perdu
    param(6) = param(3)*ener
    deallocate(grad)

    ! les cinq derniers parametres
    param(7) = 0.5d0*u*product(dx)*sum(abs(psi0)**4)
    param(8) = u*product(dx)*sum(abs(psi0)**2*(conjg(psi0)*dpsi + conjg(dpsi)*psi0))
    param(9) = 0.5d0*u*product(dx)*sum((conjg(psi0)*dpsi+conjg(dpsi)*psi0)**2 + 2.d0*abs(psi0*dpsi)**2)
    param(10) = u*product(dx)*sum(abs(dpsi)**2*(conjg(psi0)*dpsi + conjg(dpsi)*psi0))
    param(11) = 0.5d0*u*product(dx)*sum(abs(dpsi)**4)

    ! parametrs for dE/dlambda

    q(1) = param(1)*param(4)+param(7)
    q(2) = param(2)*param(4)+param(1)*param(5)+param(8)
    q(3) = param(3)*param(4)+param(2)*param(5)+param(1)*param(6)+param(9)
    q(4) = param(3)*param(5)+param(2)*param(6)+param(10)
    q(5) = param(3)*param(6)+param(11)

    r(1) = param(1)**2
    r(2) = 2.*param(1)*param(2)
    r(3) = param(2)**2+2.*param(1)*param(3)
    r(4) = 2.*param(2)*param(3)
    r(5) = param(3)**2

    p(1) = q(2)*r(1)-q(1)*r(2)
    p(2) = 2.*(q(3)*r(1)-q(1)*r(3))
    p(3) = (q(3)*r(2)-q(2)*r(3))+3.*(q(4)*r(1)-q(1)*r(4))
    p(4) = 2.*(q(4)*r(2)-q(2)*r(4))+4.*(q(5)*r(1)-q(1)*r(5))
    p(5) = 3.*(q(5)*r(2)-q(2)*r(5))+(q(4)*r(3)-q(3)*r(4))
    p(6) = 2.*(q(5)*r(3)-q(3)*r(5))
    p(7) = q(5)*r(4)-q(4)*r(5)

    ! finds the roots (stationary points of E(lambda))
    call zroots(p,degree,roots,polish)

    if(ederiv(0.d0,p,degree+1).gt.0.d0)then
       segno=-1
       i_min=degree
       i_max=1
    else
       segno=+1
       i_min=1
       i_max=degree
    endif

    found = .false.
    lambda_min = 0.d0
    do i = i_min, i_max,segno
       if((abs(imag(roots(i))).eq.0.).and.(segno*real(roots(i)).ge.0.)) then
          root = real(roots(i))
          if( (ederiv(root*(1-tol),p,degree+1).lt.0.).and.(ederiv(root*(1+tol),p,degree+1).gt.0.) ) then
             write(2,*) 'found',root,i
             found = .true.
             lambda_min = root
             exit
          else
             write(2,*) 'not minimum',root,i
          end if
       end if
    end do

    if(.not.found) then
       write(2,*) 'not found',i
       write(3,*) p(:)
    end if

    ! amener psi sur le minimum
    psi0 = psi0 + lambda_min*dpsi

  end subroutine linmin


  subroutine one_step_evolution
    implicit none
    real(8), dimension(:,:), allocatable :: meanfield

    allocate(meanfield(n1,n2))

    meanfield = u*abs(psi(:,:))**2.

    ! x_evolution
    psi = exp(-ci*(uext(:,:) + meanfield(:,:))*dwt/2.)*psi

    ! p_evolution
    in = psi
    call fft_transform(forth)
    in = exp(-ci*dwt*0.5*psq(:,:))*out
    call fft_transform(back)
    psi = out

    ! x_evolution
    psi = exp(-ci*(uext(:,:)  + meanfield(:,:))*dwt/2.)*psi


    deallocate(meanfield)

  end subroutine one_step_evolution


  subroutine ener(etot)
    implicit none
    real(8) :: etot,ek,ep,eint

    ! kinetic and potential energy the two condensates
    in = psi(:,:)
    call fft_transform(forth)
    in = 0.5d0*psq*out
    call fft_transform(back)

    ek = product(dx)*sum(real(conjg(psi(:,:))*out))
    ep = product(dx)*sum(uext(:,:)*abs(psi(:,:))**2)

    ! mean-field energy
    eint = 0.5*u*product(dx)*sum(abs(psi(:,:))**4)

    ! total energy
    etot = ek + ep + eint

    write(10,'(5(g14.6,1x))') t, etot, ek, ep, eint

  end subroutine ener


  subroutine fft_transform(sign)
    implicit none
    integer :: sign ! FFTW_FORWARD=-1, FFTW_BACKWARD=+1

    if(sign.eq.1) then
       call dfftw_execute_dft(plan_b,in,out)
       out=out/product(dim)
    else
       call dfftw_execute_dft(plan_f,in,out)
    end if

  end subroutine fft_transform


  subroutine normalize(psi)
    implicit none
    complex(8), dimension(n1,n2) :: psi

    psi(:,:) = psi(:,:)/sqrt(product(dx)*sum(abs(psi(:,:))**2))

  end subroutine normalize


  ! adds opposite phases in the two channels
  subroutine adding_phase
    implicit none

    do i1 = 1, n1
       do i2 = 1, n2
          if (x2(i2).ge.-r2.and.x2(i2).le.0) then
             psi(i1,i2) = psi(i1,i2) * exp(ci*k0*x1(i1))
          else if(x2(i2).ge.0.and.x2(i2).le.r2) then
             psi(i1,i2) = psi(i1,i2) * exp(-ci*k0*x1(i1))
          else
             psi(i1,i2) = psi(i1,i2)
          end if
       end do
    end do
    call normalize(psi)

  end subroutine adding_phase


  subroutine wrt_dat(number)
    implicit none
    character(3) :: number
    real(8) :: phase
    real(8) :: maxpsi

    ! density
    open(UNIT=23,FILE="data/dens-"//number//".dat",STATUS='unknown')
    do i1 =1, n1
       do i2 = 1, n2
          write(23,'(2(2x,f10.4),2x,g16.4E3)') x1(i1), x2(i2), abs(psi(i1,i2))**2.
       end do
       write(23,*)
    end do
    close(23)

    ! phase
    open(UNIT=24,FILE="data/phase-"//number//".dat",STATUS='unknown')
    
    maxpsi = maxval(abs(psi))
    do i1 = 1, n1
       do i2 = 1, n2
          phase = atan(dimag(psi(i1,i2))/real(psi(i1,i2)))
          if(real(psi(i1,i2)).ge.0) then
             phase = phase
          else
             phase = phase + pi
          end if
          if (phase.lt.0.d0) then
             phase = phase + 2.d0*pi
          end if
          if(abs(psi(i1,i2)).lt.0.02d0*maxpsi) then
             phase=0.d0
          end if
          write(24,'(2(2x,f10.4),2x,g16.4E3)') x1(i1), x2(i2), phase
       end do
       write(24,*)
    end do
    close(24)

    ! FT along the longitudinal and transverse direction
    in = psi
    call fft_transform(forth)
    
    do i1 =1, n1
       write(41,'(2(2x,f10.4),2x,g16.4E3)') t, p1(i1), sum(abs(out(i1,:))**2.)
    end do
    write(41,*)

     do i2 =1, n2
       write(42,'(2(2x,f10.4),2x,g16.4E3)') t, p2(i2), sum(abs(out(:,i2))**2.)
    end do
    write(42,*)

  end subroutine wrt_dat


  subroutine wrt_pot
    implicit none

    open(UNIT=23,FILE="data/pot-x.dat",STATUS='unknown')
    do i1 = 1, n1
       write(23,'(2x,f10.4,4(2x,g16.8E3))') x1(i1),uext(i1,n2/2)
    end do
    close(23)

    open(UNIT=23,FILE="data/pot-y.dat",STATUS='unknown')
    do i2 = 1, n2
       write(23,'(2x,f10.4,4(2x,g16.8E3))') x2(i2),uext(n1/2,i2)
    end do
    close(23)

  end subroutine wrt_pot

end program main


function ederiv(lambda,param,nparam)
  ! derivative of E(lambda)
  implicit none
  integer :: nparam
  real(8), dimension(nparam) :: param
  real(8) :: ederiv,lambda

  ederiv = param(1) + lambda*param(2) + lambda**2*param(3) +lambda**3*param(4) + &
       & lambda**4*param(5) + lambda**5*param(6) + lambda**6*param(7)

end function ederiv

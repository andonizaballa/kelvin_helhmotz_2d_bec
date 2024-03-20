! (C) Copr. 1986-92 Numerical Recipes Software 3^03.
! rewritten in fortran90

subroutine zroots(ar,m,roots,polish)
  ! modified for real coefficients in input only
  implicit none
  integer :: m
  real(8), dimension(m+1) :: ar
  complex(8), dimension(m) :: roots
  complex(8), dimension(m+1) :: a
  real(8), parameter :: eps=5.d-11
  logical polish
  ! uses laguer
  complex(8), dimension(m+1) :: ad
  complex(8) :: x,b,c
  integer :: i,j,jj,its

  a = ar
  ad = a

  do j=m,1,-1
     x=cmplx(0.d0,0.d0)
     call laguer(ad,j,x,its)
     if(abs(imag(x)).le.2.d0*eps**2*abs(real(x))) then
        x=cmplx(real(x),0.d0)
     end if
     roots(j)=x
     b=ad(j+1)
     do jj=j,1,-1
        c=ad(jj)
        ad(jj)=b
        b=x*b+c
     end do
  end do

  if (polish) then
     do j=1,m
        call laguer(a,m,roots(j),its)
     end do
  endif

  do j=2,m
     x=roots(j)
     do i=j-1,1,-1
        if(real(roots(i)).le.real(x)) goto 10
        roots(i+1)=roots(i)
     end do
     i=0
10   roots(i+1)=x
  end do

end subroutine zroots

subroutine laguer(a,m,x,its)
  implicit none
  integer :: m,its
  integer, parameter :: MR=8
  integer, parameter :: MT=10
  integer, parameter :: MAXIT=MT*MR
  real(8), parameter, dimension (MR) :: frac=(/.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/)
  real(8), parameter :: EPSS=5.d-12
  complex(8), dimension(m+1) ::a 
  complex(8) :: x,dx,x1,b,d,f,g,h,sq,gp,gm,g2
  real(8) :: abx,abp,abm,err
  integer iter,j

  do iter=1,MAXIT
     its=iter
     b=a(m+1)
     err=abs(b)
     d=cmplx(0.,0.)
     f=cmplx(0.,0.)
     abx=abs(x)
     do j=m,1,-1
        f=x*f+d
        d=x*d+b
        b=x*b+a(j)
        err=abs(b)+abx*err
     end do
     err=EPSS*err
     if(abs(b).le.err) then
        return
     else
        g=d/b
        g2=g*g
        h=g2-2.d0*f/b
        sq=sqrt((m-1)*(m*h-g2))
        gp=g+sq
        gm=g-sq
        abp=abs(gp)
        abm=abs(gm)
        if(abp.lt.abm) gp=gm
        if (max(abp,abm).gt.0.) then
           dx=m/gp
        else
           dx=exp(cmplx(log(1.d0+abx),dfloat(iter)))
        endif
     endif
     x1=x-dx
     if(x.eq.x1)return
     if (mod(iter,MT).ne.0) then
        x=x1
     else
        x=x-dx*frac(iter/MT)
     endif
  end do

  write(*,*) 'too many iterations in laguer'
  stop

end subroutine laguer



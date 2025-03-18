parameter(jmax=5000)
implicit real*8 (a-h,o-z)
real*8, dimension(jmax) :: r(jmax),rr(jmax),vol(jmax),mnfw(jmax),&
        rhost(jmax),rho(jmax),mhern(jmax),rhonfw(jmax),mdark(jmax),&
        grvnfw(jmax),lnd(jmax),tem(jmax),tem2(jmax),mgas(jmax),&
        fbarr(jmax),fgasr(jmax),rhoBCG(jmax),lndBCG(jmax),mgasBCG(jmax),&
        grvnfwBCG(jmax),mgasBCG_grad(jmax),lndBCG_grad(jmax),rhoBCG_grad(jmax),&
        rho_rebusco(jmax), tem_rebusco(jmax)

real*8 :: msol,mu,mp,rmin,rmax,mvir,rvir,mbcg,ahern,lsol,h,me,&
          ne0,alphast,alphasn,zfesn

real*8, dimension(jmax) :: u(jmax),flux(jmax),ne(jmax),zfe(jmax),& 
        kappa(jmax),lturb,rhofedot(jmax),rhofe(jmax),zfest(jmax),&
        amfeiniz(jmax),amfe(jmax),gradzfe(jmax),zfeobs(jmax),&
        amfeobs(jmax),rhofeobs(jmax)

pi=3.14159265359
msol = 1.989d33
cmkpc = 3.084e21
years=3.156d7
mu=0.61
boltz=1.38066e-16
guniv=6.6720e-8
mp=1.67265e-24
zfesol=1.8e-3
zfesn=0.744/1.4
snu=0.5

tnow=13.7*1.e9*years
time0=tnow-5.*1.e9*years
time=time0      

!    set the grid

rmin = 0.1*cmkpc
rmax = 2800.*cmkpc !! or 3000 ?
do j=1,jmax
   r(j)=rmin+(j-1)*rmax/(jmax-1)
enddo
do j=1,jmax-1
   rr(j)=r(j)+0.5*(r(j+1)-r(j))
enddo
rr(jmax)=rr(jmax-1)+(rr(jmax-1)-rr(jmax-2))
open(10,file='grid.dat',status='unknown')
do j=1,jmax
   write(10,*)j,real(r(j)/cmkpc),real(rr(j)/cmkpc)
enddo
close(10)

vol(1)=4.1888*r(1)**3
do j=2,jmax
   vol(j)=4.1888*(r(j)**3-r(j-1)**3)    !! centered in rr(j-1) !!
enddo

rho0nfw=7.35d-26 !! Dark Matter
rs=435.7*cmkpc
rho0=4.d-26   !! 5000 points, without BCG, ISO
rho0BCG=8.d-26   !! 5000 points, with BCG , ISO!! !!in origin was 9.d-26
rho0BCG_grad=1.5d-25   !! 5000 points, with BCG and dT/dr !!
ticm=8.9e7

rvir=2797.*cmkpc
r500=rvir/2.
fc=1.138799
mvir=1.3e15*msol
mbcg=1.d12*msol
ahern=12.*cmkpc/(1.+2.**(0.5))
aml=7.5    !! this is the mass to light ratio 

do j=1,jmax
   x=rr(j)/rs
   y=rr(j)/cmkpc
   rhonfw(j)=rho0nfw/(x*(1.+x)**2) !! dark matter density profile NFW
   rhost(j)=mbcg/(2.*pi)*(ahern/rr(j))/(ahern+rr(j))**3 !!Stellar density profile
   rho_rebusco(j)=1.937d-24 * ((4.6d-2/(1+(y/57)**2)**1.8) + (4.8d-3/(1+(y/200)**2)**0.87)) !!Rebusco profile
enddo

open(20,file='masse.dat')
mnfw(1)=0.
mhern(1)=0.
do j=2,jmax
   x=r(j)/rs
   mnfw(j)=mnfw(j-1)+rhonfw(j-1)*vol(j) !!DM mass profile 
   mdark(j)=mvir*(log(1.+x)-x/(1.+x))/fc !!DM analytical mass profile
   mhern(j)=mbcg*r(j)**2/(r(j)+ahern)**2 !!Hernquist mass profile for the BCG
   write(20,1101)r(j)/cmkpc,mnfw(j)/msol,mdark(j)/msol,mhern(j)/msol
enddo
1101 format(4(1pe12.4))
close(20)

grvnfw(1)=0.
do j=2,jmax
   grvnfw(j)=guniv*(mnfw(j))/r(j)**2   !! without BCG !!
   grvnfwBCG(j)=guniv*(mnfw(j)+ mhern(j))/r(j)**2   !! with BCG !!
   write(20,1002)r(j)/cmkpc,grvnfw(j)/msol
enddo
1002 format(2(1pe12.4))
close(20)

!
! Temperature profile
!
  temp0=8.12e7
  rtemp1=71.
  rtemp2=71.
  r0bill=40.
  t0bill=0.9
  tmaxbill=4.8
  tcoeffbill=tmaxbill-t0bill
  pbill=1.8
  qbill=0.2  !!0.15
  sbill=1.6
  qplusp = qbill + pbill
  sinv=1./sbill
!
open(20,file='temperature.dat',status='unknown')
 do j=1,jmax
    rkpc=rr(j)/cmkpc
    x=rr(j)/r500
    y=rr(j)/cmkpc
    xx=x/0.045
    tem(j)=ticm*1.35*(xx**1.9+0.45)/(xx**1.9+1.)* &   !! this is for Perseus !!
        1./(1.+(x/0.6)**2)**0.45
   tem_rebusco(j)=7*1.16d7 * ((1+(y/71)**3)/(2.3+(y/71)**3))
    write(20,1003)rr(j)/cmkpc,tem(j),ticm,tem_rebusco(j)
 enddo
close(20)
1003 format(4(1pe12.4))

!     calculate the gas density, assuming ticm !! mette il gas in eq. con il potenziale

lnd(1)=log(rho0) !! NO BCG
lndBCG(1)=log(rho0BCG)         !! BCG
lndBCG_grad(1)=log(rho0BCG_grad)         !! BCG
do j=2,jmax
   gg=grvnfw(j)
   ggBCG=grvnfwBCG(j)
   temmed=0.5*(tem(j)+tem(j-1))
   lnd(j)=lnd(j-1)-gg*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm) !!NO BCG
   lndBCG(j)=lndBCG(j-1)-ggBCG*(mu*mp)*(rr(j)-rr(j-1))/(boltz*ticm) !! BCG
   lndBCG_grad(j)=lndBCG_grad(j-1)-ggBCG*(mu*mp)*(rr(j)-rr(j-1))/(boltz*temmed) & 
   - (log(tem(j)) - log(tem(j-1))) !! dT/dr
enddo

do j=1,jmax
   rho(j)=exp(lnd(j)) !! Gas density NO BCG
   rhoBCG(j)=exp(lndBCG(j)) !! Gas density BCG
   rhoBCG_grad(j)=exp(lndBCG_grad(j)) !!Gas density BCG dT/dr
   
enddo

open(20,file='density.dat',status='unknown')
do j=1,jmax
   write(20,1111)rr(j)/cmkpc,rho(j),rhonfw(j),rhoBCG(j),rhost(j),rhoBCG_grad(j),rho_rebusco(j)
enddo
1111 format(7(1pe12.4))
close(20)

open(20,file='mgas.dat',status='unknown')
mgas(1)=rho(j)*4.188*r(1)**3 !!NO BCG
mgasBCG(1)=rhoBCG(j)*4.188*r(1)**3 !! BCG
do j=2,jmax
   mgas(j)=mgas(j-1)+rho(j-1)*vol(j) !! NO BCG
   mgasBCG(j)=mgasBCG(j-1)+rhoBCG(j-1)*vol(j) !! BCG
   mgasBCG_grad(j)=mgasBCG_grad(j-1)+rhoBCG_grad(j-1)*vol(j) !! BCG
   write(20,1100)r(j)/cmkpc,mgas(j)/msol,mnfw(j)/msol
enddo
close(20)

!! fbar=mgas(jmax-1)/(mnfw(jmax-1)+mgas(jmax-1)) !!without BCG
!! fbar=(mhern(jmax-1)+mgasBCG(jmax-1))/(mnfw(jmax-1)+mgasBCG(jmax-1)+mhern(jmax-1)) !!with BCG
fbar=(mhern(jmax-1)+mgasBCG_grad(jmax-1))/(mnfw(jmax-1)+mgasBCG_grad(jmax-1)+mhern(jmax-1)) !!with BCG and dT/dr
fgas=(mgasBCG_grad(jmax-1))/(mnfw(jmax-1)+mgasBCG_grad(jmax-1)+mhern(jmax-1))
print*,'fbari, fgas = ',real(fbar), real(fgas)

open(20,file='barfrac.dat')
do j=2,jmax-1
   fbarr(j)=(mhern(j)+mgasBCG(j))/(mnfw(j)+mgasBCG(j)+mhern(j))
   fgasr(j)=(mgasBCG(j))/(mnfw(j)+mgasBCG(j)+mhern(j))
   write(20,1100)r(j)/cmkpc,fbarr(j),fgasr(j)
enddo
close(20)
1100 format(4(1pe12.4))

!***********************************************************************
!! At this point we have the gas density profile and we can proceed
!! with the integration of the diffusion equation for rhofe
!***********************************************************************

!! Set the initial abundance profile

zfeout=0.4*zfesol   !! this is the background abundance !!

do j=1,jmax
   x=rr(j)/(80.*cmkpc)
   zfeobs(j)=zfesol*0.3*(1.4*1.15*((2.2+x**3))/(1+x**3)/1.15)  !Perseus!
   zfe(j)=zfeobs(j) - zfeout !! subtract z_Fe,out  WITHOUT BACKGROUND
   !!zfe(j)=zfeout 
   !!zfe(j)=0.
   zfe(j)=max(zfe(j),0.)
   rhofe(j)=rhoBCG_grad(j)*zfe(j)/1.4   !! -Zfe_out
   rhofeobs(j)=rhoBCG_grad(j)*zfeobs(j)/1.4
enddo

 do j=1,jmax
    zfest(j)=1.*zfesol    !! set the stellar abundance !!
 enddo

 !! Calculate the initial excess of iron mass

amfeiniz(1)=rhofe(1)*vol(1)
amfeiniz(1)=rhofeobs(1)*vol(1)
do j=2,jmax
   amfeiniz(j)=amfeiniz(j-1)+rhofe(j-1)*vol(j)
   amfeobs(j)=amfeobs(j-1)+rhofeobs(j-1)*vol(j)
enddo

open(20,file='zfe_initial.dat')
do j=1,jmax
   write(20,1500)rr(j)/cmkpc,zfe(j)/zfesol,zfeobs(j)/zfesol, &
                 r(j)/cmkpc,amfeiniz(j)/msol,amfeobs(j)/msol
enddo
close(20)
1500 format(6(1pe12.4))

open(20,file='initial.dat',status='unknown')
do j=1,jmax
   write(20,3001)rr(j)/cmkpc+0.001,zfe(j)/zfesol,zfeobs(j)/zfesol,ne(j)
enddo
close(20)
3001  format(4(1pe12.4))

!! boundary conditions (outflows)

zfe(1)=zfe(2)
zfe(jmax)=zfe(jmax-1)
rhofe(1)=rhoBCG_grad(1)*zfe(1)/1.4
rhofe(jmax)=rhoBCG_grad(jmax)*zfe(jmax)/1.4

!! Here start the time integration (use FTCS method)

 print*,'dai ncycle'
 read(*,*)ncycle
 tend=tnow

!!  set the diffusion coefficient kappa = C*v*l (for now constant)

 open(20,file='kappa.dat',status='unknown')
 vturb=260.e5   !! come Perseus !!
 lturb=15.*cmkpc  !! this is quite uncertain !!
 rscala=30.*cmkpc

 kappa=0.11*vturb*lturb !! costant Kappa

!!do j=1,jmax    !! for variable kappa
!!    kappa(j)=rhost(j)  !!0.333*vturb*lturb   !! constant !!
!!    kappa(j)=kappa(j)-0.6*kappa(j)*exp(-(r(j)/rscala)**2)
!!    write(20,*)real(r(j)/cmkpc),kappa(j)
!! enddo
 close(20)

      n=0
1000  continue      !! here start the main time cycle
      n=n+1

!! calculate the timestep (to be modified if the grid is non-uniform)

      dt=0.4*(r(5)-r(4))**2/(2.*kappa(5))  !! ok for Delta_r costant !!
      time=time+dt
   !! print*,'n,dt (yr),time (Gyr) = ',n,real(dt/years), &
     !!       real(time/1.e9/years)

!! write the source terms (SNIa and stellar winds)

 slope=1.1
 alphast=4.7e-20*(time/tnow)**(-1.26)
 alphasn=4.436e-20*(snu/aml)*(time/tnow)**(-slope)

!!  print*,'alphast,sn = ',alphast,alphasn

 do j=2,jmax-2
    rhofedot(j)=(alphast*zfest(j)/1.4+alphasn*zfesn)*rhost(j)
 enddo

!! the equation to be solved is d(n*zfe)/dt = div(kappa*n*grad(zfe)) + S
!! (according to Rebusco et al. 2006)
!! Use the FTCS scheme.

!!      goto776
!! source step

 do j=2,jmax-1
!!!    write(70,*)rhofe(j),dt*rhofedot(j)
!!!    if(j.eq.5)print*,'azz ',dt,rhofe(j),dt*rhofedot(j),rhofedot(j)
    rhofe(j)=rhofe(j) + dt*rhofedot(j)
    zfe(j)=rhofe(j)/rhoBCG_grad(j) * 1.4
 enddo

!! set the boundary conditions (outflows)

      zfe(1)=zfe(2)
      zfe(jmax)=zfe(jmax-1)
      rhofe(1)=rhofe(2)
      rhofe(jmax)=rhofe(jmax-1)
776   continue

!!  goto777
!  diffusive step   !  check the Fe conservation !

 do j=2,jmax-1
    gradzfe(j)=(zfe(j)-zfe(j-1))/(rr(j)-rr(j-1))  !! dZ/dr centered at "j" !!
 enddo
 gradzfe(1)=0.        !! B.C. !!
 gradzfe(jmax)=0.

 do j=2,jmax-1
    rhojp1=0.5*(rhoBCG_grad(j+1)+rhoBCG_grad(j))  !! rho centered at "j+1" !!
    rhoj=0.5*(rhoBCG_grad(j-1)+rhoBCG_grad(j))    !! rho centered at "j" !!
    rhofe(j)=rhofe(j) &
            + (dt/1.4)*(r(j+1)**2*kappa(j+1)*rhojp1*gradzfe(j+1) &
            -r(j)**2*kappa(j)*rhoj*gradzfe(j))   &
             / (0.33333333*(r(j+1)**3-r(j)**3))
         zfe(j)=1.4*rhofe(j)/rhoBCG_grad(j)  !! update Z_Fe with the new rho_Fe !!
      enddo
2000  format(3(1pe12.4))

!! set the boundary conditions (outflows)

      zfe(1)=zfe(2)
      zfe(jmax)=zfe(jmax-1)
      rhofe(1)=rhofe(2)
      rhofe(jmax)=rhofe(jmax-1)
777   continue

      if (time.ge.tend) goto1001
      if (n.ge.ncycle) goto1001

      goto 1000

1001  continue

      do j=2,jmax
         write(99,*)real(log10(r(j)/cmkpc)),real(log10(rhofedot(j)))
      enddo

!! calcola la massa di Fe al tempo finale

      amfe(1)=rhofe(1)*vol(1)
      do j=2,jmax
         amfe(j)=amfe(j-1)+rhofe(j-1)*vol(j)
      enddo

      write(6,3002)amfe(jmax)/msol,amfeiniz(jmax)/msol,amfeobs(jmax)/msol
      write(6,3003)amfe(180)/msol,amfeiniz(180)/msol,amfeobs(180)/msol
3002  format('M_Fe(tot), M_Fein(tot) (Msol) = ',3(1pe12.4))
3003  format('M_Fe(<100kpc), M_Fein(<100kpc) (Msol) = ',3(1pe12.4))

      print*,'TIME (Gyr) = ',time/3.156e16

      open(21,file='diff.dat',status='unknown')
      do j=2,jmax
         write(21,3000)rr(j)/cmkpc,zfe(j)/zfesol
      enddo
      close(21)
3000  format(2(1pe12.4))

stop
end

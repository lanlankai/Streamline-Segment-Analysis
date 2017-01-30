      program main

      integer me,ll,lh,npu

c the defination of variables
      character(len=49) path
      character(len=3) fileid
      character(len=8) fname

      double precision jumpeps,pace,err
      integer ndimx,ndimy,ndimz,nt
      double precision dx,dy,dz

      double precision tb,tu
      double precision jumpratio
      double precision dzmean,dlmean,curvmean,meanz,smean
      double precision mark(500000)
      integer up_down1,up_down2

      double precision deriv

      double precision umean,lmean
      double precision lengs_du(0:100,-100:100)
     . ,lengs_curv(0:100,-100:100),lengs(0:100)
     . ,lengs_cos(0:100,-100:100),lengsp(0:100),lengsn(0:100)
     . ,du_curv(-100:100,-100:100),cos_curv(-100:100,-100:100)
     . ,theta_curv(-100:100,-100:100),theta1_curv(-100:100,-100:100)
     . ,theta_str(-100:100,-100:100),theta1_str(-100:100,-100:100)
     . ,mpdf(-100:100,5),apdf(0:100)
      double precision lengs_du_log(-100:100,-100:100)
      double precision facti,factj,vec1(3),vec2(3),vec3(3)

      double precision temp1,temp2,temp3,zmax,zmin,fmax,fmin
      double precision temp01,temp02,temp03,velp,veln
      integer px,py,pz
      double precision sx,sy,sz
      double precision t11,t12,t13,t21,t22,t23,t31,t32,t33
     . ,nor1,nor2,nor3

      double precision temp,plen,nlen,posold(3),tp1,tp2,tp3,tp,pos0(3)
      integer endid,switch
      double precision posx,posy,posz,endx,endy,endz
      double precision subdvx(-1:1,-1:1,-1:1),subdvy(-1:1,-1:1,-1:1),
     . subdvz(-1:1,-1:1,-1:1)
      double precision subu1(-1:0,-1:1,-1:1),subu2(-1:1,-1:0,-1:1),
     . subu3(-1:1,-1:1,-1:0)
      double precision subv1(-1:0,-1:1,-1:1),subv2(-1:1,-1:0,-1:1),
     . subv3(-1:1,-1:1,-1:0)
      double precision subw1(-1:0,-1:1,-1:1),subw2(-1:1,-1:0,-1:1),
     . subw3(-1:1,-1:1,-1:0)
      integer i,j,k,ni,nj,nk,endfind,saddle,i0,j0,k0
      integer numends,findex,loop,force0
      character(len=2) title
      integer status
      double precision, dimension(:,:,:), pointer:: psi
      double precision, dimension(:,:,:), pointer:: u_t,v_t
      double precision, dimension(:,:,:), pointer:: u,v,w
      double precision, dimension(:,:,:), pointer:: curv,area
     . ,traj,trajs,strain
      double precision, dimension(:,:,:,:),pointer:: endlocp,endlocn
      double precision, dimension(:,:,:,:), pointer:: omg,pomg,domg
      double precision, dimension(:,:,:,:), pointer:: trajsbu

      integer, dimension(:,:,:), pointer:: cell
      double precision, dimension(:,:), pointer:: points
      integer, dimension(:,:,:,:), pointer:: numbering
      integer, dimension(:,:,:,:), pointer:: identity
      double precision, dimension(:,:,:,:), pointer:: counter
      real :: start, finish,finish_old,Ttime

      common/large0/ psi
      common/large1/ numbering,counter,identity
      common/large2/ points
      common/large3/ curv,area,traj,trajs,strain
      common/large4/ cell
      common/large5/ u,v,w

      common/global1/ jumpeps,pace,err
      common/global2/ ndimx,ndimy,ndimz,nt
      common/global3/ dx,dy,dz

      open(unit=16,file='param',status='unknown')
      read(16,*) ndimx,ndimy
      read(16,*) nt
      read(16,*) dx,dy
      read(16,*) dzmean,dlmean
      read(16,"(A)") fname
      close(16)

      write(*,*) 'ndimx is',ndimx
      write(*,*) 'ndimy is',ndimy
      write(*,*) 'Total number of snaphost is',nt
      write(*,'((A32,es14.4))') 'Maximum velocity difference is ',dzmean
      write(*,'((A25,es14.4))') 'Maximum length scale is ',dlmean
      write(*,*) 'File name is:',fname
      write(*,*) 'Please confirm above parameters by press ENTER'
      write(*,*) 'Or press CTRL+C to cancel the calculation!'

      read(*,*)

      write(*,*) 'we are begining to do the calculation'
      call cpu_time(start)


      loop=0
c      dzmean=0.1d0
c      dlmean=30.d0

      numends=0
      jumpratio=0.d0

      jumpeps=2.d-5
      pace=0.02
      err=0.1d0
      dz=dx
      ndimz=3
      allocate(u_t(ndimx,ndimy,nt),STAT=status)
      allocate(v_t(ndimx,ndimy,nt),STAT=status)
      allocate(u(ndimx,ndimy,ndimz),STAT=status)
      allocate(v(ndimx,ndimy,ndimz),STAT=status)
      allocate(w(ndimx,ndimy,ndimz),STAT=status)

      allocate(curv(ndimx,ndimy,ndimz),STAT=status)
      allocate(area(ndimx,ndimy,ndimz),STAT=status)
      allocate(traj(ndimx,ndimy,ndimz),STAT=status)
      allocate(trajs(ndimx,ndimy,ndimz),STAT=status)
      allocate(trajsbu(ndimx,ndimy,ndimz,6),STAT=status)
      allocate(strain(ndimx,ndimy,ndimz),STAT=status)

      allocate(cell(ndimx,ndimy,ndimz),STAT=status)
      allocate(points(500000,3),STAT=status)
      allocate(numbering(ndimx,ndimy,ndimz,2),STAT=status)
      allocate(counter(0:ndimx+1,0:ndimy+1,0:ndimz+1,4),STAT=status)
      allocate(identity(0:ndimx+1,0:ndimy+1,0:ndimz+1,2),STAT=status)

      allocate(endlocp(ndimx,ndimy,ndimz,4),STAT=status)
      allocate(endlocn(ndimx,ndimy,ndimz,4),STAT=status)

      allocate(omg(ndimx,ndimy,ndimz,3),STAT=status)
      allocate(pomg(ndimx,ndimy,ndimz,3),STAT=status)
      allocate(domg(ndimx,ndimy,ndimz,3),STAT=status)

      do k=1,ndimz
      do j=1,ndimy
      do i=1,ndimx
      numbering(i,j,k,1)=0
      numbering(i,j,k,2)=0
      enddo
      enddo
      enddo

      do i=0,ndimx+1
      do j=0,ndimy+1
      do k=0,ndimz+1
      counter(i,j,k,1)=0.0
      counter(i,j,k,2)=0.0
      counter(i,j,k,3)=0.0
      counter(i,j,k,4)=0.0
      identity(i,j,k,1)=0
      identity(i,j,k,2)=0
      enddo
      enddo
      enddo

      loop=loop+1
      allocate(psi(ndimx,ndimy,ndimz),STAT=status)

! =====================================
      open(unit=16,file=fname,status='unknown')
      do me=1,nt
      do j=1,ndimy
      do i=1,ndimx
      read(16,*) temp,temp,u_t(i,j,me),v_t(i,j,me)
      end do
      end do
      enddo
      close(16)

c      open(60,file='./test.dat',action='write')
c      write(60,*) '"x", "y", "u", "v"'
c      write(60,*) 'ZONE I=128, J=128, F=POINT'
c      do j=1,ndimy
c      do i=1,ndimx
c      write(60,'(4(es14.4))') i*1.d0,j*1.d0,u_t(i,j,1),v_t(i,j,1)
c      enddo
c      enddo
c      close(60)

      do k=1,ndimz
      do j=1,ndimy
      do i=1,ndimx
      curv(i,j,k)=0.d0
      enddo
      enddo
      enddo

      do j=-100,100
      mpdf(j,1)=0.d0
      mpdf(j,2)=0.d0
      mpdf(j,3)=0.d0
      mpdf(j,4)=0.d0
      mpdf(j,5)=0.d0
      enddo

      do j=-100,100
      do i=0,100
      lengs_du(i,j)=0.d0
      enddo
      do i=-100,100
      lengs_du_log(i,j)=0.d0
      enddo
      enddo
      do i=0,100
      lengs(i)=0.d0
      lengsp(i)=0.d0
      lengsn(i)=0.d0
      enddo

      call cpu_time(finish)
      write(*,*) 'It takes ',finish-start, ' second to loading data!'

      call cpu_time(finish_old)
      do me=1,nt
      write(*,*) 'Now we are doing snapshot:',me
      call cpu_time(finish)
      write(*,*) 'Time = ',finish-start
      if(me==1)finish_old=finish
      if(me==2)finish_old=finish-finish_old
c      if(me==2) write(*,*) (nt-me)*finish_old
      if(me>1)Ttime=real(nt-me)*finish_old
      if(me>1.and.me<nt)then
         write(*,*) 'The rest still need  ',Ttime, ' second!'
       endif

      do i=1,ndimx
      do j=1,ndimy
      u(i,j,1)=u_t(i,j,me)
      u(i,j,2)=u(i,j,1)
      u(i,j,3)=u(i,j,1)
      v(i,j,1)=v_t(i,j,me)
      v(i,j,2)=v(i,j,1)
      v(i,j,3)=v(i,j,1)
      w(i,j,1)=0.d0
      w(i,j,2)=0.d0
      w(i,j,3)=0.d0
      end do
      end do

      do k=2,ndimz-1
      do j=2,ndimy-1
      do i=2,ndimx-1
      trajs(i,j,k)=0.d0
      trajsbu(i,j,k,1)=0.d0
      trajsbu(i,j,k,2)=0.d0
      trajsbu(i,j,k,3)=0.d0
      trajsbu(i,j,k,4)=0.d0
      trajsbu(i,j,k,5)=0.d0
      trajsbu(i,j,k,6)=0.d0
      enddo
      enddo
      enddo

      do k=2,ndimz-1
!      write(*,*) k
      do j=2,ndimy-1
      do i=2,ndimx-1

       pos0(1)=i
       pos0(2)=j
       pos0(3)=k

       endid=0
c switch=1 means we do calculation in positive direction
       switch=1
       force0=0
       endfind=0
       posold(1)=pos0(1)
       posold(2)=pos0(2)
       posold(3)=pos0(3)
       ni=posold(1)
       nj=posold(2)
       nk=posold(3)
       posx=posold(1)-ni
       posy=posold(2)-nj
       posz=posold(3)-nk
       plen=0.d0
       up_down1=0

       do while(endfind.eq.0.and.(ni.gt.1).and.(ni.lt.ndimx).and.
     .  (nj.gt.1).and.(nj.lt.ndimy).and.(nk.gt.1).and.(nk.lt.ndimz))

c give the velocity value to sub arrays
       do i0=-1,1
       do j0=-1,1
       do k0=-1,1
       subdvx(i0,j0,k0)=u(ni+i0,nj+j0,nk+k0)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,1
       subdvy(i0,j0,k0)=v(ni+i0,nj+j0,nk+k0)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,1
       subdvz(i0,j0,k0)=0.d0
       enddo
       enddo
       enddo

c give the derivative value to sub arrays
       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subu1(i0,j0,k0)=
     . (u(ni+i0+1,nj+j0,nk+k0)-u(ni+i0,nj+j0,nk+k0))/dx
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subu2(i0,j0,k0)=
     . (u(ni+i0,nj+j0+1,nk+k0)-u(ni+i0,nj+j0,nk+k0))/dy
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subu3(i0,j0,k0)=
     . (u(ni+i0,nj+j0,nk+k0+1)-u(ni+i0,nj+j0,nk+k0))/dz
       enddo
       enddo
       enddo

       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subv1(i0,j0,k0)=
     . (v(ni+i0+1,nj+j0,nk+k0)-v(ni+i0,nj+j0,nk+k0))/dx
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subv2(i0,j0,k0)=
     . (v(ni+i0,nj+j0+1,nk+k0)-v(ni+i0,nj+j0,nk+k0))/dy
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subv3(i0,j0,k0)=
     . (v(ni+i0,nj+j0,nk+k0+1)-v(ni+i0,nj+j0,nk+k0))/dz
       enddo
       enddo
       enddo

       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subw1(i0,j0,k0)=
     . (w(ni+i0+1,nj+j0,nk+k0)-w(ni+i0,nj+j0,nk+k0))/dx
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subw2(i0,j0,k0)=
     . (w(ni+i0,nj+j0+1,nk+k0)-w(ni+i0,nj+j0,nk+k0))/dy
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subw3(i0,j0,k0)=
     . (w(ni+i0,nj+j0,nk+k0+1)-w(ni+i0,nj+j0,nk+k0))/dz
       enddo
       enddo
       enddo

       up_down1=up_down1+1
       call endpoints(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,subu1,subu2,subu3,subv1,subv2,subv3,
     . subw1,subw2,subw3,endfind,jumpratio,force0,up_down1,deriv)

       tp1=(ni+endx-posold(1))*dx
       tp2=(nj+endy-posold(2))*dy
       tp3=(nk+endz-posold(3))*dz
       posold(1)=ni+endx
       posold(2)=nj+endy
       posold(3)=nk+endz
       plen=plen+dsqrt(tp1**2+tp2**2+tp3**2)

       call tellnews(ni,nj,nk,endx,endy,endz,posx,posy,posz)

       force0=force0+1
       if(force0.gt.500) endfind=1

       if(endfind.eq.1) then
       endid=endid+1

       px=ni+posx
       py=nj+posy
       pz=nk+posz
       sx=ni+posx-px
       sy=nj+posy-py
       sz=nk+posz-pz

       if(px.eq.0.or.py.eq.0.or.pz.eq.0) write(*,*) i,j,k

       t11=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*u(px,py,pz)+
     . sx*(1.d0-sy)*(1.d0-sz)*u(px+1,py,pz)+(1.d0-sx)*sy*(1.d0-sz)
     . *u(px,py+1,pz)+(1.d0-sx)*(1.d0-sy)*sz*u(px,py,pz+1)+
     . (1.d0-sx)*sy*sz*u(px,py+1,pz+1)+sx*(1.d0-sy)*sz
     . *u(px+1,py,pz+1)+sx*sy*(1.d0-sz)*u(px+1,py+1,pz)+
     . sx*sy*sz*u(px+1,py+1,pz+1)
       t12=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*v(px,py,pz)+
     . sx*(1.d0-sy)*(1.d0-sz)*v(px+1,py,pz)+(1.d0-sx)*sy*(1.d0-sz)
     . *v(px,py+1,pz)+(1.d0-sx)*(1.d0-sy)*sz*v(px,py,pz+1)+
     . (1.d0-sx)*sy*sz*v(px,py+1,pz+1)+sx*(1.d0-sy)*sz
     . *v(px+1,py,pz+1)+sx*sy*(1.d0-sz)*v(px+1,py+1,pz)+
     . sx*sy*sz*v(px+1,py+1,pz+1)
       t13=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*w(px,py,pz)+
     . sx*(1.d0-sy)*(1.d0-sz)*w(px+1,py,pz)+(1.d0-sx)*sy*(1.d0-sz)
     . *w(px,py+1,pz)+(1.d0-sx)*(1.d0-sy)*sz*w(px,py,pz+1)+
     . (1.d0-sx)*sy*sz*w(px,py+1,pz+1)+sx*(1.d0-sy)*sz
     . *w(px+1,py,pz+1)+sx*sy*(1.d0-sz)*w(px+1,py+1,pz)+
     . sx*sy*sz*w(px+1,py+1,pz+1)
       velp=dsqrt(t11**2+t12**2+t13**2)

!$omp critical
       endlocp(i,j,k,1)=ni+posx
       endlocp(i,j,k,2)=nj+posy
       endlocp(i,j,k,3)=nk+posz
       endlocp(i,j,k,4)=velp
       trajs(i,j,k)=plen
!$omp end critical
       endif

       enddo

c switch=-1 means we do calculation in negative direction
       switch=-1

       force0=0
       endfind=0
       posold(1)=pos0(1)
       posold(2)=pos0(2)
       posold(3)=pos0(3)
       ni=posold(1)
       nj=posold(2)
       nk=posold(3)
       posx=posold(1)-ni
       posy=posold(2)-nj
       posz=posold(3)-nk
       nlen=0.d0
       up_down2=0

       do while(endfind.eq.0.and.(ni.gt.1).and.(ni.lt.ndimx).and.
     .  (nj.gt.1).and.(nj.lt.ndimy).and.(nk.gt.1).and.(nk.lt.ndimz))

       do i0=-1,1
       do j0=-1,1
       do k0=-1,1
       subdvx(i0,j0,k0)=u(ni+i0,nj+j0,nk+k0)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,1
       subdvy(i0,j0,k0)=v(ni+i0,nj+j0,nk+k0)
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,1
       subdvz(i0,j0,k0)=w(ni+i0,nj+j0,nk+k0)
       enddo
       enddo
       enddo

c give the derivative value to sub arrays
       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subu1(i0,j0,k0)=
     . (u(ni+i0+1,nj+j0,nk+k0)-u(ni+i0,nj+j0,nk+k0))/dx
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subu2(i0,j0,k0)=
     . (u(ni+i0,nj+j0+1,nk+k0)-u(ni+i0,nj+j0,nk+k0))/dy
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subu3(i0,j0,k0)=
     . (u(ni+i0,nj+j0,nk+k0+1)-u(ni+i0,nj+j0,nk+k0))/dz
       enddo
       enddo
       enddo

       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subv1(i0,j0,k0)=
     . (v(ni+i0+1,nj+j0,nk+k0)-v(ni+i0,nj+j0,nk+k0))/dx
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subv2(i0,j0,k0)=
     . (v(ni+i0,nj+j0+1,nk+k0)-v(ni+i0,nj+j0,nk+k0))/dy
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subv3(i0,j0,k0)=
     . (v(ni+i0,nj+j0,nk+k0+1)-v(ni+i0,nj+j0,nk+k0))/dz
       enddo
       enddo
       enddo

       do i0=-1,0
       do j0=-1,1
       do k0=-1,1
       subw1(i0,j0,k0)=
     . (w(ni+i0+1,nj+j0,nk+k0)-w(ni+i0,nj+j0,nk+k0))/dx
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,0
       do k0=-1,1
       subw2(i0,j0,k0)=
     . (w(ni+i0,nj+j0+1,nk+k0)-w(ni+i0,nj+j0,nk+k0))/dy
       enddo
       enddo
       enddo

       do i0=-1,1
       do j0=-1,1
       do k0=-1,0
       subw3(i0,j0,k0)=
     . (w(ni+i0,nj+j0,nk+k0+1)-w(ni+i0,nj+j0,nk+k0))/dz
       enddo
       enddo
       enddo

       up_down2=up_down2+1
       call endpoints(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,subu1,subu2,subu3,subv1,subv2,subv3,
     . subw1,subw2,subw3,endfind,jumpratio,force0,up_down2,deriv)

       tp1=(ni+endx-posold(1))*dx
       tp2=(nj+endy-posold(2))*dy
       tp3=(nk+endz-posold(3))*dz
       posold(1)=ni+endx
       posold(2)=nj+endy
       posold(3)=nk+endz
       nlen=nlen+dsqrt(tp1**2+tp2**2+tp3**2)

       call tellnews(ni,nj,nk,endx,endy,endz,posx,posy,posz)

       force0=force0+1
       if(force0.gt.500) endfind=1

       if(endfind.eq.1) then
       endid=endid+1

       px=ni+posx
       py=nj+posy
       pz=nk+posz
       sx=ni+posx-px
       sy=nj+posy-py
       sz=nk+posz-pz

       t11=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*u(px,py,pz)+
     . sx*(1.d0-sy)*(1.d0-sz)*u(px+1,py,pz)+(1.d0-sx)*sy*(1.d0-sz)
     . *u(px,py+1,pz)+(1.d0-sx)*(1.d0-sy)*sz*u(px,py,pz+1)+
     . (1.d0-sx)*sy*sz*u(px,py+1,pz+1)+sx*(1.d0-sy)*sz
     . *u(px+1,py,pz+1)+sx*sy*(1.d0-sz)*u(px+1,py+1,pz)+
     . sx*sy*sz*u(px+1,py+1,pz+1)
       t12=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*v(px,py,pz)+
     . sx*(1.d0-sy)*(1.d0-sz)*v(px+1,py,pz)+(1.d0-sx)*sy*(1.d0-sz)
     . *v(px,py+1,pz)+(1.d0-sx)*(1.d0-sy)*sz*v(px,py,pz+1)+
     . (1.d0-sx)*sy*sz*v(px,py+1,pz+1)+sx*(1.d0-sy)*sz
     . *v(px+1,py,pz+1)+sx*sy*(1.d0-sz)*v(px+1,py+1,pz)+
     . sx*sy*sz*v(px+1,py+1,pz+1)
       t13=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*w(px,py,pz)+
     . sx*(1.d0-sy)*(1.d0-sz)*w(px+1,py,pz)+(1.d0-sx)*sy*(1.d0-sz)
     . *w(px,py+1,pz)+(1.d0-sx)*(1.d0-sy)*sz*w(px,py,pz+1)+
     . (1.d0-sx)*sy*sz*w(px,py+1,pz+1)+sx*(1.d0-sy)*sz
     . *w(px+1,py,pz+1)+sx*sy*(1.d0-sz)*w(px+1,py+1,pz)+
     . sx*sy*sz*w(px+1,py+1,pz+1)

       veln=dsqrt(t11**2+t12**2+t13**2)
!$omp critical
       endlocn(i,j,k,1)=ni+posx
       endlocn(i,j,k,2)=nj+posy
       endlocn(i,j,k,3)=nk+posz
       endlocn(i,j,k,4)=veln
       trajs(i,j,k)=plen+nlen
!$omp end critical
       endif

       enddo

       if(endid.lt.1.5) trajs(i,j,k)=0.d0

 300   temp=0

       enddo
       enddo

      enddo

      do k=2,ndimz-1
      do j=2,ndimy-1
      do i=2,ndimx-1
      if(trajs(i,j,k).gt.1.d-8) then
      ni=trajs(i,j,k)/dlmean*20
      ni=min(ni,100)
c every quantity should be weighted by the flame surface area
      lengs(ni)=lengs(ni)+1

      nj=(endlocp(i,j,k,4)-endlocn(i,j,k,4))/dzmean*30
      nj=min(nj,100)
      nj=max(nj,-100)
      if(nj.eq.0) then
      factj=0.5d0
      else
      factj=1.d0
      endif
      lengs_du(ni,nj)=lengs_du(ni,nj)+factj

      ni=trajs(i,j,k)/dlmean*20
      ni=min(ni,100)
      if((endlocp(i,j,k,4)-endlocn(i,j,k,4)).ge.0.d0) then
      lengsp(ni)=lengsp(ni)+1
      else
      lengsn(ni)=lengsn(ni)+1
      endif


      ni=log(trajs(i,j,k))*20
      ni=min(ni,100)
      ni=max(ni,-100)
      nj=log(dabs(endlocp(i,j,k,4)-endlocn(i,j,k,4)))*20
      nj=min(nj,100)
      nj=max(nj,-100)
      if(ni.eq.0) then
      facti=0.5d0
      else
      facti=1.d0
      endif
      if(nj.eq.0) then
      factj=0.5d0
      else
      factj=1.d0
      endif
      lengs_du_log(ni,nj)=lengs_du_log(ni,nj)+facti*factj
      endif

      enddo
      enddo
      enddo

! this is the loop for different files
      enddo

      open(60,file='./length.dat',action='write')
      write(60,*) '%"lengs", "pdf"'
      write(60,*) '%ZONE I=101, F=POINT'
      do i=0,100
      write(60,'(4(es14.4))') i/20.d0*dlmean,lengs(i),
     . lengsp(i),lengsn(i)
      enddo
      close(60)

      open(60,file='./len_dv.dat',action='write')
      write(60,*) '%"lengs", "du", "pdf"'
      write(60,*) '%ZONE I=101, J=201, F=POINT'
      do j=-100,100
      do i=0,100
      write(60,'(3(es14.4))') i/20.d0*dlmean,j/30.d0*dzmean,
     . lengs_du(i,j)
      enddo
      enddo
      close(60)

      open(60,file='./len_dv_log.dat',action='write')
      write(60,*) '%"lengs", "du", "pdf"'
      write(60,*) '%ZONE I=201, J=201, F=POINT'
      do j=-100,100
      do i=-100,100
      write(60,'(3(es14.4))') i/20.d0,j/20.d0,
     . lengs_du_log(i,j)
      enddo
      enddo
      close(60)

      call cpu_time(finish)
      write(*,*) 'Total time is: ',finish-start

      write(*,*) 'all things are done'

      deallocate(psi)
      deallocate(points)
      deallocate(numbering)

!      enddo

      end


      subroutine tellnew(ni,nj,nk,endx,endy,endz,posx,posy,posz)
      integer ni,nj,nk
      double precision endx,endy,endz,posx,posy,posz

      if(abs(endx).le.0.5.and.abs(endy).le.0.5.and.abs(endz).le.0.5)
     . then
      ni=ni
      nj=nj
      nk=nk
      posx=endx
      posy=endy
      posz=endz
      elseif(endx.gt.0.5.and.abs(endy).le.0.5.and.abs(endz).le.0.5)then
      ni=ni+1
      nj=nj
      nk=nk
      posx=endx-1.d0
      posy=endy
      posz=endz
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5)
     .  .and.(abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj
      nk=nk
      posx=endx+1.d0
      posy=endy
      posz=endz
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni
      nj=nj+1
      nk=nk
      posx=endx
      posy=endy-1.d0
      posz=endz
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni
      nj=nj-1
      nk=nk
      posx=endx
      posy=endy+1.d0
      posz=endz
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz
      elseif((abs(endx).le.0.5).and.(abs(endy).le.0.5).and.
     . (endz.gt.0.5)) then
      ni=ni
      nj=nj
      nk=nk+1
      posx=endx
      posy=endy
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(abs(endy).le.0.5).and.(endz.gt.0.5))
     . then
      ni=ni+1
      nj=nj
      nk=nk+1
      posx=endx-1.d0
      posy=endy
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5).and.
     . (endz.gt.0.5))then
      ni=ni-1
      nj=nj
      nk=nk+1
      posx=endx+1.d0
      posy=endy
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.(endz.gt.0.5))
     . then
      ni=ni
      nj=nj+1
      nk=nk+1
      posx=endx
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.
     . (endz.gt.0.5)) then
      ni=ni
      nj=nj-1
      nk=nk+1
      posx=endx
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.(endz.gt.0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk+1
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.(endz.gt.0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk+1
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.(endz.gt.0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk+1
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.(endz.gt.0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk+1
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz-1.d0

      elseif((abs(endx).le.0.5).and.(abs(endy).le.0.5).and.
     . (endz.lt.-0.5)) then
      ni=ni
      nj=nj
      nk=nk-1
      posx=endx
      posy=endy
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(abs(endy).le.0.5).and.
     . (endz.lt.-0.5)) then
      ni=ni+1
      nj=nj
      nk=nk-1
      posx=endx-1.d0
      posy=endy
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni-1
      nj=nj
      nk=nk-1
      posx=endx+1.d0
      posy=endy
      posz=endz+1.d0
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni
      nj=nj+1
      nk=nk-1
      posx=endx
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni
      nj=nj-1
      nk=nk-1
      posx=endx
      posy=endy+1.d0
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk-1
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk-1
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk-1
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk-1
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz+1.d0

      endif

      return
      end

      subroutine endpoint(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,endfind,saddle,jumpratio)

      double precision jumpratio
      double precision posx,posy,posz,nposx,nposy,nposz,endx,endy,endz
      double precision angle(3),temp,temp1,temp2,temp3
      integer switch
      double precision subdvx(-1:1,-1:1,-1:1),subdvy(-1:1,-1:1,-1:1),
     . subdvz(-1:1,-1:1,-1:1)
      double precision normx,normy,normz,normt,jump
      integer ix,iy,iz,inner,endfind,saddle,force
      double precision xp,xm,yp,ym,zp,zm
      double precision vxp(3),vxm(3),vyp(3),vym(3),vzp(3),vzm(3)
      double precision diverg,totv,vsum(3)

      double precision jumpeps,pace,err
      integer ndimx,ndimy,ndimz,nt
      double precision dx,dy,dz

      common/global1/ jumpeps,pace,err
      common/global2/ ndimx,ndimy,ndimz,nt
      common/global3/ dx,dy,dz

      inner=1
      force=0

      nposx=posx
      nposy=posy
      nposz=posz

      normx=0.d0
      normy=0.d0
      normz=0.d0
      do ix=-1,0
      do iy=-1,1
      do iz=-1,1
      normx=max(dabs(subdvx(ix,iy,iz)),normx)
      enddo
      enddo
      enddo
      do ix=-1,1
      do iy=-1,0
      do iz=-1,1
      normy=max(dabs(subdvy(ix,iy,iz)),normy)
      enddo
      enddo
      enddo
      do ix=-1,1
      do iy=-1,1
      do iz=-1,0
      normz=max(dabs(subdvz(ix,iy,iz)),normz)
      enddo
      enddo
      enddo
      normt=dsqrt(normx*normx+normy*normy+normz*normz)+1.d-20

      do while(inner.eq.1)

      if((abs(nposx).gt.0.5d0).or.(abs(nposy).gt.0.5d0).or.
     . (abs(nposz).gt.0.5d0)) then
      inner=0
      endx=nposx
      endy=nposy
      endz=nposz
      goto 54
      endif

      call findangle(nposx,nposy,nposz,switch,
     . subdvx,subdvy,subdvz,angle,totv)

      jump=min(pace,totv/normt/4.d0)

      if(jump.lt.jumpeps) then

       xp=nposx+pace
       xm=nposx-pace
       yp=nposy+pace
       ym=nposy-pace
       zp=nposz+pace
       zm=nposz-pace
       call findangle(xp,nposy,nposz,switch,subdvx,subdvy,subdvz,
     .  vxp,temp)
       call findangle(xm,nposy,nposz,switch,subdvx,subdvy,subdvz,
     .  vxm,temp)
       call findangle(nposx,yp,nposz,switch,subdvx,subdvy,subdvz,
     .  vyp,temp)
       call findangle(nposx,ym,nposz,switch,subdvx,subdvy,subdvz,
     .  vym,temp)
       call findangle(nposx,nposy,zp,switch,subdvx,subdvy,subdvz,
     .  vzp,temp)
       call findangle(nposx,nposy,zm,switch,subdvx,subdvy,subdvz,
     .  vzm,temp)
       diverg=vxp(1)-vxm(1)+vyp(2)-vym(2)+vzp(3)-vzm(3)
       if(diverg.gt.-3.5) then

        if(totv.lt.1.d-10) then
        saddle=1
        endif

       call random_number(vsum(1))
       call random_number(vsum(2))
       call random_number(vsum(3))

        if(abs(angle(1)).ge.abs(angle(2)).
     .  and.abs(angle(1)).ge.abs(angle(3))) then
        temp2=vsum(2)**2*angle(2)
        temp3=vsum(3)**2*angle(3)
        temp1=10*(temp2*angle(2)+temp3*angle(3))/(angle(1)+1.d-20)
        temp=1.d-20+dsqrt(temp1**2+temp2**2+temp3**2)
        nposx=nposx-temp1/temp*pace
        nposy=nposy+temp2/temp*pace
        nposz=nposz+temp3/temp*pace
        goto 55
        elseif(abs(angle(2)).ge.abs(angle(1)).
     .  and.abs(angle(2)).ge.abs(angle(3))) then
        temp1=vsum(1)**2*angle(1)
        temp3=vsum(3)**2*angle(3)
        temp2=10*(temp1*angle(1)+temp3*angle(3))/(angle(2)+1.d-20)
        temp=1.d-20+dsqrt(temp1**2+temp2**2+temp3**2)
        nposx=nposx+temp1/temp*pace
        nposy=nposy-temp2/temp*pace
        nposz=nposz+temp3/temp*pace
        goto 55
        elseif(abs(angle(3)).ge.abs(angle(1)).
     .  and.abs(angle(3)).ge.abs(angle(2))) then
        temp1=vsum(1)**2*angle(1)
        temp2=vsum(2)**2*angle(2)
        temp3=10*(temp1*angle(1)+temp2*angle(2))/(angle(3)+1.d-20)
        temp=1.d-20+dsqrt(temp1**2+temp2**2+temp3**2)
        nposx=nposx+temp1/temp*pace
        nposy=nposy+temp2/temp*pace
        nposz=nposz-temp3/temp*pace
        goto 55
        endif
       endif

       endfind=1
       inner=0
       endx=nposx
       endy=nposy
       endz=nposz
       goto 54

      else
       nposx=nposx+angle(1)*jump
       nposy=nposy+angle(2)*jump
       nposz=nposz+angle(3)*jump
      endif

  55  force=force+1
      if(force.gt.10000) then
c     write(*,*) 'this is a forced ending'
      jumpratio=jumpratio+1
      endfind=1
      inner=0
      endx=nposx
      endy=nposy
      endz=nposz
      goto 54
      endif

      enddo

  54  return

      end

      subroutine findangle(nposx,nposy,nposz,switch,
     . subdvx,subdvy,subdvz,angl,totv)

      double precision nposx,nposy,nposz,totv
      integer switch
      double precision angl(3)
      double precision subdvx(-1:1,-1:1,-1:1),subdvy(-1:1,-1:1,-1:1),
     . subdvz(-1:1,-1:1,-1:1)
      double precision sx,sy,sz
      integer ix,iy,iz

      double precision jumpeps,pace,err
      integer ndimx,ndimy,ndimz,nt
      double precision dx,dy,dz

      common/global1/ jumpeps,pace,err
      common/global2/ ndimx,ndimy,ndimz,nt
      common/global3/ dx,dy,dz

      ix=-1
      iy=0
      iz=0
      sx=nposx+0.5d0
      sy=nposy
      sz=nposz
      if(nposy.lt.0.d0) then
       sy=1.d0+nposy
       iy=-1
      endif
      if(nposz.lt.0.d0) then
       sz=1.d0+nposz
       iz=-1
      endif
      angl(1)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvx(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvx(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvx(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvx(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvx(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvx(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvx(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvx(ix+1,iy+1,iz+1)
      angl(1)=angl(1)/dx*switch

      ix=0
      iy=-1
      iz=0
      sx=nposx
      sy=nposy+0.5d0
      sz=nposz
      if(nposx.lt.0.d0) then
       sx=1.d0+nposx
       ix=-1
      endif
      if(nposz.lt.0.d0) then
       sz=1.d0+nposz
       iz=-1
      endif
      angl(2)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvy(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvy(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvy(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvy(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvy(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvy(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvy(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvy(ix+1,iy+1,iz+1)
      angl(2)=angl(2)/dy*switch

      ix=0
      iy=0
      iz=-1
      sx=nposx
      sy=nposy
      sz=nposz+0.5d0
      if(nposx.lt.0.d0) then
       sx=1.d0+nposx
       ix=-1
      endif
      if(nposy.lt.0.d0) then
       sy=1.d0+nposy
       iy=-1
      endif
      angl(3)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvz(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvz(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvz(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvz(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvz(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvz(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvz(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvz(ix+1,iy+1,iz+1)
      angl(3)=angl(3)/dz*switch

      totv=1.d-30+dsqrt(angl(1)**2+angl(2)**2+angl(3)**2)
      angl(1)=angl(1)/totv
      angl(2)=angl(2)/totv
      angl(3)=angl(3)/totv
      totv=totv*dx

      return
      end

      subroutine endpoints(posx,posy,posz,endx,endy,endz,switch,
     . subdvx,subdvy,subdvz,subu1,subu2,subu3,subv1,subv2,subv3,
     . subw1,subw2,subw3,endfind,jumpratio,force0,up_down,deriv)

      double precision jumpratio
      double precision posx,posy,posz,nposx,nposy,nposz,endx,endy,endz
      double precision angle(3),temp
      double precision subdvx(-1:1,-1:1,-1:1),subdvy(-1:1,-1:1,-1:1),
     . subdvz(-1:1,-1:1,-1:1)

      double precision subu1(-1:0,-1:1,-1:1),subu2(-1:1,-1:0,-1:1),
     . subu3(-1:1,-1:1,-1:0)
      double precision subv1(-1:0,-1:1,-1:1),subv2(-1:1,-1:0,-1:1),
     . subv3(-1:1,-1:1,-1:0)
      double precision subw1(-1:0,-1:1,-1:1),subw2(-1:1,-1:0,-1:1),
     . subw3(-1:1,-1:1,-1:0)

      double precision u1,u2,u3,v1,v2,v3,w1,w2,w3

      double precision normx,normy,normz,normt,jump
      integer ix,iy,iz,endfind,force,switch
      double precision tu0,tv0,tw0,sx0,sy0,sz0,tu1,tv1,tw1,sx1,sy1,sz1
      double precision dfaids0,dfaids1,sfx,sfy,sfz
      double precision deriv
      integer up_down,force0

      double precision jumpeps,pace,err
      integer ndimx,ndimy,ndimz,nt
      double precision dx,dy,dz

      common/global1/ jumpeps,pace,err
      common/global2/ ndimx,ndimy,ndimz,nt
      common/global3/ dx,dy,dz

      inner=1
      force=0

      nposx=posx
      nposy=posy
      nposz=posz
      sx0=nposx
      sy0=nposy
      sz0=nposz

      endx=sx0
      endy=sy0
      endz=sz0

      ix=0
      iy=0
      iz=0
      sfx=sx0
      sfy=sy0
      sfz=sz0
      if(sx0.lt.0.d0) then
      ix=-1
      sfx=sx0+1
      endif
      if(sy0.lt.0.d0) then
      iy=-1
      sfy=sy0+1
      endif
      if(sz0.lt.0.d0) then
      iz=-1
      sfz=sz0+1
      endif
      tu0=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subdvx(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subdvx(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subdvx(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subdvx(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subdvx(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subdvx(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subdvx(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subdvx(ix+1,iy+1,iz+1)
      tv0=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subdvy(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subdvy(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subdvy(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subdvy(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subdvy(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subdvy(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subdvy(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subdvy(ix+1,iy+1,iz+1)
      tw0=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subdvz(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subdvz(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subdvz(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subdvz(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subdvz(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subdvz(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subdvz(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subdvz(ix+1,iy+1,iz+1)

      ix=-1
      iy=0
      iz=0
      sfx=sx0+0.5d0
      sfy=sy0
      sfz=sz0
      if(sy0.lt.0.d0) then
       sfy=1.d0+sy0
       iy=-1
      endif
      if(sz0.lt.0.d0) then
       sfz=1.d0+sz0
       iz=-1
      endif
      u1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subu1(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subu1(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subu1(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subu1(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subu1(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subu1(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subu1(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subu1(ix+1,iy+1,iz+1)
      v1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subv1(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subv1(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subv1(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subv1(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subv1(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subv1(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subv1(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subv1(ix+1,iy+1,iz+1)
      w1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subw1(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subw1(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subw1(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subw1(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subw1(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subw1(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subw1(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subw1(ix+1,iy+1,iz+1)

      ix=0
      iy=-1
      iz=0
      sfx=sx0
      sfy=sy0+0.5d0
      sfz=sz0
      if(sx0.lt.0.d0) then
       sfx=1.d0+sx0
       ix=-1
      endif
      if(sz0.lt.0.d0) then
       sfz=1.d0+sz0
       iz=-1
      endif
      u2=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subu2(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subu2(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subu2(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subu2(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subu2(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subu2(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subu2(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subu2(ix+1,iy+1,iz+1)
      v2=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subv2(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subv2(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subv2(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subv2(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subv2(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subv2(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subv2(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subv2(ix+1,iy+1,iz+1)
      w2=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subw2(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subw2(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subw2(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subw2(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subw2(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subw2(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subw2(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subw2(ix+1,iy+1,iz+1)

      ix=0
      iy=0
      iz=-1
      sfx=sx0
      sfy=sy0
      sfz=sz0+0.5d0
      if(sx0.lt.0.d0) then
       sfx=1.d0+sx0
       ix=-1
      endif
      if(sy0.lt.0.d0) then
       sfy=1.d0+sy0
       iy=-1
      endif
      u3=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subu3(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subu3(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subu3(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subu3(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subu3(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subu3(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subu3(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subu3(ix+1,iy+1,iz+1)
      v3=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subv3(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subv3(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subv3(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subv3(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subv3(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subv3(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subv3(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subv3(ix+1,iy+1,iz+1)
      w3=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subw3(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subw3(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subw3(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subw3(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subw3(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subw3(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subw3(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subw3(ix+1,iy+1,iz+1)

      dfaids0=tu0*u1*tu0+tu0*u2*tv0+tu0*u3*tw0
     . +tv0*v1*tu0+tv0*v2*tv0+tv0*v3*tw0
     . +tw0*w1*tu0+tw0*w2*tv0+tw0*w3*tw0
c this is used 'remeber' the value of last dfaids0 in the iteration
c process (up_down.gt.1). Because of the cuttoof error of computer,
c sometimes, the first deriv is not equal to the last one, and then
c error may be introduced.
      if(up_down.eq.1) then
      deriv=dfaids0
      endif

      do while(endfind.eq.0.and.(abs(endx).le.0.5d0).and.
     . (abs(endy).le.0.5d0).and.(abs(endz).le.0.5d0))

      call findangles(sx0,sy0,sz0,switch,subdvx,subdvy,subdvz,angle)
      sx1=sx0+angle(1)*pace
      sy1=sy0+angle(2)*pace
      sz1=sz0+angle(3)*pace
      force=force+1

      ix=0
      iy=0
      iz=0
      sfx=sx1
      sfy=sy1
      sfz=sz1
      if(sx1.lt.0.d0) then
      ix=-1
      sfx=sx1+1
      endif
      if(sy1.lt.0.d0) then
      iy=-1
      sfy=sy1+1
      endif
      if(sz1.lt.0.d0) then
      iz=-1
      sfz=sz1+1
      endif
      tu1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subdvx(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subdvx(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subdvx(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subdvx(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subdvx(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subdvx(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subdvx(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subdvx(ix+1,iy+1,iz+1)
      tv1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subdvy(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subdvy(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subdvy(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subdvy(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subdvy(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subdvy(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subdvy(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subdvy(ix+1,iy+1,iz+1)
      tw1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subdvz(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subdvz(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subdvz(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subdvz(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subdvz(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subdvz(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subdvz(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subdvz(ix+1,iy+1,iz+1)

      ix=-1
      iy=0
      iz=0
      sfx=sx1+0.5d0
      sfy=sy1
      sfz=sz1
      if(sy1.lt.0.d0) then
       sfy=1.d0+sy1
       iy=-1
      endif
      if(sz1.lt.0.d0) then
       sfz=1.d0+sz1
       iz=-1
      endif
      u1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subu1(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subu1(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subu1(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subu1(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subu1(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subu1(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subu1(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subu1(ix+1,iy+1,iz+1)
      v1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subv1(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subv1(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subv1(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subv1(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subv1(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subv1(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subv1(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subv1(ix+1,iy+1,iz+1)
      w1=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subw1(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subw1(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subw1(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subw1(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subw1(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subw1(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subw1(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subw1(ix+1,iy+1,iz+1)

      ix=0
      iy=-1
      iz=0
      sfx=sx1
      sfy=sy1+0.5d0
      sfz=sz1
      if(sx1.lt.0.d0) then
       sfx=1.d0+sx1
       ix=-1
      endif
      if(sz1.lt.0.d0) then
       sfz=1.d0+sz1
       iz=-1
      endif
      u2=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subu2(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subu2(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subu2(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subu2(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subu2(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subu2(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subu2(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subu2(ix+1,iy+1,iz+1)
      v2=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subv2(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subv2(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subv2(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subv2(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subv2(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subv2(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subv2(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subv2(ix+1,iy+1,iz+1)
      w2=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subw2(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subw2(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subw2(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subw2(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subw2(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subw2(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subw2(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subw2(ix+1,iy+1,iz+1)

      ix=0
      iy=0
      iz=-1
      sfx=sx1
      sfy=sy1
      sfz=sz1+0.5d0
      if(sx1.lt.0.d0) then
       sfx=1.d0+sx1
       ix=-1
      endif
      if(sy1.lt.0.d0) then
       sfy=1.d0+sy1
       iy=-1
      endif
      u3=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subu3(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subu3(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subu3(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subu3(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subu3(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subu3(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subu3(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subu3(ix+1,iy+1,iz+1)
      v3=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subv3(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subv3(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subv3(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subv3(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subv3(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subv3(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subv3(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subv3(ix+1,iy+1,iz+1)
      w3=(1.d0-sfx)*(1.d0-sfy)*(1.d0-sfz)*subw3(ix,iy,iz)+
     . sfx*(1.d0-sfy)*(1.d0-sfz)*subw3(ix+1,iy,iz)+
     . (1.d0-sfx)*sfy*(1.d0-sfz)*subw3(ix,iy+1,iz)+
     . (1.d0-sfx)*(1.d0-sfy)*sfz*subw3(ix,iy,iz+1)+
     . (1.d0-sfx)*sfy*sfz*subw3(ix,iy+1,iz+1)+
     . sfx*(1.d0-sfy)*sfz*subw3(ix+1,iy,iz+1)+
     . sfx*sfy*(1.d0-sfz)*subw3(ix+1,iy+1,iz)+
     . sfx*sfy*sfz*subw3(ix+1,iy+1,iz+1)

      dfaids1=tu1*u1*tu1+tu1*u2*tv1+tu1*u3*tw1
     . +tv1*v1*tu1+tv1*v2*tv1+tv1*v3*tw1
     . +tw1*w1*tu1+tw1*w2*tv1+tw1*w3*tw1

c this is from the consideration of the case of stagnation point
c which means stagnation points are special extremal points.
      if(tu0*tu1.le.0.d0.and.tv0*tv1.le.0.d0.and.tw0*tw1.le.0.d0)
     .  endfind=1
      if(deriv*dfaids1.lt.0.d0) endfind=1
      if(force.gt.500) endfind=1

      sx0=sx1
      sy0=sy1
      sz0=sz1
      tu0=tu1
      tv0=tv1
      dfaids0=dfaids1
      deriv=dfaids0

      endx=sx0
      endy=sy0
      endz=sz0

      enddo

      end

      subroutine findangles(nposx,nposy,nposz,switch,
     . subdvx,subdvy,subdvz,angl)

      double precision nposx,nposy,nposz
      double precision angl(3)
      double precision subdvx(-1:1,-1:1,-1:1),subdvy(-1:1,-1:1,-1:1),
     . subdvz(-1:1,-1:1,-1:1)
      double precision sx,sy,sz,totv
      double precision valuek1,valuek2,ghost(3),tap
      integer ix,iy,iz,switch

      double precision jumpeps,pace,err
      integer ndimx,ndimy,ndimz,nt
      double precision dx,dy,dz

      common/global1/ jumpeps,pace,err
      common/global2/ ndimx,ndimy,ndimz,nt
      common/global3/ dx,dy,dz

      ix=0
      iy=0
      iz=0
      sx=nposx
      sy=nposy
      sz=nposz
      if(sx.lt.0.d0) then
      ix=-1
      sx=sx+1
      endif
      if(sy.lt.0.d0) then
      iy=-1
      sy=sy+1
      endif
      if(sz.lt.0.d0) then
      iz=-1
      sz=sz+1
      endif

      angl(1)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvx(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvx(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvx(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvx(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvx(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvx(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvx(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvx(ix+1,iy+1,iz+1)
      angl(1)=angl(1)/dx*switch

      angl(2)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvy(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvy(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvy(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvy(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvy(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvy(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvy(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvy(ix+1,iy+1,iz+1)
      angl(2)=angl(2)/dy*switch

      angl(3)=(1.d0-sx)*(1.d0-sy)*(1.d0-sz)*subdvz(ix,iy,iz)+
     . sx*(1.d0-sy)*(1.d0-sz)*subdvz(ix+1,iy,iz)+
     . (1.d0-sx)*sy*(1.d0-sz)*subdvz(ix,iy+1,iz)+
     . (1.d0-sx)*(1.d0-sy)*sz*subdvz(ix,iy,iz+1)+
     . (1.d0-sx)*sy*sz*subdvz(ix,iy+1,iz+1)+
     . sx*(1.d0-sy)*sz*subdvz(ix+1,iy,iz+1)+
     . sx*sy*(1.d0-sz)*subdvz(ix+1,iy+1,iz)+
     . sx*sy*sz*subdvz(ix+1,iy+1,iz+1)
      angl(3)=angl(3)/dz*switch

      totv=1.d-20+dsqrt(angl(1)**2+angl(2)**2+angl(3)**2)
      angl(1)=angl(1)/totv
      angl(2)=angl(2)/totv
      angl(3)=angl(3)/totv
      totv=totv*dx

      return
      end

      subroutine tellnews(ni,nj,nk,endx,endy,endz,posx,posy,posz)
      integer ni,nj,nk
      double precision endx,endy,endz,posx,posy,posz

      if(abs(endx).le.0.5.and.abs(endy).le.0.5.and.abs(endz).le.0.5)
     . then
      ni=ni
      nj=nj
      nk=nk
      posx=endx
      posy=endy
      posz=endz
      elseif(endx.gt.0.5.and.abs(endy).le.0.5.and.abs(endz).le.0.5)then
      ni=ni+1
      nj=nj
      nk=nk
      posx=endx-1.d0
      posy=endy
      posz=endz
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5)
     .  .and.(abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj
      nk=nk
      posx=endx+1.d0
      posy=endy
      posz=endz
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni
      nj=nj+1
      nk=nk
      posx=endx
      posy=endy-1.d0
      posz=endz
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni
      nj=nj-1
      nk=nk
      posx=endx
      posy=endy+1.d0
      posz=endz
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.
     . (abs(endz).le.0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz
      elseif((abs(endx).le.0.5).and.(abs(endy).le.0.5).and.
     . (endz.gt.0.5)) then
      ni=ni
      nj=nj
      nk=nk+1
      posx=endx
      posy=endy
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(abs(endy).le.0.5).and.(endz.gt.0.5))
     . then
      ni=ni+1
      nj=nj
      nk=nk+1
      posx=endx-1.d0
      posy=endy
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5).and.
     . (endz.gt.0.5))then
      ni=ni-1
      nj=nj
      nk=nk+1
      posx=endx+1.d0
      posy=endy
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.(endz.gt.0.5))
     . then
      ni=ni
      nj=nj+1
      nk=nk+1
      posx=endx
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.
     . (endz.gt.0.5)) then
      ni=ni
      nj=nj-1
      nk=nk+1
      posx=endx
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.(endz.gt.0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk+1
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.(endz.gt.0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk+1
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz-1.d0
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.(endz.gt.0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk+1
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz-1.d0
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.(endz.gt.0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk+1
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz-1.d0

      elseif((abs(endx).le.0.5).and.(abs(endy).le.0.5).and.
     . (endz.lt.-0.5)) then
      ni=ni
      nj=nj
      nk=nk-1
      posx=endx
      posy=endy
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(abs(endy).le.0.5).and.
     . (endz.lt.-0.5)) then
      ni=ni+1
      nj=nj
      nk=nk-1
      posx=endx-1.d0
      posy=endy
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(abs(endy).le.0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni-1
      nj=nj
      nk=nk-1
      posx=endx+1.d0
      posy=endy
      posz=endz+1.d0
      elseif((abs(endx).le.0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni
      nj=nj+1
      nk=nk-1
      posx=endx
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((abs(endx).le.0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5))
     .  then
      ni=ni
      nj=nj-1
      nk=nk-1
      posx=endx
      posy=endy+1.d0
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5)) then
      ni=ni+1
      nj=nj+1
      nk=nk-1
      posx=endx-1.d0
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(endy.gt.0.5).and.(endz.lt.-0.5)) then
      ni=ni-1
      nj=nj+1
      nk=nk-1
      posx=endx+1.d0
      posy=endy-1.d0
      posz=endz+1.d0
      elseif((endx.lt.-0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5)) then
      ni=ni-1
      nj=nj-1
      nk=nk-1
      posx=endx+1.d0
      posy=endy+1.d0
      posz=endz+1.d0
      elseif((endx.gt.0.5).and.(endy.lt.-0.5).and.(endz.lt.-0.5)) then
      ni=ni+1
      nj=nj-1
      nk=nk-1
      posx=endx-1.d0
      posy=endy+1.d0
      posz=endz+1.d0

      endif

      return
      end

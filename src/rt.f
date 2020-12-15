!--------------------------------------------------------------------
          subroutine int3d(i,j,k,nx,ny,nz,xp,yp,zp,xyzb,nxyz,nb,ilb,ile,ibcc,
     &       xm,ym,zm,ni,nj,nk,io,icount)
!--------------------------------------------------------------------
c            implicit none
            INCLUDE 'common.h'
            INCLUDE 'mpif.h'

            integer i,j,k,nx,ny,nz,nb,ie,n,ntr(6),ntrp,ntrq
            integer io,ilb,ile,icount
            integer im,ip,jm,jp,km,kp
            real xp(nx,ny),yp(nx,ny),zp(nz),w(7),wp
            real xyzb(3,3,nb),nxyz(3,nb)
            real xcc,xccm,xccp,ycc,yccm,yccp,zcc,zccp,zccm,xccq,yccq,zccq
!            real w1,w2,w3,w4,w5,w6,w7
            real di,dj,dk,dm,wmax,wmin
            real d,wfip,wfjp,wfim,wfjm,wsum,wwal
            real xm,ym,zm,xmt,ymt,zmt
            real dist1,dist2p,dist2m,distmp,distmt,distpt
            integer iflag,jflag,kflag,iflagt,nckint,invtri,ib
            integer ibcc(nx,ny,nz),ibv0,ibv1,ibv2,ibv3,ibv4,ibv5,ibv6,ibsum
            real xa,ya,za,xb,yb,zb,xc,yc,zc
            real minxaxbxc,minyaybyc,minzazbzc
            real maxxaxbxc,maxyaybyc,maxzazbzc
            real ycen,zcen,omegat,erre,theta,pig,nerre,ntheta,grar,grarn
            real vyi,vzi,vri,vthi
            real vyj,vzj,vrj,vthj
            real vyk,vzk,vrk,vthk

            real ni,nj,nk
!
            INTEGER i1,j1,k1

            i1=4
            j1=2
            k1=12
!
            ntr = 0
            ntrp = 0
            ntrq = 0
            w=0.
            wmax=0.
            wmin=0.
            pig=2.*asin(1.)
!
            im = max(1 ,i-1)
            ip = min(nx,i+1)
            jm = max(1 ,j-1)
            jp = min(ny,j+1)
            km = max(1 ,k-1)
            kp = min(nz,k+1)
!            if (i.eq.1 ) im = nx  
!            if (i.eq.nx) ip = 1
!
            ibv0  = ibcc (i ,j ,k )
            ibv1  = ibcc (im,j ,k )
            ibv2  = ibcc (ip,j ,k )
            ibv3  = ibcc (i ,jm,k )
            ibv4  = ibcc (i ,jp,k )
            ibv5  = ibcc (i ,j ,km)
            ibv6  = ibcc (i ,j ,kp)
            ibsum = ibv1+ibv2+ibv3+ibv4+ibv5+ibv6
!
c            if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c              write(6,*) 'ibv0=',ibv0
c              write(6,*) 'ibv1=',ibv1
c              write(6,*) 'ibv2=',ibv2
c              write(6,*) 'ibv3=',ibv3
c              write(6,*) 'ibv4=',ibv4
c              write(6,*) 'ibv5=',ibv5
c              write(6,*) 'ibv6=',ibv6
c            endif
            if (ibv0 .ne. 1) write(8,*)'WARNING IBV0 ',ibv0
            if (abs(ibsum) .eq. 6) write(8,*)'WARNING IBVSUM ',ibsum,i,j,k
!
            iflag = 0
            jflag = 0
            kflag = 0
!
            grar = sqrt((xp(ip,jp)-xp(i,j))**2+(yp(ip,jp)-yp(i,j))**2+(zp(kp)-zp(k))**2)
!
            if (ibv0+ibv5==0.or.ibv0+ibv6==0) then
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '1. if (ibv0+ibv5==0.or.ibv0+ibv6==0) ',ilb,ile
c              endif
!            if (ibv5.ne.ibv6) then
               zcc  = zp(k)
               ycc  = yp(i,j)
               xcc  = xp(i,j)
               zccm = zp(km)
               zccp = zp(kp)
               do n=ilb,ile
!               do n=1,nb
                  xa = xyzb(1,1,n)
                  ya = xyzb(2,1,n)
                  za = xyzb(3,1,n)
                  xb = xyzb(1,2,n)
                  yb = xyzb(2,2,n)
                  zb = xyzb(3,2,n)
                  xc = xyzb(1,3,n)
                  yc = xyzb(2,3,n)
                  zc = xyzb(3,3,n)
!                  
                  minxaxbxc = min(xa,xb,xc)
                  minyaybyc = min(ya,yb,yc)
                  minzazbzc = min(za,zb,zc)
                  maxxaxbxc = max(xa,xb,xc)
                  maxyaybyc = max(ya,yb,yc)
                  maxzazbzc = max(za,zb,zc)
!
                  if( (ycc.ge.minyaybyc.and.ycc.le.maxyaybyc) .and.  
     &                  (xcc.ge.minxaxbxc.and.xcc.le.maxxaxbxc))then
                     call segtriint(xcc,ycc,zccm,xcc,ycc,zccp,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
!
                     if (iflagt.eq.1) then
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '1. found inters. ',xmt,ymt,zmt
c              endif

                        if(zmt.ge.zccm.and.zmt.le.zccp)then
                           kflag = 1
                           xm = xmt
                           ym = ymt
                           zm = zmt
                           ntr(5)=n
                           ntr(6)=n 
                           if (ibv5.eq.1) w(5) = (zcc-zm)/(zccm-zm)
                           if (ibv6.eq.1) w(6) = (zcc-zm)/(zccp-zm)
                           if (w(5).ne.0.0.and.w(6).ne.0.0) write(8,*) 'REALLY WARNING:kflag'
                        endif
                     endif
                  endif
               enddo
            endif
            if (ibv0+ibv3==0.or.ibv0+ibv4==0) then 
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '2. if (ibv0+ibv3==0.or.ibv0+ibv4==0)',ilb,ile
c              endif
!            if (ibv3.ne.ibv4) then 
               zcc  = zp(k)
               ycc  = yp(i,j)
               xcc  = xp(i,j)
               yccm = yp(i,jm)
               yccp = yp(i,jp)
               xccm = xp(i,jm)
               xccp = xp(i,jp)
               do n=ilb,ile
!               do n=1,nb
                  xa = xyzb(1,1,n)
                  ya = xyzb(2,1,n)
                  za = xyzb(3,1,n)
                  xb = xyzb(1,2,n)
                  yb = xyzb(2,2,n)
                  zb = xyzb(3,2,n)
                  xc = xyzb(1,3,n)
                  yc = xyzb(2,3,n)
                  zc = xyzb(3,3,n)
!                  
                  minxaxbxc = min(xa,xb,xc)
                  minyaybyc = min(ya,yb,yc)
                  minzazbzc = min(za,zb,zc)
                  maxxaxbxc = max(xa,xb,xc)
                  maxyaybyc = max(ya,yb,yc)
                  maxzazbzc = max(za,zb,zc)
!
                  if  (zcc.ge.minzazbzc.and.zcc.le.maxzazbzc) then
                     call segtriint(xccm,yccm,zcc,xccp,yccp,zcc,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
!
                     if (iflagt.eq.1) then
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '2. found inters. ',xmt,ymt,zmt
c              endif

                        distmp=sqrt((yccm-yccp)**2.+(xccm-xccp)**2.)
                        distmt=sqrt((yccm-ymt )**2.+(xccm-xmt )**2.)
                        distpt=sqrt((yccp-ymt )**2.+(xccp-xmt )**2.)
                        if(1.001*distmp.ge.(distmt+distpt))then
                           jflag = 1
                           xm = xmt
                           ym = ymt
                           zm = zmt
                           ntr(3)=n
                           ntr(4)=n 
                           dist1 =sqrt((ycc -ym)**2.+(xcc -xm)**2.)
                           dist2p=sqrt((yccp-ym)**2.+(xccp-xm)**2.)
                           dist2m=sqrt((yccm-ym)**2.+(xccm-xm)**2.)
                           if (ibv3.eq.1) w(3) = dist1/dist2m
                           if (ibv4.eq.1) w(4) = dist1/dist2p

!                           if(all(w(3:4)==0).and.ibv0+ibv3==0.and.ibv0+ibv4==0) then 
!                              w(4) = dist1/dist2m
!                              w(3) = dist1/dist2p
!                              write(*,*) w(3),w(4)
!                           endif

                           if (w(3).ne.0.0.and.w(4).ne.0.0) write(8,*) 'REALLY WARNING:jflag'
                        endif
                     endif
                  endif
               enddo
            endif
            if (ibv0+ibv1==0.or.ibv1+ibv2==0) then
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '3. if (ibv0+ibv1==0.or.ibv1+ibv2==0)',ilb,ile
c              endif
!            if (ibv1.ne.ibv2) then
               zcc  = zp(k)
               ycc  = yp(i,j)
               xcc  = xp(i,j)
               yccm = yp(im,j)
               yccp = yp(ip,j)
               xccm = xp(im,j)
               xccp = xp(ip,j)
               do n=ilb,ile
!               do n=1,nb
                  xa = xyzb(1,1,n)
                  ya = xyzb(2,1,n)
                  za = xyzb(3,1,n)
                  xb = xyzb(1,2,n)
                  yb = xyzb(2,2,n)
                  zb = xyzb(3,2,n)
                  xc = xyzb(1,3,n)
                  yc = xyzb(2,3,n)
                  zc = xyzb(3,3,n)
!                  
                  minxaxbxc = min(xa,xb,xc)
                  minyaybyc = min(ya,yb,yc)
                  minzazbzc = min(za,zb,zc)
                  maxxaxbxc = max(xa,xb,xc)
                  maxyaybyc = max(ya,yb,yc)
                  maxzazbzc = max(za,zb,zc)
!
                  if  (zcc.ge.minzazbzc.and.zcc.le.maxzazbzc) then
                     call segtriint(xccm,yccm,zcc,xccp,yccp,zcc,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
!
                     if (iflagt.eq.1) then
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '3. found inters. ',xmt,ymt,zmt
c              endif

                        distmp=sqrt((yccm-yccp)**2.+(xccm-xccp)**2.)
                        distmt=sqrt((yccm-ymt )**2.+(xccm-xmt )**2.)
                        distpt=sqrt((yccp-ymt )**2.+(xccp-xmt )**2.)
!                       if(i==i1.AND.j==j1.AND.k==k1) then
!                         WRITE(6,'(4(I4,1X),20(F14.8,1X))') I,J,K,IFLAGT 
!    &                      ,XCCM,YCCM,ZCC,XCCP,YCCP,ZCC,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XMT,YMT,ZMT
!                         WRITE(6,'(3(F14.8,1X))') distmp,distmt,distpt
!                       endif
                        if(1.001*distmp.ge.(distmt+distpt))then
                           iflag = 1
                           xm = xmt
                           ym = ymt
                           zm = zmt
                           ntr(1)=n
                           ntr(2)=n 
                           dist1 =sqrt((ycc -ym)**2.+(xcc -xm)**2.)
                           dist2p=sqrt((yccp-ym)**2.+(xccp-xm)**2.)
                           dist2m=sqrt((yccm-ym)**2.+(xccm-xm)**2.)
                           if (ibv1.eq.1) w(1) = dist1/dist2m
                           if (ibv2.eq.1) w(2) = dist1/dist2p
                           if (w(1).ne.0.0.and.w(2).ne.0.0)
     &                 write(8,*) 'REALLY WARNING:iflag'
                        endif
                     endif
                  endif
               enddo
            endif
!
            wmax = maxval(w)
            wmin = minval(w)
!            if (wmax.gt.1.0) write(*,*) 'REALLY WARNING:w>1.0'
!            if (wmax.eq.0.0) write(*,*) 'REALLY WARNING:w=0.0'
!            if (wmin.lt.0.0) write(*,*) 'REALLY WARNING:w<0.0'

            ie=0
            wp=1.0
            do n=1,6
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) 'n=',n,' ntr(n)=',ntr(n),',w(n)=',w(n)
c              endif
               if(w(n).gt.0.000001.and.w(n).lt.wp) then
                  wp=w(n)
                  ie=n
                  ntrp=ntr(n)
               endif
            enddo
!            if(ie.eq.0) then 
!               write(*,*)'error'
!               write(*,*) ibv0,ibv1,ibv2,ibv3,ibv4,ibv5,ibv6
!               write(*,*) wmax,wmin
!            endif


!            if(ie.ne.0)write(*,*)'ok'

!-----------NORMAL DISTANCE------------------------
!
!            goto 825
!
            invtri=0
            ntrq=ntrp
!            WRITE(6,*) NTRQ,'/',NB
            IF(NTRQ>NB.OR.NTRQ<1) THEN
              icount = 1
!              WRITE(6,*) NTRQ,'/',NB
!              WRITE(6,*) I,J,K
              GRAR = 100000.
              DO IB=1,NB
                DO N=1,3
                xa = xyzb(1,N,IB)
                ya = xyzb(2,N,IB)
                za = xyzb(3,N,IB)
                grarn = sqrt((xa-xp(i,j))**2.+(ya-yp(i,j))**2.+(za-zp(k))**2.)
                IF(grarn<grar) THEN
                  grar=grarn
                  xm = xa
                  ym = ya
                  zm = za
                  ni = (xp(i,j)-xa)/grar
                  nj = (yp(i,j)-ya)/grar
                  nk = (zp(k)-za)/grar
                ENDIF
               ENDDO
              ENDDO
!              WRITE(6,*) i,j,k,xp(i,j),yp(i,j),zp(k),xm,ym,zm,ni,nj,nk
              GO TO 825
            ENDIF
       
            zcc  = zp(k)
            ycc  = yp(i,j)
            xcc  = xp(i,j)
100         xccq=xcc-nxyz(1,ntrq)*grar*io
            yccq=ycc-nxyz(2,ntrq)*grar*io
            zccq=zcc-nxyz(3,ntrq)*grar*io
            nckint=0
            do n=ilb,ile
!            do n=1,nb
                  xa = xyzb(1,1,n)
                  ya = xyzb(2,1,n)
                  za = xyzb(3,1,n)
                  xb = xyzb(1,2,n)
                  yb = xyzb(2,2,n)
                  zb = xyzb(3,2,n)
                  xc = xyzb(1,3,n)
                  yc = xyzb(2,3,n)
                  zc = xyzb(3,3,n)
! 
               call segtriint(xcc,ycc,zcc,xccq,yccq,zccq,xa,ya,za,xb,yb,zb,xc,yc,zc,xmt,ymt,zmt,iflagt)
           !    write(*,*)iflagt
!
               if (iflagt.eq.1) then
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '4. found inters. ',xmt,ymt,zmt,n
c              endif

      !            write(*,*)'found intersection',n
                        distmp=sqrt((ycc-yccq)**2.+(xcc-xccq)**2.)
                        distmt=sqrt((ycc-ymt )**2.+(xcc-xmt )**2.)
                        distpt=sqrt((yccq-ymt )**2.+(xccq-xmt )**2.)
                        if(1.001*distmp.ge.(distmt+distpt))then
!                  if(xmt.ge.(min(xcc,xccq)).and.xmt.le.(max(xcc,xccq)))then
!                     if(ymt.ge.(min(ycc,yccq)).and.ymt.le.(max(ycc,yccq)))then
                        if(zmt.ge.(min(zcc,zccq)).and.zmt.le.(max(zcc,zccq)))then
                           nckint=nckint+1
                           xm = xmt
                           ym = ymt
                           zm = zmt
             !              write(*,*)'right intersection',n
                           if(n.ne.ntrq)then 
                  !            write(*,*)'diff trgl.ind',i,j,k,ntrq,n
                  !            write(*,*)'diff trgl.use',nxyz(ntrq,1),nxyz(ntrq,2),nxyz(ntrq,3)
                  !            write(*,*)'diff trgl.fou',nxyz(n,1),nxyz(n,2),nxyz(n,3)
                              ntrq=n
                              if(invtri.eq.3)goto 200
                              invtri=invtri+1
                              goto 100
200                           continue
!                              write(*,*)'skip invtri'
!                           endif
                        endif
                     endif
                  endif
               endif
            enddo
!
            if(nckint.ne.1)then 
!               write(*,*)'error intersection',nckint
               ntrq=ntrp
!
               if(ie.eq.1.or.ie.eq.2)then
               zcc  = zp(k)
               ycc  = yp(i,j)
               xcc  = xp(i,j)
               yccm = yp(im,j)
               yccp = yp(ip,j)
               xccm = xp(im,j)
               xccp = xp(ip,j)
                  xa = xyzb(1,1,ntrq)
                  ya = xyzb(2,1,ntrq)
                  za = xyzb(3,1,ntrq)
                  xb = xyzb(1,2,ntrq)
                  yb = xyzb(2,2,ntrq)
                  zb = xyzb(3,2,ntrq)
                  xc = xyzb(1,3,ntrq)
                  yc = xyzb(2,3,ntrq)
                  zc = xyzb(3,3,ntrq)
               call segtriint(xccm,yccm,zcc,xccp,yccp,zcc,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflagt)
               if (iflagt.ne.1)  pause 'rt.f:iflagt.ne.1'
               endif
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '5. found inters. ',xm,ym,zm,n
c              endif

               if(ie.eq.3.or.ie.eq.4)then
               zcc  = zp(k)
               ycc  = yp(i,j)
               xcc  = xp(i,j)
               yccm = yp(i,jm)
               yccp = yp(i,jp)
               xccm = xp(i,jm)
               xccp = xp(i,jp)
                  xa = xyzb(1,1,ntrq)
                  ya = xyzb(2,1,ntrq)
                  za = xyzb(3,1,ntrq)
                  xb = xyzb(1,2,ntrq)
                  yb = xyzb(2,2,ntrq)
                  zb = xyzb(3,2,ntrq)
                  xc = xyzb(1,3,ntrq)
                  yc = xyzb(2,3,ntrq)
                  zc = xyzb(3,3,ntrq)
               call segtriint(xccm,yccm,zcc,xccp,yccp,zcc,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflagt)
               if (iflagt.ne.1) pause 'rt.f:iflagt.ne.1'
               endif
c              if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '6. found inters. ',xm,ym,zm,n
c              endif
               if(ie.eq.5.or.ie.eq.6)then
               zcc  = zp(k)
               ycc  = yp(i,j)
               xcc  = xp(i,j)
               zccm = zp(km)
               zccp = zp(kp)
                  xa = xyzb(1,1,ntrq)
                  ya = xyzb(2,1,ntrq)
                  za = xyzb(3,1,ntrq)
                  xb = xyzb(1,2,ntrq)
                  yb = xyzb(2,2,ntrq)
                  zb = xyzb(3,2,ntrq)
                  xc = xyzb(1,3,ntrq)
                  yc = xyzb(2,3,ntrq)
                  zc = xyzb(3,3,ntrq)
               call segtriint(xcc,ycc,zccm,xcc,ycc,zccp,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflagt)
               if (iflagt.ne.1)  pause
               endif
c             if(i==i1.AND.j==j1.AND.k+MYRANK*(NZ-2)==k1) then
c                write(6,*) '7. found inters. ',xm,ym,zm,n
c              endif
            endif
!
!            erre=sqrt(xm**2.+ym**2.)
!            theta=atan(ym/xm)
!            if (xm.lt.0.0) theta=theta+pig
!            if (xm.eq.0.0) then
!               if (ym.eq.0.0) theta=0.0
!               if (ym.gt.0.0) theta=+pig/2.0
!               if (ym.lt.0.0) theta=-pig/2.0
!            endif
!            if(theta.lt.0)theta=theta+2.*pig
!            ntheta = nxyz(2,ntrq)*cos(theta)-nxyz(1,ntrq)*sin(theta)
!            nerre  = nxyz(2,ntrq)*sin(theta)+nxyz(1,ntrq)*cos(theta)
!            write(87,257)i,j,k,theta,erre,zm,ntheta,nerre,nxyz(3,ntrq)
!257         format(3(2x,i5),6(2x,e14.7))
        !    write(*,*)'end int3d'

            ni = nxyz(1,ntrq)
            nj = nxyz(2,ntrq)
            nk = nxyz(3,ntrq)
!
825         continue
            return
          end subroutine int3d!

c---------------------------------------------------------------------
C
C     PURPOSE: Find points that are outside or inside the immersed object
C
C---------------------------------------------------------------------
C
        SUBROUTINE TAG3D(X,Y,Z,NI,NJ,NK
     &                  ,XYZB,XYZCB,NXYZB,NB,ILB,ILE,S
     &                  ,IBS,IBD,ICOM)

        INCLUDE 'common.h'
        INCLUDE 'immersed.h'
        INCLUDE 'mpif.h'

        INTEGER NCR
        PARAMETER (NCR = 10)   !!!!!!

        INTEGER NI,NJ,NK,NB,ILB,ILE,IBD,ICOM
        INTEGER IBS(NI,NJ,NK)
        REAL S(NI), DS
        REAL X(NI,NJ),Y(NI,NJ),Z(NK)
        REAL XYZB(3,3,NB),XYZCB(3,NB),NXYZB(3,NB)
c     
c...  Local arrays
        INTEGER ICROSS,IFLAG,I1,I2,ISAVE,N,NINS,NEXT,NINSG,ISGN,II,IMN,IMX
        INTEGER I,J,K,IC,JN,ISGNT,NNMAX,NNMAX2,NOUT,NQ,NN,NN2,IN,NC,NDP,N1,N2
        INTEGER JJ1,JJ2,J1,J2,K1,K2,NDJ,IMIN,IMAX,IBMAXG,NCROSSMAXG
        INTEGER FIRST,LAST,MYFIRST,MYLAST
        INTEGER STATUS(MPI_STATUS_SIZE)

        REAL    QUERYXY(2)!,XYB(2,NB)
        REAL    DOTPROD
        REAL    SGN,DIST,SGN1

        INTEGER NCROSS(NI,NJ),NCROSSG(NI,NJ),NCRPRV(NJ)
        REAL    ZCROSS(NCR,NJ),ZCROSSO(NCR,NI,NJ)
        REAL    ZCROSSG(NCR,NJ),ZCROSSOG(NCR,NI,NJ)

        REAL    xcc,ycc,zcc,xm,ym,zm,ztmin,rzmin,rzmax,tol,minz,maxz,minzg,maxzg
        REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,rayzmin,rayzmax
        REAL    pcc(3),va(3),vb(3),vc(3),qvec(3),rvec(3),pint(3)

        INTEGER, DIMENSION(:), ALLOCATABLE :: INN,INN2,NNTRI

        INTEGER COUNTS,RATE,NCLOCKS
        REAL    CLOCKTEMP1,CLOCKTEMP2,CLOCKTEMP3,TCLOCK

        REAL, DIMENSION(:), ALLOCATABLE :: CLOCK,CLOCKG,CLOCKGMIN,CLOCKGMAX
c
c..... Functions
        REAL    anglerad

        NCLOCKS = 14
        ALLOCATE(CLOCK(NCLOCKS),CLOCKG(NCLOCKS),CLOCKGMIN(NCLOCKS),CLOCKGMAX(NCLOCKS))
        CLOCK = 0.0
        CLOCKG = 0.0
        CLOCKGMIN = 0.0
        CLOCKGMAX = 0.0

        rayzmin = min(zmin-zlen,-1.0e8)
        rayzmax = max(zmax+zlen,1.0e8)

        CLOCK(1) = tclock()
        
        NNMAX = 1000
        NNMAX2 = -1
        N1 = 0
        N2 = 0

        ALLOCATE(INN(NNMAX))
        ALLOCATE(NNTRI(NNMAX))
        NNTRI = 0

        IF(IMBDCMP==0) THEN
c
c in the case the immersed-boundary is not decomposed
c   
          CLOCK(2) = tclock()
          IF(ICOM==0) THEN
c
c initialization call (ICOM=0)
c
            CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 0, MYRANK, DS)
          ENDIF
          CLOCK(2) = tclock()-CLOCK(2)
        
          IBS = 1

          if(idomy==0) then

            IF(MYSIZE<=NJ-2) THEN
              JJ1 = 2
              JJ2 = (NJ-2)/MYSIZE
              J1=JJ1 + JJ2*MYRANK
              J2=J1 + JJ2 -1
              NDJ = J2-J1+1
              NQ = NI*NDJ
              IMN = IBMIN(IBD)
              IMX = IBMAX(IBD)
            ELSE
              JJ1 = 2
              JJ2 = MYSIZE/(NJ-2)
c              JJ2 = (NJ-2)/MYSIZE
              J1 = JJ1 + REAL(NJ-2)/REAL(MYSIZE)*MYRANK
              J2 = J1 + REAL(NJ-2)/REAL(MYSIZE)

              CALL MPI_ALLREDUCE(IBMAX(IBD),IMX,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
              IMN = MOD(MYRANK,JJ2)*REAL(IMX-IBMIN(IBD)+1)/REAL(JJ2)+IBMIN(IBD)
              IMX = IMN+REAL(IMX-IBMIN(IBD)+1)/REAL(JJ2)-1              
            ENDIF

          else
            J1 = 2
            J2 = NJ-1
          endif

          nins=0
          next=0
          
          IMIN = IBMIN(ibd)
          CALL MPI_ALLREDUCE(IBMAX(IBD),IMAX,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

c          DO I=IMIN,IBMAX(IBD)
          NCROSS = 0
          ZCROSSO = 0.0

          DO I=IMN,IMX

            ZCROSS = 0.0
            
            CLOCKTEMP1 = tclock()

c
c DS is a function of the X position
c
            DS = 1.0*S(I)**2.            
c            DS = 1.0*MAXVAL(S)**2.            
c            write(6,*) i,s(i)
c            DS = 2.0*S(I)**2.            
          
            DO J = J1,J2

              xcc = X(I,J)
              ycc = Y(I,J)
              QUERYXY(1) = XCC
              QUERYXY(2) = YCC
              ICROSS = 0
              CLOCKTEMP2 = tclock()
c
c committee4 finds the centers of the triangles closest to the grid point
c
              CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNTR2D, 2, 1, MYRANK, DS)
              CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2

              IF(NN>NNMAX2) NNMAX2=NN

              IF(NN>NNTR2D) THEN
                N1 = N1+1
              ELSE
                N2 = N2+1
              ENDIF

              IF(NN>NNTR2D) THEN                     
                IF(NN>NNMAX) THEN
                  clocktemp2 = tclock()
                  DEALLOCATE(INN,NNTRI)
                  NNMAX = NN
                  ALLOCATE(INN(NNMAX),NNTRI(NNMAX))
                  clock(12) = clock(12) + tclock() - clocktemp2
                ENDIF
                CLOCKTEMP2 = tclock()
                CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NN, 2, 1, MYRANK, DS)
                CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2
              ENDIF

              clocktemp3 = tclock()
              DO IN=1,MIN(NN,NNMAX)     
                N = INN(IN)
                xa = xyzb(1,1,n)
                ya = xyzb(2,1,n)
                za = xyzb(3,1,n)
                xb = xyzb(1,2,n)
                yb = xyzb(2,2,n)
                zb = xyzb(3,2,n)
                xc = xyzb(1,3,n)
                yc = xyzb(2,3,n)
                zc = xyzb(3,3,n)
                IF ((ycc.ge.min(ya,yb,yc).and.ycc.le.max(ya,yb,yc)) .and.  
     &             (xcc.ge.min(xa,xb,xc).and.xcc.le.max(xa,xb,xc)) ) THEN
c
c segtriint looks for an intersection of the triangle N with the control
c ray with equations X=XCC and Y=YCC 
c
                  CALL segtriint(xcc,ycc,rayzmin,xcc,ycc,rayzmax
     &               ,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
                  IF (iflag.eq.1) then
                    NNTRI(IN) = NNTRI(IN)+1
                    ICROSS = ICROSS+1
                    IF(ICROSS>NCR) THEN
                      WRITE(6,*) '1. WARNING MYRANK=',MYRANK,J,I,ICROSS,XCC,YCC,icross,IBD
c                      do ii=1,ncr
c                        write(6,*) ii,zcross(ii,j)
c                      enddo
                    ENDIF
                    ZCROSS(ICROSS,J) = ZM

                  ENDIF

                ENDIF
              ENDDO
              clock(9) = clock(9) + tclock() - clocktemp3

              NCROSS(I,J) = ICROSS
              IF(NCROSS(I,J)>NCR) THEN
                WRITE(6,*) '2. WARNING MYRANK=',MYRANK,J,I,NCROSS(I,J),IBD
              ENDIF

              clocktemp3 = tclock()
              !Sort z-coord. crossings in ascending order
              DO I1=1,NCROSS(I,J)
                ztmin = 1.0e10
                do i2=1,NCROSS(I,J)
                  if (ZCROSS(i2,J).le.ztmin) then
                    ztmin = ZCROSS(i2,J)
                    isave = i2;
                  endif
                enddo
c
c in ZCROSSO the Z coordinates of the intersection points are organized
c in ascending order
c
                ZCROSSO(i1,i,J) = ztmin
                ZCROSS (isave,J) = 2.0e10
              ENDDO
              clock(10) = clock(10) + tclock() - clocktemp3

              clocktemp3 = tclock()
              !Remove duplicate z-coord. crossings
                zcross(:,j) = zcrosso(:,i,j)
                nc=0
                DO i1=1,ncross(i,j)
                  ndp=0
                  DO i2=1,ncross(i,j)
!!!!!!                    if(zcross(i2,j)==zcross(i1,j).AND.i1>i2) then
                    if(abs(zcross(i2,j)-zcross(i1,j)).lt.1.e-06.AND.i1>i2) then
                      ndp=ndp+1
                    endif
                  ENDDO

!!!!!!                  if(ndp<=1) then
                  if(ndp==0) then
                    nc=nc+1
                    zcrosso(nc,i,j) = zcross(i1,j)
                  endif
                ENDDO
              clock(11) = clock(11) + tclock() - clocktemp3
c
c actual number of intersection points
c
              ncross(i,j) = nc

c
c warning if the number of intersection points is odd
c
              IF(mod(ncross(i,j),2)/=0) THEN
                write(6,*) 'WARNING:',i,j,ncross(i,j),zcrosso(1:ncross(i,j),i,j)
              ENDIF

            ENDDO

            CLOCK(4) = CLOCK(4) + tclock() - CLOCKTEMP1
          ENDDO

          if(mysize>1 .AND. idomy==0) then
            clocktemp1 = tclock()
c            if(mysize>nj-2) then
            CALL MPI_ALLREDUCE(MAXVAL(NCROSS(IMIN:IMAX,2:NJ-1)),NCROSSMAXG,1,MPI_INTEGER
     &           ,MPI_MAX,MPI_COMM_EDDY,IERR)
              CALL MPI_ALLREDUCE(NCROSS(IMIN:IMAX,2:NJ-1),NCROSSG(IMIN:IMAX,2:NJ-1)
     &           ,(IMAX-IMIN+1)*(NJ-2),MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
              CALL MPI_ALLREDUCE(ZCROSSO(1:NCROSSMAXG,IMIN:IMAX,2:NJ-1),ZCROSSOG(1:NCROSSMAXG,IMIN:IMAX,2:NJ-1)
     &             ,NCROSSMAXG*(IMAX-IMIN+1)*(NJ-2),MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
c            else
c              CALL MPI_ALLREDUCE(NCROSS,NCROSSG,NI*NJ,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c              CALL MPI_ALLREDUCE(ZCROSSO(1:NCROSSMAXG,1:NI,1:NJ),ZCROSSOG(1:NCROSSMAXG,1:NI,1:NJ)
c     &              ,NCROSSMAXG*NI*NJ,MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
c            endif

            NCROSS = NCROSSG
            ZCROSSO = ZCROSSOG
            clock(13) = clock(13) + tclock() - clocktemp1
          endif

          CLOCKTEMP1 = tclock()
c
c if the number of the intersections upstream the grid point is odd
c IBS=-1 and the point is tagged as interior
c
          DO I=IBMIN(IBD),IBMAX(IBD)
            DO J=2,NJ-1
              DO ICROSS=1,NCROSS(I,J),2
                CALL LOCATE(z,nk,zcrosso(icross,i,j),k1)
                CALL LOCATE(z,nk,zcrosso(icross+1,i,j),k2)
                if(k2>k1) then
                  k1 = k1+1
                  ibs(i,j,k1:k2) = -1
                  nins = nins+(k2-k1+1)
                  next = next+(nk-2)-(k2-k1+1)
                endif
              ENDDO
            ENDDO

            CLOCK(6) = CLOCK(6) + tclock() - CLOCKTEMP1

c            call mpi_finalize(ierr)
c            stop

          ENDDO

          CLOCK(14) = tclock()
          if(idomy==0) then
            IBS(IMIN:IBMAX(IBD),1 ,KBMIN(IBD):KBMAX(IBD)) = IBS(IMIN:IBMAX(IBD),NJ-1,KBMIN(IBD):KBMAX(IBD))
            IBS(IMIN:IBMAX(IBD),NJ,KBMIN(IBD):KBMAX(IBD)) = IBS(IMIN:IBMAX(IBD), 2 ,KBMIN(IBD):KBMAX(IBD))
          endif
c
c update of IBS for the ghost layers between processors
c
          call refreshflag(ibs,ni,nj,nk)

          CLOCK(14) = tclock() - clock(14)

          CLOCK(7) = tclock()
c
c deallocation of the variables for COMMITTEE4
c
          IF(ICOM==2) CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 2, MYRANK, DS)
          clock(7) = tclock()-clock(7)

        ELSE
c
c in the case the immersed-boundary is decomposed
c
c          MINZ = MINVAL(XYZB(3,:,:))
c          MAXZ = MAXVAL(XYZB(3,:,:))

c          IF(NB==0) MINZ=100000.0
c          CALL MPI_ALLREDUCE(MINZ,MINZG,1,MTYPE,MPI_MIN,MPI_COMM_EDDY,IERR)
c          IF(NB==0) MAXZ=-100000.0
c          CALL MPI_ALLREDUCE(MAXZ,MAXZG,1,MTYPE,MPI_MAX,MPI_COMM_EDDY,IERR)

c          write(6,*) 'myrank=',myrank,'minz=',minz,'minzg=',minzg
c     &         ,',maxz=',maxz,'maxzg=',maxzg,'z:',z(1),z(nk-1)
c          call mpi_finalize(ierr)
c          stop

          FIRST=MYSIZE+1
          LAST=-1
          IF(NB>0) THEN
            FIRST=MYRANK
            LAST=MYRANK           
          ENDIF

c          IF(MINZG>Z(1) .AND. MINZG<=Z(NK-1)) THEN
c            FIRST=MYRANK
c          ELSE
c            FIRST=-1
c          ENDIF

          CALL MPI_ALLREDUCE(FIRST,MYFIRST,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)

c          IF(MAXZG>Z(1) .AND. MAXZG<=Z(NK-1)) THEN
c            LAST=MYRANK
c          ELSE
c            LAST=-1
c          ENDIF

          CALL MPI_ALLREDUCE(LAST,MYLAST,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

c          write(6,*) 'myrank=',myrank,myfirst,mylast
c          call mpi_finalize(ierr)
c          stop

          CALL MPI_ALLREDUCE(IBMAX(IBD),IMAX,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)
c          CALL MPI_ALLREDUCE(IBMIN(IBD),IMIN,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)

          
          eps=1.0e-11
          rzmin=z(1)+eps
          rzmax=z(nk-1)
          if(myrank==myfirst) rzmin=-100.0

          CLOCK(2) = tclock()
          IF(ICOM==0 .AND. NB>0) CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 0, MYRANK, DS)
          CLOCK(2) = tclock()-CLOCK(2)
        
          IBS = 1

          J1 = 2
          J2 = NJ-1

          nins=0
          next=0

          IMIN = IBMIN(ibd)

          IF(NB>0) THEN

            NCROSS = 0
            ZCROSSO = 0.0
            DO I=IMIN,IMAX
              ZCROSS = 0.0
            
              CLOCKTEMP1 = tclock()

              DS = 0.5*S(I)**2.            
          
              DO J = J1,J2

                xcc = X(I,J)
                ycc = Y(I,J)
                QUERYXY(1) = XCC
                QUERYXY(2) = YCC
                ICROSS = 0
                CLOCKTEMP2 = tclock()
                CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 1, MYRANK, DS)
                CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2

                IF(NN>NNMAX) THEN 
                  clocktemp2 = tclock()
                  DEALLOCATE(INN)
                  NNMAX = NN
                  ALLOCATE(INN(NNMAX))
                  clock(12) = clock(12) + tclock() - clocktemp2
                  CLOCKTEMP2 = tclock()
                  CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 1, MYRANK, DS)
                  CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2
                ENDIF

                clocktemp3 = tclock()
                DO IN=1,MIN(NN,NNMAX)
                  N = INN(IN)
                  xa = xyzb(1,1,n)
                  ya = xyzb(2,1,n)
                  za = xyzb(3,1,n)
                  xb = xyzb(1,2,n)
                  yb = xyzb(2,2,n)
                  zb = xyzb(3,2,n)
                  xc = xyzb(1,3,n)
                  yc = xyzb(2,3,n)
                  zc = xyzb(3,3,n)
                  IF ((ycc.ge.min(ya,yb,yc).and.ycc.le.max(ya,yb,yc)) .and.  
     &              (xcc.ge.min(xa,xb,xc).and.xcc.le.max(xa,xb,xc)) ) THEN
                    CALL segtriint(xcc,ycc,rzmin,xcc,ycc,rzmax
     &               ,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
                    IF (iflag.eq.1) then
                      ICROSS = ICROSS+1
                      IF(ICROSS>NCR) THEN
                        WRITE(6,*) '1. WARNING MYRANK=',MYRANK,J,I,ICROSS,XCC,YCC,IBD
c                        do ii=1,ncr
c                          write(6,*) ii,zcross(ii,j)
c                        enddo
                      ENDIF
                      ZCROSS(ICROSS,J) = ZM

                    ENDIF

                  ENDIF
                ENDDO
                clock(9) = clock(9) + tclock() - clocktemp3

                NCROSS(I,J) = ICROSS
                IF(NCROSS(I,J)>NCR) THEN
                  WRITE(6,*) '2. WARNING MYRANK=',MYRANK,J,I,NCROSS(I,J),IBD
                ENDIF

                clocktemp3 = tclock()
                !Sort z-coord. crossings in ascending order
                DO I1=1,NCROSS(I,J)
                  ztmin = 1.0e10
                  do i2=1,NCROSS(I,J)
                    if (ZCROSS(i2,J).le.ztmin) then
                      ztmin = ZCROSS(i2,J)
                      isave = i2;
                    endif
                  enddo
                  ZCROSSO(i1,I,J) = ztmin
                  ZCROSS (isave,J) = 2.0e10
                ENDDO
                clock(10) = clock(10) + tclock() - clocktemp3

                clocktemp3 = tclock()
               !Remove duplicate z-coord. crossings
                zcross(:,j) = zcrosso(:,i,j)
                nc=0
                DO i1=1,ncross(i,j)
                  ndp=0
                  DO i2=1,ncross(i,j)
!!!!!!                    if(zcross(i2,j)==zcross(i1,j).AND.i1>i2) then
                    if(abs(zcross(i2,j)-zcross(i1,j)).lt.1.e-06.AND.i1>i2) then
                      ndp=ndp+1
                    endif
                  ENDDO

                  if(ndp==0) then
                    nc=nc+1
                    zcrosso(nc,i,j) = zcross(i1,j)
                  endif
                ENDDO
                clock(11) = clock(11) + tclock() - clocktemp3

                ncross(i,j) = nc
              ENDDO

              CLOCK(4) = CLOCK(4) + tclock() - CLOCKTEMP1

 
              CLOCKTEMP1 = tclock()

c              IF(MYRANK==MYFIRST .AND. MYFIRST/=MYLAST) THEN
c                CALL MPI_SEND(NCROSS,NY,MPI_INTEGER,JP,JP,MPI_COMM_EDDY,IERR)
c              ENDIF

              ncrprv=0

              DO J=2,NJ-1

                DO k=KBMIN(IBD),KBMAX(IBD)
                  zcc = z(k)
                  DO ICROSS=1,NCROSS(I,J)
                    IF (ZCROSSO(ICROSS,I,J) .LE. zcc) THEN
                      ibs(i,j,k) = -ibs(i,j,k) 
                    ENDIF
                  ENDDO

                  IF(ibs(i,j,k).eq.-1) nins=nins+1
                  IF(ibs(i,j,k).eq.1)  next=next+1
                ENDDO

              ENDDO

              IF(MYRANK>MYFIRST) THEN
                CALL MPI_RECV(NCRPRV,NJ,MPI_INTEGER,MYRANK-1,MYRANK-1,MPI_COMM_EDDY,STATUS,IERR)
              ENDIF

              IF(MYRANK<MYLAST) THEN
                CALL MPI_SEND(NCRPRV+NCROSS(I,:),NJ,MPI_INTEGER,MYRANK+1,MYRANK,MPI_COMM_EDDY,IERR)
              ENDIF

              IF(MYRANK>MYFIRST) THEN
                DO J=2,NJ-1
                  IF(MOD(NCRPRV(J),2)==1) THEN
c                  IBS(I,J,KBMIN(IBD):KBMAX(IBD))=-1
                    IBS(I,J,KBMIN(IBD):KBMAX(IBD))=-IBS(I,J,KBMIN(IBD):KBMAX(IBD))
                  ENDIF
                ENDDO
 
              ENDIF

            CLOCK(6) = CLOCK(6) + tclock() - CLOCKTEMP1

            ENDDO

          ENDIF

          CLOCK(14) = tclock()
          if(idomy==0) then
            IBS(IMIN:IBMAX(IBD),1,KBMIN(IBD):KBMAX(IBD))  = IBS(IMIN:IBMAX(IBD),NJ-1,KBMIN(IBD):KBMAX(IBD))
            IBS(IMIN:IBMAX(IBD),NJ,KBMIN(IBD):KBMAX(IBD)) = IBS(IMIN:IBMAX(IBD),2,KBMIN(IBD):KBMAX(IBD))
          endif
         
          call refreshflag(ibs,ni,nj,nk)

          CLOCK(14) = tclock() - clock(14)

          CLOCK(7) = tclock()
          IF(ICOM==2 .AND. NB>0) CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 2, MYRANK, DS)
          clock(7) = tclock()-clock(7)

       ENDIF

       CLOCK(1) = tclock() - CLOCK(1)


       IF(IOLVL>0) THEN
         CALL MPI_REDUCE(NINS,NINSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)

         IF(MYRANK.EQ.0) THEN
           OPEN(UNIT=16, FILE='stats_imb.dat', FORM='FORMATTED'
     &          ,POSITION='APPEND')
           N = ((NI-2)*(NJ-2)*(NK-2)*MYSIZE)
           WRITE(16,'(2(A,1X,I12),A,F6.2,A)') 'TAG3D: No. INSIDE='
     &          ,NINSG,'/',N,', (~',100.*REAL(NINSG)/REAL(N),'/100)'
           CLOSE(16)
         ENDIF
       ENDIF

       IF(IOCLOCK>0) THEN
         CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
         CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
         CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
          
         CLOCKG = CLOCKG/REAL(MYSIZE)


         IF(MYRANK.EQ.0) THEN
c           WRITE(6,*) 'N1=',N1,REAL(N1)/REAL(N1+N2)*100,', N2=',N2,REAL(n2)/REAL(N1+N2)*100
c           WRITE(6,*) 'MAX NN=',NNMAX2
c           DO I=1,NNMAX
c             IF(NNTRI(I)>0) THEN
c               WRITE(6,*) I,NNTRI(I)
c             ENDIF
c           ENDDO
           OPEN(UNIT=16, FILE='clock.dat',FORM='FORMATTED'
     &           ,POSITION='APPEND')
           WRITE(16,'(A)') '---- TAG3D: ----------------------------'
           WRITE(16,'(A)') '   Task/Time            Ave.              Max.      
     &           Min.'
           WRITE(16,'(A,3(5x,F8.4,5x))') 'Total               :',CLOCKG(1),CLOCKGMAX(1),CLOCKGMIN(1)
           i=3
           write(16,905) 'Build query array  :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=2
           write(16,905) 'Build ANN structure:'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=4
           write(16,905) 'Find z crossings   :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=8
           write(16,905) ' -Find NN          :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=9
           write(16,905) ' -Find intrsctns.  :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=10
           write(16,905) ' -Order z crossings:'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=11
           write(16,905) ' -Remove duplic.   :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=12
           write(16,905) ' -Dealloc. & alloc.:'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=13
           write(16,905) 'MPI_ALLreduces.    :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=6
           write(16,905) 'Classify points    :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=14
           write(16,905) 'Refresh ibs y-dir  :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=7
           write(16,905) 'Delete ANN struct. :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           WRITE(16,'(A)') '----------------------------------------'
           CLOSE(16)
                      
         ENDIF
       ENDIF 

       DEALLOCATE(CLOCK,CLOCKG,CLOCKGMIN,CLOCKGMAX)
       DEALLOCATE(INN)
       DEALLOCATE(NNTRI)

       RETURN

 905  format(A,3(1x,F8.4,A,F6.2,A,2x))

       END
C---------------------------------------------------------------------

!--------------------------------------------------------------------
          subroutine segtriint(x1,y1,z1,x2,y2,z2,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
!--------------------------------------------------------------------
            implicit none
            real x1,y1,z1,x2,y2,z2,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm
            real q(3),r(3),va(3),vb(3),vc(3),t,u,v
            integer iflag,code,itemp

!
!...redefine ray as origin plus lenght vector
!
            q(1) = x1
            q(2) = y1
            q(3) = z1
            r(1) = x2-x1
            r(2) = y2-y1
            r(3) = z2-z1
            va(1) = xa
            va(2) = ya
            va(3) = za
            vb(1) = xb
            vb(2) = yb
            vb(3) = zb
            vc(1) = xc
            vc(2) = yc
            vc(3) = zc
!
            ITEMP=0
c            IF(IFLAG==-1) THEN
c              ITEMP=1
c              WRITE(6,'(A,15(1X,F14.8))') '1. TEST SEGTRI',q(1),q(2),q(3)
c     &         ,r(1),r(2),r(3)
c     &         ,xa,xb,xc,ya,yb,yc,za,zb,zc
c            ENDIF
            iflag = 0.
            call intersect_triangle(q,r,va,vb,vc,t,u,v,code)
            IF(code.eq.1.AND.(t>1.0.OR.t<0.0)) THEN
              iflag = 0
              return
            ENDIF
c            IF(ITEMP==1) THEN
c              WRITE(6,'(A,I4,3(1X,F14.8))') '2. TEST SEGTRI',CODE,T,U,V
c            ENDIF
            if (code.eq.2) call intersect_triangle(q,r,va,vc,vb,t,u,v,code)
!
            xm = q(1)+t*r(1)
            ym = q(2)+t*r(2)
            zm = q(3)+t*r(3)
!
            if (code.eq.1) iflag = 1
!
            return
          end subroutine segtriint
!
!--------------------------------------------------------------------
          subroutine dot(xy,x,y)
!--------------------------------------------------------------------
            implicit none
            real x(3),y(3),xy
!
            xy = x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
            return 
          end subroutine dot
!--------------------------------------------------------------------
          subroutine  sub(xmy,x,y)
!--------------------------------------------------------------------
            implicit none
            real x(3),y(3),xmy(3)
!
            xmy = x-y
            return 
          end subroutine sub
!--------------------------------------------------------------------
          subroutine  cross(xcy,x,y)
!--------------------------------------------------------------------
            implicit none
            real x(3),y(3),xcy(3)
!
            xcy(1) = x(2)*y(3)-y(2)*x(3)
            xcy(2) = x(3)*y(1)-y(3)*x(1)
            xcy(3) = x(1)*y(2)-y(1)*x(2)
            return 
          end subroutine cross
!
!--------------------------------------------------------------------
          subroutine intersect_triangle(orig,dir,vert0,vert1,vert2,t,u,v,intersect)
!--------------------------------------------------------------------
            implicit none
            integer intersect
            real u,v,t,det,inv_det
            real, dimension (3) ::  edge1,edge2,tvec,pvec,qvec,orig,dir,vert0,vert1,vert2
            real epsilon
            epsilon = 1.0e-18
!
            t = 0.
            u = 0.
            v = 0. 

!.....find vectors for two edges sharing vert0 
            call sub(edge1, vert1, vert0)
            call sub(edge2, vert2, vert0)

!.....begin calculating determinant - also used to calculate U parameter 
            call cross(pvec, dir, edge2)

!.....if determinant is near zero, ray lies in plane of triangle
c PVEC is the first row of the inverse matrix * its determinant
            call dot(det, edge1, pvec)

            if (det.gt.-epsilon.and.det.lt.epsilon) then
               intersect = 2
               return
            endif
            inv_det = 1.0 / det

!......calculate distance from vert0 to ray origin 
            call sub(tvec, orig, vert0)

!.....calculate U parameter and test bounds
c used to establish if the intersection DIR-PLANE is inside the triangle 
            call dot(u, tvec, pvec)
            u = u * inv_det

            if (u.lt.0.0.or.u.gt.1.0) then
               intersect = 3
               return
            endif

!.....prepare to test V parameter 
            call cross(qvec, tvec, edge1)

!.....calculate V parameter and test bounds 
c used to establish if the intersection DIR-PLANE is inside the triangle
            call dot(v, dir, qvec)
            v = v * inv_det

!!!!!!            if (v.lt.0.0.or.(u+v).gt.1.0) then
            if (v.lt.0.0.or.(u+v-1.).gt.1.e-07) then
               intersect = 3
               return
            endif

!.....calculate t, ray intersects triangle 
c t is measured along DIR
            call dot(t, edge2, qvec)
            t = t * inv_det

            intersect = 1
!
            return
          end subroutine intersect_triangle

c-----------------------------------------------------------------------
          integer function pointintriangle(p,a,b,c)
       
          integer sameside
          real p(2),a(2),b(2),c(2)

          pointintriangle = 0

          IF (sameside(p,a,b,c)>0 .AND. 
     &        sameside(p,b,a,c)>0 .AND. sameside (p,c,a,b)>0 ) THEN
            pointintriangle = 1
          ENDIF

          return

          end
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
          integer function sameside(p1,p2,a,b)

          real crossproduct2
          real p1(2),p2(2),a(2),b(2)
          real cp1,cp2
          
          cp1 = crossproduct2(a,b,p1)
          cp2 = crossproduct2(a,b,p2)

          sameside = 0
          if(cp1*cp2>=0.) sameside=1

          return
          end
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
          real function crossproduct2(a,b,c)
c
c         PURPOSE:  Given 2D points a,b, and c it calculates the 
c         cross product of vectors b-a and c-a and stores it into
c         vector v.

          real a(2),b(2),c(2)

          integer i
          real va(2),vb(2),v(2)

          do i=1,2
            va(i) = b(i)-a(i)
            vb(i) = c(i)-a(i)
          enddo

          crossproduct2 = va(1)*vb(2) - vb(1)*va(2)

          return

          end
c-----------------------------------------------------------------------
   

C---- function ray_face_int-----------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function ray_face_int(q,r,xu,yu,zu,nx,ny,nz,i,j,k,dir,xint,yint,zint,icyl)
C     
C     PURPOSE: Find intersection of a ray q+r with a grid face in 
C     direction dir
C
C---------------------------------------------------------------------------
      implicit none

      integer nx,ny,nz,i,j,k,dir,icyl
      real    xint,yint,zint
c
c.... Input/Output Arrays
      real    q(3),r(3),xu(nx,ny),yu(nx,ny),zu(nz)
      real    va(3),vb(3),vc(3),vd(3)
      real    vertex(4,3)
c
c.... Functions
      integer ray_poly_int,ray_xycylface_int
      integer ray_xzface_int,ray_yzcylface_int

      if(i<1.OR.i>nx.OR.j<1.OR.j>ny) then
        ray_face_int=0
        return
      endif

      IF(dir==1) THEN

        vertex(1,1)=xu(i,j)
        vertex(1,2)=yu(i,j)
        vertex(1,3)=zu(k)
        vertex(2,1)=xu(i,j+1)
        vertex(2,2)=yu(i,j+1)
        vertex(2,3)=zu(k)
        vertex(3,1)=xu(i,j+1)
        vertex(3,2)=yu(i,j+1)
        vertex(3,3)=zu(k+1)
        vertex(4,1)=xu(i,j)
        vertex(4,2)=yu(i,j)
        vertex(4,3)=zu(k+1)
        if(icyl==0) then
          ray_face_int = ray_poly_int(q,r,vertex,4,xint,yint,zint)
        else
          if(i==1) then
c            vertex(1:4,1) = 0.0
c            vertex(1:4,2) = 0.0
            ray_face_int = 0
          else             
            ray_face_int = ray_yzcylface_int(q,r,vertex,4,xint,yint,zint)
          endif
        endif

      ELSEIF(dir==2) THEN

        vertex(1,1)=xu(i,j)
        vertex(1,2)=yu(i,j)
        vertex(1,3)=zu(k)
        vertex(2,1)=xu(i+1,j)
        vertex(2,2)=yu(i+1,j)
        vertex(2,3)=zu(k)
        vertex(3,1)=xu(i+1,j)
        vertex(3,2)=yu(i+1,j)
        vertex(3,3)=zu(k+1)
        vertex(4,1)=xu(i,j)
        vertex(4,2)=yu(i,j)
        vertex(4,3)=zu(k+1)

        if(icyl==1.AND.i==1) then
          vertex(1,1) = 0.0
          vertex(1,2) = 0.0
          vertex(4,1) = 0.0
          vertex(4,2) = 0.0
        endif
        ray_face_int = ray_xzface_int(q,r,vertex,4,xint,yint,zint,icyl)

      ELSE

        vertex(1,1)=xu(i,j)
        vertex(1,2)=yu(i,j)
        vertex(1,3)=zu(k)
        vertex(2,1)=xu(i+1,j)
        vertex(2,2)=yu(i+1,j)
        vertex(2,3)=zu(k)
        vertex(3,1)=xu(i+1,j+1)
        vertex(3,2)=yu(i+1,j+1)
        vertex(3,3)=zu(k)
        vertex(4,1)=xu(i,j+1)
        vertex(4,2)=yu(i,j+1)
        vertex(4,3)=zu(k)

        if(icyl==0) then
          ray_face_int = ray_poly_int(q,r,vertex,4,xint,yint,zint)
        else
          if(i==1) then
            vertex(1,1) = 0.0
            vertex(1,2) = 0.0
            vertex(4,1) = 0.0
            vertex(4,2) = 0.0
          endif
          ray_face_int = ray_xycylface_int(q,r,vertex,4,xint,yint,zint)
        endif

      ENDIF

      RETURN
      END
C------------------------------------------------------------------------




C---- function ray_poly_int-------------------N. Beratlis-31 Dec. 2008---
C
      INTEGER function ray_poly_int(q,r,vertex,nv,xint,yint,zint)
C
C     PURPOSE: Find intersection of ray q+r with polygon defined by vertex
C     Algorithm found at:
c     http://www.siggraph.org/education/materials/HyperGraph/raytrace/rayplane_intersection.htm
C
C------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER nv
      REAL    xint,yint,zint
      REAL    q(3),r(3),vertex(nv,3)

      INTEGER dir,i
      REAL    a,b,c,d
      REAL    vecmag,dotprod,vd,v1,t
      INTEGER point_poly_clasf
      REAL    va(3),vb(3),vc(3),pn(3),rn(3)
      REAL    vp(nv,2)

c      write(6,*) q,r,vertex,nv
      !Define plane of polygon from the first 3 vertices
      va(1) = vertex(1,1)
      va(2) = vertex(1,2)
      va(3) = vertex(1,3)
      vb(1) = vertex(2,1)
      vb(2) = vertex(2,2)
      vb(3) = vertex(2,3)
      vc(1) = vertex(3,1)
      vc(2) = vertex(3,2)
      vc(3) = vertex(3,3)
      call plane_tri(va,vb,vc,a,b,c,d)
      pn(1)=a
      pn(2)=b
      pn(3)=c

      vd = dotprod(pn,r)


      ray_poly_int = 0
      IF(vd/=0) THEN
        v1 = -(dotprod(pn,q)+d)
        t=v1/vd
        IF(t>=0) THEN
          xint = q(1)+t*r(1)
          yint = q(2)+t*r(2)
          zint = q(3)+t*r(3)
          dir=maxloc(pn,1)
          !Project polygon to 2D plane
          IF(dir==1) THEN
            vp(1:nv,1) = vertex(1:nv,2)
            vp(1:nv,2) = vertex(1:nv,3)
            if(point_poly_clasf(yint,zint,vp,nv)==1) then
              ray_poly_int=1
            endif
          ELSEIF(dir==2) THEN
            vp(1:nv,1) = vertex(1:nv,1)
            vp(1:nv,2) = vertex(1:nv,3)
            if(point_poly_clasf(xint,zint,vp,nv)==1) then
              ray_poly_int=1
            endif
          ELSE
            vp(1:nv,1) = vertex(1:nv,1)
            vp(1:nv,2) = vertex(1:nv,2)            
            if(point_poly_clasf(xint,yint,vp,nv)==1) then
              ray_poly_int=1
            endif
          ENDIF
          
        ENDIF
      ENDIF

c      write(6,*) ray_poly_int

      RETURN
      END
C------------------------------------------------------------------------



C---- function point_poly----------------------N. Beratlis-31 Dec. 2008--
C
      INTEGER function point_poly_clasf(x,y,vertices,nv)
C
C     PURPOSE: Determine if point x,y is inside or outside polygon defined
C     by nv vertices. 
C------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER nv
      REAL    x,y
      REAL    vertices(nv,2)

      INTEGER i,cross
      INTEGER ray_seg_int
      REAL    xint,yint
      REAL    p(2),q(2),r(2),s(2)

      p(1) = x
      p(2) = y
      q(1) = 1000.
      q(2) = 0.

      cross=0
      do i=1,nv-1
        r(1) = vertices(i,1)
        r(2) = vertices(i,2)
        s(1) = vertices(i+1,1)-vertices(i,1)
        s(2) = vertices(i+1,2)-vertices(i,2)

        cross = cross+ray_seg_int(p,q,r,s,xint,yint)
      enddo

c      write(6,*) 'cross=',cross
      if(mod(cross,2)==0) then
        point_poly_clasf=0
      else
        point_poly_clasf=1
      endif

      RETURN
      END
C-----------------------------------------------------------------------


c---- function ray_xycylface_int--------------N. Beratlis-30 May 2009---
C
c     PURPOSE: Find interesection of ray q+r with xy face of cylindical
C     cell defined by 4 vertices.
C
C-----------------------------------------------------------------------
 
      integer function ray_xycylface_int(q,r,vertex,nv,xint,yint,zint)

      implicit none
c
c.... Input/Output arrays
      integer nv
      real    xint,yint,zint
      real    q(3),r(3),vertex(nv,3)
c
c.... Local arrays
      real    a,b,c,d,vd,v1,t
      real    rmin,rmax,thetamin,thetamax,thetap,rp,theta1,xp,yp,x2,y2
      real    va(3),vb(3),vc(3),pn(3),rn(3)
c
c.... Functions
      real    dotprod,anglerad

      va(1) = vertex(1,1)
      va(2) = vertex(1,2)
      va(3) = vertex(1,3)
      vb(1) = vertex(2,1)
      vb(2) = vertex(2,2)
      vb(3) = vertex(2,3)
      vc(1) = vertex(3,1)
      vc(2) = vertex(3,2)
      vc(3) = vertex(3,3)
c
c.... Find coefficients defining plane and its normal (Ax+By+Cz+D=0)      
      call plane_tri(va,vb,vc,a,b,c,d)
      pn(1)=a
      pn(2)=b
      pn(3)=c

c      write(6,*) va,vb,vc
c      write(6,*) pn,r
      vd = dotprod(pn,r)

      ray_xycylface_int = 1

c      write(6,*) 'vd=',vd
      if(vd==0) then !Ray is parallel to face
        ray_xycylface_int = 0
        return
      else
        v1 = -(dotprod(pn,q)+d)
        t=v1/vd
        if(t<0) then !Ray is pointing away from face
          ray_xycylface_int = 0
          return
        else
          xint = q(1)+t*r(1)
          yint = q(2)+t*r(2)
          zint = q(3)+t*r(3)

          rp     = sqrt(yint**2. + xint**2.)

          rmin = sqrt(vertex(1,1)**2.+vertex(1,2)**2.)
          rmax = sqrt(vertex(2,1)**2.+vertex(2,2)**2.)

          if(rp<rmin) then
            ray_xycylface_int = 0
            return
          endif

          if(rp>rmax) then
            ray_xycylface_int = 0
            return
          endif
          
          theta1 = anglerad(vertex(2,1),vertex(2,2))
c.... Perform coordinate rotation to aling vertex(1,1) with x-axis (theta1->0) 
          call cord_2drot(xint,yint,theta1,xp,yp)
          thetap = anglerad(xp,yp)
          if(thetap<0) then
            ray_xycylface_int = 0
            return
          endif
          
          call cord_2drot(vertex(3,1),vertex(3,2),theta1,x2,y2)
          thetamax = anglerad(x2,y2)
          if(thetap>thetamax) then
            ray_xycylface_int = 0
            return
          endif

        endif
      endif

      return

      end
C-----------------------------------------------------------------------


c---- function ray_xzface_int ----------------N. Beratlis-30 May 2009---
C
c     PURPOSE: Find interesection of ray q+r with xy face of cell defined
C     by 4 vertices.
C
C-----------------------------------------------------------------------
 
      integer function ray_xzface_int(q,r,vertex,nv,xint,yint,zint,icyl)

      implicit none
c
c.... Input/Output arrays
      integer nv,icyl
      real    xint,yint,zint
      real    q(3),r(3),vertex(nv,3)
c
c.... Local arrays
      real    a,b,c,d,vd,v1,t
      real    rmin,rmax,zmin,zmax,zp,rp
      real    va(3),vb(3),vc(3),pn(3),rn(3)
c
c.... Functions
      real    dotprod,anglerad

      va(1) = vertex(1,1)
      va(2) = vertex(1,2)
      va(3) = vertex(1,3)
      vb(1) = vertex(2,1)
      vb(2) = vertex(2,2)
      vb(3) = vertex(2,3)
      vc(1) = vertex(3,1)
      vc(2) = vertex(3,2)
      vc(3) = vertex(3,3)
c
c.... Find coefficiens defining plane and its normal (Ax+By+Cz+D=0)      
      call plane_tri(va,vb,vc,a,b,c,d)
      pn(1)=a
      pn(2)=b
      pn(3)=c

      vd = dotprod(pn,r)

      ray_xzface_int = 1

      if(vd==0) then !Ray is parallel to face
        ray_xzface_int = 0
        return
      else
        v1 = -(dotprod(pn,q)+d)
        t=v1/vd
c        write(6,*) 't=',t
        if(t<0) then !Ray is pointing away from face
          ray_xzface_int = 0
          return
        else
          xint = q(1)+t*r(1)
          yint = q(2)+t*r(2)
          zint = q(3)+t*r(3)
          
          if(icyl==0) then
            rmin = vertex(1,1)
            rmax = vertex(2,1)
          else
            rmin = sqrt(vertex(1,1)**2.+vertex(1,2)**2.)
            rmax = sqrt(vertex(2,1)**2.+vertex(2,2)**2.)
          endif

          if(rmin>rmax) then
            a = rmin
            rmin = rmax
            rmax = a
          endif

          zmin = vertex(1,3)
          zmax = vertex(4,3)
          if(zmin>zmax) then
            a = zmin
            zmin = zmax
            zmax = a
          endif

          if(icyl==0) then
            rp = xint
          else
            rp = sqrt(xint**2. + yint**2.)
          endif
          zp = zint

c.... Determine if points is outside radius range          
          if(rp<rmin) then
            ray_xzface_int = 0
            return
          elseif(rp>rmax) then
            ray_xzface_int = 0
            return
          endif

c.... Determine if point is outside angle range
          if(zp<zmin) then
            ray_xzface_int = 0
            return
          elseif(zp>zmax) then
            ray_xzface_int = 0
            return
          endif

        endif
      endif

      return

      end
C-----------------------------------------------------------------------


c---- function ray_yzcylface_int--------------N. Beratlis-30 May 2009---
C
c     PURPOSE: Find interesection of ray q+r with segment of cylindical
C     surface defined by 4 vertices.
C
C-----------------------------------------------------------------------
 
      integer function ray_yzcylface_int(q,r,vertex,nv,xint,yint,zint)

      implicit none
c
c.... Input/Output arrays
      integer nv
      real    xint,yint,zint
      real    q(3),r(3)
      real    vertex(nv,3)
c
c.... Local arrays
      integer bn
      real    a,b,c,R2,det,t1,t2,theta1,theta2,z1,z2,rot_rad
      real    p(3),prot(3)
c
c.... Functions
      real    anglerad
      integer bounded_cyl_seg

      theta1 = anglerad(vertex(1,1),vertex(1,2))
      theta2 = anglerad(vertex(3,1),vertex(3,2))

      rot_rad = theta1
      call cord_2drot(vertex(3,1),vertex(3,2),rot_rad,a,b)

      theta1 = 0.0
      theta2 = anglerad(a,b)

      z1 = vertex(1,3)
      z2 = vertex(3,3)

      if(z1>z2) then
        a = z1
        z1 = z2
        z2 = a
      endif

      R2 = vertex(1,1)**2.+vertex(1,2)**2.


      a = r(1)**2.+r(2)**2.
      b = 2.*(q(1)*r(1)+q(2)*r(2))
      c = q(1)**2+q(2)**2.- R2

      det = b**2-4*a*c

      if(det<0 .OR. a==0.0) then
        ray_yzcylface_int=0
      elseif(det==0) then
        t1 = -b/(2.0*a)

        p(1) = q(1)+t1*r(1)
        p(2) = q(2)+t1*r(2)
        p(3) = q(3)+t1*r(3)
        prot = p

        call cord_2drot(p(1),p(2),rot_rad,prot(1),prot(2))

        bn = bounded_cyl_seg(prot,theta1,theta2,z1,z2)
        if(bn==1) then
          xint = p(1)
          yint = p(2)
          zint = p(3)
          ray_yzcylface_int = 1
          return
        else
          ray_yzcylface_int = 0
        endif
      else
        t1 = (-b+sqrt(det))/(2*a)        
        t2 = (-b-sqrt(det))/(2*a)
       
        p(1) = q(1)+t1*r(1)
        p(2) = q(2)+t1*r(2)
        p(3) = q(3)+t1*r(3)
        prot = p

        call cord_2drot(p(1),p(2),rot_rad,prot(1),prot(2))

        bn = bounded_cyl_seg(prot,theta1,theta2,z1,z2)
        if(bn==1) then
          xint = p(1)
          yint = p(2)
          zint = p(3)
          ray_yzcylface_int = 1
          return
        elseif(t2>0.0) then
          p(1) = q(1)+t2*r(1)
          p(2) = q(2)+t2*r(2)
          p(3) = q(3)+t2*r(3)
          prot = p

          call cord_2drot(p(1),p(2),rot_rad,prot(1),prot(2))

          bn = bounded_cyl_seg(prot,theta1,theta2,z1,z2)
          if(bn==1) then
            xint = p(1)
            yint = p(2)
            zint = p(3)
            ray_yzcylface_int = 1
            return
          else
            ray_yzcylface_int = 0
          endif

        else
          ray_yzcylface_int = 0
        endif

      endif

      return

      end
C-----------------------------------------------------------------------


c---- function bounded_cyl_seg----------------N. Beratlis-30 May 2009---
C
C     PURPOSE: Check if point p is bounded by segment of cylindrical 
C     surface defined by theta1-theta2 and z1,z2
C
C-----------------------------------------------------------------------
      integer function bounded_cyl_seg(p,theta1,theta2,z1,z2)
      
      implicit none
c
c.... Input/Output arrays
      real z1,z2
      real theta1,theta2
      real p(3)
c
c.... Local arrays
      integer b
      real theta,z
c
c.... Functions
      real anglerad

      theta = anglerad(p(1),p(2))
      z     = p(3)

      b=1

      if(theta<theta1) b=0
      if(theta>theta2) b=0
      if(z<z1) b=0
      if(z>z2) b=0

      bounded_cyl_seg = b

      return

      end
C-----------------------------------------------------------------------





C---- function seg_seg_int----------------------N. Beratlis-31 Dec. 2008-
C
      INTEGER function seg_seg_int(va1,va2,vb1,vb2,x,y)
C
C     PURPOSE: Find if line segment defined by endpoints va1,va2 intersects
C     line segment defined by endpoints vb1,vb2
C
C------------------------------------------------------------------------
      IMPLICIT NONE

      REAL    x,y
      REAL    va1(2),va2(2),vb1(2),vb2(2)

      REAL    a1,b1,a2,b2

      seg_seg_int = 0

      a1 = (va2(2)-va1(2))/(va2(1)-va1(1))
      b1 = va1(2)-a1*va1(1)

      a2 = (vb2(2)-vb1(2))/(vb2(1)-vb1(1))
      b2 = vb1(2)-a2*vb1(1)
      
      if(a1/=a2) then
        x = (b2-b1)/(a1-a2)
        y = a1/(a1-a2)*(b2-b1)+b1
        seg_seg_int = 1
      endif

      RETURN
      END
C------------------------------------------------------------------------





C---- function ray_seg_int----------------------N. Beratlis-31 Dec. 2008-
C
      INTEGER function ray_seg_int(p,q,r,s,x,y)
C
C     PURPOSE: Find if ray defined by P+a*Q intersects line segment defined
C     by R+b*S, where 0<a<1, 0<b<1.
C
C------------------------------------------------------------------------
      IMPLICIT NONE

      REAL    x,y
      REAL    p(2),q(2),r(2),s(2)

      REAL    d0,a,b
      REAL    dotproduct
      REAL    f(2),g(2)

      ray_seg_int = 0

      f(1) = s(2)
      f(2) =-s(1) 
      d0 = dotproduct(f,q,2)
c      write(6,*) 'f=',f,',q=',q
c      write(6,*) 'd0=',d0
      if(d0/=0) then
        g(1) = r(1)-p(1)
        g(2) = r(2)-p(2)
c        write(6,*) 'f=',f,',g=',g
        a = dotproduct(f,g,2)/d0
        if(a>=0) then
          f(1) = q(2)
          f(2) =-q(1)
          b = dotproduct(f,g,2)/d0
c          WRITE(6,*) a,b
          if(b>=0.0 .AND. b<=1.0) then
            ray_seg_int=1
            x=p(1)+a*q(1)
            y=p(2)+a*q(2)
         endif
        endif
      endif

      RETURN
      END
C------------------------------------------------------------------------


C---- function dotproduct----------------------N. Beratlis-31 Dec. 2008--
C
      REAL FUNCTION DOTPRODUCT(a,b,n)
C
C     PURPOSE: Calculate dot product between n-th order vectors a and b
C
C------------------------------------------------------------------------
c
      IMPLICIT NONE
c
      INTEGER n
      REAL    a(n),b(n)
c
c.....Local arrays
      INTEGER I
c
      dotproduct = 0.
      DO i=1,n
        dotproduct = dotproduct + a(i)*b(i)
      ENDDO

      RETURN

      END
C------------------------------------------------------------------------




C---- subroutine plane_tri----------------------N. Beratlis-31 Dec. 2008-
C
      subroutine plane_tri(va,vb,vc,a,b,c,d)
C
C     PURPOSE: Get coefficients for plane defined by A*x+B*y+C*z+D=0 
C     using 3 points, va,vb,vc
C
C------------------------------------------------------------------------
      REAL a,b,c,d
      REAL va(3),vb(3),vc(3)

      REAL s1(3),s2(3),un(3)

      s1(1) = vb(1)-va(1)
      s1(2) = vb(2)-va(2)
      s1(3) = vb(3)-va(3)
      s2(1) = vc(1)-va(1)
      s2(2) = vc(2)-va(2)
      s2(3) = vc(3)-va(3)

      CALL cross(un,s1,s2)
      A = un(1)
      B = un(2)
      C = un(3)
      D = -(A*va(1)+B*va(2)+C*va(3))

      RETURN
      END
C------------------------------------------------------------------------


C---- function vecmag ------------------------N. Beratlis-31 Dec. 2008--
C
C     PURPOSE: Calculates magnitude of a vector
C
C-----------------------------------------------------------------------
      real function vecmag(a,n)

      integer n
      real    a(n)

      vecmag=0.0
      do i=1,n
        vecmag=vecmag+a(i)**2. 
      enddo
      vecmag=sqrt(vecmag)

      return

      end
C-----------------------------------------------------------------------




C---- function grid2grid_intr------------------N. Beratlis-16 Jan. 2009-
C
C     PURPOSE: From point xp,yp,zp lying on cell face aligned with plane 
C     "dir" find intersection with next cell face along vector "rvec"
C
C-----------------------------------------------------------------------
      integer function grid2grid_intr(xp,yp,zp,xinp,yinp,zinp,r
     &     ,xc_car,yc_car,zc,nx,ny,nz,i,j,k,i1,j1,k1,dir,dir1,icyl,ifacet)

      IMPLICIT NONE
c      include 'common.h'
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output
      integer nx,ny,nz,i,j,k,i1,j1,k1,icyl,dir,dir1,ifacet
      real    xp,yp,zp,xinp,yinp,zinp
      real    xc_car(nx,ny),yc_car(nx,ny),zc(nz)
      real    r(3),rext(3)
c
c.... Local arrays      
      integer i2,j2,k2
      integer iext,jext,kext,iext1,jext1,kext1,iext2,jext2,kext2
      integer iflagt,iflag1,iflag2,iflag3
      integer iexta,jexta,kexta
      real    q(3),theta,dely,jy,dz,rm
c
c.... Functions
      integer ray_face_int
      real     vecmag,anglerad

      q(1) = xp
      q(2) = yp
      q(3) = zp

      iflag1 = 0
      iflag2 = 0
      iflag3 = 0

      rext = r/vecmag(r,3)*0.01

      jy = -1.5
      if(yc_car(2,1)==0.0) jy=-1.0

      dely = 2.*pi/(ny-2)
      theta = dely*(real(j)+jy)


      iext=0
      jext=0
      kext=0

      iext2 = -1
      jext2 = -1
      kext2 = -1

      iflagt=0

      if(dir==1) then
    
        rm = sqrt(xc_car(i-1,j)**2. + yc_car(i-1,j)**2.)
        if(icyl==1 .AND. i==2) rm=0.0
        call vec_ijkext_face(q,r,iext,jext,kext,icyl,1,rm,dely,dir)

        i2 = i+iext
        j2 = j
        k2 = k
        if(icyl==1) call index_bnd(i+iext,j,k,i2,j2,k2,nx,ny,nz)
        iflag1 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &        ,i2,j2,k2,1,xinp,yinp,zinp,icyl)

        if(iflag1==1) then
          iflagt=1
          dir1=1
          i1 = i2
          j1 = j2
          k1 = k2
        else

          call vec_ijkext_face(q,r,iext,jext,kext,icyl,2,theta,dely,dir) 
          i2 = i+iext
          j2 = j+jext
          k2 = k
          if(icyl==1) call index_bnd(i+iext,j+jext,k,i2,j2,k2,nx,ny,nz)
          iflag2 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &         ,i2,j2,k2,2,xinp,yinp,zinp,icyl)

          if(iflag2==1) then
            iflagt=1
            dir1=2
            i1 = i2
            j1 = j2
            k1 = k2
          else

            dz = zc(k)-q(3)
            if(r(3)>=0.0) then
              dz = zc(k+1)-q(3)
            endif
            call vec_ijkext_face(q,r,iext,jext,kext,icyl,3,dz,dely,dir)

            i2 = i+iext
            j2 = j
            k2 = k+kext
            if(icyl==1) call index_bnd(i+iext,j,k+kext,i2,j2,k2,nx,ny,nz)
            iflag3 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &         ,i2,j2,k2,3,xinp,yinp,zinp,icyl)

            if(iflag3==1) then
              iflagt=1
              dir1=3
              i1 = i2
              j1 = j2
              k1 = k2
            endif
          endif
        endif
      elseif(dir==2) then
        
        rm = sqrt(xc_car(i,j)**2. + yc_car(i,j)**2.)
        if(icyl==1 .AND. i==1) rm=0.0
        call vec_ijkext_face(q,r,iext,jext,kext,icyl,1,rm,dely,dir) 
        i2 = i+iext
        j2 = j+jext
        k2 = k
        if(icyl==1) call index_bnd(i+iext,j+jext,k,i2,j2,k2,nx,ny,nz)

        iflag1 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &        ,i2,j2,k2,1,xinp,yinp,zinp,icyl)

        if(iflag1==1) then
          iflagt=1
          dir1=1
          i1 = i2
          j1 = j2
          k1 = k2
        else

          call vec_ijkext_face(q,r,iext,jext,kext,icyl,2,theta,dely,dir) 
          i2 = i
          j2 = j+jext
          k2 = k
          if(icyl==1) call index_bnd(i,j+jext,k,i2,j2,k2,nx,ny,nz)

          iflag2 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &         ,i2,j2,k2,2,xinp,yinp,zinp,icyl)
          if(iflag2==1) then
            iflagt=1
            dir1=2
            i1 = i2
            j1 = j2
            k1 = k2
          else

            dz = q(3)-zc(k)
            if(r(3)>=0.0) then
              dz = zc(k+1)-q(3)
            endif
            call vec_ijkext_face(q,r,iext,jext,kext,icyl,3,dz,dely,dir)

            i2 = i
            j2 = j+jext
            k2 = k+kext
            if(icyl==1) call index_bnd(i,j+jext,k+kext,i2,j2,k2,nx,ny,nz)
            iflag3 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &            ,i2,j2,k2,3,xinp,yinp,zinp,icyl)
            if(iflag3==1) then
              iflagt=1
              dir1=3
              i1 = i2
              j1 = j2
              k1 = k2
            endif
          endif
        endif
      elseif(dir==3) then

        rm = sqrt(xc_car(i,j)**2. + yc_car(i,j)**2.)
        if(icyl==1 .AND. i==1) rm=0.0
        call vec_ijkext_face(q,r,iext,jext,kext,icyl,1,rm,dely,dir)

        i2 = i+iext
        j2 = j
        k2 = k+kext
        if(icyl==1) call index_bnd(i+iext,j,k+kext,i2,j2,k2,nx,ny,nz)

        iflag1 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &        ,i2,j2,k2,1,xinp,yinp,zinp,icyl)

        if(iflag1==1) then
          iflagt=1
          dir1=1
          i1 = i2
          j1 = j2
          k1 = k2
        else

          call vec_ijkext_face(q,r,iext,jext,kext,icyl,2,theta,dely,dir)

          i2 = i
          j2 = j+jext
          k2 = k+kext
          if(icyl==1) call index_bnd(i,j+jext,k+kext,i2,j2,k2,nx,ny,nz)
          iflag2 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &         ,i2,j2,k2,2,xinp,yinp,zinp,icyl)

          if(iflag2==1) then
            iflagt=1
            dir1=2
            i1 = i2
            j1 = j2
            k1 = k2
          else

            dz = zc(k + int(sign(1.0,r(3))))-zc(k)
            call vec_ijkext_face(q,r,iext,jext,kext,icyl,3,dz,dely,dir)

            i2 = i
            j2 = j
            k2 = k+kext
            if(icyl==1) call index_bnd(i,j,k+kext,i2,j2,k2,nx,ny,nz)

            iflag3 = ray_face_int(q,r,xc_car,yc_car,zc,nx,ny,nz
     &         ,i2,j2,k2,3,xinp,yinp,zinp,icyl)

            if(iflag3==1) then
              iflagt=1
              dir1=3
              i1 = i2
              j1 = j2
              k1 = k2
            endif
          endif
        endif
      endif


      if(iflag1==0.AND.iflag2==0.AND.iflag3==0) then
        write(6,*) 'WARNING: grid2grid_intr no face intersection',ifacet,i,j,k,dir,q,r
        write(6,*) '     q(1)=',q(1)
        write(6,*) '     q(2)=',q(2)
        write(6,*) '     q(3)=',q(3)
        write(6,*) '     r(1)=',r(1)
        write(6,*) '     r(2)=',r(2)
        write(6,*) '     r(3)=',r(3)

        if(dir==1) then

          rm = sqrt(xc_car(i-1,j)**2. + yc_car(i-1,j)**2.)
          if(icyl==1 .AND. i==1) rm=0.0
          call vec_ijkext_face(q,r,iext,jext,kext,icyl,1,rm,dely,dir) 
          call index_bnd(i+iext,j+jext,k,i2,j2,k2,nx,ny,nz)
          write(6,*) '1.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,rm

          call vec_ijkext_face(q,r,iext,jext,kext,icyl,2,theta,dely,dir) 
          call index_bnd(i,j+jext,k,i2,j2,k2,nx,ny,nz)
          write(6,*) '2.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,theta

          dz = q(3)-zc(k)
          if(r(3)>=0.0) then
            dz = zc(k+1)-q(3)
          endif
          call vec_ijkext_face(q,r,iext,jext,kext,icyl,3,dz,dely,dir)
          write(6,*) '3.',ifacet,i,j,k,iext,jext,kext,dz

        elseif(dir==2) then

          rm = sqrt(xc_car(i,j)**2. + yc_car(i,j)**2.)
          if(icyl==1 .AND. i==1) rm=0.0
          call vec_ijkext_face(q,r,iext,jext,kext,icyl,1,rm,dely,dir) 
          call index_bnd(i+iext,j+jext,k,i2,j2,k2,nx,ny,nz)
          write(6,*) '1.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,rm

          call vec_ijkext_face(q,r,iext,jext,kext,icyl,2,theta,dely,dir) 
          call index_bnd(i,j+jext,k,i2,j2,k2,nx,ny,nz)
          write(6,*) '2.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,theta

          dz = q(3)-zc(k)
          if(r(3)>=0.0) then
            dz = zc(k+1)-q(3)
          endif
          call vec_ijkext_face(q,r,iext,jext,kext,icyl,3,dz,dely,dir)
          call index_bnd(i,j+jext,k+kext,i2,j2,k2,nx,ny,nz)
          write(6,*) '3.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,dz

        else

          rm = sqrt(xc_car(i,j)**2. + yc_car(i,j)**2.)
          if(icyl==1 .AND. i==1) rm=0.0
          call vec_ijkext_face(q,r,iext,jext,kext,icyl,1,rm,dely,dir)
          call index_bnd(i+iext,j,k+kext,i2,j2,k2,nx,ny,nz)
          write(6,*) '1.',ifacet,i,j,k,iext,jext,kext,rm

          call vec_ijkext_face(q,r,iext,jext,kext,icyl,2,theta,dely,dir)
          call index_bnd(i,j+jext,k+kext,i2,j2,k2,nx,ny,nz)
          write(6,*) '2.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,theta

          dz = zc(k + int(sign(1.0,r(3))))-zc(k)
          call vec_ijkext_face(q,r,iext,jext,kext,icyl,3,dz,dely,dir)
          call index_bnd(i,j,k+kext,i2,j2,k2,nx,ny,nz)
          write(6,*) '3.',ifacet,i,j,k,iext,jext,kext,i2,j2,k2,dz

        endif

c        stop
      endif

      grid2grid_intr=iflagt

      return
      end
C-----------------------------------------------------------------------






c---- subroutine cord_2drot------------------N. Beratlis-31 May 2009---
C
C     PURPOSE: Perform coordinate axis rotation of point x,y by angle theta
C
C-----------------------------------------------------------------------
      subroutine cord_2drot(x,y,theta,xp,yp)

      implicit none

      real x,y,xp,yp,theta

      xp = y*sin(theta) + x*cos(theta)
      yp = y*cos(theta) - x*sin(theta)

      return

      end
C-----------------------------------------------------------------------

c---- subroutine tri_closest_pt----------------------N. Beratlis-18 Jun. 2011---
C
C     PURPOSE: Find closest point to triangle. Original algorithm from:
C     http://www.geometrictools.com/LibMathematics/Distance/Distance.html     
C
c-------------------------------------------------------------------------------
      subroutine tri_closest_pt(vertex,p,q)
c
      implicit none
c     
c.... Input/Output arrays
      real p(3),vertex(3,3),q(3)
c
c.... Local arrays
      real a00,a01,a11,b0,b1,c,s,t,det,sd,invdet,tmp0,tmp1,numer,denom
      real E0(3),E1(3),D0(3),diff(3)
c
c.... Function
      real dotproduct

      E0 = vertex(:,2)-vertex(:,1)
      E1 = vertex(:,3)-vertex(:,1)
      diff = vertex(:,1)-p
c
      a00 = dotproduct(E0,E0,3)
      a01 = dotproduct(E0,E1,3)
      a11 = dotproduct(E1,E1,3)
      b0 = dotproduct(E0,diff,3)
      b1 = dotproduct(E1,diff,3)
      c = dotproduct(diff,diff,3)

      det = abs(a00*a11-a01*a01)
      s = a01*b1 - a11*b0
      t = a01*b0 - a00*b1

      if(s+t<=det) then
        if(s<0) then
          if(t<0) then !Region 4
            if(b0<0.0) then
              t = 0.0
              if(-b0 >= a00) then
                s = 1.0
                sd = a00 + 2.0*b0 + c
              else
                s = -b0/a00
                sd = b0*s + c
              endif
            else
              s = 0.0
              if(b1>0.0) then
                t = 0.0
                sd = c
              elseif(-b1 >= a11) then
                t = 1.0
                sd = a11 + 2.0*b1 + c
              else
                t = -b1/a11
                sd = b1*t + c
              endif
            endif
          else !Region 3
            s = 0.0
            if(b1>=0.0) then
              t = 0.0
              sd = c
            elseif(-b1>a11) then
              t = 1.0
              sd = a11 + 2.0*b1 + c
            else
              t = -b1/a11
              sd = b1*t + c
            endif
          endif
        elseif(t<0.0) then !Region 5
          t = 0.0
          if(b0 >= 0.0) then
            s = 0.0
            sd = c
          elseif (-b0 >= a00) then
            s = 1.0
            sd = a00 + 2.0*b0 + c
          else
            s = -b0/a00
            sd = b0*s + c
          endif
        else !Region 0
          invDet = 1.0/det
          s = s*invDet
          t = t*invDet
          sd = s*(a00*s + a01*t + 2.0*b0) + t*(a01*s + a11*t + 2.0*b1) + c
        endif
      else
        if(s<0.0) then !Region 2
          tmp0 = a01 + b0
          tmp1 = a11 + b1
          if(tmp1>tmp0) then
            numer = tmp1-tmp0
            denom = a00 - 2.0*a01 + a11
            if(numer>=denom) then
              s = 1.0
              t = 0.0
              sd = a00 + 2.0*b0 + c
            else
              s = numer/denom
              t = 1.0-s
              sd = s*(a00*s + a01*t + 2.0*b0) + t*(a01*s + a11*t + 2.0*b1) + c
            endif
          else
            s = 0.0
            if(tmp1<=0.0) then
              t = 1.0
              sd = a11 + 2.0*b1 + c
            elseif(b1>=0.0) then
              t = 0.0
              sd = c
            else
              t = -b1/a11
              sd = b1*t + c
            endif
          endif
        elseif(t<0.0) then !Region 6
          tmp0 = a01 + b1
          tmp1 = a00 + b0
          if(tmp1>tmp0) then
            numer = tmp1 - tmp0
            denom = a00 - 2.0*a01 + a11
            if(numer>=denom) then
              t = 1.0
              s = 0.0
              sd = a11 + 2.0*b1 + c
            else
              t = numer/denom
              s = 1.0 - t
              sd = s*(a00*s + a01*t + 2.0*b0) + t*(a01*s + a11*t + 2.0*b1) + c
            endif
          else
            t = 0.0
            if(tmp1<0.0) then
              s = 1.0
              sd = a00 + 2.0*b0 + c
            elseif(b0 >= 0.0) then
              s = 0.0
              sd = c
            else
              s = -b0/a00
              sd = b0*s + c
            endif
          endif
        else !region 1
          numer = a11 + b1 - a01 - b0
          if(numer<=0.0) then
            s = 0.0
            t = 1.0
            sd = a11 + 2.0*b1 + c
          else
            denom = a00 - 2.0*a01 + a11
            if(numer >=denom) then
              s = 1.0
              t = 0.0
              sd = a00 + 2.0*b0 + c
            else
              s = numer/denom
              t = 1.0-s
              sd = s*(a00*s + a01*t +2.0*b0) + t*(a01*s + a11*t + 2.0*b1) + c
            endif
          endif
        endif
      endif
      
      if(sd<0.0) sd = 0.0
      
      q = vertex(:,1) + s*E0 + t*E1

      return

      end
C-------------------------------------------------------------------------------


c---------------------------------------------------------------------
C
C     PURPOSE: Find points that are outside or inside the immersed object
C
C---------------------------------------------------------------------
C
        SUBROUTINE TAG3D_MOD(X,Y,Z,NI,NJ,NK
     &                  ,XYZB,XYZCB,NXYZB,NB,S
     &                  ,IBS,IBD,ICOM)

        INCLUDE 'common.h'
        INCLUDE 'immersed.h'
        INCLUDE 'mpif.h'

        INTEGER NCR
        PARAMETER (NCR = 10)   !!!!!!

        INTEGER NI,NJ,NK,NB,IBD,ICOM
        INTEGER IBS(NI,NJ,NK)
        REAL S(NI), DS
        REAL X(NI,NJ),Y(NI,NJ),Z(NK)
        REAL XYZB(3,3,NB),XYZCB(3,NB),NXYZB(3,NB)
c     
c...  Local arrays
        INTEGER ICROSS,IFLAG,I1,I2,ISAVE,N,NINS,NEXT,NINSG,ISGN,II,IMN,IMX
        INTEGER I,J,K,IC,JN,ISGNT,NNMAX,NNMAX2,NOUT,NQ,NN,NN2,IN,NC,NDP,N1,N2
        INTEGER JJ1,JJ2,J1,J2,K1,K2,NDJ,IMIN,IMAX,IBMAXG,NCROSSMAXG
        INTEGER FIRST,LAST,MYFIRST,MYLAST
        INTEGER STATUS(MPI_STATUS_SIZE)

        REAL    QUERYXY(2)!,XYB(2,NB)
        REAL    DOTPROD
        REAL    SGN,DIST,SGN1

        INTEGER NCROSS(NI,NJ),NCROSSG(NI,NJ),NCRPRV(NJ)
        REAL    ZCROSS(NCR,NJ),ZCROSSO(NCR,NI,NJ)
        REAL    ZCROSSG(NCR,NJ),ZCROSSOG(NCR,NI,NJ)

        REAL    xcc,ycc,zcc,xm,ym,zm,ztmin,rzmin,rzmax,tol,minz,maxz,minzg,maxzg
        REAL    xa,ya,za,xb,yb,zb,xc,yc,zc,rayzmin,rayzmax
        REAL    pcc(3),va(3),vb(3),vc(3),qvec(3),rvec(3),pint(3)

        INTEGER, DIMENSION(:), ALLOCATABLE :: INN,INN2,NNTRI

        INTEGER COUNTS,RATE,NCLOCKS
        REAL    CLOCKTEMP1,CLOCKTEMP2,CLOCKTEMP3,TCLOCK

        REAL, DIMENSION(:), ALLOCATABLE :: CLOCK,CLOCKG,CLOCKGMIN,CLOCKGMAX
c
c..... Functions
        REAL    anglerad

        NCLOCKS = 14
        ALLOCATE(CLOCK(NCLOCKS),CLOCKG(NCLOCKS),CLOCKGMIN(NCLOCKS),CLOCKGMAX(NCLOCKS))
        CLOCK = 0.0
        CLOCKG = 0.0
        CLOCKGMIN = 0.0
        CLOCKGMAX = 0.0

        rayzmin = min(zmin-zlen,-1.0e8)
        rayzmax = max(zmax+zlen,1.0e8)

        CLOCK(1) = tclock()
        
        NNMAX = 1000
        NNMAX2 = -1
        N1 = 0
        N2 = 0

        ALLOCATE(INN(NNMAX))
        ALLOCATE(NNTRI(NNMAX))
        NNTRI = 0

        IF(IMBDCMP==0) THEN
c
c in the case the immersed-boundary is not decomposed
c   
          CLOCK(2) = tclock()
          IF(ICOM==0) THEN
c
c initialization call (ICOM=0)
c
            CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 0, MYRANK, DS)
          ENDIF
          CLOCK(2) = tclock()-CLOCK(2)
        
          IBS = 1

          if(idomy==0) then

            IF(MYSIZE<=NJ-2) THEN
              JJ1 = 2
              JJ2 = (NJ-2)/MYSIZE
              J1=JJ1 + JJ2*MYRANK
              J2=J1 + JJ2 -1
              NDJ = J2-J1+1
              NQ = NI*NDJ
              IMN = IBMIN(IBD)
              IMX = IBMAX(IBD)
            ELSE
              JJ1 = 2
              JJ2 = MYSIZE/(NJ-2)
c              JJ2 = (NJ-2)/MYSIZE
              J1 = JJ1 + REAL(NJ-2)/REAL(MYSIZE)*MYRANK
              J2 = J1 + REAL(NJ-2)/REAL(MYSIZE)

              CALL MPI_ALLREDUCE(IBMAX(IBD),IMX,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)
              IMX = IMX + JJ2 - MOD(IMX-IBMIN(IBD)+1,JJ2)
              IMN = MOD(MYRANK,JJ2)*(IMX-IBMIN(IBD)+1)/JJ2+IBMIN(IBD)
              IMX = IMN+(IMX-IBMIN(IBD)+1)/JJ2-1              
!              IMN = MOD(MYRANK,JJ2)*REAL(IMX-IBMIN(IBD)+1)/REAL(JJ2)+IBMIN(IBD)
!              IMX = IMN+REAL(IMX-IBMIN(IBD)+1)/REAL(JJ2)-1              
            ENDIF

          else
            J1 = 2
            J2 = NJ-1
          endif

          nins=0
          next=0
          
          IMIN = IBMIN(ibd)
          CALL MPI_ALLREDUCE(IBMAX(IBD),IMAX,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

c          DO I=IMIN,IBMAX(IBD)
          NCROSS = 0
          ZCROSSO = 0.0

          DO I=IMN,IMX

            ZCROSS = 0.0
            
            CLOCKTEMP1 = tclock()

c
c DS is a function of the X position
c
            DS = 1.0*S(I)**2.            
c            DS = 1.0*MAXVAL(S)**2.            
c            write(6,*) i,s(i)
c            DS = 2.0*S(I)**2.            
          
            DO J = J1,J2

              xcc = X(I,J)
              ycc = Y(I,J)
              QUERYXY(1) = XCC
              QUERYXY(2) = YCC
              ICROSS = 0
              CLOCKTEMP2 = tclock()
c
c committee4 finds the centers of the triangles closest to the grid point
c
              CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNTR2D, 2, 1, MYRANK, DS)
              CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2

              IF(NN>NNMAX2) NNMAX2=NN

              IF(NN>NNTR2D) THEN
                N1 = N1+1
              ELSE
                N2 = N2+1
              ENDIF

              IF(NN>NNTR2D) THEN                     
                IF(NN>NNMAX) THEN
                  clocktemp2 = tclock()
                  DEALLOCATE(INN,NNTRI)
                  NNMAX = NN
                  ALLOCATE(INN(NNMAX),NNTRI(NNMAX))
                  clock(12) = clock(12) + tclock() - clocktemp2
                ENDIF
                CLOCKTEMP2 = tclock()
                CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NN, 2, 1, MYRANK, DS)
                CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2
              ENDIF

              clocktemp3 = tclock()
              DO IN=1,MIN(NN,NNMAX)     
                N = INN(IN)
                xa = xyzb(1,1,n)
                ya = xyzb(2,1,n)
                za = xyzb(3,1,n)
                xb = xyzb(1,2,n)
                yb = xyzb(2,2,n)
                zb = xyzb(3,2,n)
                xc = xyzb(1,3,n)
                yc = xyzb(2,3,n)
                zc = xyzb(3,3,n)
                IF ((ycc.ge.min(ya,yb,yc).and.ycc.le.max(ya,yb,yc)) .and.  
     &             (xcc.ge.min(xa,xb,xc).and.xcc.le.max(xa,xb,xc)) ) THEN
c
c segtriint looks for an intersection of the triangle N with the control
c ray with equations X=XCC and Y=YCC 
c
                  CALL segtriint(xcc,ycc,rayzmin,xcc,ycc,rayzmax
     &               ,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
                  IF (iflag.eq.1) then
                    NNTRI(IN) = NNTRI(IN)+1
                    ICROSS = ICROSS+1
                    IF(ICROSS>NCR) THEN
                      WRITE(6,*) '1. WARNING MYRANK=',MYRANK,J,I,ICROSS,XCC,YCC,icross,IBD
c                      do ii=1,ncr
c                        write(6,*) ii,zcross(ii,j)
c                      enddo
                    ENDIF
                    ZCROSS(ICROSS,J) = ZM

                  ENDIF

                ENDIF
              ENDDO
              clock(9) = clock(9) + tclock() - clocktemp3

              NCROSS(I,J) = ICROSS
              IF(NCROSS(I,J)>NCR) THEN
                WRITE(6,*) '2. WARNING MYRANK=',MYRANK,J,I,NCROSS(I,J),IBD
              ENDIF

              clocktemp3 = tclock()
              !Sort z-coord. crossings in ascending order
              DO I1=1,NCROSS(I,J)
                ztmin = 1.0e10
                do i2=1,NCROSS(I,J)
                  if (ZCROSS(i2,J).le.ztmin) then
                    ztmin = ZCROSS(i2,J)
                    isave = i2;
                  endif
                enddo
c
c in ZCROSSO the Z coordinates of the intersection points are organized
c in ascending order
c
                ZCROSSO(i1,i,J) = ztmin
                ZCROSS (isave,J) = 2.0e10
              ENDDO
              clock(10) = clock(10) + tclock() - clocktemp3

              clocktemp3 = tclock()
              !Remove duplicate z-coord. crossings
                zcross(:,j) = zcrosso(:,i,j)
                nc=0
                DO i1=1,ncross(i,j)
                  ndp=0
                  DO i2=1,ncross(i,j)
!!!!!!                    if(zcross(i2,j)==zcross(i1,j).AND.i1>i2) then
                    if(abs(zcross(i2,j)-zcross(i1,j)).lt.1.e-06.AND.i1>i2) then
                      ndp=ndp+1
                    endif
                  ENDDO

!!!!!!                  if(ndp<=1) then
                  if(ndp==0) then
                    nc=nc+1
                    zcrosso(nc,i,j) = zcross(i1,j)
                  endif
                ENDDO
              clock(11) = clock(11) + tclock() - clocktemp3
c
c actual number of intersection points
c
              ncross(i,j) = nc
              if(ncross(i,j).eq.1) ncross(i,j)=0

c
c warning if the number of intersection points is odd
c
              IF(mod(ncross(i,j),2)/=0) THEN
                write(6,*) 'WARNING:',i,j,ncross(i,j),zcrosso(1:ncross(i,j),i,j)
              ENDIF

            ENDDO

            CLOCK(4) = CLOCK(4) + tclock() - CLOCKTEMP1
          ENDDO

          if(mysize>1 .AND. idomy==0) then
            clocktemp1 = tclock()
c            if(mysize>nj-2) then
            CALL MPI_ALLREDUCE(MAXVAL(NCROSS(IMIN:IMAX,2:NJ-1)),NCROSSMAXG,1,MPI_INTEGER
     &           ,MPI_MAX,MPI_COMM_EDDY,IERR)
              CALL MPI_ALLREDUCE(NCROSS(IMIN:IMAX,2:NJ-1),NCROSSG(IMIN:IMAX,2:NJ-1)
     &           ,(IMAX-IMIN+1)*(NJ-2),MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
              CALL MPI_ALLREDUCE(ZCROSSO(1:NCROSSMAXG,IMIN:IMAX,2:NJ-1),ZCROSSOG(1:NCROSSMAXG,IMIN:IMAX,2:NJ-1)
     &             ,NCROSSMAXG*(IMAX-IMIN+1)*(NJ-2),MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
c            else
c              CALL MPI_ALLREDUCE(NCROSS,NCROSSG,NI*NJ,MPI_INTEGER,MPI_SUM,MPI_COMM_EDDY,IERR)
c              CALL MPI_ALLREDUCE(ZCROSSO(1:NCROSSMAXG,1:NI,1:NJ),ZCROSSOG(1:NCROSSMAXG,1:NI,1:NJ)
c     &              ,NCROSSMAXG*NI*NJ,MTYPE,MPI_SUM,MPI_COMM_EDDY,IERR)
c            endif

            NCROSS = NCROSSG
            ZCROSSO = ZCROSSOG
            clock(13) = clock(13) + tclock() - clocktemp1
          endif

          CLOCKTEMP1 = tclock()
c
c if the number of the intersections upstream the grid point is odd
c IBS=-1 and the point is tagged as interior
c
          DO I=IBMIN(IBD),IBMAX(IBD)
            DO J=2,NJ-1
              DO ICROSS=1,NCROSS(I,J),2
                CALL LOCATE(z,nk,zcrosso(icross,i,j),k1)
                CALL LOCATE(z,nk,zcrosso(icross+1,i,j),k2)
                if(k2>k1) then
                  k1 = k1+1
                  ibs(i,j,k1:k2) = -1
                  nins = nins+(k2-k1+1)
                  next = next+(nk-2)-(k2-k1+1)
                endif
              ENDDO
            ENDDO

            CLOCK(6) = CLOCK(6) + tclock() - CLOCKTEMP1

c            call mpi_finalize(ierr)
c            stop

          ENDDO

          CLOCK(14) = tclock()
          if(idomy==0) then
            IBS(IMIN:IBMAX(IBD),1 ,KBMIN(IBD):KBMAX(IBD)) = IBS(IMIN:IBMAX(IBD),NJ-1,KBMIN(IBD):KBMAX(IBD))
            IBS(IMIN:IBMAX(IBD),NJ,KBMIN(IBD):KBMAX(IBD)) = IBS(IMIN:IBMAX(IBD), 2 ,KBMIN(IBD):KBMAX(IBD))
          endif
c
c update of IBS for the ghost layers between processors
c
          call refreshflag(ibs,ni,nj,nk)

          CLOCK(14) = tclock() - clock(14)

          CLOCK(7) = tclock()
c
c deallocation of the variables for COMMITTEE4
c
          IF(ICOM==2) CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 2, MYRANK, DS)
          clock(7) = tclock()-clock(7)

        ELSE
c
c in the case the immersed-boundary is decomposed
c
c          MINZ = MINVAL(XYZB(3,:,:))
c          MAXZ = MAXVAL(XYZB(3,:,:))

c          IF(NB==0) MINZ=100000.0
c          CALL MPI_ALLREDUCE(MINZ,MINZG,1,MTYPE,MPI_MIN,MPI_COMM_EDDY,IERR)
c          IF(NB==0) MAXZ=-100000.0
c          CALL MPI_ALLREDUCE(MAXZ,MAXZG,1,MTYPE,MPI_MAX,MPI_COMM_EDDY,IERR)

c          write(6,*) 'myrank=',myrank,'minz=',minz,'minzg=',minzg
c     &         ,',maxz=',maxz,'maxzg=',maxzg,'z:',z(1),z(nk-1)
c          call mpi_finalize(ierr)
c          stop

          FIRST=MYSIZE+1
          LAST=-1
          IF(NB>0) THEN
            FIRST=MYRANK
            LAST=MYRANK           
          ENDIF

c          IF(MINZG>Z(1) .AND. MINZG<=Z(NK-1)) THEN
c            FIRST=MYRANK
c          ELSE
c            FIRST=-1
c          ENDIF

          CALL MPI_ALLREDUCE(FIRST,MYFIRST,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)

c          IF(MAXZG>Z(1) .AND. MAXZG<=Z(NK-1)) THEN
c            LAST=MYRANK
c          ELSE
c            LAST=-1
c          ENDIF

          CALL MPI_ALLREDUCE(LAST,MYLAST,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)

c          write(6,*) 'myrank=',myrank,myfirst,mylast
c          call mpi_finalize(ierr)
c          stop

          CALL MPI_ALLREDUCE(IBMAX(IBD),IMAX,1,MPI_INTEGER,MPI_MAX,MPI_COMM_EDDY,IERR)
c          CALL MPI_ALLREDUCE(IBMIN(IBD),IMIN,1,MPI_INTEGER,MPI_MIN,MPI_COMM_EDDY,IERR)

          
          eps=1.0e-11
          rzmin=z(1)+eps
          rzmax=z(nk-1)
          if(myrank==myfirst) rzmin=-100.0

          CLOCK(2) = tclock()
          IF(ICOM==0 .AND. NB>0) CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 0, MYRANK, DS)
          CLOCK(2) = tclock()-CLOCK(2)
        
          IBS = 1

          J1 = 2
          J2 = NJ-1

          nins=0
          next=0

          IMIN = IBMIN(ibd)

          IF(NB>0) THEN

            NCROSS = 0
            ZCROSSO = 0.0
            DO I=IMIN,IMAX
              ZCROSS = 0.0
            
              CLOCKTEMP1 = tclock()

              DS = 0.5*S(I)**2.            
          
              DO J = J1,J2

                xcc = X(I,J)
                ycc = Y(I,J)
                QUERYXY(1) = XCC
                QUERYXY(2) = YCC
                ICROSS = 0
                CLOCKTEMP2 = tclock()
                CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 1, MYRANK, DS)
                CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2

                IF(NN>NNMAX) THEN 
                  clocktemp2 = tclock()
                  DEALLOCATE(INN)
                  NNMAX = NN
                  ALLOCATE(INN(NNMAX))
                  clock(12) = clock(12) + tclock() - clocktemp2
                  CLOCKTEMP2 = tclock()
                  CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 1, MYRANK, DS)
                  CLOCK(8) = CLOCK(8) + tclock() - CLOCKTEMP2
                ENDIF

                clocktemp3 = tclock()
                DO IN=1,MIN(NN,NNMAX)
                  N = INN(IN)
                  xa = xyzb(1,1,n)
                  ya = xyzb(2,1,n)
                  za = xyzb(3,1,n)
                  xb = xyzb(1,2,n)
                  yb = xyzb(2,2,n)
                  zb = xyzb(3,2,n)
                  xc = xyzb(1,3,n)
                  yc = xyzb(2,3,n)
                  zc = xyzb(3,3,n)
                  IF ((ycc.ge.min(ya,yb,yc).and.ycc.le.max(ya,yb,yc)) .and.  
     &              (xcc.ge.min(xa,xb,xc).and.xcc.le.max(xa,xb,xc)) ) THEN
                    CALL segtriint(xcc,ycc,rzmin,xcc,ycc,rzmax
     &               ,xa,ya,za,xb,yb,zb,xc,yc,zc,xm,ym,zm,iflag)
                    IF (iflag.eq.1) then
                      ICROSS = ICROSS+1
                      IF(ICROSS>NCR) THEN
                        WRITE(6,*) '1. WARNING MYRANK=',MYRANK,J,I,ICROSS,XCC,YCC,IBD
c                        do ii=1,ncr
c                          write(6,*) ii,zcross(ii,j)
c                        enddo
                      ENDIF
                      ZCROSS(ICROSS,J) = ZM

                    ENDIF

                  ENDIF
                ENDDO
                clock(9) = clock(9) + tclock() - clocktemp3

                NCROSS(I,J) = ICROSS
                IF(NCROSS(I,J)>NCR) THEN
                  WRITE(6,*) '2. WARNING MYRANK=',MYRANK,J,I,NCROSS(I,J),IBD
                ENDIF

                clocktemp3 = tclock()
                !Sort z-coord. crossings in ascending order
                DO I1=1,NCROSS(I,J)
                  ztmin = 1.0e10
                  do i2=1,NCROSS(I,J)
                    if (ZCROSS(i2,J).le.ztmin) then
                      ztmin = ZCROSS(i2,J)
                      isave = i2;
                    endif
                  enddo
                  ZCROSSO(i1,I,J) = ztmin
                  ZCROSS (isave,J) = 2.0e10
                ENDDO
                clock(10) = clock(10) + tclock() - clocktemp3

                clocktemp3 = tclock()
               !Remove duplicate z-coord. crossings
                zcross(:,j) = zcrosso(:,i,j)
                nc=0
                DO i1=1,ncross(i,j)
                  ndp=0
                  DO i2=1,ncross(i,j)
!!!!!!                    if(zcross(i2,j)==zcross(i1,j).AND.i1>i2) then
                    if(abs(zcross(i2,j)-zcross(i1,j)).lt.1.e-06.AND.i1>i2) then
                      ndp=ndp+1
                    endif
                  ENDDO

                  if(ndp==0) then
                    nc=nc+1
                    zcrosso(nc,i,j) = zcross(i1,j)
                  endif
                ENDDO
                clock(11) = clock(11) + tclock() - clocktemp3

                ncross(i,j) = nc
              ENDDO

              CLOCK(4) = CLOCK(4) + tclock() - CLOCKTEMP1

 
              CLOCKTEMP1 = tclock()

c              IF(MYRANK==MYFIRST .AND. MYFIRST/=MYLAST) THEN
c                CALL MPI_SEND(NCROSS,NY,MPI_INTEGER,JP,JP,MPI_COMM_EDDY,IERR)
c              ENDIF

              ncrprv=0

              DO J=2,NJ-1

                DO k=KBMIN(IBD),KBMAX(IBD)
                  zcc = z(k)
                  DO ICROSS=1,NCROSS(I,J)
                    IF (ZCROSSO(ICROSS,I,J) .LE. zcc) THEN
                      ibs(i,j,k) = -ibs(i,j,k) 
                    ENDIF
                  ENDDO

                  IF(ibs(i,j,k).eq.-1) nins=nins+1
                  IF(ibs(i,j,k).eq.1)  next=next+1
                ENDDO

              ENDDO

              IF(MYRANK>MYFIRST) THEN
                CALL MPI_RECV(NCRPRV,NJ,MPI_INTEGER,MYRANK-1,MYRANK-1,MPI_COMM_EDDY,STATUS,IERR)
              ENDIF

              IF(MYRANK<MYLAST) THEN
                CALL MPI_SEND(NCRPRV+NCROSS(I,:),NJ,MPI_INTEGER,MYRANK+1,MYRANK,MPI_COMM_EDDY,IERR)
              ENDIF

              IF(MYRANK>MYFIRST) THEN
                DO J=2,NJ-1
                  IF(MOD(NCRPRV(J),2)==1) THEN
c                  IBS(I,J,KBMIN(IBD):KBMAX(IBD))=-1
                    IBS(I,J,KBMIN(IBD):KBMAX(IBD))=-IBS(I,J,KBMIN(IBD):KBMAX(IBD))
                  ENDIF
                ENDDO
 
              ENDIF

            CLOCK(6) = CLOCK(6) + tclock() - CLOCKTEMP1

            ENDDO

          ENDIF

          CLOCK(14) = tclock()
          if(idomy==0) then
            IBS(IMIN:IBMAX(IBD),1,KBMIN(IBD):KBMAX(IBD))  = IBS(IMIN:IBMAX(IBD),NJ-1,KBMIN(IBD):KBMAX(IBD))
            IBS(IMIN:IBMAX(IBD),NJ,KBMIN(IBD):KBMAX(IBD)) = IBS(IMIN:IBMAX(IBD),2,KBMIN(IBD):KBMAX(IBD))
          endif
         
          call refreshflag(ibs,ni,nj,nk)

          CLOCK(14) = tclock() - clock(14)

          CLOCK(7) = tclock()
          IF(ICOM==2 .AND. NB>0) CALL COMMITTEE4(QUERYXY, XYZCB, NB, INN, NN, NNMAX, 2, 2, MYRANK, DS)
          clock(7) = tclock()-clock(7)

       ENDIF

       CLOCK(1) = tclock() - CLOCK(1)


       IF(IOLVL>0) THEN
         CALL MPI_REDUCE(NINS,NINSG,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_EDDY,IERR)

         IF(MYRANK.EQ.0) THEN
           OPEN(UNIT=16, FILE='stats_imb.dat', FORM='FORMATTED'
     &          ,POSITION='APPEND')
           N = ((NI-2)*(NJ-2)*(NK-2)*MYSIZE)
           WRITE(16,'(2(A,1X,I12),A,F6.2,A)') 'TAG3D: No. INSIDE='
     &          ,NINSG,'/',N,', (~',100.*REAL(NINSG)/REAL(N),'/100)'
           CLOSE(16)
         ENDIF
       ENDIF

       IF(IOCLOCK>0) THEN
         CALL MPI_REDUCE(CLOCK,CLOCKG,NCLOCKS,MTYPE,MPI_SUM,0,MPI_COMM_EDDY,IERR)
         CALL MPI_REDUCE(CLOCK,CLOCKGMIN,NCLOCKS,MTYPE,MPI_MIN,0,MPI_COMM_EDDY,IERR)
         CALL MPI_REDUCE(CLOCK,CLOCKGMAX,NCLOCKS,MTYPE,MPI_MAX,0,MPI_COMM_EDDY,IERR)
          
         CLOCKG = CLOCKG/REAL(MYSIZE)


         IF(MYRANK.EQ.0) THEN
c           WRITE(6,*) 'N1=',N1,REAL(N1)/REAL(N1+N2)*100,', N2=',N2,REAL(n2)/REAL(N1+N2)*100
c           WRITE(6,*) 'MAX NN=',NNMAX2
c           DO I=1,NNMAX
c             IF(NNTRI(I)>0) THEN
c               WRITE(6,*) I,NNTRI(I)
c             ENDIF
c           ENDDO
           OPEN(UNIT=16, FILE='clock.dat',FORM='FORMATTED'
     &           ,POSITION='APPEND')
           WRITE(16,'(A)') '---- TAG3D: ----------------------------'
           WRITE(16,'(A)') '   Task/Time            Ave.              Max.      
     &           Min.'
           WRITE(16,'(A,3(5x,F8.4,5x))') 'Total               :',CLOCKG(1),CLOCKGMAX(1),CLOCKGMIN(1)
           i=3
           write(16,905) 'Build query array  :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=2
           write(16,905) 'Build ANN structure:'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=4
           write(16,905) 'Find z crossings   :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=8
           write(16,905) ' -Find NN          :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=9
           write(16,905) ' -Find intrsctns.  :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=10
           write(16,905) ' -Order z crossings:'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=11
           write(16,905) ' -Remove duplic.   :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=12
           write(16,905) ' -Dealloc. & alloc.:'
     &          ,clockg(i),'(',100*clockg(i)/clockg(4),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(4),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(4),'%)'
           i=13
           write(16,905) 'MPI_ALLreduces.    :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=6
           write(16,905) 'Classify points    :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=14
           write(16,905) 'Refresh ibs y-dir  :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           i=7
           write(16,905) 'Delete ANN struct. :'
     &          ,clockg(i),'(',100*clockg(i)/clockg(1),'%)'
     &          ,clockgmax(i),'(',100*clockgmax(i)/clockgmax(1),'%)'
     &          ,clockgmin(i),'(',100*clockgmin(i)/clockgmin(1),'%)'
           WRITE(16,'(A)') '----------------------------------------'
           CLOSE(16)
                      
         ENDIF
       ENDIF 

       DEALLOCATE(CLOCK,CLOCKG,CLOCKGMIN,CLOCKGMAX)
       DEALLOCATE(INN)
       DEALLOCATE(NNTRI)

       RETURN

 905  format(A,3(1x,F8.4,A,F6.2,A,2x))

       END
C---------------------------------------------------------------------

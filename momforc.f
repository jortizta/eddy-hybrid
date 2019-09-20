C---- subroutine momforc------------------------N. Beratlis-25 Apr. 2009--
C
C     PURPOSE: Compute momentum forcing and correct velocity near immersed
C     body surface.
C
C-------------------------------------------------------------------------
      subroutine momforc(us,ub,xu,yu,zu,nx,ny,nz,iu,ju,ku,umtrx,uindx
     &     ,uim,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,impl,ivel)
      
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,impl,nflumax,ivel
      integer limu(nflu),mimu(nflu)
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax),mrku(nfcmax),diru(nfcmax)
      integer uindx(nsp,nfcmax)
      real    xu(nx),yu(ny),zu(nz)
      real    umtrx(nsp,nsp,nfcmax)
      real    nxu(nfcmax),nyu(nfcmax),nzu(nfcmax),uim(nfcmax)
      real    us(nx,ny,nz),ub(nx,ny,nz)
c
c.... Local var. and arrays
      integer im,ism,i,j,k,ilp,i1,j1,k1
      real    xint,yint,zint,xp,yp,zp
      integer indx(nsp)
      real    a(nsp,nsp),b(nsp)
c
C.... Functions
      real    interp_cellface,anglerad

      IF(impl==1) THEN
c implicit treatment

        DO ilp=1,nflumax

          DO im=limu(ilp)+1,limu(ilp)+mimu(ilp)
            a(:,:)=umtrx(:,:,im)
            indx(:)=uindx(:,im)
            b(1)=uim(im)
c b(1) is the velocity on the immersed-boundary
            if(mrku(im)==1 .OR. mrku(im)==3) then
c cases for which the interpolation stencil is found along the normal direction
c or using the closest point criterion: the projection point has been found as
c intersection between the outward direction and a grid face
c N.B. the explicit velocity UB is used for this interpolation
              b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &            ,ub,nx,ny,nz,icyl,diru(im))
c b(2) is the interpolated velocity at the extension point
            elseif(mrku(im)==2) then
c case for which the fluid point is found along a grid line
              do ism=1,nsm
                i1 = iu(im)+int(nxu(im))
                j1 = ju(im)+int(nyu(im))
                k1 = ku(im)+int(nzu(im))
                call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
                b(ism+1)=ub(i1,j1,k1)
              enddo
            endif

            call lubksb(a,nsp,nsp,indx,b)   !!!!!!
c the interpolated velocity at the interface point is stored in b(1)

            us(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))
     &           +  b(1)-ub(iu(im),ju(im),ku(im))
c the above row allows to eliminate one half of the terms treated by C.N.
c N.B. the velocity UB used to define the boundary condition is explicit

c the value of the boundary condition is stored in UB
            ub(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))

          ENDDO

c          CALL REFRESHBC(US,NX*NY,NZ)
          CALL REFRESHBC(UB,NX*NY,NZ)

          !Periodic y bc
          us(:,1 ,:) = us(:,ny-1,:)
          us(:,ny,:) = us(:,  2 ,:)

          !Centerline bc
          if(icyl==1) then
            if(ivel==1) then
              do j=1,ny
                us(1,j,:) = 0.5*(us(2,j,:)-us(2,jsym(j),:))    !!!!!! - instead of +
              enddo
            elseif(ivel==2) then
              do j=1,ny
                us(1,j,:) =-us(2,jsym(j),:)
              enddo
            elseif(ivel==3) then     !!!!!! 3 instead of 2
              do j=1,ny
                us(1,j,:) = us(2,jsym(j),:)
              enddo
            endif
          endif

        ENDDO

c the value of the boundary condition is stored in UB   !!!!!!
!!!!!!        DO ilp=1,nflumax
!!!!!!          DO im=limu(ilp)+1,limu(ilp)+mimu(ilp)
!!!!!!            ub(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))
!!!!!!          ENDDO
!!!!!!        ENDDO

        CALL REFRESHBC(US,NX*NY,NZ)
!!!!!!        CALL REFRESHBC(UB,NX*NY,NZ)

      ELSE
c explicit treatment

        DO ilp=1,nflumax     !!!!!! instead of nflu

          DO im=limu(ilp)+1,limu(ilp)+mimu(ilp)
            a(:,:)=umtrx(:,:,im)
            indx(:)=uindx(:,im)
            b(1)=uim(im)
            if(mrku(im)==1 .OR. mrku(im)==3) then
              b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &            ,us,nx,ny,nz,icyl,diru(im))
            elseif(mrku(im)==2) then
              do ism=1,nsm
                i1 = iu(im)+int(nxu(im))
                j1 = ju(im)+int(nyu(im))
                k1 = ku(im)+int(nzu(im))
                call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
                b(ism+1)=us(i1,j1,k1)
              enddo
            endif

            call lubksb(a,nsp,nsp,indx,b)

            us(iu(im),ju(im),ku(im))=b(1)

          ENDDO

          CALL REFRESHBC(US,NX*NY,NZ)

          !Periodic y bc
          us(:,1 ,:) = us(:,ny-1,:)
          us(:,ny,:) = us(:,  2 ,:)

          !Centerline bc
          if(icyl==1) then
            if(ivel==1) then
              do j=1,ny
                us(1,j,:) = 0.5*(us(2,j,:)-us(2,jsym(j),:))
              enddo
            elseif(ivel==2) then
              do j=1,ny
                us(1,j,:) =-us(2,jsym(j),:)
              enddo
            elseif(ivel==3) then     !!!!!! 3 instead of 2
              do j=1,ny
                us(1,j,:) = us(2,jsym(j),:)
              enddo
            endif
          endif

        ENDDO

      ENDIF


      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          us(1,j,:) = us(2,jsym(j),:)
c        enddo
c      endif


      return

      end
C-------------------------------------------------------------------------


C---- subroutine velinterior-------------------N. Beratlis-26 Apr. 2009---
C
C     PURPOSE: Set the velocity inside immersed body 
C
C-------------------------------------------------------------------------
      subroutine velinterior(us,vs,ws,ub,vb,wb,flagu,flagv,flagw,xu,yv,zw
     &     ,xc,yc,zc,nx,ny,nz,nbd,tlevel,ibd,mov,impl)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... input/output var. and arrays
      integer nx,ny,nz,ibd,nbd,mov,impl
      real    tlevel
      integer flagu(nx,ny,nz,nbd),flagv(nx,ny,nz,nbd),flagw(nx,ny,nz,nbd)
      real    xu(nx),xc(nx),yv(ny),yc(ny),zw(nz),zc(nz)
      real    us(nx,ny,nz),vs(nx,ny,nz),ws(nx,ny,nz)
      real    ub(nx,ny,nz),vb(nx,ny,nz),wb(nx,ny,nz)
c
C.... local var. and arrays
      integer i,j,k
c
c.... function
      real    ubd,vbd,wbd

      if(mov==0) then
c Stationary boundaries
        DO k=kbmin(ibd),kbmax(ibd)
        DO j=jbmin(ibd),jbmax(ibd)
        DO i=ibmin(ibd),ibmax(ibd)
          IF(flagu(i,j,k,ibd)==ibd) us(i,j,k) = 0.0
          IF(flagv(i,j,k,ibd)==ibd) vs(i,j,k) = 0.0
          if(flagw(i,j,k,ibd)==ibd) ws(i,j,k) = 0.0
        ENDDO
        ENDDO
        ENDDO
      else

c Moving boundaries
        if(icyl==0) then

          if(impl==0) then
c Case of explicit treatment
            DO k=kbmin(ibd),kbmax(ibd)
            DO j=jbmin(ibd),jbmax(ibd)
            DO i=ibmin(ibd),ibmax(ibd)
              if(flagu(i,j,k,ibd)==ibd) us(i,j,k) = ubd(xu(i),yc(j),zc(k),tlevel,ibd)
              if(flagv(i,j,k,ibd)==ibd) vs(i,j,k) = vbd(xc(i),yv(j),zc(k),tlevel,ibd)
              if(flagw(i,j,k,ibd)==ibd) ws(i,j,k) = wbd(xc(i),yc(j),zw(k),tlevel,ibd)
            ENDDO
            ENDDO
            ENDDO
          else
c Case of implicit treatment
            DO k=kbmin(ibd),kbmax(ibd)
            DO j=jbmin(ibd),jbmax(ibd)
            DO i=ibmin(ibd),ibmax(ibd)
              if(flagu(i,j,k,ibd)==ibd) us(i,j,k) = us(i,j,k)-ub(i,j,k)
     &           + ubd(xu(i),yc(j),zc(k),tlevel,ibd)
              if(flagv(i,j,k,ibd)==ibd) vs(i,j,k) = vs(i,j,k)-vb(i,j,k)
     &           + vbd(xc(i),yv(j),zc(k),tlevel,ibd)
              if(flagw(i,j,k,ibd)==ibd) ws(i,j,k) = ws(i,j,k)-wb(i,j,k)
     &           + wbd(xc(i),yc(j),zw(k),tlevel,ibd)
            ENDDO
            ENDDO
            ENDDO
          endif

        else

          if(impl==0) then
            DO k=kbmin(ibd),kbmax(ibd)
            DO j=jbmin(ibd),jbmax(ibd)
            DO i=ibmin(ibd),ibmax(ibd)
              if(flagu(i,j,k,ibd)==ibd) us(i,j,k) = ubd(xu(i)*cos(yc(j)),xu(i)*sin(yc(j)),zc(k),tlevel,ibd)
              if(flagv(i,j,k,ibd)==ibd) vs(i,j,k) = vbd(xc(i)*cos(yv(j)),xc(i)*sin(yv(j)),zc(k),tlevel,ibd)
              if(flagw(i,j,k,ibd)==ibd) ws(i,j,k) = wbd(xc(i)*cos(yc(j)),xc(i)*sin(yc(j)),zw(k),tlevel,ibd)
            ENDDO
            ENDDO
            ENDDO
          else
            DO k=kbmin(ibd),kbmax(ibd)
            DO j=jbmin(ibd),jbmax(ibd)
            DO i=ibmin(ibd),ibmax(ibd)
              if(flagu(i,j,k,ibd)==ibd) us(i,j,k) = us(i,j,k)-ub(i,j,k)
     &           + ubd(xu(i)*cos(yc(j)),xu(i)*sin(yc(j)),zc(k),tlevel,ibd)
              if(flagv(i,j,k,ibd)==ibd) vs(i,j,k) = vs(i,j,k)-vb(i,j,k)
     &           + vbd(xc(i)*cos(yv(j)),xc(i)*sin(yv(j)),zc(k),tlevel,ibd)
              if(flagw(i,j,k,ibd)==ibd) ws(i,j,k) = ws(i,j,k)-wb(i,j,k)
     &           + wbd(xc(i)*cos(yc(j)),xc(i)*sin(yc(j)),zw(k),tlevel,ibd)
            ENDDO
            ENDDO
            ENDDO
          endif

        endif
      endif

      us(:,1,:)  = us(:,ny-1,:)
      us(:,ny,:) = us(:,2,:)
      vs(:,1,:)  = vs(:,ny-1,:)
      vs(:,ny,:) = vs(:,2,:)
      ws(:,1,:)  = ws(:,ny-1,:)
      ws(:,ny,:) = ws(:,2,:)

      CALL REFRESHBC(US,NX*NY,NZ)
      CALL REFRESHBC(VS,NX*NY,NZ)
      CALL REFRESHBC(WS,NX*NY,NZ)

      return

      end
C-------------------------------------------------------------------------


C---- subroutine tvinterior-------------------N. Beratlis-26 Apr. 2009---
C
C     PURPOSE: Set the eddy viscosity inside immersed body 
C
C-------------------------------------------------------------------------
      subroutine tvinterior(tv,flagtv,nx,ny,nz,nbd,ibd)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... input/output var. and arrays
      integer nx,ny,nz,nbd,ibd
      integer flagtv(nx,ny,nz,nbd)
      real    tv(nx,ny,nz)
c
C.... local var. and arrays
      integer i,j,k

      DO k=kbmin(ibd),kbmax(ibd)
      DO j=jbmin(ibd),jbmax(ibd)
      DO i=ibmin(ibd),ibmax(ibd)
        IF(flagtv(i,j,k,ibd)==ibd) tv(i,j,k) = 0.0
      ENDDO
      ENDDO
      ENDDO

      tv(:,1,:)  = tv(:,ny-1,:)
      tv(:,ny,:) = tv(:,2,:)

      CALL REFRESHBC(TV,NX*NY,NZ)

      return

      end
C-------------------------------------------------------------------------



c---- subroutine correctpres----------------N. Beratlis-18 May 2009---
c------------------------------------------- Modified by A. Posa -----
C
C     PURPOSE: Correct pressure
C
C---------------------------------------------------------------------
      subroutine correctpres(p,nx,ny,nz,ip,jp,kp,pmtrx
     &     ,pim,nxp,nyp,nzp,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp
     &     ,nfacet,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax,tlevel,dt,ibd,mbd)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nfacet,ibd,mbd,nflpmax
      integer limp(nflp),mimp(nflp)
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      integer mrkp(nfcmax),dirp(nfcmax)
      integer pindx(nsp,nfcmax)
      integer fp(nfcmax)
      real    tlevel,dt
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    pim(nfcmax),xnp(nfcmax),ynp(nfcmax),znp(nfcmax),dpdnn(nfcmax)
      real    unvect(3,nfacet)
      real    nxp(nfcmax),nyp(nfcmax),nzp(nfcmax)
      real    pmtrx(nsp,nsp,nfcmax)
c
c.... Local arrays
      integer i,j,ism,ii,jj,kk,in,jn,kn,in1,jn1,kn1,in2,jn2,kn2,ilp
      integer i1,j1,k1
      integer i2,j2,k2
      integer i3,j3,k3
      integer indx(nsp)
      real    xpcar,ypcar,zpcar,s1,s2,xp1car,yp1car,zp1car,xp2car,yp2car,zp2car
      real    v1(3),v2(3)
      real    a(nsp,nsp),b(nsp)
c
c.... Function
      real    interp_cellface
      real    dotproduct,anglerad

c The normal pressure gradient is equal to 0 in the case of stationary
c boundaries
      dpdnn(limp(1)+1:limp(1)+sum(mimp(:))) = 0.0   !!!!!! instead of dpdnn = 0.0

      do ilp=1,nflpmax


        if(ibd>=mbd) then
c For the moving boundaries the normal pressure gradient is evaluated
          call presgrad_flagp(dpdnn,mrkp,xnp,ynp,znp,nxp,nyp,nzp,ip,jp,kp
     &          ,xc_car,yc_car,zc,nx,ny,nz,unvect,fp
     &          ,nfacet,limp(ilp),mimp(ilp),tlevel,dt,ibd)
        endif


        do i=limp(ilp)+1,limp(ilp)+mimp(ilp)

          a(:,:) = pmtrx(:,:,i)
          indx(:) = pindx(:,i)

!	      in = ip(i)	  !!!!!!!
!	      jn = jp(i)	  !!!!!!!
!	      kn = kp(i)	  !!!!!!!

          b(1) = dpdnn(i)

          if(mrkp(i)==1 .OR. mrkp(i)==3) then
c the pressure at the extension fluid point is evaluated by interpolation
            pim(i) = interp_cellface(nxp(i),nyp(i),nzp(i)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dirp(i))
          else
            i1 = ip(i)+int(nxp(i))
            j1 = jp(i)+int(nyp(i))
            k1 = kp(i)+int(nzp(i))
            call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
c the pressure at the exterior grid point is determined
            pim(i) = p(i1,j1,k1)
          endif
          b(2) = pim(i)

c the local system for the pressure stencil is solved
          call lubksb(a,nsp,nsp,indx,b)
c the pressure condition at the interface point is enforced
          p(ip(i),jp(i),kp(i))=b(1)

        enddo


        call refreshbc(p,nx*ny,nz)
        
        !Periodic y bc
        p(:,1,:) = p(:,ny-1,:)
        p(:,ny,:) = p(:,2,:)

        !Centerline bc
        if(icyl==1) then
          do j=1,ny
            p(1,j,:) = p(2,jsym(j),:)
          enddo
        endif

      enddo

      !Periodic y bc
c      p(:,1,:) = p(:,ny-1,:)
c      p(:,ny,:) = p(:,2,:)

      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          p(1,j,:) = p(2,jsym(j),:)
c        enddo
c      endif

      return

      end
c---------------------------------------------------------------------



C---- subroutine momforc1-----------------------N. Beratlis-25 Apr. 2009--
C
C     PURPOSE: Compute momentum forcing and correct velocity near immersed
C     body surface.
C
C-------------------------------------------------------------------------
      subroutine momforc1(us,ub,xu,yu,zu,nx,ny,nz,iu,ju,ku,iumtrx,jumtrx,kumtrx
     &       ,umtrx,uindx,uim,mrku,diru,nxu,nyu,nzu,xnu,ynu,znu,limu,mimu,impl)
      
c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,impl,limu,mimu
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax),mrku(nfcmax),diru(nfcmax)
      integer uindx(nsp,nfcmax)
      integer iumtrx(nsm,nfcmax),jumtrx(nsm,nfcmax),kumtrx(nsm,nfcmax)
      real    xu(nx),yu(ny),zu(nz)
      real    umtrx(nsp,nsp,nfcmax)
      real    nxu(nfcmax),nyu(nfcmax),nzu(nfcmax),uim(nfcmax)
      real    xnu(nfcmax),ynu(nfcmax),znu(nfcmax)
      real    us(nx,ny,nz),ub(nx,ny,nz)
c
c.... Local var. and arrays
      integer im,ism,i,j,k
      real    s1,s2,u1,u2,a1
      integer indx(nsp)
      real    a(nsp,nsp),b(nsp)
c
C.... Functions
      real    interp_cellface

      IF(impl==1) THEN

        DO im=limu+1,mimu
          a(:,:)=umtrx(:,:,im)
          indx(:)=uindx(:,im)
          b(1)=uim(im)
          if(mrku(im)==1) then
            b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &          ,ub,nx,ny,nz,icyl,diru(im))
          elseif(mrku(im)==2) then
            do ism=1,nsm
              b(ism+1)=ub(iumtrx(ism,im),jumtrx(ism,im),kumtrx(ism,im))
            enddo
          endif
          call lubksb(a,nsp,nsp,indx,b)
          us(iu(im),ju(im),ku(im))=us(iu(im),ju(im),ku(im))
     &         +  b(1)-ub(iu(im),ju(im),ku(im))
        ENDDO

      ELSE

        DO im=limu+1,mimu
          a(:,:)=umtrx(:,:,im)
          indx(:)=uindx(:,im)
          b(1)=uim(im)
          if(mrku(im)==1) then
            b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &          ,us,nx,ny,nz,icyl,diru(im))
          elseif(mrku(im)==2) then
            do ism=1,nsm
              b(ism+1)=us(iumtrx(ism,im),jumtrx(ism,im),kumtrx(ism,im))
            enddo
          endif

          u2 = b(2)
          call lubksb(a,nsp,nsp,indx,b)
          us(iu(im),ju(im),ku(im))=b(1)

          if(icyl==0) then
            s1 = sqrt( (xnu(im)-xu(iu(im)))**2.
     &            +  (ynu(im)-yu(ju(im)))**2.
     &            +  (znu(im)-zu(ku(im)))**2.)
            if(mrku(im)==1) then
              s2 = sqrt( (xnu(im)-nxu(im))**2.
     &              + (ynu(im)-nyu(im))**2.
     &              + (znu(im)-nzu(im))**2. )
            elseif(mrku(im)==2) then
              s2 = sqrt( (xnu(im)-xu(iumtrx(1,im)))**2.
     &              + (ynu(im)-yu(jumtrx(1,im)))**2.
     &              + (znu(im)-zu(kumtrx(1,im)))**2. )
            endif
          else
            s1 = sqrt( (xnu(im)-xu(iu(im))*cos(yu(ju(im))))**2.
     &            +  (ynu(im)-xu(iu(im))*sin(yu(ju(im))) )**2.
     &            +  (znu(im)-zu(ku(im)))**2.)
            if(mrku(im)==1) then
              s2 = sqrt( (nxu(im)-xu(iu(im))*cos(yu(ju(im))))**2.
     &              + (nyu(im)-xu(iu(im))*sin(yu(ju(im))))**2.
     &              + (nzu(im)-zu(ku(im)))**2. )
            elseif(mrku(im)==2) then
              s2 = sqrt( (xu(iu(im))*cos(yu(ju(im)))-xu(iumtrx(1,im))*cos(yu(jumtrx(1,im))))**2.
     &              + (xu(iu(im))*sin(yu(ju(im)))-xu(iumtrx(1,im))*sin(yu(jumtrx(1,im))))**2.
     &              + (zu(ku(im))-zu(kumtrx(1,im)))**2. )              
            endif
          endif

          a1=abs(s2/s1)/(abs(u2-b(1))/abs(b(1)-uim(im)))
        ENDDO

      ENDIF

      CALL REFRESHBC(US,NX*NY,NZ)

      return

      end
C-------------------------------------------------------------------------

      subroutine densinterior(dens,nx,ny,nz,flagp,ibd,nbd)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,ibd,nbd
      integer flagp(nx,ny,nz,nbd)
      real    dens(nx,ny,nz)
      
      integer i,j,k

c      write(6,*) 'presinterior:',myrank,ibmin(ibd),ibmax(ibd),jbmin(ibd),jbmax(ibd),kbmin(ibd),kbmax(ibd)
c      write(6,*) 'flagpo:',myrank,flagp(4,2,126,ibd)
      do k=kbmin(ibd),kbmax(ibd)
      do j=jbmin(ibd),jbmax(ibd)
      do i=ibmin(ibd),ibmax(ibd)
        if((flagp(i,j,k,ibd)==ibd)) dens(i,j,k) = 0.0
c        if(i==4 .AND. k+myrank*(nz-2)==126) write(6,*) 'myrank=',myrank,flagp(i,j,k,ibd),p(i,j,k)
      enddo
      enddo
      enddo

      CALL REFRESHBC(DENS,NX*NY,NZ)
      
      return

      end


C---- subroutine velinterior1 ---------------N. Beratlis-26 Apr. 2009---
C
C     PURPOSE: Set the velocity inside immersed body 
C
C-----------------------------------------------------------------------
      subroutine velinterior1(us,vs,ws,ub,vb,wb,iu,ju,ku,iv,jv,kv,iw,jw,kw
     &     ,xu,yv,zw,xc,yc,zc,nx,ny,nz,limu,mimu,limv,mimv,limw,mimw,tlevel,ibd,mov,impl)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... input/output var. and arrays
      integer nx,ny,nz,ibd,nbd,mov,impl
      integer limu,limv,limw,mimu,mimv,mimw
      real    tlevel
      integer iu(nintmax),ju(nintmax),ku(nintmax)
      integer iv(nintmax),jv(nintmax),kv(nintmax)
      integer iw(nintmax),jw(nintmax),kw(nintmax)
      real    xu(nx),xc(nx),yv(ny),yc(ny),zw(nz),zc(nz)
      real    us(nx,ny,nz),vs(nx,ny,nz),ws(nx,ny,nz)
      real    ub(nx,ny,nz),vb(nx,ny,nz),wb(nx,ny,nz)
c
C.... local var. and arrays
      integer i,j,k,im
c
c.... function
      real    ubd,vbd,wbd

      if(mov==0) then
        do i=limu+1,limu+mimu
          us(iu(i),ju(i),ku(i))=0.0
        enddo
        do i=limv+1,limv+mimv
          vs(iv(i),jv(i),kv(i))=0.0
        enddo
        do i=limw+1,limw+mimw
          ws(iw(i),jw(i),kw(i))=0.0
        enddo
      else

        if(icyl==0) then

          if(impl==0) then
            do im=limu+1,limu+mimu
              i=iu(im)
              j=ju(im)
              k=ku(im)
              us(i,j,k) = ubd(xu(i),yc(j),zc(k),tlevel,ibd)
            enddo
            do im=limv+1,limv+mimv
              i=iv(im)
              j=jv(im)
              k=kv(im)
              vs(i,j,k) = vbd(xc(i),yv(j),zc(k),tlevel,ibd)
            enddo
            do im=limw+1,limw+mimw
              i=iw(im)
              j=jw(im)
              k=kw(im)
              ws(i,j,k) = wbd(xc(i),yc(j),zw(k),tlevel,ibd)
            enddo
          else
            do im=limu+1,limu+mimu
              i=iu(im)
              j=ju(im)
              k=ku(im)
              us(i,j,k) = us(i,j,k)-ub(i,j,k)
     &           + ubd(xu(i),yc(j),zc(k),tlevel,ibd)
            enddo
            do im=limv+1,limv+mimv
              i=iv(im)
              j=jv(im)
              k=kv(im)
              vs(i,j,k) = vs(i,j,k)-vb(i,j,k)
     &           + vbd(xc(i),yv(j),zc(k),tlevel,ibd)
            enddo
            do im=limw+1,limw+mimw
              i=iw(im)
              j=jw(im)
              k=kw(im)
              ws(i,j,k) = ws(i,j,k)-wb(i,j,k)
     &           + wbd(xc(i),yc(j),zw(k),tlevel,ibd)
            enddo

          endif

        else

          if(impl==0) then

            do im=limu+1,limu+mimu
              i=iu(im)
              j=ju(im)
              k=ku(im)
              us(i,j,k) = ubd(xu(i)*cos(yc(j)),xu(i)*sin(yc(j)),zc(k),tlevel,ibd)
            enddo
            do im=limv+1,limv+mimv
              i=iv(im)
              j=jv(im)
              k=kv(im)
              vs(i,j,k) = vbd(xc(i)*cos(yv(j)),xc(i)*sin(yv(j)),zc(k),tlevel,ibd)
            enddo
            do im=limw+1,limw+mimw
              i=iw(im)
              j=jw(im)
              k=kw(im)
              ws(i,j,k) = wbd(xc(i)*cos(yc(j)),xc(i)*sin(yc(j)),zw(k),tlevel,ibd)
            enddo

          else

            do im=limu+1,limu+mimu
              i=iu(im)
              j=ju(im)
              k=ku(im)
              us(i,j,k) = us(i,j,k)-ub(i,j,k)
     &           + ubd(xu(i)*cos(yc(j)),xu(i)*sin(yc(j)),zc(k),tlevel,ibd)
            enddo
            do im=limv+1,limv+mimv
              i=iv(im)
              j=jv(im)
              k=kv(im)
              vs(i,j,k) = vs(i,j,k)-vb(i,j,k)
     &           + vbd(xc(i)*cos(yv(j)),xc(i)*sin(yv(j)),zc(k),tlevel,ibd)
            enddo
            do im=limw+1,limw+mimw
              i=iw(im)
              j=jw(im)
              k=kw(im)
              ws(i,j,k) = ws(i,j,k)-wb(i,j,k)
     &           + wbd(xc(i)*cos(yc(j)),xc(i)*sin(yc(j)),zw(k),tlevel,ibd)
            enddo

          endif

        endif
      endif

      us(:,1,:)  = us(:,ny-1,:)
      us(:,ny,:) = us(:,2,:)
      vs(:,1,:)  = vs(:,ny-1,:)
      vs(:,ny,:) = vs(:,2,:)
      ws(:,1,:)  = ws(:,ny-1,:)
      ws(:,ny,:) = ws(:,2,:)

      CALL REFRESHBC(US,NX*NY,NZ)
      CALL REFRESHBC(VS,NX*NY,NZ)
      CALL REFRESHBC(WS,NX*NY,NZ)

      return

      end
C-----------------------------------------------------------------------


C---- subroutine tvinterior1 ----------------N. Beratlis-26 Apr. 2009---
C
C     PURPOSE: Set the eddy viscosity inside immersed body 
C
C-----------------------------------------------------------------------
      subroutine tvinterior1(tv,itv,jtv,ktv,xc,yc,zc,nx,ny,nz,limtv,mimtv)   !!!!!! ibd

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... input/output var. and arrays
      integer nx,ny,nz    !!!!!! ,ibd
      integer limtv,mimtv
      integer itv(nintmax),jtv(nintmax),ktv(nintmax)
      real    xc(nx),yc(ny),zc(nz)
      real    tv(nx,ny,nz)
c
C.... local var. and arrays
      integer i,j,k,im

      do im=limtv+1,limtv+mimtv
        i=itv(im)
        j=jtv(im)
        k=ktv(im)
        tv(i,j,k) = 0.0
      enddo

      tv(:,1,:) = tv(:,ny-1,:)
      tv(:,ny,:) = tv(:,2,:)

      CALL REFRESHBC(TV,NX*NY,NZ)

      return

      end
C-----------------------------------------------------------------------



C---- subroutine presinterior1---------------N. Beratlis-18 May 2009--
C
C     PURPOSE: Set the pressure in the interior to zero
C
c---------------------------------------------------------------------
      subroutine presinterior1(p,nx,ny,nz,ip,jp,kp,limp,mimp)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,limp,mimp
      integer ip(nintmax),jp(nintmax),kp(nintmax)
      real    p(nx,ny,nz)
      
      integer i

      do i=limp+1,limp+mimp
        p(ip(i),jp(i),kp(i))=0
      enddo

      CALL REFRESHBC(P,NX*NY,NZ)
      
      return

      end


      subroutine densinterior1(dens,nx,ny,nz,ip,jp,kp,limp,mimp)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,limp,mimp
      integer ip(nintmax),jp(nintmax),kp(nintmax)
      real    dens(nx,ny,nz)
      
      integer i

      do i=limp+1,limp+mimp
        dens(ip(i),jp(i),kp(i))=0.0
      enddo

      CALL REFRESHBC(DENS,NX*NY,NZ)
      
      return

      end
c---------------------------------------------------------------------

C---- subroutine convert_velcyl2car -----------N. Beratlis-14 Nov.
C2010---
C
C     PURPOSE: Convert velocity from cylindrical coordinates to
C     cartesian
C
C-------------------------------------------------------------------------
      subroutine convert_velcyl2car(ur,ut,ux,uy,x,y)
c
      implicit none
c
c.... Input/Output arrays
      real x,y,ur,ut,ux,uy
c
c.... Local arrays
      real u1,u2,theta
c
c.... Functions
      real anglerad

      theta = anglerad(x,y)

      u1 = ur*cos(theta) - uy*sin(theta)
      u2 = ur*sin(theta) + uy*cos(theta)

      ux = u1
      uy = u2

      return

      end
C-------------------------------------------------------------------------


C---- subroutine pfilter-----------------------N. Beratlis-18 Nov. 2010---
C
C     PURPOSE: Smooth pressure
C
C-------------------------------------------------------------------------
      subroutine pfilter(p,dp,flagp,nx,ny,nz,nbd,ip,jp,kp,limp,mimp
     $     ,nflpmax,ibd,npass)
C
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nbd,ibd,nflpmax,npass
      integer limp(nflp),mimp(nflp)
      integer flagp(nx,ny,nz,nbd)
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      real    xc(nx),yc(ny),zc(nz)
      real    p(nx,ny,nz),dp(nx,ny,nz)
c
c.... Local arrays
      integer i,j,k,il,ii,im
      real    dsum,psum,d

      do ii=1,npass
        
        dp = p

        do il=1,nflpmax

          do im=limp(ibd)+1,limp(ibd)+mimp(ibd)
            i = ip(im)
            j = jp(im)
            k = kp(im)

            psum = 0.0
            dsum = 0.0

            if(flagp(i+1,j,k,ibd)<=0) then
              d = au(i)
              dsum = dsum + d
              psum = psum + d*dp(i+1,j,k)
            endif

            if(flagp(i-1,j,k,ibd)<=0) then
              d = au(i-1)
              dsum = dsum + d
              psum = psum + d*dp(i-1,j,k)
            endif

c            if(flagp(i,j+1,k,ibd)<=0) then
c              d = bv(j)/rp(i)
c              dsum = dsum + d
c              psum = psum + d*dp(i,j+1,k)
c            endif

c            if(flagp(i,j-1,k,ibd)<=0) then
c              d = bv(j-1)/rp(i)            
c              dsum = dsum + d
c              psum = psum + d*dp(i,j-1,k)
c            endif

            if(flagp(i,j,k+1,ibd)<=0) then
              d = cw(k)
              dsum = dsum + d
              psum = psum + d*dp(i,j,k+1)
            endif

            if(flagp(i,j,k-1,ibd)<=0) then
              d = cw(k-1)
              dsum = dsum + d
              psum = psum + d*dp(i,j,k-1)
            endif

            if(dsum>0.0) then
              p(i,j,k) = (dp(i,j,k)*dsum + psum)/(2.0*dsum)
            endif

          enddo

          call refreshbc(p,nx*ny,nz)
         
          !Periodic y bc
          p(:,1,:) = p(:,ny-1,:)
          p(:,ny,:) = p(:,2,:)

          !Centerline bc
          if(icyl==1) then
            do j=1,ny
              p(1,j,:) = p(2,jsym(j),:)
            enddo
          endif

        enddo

      enddo

      return

      end
C-------------------------------------------------------------------------



C---- subroutine extend_zgrid------------------N. Beratlis-20 Nov. 2010---
C
C     PURPOSE: Extend local zgrid
C
C-------------------------------------------------------------------------
      subroutine extend_zgrid(zwg,zcg,nzg,zw,zc,nzl)
c
c      implicit none
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output array
      integer nzl,nzg
      real    zw(nzl),zc(nzl),zwg(nzg),zcg(nzg)
c
c.... Local arrays
      integer k,kl,nz,ext,kg

      nz = (nzg-2)/mysize + 2
      ext = (nzl-nz)/2

      do k=1,nzl
        kg = myrank*(nz-2)+k-ext
c        write(myrank+10,*) myrank,k,kg,ext,nz,nzl
        if(kg<=0) then
          zw(k) = zwg(1)-(zwg(2)-zwg(1))*abs(kg)
          zc(k) = zcg(1)-(zcg(2)-zcg(1))*abs(kg)
        elseif(kg>nzg) then
          zw(k) = zwg(nzg)+(zwg(nzg)-zwg(nzg-1))*(kg-nzg)
          zc(k) = zcg(nzg)+(zcg(nzg)-zcg(nzg-1))*(kg-nzg)
        else
          zw(k) = zwg(kg)
          zc(k) = zcg(kg)
        endif
      enddo

      return

      end
C-------------------------------------------------------------------------


C---- subroutine extend_flag_z-----------------N. Beratlis-20 Nov. 2010---
C
C     PURPOSE: Extend flag array in z directions
C
C-------------------------------------------------------------------------
      subroutine extend_flag_z(flag1,flag2,nx,ny,nzl,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nzl
      integer flag1(nx,ny,nzl),flag2(nx,ny,nz)
c
c.... Local arrays
      integer k,kg,ext,nzg
      integer status(MPI_STATUS_SIZE)

      ext = (nz-nzl)/2
      nzg = mysize*(nzl-2)+2

      do k=1,nzl
        flag2(:,:,k+ext) = flag1(:,:,k)
      enddo

      DO k=1,ext
        CALL MPI_SENDRECV(FLAG1(1,1,NZL-EXT-2+K),NX*NY,MPI_INTEGER,MYRIGHT,10
     &         ,          FLAG2(1,1,K),NX*NY,MPI_INTEGER,MYLEFT ,10
     &         ,            MPI_COMM_EDDY,STATUS,IERR)

        CALL MPI_SENDRECV(FLAG1(1,1,2+K),NX*NY,MPI_INTEGER,MYLEFT ,11
     &       ,            FLAG2(1,1,NZ-EXT+K),NX*NY,MPI_INTEGER,MYRIGHT,11
     &       ,            MPI_COMM_EDDY,STATUS,IERR)
      ENDDO


      return

      end
C-------------------------------------------------------------------------



C---- subroutine extend_var_z------------------N. Beratlis-20 Nov. 2010---
C
C     PURPOSE: Extend flag array in z directions
C
C-------------------------------------------------------------------------
      subroutine extend_var_z(var1,var2,nx,ny,nzl,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nzl
      real    var1(nx,ny,nzl),var2(nx,ny,nz)
c
c.... Local arrays
      integer k,kg,ext,nzg
      integer status(MPI_STATUS_SIZE)

      ext = (nz-nzl)/2
      nzg = mysize*(nzl-2)+2

      do k=1,nzl
        var2(:,:,k+ext) = var1(:,:,k)
      enddo

      DO k=1,ext
c        write(6,*) 'myrank=',myrank,k,', send',nzl-ext-2+k,',recv at',k
c     $        ,', send',2+k,',rcv at',nzl+k
        CALL MPI_SENDRECV(VAR1(1,1,NZL-EXT-2+K),NX*NY,MTYPE,MYRIGHT,10
     &         ,          VAR2(1,1,K),NX*NY,MTYPE,MYLEFT ,10
     &         ,            MPI_COMM_EDDY,STATUS,IERR)
c        write(6,*) 'var: myrank=',myrank,var1(109,4,128),var2(109,4,2)

        CALL MPI_SENDRECV(VAR1(1,1,2+K),NX*NY,MTYPE,MYLEFT ,11
     &       ,            VAR2(1,1,NZ-EXT+K),NX*NY,MTYPE,MYRIGHT,11
     &       ,            MPI_COMM_EDDY,STATUS,IERR)
      ENDDO


      return

      end
C-------------------------------------------------------------------------


C---- subroutine mpi_var_interp ---------------N. Beratlis-30 Nov.2010---
C
C     PURPOSE: Interpolate and communicate variable from left and right
C     proc.
C
C-------------------------------------------------------------------------
      subroutine mpi_var_interp(cordl,ureql,nl,cordr,ureqr,nr,u,xu,yu,zu
     $     ,nx,ny,nz)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer n,nx,ny,nz,nl,nr
      real    xu(nx),yu(ny),zu(nz)
      real    cordl(3,nl),cordr(3,nr),ureql(nl),ureqr(nr)
      real    u(nx,ny,nz)
c
c.... Local arrays
      integer nl_req,nl_sup,nr_req,nr_sup,nmsg,ii
      real, dimension(:,:), allocatable :: cord
      real, dimension(:), allocatable :: u_sup
      integer status(mpi_status_size)
c
c.... Functions
      real    interp3d

      nmsg = 10

      nl_req = nl
      nr_req = nr

      CALL MPI_SENDRECV(NL_REQ,1,MPI_INTEGER,MYLEFT,1
     &                   ,NR_SUP,1,MPI_INTEGER,MYRIGHT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

      CALL MPI_SENDRECV(NR_REQ,1,MPI_INTEGER,MYRIGHT,1
     &                   ,NL_SUP,1,MPI_INTEGER,MYLEFT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

      IF(MOD(MYRANK,2)==0) THEN

        IF(NL_REQ>0) THEN
          CALL MPI_SEND(CORDL,3*NL_REQ,MTYPE,MYLEFT,2,MPI_COMM_EDDY,IERR)
          CALL MPI_RECV(UREQL,NL_REQ,MTYPE,MYLEFT,6,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

        IF(NR_SUP>0) THEN
          ALLOCATE(CORD(3,NR_SUP),U_SUP(NR_SUP))
          CALL MPI_RECV(CORD,NR_SUP*3,MTYPE,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NR_SUP
            u_sup(ii) = interp3d(cord(1,ii),cord(2,ii),cord(3,ii),xu,yu,zu,u,nx,ny,nz,icyl)
          ENDDO
          CALL MPI_SEND(U_SUP,NR_SUP,MTYPE,MYRIGHT,6,MPI_COMM_EDDY,IERR)
          DEALLOCATE(CORD,U_SUP)
        ENDIF

      ELSE

        IF(NR_SUP>0) THEN
          ALLOCATE(CORD(3,NR_SUP),U_SUP(NR_SUP))
          CALL MPI_RECV(CORD,NR_SUP*3,MTYPE,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NR_SUP
            u_sup(ii) = interp3d(cord(1,ii),cord(2,ii),cord(3,ii),xu,yu,zu,u,nx,ny,nz,icyl)
          ENDDO
          CALL MPI_SEND(U_SUP,NR_SUP,MTYPE,MYRIGHT,6,MPI_COMM_EDDY,IERR)
          DEALLOCATE(CORD,U_SUP)
        ENDIF

        IF(NL_REQ>0) THEN
          CALL MPI_SEND(CORDL,3*NL_REQ,MTYPE,MYLEFT,2,MPI_COMM_EDDY,IERR)
          CALL MPI_RECV(UREQL,NL_REQ,MTYPE,MYLEFT,6,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

      ENDIF

      IF(MOD(MYRANK,2)==0) THEN

        IF(NR_REQ>0) THEN
          CALL MPI_SEND(CORDR,3*NR_REQ,MTYPE,MYRIGHT,3,MPI_COMM_EDDY,IERR)
          CALL MPI_RECV(UREQR,NR_REQ,MTYPE,MYRIGHT,5,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

        IF(NL_SUP>0) THEN
          ALLOCATE(CORD(3,NL_SUP),U_SUP(NL_SUP))
          CALL MPI_RECV(CORD,NL_SUP*3,MTYPE,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NL_SUP
            u_sup(ii) = interp3d(cord(1,ii),cord(2,ii),cord(3,ii),xu,yu,zu,u,nx,ny,nz,icyl)
          ENDDO
          CALL MPI_SEND(U_SUP,NL_SUP,MTYPE,MYLEFT,5,MPI_COMM_EDDY,IERR)
          DEALLOCATE(CORD,U_SUP)
        ENDIF

      ELSE

        IF(NL_SUP>0) THEN
          ALLOCATE(CORD(3,NL_SUP),U_SUP(NL_SUP))
          CALL MPI_RECV(CORD,NL_SUP*3,MTYPE,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NL_SUP
            u_sup(ii) = interp3d(cord(1,ii),cord(2,ii),cord(3,ii),xu,yu,zu,u,nx,ny,nz,icyl)
          ENDDO
          CALL MPI_SEND(U_SUP,NL_SUP,MTYPE,MYLEFT,5,MPI_COMM_EDDY,IERR)
          DEALLOCATE(CORD,U_SUP)
        ENDIF

        IF(NR_REQ>0) THEN
          CALL MPI_SEND(CORDR,3*NR_REQ,MTYPE,MYRIGHT,3,MPI_COMM_EDDY,IERR)
          CALL MPI_RECV(UREQR,NR_REQ,MTYPE,MYRIGHT,5,MPI_COMM_EDDY,STATUS,IERR)
        ENDIF

      ENDIF

      return

      end
C-------------------------------------------------------------------------



C---- subroutine momforc_mod2-------------------N. Beratlis-25 Apr. 2009--
C--------------------------------------------------A. Posa - 02 Nov 2012--
C
C     PURPOSE: Compute momentum forcing and correct velocity near immersed
C     body surface.
C
C-------------------------------------------------------------------------
      subroutine momforc_mod2(us,ub,xu,yu,zu,nx,ny,nz,iu,ju,ku,umtrx,uindx
     &     ,uim,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,impl,ivel,nbd,nmax)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,impl,nflumax(nbd),ivel,nbd
      integer limu(nbd,nflu),mimu(nbd,nflu)
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax),mrku(nfcmax),diru(nfcmax)
      integer uindx(nsp,nfcmax)
      real    xu(nx),yu(ny),zu(nz)
      real    umtrx(nsp,nsp,nfcmax)
      real    nxu(nfcmax),nyu(nfcmax),nzu(nfcmax),uim(nfcmax)
      real    us(nx,ny,nz),ub(nx,ny,nz)
      integer nmax
c
c.... Local var. and arrays
      integer im,ism,i,j,k,ilp,i1,j1,k1,nloopmax,ibd,nloopg
      real    xint,yint,zint,xp,yp,zp
      integer indx(nsp)
      real    a(nsp,nsp),b(nsp)
      integer ilu,iru,indexl(nmax,2),indexr(nmax,2)
      integer ilut,irut,indexlt(nmax,2),indexrt(nmax,2)
c
C.... Functions
      real    interp_cellface,anglerad

      nloopmax = maxval(nflumax(1:nbd))
      CALL MPI_ALLREDUCE(nloopmax,nloopg,1,mpi_integer,MPI_MAX,MPI_COMM_EDDY,IERR)

      IF(impl==1) THEN
c implicit treatment

        ilut = 0
        irut = 0

        DO ilp=1,nloopg

          ilu = 0
          iru = 0

          DO ibd = 1, nbd

            if(ilp.gt.nflumax(ibd)) cycle

            DO im=limu(ibd,ilp)+1,limu(ibd,ilp)+mimu(ibd,ilp)
              a(:,:)=umtrx(:,:,im)
              indx(:)=uindx(:,im)
              b(1)=uim(im)
c b(1) is the velocity on the immersed-boundary
              if(mrku(im)==1 .OR. mrku(im)==3) then
c cases for which the interpolation stencil is found along the normal direction
c or using the closest point criterion: the projection point has been found as
c intersection between the outward direction and a grid face
c N.B. the explicit velocity UB is used for this interpolation
                b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &            ,ub,nx,ny,nz,icyl,diru(im))
c b(2) is the interpolated velocity at the extension point
              elseif(mrku(im)==2) then
c case for which the fluid point is found along a grid line
                do ism=1,nsm
                  i1 = iu(im)+int(nxu(im))
                  j1 = ju(im)+int(nyu(im))
                  k1 = ku(im)+int(nzu(im))
                  call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
                  b(ism+1)=ub(i1,j1,k1)
                enddo
              endif

              call lubksb(a,nsp,nsp,indx,b)   !!!!!!
c the interpolated velocity at the interface point is stored in b(1)

              us(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))
     &           +  b(1)-ub(iu(im),ju(im),ku(im))
c the above row allows to eliminate one half of the terms treated by C.N.
c N.B. the velocity UB used to define the boundary condition is explicit

c the value of the boundary condition is stored in UB
              ub(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))

              if(ku(im).eq.kz1) then
                 ilu = ilu+1
                 indexl(ilu,1) = iu(im)
                 indexl(ilu,2) = ju(im)
              endif
              if(ku(im).eq.kz2) then
                 iru = iru+1
                 indexr(iru,1) = iu(im)
                 indexr(iru,2) = ju(im)
              endif

            ENDDO

          ENDDO

          call local_refreshbc(indexl(1:ilu,:),ilu,indexr(1:iru,:),iru,ub,nx,ny,nz)

c          CALL REFRESHBC(US,NX*NY,NZ)
!          CALL REFRESHBC(UB,NX*NY,NZ)

          !Periodic y bc
          us(:,1 ,:) = us(:,ny-1,:)
          us(:,ny,:) = us(:,  2 ,:)

          !Centerline bc
          if(icyl==1) then
            if(ivel==1) then
              do j=1,ny
                us(1,j,:) = 0.5*(us(2,j,:)-us(2,jsym(j),:))   !!!!!! - instead of +
              enddo
            elseif(ivel==2) then
              do j=1,ny
                us(1,j,:) =-us(2,jsym(j),:)
              enddo
            elseif(ivel==3) then     !!!!!! 3 instead of 2
              do j=1,ny
                us(1,j,:) = us(2,jsym(j),:)
              enddo
            endif
          endif

          indexlt(ilut+1:ilut+ilu,:) = indexl(1:ilu,:)
          indexrt(irut+1:irut+iru,:) = indexr(1:iru,:)
          ilut = ilut + ilu
          irut = irut + iru

        ENDDO

c the value of the boundary condition (with the implicit treatment) is stored in UB   !!!!!!
!!!!!!        DO ilp=1,nflumax
!!!!!!          DO im=limu(ilp)+1,limu(ilp)+mimu(ilp)
!!!!!!            ub(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))
!!!!!!          ENDDO
!!!!!!        ENDDO

        call local_refreshbc(indexlt(1:ilut,:),ilut,indexrt(1:irut,:),irut,us,nx,ny,nz)

!        CALL REFRESHBC(US,NX*NY,NZ)
!!!!!!        CALL REFRESHBC(UB,NX*NY,NZ)

      ELSE
c explicit treatment

        DO ilp=1,nloopg     !!!!!! instead of nflu

          ilu = 0
          iru = 0

          DO ibd = 1, nbd

            if(ilp.gt.nflumax(ibd)) cycle

            DO im=limu(ibd,ilp)+1,limu(ibd,ilp)+mimu(ibd,ilp)
              a(:,:)=umtrx(:,:,im)
              indx(:)=uindx(:,im)
              b(1)=uim(im)
              if(mrku(im)==1 .OR. mrku(im)==3) then
                b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &            ,us,nx,ny,nz,icyl,diru(im))
              elseif(mrku(im)==2) then
                do ism=1,nsm
                  i1 = iu(im)+int(nxu(im))
                  j1 = ju(im)+int(nyu(im))
                  k1 = ku(im)+int(nzu(im))
                  call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
                  b(ism+1)=us(i1,j1,k1)
                enddo
              endif

              call lubksb(a,nsp,nsp,indx,b)

              us(iu(im),ju(im),ku(im))=b(1)

              if(ku(im).eq.kz1) then
                 ilu = ilu+1
                 indexl(ilu,1) = iu(im)
                 indexl(ilu,2) = ju(im)
              endif
              if(ku(im).eq.kz2) then
                 iru = iru+1
                 indexr(iru,1) = iu(im)
                 indexr(iru,2) = ju(im)
              endif

            ENDDO

          ENDDO

          call local_refreshbc(indexl(1:ilu,:),ilu,indexr(1:iru,:),iru,us,nx,ny,nz)

!          CALL REFRESHBC(US,NX*NY,NZ)

          !Periodic y bc
          us(:,1 ,:) = us(:,ny-1,:)
          us(:,ny,:) = us(:,  2 ,:)

          !Centerline bc
          if(icyl==1) then
            if(ivel==1) then
              do j=1,ny
                us(1,j,:) = 0.5*(us(2,j,:)-us(2,jsym(j),:))
              enddo
            elseif(ivel==2) then
              do j=1,ny
                us(1,j,:) =-us(2,jsym(j),:)
              enddo
            elseif(ivel==3) then     !!!!!! 3 instead of 2
              do j=1,ny
                us(1,j,:) = us(2,jsym(j),:)
              enddo
            endif
          endif

        ENDDO

      ENDIF


      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          us(1,j,:) = us(2,jsym(j),:)
c        enddo
c      endif


      return

      end
C-------------------------------------------------------------------------

C---- subroutine local_refreshbc --------------N. Beratlis-30 Nov. 2010---
C------------------------------------------------A. Posa - 2 Nov 2012-----
C
C     Refresh BC at the ghost layers
C
C-------------------------------------------------------------------------
      subroutine local_refreshbc(indexl,nl_sup,indexr,nr_sup,var,nx,ny,nz)
c
      include 'common.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nl_sup,nr_sup
      integer indexl(nl_sup,2),indexr(nr_sup,2)
      real    var(nx,ny,nz)
c
c.... Local arrays
      integer ii,nl_req,nr_req
      real usupl(nl_sup),usupr(nr_sup)
      integer, dimension(:,:), allocatable :: index_req
      real, dimension(:), allocatable :: u_req
      integer status(mpi_status_size)
      real buffer(nx,ny)
c

      CALL MPI_SENDRECV(NL_SUP,1,MPI_INTEGER,MYLEFT,1
     &                   ,NR_REQ,1,MPI_INTEGER,MYRIGHT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

      CALL MPI_SENDRECV(NR_SUP,1,MPI_INTEGER,MYRIGHT,1
     &                   ,NL_REQ,1,MPI_INTEGER,MYLEFT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

      IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) buffer(1:nx,1:ny) = var(1:nx,1:ny,nz)

      IF(MOD(MYRANK,2)==0) THEN

        IF(NL_SUP>0) THEN
          CALL MPI_SEND(INDEXL,2*NL_SUP,MPI_INTEGER,MYLEFT,2,MPI_COMM_EDDY,IERR)
          DO II=1,NL_SUP
            usupl(ii) = var(indexl(ii,1),indexl(ii,2),kz1)
          ENDDO
          CALL MPI_SEND(USUPL,NL_SUP,MTYPE,MYLEFT,6,MPI_COMM_EDDY,IERR)
        ENDIF

        IF(NR_REQ>0) THEN
          ALLOCATE(INDEX_REQ(NR_REQ,2),U_REQ(NR_REQ))
          CALL MPI_RECV(INDEX_REQ,NR_REQ*2,MPI_INTEGER,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(U_REQ,NR_REQ,MTYPE,MYRIGHT,6,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NR_REQ
            var(index_req(ii,1),index_req(ii,2),nz) = u_req(ii)
          ENDDO
          DEALLOCATE(INDEX_REQ,U_REQ)
        ENDIF

      ELSE

        IF(NR_REQ>0) THEN
          ALLOCATE(INDEX_REQ(NR_REQ,2),U_REQ(NR_REQ))
          CALL MPI_RECV(INDEX_REQ,NR_REQ*2,MPI_INTEGER,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(U_REQ,NR_REQ,MTYPE,MYRIGHT,6,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NR_REQ
            var(index_req(ii,1),index_req(ii,2),nz) = u_req(ii)
          ENDDO
          DEALLOCATE(INDEX_REQ,U_REQ)
        ENDIF

        IF(NL_SUP>0) THEN
          CALL MPI_SEND(INDEXL,2*NL_SUP,MPI_INTEGER,MYLEFT,2,MPI_COMM_EDDY,IERR)
          DO II=1,NL_SUP
            usupl(ii) = var(indexl(ii,1),indexl(ii,2),kz1)
          ENDDO
          CALL MPI_SEND(USUPL,NL_SUP,MTYPE,MYLEFT,6,MPI_COMM_EDDY,IERR)
        ENDIF

      ENDIF

      IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) var(1:nx,1:ny,nz) = buffer(1:nx,1:ny)

      IF(MYRANK==0 .AND. ITYPE(5)/=500) buffer(1:nx,1:ny) = var(1:nx,1:ny,1)

      IF(MOD(MYRANK,2)==0) THEN

        IF(NR_SUP>0) THEN
          CALL MPI_SEND(INDEXR,2*NR_SUP,MPI_INTEGER,MYRIGHT,3,MPI_COMM_EDDY,IERR)
          DO II=1,NR_SUP
            usupr(ii) = var(indexr(ii,1),indexr(ii,2),kz2)
          ENDDO
          CALL MPI_SEND(USUPR,NR_SUP,MTYPE,MYRIGHT,5,MPI_COMM_EDDY,IERR)
        ENDIF

        IF(NL_REQ>0) THEN
          ALLOCATE(INDEX_REQ(NL_REQ,2),U_REQ(NL_REQ))
          CALL MPI_RECV(INDEX_REQ,NL_REQ*2,MPI_INTEGER,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(U_REQ,NL_REQ,MTYPE,MYLEFT,5,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NL_REQ
            var(index_req(ii,1),index_req(ii,2),1) = u_req(ii)
          ENDDO
          DEALLOCATE(INDEX_REQ,U_REQ)
        ENDIF

      ELSE

        IF(NL_REQ>0) THEN
          ALLOCATE(INDEX_REQ(NL_REQ,2),U_REQ(NL_REQ))
          CALL MPI_RECV(INDEX_REQ,NL_REQ*2,MPI_INTEGER,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
          CALL MPI_RECV(U_REQ,NL_REQ,MTYPE,MYLEFT,5,MPI_COMM_EDDY,STATUS,IERR)
          DO II=1,NL_REQ
            var(index_req(ii,1),index_req(ii,2),1) = u_req(ii)
          ENDDO
          DEALLOCATE(INDEX_REQ,U_REQ)
        ENDIF

        IF(NR_SUP>0) THEN
          CALL MPI_SEND(INDEXR,2*NR_SUP,MPI_INTEGER,MYRIGHT,3,MPI_COMM_EDDY,IERR)
          DO II=1,NR_SUP
            usupr(ii) = var(indexr(ii,1),indexr(ii,2),kz2)
          ENDDO
          CALL MPI_SEND(USUPR,NR_SUP,MTYPE,MYRIGHT,5,MPI_COMM_EDDY,IERR)
        ENDIF

      ENDIF

      IF(MYRANK==0 .AND. ITYPE(5)/=500) var(1:nx,1:ny,1) = buffer(1:nx,1:ny)

      return

      end
C-------------------------------------------------------------------------

c---- subroutine correctpres_mod2 ----------N. Beratlis-18 May 2009---
c---------------------------------------------A. Posa - 3 Nov 2012----
C
C     PURPOSE: Correct pressure
C
C---------------------------------------------------------------------
      subroutine correctpres_mod2(p,nx,ny,nz,ip,jp,kp,pmtrx
     &     ,pim,nxp,nyp,nzp,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp
     &     ,nfacet,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax,tlevel,dt,mbd,nbd,nmax)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nfacet,mbd,nflpmax(nbd),nbd
      integer limp(nbd,nflp),mimp(nbd,nflp)
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      integer mrkp(nfcmax),dirp(nfcmax)
      integer pindx(nsp,nfcmax)
      integer fp(nfcmax)
      real    tlevel,dt
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    pim(nfcmax),xnp(nfcmax),ynp(nfcmax),znp(nfcmax),dpdnn(nfcmax)
      real    unvect(3,nfacet)
      real    nxp(nfcmax),nyp(nfcmax),nzp(nfcmax)
      real    pmtrx(nsp,nsp,nfcmax)
      integer nmax
c
c.... Local arrays
      integer i,j,ism,ii,jj,kk,in,jn,kn,in1,jn1,kn1,in2,jn2,kn2,ilp
      integer i1,j1,k1
      integer i2,j2,k2
      integer i3,j3,k3
      integer indx(nsp)
      real    xpcar,ypcar,zpcar,s1,s2,xp1car,yp1car,zp1car,xp2car,yp2car,zp2car
      real    v1(3),v2(3)
      real    a(nsp,nsp),b(nsp)
      integer nloopmax,ibd,ilb,ile,nloopg
      integer ilpr,irpr,indexl(nmax,2),indexr(nmax,2)
c
c.... Function
      real    interp_cellface
      real    dotproduct,anglerad

c The normal pressure gradient is equal to 0 in the case of stationary
c boundaries
      dpdnn = 0.0

      nloopmax = maxval(nflpmax(1:nbd))
      CALL MPI_ALLREDUCE(nloopmax,nloopg,1,mpi_integer,MPI_MAX,MPI_COMM_EDDY,IERR)

      do ilp=1,nloopg

        ilpr = 0
        irpr = 0

        do ibd=1,nbd

          if(ilp.gt.nflpmax(ibd)) cycle

          ilb = lb(ibd)+1
          ile = lb(ibd)+mb(ibd)

          if(ibd>=mbd) then
c For the moving boundaries the normal pressure gradient is evaluated
            call presgrad_flagp(dpdnn,mrkp,xnp,ynp,znp,nxp,nyp,nzp,ip,jp,kp
     &          ,xc_car,yc_car,zc,nx,ny,nz,unvect(:,ilb:ile),fp
     &          ,mb(ibd),limp(ibd,ilp),mimp(ibd,ilp),tlevel,dt,ibd)
          endif

          do i=limp(ibd,ilp)+1,limp(ibd,ilp)+mimp(ibd,ilp)

            a(:,:) = pmtrx(:,:,i)
            indx(:) = pindx(:,i)

!             in = ip(i)          !!!!!!!
!             jn = jp(i)          !!!!!!!
!             kn = kp(i)          !!!!!!!

            b(1) = dpdnn(i)

            if(mrkp(i)==1 .OR. mrkp(i)==3) then
c the pressure at the extension fluid point is evaluated by interpolation
              pim(i) = interp_cellface(nxp(i),nyp(i),nzp(i)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dirp(i))
            else
              i1 = ip(i)+int(nxp(i))
              j1 = jp(i)+int(nyp(i))
              k1 = kp(i)+int(nzp(i))
              call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
c the pressure at the exterior grid point is determined
              pim(i) = p(i1,j1,k1)
            endif
            b(2) = pim(i)

c the local system for the pressure stencil is solved
            call lubksb(a,nsp,nsp,indx,b)
c the pressure condition at the interface point is enforced
            p(ip(i),jp(i),kp(i))=b(1)

            if(kp(i).eq.kz1) then
              ilpr = ilpr + 1
              indexl(ilpr,1) = ip(i)
              indexl(ilpr,2) = jp(i)
            endif
            if(kp(i).eq.kz2) then
              irpr = irpr + 1
              indexr(irpr,1) = ip(i)
              indexr(irpr,2) = jp(i)
            endif

          enddo

        enddo

        call local_refreshbc(indexl(1:ilpr,:),ilpr,indexr(1:irpr,:),irpr,p,nx,ny,nz)

!        call refreshbc(p,nx*ny,nz)

        !Periodic y bc
        p(:,1,:) = p(:,ny-1,:)
        p(:,ny,:) = p(:,2,:)

        !Centerline bc
        if(icyl==1) then
          do j=1,ny
            p(1,j,:) = p(2,jsym(j),:)
          enddo
        endif

      enddo

      !Periodic y bc
c      p(:,1,:) = p(:,ny-1,:)
c      p(:,ny,:) = p(:,2,:)

      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          p(1,j,:) = p(2,jsym(j),:)
c        enddo
c      endif

      return

      end


      subroutine correctdens_mod2(dens,nx,ny,nz,ip,jp,kp,pmtrx
     &     ,rhoim,nxp,nyp,nzp,xnp,ynp,znp,dirp,mrkp,pindx,drhodnn,unvect,fp
     &     ,nfacet,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax,tlevel,dt,mbd,nbd,nmax)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nfacet,mbd,nflpmax(nbd),nbd
      integer limp(nbd,nflp),mimp(nbd,nflp)
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      integer mrkp(nfcmax),dirp(nfcmax)
      integer pindx(nsp,nfcmax)
      integer fp(nfcmax)
      real    tlevel,dt
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    dens(nx,ny,nz)
      real    rhoim(nfcmax),xnp(nfcmax),ynp(nfcmax),znp(nfcmax),drhodnn(nfcmax)
      real    unvect(3,nfacet)
      real    nxp(nfcmax),nyp(nfcmax),nzp(nfcmax)
      real    pmtrx(nsp,nsp,nfcmax)
      integer nmax
c
c.... Local arrays
      integer i,j,ism,ii,jj,kk,in,jn,kn,in1,jn1,kn1,in2,jn2,kn2,ilp
      integer i1,j1,k1
      integer i2,j2,k2
      integer i3,j3,k3
      integer indx(nsp)
      real    xpcar,ypcar,zpcar,s1,s2,xp1car,yp1car,zp1car,xp2car,yp2car,zp2car
      real    v1(3),v2(3)
      real    a(nsp,nsp),b(nsp)
      integer nloopmax,ibd,ilb,ile,nloopg
      integer ilpr,irpr,indexl(nmax,2),indexr(nmax,2)
c
c.... Function
      real    interp_cellface
      real    dotproduct,anglerad

c The normal pressure gradient is equal to 0 in the case of stationary
c boundaries
      drhodnn = 0.0

      nloopmax = maxval(nflpmax(1:nbd))
      CALL MPI_ALLREDUCE(nloopmax,nloopg,1,mpi_integer,MPI_MAX,MPI_COMM_EDDY,IERR)

      do ilp=1,nloopg

        ilpr = 0
        irpr = 0

        do ibd=1,nbd

          if(ilp.gt.nflpmax(ibd)) cycle

          ilb = lb(ibd)+1
          ile = lb(ibd)+mb(ibd)

!           if(ibd>=mbd) then
! c For the moving boundaries the normal pressure gradient is evaluated
!             call presgrad_flagp(dpdnn,mrkp,xnp,ynp,znp,nxp,nyp,nzp,ip,jp,kp
!      &          ,xc_car,yc_car,zc,nx,ny,nz,unvect(:,ilb:ile),fp
!      &          ,mb(ibd),limp(ibd,ilp),mimp(ibd,ilp),tlevel,dt,ibd)
!           endif

          do i=limp(ibd,ilp)+1,limp(ibd,ilp)+mimp(ibd,ilp)

            a(:,:) = pmtrx(:,:,i)
            indx(:) = pindx(:,i)

!             in = ip(i)          !!!!!!!
!             jn = jp(i)          !!!!!!!
!             kn = kp(i)          !!!!!!!

            b(1) = drhodnn(i)

            if(mrkp(i)==1 .OR. mrkp(i)==3) then
            
c the pressure at the extension fluid point is evaluated by interpolation
              rhoim(i) = interp_cellface(nxp(i),nyp(i),nzp(i)
     &            ,xc,yc,zc,dens,nx,ny,nz,icyl,dirp(i))
            else
              i1 = ip(i)+int(nxp(i))
              j1 = jp(i)+int(nyp(i))
              k1 = kp(i)+int(nzp(i))
              call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
c the pressure at the exterior grid point is determined
              rhoim(i) = dens(i1,j1,k1)
            endif
            b(2) = rhoim(i)

c the local system for the pressure stencil is solved
            call lubksb(a,nsp,nsp,indx,b)
c the pressure condition at the interface point is enforced
            dens(ip(i),jp(i),kp(i))=b(1)

            if(kp(i).eq.kz1) then
              ilpr = ilpr + 1
              indexl(ilpr,1) = ip(i)
              indexl(ilpr,2) = jp(i)
            endif
            if(kp(i).eq.kz2) then
              irpr = irpr + 1
              indexr(irpr,1) = ip(i)
              indexr(irpr,2) = jp(i)
            endif

          enddo

        enddo

        call local_refreshbc(indexl(1:ilpr,:),ilpr,indexr(1:irpr,:),irpr,dens,nx,ny,nz)

!        call refreshbc(p,nx*ny,nz)

!         !Periodic y bc
!         p(:,1,:) = p(:,ny-1,:)
!         p(:,ny,:) = p(:,2,:)
! 
!         !Centerline bc
!         if(icyl==1) then
!           do j=1,ny
!             p(1,j,:) = p(2,jsym(j),:)
!           enddo
!         endif

      enddo

      !Periodic y bc
c      p(:,1,:) = p(:,ny-1,:)
c      p(:,ny,:) = p(:,2,:)

      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          p(1,j,:) = p(2,jsym(j),:)
c        enddo
c      endif

      return

      end
c
c---------------------------------------------------------------------

C------------------------------------------------A. Posa - 5 Nov 2012-----
C
C     PURPOSE: Refresh ghost cells of flag array between processors.
C
C-------------------------------------------------------------------------
      subroutine local_refreshflag(indexl,nl_sup,indexr,nr_sup,flag,nx,ny,nz)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nl_sup,nr_sup
      integer indexl(nl_sup,2),indexr(nr_sup,2)
      integer flag(nx,ny,nz)
c
c.... Local arrays
      integer ii,nl_req,nr_req
      integer fsupl(nl_sup),fsupr(nr_sup)
      integer, dimension(:,:), allocatable :: index_req
      integer, dimension(:), allocatable :: flag_req
      integer status(mpi_status_size)
      integer buffer(nx,ny)
c
      if(idomy==0) then

        CALL MPI_SENDRECV(NL_SUP,1,MPI_INTEGER,MYLEFT,1
     &                   ,NR_REQ,1,MPI_INTEGER,MYRIGHT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

        CALL MPI_SENDRECV(NR_SUP,1,MPI_INTEGER,MYRIGHT,1
     &                   ,NL_REQ,1,MPI_INTEGER,MYLEFT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

        IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) buffer(1:nx,1:ny) = flag(1:nx,1:ny,nz)

        IF(MOD(MYRANK,2)==0) THEN

          IF(NL_SUP>0) THEN
            CALL MPI_SEND(INDEXL,2*NL_SUP,MPI_INTEGER,MYLEFT,2,MPI_COMM_EDDY,IERR)
            DO II=1,NL_SUP
              fsupl(ii) = flag(indexl(ii,1),indexl(ii,2),kz1)
            ENDDO
            CALL MPI_SEND(FSUPL,NL_SUP,MPI_INTEGER,MYLEFT,6,MPI_COMM_EDDY,IERR)
          ENDIF

          IF(NR_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NR_REQ,2),FLAG_REQ(NR_REQ))
            CALL MPI_RECV(INDEX_REQ,NR_REQ*2,MPI_INTEGER,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(FLAG_REQ,NR_REQ,MPI_INTEGER,MYRIGHT,6,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              flag(index_req(ii,1),index_req(ii,2),nz) = flag_req(ii)
            ENDDO
            DEALLOCATE(INDEX_REQ,FLAG_REQ)
          ENDIF

        ELSE

          IF(NR_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NR_REQ,2),FLAG_REQ(NR_REQ))
            CALL MPI_RECV(INDEX_REQ,NR_REQ*2,MPI_INTEGER,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(FLAG_REQ,NR_REQ,MPI_INTEGER,MYRIGHT,6,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              flag(index_req(ii,1),index_req(ii,2),nz) = flag_req(ii)
            ENDDO
            DEALLOCATE(INDEX_REQ,FLAG_REQ)
          ENDIF

          IF(NL_SUP>0) THEN
            CALL MPI_SEND(INDEXL,2*NL_SUP,MPI_INTEGER,MYLEFT,2,MPI_COMM_EDDY,IERR)
            DO II=1,NL_SUP
              fsupl(ii) = flag(indexl(ii,1),indexl(ii,2),kz1)
            ENDDO
            CALL MPI_SEND(FSUPL,NL_SUP,MPI_INTEGER,MYLEFT,6,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

        IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) flag(1:nx,1:ny,nz) = buffer(1:nx,1:ny)

        IF(MYRANK==0 .AND. ITYPE(5)/=500) buffer(1:nx,1:ny) = flag(1:nx,1:ny,1)

        IF(MOD(MYRANK,2)==0) THEN

          IF(NR_SUP>0) THEN
            CALL MPI_SEND(INDEXR,2*NR_SUP,MPI_INTEGER,MYRIGHT,3,MPI_COMM_EDDY,IERR)
            DO II=1,NR_SUP
              fsupr(ii) = flag(indexr(ii,1),indexr(ii,2),kz2)
            ENDDO
            CALL MPI_SEND(FSUPR,NR_SUP,MPI_INTEGER,MYRIGHT,5,MPI_COMM_EDDY,IERR)
          ENDIF

          IF(NL_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NL_REQ,2),FLAG_REQ(NL_REQ))
            CALL MPI_RECV(INDEX_REQ,NL_REQ*2,MPI_INTEGER,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(FLAG_REQ,NL_REQ,MPI_INTEGER,MYLEFT,5,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              flag(index_req(ii,1),index_req(ii,2),1) = flag_req(ii)
            ENDDO
            DEALLOCATE(INDEX_REQ,FLAG_REQ)
          ENDIF

        ELSE

          IF(NL_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NL_REQ,2),FLAG_REQ(NL_REQ))
            CALL MPI_RECV(INDEX_REQ,NL_REQ*2,MPI_INTEGER,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
            CALL MPI_RECV(FLAG_REQ,NL_REQ,MPI_INTEGER,MYLEFT,5,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              flag(index_req(ii,1),index_req(ii,2),1) = flag_req(ii)
            ENDDO
            DEALLOCATE(INDEX_REQ,FLAG_REQ)
          ENDIF

          IF(NR_SUP>0) THEN
            CALL MPI_SEND(INDEXR,2*NR_SUP,MPI_INTEGER,MYRIGHT,3,MPI_COMM_EDDY,IERR)
            DO II=1,NR_SUP
              fsupr(ii) = flag(indexr(ii,1),indexr(ii,2),kz2)
            ENDDO
            CALL MPI_SEND(FSUPR,NR_SUP,MPI_INTEGER,MYRIGHT,5,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

        IF(MYRANK==0 .AND. ITYPE(5)/=500) flag(1:nx,1:ny,1) = buffer(1:nx,1:ny)

      else

        call mpi_sendrecv(flag(:,ny-1,:),nx*nz,mpi_integer,myright,10
     &        ,           flag(:,1,:)   ,nx*nz,mpi_integer,myleft ,10
     &        ,           mpi_comm_eddy,status,ierr)

        call mpi_sendrecv(flag(:,2,:   ),nx*nz,mpi_integer,myleft,11
     &        ,           flag(:,ny,:  ),nx*nz,mpi_integer,myright ,11
     &        ,           mpi_comm_eddy,status,ierr)

      endif

      return

      end
C-------------------------------------------------------------------------

C---- subroutine local_refreshflag_smp --------N. Beratlis-30 Nov. 2010---
C------------------------------------------------A. Posa - 5 Nov 2012-----
C
C     PURPOSE: Refresh ghost cells of flag array between processors.
C
C-------------------------------------------------------------------------
      subroutine local_refreshflag_smp(indexl,nl_sup,indexr,nr_sup,flag,nx,ny,nz)
c
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output arrays
      integer nx,ny,nz,nl_sup,nr_sup
      integer indexl(nl_sup,2),indexr(nr_sup,2)
      integer flag(nx,ny,nz)
c
c.... Local arrays
      integer ii,nl_req,nr_req
      integer, dimension(:,:), allocatable :: index_req
      integer status(mpi_status_size)
      integer buffer(nx,ny)
c
      if(idomy==0) then

        CALL MPI_SENDRECV(NL_SUP,1,MPI_INTEGER,MYLEFT,1
     &                   ,NR_REQ,1,MPI_INTEGER,MYRIGHT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

        CALL MPI_SENDRECV(NR_SUP,1,MPI_INTEGER,MYRIGHT,1
     &                   ,NL_REQ,1,MPI_INTEGER,MYLEFT,1
     &                   ,MPI_COMM_EDDY,STATUS,IERR)

        IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) buffer(1:nx,1:ny) = flag(1:nx,1:ny,nz)

        IF(MOD(MYRANK,2)==0) THEN

          IF(NL_SUP>0) THEN
            CALL MPI_SEND(INDEXL,2*NL_SUP,MPI_INTEGER,MYLEFT,2,MPI_COMM_EDDY,IERR)
          ENDIF

          IF(NR_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NR_REQ,2))
            CALL MPI_RECV(INDEX_REQ,NR_REQ*2,MPI_INTEGER,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              flag(index_req(ii,1),index_req(ii,2),nz) = 0
            ENDDO
            DEALLOCATE(INDEX_REQ)
          ENDIF

        ELSE

          IF(NR_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NR_REQ,2))
            CALL MPI_RECV(INDEX_REQ,NR_REQ*2,MPI_INTEGER,MYRIGHT,2,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NR_REQ
              flag(index_req(ii,1),index_req(ii,2),nz) = 0
            ENDDO
            DEALLOCATE(INDEX_REQ)
          ENDIF

          IF(NL_SUP>0) THEN
            CALL MPI_SEND(INDEXL,2*NL_SUP,MPI_INTEGER,MYLEFT,2,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

        IF(MYRANK==MYSIZE-1 .AND. ITYPE(6)/=500) flag(1:nx,1:ny,nz) = buffer(1:nx,1:ny)

        IF(MYRANK==0 .AND. ITYPE(5)/=500) buffer(1:nx,1:ny) = flag(1:nx,1:ny,1)

        IF(MOD(MYRANK,2)==0) THEN

          IF(NR_SUP>0) THEN
            CALL MPI_SEND(INDEXR,2*NR_SUP,MPI_INTEGER,MYRIGHT,3,MPI_COMM_EDDY,IERR)
          ENDIF

          IF(NL_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NL_REQ,2))
            CALL MPI_RECV(INDEX_REQ,NL_REQ*2,MPI_INTEGER,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              flag(index_req(ii,1),index_req(ii,2),1) = 0
            ENDDO
            DEALLOCATE(INDEX_REQ)
          ENDIF

        ELSE

          IF(NL_REQ>0) THEN
            ALLOCATE(INDEX_REQ(NL_REQ,2))
            CALL MPI_RECV(INDEX_REQ,NL_REQ*2,MPI_INTEGER,MYLEFT,3,MPI_COMM_EDDY,STATUS,IERR)
            DO II=1,NL_REQ
              flag(index_req(ii,1),index_req(ii,2),1) = 0
            ENDDO
            DEALLOCATE(INDEX_REQ)
          ENDIF

          IF(NR_SUP>0) THEN
            CALL MPI_SEND(INDEXR,2*NR_SUP,MPI_INTEGER,MYRIGHT,3,MPI_COMM_EDDY,IERR)
          ENDIF

        ENDIF

        IF(MYRANK==0 .AND. ITYPE(5)/=500) flag(1:nx,1:ny,1) = buffer(1:nx,1:ny)

      else

        call mpi_sendrecv(flag(:,ny-1,:),nx*nz,mpi_integer,myright,10
     &        ,           flag(:,1,:)   ,nx*nz,mpi_integer,myleft ,10
     &        ,           mpi_comm_eddy,status,ierr)

        call mpi_sendrecv(flag(:,2,:   ),nx*nz,mpi_integer,myleft,11
     &        ,           flag(:,ny,:  ),nx*nz,mpi_integer,myright ,11
     &        ,           mpi_comm_eddy,status,ierr)

      endif

      return

      end
C-------------------------------------------------------------------------

C---- subroutine momforc_mod3 ------------------N. Beratlis-25 Apr. 2009--
C--------------------------------------------------A. Posa - 18 Nov 2012--
C
C     PURPOSE: Compute momentum forcing and correct velocity near immersed
C     body surface.
C
C-------------------------------------------------------------------------
      subroutine momforc_mod3(us,ub,xu,yu,zu,nx,ny,nz,iu,ju,ku,umtrx,uindx
     &     ,uim,mrku,diru,nxu,nyu,nzu,nflumax,limu,mimu,impl,ivel,nbd,nmax)

c      IMPLICIT NONE
      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'

      integer nx,ny,nz,impl,nflumax(nbd),ivel,nbd
      integer limu(nbd,nflu),mimu(nbd,nflu)
      integer iu(nfcmax),ju(nfcmax),ku(nfcmax),mrku(nfcmax),diru(nfcmax)
      integer uindx(nsp,nfcmax)
      real    xu(nx),yu(ny),zu(nz)
      real    umtrx(nsp,nsp,nfcmax)
      real    nxu(nfcmax),nyu(nfcmax),nzu(nfcmax),uim(nfcmax)
      real    us(nx,ny,nz),ub(nx,ny,nz)
      integer nmax
c
c.... Local var. and arrays
      integer im,ism,i,j,k,ilp,i1,j1,k1,nloopmax,ibd,nloopg
      real    xint,yint,zint,xp,yp,zp
      integer indx(nsp)
      real    a(nsp,nsp),b(nsp)
      integer ilu,iru,indexl(nmax,2),indexr(nmax,2)
      integer ilut,irut,indexlt(nmax,2),indexrt(nmax,2)
      integer ibdf,ibdl
      parameter (ibdf=1,ibdl=1)
c
C.... Functions
      real    interp_cellface,anglerad

      nloopmax = maxval(nflumax(ibdf:ibdl))
      CALL MPI_ALLREDUCE(nloopmax,nloopg,1,mpi_integer,MPI_MAX,MPI_COMM_EDDY,IERR)

      IF(impl==1) THEN
c implicit treatment

        ilut = 0
        irut = 0

        DO ilp=1,nloopg

          ilu = 0
          iru = 0

          DO ibd = ibdf, ibdl

            if(ilp.gt.nflumax(ibd)) cycle

            DO im=limu(ibd,ilp)+1,limu(ibd,ilp)+mimu(ibd,ilp)
              a(:,:)=umtrx(:,:,im)
              indx(:)=uindx(:,im)
              b(1)=uim(im)
c b(1) is the velocity on the immersed-boundary
              if(mrku(im)==1 .OR. mrku(im)==3) then
c cases for which the interpolation stencil is found along the normal direction
c or using the closest point criterion: the projection point has been found as
c intersection between the outward direction and a grid face
c N.B. the explicit velocity UB is used for this interpolation
                b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &            ,ub,nx,ny,nz,icyl,diru(im))
c b(2) is the interpolated velocity at the extension point
              elseif(mrku(im)==2) then
c case for which the fluid point is found along a grid line
                do ism=1,nsm
                  i1 = iu(im)+int(nxu(im))
                  j1 = ju(im)+int(nyu(im))
                  k1 = ku(im)+int(nzu(im))
                  call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
                  b(ism+1)=ub(i1,j1,k1)
                enddo
              endif

              call lubksb(a,nsp,nsp,indx,b)   !!!!!!
c the interpolated velocity at the interface point is stored in b(1)

              us(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))
     &           +  b(1)-ub(iu(im),ju(im),ku(im))
c the above row allows to eliminate one half of the terms treated by C.N.
c N.B. the velocity UB used to define the boundary condition is explicit

c the value of the boundary condition is stored in UB
              ub(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))

              if(ku(im).eq.kz1) then
                 ilu = ilu+1
                 indexl(ilu,1) = iu(im)
                 indexl(ilu,2) = ju(im)
              endif
              if(ku(im).eq.kz2) then
                 iru = iru+1
                 indexr(iru,1) = iu(im)
                 indexr(iru,2) = ju(im)
              endif

            ENDDO

          ENDDO

          call local_refreshbc(indexl(1:ilu,:),ilu,indexr(1:iru,:),iru,ub,nx,ny,nz)

c          CALL REFRESHBC(US,NX*NY,NZ)
!          CALL REFRESHBC(UB,NX*NY,NZ)

          !Periodic y bc
          us(:,1 ,:) = us(:,ny-1,:)
          us(:,ny,:) = us(:,  2 ,:)

          !Centerline bc
          if(icyl==1) then
            if(ivel==1) then
              do j=1,ny
                us(1,j,:) = 0.5*(us(2,j,:)-us(2,jsym(j),:))    !!!!!! - instead of +
              enddo
            elseif(ivel==2) then
              do j=1,ny
                us(1,j,:) =-us(2,jsym(j),:)
              enddo
            elseif(ivel==3) then     !!!!!! 3 instead of 2
              do j=1,ny
                us(1,j,:) = us(2,jsym(j),:)
              enddo
            endif
          endif

          indexlt(ilut+1:ilut+ilu,:) = indexl(1:ilu,:)
          indexrt(irut+1:irut+iru,:) = indexr(1:iru,:)
          ilut = ilut + ilu
          irut = irut + iru

        ENDDO

c the value of the boundary condition (with the implicit treatment) is stored in UB   !!!!!!
!!!!!!        DO ilp=1,nflumax
!!!!!!          DO im=limu(ilp)+1,limu(ilp)+mimu(ilp)
!!!!!!            ub(iu(im),ju(im),ku(im)) = us(iu(im),ju(im),ku(im))
!!!!!!          ENDDO
!!!!!!        ENDDO

        call local_refreshbc(indexlt(1:ilut,:),ilut,indexrt(1:irut,:),irut,us,nx,ny,nz)

!        CALL REFRESHBC(US,NX*NY,NZ)
!!!!!!        CALL REFRESHBC(UB,NX*NY,NZ)

      ELSE
c explicit treatment

        DO ilp=1,nloopg     !!!!!! instead of nflu

          ilu = 0
          iru = 0

          DO ibd = ibdf, ibdl

            if(ilp.gt.nflumax(ibd)) cycle

            DO im=limu(ibd,ilp)+1,limu(ibd,ilp)+mimu(ibd,ilp)
              a(:,:)=umtrx(:,:,im)
              indx(:)=uindx(:,im)
              b(1)=uim(im)
              if(mrku(im)==1 .OR. mrku(im)==3) then
                b(2) = interp_cellface(nxu(im),nyu(im),nzu(im),xu,yu,zu
     &            ,us,nx,ny,nz,icyl,diru(im))
              elseif(mrku(im)==2) then
                do ism=1,nsm
                  i1 = iu(im)+int(nxu(im))
                  j1 = ju(im)+int(nyu(im))
                  k1 = ku(im)+int(nzu(im))
                  call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
                  b(ism+1)=us(i1,j1,k1)
                enddo
              endif

              call lubksb(a,nsp,nsp,indx,b)

              us(iu(im),ju(im),ku(im))=b(1)

              if(ku(im).eq.kz1) then
                 ilu = ilu+1
                 indexl(ilu,1) = iu(im)
                 indexl(ilu,2) = ju(im)
              endif
              if(ku(im).eq.kz2) then
                 iru = iru+1
                 indexr(iru,1) = iu(im)
                 indexr(iru,2) = ju(im)
              endif

            ENDDO

          ENDDO

          call local_refreshbc(indexl(1:ilu,:),ilu,indexr(1:iru,:),iru,us,nx,ny,nz)

!          CALL REFRESHBC(US,NX*NY,NZ)

          !Periodic y bc
          us(:,1 ,:) = us(:,ny-1,:)
          us(:,ny,:) = us(:,  2 ,:)

          !Centerline bc
          if(icyl==1) then
            if(ivel==1) then
              do j=1,ny
                us(1,j,:) = 0.5*(us(2,j,:)-us(2,jsym(j),:))
              enddo
            elseif(ivel==2) then
              do j=1,ny
                us(1,j,:) =-us(2,jsym(j),:)
              enddo
            elseif(ivel==3) then     !!!!!! 3 instead of 2
              do j=1,ny
                us(1,j,:) = us(2,jsym(j),:)
              enddo
            endif
          endif

        ENDDO

      ENDIF


      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          us(1,j,:) = us(2,jsym(j),:)
c        enddo
c      endif


      return

      end
C-------------------------------------------------------------------------

c---- subroutine correctpres_mod3 ----------N. Beratlis-18 May 2009---
c---------------------------------------------A. Posa - 18 Nov 2012---
C
C     PURPOSE: Correct pressure
C
C---------------------------------------------------------------------
      subroutine correctpres_mod3(p,nx,ny,nz,ip,jp,kp,pmtrx
     &     ,pim,nxp,nyp,nzp,xnp,ynp,znp,dirp,mrkp,pindx,dpdnn,unvect,fp
     &     ,nfacet,xc_car,yc_car,xc,yc,zc,limp,mimp,nflpmax,tlevel,dt,mbd,nbd,nmax)

      include 'common.h'
      include 'immersed.h'
      include 'mpif.h'
c
c.... Input/Output Arrays
      integer nx,ny,nz,nfacet,mbd,nflpmax(nbd),nbd
      integer limp(nbd,nflp),mimp(nbd,nflp)
      integer ip(nfcmax),jp(nfcmax),kp(nfcmax)
      integer mrkp(nfcmax),dirp(nfcmax)
      integer pindx(nsp,nfcmax)
      integer fp(nfcmax)
      real    tlevel,dt
      real    xc(nx),yc(ny),zc(nz)
      real    xc_car(nx,ny),yc_car(nx,ny)
      real    p(nx,ny,nz)
      real    pim(nfcmax),xnp(nfcmax),ynp(nfcmax),znp(nfcmax),dpdnn(nfcmax)
      real    unvect(3,nfacet)
      real    nxp(nfcmax),nyp(nfcmax),nzp(nfcmax)
      real    pmtrx(nsp,nsp,nfcmax)
      integer nmax
c
c.... Local arrays
      integer i,j,ism,ii,jj,kk,in,jn,kn,in1,jn1,kn1,in2,jn2,kn2,ilp
      integer i1,j1,k1
      integer i2,j2,k2
      integer i3,j3,k3
      integer indx(nsp)
      real    xpcar,ypcar,zpcar,s1,s2,xp1car,yp1car,zp1car,xp2car,yp2car,zp2car
      real    v1(3),v2(3)
      real    a(nsp,nsp),b(nsp)
      integer nloopmax,ibd,ilb,ile,nloopg
      integer ilpr,irpr,indexl(nmax,2),indexr(nmax,2)
      integer ibdf,ibdl
      parameter (ibdf=1,ibdl=1)
c
c.... Function
      real    interp_cellface
      real    dotproduct,anglerad

c The normal pressure gradient is equal to 0 in the case of stationary
c boundaries
      dpdnn = 0.0   

      nloopmax = maxval(nflpmax(ibdf:ibdl))
      CALL MPI_ALLREDUCE(nloopmax,nloopg,1,mpi_integer,MPI_MAX,MPI_COMM_EDDY,IERR)

      do ilp=1,nloopg

        ilpr = 0
        irpr = 0

        do ibd=ibdf,ibdl

          if(ilp.gt.nflpmax(ibd)) cycle

          ilb = lb(ibd)+1
          ile = lb(ibd)+mb(ibd)

          if(ibd>=mbd) then
c For the moving boundaries the normal pressure gradient is evaluated
            call presgrad_flagp(dpdnn,mrkp,xnp,ynp,znp,nxp,nyp,nzp,ip,jp,kp
     &          ,xc_car,yc_car,zc,nx,ny,nz,unvect(:,ilb:ile),fp
     &          ,mb(ibd),limp(ibd,ilp),mimp(ibd,ilp),tlevel,dt,ibd)
          endif

          do i=limp(ibd,ilp)+1,limp(ibd,ilp)+mimp(ibd,ilp)

            a(:,:) = pmtrx(:,:,i)
            indx(:) = pindx(:,i)

!             in = ip(i)          !!!!!!!
!             jn = jp(i)          !!!!!!!
!             kn = kp(i)          !!!!!!!

            b(1) = dpdnn(i)

            if(mrkp(i)==1 .OR. mrkp(i)==3) then
c the pressure at the extension fluid point is evaluated by interpolation
              pim(i) = interp_cellface(nxp(i),nyp(i),nzp(i)
     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dirp(i))
c              pim(i) = interp_cellface_mod(nxp(i),nyp(i),nzp(i)
c     &            ,xc,yc,zc,p,nx,ny,nz,icyl,dirp(i))
            else
              i1 = ip(i)+int(nxp(i))
              j1 = jp(i)+int(nyp(i))
              k1 = kp(i)+int(nzp(i))
              call index_bnd(i1,j1,k1,i1,j1,k1,nx,ny,nz)
c the pressure at the exterior grid point is determined
              pim(i) = p(i1,j1,k1)
            endif
            b(2) = pim(i)

c the local system for the pressure stencil is solved
            call lubksb(a,nsp,nsp,indx,b)
c the pressure condition at the interface point is enforced
            p(ip(i),jp(i),kp(i))=b(1)
            
            if(kp(i).eq.kz1) then
              ilpr = ilpr + 1
              indexl(ilpr,1) = ip(i)
              indexl(ilpr,2) = jp(i)
            endif
            if(kp(i).eq.kz2) then
              irpr = irpr + 1
              indexr(irpr,1) = ip(i)
              indexr(irpr,2) = jp(i)
            endif

          enddo

        enddo

        call local_refreshbc(indexl(1:ilpr,:),ilpr,indexr(1:irpr,:),irpr,p,nx,ny,nz)

!        call refreshbc(p,nx*ny,nz)

        !Periodic y bc
        p(:,1,:) = p(:,ny-1,:)
        p(:,ny,:) = p(:,2,:)

        !Centerline bc
        if(icyl==1) then
          do j=1,ny
            p(1,j,:) = p(2,jsym(j),:)
          enddo
        endif

      enddo

      !Periodic y bc
c      p(:,1,:) = p(:,ny-1,:)
c      p(:,ny,:) = p(:,2,:)

      !Centerline bc
c      if(icyl==1) then
c        do j=1,ny
c          p(1,j,:) = p(2,jsym(j),:)
c        enddo
c      endif

      return

      end
c---------------------------------------------------------------------

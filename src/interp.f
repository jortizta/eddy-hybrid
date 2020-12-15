C---- function interp_face-------------------N. Beratlis-05 Dec. 2008-
C
      REAL function interp_cellface(xp,yp,zp,x,y,z,p,nx,ny,nz,icyl,dir)
C
C     PURPOSE: Interpolate on a plane parallel to face of cell
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
c      include 'common.h'
c      include 'mpif.h'

      INTEGER nx,ny,nz,icyl,dir
      REAL    xp,yp,zp
      REAL    x(nx),y(ny),z(nz)
      REAL    p(nx,ny,nz)

      INTEGER i,j,k,plane,i1,j1,k1,i2,j2,k2,i3,j3,k3,i4,j4,k4,s
      REAL    pint,xint,yint,zint,a1,b1,c1
      REAL    interp2d,anglerad
     

      IF(icyl==1) THEN
        xint = sqrt(xp**2. + yp**2.)
        yint = anglerad(xp,yp)
      ELSE
        xint = xp
        yint = yp
      ENDIF
      zint = zp


      IF(abs(dir)==1) THEN
        call closest(x,nx,xint,i)
        call locate(y,ny,yint,j)
        call locate(z,nz,zint,k)
        b1 = (yint-y(j))/(y(j+1)-y(j))
        c1 = (zint-z(k))/(z(k+1)-z(k))
        call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
        call index_bnd(i  ,j+1,k  ,i2,j2,k2,nx,ny,nz)
        call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
        call index_bnd(i  ,j+1,k+1,i4,j4,k4,nx,ny,nz)
        pint = (1-b1)*(1-c1)*p(i1,j1,k1) + b1*(1-c1)*p(i2,j2,k2)
     &       + (1-b1)*  c1  *p(i3,j3,k3) + b1*  c1  *p(i4,j4,k4)
      ELSEIF(abs(dir)==2) THEN
        call closest(y,ny,yint,j)
        call locate(x,nx,xint,i)
        call locate(z,nz,zint,k)
        a1 = (xint-x(i))/(x(i+1)-x(i))
        c1 = (zint-z(k))/(z(k+1)-z(k))
        call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
        call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
        call index_bnd(i  ,j  ,k+1,i3,j3,k3,nx,ny,nz)
        call index_bnd(i+1,j  ,k+1,i4,j4,k4,nx,ny,nz)
        pint = (1-a1)*(1-c1)*p(i1,j1,k1) + a1*(1-c1)*p(i2,j2,k2)
     &       + (1-a1)*  c1  *p(i3,j3,k3) + a1*  c1  *p(i4,j4,k4)
      ELSEIF(abs(dir)==3) THEN
        call closest(z,nz,zint,k)
        call locate(x,nx,xint,i)
        call locate(y,ny,yint,j)
        a1 = (xint-x(i))/(x(i+1)-x(i))
        b1 = (yint-y(j))/(y(j+1)-y(j))
        call index_bnd(i  ,j  ,k  ,i1,j1,k1,nx,ny,nz)
        call index_bnd(i+1,j  ,k  ,i2,j2,k2,nx,ny,nz)
        call index_bnd(i  ,j+1,k  ,i3,j3,k3,nx,ny,nz)
        call index_bnd(i+1,j+1,k  ,i4,j4,k4,nx,ny,nz)
        pint = (1-a1)*(1-b1)*p(i1,j1,k1) + a1*(1-b1)*p(i2,j2,k2)
     &       + (1-a1)*  b1  *p(i3,j3,k3) + a1*  b1  *p(i4,j4,k4)
      ENDIF

      interp_cellface = pint

      RETURN

      END
C-----------------------------------------------------------------------






C---- function interp2d------------------------N. Beratlis-26 Dec. 2008-
C
      REAL function interp2d(xp,yp,x,y,p,nx,ny)
C
C     PURPOSE: Interpolate in 2 dimensions
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER nx,ny
      REAL    xp,yp
      REAL    x(nx),y(ny),p(nx,ny)

      INTEGER i,j
      REAL    a1,b1

      CALL LOCATE(x,nx,xp,i)
      CALL LOCATE(y,ny,yp,j)

      a1 = (xp-x(i))/(x(i+1)-x(i))
      b1 = (yp-y(j))/(y(j+1)-y(j))
      
      interp2d = (1.-a1)*(1.-b1)*p(i  ,j  )
     &         + (1.-a1)*(  b1 )*p(i  ,j+1)
     &         + (  a1 )*(1.-b1)*p(i+1,j  )
     &         + (  a1 )*(  b1 )*p(i+1,j+1)   

      RETURN
      END
C-----------------------------------------------------------------------


c-----function interp3d------------------------N. Beratlis 19 Dec. 2008--
C
      REAL function interp3d(x,y,z,xc,yc,zc,p,nx,ny,nz,icyl)
c
C     PURPOSE: Interpolate value in a point within cell i-j-k
c
C------------------------------------------------------------------------
c
      IMPLICIT NONE
c
c.... Input/Output Array
      INTEGER nx,ny,nz,icyl
      REAL    x,y,z
      REAL    xc(nx),yc(ny),zc(nz)
      REAL    p(nx,ny,nz)
c
c.... Local Arrays
      INTEGER i,j,k
      REAL    a1,b1,c1,a2,b2,c2,dx,dy,dz,pint,r,theta
c
c.... Function
      real    anglerad

      interp3d=0.0

      CALL ijk_xyz(x,y,z,xc,yc,zc,nx,ny,nz,i,j,k,icyl)
      dx = xc(i+1)-xc(i)
      dy = yc(j+1)-yc(j)
      dz = zc(k+1)-zc(k)

      if(icyl==0) then
        a1 = (x-xc(i))/dx
        a2 = 1.-a1
        b1 = (y-yc(j))/dy
        b2 = 1.-b1
      else
        r = sqrt(x**2.+y**2.)
        theta = anglerad(x,y)
        a1 = (r-xc(i))/dx
        a2 = 1.-a1
        b1 = (theta-yc(j))/dy
        b2 = 1.-b1
        if(a1<0.0 .OR. a1>1.0) then
          write(6,*) 'WARNING: a1 out of bounds:',a1,x,y,z,r,theta,xc(i),yc(j),i,j,k
        endif
        if(b1<0.0 .OR. b1>1.0) then
          write(6,*) 'WARNING: b1 out of bounds:',b1,x,y,z,r,theta,xc(i),yc(j),i,j,k
         endif
      endif
      c1 = (z-zc(k))/dz
      c2 = 1.-c1

      pint = a2*b2*c2*p(i  ,j  ,k  )
     &     + a1*b2*c2*p(i+1,j  ,k  )
     &     + a2*b1*c2*p(i  ,j+1,k  )
     &     + a1*b1*c2*p(i+1,j+1,k  )
     &     + a2*b2*c1*p(i  ,j  ,k+1)
     &     + a1*b2*c1*p(i+1,j  ,k+1)
     &     + a2*b1*c1*p(i  ,j+1,k+1)
     &     + a1*b1*c1*p(i+1,j+1,k+1)

c      write(6,*) i,j,k,a1,b1,c1,a2,b2,c2,r-xc(i),dx,pint

      interp3d = pint

      return

      end
C-----------------------------------------------------------------------


c-----function interp3d1 ----------------------N. Beratlis 19 Dec. 2008--
C
      REAL function interp3d1(x,y,z,xc,yc,zc,p,nx,ny,nz,icyl,ie)
c
C     PURPOSE: Interpolate value in a point within cell i-j-k
c
C------------------------------------------------------------------------
c
      IMPLICIT NONE
c
c.... Input/Output Array
      INTEGER nx,ny,nz,icyl,ie
      REAL    x,y,z
      REAL    xc(nx),yc(ny),zc(nz)
      REAL    p(nx,ny,nz)
c
c.... Local Arrays
      INTEGER i,j,k
      REAL    a1,b1,c1,a2,b2,c2,dx,dy,dz,pint,r,theta
c
c.... Function
      real    anglerad

      interp3d1=0.0

      CALL ijk_xyz(x,y,z,xc,yc,zc,nx,ny,nz,i,j,k,icyl)

      dx = xc(i+1)-xc(i)
      dy = yc(j+1)-yc(j)
      dz = zc(k+1)-zc(k)

      if(icyl==0) then
        a1 = (x-xc(i))/dx
        a2 = 1.-a1
        b1 = (y-yc(j))/dy
        b2 = 1.-b1
      else
        r = sqrt(x**2.+y**2.)
        theta = anglerad(x,y)
        a1 = (r-xc(i))/dx
        a2 = 1.-a1
        b1 = (theta-yc(j))/dy
        b2 = 1.-b1
        if(a1<0.0 .OR. a1>1.0) then
          write(6,*) 'WARNING: a1 out of bounds:',a1,x,y,z,r,theta,xc(i),yc(j)
        endif
        if(b1<0.0 .OR. b1>1.0) then
          write(6,*) 'WARNING: b1 out of bounds:',b1,x,y,z,r,theta,xc(i),yc(j)
        endif
      endif
      c1 = (z-zc(k))/dz
      c2 = 1.-c1

      pint = (1.-a1)*(1.-b1)*(1.-c1)*p(i  ,j  ,k  )
     &     + (  a1 )*(1.-b1)*(1.-c1)*p(i+1,j  ,k  )
     &     + (1.-a1)*(  b1 )*(1.-c1)*p(i  ,j+1,k  )
     &     + (  a1 )*(  b1 )*(1.-c1)*p(i+1,j+1,k  )
     &     + (1.-a1)*(1.-b1)*(  c1 )*p(i  ,j  ,k+1)
     &     + (  a1 )*(1.-b1)*(  c1 )*p(i+1,j  ,k+1)
     &     + (1.-a1)*(  b1 )*(  c1 )*p(i  ,j+1,k+1)
     &     + (  a1 )*(  b1 )*(  c1 )*p(i+1,j+1,k+1)

      interp3d1 = pint

      return

      end
C-----------------------------------------------------------------------



C---- function interp_face-------------------N. Beratlis-05 Dec. 2008-
C
      REAL function interp_cellface1(xp,yp,zp,x,y,z,p,nx,ny,nz,icyl,dir)
C
C     PURPOSE: Interpolate on a plane parallel to face of cell
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
c      include 'common.h'
c      include 'mpif.h'

      INTEGER nx,ny,nz,icyl,dir
      REAL    xp,yp,zp
      REAL    x(nx),y(ny),z(nz)
      REAL    p(nx,ny,nz)

      INTEGER i,j,k,plane
      REAL    pint,xint,yint,zint
      REAL    interp2d,anglerad
      REAL, DIMENSION(:,:), ALLOCATABLE :: p2d

      IF(icyl==1) THEN
        xint = sqrt(xp**2. + yp**2.)
        yint = anglerad(xp,yp)
      ELSE
        xint = xp
        yint = yp
      ENDIF
      zint = zp

      IF(abs(dir)==1) THEN
        call closest(x,nx,xint,plane)
        ALLOCATE(p2d(ny,nz))
        DO j=1,ny
        DO k=1,nz
          p2d(j,k) = p(plane,j,k)
        ENDDO
        ENDDO
        pint = interp2d(yint,zint,y,z,p2d,ny,nz)
      ELSEIF(abs(dir)==2) THEN
        call closest(y,ny,yint,plane)
        ALLOCATE(p2d(nx,nz))
        DO i=1,nx
        DO k=1,nz
          p2d(i,k) = p(i,plane,k)
        ENDDO
        ENDDO
        pint = interp2d(xint,zint,x,z,p2d,nx,nz)
      ELSEIF(abs(dir)==3) THEN
        call closest(z,nz,zint,plane)
        ALLOCATE(p2d(nx,ny))
        DO i=1,nx
        DO j=1,ny
          p2d(i,j) = p(i,j,plane)
        ENDDO
        ENDDO
        pint = interp2d(xint,yint,x,y,p2d,nx,ny)
      ENDIF

      DEALLOCATE(p2d)

      interp_cellface1 = pint

c      write(6,*) xp,yp,zp,dir,pint

      RETURN

      END
C-----------------------------------------------------------------------



C---- function triangle_interp ---------------N. Beratlis-26 Jun 2009---
C
C     PURPOSE: Triangle interpolation
C
C-----------------------------------------------------------------------
      real function triangle_interp(x,y,z,node,var,ie)
c
      implicit none
c
c.... Input/Output Arrays
      real    x,y,z
      real    node(3,3)
      real    var(3)
c
c.... Local arrays
      integer dir,ie
      real    a,b,c,d,xp,yp,p
      real    va(2),vb(2),vc(2),un(3)
      real    amtrx(3,3),indx(3),rmtrx(3)
c
c.... Function
      real    vecmag,dotproduct

      if(.false.) then
      rmtrx(1) = sqrt( (node(1,1)-x)**2. + (node(1,2)-y)**2. + (node(1,3)-z)**2.)
      rmtrx(2) = sqrt( (node(2,1)-x)**2. + (node(2,2)-y)**2. + (node(2,3)-z)**2.)
      rmtrx(3) = sqrt( (node(3,1)-x)**2. + (node(3,2)-y)**2. + (node(3,3)-z)**2.)

c      d = sum(rmtrx)
c      a = 2.0 - (rmtrx(1))/d
c      b = 2.0 - (rmtrx(2))/d
c      c = 2.0 - (rmtrx(3))/d
      a = rmtrx(2)*rmtrx(3)
      b = rmtrx(1)*rmtrx(3)
      c = rmtrx(1)*rmtrx(2)
      d = a+b+c

      a = a/d
      b = b/d
      c = c/d
      
      p = a*var(1) + b*var(2) + c*var(3)
c      p = (var(1)+var(2)+var(3))/3.0

c      if(ie==103) then
c        write(6,*) ie,a,b,c,rmtrx,var,p
c      endif

      triangle_interp = p
      return
      endif

c      write(6,*) rmtrx,var,p
c      stop

      !Compute normal vector
      call plane_tri(node(:,1),node(:,2),node(:,3),a,b,c,d)
      un(1) = a
      un(2) = b
      un(3) = c
      d = vecmag(un,3)
c      write(6,*) 'un=',un,', vecmag=',d      
      un = un/d

      dir = maxloc(abs(un),1)
c      write(6,*) 'dir=',dir,abs(un)

      !Project triangle along direction most aligned with its normal
      if(dir==1) then
        va(1) = node(1,2)
        va(2) = node(1,3)
        vb(1) = node(2,2)
        vb(2) = node(2,3)
        vc(1) = node(3,2)
        vc(2) = node(3,3)
        xp = y
        yp = z
      elseif(dir==2) then
        va(1) = node(1,1)
        va(2) = node(1,3)
        vb(1) = node(2,1)
        vb(2) = node(2,3)
        vc(1) = node(3,1)
        vc(2) = node(3,3)
        xp = x
        yp = z
      else
        va(1) = node(1,1)
        va(2) = node(1,2)
        vb(1) = node(2,1)
        vb(2) = node(2,2)
        vc(1) = node(3,1)
        vc(2) = node(3,2)
        xp = x
        yp = y
      endif

      amtrx(1,1) = va(1)
      amtrx(2,1) = vb(1)
      amtrx(3,1) = vc(1)
      amtrx(1,2) = va(2)
      amtrx(2,2) = vb(2)
      amtrx(3,2) = vc(2)
      amtrx(:,3) = 1.0

      call ludcmp(amtrx,3,3,indx,d)
      rmtrx(1:3) = var(1:3)
      call lubksb(amtrx,3,3,indx,rmtrx)
      
      p = rmtrx(1)*xp + rmtrx(2)*yp + rmtrx(3)
c      if(ie==255) then
c      write(6,*) dir,rmtrx,xp,yp,va,vb,vc,var,p
c      write(6,*) node,var
c      write(6,*) rmtrx(1)*va(1)+rmtrx(2)*va(2)+rmtrx(3),var(1)
c      write(6,*) rmtrx(1)*vb(1)+rmtrx(2)*vb(2)+rmtrx(3),var(2)
c      write(6,*) rmtrx(1)*vc(1)+rmtrx(2)*vc(2)+rmtrx(3),var(3)

c      call mpi_finalize(ierr)
c      stop
c      endif

      triangle_interp = p

      return

      end
C-----------------------------------------------------------------------

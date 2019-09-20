c A rigid body rotation around Z is considered !!!!!!!
      SUBROUTINE RBM(UNVECT,VERTEX,VERTEXC,NFACET,MBD,NBD,TMOVE,TLEVEL)

      implicit none
      include 'immersed.h'
      REAL, PARAMETER :: pi=3.141592653589793
      real    tmove,tlevel
      integer nfacet,nbd,mbd,ibd
      real    unvect(3,nfacet),vertex(3,3,nfacet),vertexc(3,nfacet)
c
c.... Local arrays
      integer i,ilb,ile,ivtx
      real    dtheta,x,y,z,omega

c.... Functions
      real    rotspeed
!      real    osc_freq

!     no change in Z direction
!      OMEGA = rotspeed(tlevel)
!      DTHETA = OMEGA*TMOVE

!     OMEGA = osc_freq(tlevel)
      DTHETA = 0.8*sin(0.8*tlevel)
!      write(*,*) DTHETA,OMEGA,tlevel,mbd,nbd
      DO ibd=1,nbd
        ilb = lb(ibd)+1
        ile = lb(ibd)+mb(ibd)       
!        write(*,*)ilb,ile 
        DO I=ILB,ILE
c...Normal vector
          X = UNVECT(1,I)
!          Y = UNVECT(2,I)
!          Z = UNVECT(3,I)
c          UNVECT(1,I) = X*COS(DTHETA)+Z*SIN(DTHETA)
c          UNVECT(3,I) = Z*COS(DTHETA)-X*SIN(DTHETA)
!          UNVECT(1,I) = X*COS(DTHETA)-Y*SIN(DTHETA)
!          UNVECT(2,I) = Y*COS(DTHETA)+X*SIN(DTHETA)
!           UNVECT(3,I) = Z + DTHETA
           UNVECT(1,I) = X+DTHETA
!           WRITE(*,*) X,UNVECT(1,I)
c...Vertices
          DO IVTX=1,3
            X = VERTEX(1,IVTX,I)
            Y = VERTEX(2,IVTX,I)
            Z = VERTEX(3,IVTX,I)
c            VERTEX(1,IVTX,I) = X*COS(DTHETA)+Z*SIN(DTHETA)
c            VERTEX(3,IVTX,I) = Z*COS(DTHETA)-X*SIN(DTHETA)
!            VERTEX(1,IVTX,I) = X*COS(DTHETA)-Y*SIN(DTHETA)
!            VERTEX(2,IVTX,I) = Y*COS(DTHETA)+X*SIN(DTHETA)
             VERTEX(1,IVTX,I) = X+DTHETA
!            VERTEX(3,IVTX,I) = Z+DTHETA

          ENDDO

          VERTEXC(:,I)=SUM(VERTEX(:,:,I),2)/3.0

        ENDDO

c        write(6,*) '1. rbm:',vertex(:,:,1),vertexc(:,1)
c        write(6,*) '2. rbm:',vertex(:,:,1),sum(vertex(:,:,1),2)/3.0

      ENDDO

      RETURN

      END
c---------------------------------------------------------------------


c A rigid body rotation around Z is considered !!!!!!!
C---- SUBROUTINE RBM_NODES -------------- N. Beratlis-01 Nov. 2010 ---
C
C     PURPOSE: Move nodes and normals of triangles
C
c---------------------------------------------------------------------
      subroutine rbm_nodes(node,nd_nrm,np,mbd,nbd,dt,tlevel)
c
      implicit none
      include 'immersed.h'
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output arrays
      integer np,mbd,nbd
      real    dt,tlevel
      real    node(3,np),nd_nrm(3,np)
c
c.... Local arrays
      integer ibd,i,ilb,ile
      real    omega,dtheta,x,y,z
c
c.... Functions
      real    rotspeed

      OMEGA = rotspeed(tlevel)
      DTHETA = OMEGA*DT

      DO ibd=mbd,nbd
        ilb = lv(ibd)+1
        ile = lv(ibd)+mv(ibd)        
        DO I=ILB,ILE
c...Normal vector
          X = ND_NRM(1,I)
          Y = ND_NRM(2,I)          
c          Z = ND_NRM(3,I)
c          ND_NRM(1,I) = X*COS(DTHETA)+Z*SIN(DTHETA)
c          ND_NRM(3,I) = Z*COS(DTHETA)-X*SIN(DTHETA)
          ND_NRM(1,I) = X*COS(DTHETA)-Y*SIN(DTHETA)
          ND_NRM(2,I) = Y*COS(DTHETA)+X*SIN(DTHETA)          
c... Vertices
          X = NODE(1,I)
          Y = NODE(2,I)
c          Z = NODE(3,I)
c          NODE(1,I) = X*COS(DTHETA)+Z*SIN(DTHETA)
c          NODE(3,I) = Z*COS(DTHETA)-X*SIN(DTHETA)
          NODE(1,I) = X*COS(DTHETA)-Y*SIN(DTHETA)
          NODE(2,I) = Y*COS(DTHETA)+X*SIN(DTHETA)          
        ENDDO
      ENDDO

      return

      end
c---------------------------------------------------------------------



c---- FUNCTION OMEGA -------------------------------------------------
c
c     PURPOSE: Compute rotational speed of a spinning sphere
c
c---------------------------------------------------------------------
c
      REAL function rotspeed(t)

      INCLUDE 'common.h'
      
      REAL t,mag

      !Golfball a=omega*R/U
      !omega = 1500rpm*2*pi/60
      !R=D/2=43.0e-3/2
      !U=40m/s ->Re=110k

c      mag = 0.1464 !2500rpm, Re=110k
c      mag = 0.19646 !3500rpm, Re=110k
      mag = amp !3500rpm, Re=17000

      rotspeed = mag

      RETURN

      END
c
c---------------------------------------------------------------------


c---- function ubd-------------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set u component of immersed body velocity. 
C     ATTENTION: x,y,z are assumed in cartesian coordinates, but ubd
C     can be cylindrical or cartesian velocity.
c
c     A rigid rotation around Z is considered:
c     the radial velocity is equal to 0 !!!!!!!
C
c---------------------------------------------------------------------
      real function ubd(x,y,z,t,ibd)

c      include 'common.h'
c      include 'immersed.h'
      implicit none
c
c.... Input/Output arrays
      integer ibd
      real x,y,z,t
c
c.... Local arrays
c      real xcar,ycar,r,a
c      real omega,phi,theta,U,ux
c
c.... Function
c      real rotspeed,anglerad

!      ubd = 0.0
       ubd = 0.8*sin(0.8*t)

c      omega = rotspeed(t)

c      r = sqrt(x**2. + z**2.)
c      U = omega*r

c      if(r>0.0) then
c        ux = U*z/r
c        r = sqrt(x**2. + y**2.)
c        if(r>0.0) ubd = ux*x/r
c      endif

      RETURN

      END
c---------------------------------------------------------------------



c---- function vbd-------------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set v component of immersed body velocity. 
C     ATTENTION: x,y,z are assumed in cartesian coordinates, but vbd
C     can be cylindrical or cartesian velocity.
c
c     A rigid rotation around Z is considered:
c     the azimuthal velocity is equal to r*omega !!!!!!!
C
c---------------------------------------------------------------------
      real function vbd(x,y,z,t,ibd)

c      include 'common.h'
c      include 'immersed.h'
      implicit none
c
c.... Input/Output arrays
      integer ibd
      real x,y,z,t
c
c.... Local arrays
c      real xcar,ycar,r,a
c      real omega,phi,theta,U,ux
      real r,omega
c
c.... Function
c      real rotspeed,anglerad
      real rotspeed

c      vbd = 0.0

      omega = rotspeed(t)

      r = sqrt(x**2. + y**2.)
      vbd = omega*r

c      if(r>0.0) then
c        ux = U*z/r
c        r = sqrt(x**2. + y**2.)
c        if(r>0.0) vbd =-ux*y/r
c      endif


      RETURN

      END
c---------------------------------------------------------------------



c---- function wbd-------------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set w component of immersed body velocity. 
C     ATTENTION: x,y,z are assumed in cartesian coordinates, but wbd
C     can be cylindrical or cartesian velocity.
c
c     A rigid rotation around Z is considered:
c     the axial velocity is equal to 0 !!!!!!!
C
c---------------------------------------------------------------------
      real function wbd(x,y,z,t,ibd)

c      include 'common.h'
c      include 'immersed.h'
      implicit none
c
c.... Input/Output arrays
      integer ibd
      real x,y,z,t
c
c.... Local arrays
c      real xcar,ycar,r,a
c      real omega,phi,theta,U,ux
c
c.... Function
c      real rotspeed,anglerad

      wbd = 0.0
!       wbd  = 0.8*sin(0.8*t)
c      omega = rotspeed(t)

c      r = sqrt(x**2. + z**2.)
c      U = omega*r

c      if(r>0.0) then
c        wbd =-U*x/r
c      endif

      RETURN

      END
c---------------------------------------------------------------------



c---- function dpdn------------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set pressure gradient along normal to immersed body. 
c     ATTENTION: x,y,z are assumed in cartesian coordinates.
c
c     ATTENTION: a rigid body rotation around Z is considered !!!!!!
c
c---------------------------------------------------------------------
      real function dpdn(x,y,z,nx,ny,nz,t,dt,ibd)

      implicit none
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output arrays      
      integer ibd
      real    x,y,z,nx,ny,nz,t,dt
c
c.... Local arrays
      real    omega,dtheta,un1,un2,dundt,x2,y2,r,U
      real    u1(3),u2(3),nv1(3),nv2(3)
      real    un(3),dp(3)
c
c.... Functions
      real    rotspeed,dotproduct

      dpdn = 0.0

      !Find velocities in cartesian coordinates
      r = sqrt(x**2. + y**2.)
      omega = rotspeed(t)
      U = omega*r

      if(r>0.0) then
        u1(1) =-U*y/r
        u1(2) = U*x/r
        u1(3) = 0.0
        
        nv1(1) = nx
        nv1(2) = ny
        nv1(3) = nz

        un1 = dotproduct(u1,nv1,3)
c UN1 is the dot product of the normal vector and the velocity of
c the body before the rotation
      endif

      !Perform rigid body motion
      dtheta = omega*dt
      x2 = x*cos(dtheta)-y*sin(dtheta)
      y2 = y*cos(dtheta)+x*sin(dtheta)
      nv2(1) = nx*cos(dtheta)-ny*sin(dtheta)
      nv2(2) = ny*cos(dtheta)+nx*sin(dtheta)
      nv2(3) = nz

      r = sqrt(x2**2. + y2**2.)
      omega = rotspeed(t+dt)
      U = omega*r

      if(r>0.0) then
        u2(1) =-U*y2/r
        u2(2) = U*x2/r
        u2(3) = 0.0
        
        un2 = dotproduct(u2,nv2,3)
c UN2 is the dot product of the normal vector and the velocity of
c the body after the rotation
      endif

      dundt = (un2-un1)/dt

      dpdn = -dundt

      return

      end
c---------------------------------------------------------------------




c---- function dpdx1 ----------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set pressure gradient along x direction in cartesian coordinates. 
c     ATTENTION: x,y,z are assumed in cartesian coordinates.
c
c     ATTENTION: a rigid body rotation around Z is considered !!!!!!
c
c---------------------------------------------------------------------
      real function dpdx1(x,y,z,t,dt,ibd)

      implicit none
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output arrays      
      integer ibd
      real    x,y,z,t,dt
c
c.... Local arrays
      real    omega,dtheta,ux1,ux2,duxdt,x2,y2,r,U
      real    u1(3),u2(3),nv1(3),nv2(3)
      real    un(3),dp(3)
c
c.... Functions
      real    rotspeed,dotproduct

      dpdx1 = 0.0

      !Find velocities in cartesian coordinates
      r = sqrt(x**2. + y**2.)

      if(r>0.0) then
        omega = rotspeed(t)
        U = omega*r
        ux1 =-U*y/r

        !Perform rigid body motion and advance body in time
        dtheta = omega*dt
        y2 = y*cos(dtheta) + x*sin(dtheta)

        !Compute velocity in new position
        omega = rotspeed(t+dt)
        U = omega*r
        ux2 =-U*y2/r

        duxdt = (ux2-ux1)/dt
        dpdx1 = -duxdt

      endif

      return

      end
c---------------------------------------------------------------------



c---- function dpdx2 ----------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set pressure gradient along y direction in cartesian coordinate. 
c     ATTENTION: x,y,z are assumed in cartesian coordinates.
c
c     ATTENTION: a rigid body rotation around Z is considered !!!!!!
c
c---------------------------------------------------------------------
      real function dpdx2(x,y,z,t,dt,ibd)

      implicit none
c
c.... Input/Output arrays      
      integer ibd
      real    x,y,z,t,dt
c
c.... Local arrays
      real    U,r,x2,uy1,uy2,duydt,dtheta,omega
c
c.... Functions
      real     rotspeed

      dpdx2 = 0.0

      !Find velocities in cartesian coordinates
      r = sqrt(x**2. + y**2.)

      if(r>0.0) then
        omega = rotspeed(t)
        U = omega*r
        uy1 = U*x/r

        !Perform rigid body motion and advance body in time
        dtheta = omega*dt
        x2 = x*cos(dtheta) - y*sin(dtheta)    !!!!!! - instead of +

        !Compute velocity in new position
        omega = rotspeed(t+dt)
        U = omega*r
        uy2 =U*x2/r

        duydt = (uy2-uy1)/dt
        dpdx2 = -duydt

      endif

      return

      end
c---------------------------------------------------------------------



c---- function dpdx3 ----------------------N. Beratlis-10 Jul. 2009---
c
c     PURPOSE: Set pressure gradient along z direction (cartesian). 
c     ATTENTION: x,y,z are assumed in cartesian coordinates.
c
c     ATTENTION: a rigid body rotation around Z is considered !!!!!!
c
c---------------------------------------------------------------------
      real function dpdx3(x,y,z,t,dt,ibd)

      implicit none
      REAL, PARAMETER :: pi=3.141592653589793
c
c.... Input/Output arrays      
      integer ibd
      real    x,y,z,t,dt

      dpdx3 = 0.0

      return

      end
c---------------------------------------------------------------------


c---- function osc_freq --------------------N. Beratlis-3 Aug. 2009---
c
c     PURPOSE: Compute oscillation frequency.
c
c---------------------------------------------------------------------
c
      real function osc_freq(t)

      implicit none

      real t

      osc_freq = 0.8

      return

      end
c---------------------------------------------------------------------



c---- function osc_amp --------------------N. Beratlis-3 Aug. 2009---
c
c     PURPOSE: Compute oscillation amplitude.
c
c---------------------------------------------------------------------
c
      real function osc_amp(t)

      implicit none

      real t

      osc_amp = 0.8

      return

      end
c---------------------------------------------------------------------

c---- function vbdm------------------------N. Beratlis-10 Jul. 2009---
c---------------------------------------------A. Posa - August 2012---
c
c     PURPOSE: Set v component of immersed body velocity. 
C     ATTENTION: x,y,z are assumed in cartesian coordinates, but vbd
C     can be cylindrical or cartesian velocity.
c
c     A rigid rotation around Z is considered:
c     the azimuthal velocity increases linearly  
C
c---------------------------------------------------------------------
      real function vbdm(x,y,z,t,ibd,myrank,ind)

c      include 'common.h'
c      include 'immersed.h'
      implicit none
c
c.... Input/Output arrays
      integer ibd,myrank,ind
      real x,y,z,t,t0,t1
      parameter(t0=83.78,t1=97.22)
c
c.... Local arrays
c      real xcar,ycar,r,a
c      real omega,phi,theta,U,ux
      real r,omega
c
c.... Function
c      real rotspeed,anglerad
      real rotspeed

c      vbdm = 0.0

      omega = min(1.,(t-t0)/(t1-t0))*rotspeed(t)

      r = sqrt(x**2. + y**2.)
      vbdm = omega*r
      if((myrank.eq.0).and.(ind.eq.1000))write(6,*)'rot_speed [rad/s] =',omega

c      if(r>0.0) then
c        ux = U*z/r
c        r = sqrt(x**2. + y**2.)
c        if(r>0.0) vbdm =-ux*y/r
c      endif


      RETURN

      END
c---------------------------------------------------------------------

c---- function vbdm2-----------------------N. Beratlis-10 Jul. 2009---
c---------------------------------------------A. Posa - August 2012---
c
c     PURPOSE: Set v component of immersed body velocity. 
C     ATTENTION: x,y,z are assumed in cartesian coordinates, but vbd
C     can be cylindrical or cartesian velocity.
c
c     A rigid rotation around Z is considered:
c     the azimuthal velocity increases linearly  
C
c---------------------------------------------------------------------
      real function vbdm2(x,y,z,t,ibd,icycle)

c      include 'common.h'
c      include 'immersed.h'
      implicit none
c
c.... Input/Output arrays
      integer ibd,icycle,i0,i1
      real x,y,z,t
      parameter(i0=0,i1=1000)
c
c.... Local arrays
c      real xcar,ycar,r,a
c      real omega,phi,theta,U,ux
      real r,omega
c
c.... Function
c      real rotspeed,anglerad
      real rotspeed

c      vbdm2 = 0.0

      omega = min(1.,real(icycle-i0)/real(i1-i0))*rotspeed(t)

      r = sqrt(x**2. + y**2.)
      vbdm2 = omega*r

c      if(r>0.0) then
c        ux = U*z/r
c        r = sqrt(x**2. + y**2.)
c        if(r>0.0) vbdm2 =-ux*y/r
c      endif


      RETURN

      END
c---------------------------------------------------------------------


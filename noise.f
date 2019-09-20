c------------------------------------ A. Posa - July 2012 ------------

      subroutine noise(nx,ny,nz,jy1,jy2,kz1,kz2,u,v,w,myrank,iflag)

      implicit none

      integer nx,ny,nz,jy1,jy2,kz1,kz2,myrank,iflag
      integer i1,i2,i,j,k,k1g,k2g
      parameter(i1=175,i2=180,k1g=332,k2g=337)
      real u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)  
      integer, save :: idum = 1
      integer, save :: k1,k2
      real mean,stddev,gasdev

      parameter(mean=0.,stddev=0.1)
 
      if(iflag.eq.0) then

        k1=k1g-myrank*(nz-2)
        if(k1.le.1) then
          k1=kz1
        elseif(k1.ge.nz) then
!          return
          k1=nz
        endif
        k2=k2g-myrank*(nz-2)
        if(k2.ge.nz) then
          k2=kz2
        elseif(k2.le.1) then
!          return
          k2=1
        endif

      else

        do 100 k=k1,k2
        do 100 j=jy1,jy2
        do 100 i=i1,i2
 100    u(i,j,k)=u(i,j,k)+gasdev(idum,mean,stddev)      

        do 101 k=k1,k2
        do 101 j=jy1,jy2
        do 101 i=i1,i2
 101    v(i,j,k)=v(i,j,k)+gasdev(idum,mean,stddev)   

        do 102 k=k1,k2
        do 102 j=jy1,jy2
        do 102 i=i1,i2
 102    w(i,j,k)=w(i,j,k)+gasdev(idum,mean,stddev)   

      endif

      return
      end

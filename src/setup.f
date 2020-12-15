
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE SETUP
*                                                         *
*      Fills in all arrays containing information         *
*      about starting and ending indexes of each          *
*       boundary                                          *
*                                                         *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INCLUDE 'common.h'
*
**** Starting and ending indexes for each boundary
**** Starting and ending indexes for u,v,w momentum
*
***** Boundary "itype(1)", (west)
*
      ib1 = ix1-1
      ibu = ix1
c
      ibc = 1
      if(itype(1)==50.or.itype(1)==300) then
        ibc(1:3,1) =  1
      elseif(itype(1)==60) then
        ibc(1  ,1) =  1
        ibc(2:3,1) = -1 
      elseif(itype(1)==81) then
        ibc(1  ,1) =  1
        ibc(2  ,1) = -1 
        ibc(3  ,1) =  1
      endif
*
**** Boundary "itype(2)", (east)
*
      ib2 = ix2+1
      ieu = ix2-1
      if(itype(2)==500) ieu = ix2

      if(itype(2)==50) then
        ibc(1:3,2) =  1
      elseif(itype(2)==60) then
        ibc(1  ,2) =  1
        ibc(2:3,2) = -1
      elseif(itype(2)==80) then
        ibc(1  ,2) =  1
        ibc(2:3,2) = -1 
      elseif(itype(2)==81) then
        ibc(1  ,2) =  1
        ibc(2  ,2) = -1 
        ibc(3  ,2) =  1
      endif
*
**** Boundary "itype(3)", (front)
*
      jb1 = jy1-1
      jbv = jy1
*
**** Boundary "itype(4)", (back)
*
      jb2 = jy2+1
      jev = jy2-1
      if(itype(4)==500) jev = jy2
*
**** Boundary "itype(5)", (bottom)
*
      kb1 = kz1-1
      kbw = kz1
*
**** Boundary "itype(6)", (top)
*
      kb2 = kz2+1
      kew = kz2-1
      if(itype(6)==500.or.itype(6)==0) kew = kz2
	write(*,*) kb2,kew,kz2,kb1,kbw
      RETURN
      END

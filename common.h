c
c-----INCLUDE-FILE-- common.h -----------------I. Balaras-20/05/1998----
c
c-----------------------------------------------------------------------
c
c         I. BALARAS
c         LAST UPDATE      :  29/11/1998
c
c-----------------------------------------------------------------------
c                                                             parameters
C-----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER, PARAMETER :: MMX=1000,MMY=514,MMZ=7000
      INTEGER, PARAMETER :: NIMPLXMAX=2,NIMPLYMAX=2
      REAL, PARAMETER :: pi=3.141592653589793
     
      INTEGER background_time,begin_time
      LOGICAL save_res
c
      CHARACTER*120  str1,str2,str3,str4,str5,str6
      CHARACTER     index(0:9999)*4,proc(0:999)*3
      INTEGER  icyl,implx,imply,implz,implcx,implcy,implcz,icalf,icfl
      INTEGER  implxdcmp,implydcmp,implzdcmp
      INTEGER  nimplx,nimply,nimplz
      INTEGER  implxulim(nimplxmax,2),implxvlim(nimplxmax,2),implxwlim(nimplxmax,2)
     &        ,implyulim(nimplymax,2),implyvlim(nimplymax,2),implywlim(nimplymax,2)
      INTEGER  ischm,ibm,ipres,isgs,ifield,idens,inoise
      INTEGER  igrid,jgrid,kgrid,iogrid
      INTEGER  ix1,ix2,jy1,jy2,kz1,kz2
      INTEGER  ib1,ib2,jb1,jb2,kb1,kb2
      INTEGER  ibu,ieu,jbv,jev,kbw,kew
      INTEGER  itype(6)
      INTEGER  mp,lp,np
      INTEGER  itmax,itini,itscr,istat,itser,itres,itpost,it5p,kmin5p,kmax5p,
     &         itpln,ivrtx,itbdy,itcalf,it2d,resstep,res1,res2,pstride
     &        ,itVPfield,nVPfield,nwrite,nwrite2d,nwritebdy,iolvl,ioclock
      INTEGER  i2dflag,i3dflag
      INTEGER  itmprb,itke,yspctr,itmavg,iregtmavg
c
      INTEGER  MYSIZE,MYRANK,MYLEFT,MYRIGHT,MPI_COMM_EDDY,IERR,MTYPE,MTYPE2,COMM3D
c      INTEGER ,PARAMETER :: smpistatus=MPI_STATUS_SIZE
      INTEGER  rankx1x2,rankx2x3,rankx1x3,commx1x2, commx2x3, commx1x3,sizex1x2,sizex2x3,sizex1x3,
     &         sizex1,sizex2,sizex3,commx1, commx2, commx3,rankx1,rankx2,rankx3
      INTEGER  comm3dp,comm3dl,comm3dc,myid,nbrx1m1,nbrx1p1,nbrx2m1,nbrx2p1,nbrx3m1,nbrx3p1
      INTEGER  dims(3),coords(3)
      LOGICAL  nprscr,periods(3)
      INTEGER  SX,EX,SY,EY,SZ,EZ,CEX,nxprocs,nyprocs,nzprocs
c      INTEGER, PARAMETER :: realtype=MPI_DOUBLE_PRECISION
c      INTEGER, PARAMETER :: inttype=MPI_INTEGER
      REAL  ap(mmx),au(mmx),av(mmx),aw(mmx),
     &      bp(mmy),bu(mmy),bv(mmy),bw(mmy),
     &      cp(mmz),cu(mmz),cv(mmz),cw(mmz),cpg(mmz),cwg(mmz),
     &      ruimplx(mmx),rvimplx(mmx),rwimplx(mmx),ruimply(mmx),rvimply(mmx),rwimply(mmx)
      REAL  xmin,xmax,ymin,ymax,zmin,zmax,xlen,ylen,zlen
     &     ,delx,dely,delz,delxsq,delysq,delzsq
      REAL  ru1,css,eps,dpdx,dpdy,dpdz,frn,prn,denP1,grav,rho_0,gt1,gt2,g_orig
      REAL  cflc,tstep,dtsave,tini
      REAL  alf(3),rho(3),gam(3)
      REAL  whx(1:mmx,-1:1),why(1:mmy,-1:1),whz(1:mmz,-1:1)
      REAL  ibc(3,6)
      REAL  fre,amp,cphs,kwav,dummyfp1,dummyfp2,dummyfp3,dummyfp4,dummyfp5,ampl
      REAL  umov1,umov2,vmov1,vmov2,wmov1,wmov2

      REAL  RU(MMX),RP(MMX)
      INTEGER JSYM(MMY),RANKSYM

      LOGICAL stretch
      INTEGER start_feed, end_feed, stride_feed,index_feed
      CHARACTER*600 path_feed


c-----------------------------------------------------------------------
c                                                                commons
c-----------------------------------------------------------------------
c
      COMMON /SOLVER/  ischm,ibm,isgs,ipres,ifield,idens,inoise
      COMMON /CONTROL/ icyl,implx,imply,implz,implcx,implcy,implcz,icalf,icfl
     &                ,implxdcmp,implydcmp,implzdcmp
      COMMON /IMPLIM/  implxulim,implxvlim,implxwlim,implyulim,implyvlim,implywlim,nimplx,nimply,nimplz
      COMMON /INTGRD/  igrid,jgrid,kgrid,iogrid,stretch
      COMMON /INDEX/   ix1,ix2,jy1,jy2,kz1,kz2,ib1,ib2,jb1,jb2,kb1,kb2
      COMMON /CMP/     xmin,xmax,ymin,ymax,zmin,zmax,xlen,ylen,zlen
     &                ,delx,dely,delz,delxsq,delysq,delzsq
      COMMON /MOMENT/  ibu,ieu,jbv,jev,kbw,kew
      COMMON /BOUND1/  itype
      COMMON /MTCOEF/  ap,au,av,aw,bp,bu,bv,bw,cp,cu,cv,cw,cpg,cwg
      COMMON /IMPCOEF/ ruimplx,rvimplx,rwimplx,ruimply,rvimply,rwimply
      COMMON /DAT1/    itmax,itini,itscr,istat,itser,itres,itpost,it5p,kmin5p,kmax5p,
     &                 itpln,ivrtx,itbdy,itcalf,it2d,itVPfield,resstep
     &                 ,nVPfield,nwrite,nwrite2d,nwritebdy,iolvl,ioclock,background_time
      COMMON /DAT2/    ru1,css,eps,dpdx,dpdy,dpdz,frn,prn,denP1,grav,rho_0,gt1,gt2,g_orig
      COMMON /IOFILE/  str1,str2,str3,str4,str5,str6
      COMMON /tmstep/  cflc,tstep,dtsave,tini
      COMMON /scheme/  alf,rho,gam
      COMMON /bound2/  ibc
      COMMON /bound3/  mp,lp,np
      COMMON /weight/  whx,why,whz
      COMMON /dumprm/  fre,amp,cphs,kwav,dummyfp1,dummyfp2,dummyfp3,dummyfp4,dummyfp5,ampl
      COMMON /movwall/ umov1,umov2,vmov1,vmov2,wmov1,wmov2
      COMMON /pulsint/ i2dflag,i3dflag
      COMMON /chindx/  index
      COMMON /chproc/  proc
      COMMON /stats/   itmprb,itke,yspctr,itmavg,iregtmavg

      COMMON /parall/  MYSIZE,MYRANK,MYLEFT,MYRIGHT,MPI_COMM_EDDY,IERR,MTYPE,MTYPE2,COMM3D,SX,EX,SY,EY,SZ,EZ,CEX,
     &                 rankx1x2,rankx2x3,rankx1x3,commx1x2, commx2x3, commx1x3,sizex1x2,sizex2x3,sizex1x3,
     &		       sizex1,sizex2,sizex3,commx1, commx2,nxprocs,nyprocs,nzprocs,commx3,rankx1,rankx2,rankx3

      COMMON /cyl/     ru,rp
      COMMON /cylsym/  jsym,ranksym
      COMMON /hybprm/  start_feed,end_feed,stride_feed,index_feed,path_feed

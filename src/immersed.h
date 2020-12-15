      integer nbdmax,nfcmax,nintmax
      integer nsp,nsm
      parameter (nbdmax=3)
      parameter (nsp=2,nsm=1)

      integer iprsdv,ivelrcn,nflu,nflp,idomy,imb3d,imbdcmp,nntr2d,nntr3d,igrdint,itagy,imls,izero
      real    imbovlp,uext1,uext2,uext3
      integer lb(nbdmax),mb(nbdmax),lv(nbdmax),mv(nbdmax)
      integer ibmin(nbdmax),ibmax(nbdmax),jbmin(nbdmax),jbmax(nbdmax),kbmin(nbdmax),kbmax(nbdmax)
      integer mrb(nbdmax)
      real    crb(3,nbdmax),trb(3,nbdmax),erb(3,nbdmax),erb0(3,nbdmax),rrb(3,nbdmax),arb(3,nbdmax)
      character solid(nbdmax)*80

      common /iobject/ lb,mb,mrb,ibmin,ibmax,jbmin,jbmax,kbmin,kbmax,iprsdv,imbdcmp,nflu,nflp,nfcmax,nintmax
     &                ,idomy,imb3d,nntr2d,nntr3d,lv,mv,ivelrcn,igrdint,itagy,imls,izero
      common /robject/ crb,trb,erb,erb0,rrb,arb,imbovlp,uext1,uext2,uext3
      common /cobject/ solid

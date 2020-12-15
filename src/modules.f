c-------------------------------------------------------- A. Posa - Oct 2011 -----

      module sections
         integer, dimension(:), allocatable :: ksez2,ksez4,ksez6   
      end module sections

      module mediey
         real, dimension(:,:,:), allocatable :: tempo1,
     %uavg,vavg,wavg,vmodavg,vorazavg,voraz,vorravg,vorr,
     %vorzavg,vorz,pavg,ptotavg,tvavg
      end module mediey

      module rmsy
         real, dimension(:,:,:), allocatable :: vrms,urms,wrms,
     %prms,uvs,uws,vws
      end module rmsy  

      module mediex
         real, dimension(:,:,:), allocatable :: tempo2,
     %ruavg,rvavg,rwavg,rvmodavg,rvorazavg,rvoraz,rvorravg,rvorr,
     %rvorzavg,rvorz,rpavg,rptotavg,rtvavg
      end module mediex

      module rmsx
         real, dimension(:,:,:), allocatable :: rvrms,rurms,rwrms,
     %rprms,ruvs,ruws,rvws
      end module rmsx

      module mediez
         real, dimension(:,:,:), allocatable :: ztempo,
     %zuavg,zvavg,zwavg,zvmodavg,zvorazavg,zvorravg,
     %zvorzavg,zpavg,zptotavg,ztvavg
      end module mediez

      module rmsz
         real, dimension(:,:,:), allocatable :: zvrms,zurms,zwrms,
     %zprms,zuvs,zuws,zvws
      end module rmsz

      module vorticity
         real, dimension(:,:,:), allocatable :: voraz1,vorr1,vorz1,
     %voraz2,vorr2,vorz2,voraz14,vorr14,vorz14
      end module vorticity

      module mediecalc
         real, dimension(:,:,:), allocatable :: tempo4,tempo5,tempo14
      end module mediecalc

      module mediecalc1
         real, dimension(:,:,:), allocatable :: vplot,uplot,wplot,
     %pplot,vrmsplot,urmsplot,wrmsplot,prmsplot,uvplot,vwplot,uwplot,
     %vmodplot,tvplot
      end module mediecalc1

      module mediecalc2
         real, dimension(:,:,:), allocatable :: vplot2,uplot2,wplot2,
     %pplot2,vrmsplot2,urmsplot2,wrmsplot2,prmsplot2,uvplot2,vwplot2,
     %uwplot2,vmodplot2,tvplot2
      end module mediecalc2

      module mediecalc3
         real, dimension(:,:,:), allocatable :: vplot14,uplot14,
     %wplot14,pplot14,vrmsplot14,urmsplot14,wrmsplot14,prmsplot14,
     %uvplot14,vwplot14,uwplot14,vmodplot14,tvplot14
      end module mediecalc3

      module vortavg
         real, dimension(:,:,:), allocatable :: voraz1med,vorr1med,
     %vorz1med,voraz2med,vorr2med,vorz2med,voraz14med,vorr14med,
     %vorz14med
      end module vortavg

      module points
         integer, dimension(:), allocatable :: iprb2,jprb2,kprb2
         real, dimension(:), allocatable :: timepoints1
         real, dimension(:,:), allocatable :: avrpoints1,rmspoints1
      end module points

      module spnge
         integer :: vspngx1,vspngx3in,vspngx3out 
         integer :: dspngx1,dspngx3in,dspngx3out 
         integer :: idspngl,ivspngl 
         real, dimension(:,:), allocatable :: PhiX1
         real, dimension(:,:,:,:), allocatable :: X1inf
         real, dimension(:), allocatable :: dfxug,dfxwg,dfxul,dfxwl,dfxu
     &                                      ,dfxw,dfxpg,dfxpl

         real, dimension(:), allocatable :: dfxdg,dfxdl,dfxd

      end module spnge

      module density_bg
         real, dimension(:,:), allocatable :: dens_bg
      end module density_bg

         

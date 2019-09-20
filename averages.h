      integer nsez1,nsez2,nsez3,nsez4,nsez5,nsez6,nsez9,nsez91,
     %nsez14,m1,m2,m3,nsezmax19
      integer irms,nprobes,imedie
      real periodo,periodo2,periodo3,periodo4,periodo5
      parameter (nsez1=3,nsez2=3)
      parameter (nsez3=0,nsez4=0,nsez5=0,nsez6=0,nsez9=0,nsez91=0)
      parameter (nsez14=2)
      parameter (nsezmax19=1)
      parameter (m1=329,m2=329,m3=1026)
      parameter (irms=1,nprobes=18)
      parameter (periodo=1000000000.,periodo2=51.3,periodo3=51.3,
     %periodo4=51.3,periodo5=51.3,imedie=1)
      integer nsez2max,nsez4max,nsez6max,
     %iplant4,iplant5,iplant8,iplant9
      integer isez1(nsez1),isez3(nsez3),jsez5(nsez5),
     %isez9(nsez9),jsez91(nsez91),jsez14(nsez14),ksez2g(nsez2)
      integer isezindx2,isezindx4,isezindx6
      real tempo
      real timep,tempo3
      integer prbindx2,nprbmax2
      real timepoints
      integer jpv(m2),jmv(m2)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      common /cmedie/ tempo,iplant8,iplant9
      common /cmediey/ isez3,nsez4max
      common /cmediex/ jsez5,nsez6max
      common /cmediez/ isez9,jsez91
      common /cvorticity/ isez1,jsez14,jpv,jmv
      common /cplot/ ksez2g,nsez2max
      common /cmediecalc/ timep,tempo3,iplant4,iplant5
      common /csezindices/ isezindx2,isezindx4,isezindx6
      common /cpoints/ prbindx2,nprbmax2,timepoints

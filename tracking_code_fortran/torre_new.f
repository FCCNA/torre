      program torre
*
**********************************************************
*                                                        * 
*    Reconstruction program for the ATLAS                *  
*    RPC test station (torre)                            *
*                                                        *
*    v 5.2                                               *
*    v 5.1    1- 4-00                                    *
*    v 4.1    1- 2-00                                    *
*    v 3.3    5-10-99                                    *
*    v 3.2   14- 8-99                                    *   
*    v 3.1   20- 7-99                                    * 
*    v 3.0    9- 7-99                                    *  
*    v 2.4   30- 6-99                                    *
*    v 2.3   13- 5-99                                    *
*    v 2.2    6- 5-99                                    *
*    v 2.1   20- 4-99                                    *
*    v 2.0   25- 3-99                                    *
*    v 1.1   27-11-98                                    * 
*    v 1.0   24-11-98                                    * 
*                                                        * 
*                                   Authors:             *
*                          M.Alviggi,M.Caprio,G.Carlino  *
*                                                        *
**********************************************************
*                                                          
      
      INTEGER event,flagmc,nevmax,ievent,event_selected
*      REAL xvt,xvb
*     
      PARAMETER (emax = 1000.)
      PARAMETER (emin = 0.3) 
      PARAMETER (pi = 3.1415926)   
      REAL ene 
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
      DATA  eve, endofdata, good_xtrack, good_ytrack 
     >    /.true.,.false.,.false.,.false./
      LOGICAL dcl,rpcl,scil,patl,trackl
      COMMON /logic/ dcl,rpcl,scil,patl,trackl
      DATA dcl,     rpcl,    scil,   patl,   trackl 
     >  /.false., .false., .false., .false., .false./ 
*
      PARAMETER (NPAWC=10000000)
      PARAMETER (LRECL = 1024)
      COMMON /PAWC/ HMEMOR(NPAWC)
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      CHARACTER*10 TOPDS
 1    CHARACTER*80 CURDR
      DATA TOPDS /'//USRNTP'/
 
      EXTERNAL rintegral
*
* --------------------------------------------------------
*
* ** DataCards & Hbook initialization
*
      call init
*            
* ** Geometry of the test station
*
      call geometry
*
* ** Data (flagmc =1) or Monte Carlo (flagmc=-1). Two different paths
*
      if (flagmc.eq.1) then
*
* ** Drift velocity measurement 
*
        if (ivdm(1).eq.1) then
          CALL VDRMEAS
          goto 800
        endif
*
* ** Read the data file
*
        lf=index(ifile,' ')
        OPEN (unit=11,file=ifile(1:lf),status='unknown')     
*
* ** Read the global header 
*
        CALL READHEADER
*   
* ** Loop on the events
*
        event_selected = 0
        ievent = 0       
 100    if (debug(1)) print *,'   '
        ievent = ievent + 1
        if (ievent.gt.nevmax) goto 999    
        CALL READDATA
        if (.not.eve) goto 100
        if (endofdata) goto 999
*
* ** Decode the Scintillators pattern
*        
        CALL SCINTPAT
*
* ** RPC        
*
* **** Decode the TDC data 
*
        CALL RPCTDC
*
* **** RPC clusters
*
        CALL RPCCLUS
*
* ** Drift Chambers
*
* **** Decode the TDC data and find the hits position
*
        CALL DRIFTCHTDC 
*
* **** Pattern recognition and fit 
*
        CALL DRIFTCHPAT
*
* ** Matchings
*
        if (trackl) then          
*
* **** Track with scintillators
*
          CALL SCINTMATCH
*
* **** Track with RPCs
*
          CALL RPCMATCH      
*
        endif
*  
        event_selected = event_selected + 1
        if (debug(1)) print *,'selected event ',event_selected
        CALL hfnt(9)
        goto 100
*
 999    continue
        CLOSE (11)
*
* ** Final results
*
        if (debug(1)) then
          print *,'  ' 
          print *,'**** Selected events =',event_selected
        endif
        CALL FINAL        
*
      elseif (flagmc.eq. -1) then
**
* ** Monte Carlo Simulation
*  
        print *,' '
        print *,' ** MC Simulation ** '
        print *,' '
*
        nevsim = 0
        event = 0
        nevsim_trigger = 0   
*
* Cosmic Flux (expressed in Hz) determination
*  
        epsilon = 1.e-5
        summ_m = gauss(rintegral,emin,emax,epsilon)
        call flux_electron(summ_e)
        call surface(itrig_mc,surf)
        fluxmu = 2*pi/3*summ_m*surf
        fluxe  = 2*pi/4.6*summ_e*surf
        print *,"--- flux of muon --- ",fluxmu,
     >          " --- flux of electron --- ",fluxe

 500    continue
          nevsim = nevsim+1
          event = event + 1         
*
          CALL MCSIMULATION(nevsim_trigger,nevsim)
*
* ** Pattern recognition and fit 
**
          CALL DRIFTCHPAT
**
* ** RPC Cluster Simulation
**
*          CALL RPCCLUSMC
**

          CALL hfnt(9)
          if (nevsim.ne.nevmax) goto 500 
*
        print *,'   '
        print *,'Total Number of events   = ',nevsim
        print *,'Total Number of triggers = ',nevsim_trigger
        trate = (fluxmu+fluxe)*float(nevsim_trigger)/float(nevsim)
        print *,'Trigger Rate (Hz)        = ',trate
*
      endif     
*
* ** Chiusura ntupla
*
      call hcdir (curdr,'R')
      lench = lenocc(curdr)
      call hcdir(topds,' ')
*      call hprnt(9)
      call hrout(0,icycle,' ')
      call hcdir (curdr,' ')
      call hrend('usrntp')
*
 800  continue

      call hrput (0,'vd.hist',' ')
*
      END
*
* **********************************************************
*
* **********************************************************
*
      SUBROUTINE INIT
*
*     Inizializzazione datacards, Hbook and t0corrections 
*
* **********************************************************
*
      INTEGER flagmc,nevmax,nbin(2)
      REAL xvar(2)
      COMMON /t0/ t0corr(96)
*
      PARAMETER (NPAWC=10000000)
      PARAMETER (LRECL = 1024)
      COMMON /PAWC/ HMEMOR(NPAWC)
      CHARACTER*2 cnhdcm,cnfitm,cnscim,cnsvmem,cnscamm,cnmaxc      
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      CHARACTER*10 TOPDS
      CHARACTER*80 CURDR
      DATA TOPDS /'//USRNTP'/
*
* --------------------------------------------------
* 
* ** Read datacards
*
      CALL FFINIT(0)

      flagmc  = 1      ! 1 Data -(-1) MonteCarlo
      nevmax = 1000
      vd(1) = 0.04     ! mm/ns x bottom layers
      vd(2) = 0.04     ! mm/ns x top layers 
      vd(3) = 0.04     ! mm/ns y bottom layers
      vd(4) = 0.04     ! mm/ns y top layers
      ivdm(1) = 0         ! drift velocity measurement
      ivdm(2) = 1000      ! # of events for vd measurement 
      ivdm(3) = 11        ! # of vd steps
      ivdm(4) = 10        ! vd step / 10000
      tdcr(1) = 1.04        ! (ns) resolution VME TDC 
      tdcr(2) = 0.025       ! (ns) resolution CAMAC TDC
      ifile = 'ifile.txt'   ! data file
      t0fil = 't0file.txt'  ! T0cal file for the DC
      ntdc = 1              ! Number of TDCs for the RPC 
      chi2cut = 100.      ! refitting chi2cut
      flag_ref = 1        ! refitting flag
      nbin(1) = 45        ! Pattern Recognition - Binning slope 
      nbin(2) = 50        ! Pattern Recognition - Binning intercepta
      itrig_mc = 1        ! Trigger in the MC simulation
      call vzero(isele,3) ! Event selection
      debug(1) = .false.  ! Warnings      
      debug(2) = .false.  ! Global printing       
      debug(3) = .false.  ! RPC 
      debug(4) = .false.  ! DC and Scintillators
      debug(5) = .false.  ! DC pattern recognition 
      debug(6) = .false.  ! DC pattern recognition (verbose)
      debug(7) = .false.  ! Drift velocity measurement
      debug(8) = .false.  ! Simulation
      debug(9) = .false.  ! Raw data 
      debug(10) = .false. ! Geometry
*
      CALL FFKEY('DAMC',flagmc,1,'INTEGER')
      CALL FFKEY('IFIL',ifile,80,'MIXED')
      CALL FFKEY('T0FI',t0fil,80,'MIXED')
      CALL FFKEY('NEVM',nevmax,1,'INTEGER')
      CALL FFKEY('VDRI',vd,4,'REAL')
      CALL FFKEY('VDRM',ivdm,4,'INTEGER')
      CALL FFKEY('TDCR',tdcr,2,'REAL')     
      CALL FFKEY('NTDC',ntdc,1,'INTEGER')
      CALL FFKEY('CHIC',chi2cut,1,'REAL')
      CALL FFKEY('HBIN',nbin,2,'INTEGER')
      CALL FFKEY('DEBU',debug,10,'LOGICAL')  
      CALL FFKEY('REFI',flag_ref,1,'INTEGER')
      CALL FFKEY('SELE',isele,3,'INTEGER')
      CALL FFKEY('XVAR',xvar,2,'REAL')
      CALL FFKEY('TRIG',itrig_mc,1,'INTEGER')
      CALL FFGO       
      iflag = flagmc
      ns=nbin(1)
      nx=nbin(2)
      xvt=xvar(1)
      xvb=xvar(2)
      xv1t=xvar(1)
      xv1b=xvar(2)
 
      if (isele(1).eq.1) then
        do i=1,10
          debug(i) = .false.
        enddo
      endif  
*
* ** Hbook
*
      call hlimit(npawc)      
*
* ** Ntuple
*
      CALL HROPEN(9,'USRNTP','torre.ntp','N',1024,ISTAT)
      IF (ISTAT.NE.0) PRINT *, 'PROBLEMS WITH HROPEN!!!!!'
      CALL HCDIR (CURDR,'R')
      LENCH = LENOCC(CURDR)
      CALL HCDIR(TOPDS,' ')
      CALL HBNT (9,'ntupla torre',' ')
      CALL HBNAME (9,'GLOBAL',event,
     + 'nevnt,iblock,iflag,tempdc,temprpc,gasp,hv,itrig')
      write(cnscim,'(i2)') nscim        
      CALL HBNAME (9,'SCINT',nscint,
     + 'nscint[0,'//cnscim//']'//
     + ',iscint(nscint),nscintl,iscitime(nscint),iscil(nscint)') 
      write(cnsvmem,'(i2)') nsvmem
      write(cnscamm,'(i2)') nscamm
      CALL HBNAME (9,'RPC',nhrpc,
     + 'nhrpc'//
     + ',nsxvme[0,'//cnsvmem//'],isxvme(nsxvme)'//
     + ',itdcdataxvme(nsxvme),timexvme(nsxvme)'//
     + ',nsyvme[0,'//cnsvmem//'],isyvme(nsyvme)'//
     + ',itdcdatayvme(nsyvme),timeyvme(nsyvme)'//
     + ',nsxcam[0,'//cnscamm//'],isxcam(nsxcam)'//
     + ',itdcdataxcam(nsxcam),timexcam(nsxcam)'//
     + ',nsycam[0,'//cnscamm//'],isycam(nsycam)'//
     + ',itdcdataycam(nsycam),timeycam(nsycam)')
      write(cnmaxc,'(i2)') nmaxc
      CALL HBNAME (9,'RPCCLU',ncxvme,
     + ' ncxvme[0,'//cnmaxc//'],ifscxvme(ncxvme),mulcxvme(ncxvme)'//
     + ',ncyvme[0,'//cnmaxc//'],ifscyvme(ncyvme),mulcyvme(ncyvme)'//
     + ',ncxcam[0,'//cnmaxc//'],ifscxcam(ncxcam),mulcxcam(ncxcam)'//
     + ',ncycam[0,'//cnmaxc//'],ifscycam(ncycam),mulcycam(ncycam)')
      write(cnhdcm,'(i2)') nhdcm      
      CALL HBNAME (9,'DRIFTCH',nhdc,
     + 'nhdc[0,'//cnhdcm//'],vdxbottom,vdxtop,vdybottom,vdytop'//
     + ',nxlayer,nylayer,nxhdc(4),nyhdc(4)'//
     + ',ichan(nhdc),layer(nhdc),hitwire(nhdc),zhitwire(nhdc)'//
     + ',itdcdatadc(nhdc),dtime0(nhdc),dtime(nhdc),dlength(nhdc)'//
     + ',xvt,xvb')
      write(cnfitm,'(i2)') nfitm
      CALL HBNAME (9,'XTRACK',ixtrack,
     + 'ixtrack,thetax,ax,bx,chi2x,nxfit[0,'//cnfitm//']'//
     + ',layerx(nxfit),xres(nxfit),xpos(nxfit)')
      CALL HBNAME (9,'YTRACK',iytrack,
     + 'iytrack,thetay,ay,by,chi2y,nyfit[0,'//cnfitm//']'//
     + ',layery(nyfit),yres(nyfit),ypos(nyfit)')
      CALL HBNAME (9,'TRACK',itrack,
     + 'itrack,theta3D, phi3D')
      CALL HBNAME (9,'MCSIM',itrig_mc,
     + 'itrig_mc,fluxmu,fluxe,trate,theta_mc,phi_mc') 
      CALL HBNAME (9,'TRACKMC',nhdc_mc,
     + 'nhdc_mc [0,'//cnhdcm//']'//
     + ',ichanmc(nhdc_mc),layermc(nhdc_mc),posmc(nhdc_mc)'//
     + ',hitwiremc(nhdc_mc),zhitwiremc(nhdc_mc),dlengthmc(nhdc_mc)'//
     + ',thetaxmc,axmc,bxmc,thetaymc,aymc,bymc,gauss(nhdc_mc)'//
     + ',iscintmc(2),itrf,multipmc,nxstripmc(10),nxhdcmc(4)'//
     + ',nyhdcmc(4)')
      CALL HBNAME (9,'RPCM',nsxexp,
     + 'nsxexp[0,10],isxexp(nsxexp),nsyexp[0,10],isyexp(nsyexp)')
      CALL HBNAME (9,'SCINTM',nscintexp,
     + 'nscintexp[0,2],iscintexp(nscintexp)') 
*
* ** Histos for vd measurment
*
       if (ivdm(1).eq.1) then
         vd0xb = vd(1) 
         vd0xt = vd(2) 
         vd0yb = vd(3) 
         vd0yt = vd(4) 
         vdstep = ivdm(4)/10000.
         nstep = ivdm(3)
         vdxbmin = vd0xb - vdstep*(nstep/2.)
         vdxtmin = vd0xt - vdstep*(nstep/2.)
         vdybmin = vd0yb - vdstep*(nstep/2.)
         vdytmin = vd0yt - vdstep*(nstep/2.)
         do iv=0, nstep-1
           vdxb = vdxbmin + vdstep*iv
           vdxt = vdxtmin + vdstep*iv
           vdyb = vdybmin + vdstep*iv
           vdyt = vdytmin + vdstep*iv idbx = (1000+vdxb*10000)*10+1
          idbx = (1000+vdxb*10000)*10+1
          idby = (1000+vdyb*10000)*10+3
          idty = (1000+vdyt*10000)*10+4 
           call hbook1(idbx,'X bottom layers',40,-5.,5.,0.)
           call hbook1(idtx,'X top layers',40,-5.,5.,0.)
           call hbook1(idby,'Y bottom layers',40,-5.,5.,0.)
           call hbook1(idty,'Y top layers',40,-5.,5.,0.)
         enddo
       endif
*
* ** Load t0 file
*
      OPEN (unit=12,file='t0ch.txt',status='unknown') 
      READ (12,*,err=999) t0corr
*
 999  continue
*      do i=1,96 
*        print *,i,' t0corr =',t0corr(i)
*      enddo
      close (12)
*      
      RETURN 
      END
*
* **************************************************************
*
      SUBROUTINE geometry
*
*     Description of the geometry of the tower
*
*|                                              
*|              5          6         7     |     8                
*|     |--|----[]----|----[]----|----[]----|----[]----|-|  -----   3245  mm  
*|     |  |                     |          |                  |     (stsbz)    
*|      170 mm                     250 mm                     |         
*|      (sdy)                                                 | (stdz) 
*|                                                            |  130 mm
*|                                                            |         
*|     |-|-------------------------------------------------| -------------- 
*|     |---------------------------------------------------|             |
*|   y |   84                                        95    |  125 mm  -->|  
*|     |---------------------------------------------------|             |  
*|   y |  72                                        83     |             |  
*|     |---------------------------------------------------|  -----      | 
*|   x |  |60 |   |                                  |71 | |  29 mm (zc) |     
*|     |---------------------------------------------------|  -----      | 
*|   x | |48 |49 |   |                             |59 |   |             | 
*|     |---------------------------------------------------|             |
*|     |-|-------------------------------------------------|  -------------    
*|                                                            /          |
*|                                                            /          |
*|      RPC                                                   /          |
*|     layer  Y:   --------------------------               ----- +6     | 
*|            X:   - - - -          - - - - -                     +6     |
*|            Y:   --------------------------                     +6     |
*|            X:   - - - -          - - - - -               -----  ? mm  | 
*|                                                            |          |
*|     |---  ?  ---|-------- 496 mm  --------|                |          |
*|                                                            /          | 
*|                                                            |          |
*|     87.75 (cdy)                                               /       |
*|     |-|-------------------------------------------------|             |
*|     |---------------------------------------------------|             |
*|   y |    36                                        47   |             |
*|     |---------------------------------------------------|             |
*|   y |   24                                        35    |     (zz0)-> |
*|     |---------------------------------------------------|    2850 mm  |
*|   x |  |12 |   |                                  |23 | |             |
*|     |---------------------------------------------------|             |
*|   x | | 0 | 1 |   |                             | 11|   |             |
*|     |---------------------------------------------------|             |
*|     |-|-------------------------------------------------| -------------   
*|    88.75    |   |                                          |             
*|   (cdx)      92 mm                                         |    (sbdz)   
*|               (xc)                   125 mm                |    140 mm
*|              1          2         3|    |     4            |    
*|     |--|----[]----|----[]----|----[]----|----[]----|-|   -----       
*|     |  |                     |          |                  |     
*|     170 mm                      250 mm                     |     
*|     | (sdy)                      (sy)                      |  ?  (sb0z)
*| d0y |                                                      |
*|<--->|                                                      | 
*|---------------------------------------------------------------------> y
*0
*
*
*
*                                                            z |
*                                                              |   
* Local Reference Frame: (x,y,z) = (0,0,0) corrisponding       | 
* to the  origin of reference frame                            |_____ y
*                                                              / 
*                                                             /
*                                                            / x
*
* **************************************************************************

      INCLUDE "dc.inc"
*
* ----------------------------------------------------------------
*
      if (debug(10)) then
        print *,' '
        print *,' GEOMETRY '
        print *,'  '
      endif   
*
* ** Scintillator geometry
*
      CALL SCINTG
*
* ** Drift chambers geometry 
*
      CALL DRIFTCHG 
*
* ** RPC geometry
*
      CALL RPCG
*
      if (debug(10)) then
        print *,' '
        print *,' END GEOMETRY '
        print *,'  '
      endif  
*
      RETURN
      END
*
* *****************************************************************
*
      SUBROUTINE SCINTG
*
* Scintillators position
* scint 1 (73.8+125,573.8,-157)
* scint 2 (73.8+250+125,573.8,-157)
* scint 3 (73.8+250*2+125,573.8,-157) 
* scint 4 (73.8+250*3+125,573.8,-157)
* scint 5 (73.8+125,573.8,794)
* scint 6 (73.8+250+125,573.8,794)
* scint 7 (73.8+250*2+125,573.9,794) 
* scint 8 (73.8+250*3+125,573.8,794)
*
***************************************************************
*
      INCLUDE "dc.inc"

      PARAMETER (nscintil=8)
      COMMON /scintgeo/  yscint(nscintil),zscint(nscintil)
      COMMON /CONST/ sdy,sy,stsbz,sb0z,d0x,d0y,cdx,cdy
*     sdy = distanza lungo y tra scint. e DC
*     sy = dimensioni y scintillatore
*     stsbz = distanza scint. top scint. bottom
*     st0z = distanza scint top dall'origine lungo z
*     sb0z = distanza scint bottom dall'origine lungo z
*     d0x = distanza DC dall'origine lungo x
*     d0y = distanza DC dall'origine lungo y
*     cdx = distanza ibeam - bordo della camera lungo x
*     cdy = distanza ibeam - bordo della camera lungo y
*
      DATA sdy, sy , stsbz, sb0z, d0x, d0y, cdx, cdy
     >     /170., 250., 3245., 0., 0., 0., 88.75, 87.75/
*     >     /73.8, 250., 951., -157., 0., 0., 0., 0./ ! vecchia geometria
*
*-------------------------------------------------------------
*
      do j= 1,4
         yscint(j)=d0y+sdy+(j-1)*sy+sy/2
         yscint (j+4)=d0y+sdy+(j-1)*sy+sy/2
         zscint (j)=sb0z
         zscint (j+4)=stsbz + sb0z
      enddo   
*
      if (debug(10)) then
         print *,'   '
         print *,' Scintillators'
         print *,'  '
         do i=1,8
            print *,'scint,pos,zpos',i,yscint(i),zscint(i)
         enddo
      endif  
*       
      RETURN
      END
*
****************************************************************
*
*
      SUBROUTINE Surface(itrig_mc,surf)
*
*     According to the trigger (itrig_mc) used in the 
*     MC Simulation determine the scintillator surface
*
*     itrig_mc = 1 -  General Trigger Oriz. - x = 100 cm , y = 100 cm
*     itrig_mc = 2 -  Central Trigger Oriz. - x = 100 cm , y = 50 cm
*     itrig_mc = 3 -  General Trigger Vert. - x = 100 cm , y = 100 cm
*     itrig_mc = 4 -  Central Trigger Vert. - x = 100 cm , y = 50 cm
*
* *****************************************************************
*
      xsurf = 100. 
      if (itrig_mc .eq.1.or.itrig_mc.eq.3) then
        ysurf = 100.
      elseif (itrig_mc.eq.2.or.itrig_mc.eq.4) then  
        ysurf = 50.
      endif
*
      surf = xsurf*ysurf  
*
      RETURN
      END 
*
********************************************************************
*
      SUBROUTINE RPCG
*
*      RPC geometry
*
*   Strip position:
*   X layer: 
*            strip 0  (342+15,340)   
*            strip 15 (342+30*15+15,340)
*   Y layer: 
*            strip 0  (352+15,346)
*            strip 15 (352+30*15+15,346)
*
****************************************************************
*
      INCLUDE "dc.inc"

      PARAMETER (stripw = 30.)            
      PARAMETER (offx = 342.)
      PARAMETER (offy = 352.)
      PARAMETER (zxstrip = 340.)
      PARAMETER (zystrip = 346.) 
      PARAMETER (nxstrip = 16)
      PARAMETER (nystrip = 16)
* 
      COMMON /rpcgeo/ stx(0:nxstrip-1),stzx(0:nxstrip-1),
     >                 sty(0:nystrip-1),stzy(0:nystrip-1)                  
*
* --------------------------------------------------------------
*
      do is=0,nxstrip-1
        stx(is)=offx+is*stripw+stripw/2.
        stzx(is) = zxstrip
      enddo
*
      do is=0,nystrip-1
        sty(is)=offy+is*stripw+stripw/2.
        stzy(is) = zystrip
      enddo
*
      if (debug(10)) then 
        print *,'   '
        print *,' RPC strips'
        print *,'  '
        do i=0,nxstrip-1 
          print *,'strip ',i,' =',stx(i),stzx(i)
          print *,'strip ',i,' =',sty(i),stzy(i)          
        enddo
      endif  
*      
      RETURN
      END
*
****************************************************************
*
       SUBROUTINE DRIFTCHG                         
*
*      Drift chambers geometry
*
*
* Wires position: 
* cell 0  (92/2,-,29/2),          cell 11 (92/2+92*11,-,29/2)
* cell 12 (92*1,-,29+29/2),       cell 23 (92*12,-,29+29/2)
* cell 24 (-,92/2,29*2+29/2)      cell 35 (-,92/2+92*11,29*2+29/2)
* cell 36 (-,92*1 ,29*3+29/2)     cell 47 (-,92*12,29*3+29/2) 
* cell 48 (92/2*1,-,521+29/2),    cell 59 (92/2*12,-.521+29/2)
* cell 60 (92*1,-,521+29+29/2),   cell 71 (92*12,-,521+29+29/2)
* cell 72 (-,92/2,521+29*2+29/2), cell 83 (-,92/2+92*11,521+29*2+29/2)
* cell 84 (-,92*1,521+29*3+29/2), cell 95 (-,92*12,521+29*3+29/2)
*
* *******************************************************************          
*
       INCLUDE "dc.inc"

      COMMON /driftchgeo/  wire(0:95),zwire(0:95)
      COMMON /CONST/ sdy,sy,stsbz,sb0z,d0x,d0y,cdx,cdy
      PARAMETER (xc=92.) ! dimensioni cella lungo x = lungo y
      PARAMETER (zc=29.) ! dimensioni cella lungo z
      PARAMETER (zz0=2850.) !distanza fra le due camere a drift(bot-bot)
      PARAMETER (wibeam = 1.5)
*      PARAMETER (cdx=88.75)     !distanza bordo DC-bordo cella 0 
*      PARAMETER (cdy=87.75)     !distanza bordo DC-bordo cella 24 
*      PARAMETER (xvt = 0.)      !spostamento variabile DC-d0x top
*      PARAMETER (xvb = 0.)      !spostamento variabile DC-d0x bottom      
      PARAMETER (sbdz = 140.) !distanza DC-scint(bottom) lungo z
      PARAMETER (stdz = 130.) !distanza DC-scint (top) lungo z

*     
* --------------------------------------------------------------------
*
* Load the wires position
*
*** X-chamber ***
      
      do ich = 0,11
           
         wire(ich) = (xvb+cdx+d0x)+ xc/2.+wibeam/2. + (xc+wibeam)*ich
         wire(ich+12) = wire(ich)+ xc/2.+wibeam/2.
         wire(ich+48) = wire(ich)- xvb + xvt
         wire(ich+60) = wire(ich+48) + xc/2.+wibeam/2.
      enddo
           
      
*** Y-chamber ***

      do i = 0,11
         ich = i + 24
         wire(ich) = (cdy+d0y)+ xc/2.+wibeam/2. + (xc+wibeam)*i
         wire(ich+12) =wire(ich)+ xc/2.+ wibeam/2.
         wire(ich+48) = wire(ich)
         wire(ich+60) = wire(ich+48) + xc/2.+wibeam/2.
      enddo

*

      do i=1,4
         zcc = zc/2. + zc*(i-1)  
         do j=0,11
            ich=j + 12*(i-1)       
            zwire(ich) = zcc + 5. + sbdz + sb0z
***   ho 5. mm dal bordo della dc al bordo della cella lungo z ***
         enddo 
      enddo 
      do i=0,47
         ich2 = 48+i
         zwire(ich2) = zwire(i) + zz0
      enddo 
*
      if (debug(10)) then
         print *,'    '
         print *,' DC geometry'
         print *,'   '
         do i=0,95 
            print *,i,' wire =',wire(i),zwire(i)
         enddo
      endif  
*       
      RETURN 
      END   
*
* ***********************************************************
*
      SUBROUTINE VDRMEAS
*
*    Automatic measurement of the drift velocity
*    (The routine is executed before the main routines)
*
* ************************************************************
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
*      PARAMETER (np=6)
*      REAL par(np),step(np),pmin(np),pmax(np),sigpar(np)
*      REAL parx(3,4),pary(3,4),chi2x(4),chi2y(4)
*
* ---------------------------------------------------------      
*
      if (debug(7)) then   
        print *, '        '
        print *,'VDRMEAS'
        print *,'      '
      endif  
*
      vd0xb = vd(1) 
      vd0xt = vd(2)
      vd0yb = vd(3) 
      vd0yt = vd(4)
*
* ** Drift velocity variation
*
      vdstep = ivdm(4)/10000.
      nstep = ivdm(3)
      if (debug(7)) then   
        write (*,444) vd0xb,vd0xt
        write (*,445) vd0yb,vd0yt
 444    format (1x,'vd0 X layers - bottom =',f6.4,1x,'top =',f6.4)
 445    format (1x,'vd0 Y layers - bottom =',f6.4,1x,'top =',f6.4)
        write (*,443) nstep , vdstep
 443    format (1x,'# of steps =',i3,1x,'- vdstep =',f6.4)
      endif        
      vdxbmin = vd0xb - vdstep*(nstep/2.)
      vdxtmin = vd0xt - vdstep*(nstep/2.)       
      vdybmin = vd0yb - vdstep*(nstep/2.)
      vdytmin = vd0yt - vdstep*(nstep/2.)      
      do iv=0,nstep-1
        vdxb = vdxbmin + vdstep*iv
        vdxt = vdxtmin + vdstep*iv
        vdyb = vdybmin + vdstep*iv
        vdyt = vdytmin + vdstep*iv	
	if (debug(7)) then
          print *,'      '
          print *,'New drift velocities (mm/ns)'
	  write (*,446) vdxb , vdxt
	  write (*,447) vdyb , vdyt  
 446	  format (1x,'X layers - bottom =',f6.4,1x,'- top =',f6.4)
 447	  format (1x,'Y layers - bottom =',f6.4,1x,'- top =',f6.4)  
        endif   
        idxb = (1000+vdxb*10000)*10
        idxt = (1000+vdxt*10000)*10
        idyb = (1000+vdyb*10000)*10
        idyt = (1000+vdyt*10000)*10 	   
        vd(1) = vdxb
        vd(2) = vdxt
	vd(3) = vdyb
	vd(4) = vdyt
*
* ** Read the data file
*  
        lf=index(ifile,' ')
        OPEN (unit=11,file=ifile(1:lf),status='unknown')     

*
* ** Read the global header 
*
        CALL READHEADER
*   
* ** Loop on the events
*
        ivdmev = 0       
 777    continue
        ivdmev = ivdmev + 1
        if (ivdmev.gt.ivdm(2)) goto 888    
        CALL READDATA
        if (.not.eve) then
          goto 777
        endif 
        CALL DRIFTCHTDC 
        CALL DRIFTCHPAT
*     
        if (nfitx_good.eq.4.and.chi2redx_good.lt.15) then
          idhxb = idxb+1
	  idhxt = idxt+2
	  call hf1(idhxb,xres_good(1),1.)
	  call hf1(idhxb,xres_good(2),1.)
	  call hf1(idhxt,xres_good(3),1.)
	  call hf1(idhxt,xres_good(4),1.)
        endif  
        if (nfity_good.eq.4.and.chi2redy_good.lt.15) then
          idhyb = idyb+3
	  idhyt = idyt+4
	  call hf1(idhyb,yres_good(1),1.)
	  call hf1(idhyb,yres_good(2),1.)
	  call hf1(idhyt,yres_good(3),1.)
	  call hf1(idhyt,yres_good(4),1.)	
        endif  
*
        goto 777
*
 888    continue
*
        CLOSE (11)
*
* ** Histograms fit
*
*        par(1) = 100.
*        par(2) = -1.5
*        par(3) = 1.
*        par(4) = 100.
*        par(5) = 1.5
*        par(6) = 1.
*        pmin(1) = -100.
*        pmin(2) = -2.
*        pmin(3) = -4.
*        pmin(4) = -100.
*        pmin(5) = -2.
*        pmin(6) = -4.
*        pmax(1) = 100.
*        pmax(2) = 2.
*        pmax(3) = 4.
*        pmax(4) = 100.
*        pmax(5) = 2.
*        pmax(6) = 4.
*        step(1) = 1.
*        step(2) = .01
*        step(3) = .01
*        step(4) = 1.
*        step(5) = .01
*        step(6) = .01     

*        call hfithn(71,'G','Q',3,par,step,pmin,pmax,sigpar,chi2)   
*      call hfithn(71,'G+G','B',np,par,step,pmin,pmax,sigpar,chi2)
*
*        print *,'par =',par
*
*          do j=1,2
*        j=1
*        idx = 70+j
*        idy = 80+j 
*         do ip =1,3
*           parx(ip,j) = par(ip)
*         enddo
*         chi2x(j) = chi2
*         call hfithn(idy,'G','V',3,par,step,pmin,pmax,sigpar,chi2) 
*         do ip =1,3
*           pary(ip,j) = par(ip)
*          enddo            
*         chi2y(j) = chi2  
*          enddo  
*          do j=1,2
*            j=1
*            do i=1,3
*              print *,i,j,' parx =',parx(i,j)
*            enddo
*            print *,'chi2x =',chi2x(j)
*          enddo 
*          do j=1,4
*            do i=1,3
*              print *,i,j,'pary =',pary(i,j)
*            enddo
*            print *,'chi2y =',chi2y(j)
*          enddo  
*        endif  
      
*        call hreset(0,' ')
*      enddo

      enddo

      RETURN
      END
*
* ***********************************************************
*
      SUBROUTINE READHEADER
*
*    Read the global header 
*    Run, RPC serial number, gas mixture .......     
*
* ************************************************************
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
      CHARACTER*80 date
*
* --------------------------------------------------------------------
*
      if (debug(2)) then
        print *,'  '
        print *,' READHEADER'
        print *,'   '
      endif  
*
      READ (11,*,end=999) date
      READ (11,*,end=999) jrun 
      if (debug(9)) print *,'jrun =',jrun
      itwobit= jbyt(jrun,31,2)
      if (itwobit.eq.3) then
        irun = jbyt(jrun,1,30)
      else
        if (debug(1)) print *,' - Error in global header - run'
      endif
      READ (11,*,end=999) jrpcserial
      if (debug(9)) print *,'jrpcserial =',jrpcserial
      itwobit= jbyt(jrpcserial,31,2)
      if (itwobit.eq.3) then
        irpcserial = jbyt(jrpcserial,1,30)
      else
        if (debug(1)) print *,' - Error in global header - serial'
      endif
      READ (11,*,end=999) jgasm
      if (debug(9)) print *,'jgasm =',jgasm
      itwobit= jbyt(jgasm,31,2)
      if (itwobit.eq.3) then
        igasm = jbyt(jgasm,1,30)
      else 
        if (debug(1)) print *,' - Error in global header - gas mixture'
      endif
      READ (11,*,end=999) jtrig
      if (debug(9)) print *,'jtrig =',jtrig
      itwobit= jbyt(jtrig,31,2)
      if (itwobit.eq.3) then
        itrig = jbyt(jtrig,1,30)
      else 
        if (debug(1)) print *,' - Error in global header - trigger'
      endif
      READ (11,*,end=999) jrescam
      if (debug(9)) print *,'jrescam =',jrescam
      itwobit= jbyt(jrescam,31,2)
      if (itwobit.eq.3) then
        irescam = jbyt(jrescam,1,30)
      else 
         if (debug(1)) print *,
     >        ' -Error in global header - CAMAC TDC res'
      endif

      if (debug(2)) then
        print *,'Run number =',irun
        print *,'RPC serial Number =',irpcserial
        print *,'gas mixture =',igasm 
        print *,'Trigger =',itrig
        print *,'CAMAC TDC resolution =',irescam
      endif
*
 999  continue
*
      RETURN
      END
*
* ************************************************************
*
      SUBROUTINE READDATA
*
*     Read Data file
*
* *************************************************************
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      INTEGER worddc,wordrpc,word,word1,word2,word3,wordscitime
      COMMON /datadc/ dclen,tdcdc
      COMMON /datasci/ scipat,tdcsci,patlen
      COMMON /datarpc/ rpclen,tdcrpc,nvme,ncam
      PARAMETER (nhrpcm = 80)
      PARAMETER (nhdctot = 200)
      INTEGER scipat,rpclen,dclen,tdcdc(nhdctot),event,tdcrpc(nhrpcm),
     >     tdcsci(8)
      LOGICAL dcl,rpcl,scil,patl,trackl
      COMMON /logic/ dcl,rpcl,scil,patl,trackl
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
*            
* -------------------------------------------------------------
*    
      if (debug(2)) then
        print *,'  '
        print *,' READDATA'
        print *,'   '     
      endif  
*
      dcl =.false.
      rpcl =.false.
      scil =.false.
      patl =.false.
      eve = .true.
*
      READ (11,*,end=999,err=888) word
      if (debug(9)) print *,'word =',word
      if (word.eq.-16757216) then
	print *,'End Of Data' 
        endofdata = .true.
        RETURN
      endif
      
      itwobit= jbyt(word,31,2)
      if (itwobit.eq.3) then
*
* ** Read the block header (temperature, gas pressure, HV)
*
        igasp = jbyt(word,1,30)
        gasp = float(igasp)/10.
        READ (11,*,end=999,err=888) word1
        if (debug(9)) print *,'word1 =',word1
        itwobit= jbyt(word1,31,2)
        if (itwobit.eq.3) itemp = jbyt(word1,1,30)
        tempdc = float(itemp)/10. 
        READ (11,*,end=999,err=888) word2
        if (debug(9)) print *,'word2 =',word2
        itwobit= jbyt(word2,31,2)
        if (itwobit.eq.3) itemp = jbyt(word2,1,30)
        temprpc = float(itemp)/10.        
        READ (11,*,end=999,err=888) word3
        if (debug(9)) print *,'word3 =',word3
        itwobit= jbyt(word3,31,2)
        if (itwobit.eq.3) ihv = jbyt(word3,1,29)
        hv = float(ihv)/1000. 
*
        if (debug(2)) then  
          print *,'Temperature DC, RPC=',tempdc, temprpc
          print *,'Gas pressure =',igasp
          print *,'HV =',hv
        endif   
        eve = .false.
        RETURN        
*
      else if (itwobit.eq. 2) then 
*
* ** Read the data words 
*       
        event = jbyt(word,1,16) 
        iblock = jbyt(word,17,14)
*
        if (isele(1) .eq.1) then
          if (isele(2).eq.event.and.isele(3).eq.iblock) then
            do j=1,9
              debug(j) = .true.
            enddo
          endif
          if ((isele(2)+1).eq.event.and.isele(3).eq.iblock) then
            do j=1,9
              debug(j) = .false.
            enddo
          endif
        endif  

*          
        if (debug(1)) print *,'Event =',event,' Block =',iblock
*
        READ (11,*,end=999,err=888) word2
        if (debug(9)) print *,'word2 =',word2 
        itwobit2= jbyt(word2,31,2)
        if (itwobit2.eq. 2) then
          rpclen = jbyt(word2,9,8)
          dclen = jbyt(word2,17,8)
          scipat = jbyt(word2,1,8)
          patlen = jbyt(word2,25,3)
          if (debug(2)) then
            print *,'dclen =',dclen
            print *,'rpclen =',rpclen
            print *,'scipat =',scipat
            print *,'patlen =',patlen
          endif  
          if (dclen.gt.0)  dcl =.true.
          if (rpclen.gt.0) rpcl = .true.
          if (scipat.gt.0) scil = .true.
          if (patlen.gt.0) patl = .true.
        endif   
      else
        if (debug(1)) print *,' - Word error - Skip event'
        eve = .false.
        RETURN
      endif    
*
* ** DC
*
      do jj=1,dclen
        READ (11,*,end=999,err=888) worddc 
        itwobit=jbyt(worddc,29,4)
        if (debug(9)) print *,jj,' worddc =',worddc,itwobit  
        if (itwobit .eq.0) then 
          tdcdc(jj) = jbyt(worddc,1,24)
        else
          if (debug(1)) print *,' - Error in drift chamber word'          
        endif
      enddo  
*
* ** RPC
*
      nvme = 0 
      ncam = 0
*      rpclen = rpclen + 32
*
      do jj=1,rpclen
        READ (11,*,end=999,err=888) wordrpc
        if (debug(9)) print *,jj,' wordrpc =',wordrpc
        if (ntdc.eq.1) then
           itwobit=jbyt(wordrpc,29,4)
           if (itwobit .eq.4) then
              tdcrpc(jj) = jbyt(wordrpc,1,24) 
           else
              if (debug(1)) print *,' - Error in rpc word. Skip event'
           endif
        else if (ntdc.eq.2) then 
           ithreebit=jbyt(wordrpc,29,4)
           if (debug(1)) print *,'4bitsrpc ',ithreebit
           if (ithreebit .eq.4) then
              nvme = nvme + 1
              tdcrpc(jj) = jbyt(wordrpc,1,24) 
           else if (ithreebit .eq.6) then
              tdcrpc(jj) = jbyt(wordrpc,1,24)
              ncam = ncam + 1               
           else
              if (debug(1)) print *,' - Error in rpc word'
              eve = .false.
              do jjj = jj+1, rpclen
                 READ (11,*,end=999,err=888) wordrpc     
              enddo       
              RETURN  
           endif
        endif
      enddo   
            

*
*** Scintillator Time
*
      do ii=1,patlen
         READ (11,*,end=999,err=888) wordscitime
         if (debug(9)) print *,' wordscintime =',wordscitime  
         ifourbit=jbyt(wordscitime,29,4)
*         if (debug(1)) print *,'4bits ',ifourbit
         if (ifourbit.eq.7) then
            tdcsci(ii) = jbyt(wordscitime,1,19)
         else
            if (debug(1)) print *,' - Error in scitime word' 
         endif
      enddo
      
*
      

      if (dclen.eq.0) then
        if (debug(1)) print *,' - No DC Hits'  
      endif  
      if (rpclen.eq.0) then        
        if (debug(1)) print *,' - No RPC Hits'
      endif
      if (patlen.eq.0) then
        if (debug(1)) print *,' - No fired scintillators'  
      endif 
*
      if (dclen.gt.nhdcm) then
        if (debug(1)) print *,' - Too many DC hits - Skip event'
        eve = .false.
        RETURN
      endif
*
      goto 999
*
 888  if (debug(1)) print *,' - Error reading data'
 999  continue
*
      RETURN
      END       
*
* ************************************************************
*
      SUBROUTINE SCINTPAT
*
*     Decode the scintillators pattern
*
*************************************************************
*
      INTEGER scipat,isci,tdcsci(8),iscil,iscitime
      COMMON /datasci/ scipat,tdcsci,patlen
*     
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
* ------------------------------------------------------
*
      if (debug(4)) then
        print *,'  '
        print *,' SCINTPAT '
        print *,'   '
      endif  
*      
      call vzero(iscint,8)
      call vzero(iscil,8)
      call vzero(iscitime,8)
      nscint = 0
      do i=1,8
        isci=jbyt(scipat,i,1)
        if (isci.eq.1) then 
          nscint = nscint + 1
          iscint(nscint) = i
        endif   
      enddo      
      do i=1,patlen
*         print *,'tdcsci ',tdcsci(i)
         iscil(i) = jbyt(tdcsci(i),17,3)
         iscitime(i) = jbyt(tdcsci(i),1,16)
*         print *,i,'scint ',iscil(i),'time ',iscitime(i)
      enddo
      nscintl = patlen
      if (debug(4)) then
        print *,'scipat =',scipat
        print *,'nscint = ',nscint
        do i=1,nscint  
          print *,'scintillator =',iscint(i)
        enddo        
      endif  
*
      RETURN 
      END
*
* ********************************************************
*
      SUBROUTINE RPCTDC
*
* Decode the RPC TDC word.    
*
***********************************************************
*
      PARAMETER (nhrpcm = 80)
      COMMON /datarpc/ rpclen,tdcrpc,nvme,ncam
      INTEGER tdcrpc(nhrpcm),rpclen
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
* ------------------------------------------------------
*
      if (debug(3)) then
        print *,'   '
        print *,' RPCTDC'
        print *,'  ' 
      endif
*
      nsxvme = 0
      nsyvme = 0
      nsxcam = 0
      nsycam = 0
      call vzero(isxvme,nsvmem)
      call vzero(isyvme,nsvmem)
      call vzero(isxcam,nscamm)
      call vzero(isycam,nscamm)
*
      nhrpc = rpclen
      if (debug(3)) print *,'Total Number of strips =',nhrpc
      if (nhrpc.eq.0) RETURN
*
      do is=1,nhrpc  
        istrip = jbyt(tdcrpc(is),17,7)  
        idat = jbyt(tdcrpc(is),1,16)
        if (is.le.nvme) then
          if (istrip.le.15) then
            nsxvme = nsxvme + 1 
            isxvme(nsxvme)=istrip
            itdcdataxvme(nsxvme) = idat      
            timexvme(nsxvme) = itdcdataxvme(nsxvme)*tdcr(1)
          else if (istrip.gt.15) then
            nsyvme = nsyvme + 1 
            isyvme(nsyvme)=istrip
            itdcdatayvme(nsyvme) = idat      
            timeyvme(nsyvme) = itdcdatayvme(nsyvme)*tdcr(1)  
          endif         
        else
          if(idat.ne.0.and.idat.ne.4095) then 
            if (istrip.le.15) then
              nsxcam = nsxcam + 1 
              isxcam(nsxcam)=istrip
              itdcdataxcam(nsxcam) = idat      
              timexcam(nsxcam) = itdcdataxcam(nsxcam)*tdcr(2)
            else if (istrip.gt.15) then
              nsycam = nsycam + 1 
              isycam(nsycam)=istrip
              itdcdataycam(nsycam) = idat      
              timeycam(nsycam) = itdcdataycam(nsycam)*tdcr(2)           
            endif  
          endif
        endif   
      enddo
*
      if (debug(3)) then
        print *,'Strip from VME, CAMAC =',nvme, ncam
        print *,'  '
        print *,'X VME strips =',nsxvme
        do i=1,nsxvme
          print *,i,' strip ',isxvme(i),
     >    ' tdc=',itdcdataxvme(i),' time=',timexvme(i)
        enddo
        print *,'Y VME strips =',nsyvme
        do i=1,nsyvme
          print *,i,' strip ',isyvme(i),
     >    ' tdc=',itdcdatayvme(i),' time=',timeyvme(i)
        enddo 
        print *,'X CAMAC strips =',nsxcam
        do i=1,nsxcam
          print *,i,' strip ',isxcam(i),
     >    ' tdc=',itdcdataxcam(i),' time=',timexcam(i)
        enddo
        print *,'Y CAMAC strips =',nsycam
        do i=1,nsycam
          print *,i,' strip ',isycam(i),
     >    ' tdc=',itdcdataycam(i),' time=',timeycam(i)
        enddo 

      endif
*
      RETURN
      END        
*
* ****************************************************
*
      SUBROUTINE RPCCLUS
*
* Define the cluster of strips in the RPC layers
*
* *****************************************************
*
      PARAMETER (nhrpcm = 80)
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
* ------------------------------------------------------
*
      if (debug(3)) then
        print *,'   '
        print *,' RPCCLUS'
        print *,'  '
      endif
*
      ncxvme = 0 
      ncyvme = 0
      ncxcam = 0 
      ncycam = 0      
      call vzero(ifsxvme,nmaxc)
      call vzero(ifsyvme,nmaxc)
      call vzero(ifsxcam,nmaxc)
      call vzero(ifsycam,nmaxc)
      call vzero(mulcxvme,nmaxc)
      call vzero(mulcyvme,nmaxc)
      call vzero(mulcxcam,nmaxc)
      call vzero(mulcycam,nmaxc)
*
      if (nhrpc.eq.0) RETURN
*
      do j=1,4
*
* ** X VME layer
*
        if (j.eq.1 .and. nsxvme.ne.0) then
          CALL RPCDOCLUS(nsxvme,isxvme,ncxvme,ifscxvme,mulcxvme)
*
          if (debug(3)) then
            print *,'      '
            print *,'X VME layer ='
            print *,'Number of clusters =',ncxvme
            do i=1,ncxvme
              print *,i,'  First strip, Multiplicity =',
     >        ifscxvme(i),mulcxvme(i)
            enddo
          endif 
*
* ** Y VME layer
* 
        elseif (j.eq.2.and.nsyvme.ne.0) then
          CALL RPCDOCLUS(nsyvme,isyvme,ncyvme,ifscyvme,mulcyvme)
*
          if (debug(3)) then
            print *,'      '
            print *,'Y VME layer ='
            print *,'Number of clusters =',ncyvme
            do i=1,ncyvme
              print *,i,'  First strip, Multiplicity =',
     >        ifscyvme(i),mulcyvme(i)
            enddo
          endif 
*
* ** X CAMAC layer
*
        elseif (j.eq.3 .and. nsxcam.ne.0) then
          CALL RPCDOCLUS(nsxcam,isxcam,ncxcam,ifscxcam,mulcxcam)
*
          if (debug(3)) then
            print *,'      '
            print *,'X CAMAC layer ='
            print *,'Number of clusters =',ncxcam
            do i=1,ncxcam
              print *,i,'  First strip, Multiplicity =',
     >        ifscxcam(i),mulcxcam(i)
            enddo
          endif  
*
* ** Y CAMAC layer
* 
        elseif (j.eq.4.and.nsycam.ne.0) then
          CALL RPCDOCLUS(nsycam,isycam,ncycam,ifscycam,mulcycam)              
*
          if (debug(3)) then
            print *,'      '
            print *,'Y CAMAC layer ='
            print *,'Number of clusters =',ncycam
            do i=1,ncycam
              print *,i,'  First strip, Multiplicity =',
     >        ifscycam(i),mulcycam(i)
            enddo
          endif 
*
        endif
      enddo
*
      RETURN
      END
*
* ***********************************************************
*
      SUBROUTINE RPCDOCLUS(mul,istrip,nc,fsc,mulc)
*
*     Here the clusters are found
*
* *********************************************************
*
      PARAMETER (nmaxc = 10)
      INTEGER mulc(nmaxc),fsc(nmaxc),pos(64),istrip(64)
*
* -------------------------------------------------------
*
      call vzero(pos,16)
      call vzero(mulc,16)       
*
      nc = 0
      do i=1,mul
        pc = pos(i)
        if (pc.eq.0) then
          nc = nc + 1
          pos(i) = nc 
          mulc(nc) = mulc(nc) + 1 
          next = 1
          do jc = 1,mul-1 
            istripnext = istrip(i)+next
            istripprev = istrip(i)+next-1
            if (istrip(i+jc).eq.istripnext) then
              pos(i+jc)= nc
              mulc(nc) = mulc(nc)+ 1
              next = next + 1
            elseif (istrip(i+jc).eq.istripprev) then
              pos(i+jc) = nc
            else
              goto 20
            endif
          enddo
        endif
 20     continue
      enddo
*
      do kc=1,nc
        fsc(kc) = istrip(1)
        countst = 0
        nstr = 0
        do ks=1,mul
          pc = pos(ks)
          if (pc.eq.kc) then
            nstr = nstr + 1
            if (nstr.eq.1) fsc(kc)= istrip(ks)
          endif  
        enddo   
      enddo  
*
      RETURN
      END
*
*************************************************************
*
         SUBROUTINE DRIFTCHTDC
*
* Decode the TDC data finding the drift time and the 
* space position in the cell.
* Codification of the TDC word:
* ____ __ ... ______ ..... ______ ........ ___ 
* |0|0|  empty   |  channel  |    TDC data    |          
* |_|_|__ ... ___|__ ..... __|___ ........ ___|     
* |   |29 ...  23|22 ..... 16| 15 ........   0| bits
*
*************************************************************
* 
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      PARAMETER (nhdctot = 200)
      INTEGER dclen,tdcdc(nhdctot)    
      DIMENSION pos(8,2*nhdcm),zpos(8,2*nhdcm),dlen(8,2*nhdcm),nhitl(8)
      REAL dtime,dlength
      INTEGER lay(8)
      COMMON /driftchgeo/  wire(0:95),zwire(0:95)
      COMMON /hits/ xpos(4,2*nhdcm),zxpos(4,2*nhdcm),xdlen(4,2*nhdcm),
     >              ypos(4,2*nhdcm),zypos(4,2*nhdcm),ydlen(4,2*nhdcm),
     >              nxhit(4),nyhit(4)
      COMMON /datadc/ dclen,tdcdc
      COMMON /t0/ t0corr(96)
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
*
* ------------------------------------------------------------ 
*
      if (debug(4)) then
        print *,'      '
        print *,' DRIFTCHTDC  '
        print *,'      '
      endif  
*
      good_xtrack = .true.
      good_ytrack = .true.           
*
      call vzero(pos,8*2*nhdcm)
      call vzero(xpos,4*2*nhdcm)
      call vzero(ypos,4*2*nhdcm)
      call vzero(zxpos,4*2*nhdcm)
      call vzero(zypos,4*2*nhdcm)
      call vzero(nhitl,8)
      call vzero(ichan,nhdcm)
      call vzero(layer,nhdmc)
      call vzero(itdcdatadc,nhdcm)
      call vzero(hitwire,nhdcm)
      call vzero(zhitwire,nhdmc)  
      call vzero(lay,8)
      call vzero(nxhit,4)
      call vzero(nyhit,4)
      call vzero(nxhdc,4)
      call vzero(nyhdc,4)
      nxlayer = 0
      nylayer = 0 

      nhdc = dclen
      if (debug(4)) print *,'Number of hits =',nhdc      
*
      if (nhdc.eq.0) then
        good_xtrack =.false.
        good_ytrack =.false.
        RETURN
      endif   
*
      if (debug(4)) then    
        print *,'Drift velocity '
        print *,'X layers - top, bottom =', vd(2),vd(1)
        print *,'Y layer - top, bottom =',vd(4),vd(3) 
        print *,'  '
      endif

      do i=1,nhdc
        ichan(i) = jbyt(tdcdc(i),17,7)
        layer(i) = ichan(i)/12 + 1
        lay(layer(i)) = lay(layer(i)) + 1
        itdcdatadc(i) = jbyt(tdcdc(i),1,16)
        hitwire(i) = wire(ichan(i))
        zhitwire(i) = zwire(ichan(i))
        if (itdcdatadc(i) .gt. 2047) then
          if (debug(1)) 
     >    write (*,1000) i,layer(i),ichan(i),itdcdatadc(i)
 1000     format (1x,'- Hit',i2,1x,'layer',i2,1x,'cell',i3,
     >    'tdcdata',i5,' TDC data too large')
        endif
*
* ** Time vs space
*
        dtime0(i) = 1700.-tdcr(1)*itdcdatadc(i)
        dtime(i) = dtime0(i)-t0corr(ichan(i)+1)    
        if (dtime(i).ge.0) then 
          if (layer(i).le.2) then
            vdrift = vd(1)
	    vdxbottom = vdrift
	  elseif (layer(i).eq.3.or.layer(i).eq.4) then
	    vdrift = vd(3)
	    vdybottom =vdrift
          elseif (layer(i).eq.5.or. layer(i).eq.6) then
            vdrift = vd(2)
	    vdxtop = vdrift
          elseif (layer(i).ge.7) then
	    vdrift = vd(4)
	    vdytop = vdrift
          endif
          dlength(i) = dtime(i)*vdrift	  
*
          if (debug(4)) then
            print *,'** Hit number',i
            print *,'layer =',layer(i),' - cell =',ichan(i)
            print *,'tdc =',itdcdatadc(i),' - time, timet0 =',
     >               dtime0(i),dtime(i) 
            print *,'wire, zwire =',wire(ichan(i)),zwire(ichan(i)),
     >               ' - drift length =',dlength(i)     
          endif  
          if (debug(9)) print *,'tdc word =',tdcdc(i) 
*       
          nhitl(layer(i)) = nhitl(layer(i)) + 1    
          dlen(layer(i),nhitl(layer(i))) = dlength(i)  
          pos(layer(i),nhitl(layer(i))) = wire(ichan(i)) + dlength(i) 
          zpos(layer(i),nhitl(layer(i))) = zwire(ichan(i)) 
          nhitl(layer(i)) = nhitl(layer(i)) + 1   
          dlen(layer(i),nhitl(layer(i))) = dlength(i)  
          pos(layer(i),nhitl(layer(i))) = wire(ichan(i)) - dlength(i) 
          zpos(layer(i),nhitl(layer(i))) = zwire(ichan(i)) 
*
        else
          lay(layer(i)) = lay(layer(i)) - 1   
          if (debug(4)) then
            print *,'** Hit number',i
            print *,'layer =',layer(i),' - cell =',ichan(i)
            print *,'tdc =',itdcdatadc(i),' - time, timet0 =',
     >               dtime0(i),dtime(i) 
            print *,'Negative time. Hit skipped'
          endif  
          if (debug(5)) print *,'tdc word =',tdcdc(i) 
        endif  
      enddo
*
* ** require at least 3 layer in X and 3 in Y with 1 hit each
*
      if (debug(4)) print *,'     ' 
      do i=1,8
        if (debug(4)) print *,'Layer',i,' hit =',lay(i)   
        if (lay(i).eq.0) then
          if (debug(4)) print *,' ** Layer',i,' missing'
        else
          if (i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) then
            nxlayer = nxlayer + 1
          else
            nylayer = nylayer + 1
          endif    
        endif
      enddo 
      if (debug(4)) then
        print *,'  '
        print *,'X layer = ',nxlayer
        print *,'Y layer = ',nylayer
      endif   
*
      if (nxlayer.lt.3) then 
        if (debug(1)) print *,' - Less than 3 X layer. No X track'
        good_xtrack = .false.
      endif   
      if (nylayer.lt.3) then
        if (debug(1)) print *,' - Less than 3 Y layer. No Y track'
        good_ytrack = .false.
      endif
*   
      do i=1,4
        if (i.eq.1.or.i.eq.2) then
          nxhit(i) = nhitl(i)
          do j=1,nhitl(i)
            xpos(i,j) = pos(i,j)
            zxpos(i,j) = zpos(i,j)
            xdlen(i,j) = dlen(i,j)             
          enddo
          nxhit(i+2) = nhitl(i+4)
          do j=1,nhitl(i+4)
            xpos(i+2,j) = pos(i+4,j)
            zxpos(i+2,j) = zpos(i+4,j)
            xdlen(i+2,j) = dlen(i+4,j) 
          enddo      
        elseif (i.eq.3.or.i.eq.4) then
          nyhit(i-2) = nhitl(i)
          do j=1,nhitl(i)
            ypos(i-2,j) = pos(i,j)
            zypos(i-2,j) = zpos(i,j)
            ydlen(i-2,j) = dlen(i,j) 
          enddo
           nyhit(i) = nhitl(i+4)
          do j=1,nhitl(i+4)
            ypos(i,j) = pos(i+4,j)
            zypos(i,j) = zpos(i+4,j)
            ydlen(i,j) = dlen(i+4,j)
          enddo   
        endif
      enddo
*
      do i=1,4
         anx = float(nxhit(i))/2.
         any = float(nyhit(i))/2.
         nxhdc(i)=int(anx)
         nyhdc(i)=int(any)        
      enddo    
*
      if (debug(4)) then
        print *,'    '
        print *,'*** X plane'
        do j=1,4
          print *,'layer',j,' Number of hits =',nxhit(j)
          do i=1,nxhit(j)
            print *,i,'  xpos, zpos, xdlen =',
     >      xpos(j,i),zxpos(j,i),xdlen(j,i)
          enddo        
        enddo  
        print *,'  '
        print *,'*** Y plane'
        do j=1,4
          print *,'layer',j,' Number of  hits =',nyhit(j)
          do i=1,nyhit(j)
            print *,i,'  ypos, zpos, ydlen =',
     >      ypos(j,i),zypos(j,i),ydlen(j,i)
          enddo        
        enddo  
      endif  
*
      RETURN
      END 
*
***********************************************************************
*
      SUBROUTINE DRIFTCHPAT
*
*     Pattern recognition
*
* *********************************************************************
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      PARAMETER (pi = 3.1415926) 
      COMMON /hits/ xpos(4,2*nhdcm),zxpos(4,2*nhdcm),xdlen(4,2*nhdcm),
     >              ypos(4,2*nhdcm),zypos(4,2*nhdcm),ydlen(4,2*nhdcm),
     >              nxhit(4),nyhit(4)      
      COMMON /fit/ a,b,chi2,chi2red
      PARAMETER (ntrackm=150)
      COMMON /tracks/ ntrack,nlh(4),pos_track(ntrackm,4),
     >                fit_track(ntrackm,4),res_track(ntrackm,4),
     >                chi2_track(ntrackm),theta_track(ntrackm),
     >                a_track(ntrackm),b_track(ntrackm),
     >                nfit_track(ntrackm),ihhist_track(ntrackm,4),
     >                itrackflag(ntrackm)
      INTEGER nlhx_good,nlhy_good,hitx_good,hity_good
      COMMON /good/ xzpos_good(4),yzpos_good(4), 
     >              nlhx_good(4),nlhy_good(4),hitx_good(4),hity_good(4)
      DIMENSION nlayerx(4),nlayery(4),xdl(4),ydl(4)
*
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
      LOGICAL dcl,rpcl,scil,patl,trackl
      COMMON /logic/ dcl,rpcl,scil,patl,trackl
*
* --------------------------------------------------------
*
      if (debug(4)) then
        print *,'  '
        print *,' DRIFTCHPAT '
        print *,'   '
      endif  
*
* ** Put to zero the ntuples variables
*
      thetax_good = 0
      thetay_good = 0
      ax_good = 0
      ay_good = 0
      bx_good = 0
      by_good = 0
      chi2redx_good = 0
      chi2redy_good = 0
      nfitx_good = 0
      nfity_good = 0
      call vzero(xres_good,nfitm)
      call vzero(yres_good,nfitm)
      call vzero(xpos_good,nfitm)
      call vzero(ypos_good,nfitm)
      ixtrack = 0
      iytrack = 0
*      
      if (good_xtrack) ixtrack = 1
      if (good_ytrack) iytrack = 1

*              
      if (debug(4)) then
        print *,'ixtrack, iytrack =',ixtrack,iytrack
        print *,'     '
      endif
*
* ** Pattern recognition in X and Y planes (ixy =1, ixy=2)
*           
      do ixy =1, 2
        if (ixy.eq.1 .and. ixtrack.eq.1) then
*
          if (debug(4)) print *,'*** X plane ***' 
*
* ** Perform Histogram method
*
          if (debug(4)) print *,'nxhit =',nxhit

          CALL DRIFTCHISTO(nxhit,xpos,zxpos,ixy)   
*
* ** Select the best track
* 
          if (ntrack.ne.0) then
            if (debug(4)) then
              print *,'    '
              print *,'*** Select the best track ***'
              print *,'Number of tracks = ',ntrack    
              do i=1,ntrack
                print *,i,' theta, chi2red flag =',
     >                theta_track(i)*180./3.1415,chi2_track(i),
     >          itrackflag(i),' - a, b =',a_track(i), b_track(i)
              enddo     
            endif
            nsel = 1
            do while ((itrackflag(nsel).eq.0).and.(nsel.le.ntrack))
               nsel = nsel + 1
            enddo
            if (debug(4)) print *,'nsel ',nsel
            if (nsel.gt.ntrack) then
               ixtrack = 0
               print *,'No X Track cross scintillators'
               goto 332
            endif
            do i=nsel+1,ntrack
              if ((itrackflag(i).eq.1).and.
     >              (chi2_track(i).lt.chi2_track(nsel))) then
                nsel = i
              endif
            enddo   
          else 
            ixtrack = 0 
            if (debug(4)) 
     >      print *,'No X Track reconstructed by histo'
            goto 332
          endif 
*      
          nxgood_track = nsel
          nfitx_good = nfit_track(nxgood_track)
          chi2redx_good = chi2_track(nxgood_track) 
          thetax_good = theta_track(nxgood_track)
          ax_good = a_track(nxgood_track)
          bx_good = b_track(nxgood_track)
          do i=1,4
            nlhx_good(i) = nlh(i)
            hitx_good(i) = ihhist_track(nxgood_track,i)
          enddo
*
          if (debug(4)) then
            print *,'  '
            print *,'Selected track =',nxgood_track
            print *,'chi2 flag =',chi2redx_good,
     >           itrackflag(nxgood_track),' - theta =',
     >      thetax_good*180./acos(-1.),' - a, b =',ax_good,bx_good
          endif 
*
          if (debug(4)) print *,'Number of hits =',nfitx_good
          j=0
          do i=1,4
            if (nlhx_good(i).eq.1) then
              j=j+1
              layerx_good(j) = i
              xres_good(j) = res_track(nxgood_track,i)
              xpos_good(j) = pos_track(nxgood_track,i)
              xzpos_good(j) = zxpos(i,1)
              if (debug(4)) print *,'Layer =',layerx_good(j),
     >           ' - hit =',hitx_good(i),
     >           ' - xpos =',xpos_good(j),' - xres =',xres_good(j),
     >           ' - zpos =',xzpos_good(j) 
            else
              if (debug(4)) print *,'Layer =',i,' - No hits'
            endif         
          enddo
*
* ** Perform Refit if needed
*
          correctx = 0

          if (flag_ref.eq.1 .and.chi2redx_good.gt.chi2cut) then                
*         
            if (debug(4)) then
              print *,'   ' 
              print *,'Refit of the track'
              print *,' '
            endif   
*
* ** Sort the hits accordingly to drift lengths
*
            do i=1,4
              if (nlhx_good(i).eq.1) then
                xdl(i) = xdlen(i,hitx_good(i))  
cm                print *,' xdlen =  ',xdl(i),'  ',i
              else
                xdl(i) = 999.
              endif
            enddo              
cm            print *,''
            do m = 1,4
              nlayerx(m) = m
            enddo
            do i=1,3
              do j=i+1,4
                if (xdl(i).gt.xdl(j)) then
                  temp = xdl(i)
                  xdl(i) = xdl(j)
                  xdl(j) = temp
                  itemp = nlayerx(i)
                  nlayerx(i) = nlayerx(j)
                  nlayerx(j) = itemp  
                endif
              enddo
            enddo
*            
            call DRIFTCHREFIT(ixy,xdl,nlayerx)
*
            if (debug(4)) then
              print *,'  '
              print *,'Selected track =',nxgood_track
              print *,'chi2 flag =',chi2redx_good,
     >             itrackflag(nxgood_track),' - theta =',
     >        thetax_good*180./acos(-1.),' - a, b =',ax_good,bx_good
              j=0
              do i=1,4
                if (nlhx_good(i).eq.1) then
                  j=j+1
                  print *,'Layer =',i,
     >            ' - hit =',hitx_good(i),
     >            ' - xpos =',xpos_good(j),' - xres =',xres_good(j),
     >            ' - zpos =',xzpos_good(j) 
                else
                  print *,'Layer =',i,' - No hits'
                endif         
              enddo  
            endif  
*
          endif   
*
          if (chi2redx_good.gt.1000) then
            if(debug(4).or.debug(1)) 
     >      print *,' - chi2 X track =',chi2redx_good,' too large.'            
            ixtrack = 0 
          endif        
*
 332      continue

        else if (ixy.eq.2 .and. iytrack.eq.1) then
*
          if (debug(4)) then
            print *,'   '
            print *,'*** Y plane ***'
          endif  
*
* ** Perform Histogram method
*
          CALL DRIFTCHISTO(nyhit,ypos,zypos,ixy)   
*
* ** Select the best track
* 
          if (ntrack.ne.0) then
            if (debug(4)) then
              print *,'    '
              print *,'Select the best track'
              print *,'selected tracks = ',ntrack    
              do i=1,ntrack
                 print *,i,' theta, chi2red flag =',
     >         theta_track(i)*180./3.1415,chi2_track(i),
     >         itrackflag(i),' - a, b =',a_track(i), b_track(i)
              enddo     
            endif
            nsel = 1
            do while ((itrackflag(nsel).eq.0).and.(nsel.le.ntrack))
               nsel = nsel + 1
            enddo
            if (debug(4)) print *,'nsel ',nsel
            if (nsel.gt.ntrack) then
               iytrack = 0
               if (debug(4)) print *,'No Y Track cross scintillators'
               goto 333
            endif
            do i=nsel+1,ntrack
              if ((itrackflag(i).eq.1).and.
     >         (chi2_track(i).lt.chi2_track(nsel))) then
                nsel = i
              endif
            enddo  
          else 
            iytrack = 0 
            if (debug(4)) 
     >      print *,'No Y Track reconstructed by histo'
            goto 333
          endif 
*       
          nygood_track = nsel
          nfity_good = nfit_track(nygood_track)
          chi2redy_good = chi2_track(nygood_track) 
          thetay_good = theta_track(nygood_track)
          ay_good = a_track(nygood_track)
          by_good = b_track(nygood_track)
          do i=1,4
            nlhy_good(i) = nlh(i)
            hity_good(i) = ihhist_track(nygood_track,i)
          enddo
*
          if (debug(4)) then
            print *,'  '
            print *,'Selected track =',nygood_track
            print *,'chi2 flag =',chi2redy_good,
     >           itrackflag(nygood_track),' - theta =',
     >      thetay_good*180./acos(-1.),' - a, b =',ay_good, by_good
          endif 
*
          if (debug(4)) print *,'Number of hits =',nfity_good
          j=0
          do i=1,4 
            if (nlhy_good(i).eq.1) then
              j=j+1
              layery_good(j) = i
              yres_good(j) = res_track(nygood_track,i)
              ypos_good(j) = pos_track(nygood_track,i)
              yzpos_good(j) = zypos(i,1)
              if (debug(4)) print *,'Layer =',layery_good(j), 
     >        ' - hit =',hity_good(i),
     >        ' - ypos =',ypos_good(j),' - yres =',yres_good(j),
     >        ' - zpos =',yzpos_good(j) 
            else
              if (debug(4)) print *,'Layer =',i,' - No hits' 
            endif  
          enddo
*
* ** Perform Refif if needed
*
          correcty = 0

          if (flag_ref.eq.1 .and.chi2redy_good.gt.chi2cut) then         
*
            if (debug(4)) then
              print *,'   '
              print *,'Refit of the track'
              print *,'    '    
            endif
*         
* ** Sort the hits accordingly to drift lengths
*
            do i=1,4
              if (nlhy_good(i).eq.1) then
                ydl(i) = ydlen(i,hity_good(i))  
              else
                ydl(i) = 999.
              endif
            enddo              
            do m = 1,4
              nlayery(m) = m
            enddo
            do i=1,3
              do j=i+1,4
                if (ydl(i).gt.ydl(j)) then
                  temp = ydl(i)
                  ydl(i) = ydl(j)
                  ydl(j) = temp
                  itemp = nlayery(i)
                  nlayery(i) = nlayery(j)
                  nlayery(j) = itemp  
                endif
              enddo
            enddo
*           
            call DRIFTCHREFIT(ixy,ydl,nlayery)
*
            if (debug(4)) then
              print *,'  '
              print *,'Selected track =',nygood_track
              print *,'chi2 =',chi2redy_good,' - theta =',
     >        thetay_good*180./acos(-1.),' - a, b =',ay_good,by_good
              j=0
              do i=1,4
                if (nlhy_good(i).eq.1) then
                  j=j+1
                  print *,'Layer =',i,
     >            ' - hit =',hity_good(i),
     >            ' - xpos =',ypos_good(j),
     >            ' - yres =',yres_good(j),
     >            ' - zpos =',yzpos_good(j) 
                else
                  print *,'Layer =',i,' - No hits'
                endif         
              enddo
            endif           
*
          endif   
*
          if (chi2redy_good.gt.1000) then
            if(debug(4)) print *,'chi2 Y track=',chi2redy_good,
     >      ' too large' 
            iytrack = 0
          endif        
*
 333      continue

        endif
      enddo      
*
      if (debug(4)) then
        print *,'  '
        if (ixtrack.eq.1) then
          print *,' ** X track reconstructed'
        else
          print *,' ** X track NOT reconstructed'
        endif
        if (iytrack.eq.1) then
          print *,' ** Y track reconstructed'
        else
          print *,' ** Y track NOT reconstructed'
        endif
      endif  
*
      if (ixtrack.eq.0 .and. iytrack.eq.0) then
        trackl = .false.
        if (debug(1)) print *,' - No tracks reconstructed' 
      else
        trackl = .true. 
        itrack = 1       
*
        th = sqrt(tan(thetax_good**2)+tan(thetay_good**2))
        theta3D = atan(th)        
        if (tan(thetax_good).ne.0) then
           ph = tan(thetay_good)/tan(thetax_good)
*          if (tan(thetax_good).gt.0.and.tan(thetay_good).gt.0)  then
*            phi3D = atan(ph)  
*          elseif (tan(thetax_good).lt.0.and.tan(thetay_good).gt.0) then
*            phi3D = pi+atan(ph)
*          elseif (tan(thetax_good).lt.0.and.tan(thetay_good).lt.0) then 
*            phi3D = pi+atan(ph)
*          elseif (tan(thetax_good).gt.0.and.tan(thetay_good).lt.0) then
*            phi3D = 2*pi+atan(ph)          
*          endif 
           if (tan(thetax_good).gt.0.and.tan(thetay_good).gt.0)  then
              phi3D = atan(ph)  
           elseif (tan(thetax_good).lt.0.and.tan(thetay_good).gt.0) then
              phi3D = pi+atan(ph)
           elseif (tan(thetax_good).lt.0.and.tan(thetay_good).lt.0) then 
              phi3D = -pi+atan(ph)
           elseif (tan(thetax_good).gt.0.and.tan(thetay_good).lt.0) then
              phi3D = +atan(ph)          
           endif 
        endif
*  
      endif          
*
      RETURN
      END
*
* ***************************************************************
*
      SUBROUTINE DRIFTCHISTO(nhit,pos,zpos,ixy)           
*
* Perform a pattern recognition in the drift chambers
* applying the Histogram Method
*
* ****************************************************************
*
      INCLUDE "dc.inc"

*      DATA slpmin /-0.35/, x0min /0./
      DIMENSION ihist(0:ns-1,0:nx-1)
      DIMENSION diffmini(4),xyhist(4),ihhist(4),zhist(4)
      DIMENSION tfit(4),res(4)
      PARAMETER (ntrackm=150)
      PARAMETER (nhdcm = 32)
      DIMENSION nhit(4),pos(4,2*nhdcm),zpos(4,2*nhdcm)
      COMMON /CONST/ sdy,sy,stsbz,sb0z,d0x,d0y,cdx,cdy
      COMMON /tracks/ ntrack,nlh(4),pos_track(ntrackm,4),
     >                fit_track(ntrackm,4),res_track(ntrackm,4),
     >                chi2_track(ntrackm),theta_track(ntrackm),
     >                a_track(ntrackm),b_track(ntrackm),
     >                nfit_track(ntrackm),ihhist_track(ntrackm,4),
     >                itrackflag(ntrackm)
      COMMON /fit/ a,b,chi2,chi2red
*
* -------------------------------------------------------
*
      if (debug(5)) then
        print *,'  '
        print *,' DRIFTCHISTO '
        print *,'   '
      endif  

**  
      if (ixy.eq.1) then
         slpmin = ((xvt-50.)-(xvb+1000.+50.))/stsbz     
         slpmax = ((xvt+1000.+50.)-(xvb-50.))/stsbz
         x0min  = (xvb+sdy-50.) - slpmax*sb0z
         x0max  = (xvb+sdy+1000.+50.) - slpmin*sb0z
      elseif (ixy.eq.2) then
         slpmin = -(1000.+100.)/stsbz     
         slpmax =  (1000.+100.)/stsbz
         x0min  = sdy - slpmax*sb0z
         x0max  = sdy+1000. - slpmin*sb0z  
      endif

***  50. e 100. sono margini sui bordi degli scintillatori ***      
    
      diffslop = slpmax-slpmin
      diffx0 = x0max-x0min

      sstep = diffslop/ns
      xstep = diffx0/nx
*
* ** Histogram method 
*       
*      sstep = 3.6/ns
*      xstep = 1140./nx
     

*
      call vzero(ihist,ns*nx)
      if (debug(5)) then
        print *,'    '
        print *,'Histogram method'
        print *,'Number of bins =',ns, nx
        print *,'Step =',sstep, xstep
        print *,'Min Max, slop and inter =', slpmin,x0min,slpmax,x0max
      endif  
*
      if (debug(6)) then
         do i=0,ns-1
            binsmin = slpmin + i*sstep
            binsmax = slpmin + (i+1)*sstep
            print *,'slope bin, estremi =',i,binsmin, binsmax 
         enddo
         do i=0,nx-1
            binxmin = xmin + i*xstep
            binxmax = xmin + (i+1)*xstep
            print *,'x bin, estremi =',i,binxmin, binxmax 
         enddo
      endif
*
* ** Loop over the cells.
* ** Slopes and interceptas of all the segments determined by each 
* ** pair of hits are loaded in a bidimentional histogram.
*
      do icell = 1,3
        if (debug(5)) print *,'nhit',nhit(icell)
        do ih = 1,nhit(icell)
          x1 = pos(icell,ih)
          z1 = zpos(icell,ih)
          if (debug(6)) print *,'ih, x1, z1 =',ih,x1,z1
            do icell2 = icell+1,4
              if (debug(6)) print *,'** icell2 =',icell2
               do jh=1,nhit(icell2)
                x2 = pos(icell2,jh)
                z2 = zpos(icell2,jh)
                if (debug(6)) print *,'*** x1,x2,z1, z2 =',x1, x2,z1, z2       
                  slp = (x2-x1)/(z2-z1)
              x0=x1-slp*(z1)
              islp=(slp-slpmin)/sstep
              if ((slp-slpmin).lt.0.) islp = -1
              ix0 = (x0-x0min)/xstep 
              if (debug(6)) print *,'*** slp, islp  =',slp,islp 
              if (debug(6)) print *,'*** x0, ix0 =',x0,ix0
              if (islp.ge.0 .and. islp.le.(ns-1) .and.
     >        ix0.ge.0 .and.ix0.le. (nx-1)) 
     >        ihist(islp,ix0) = ihist(islp,ix0) + 1               
            enddo   
          enddo  
        enddo  
      enddo 
*  
      if (debug(6)) then
        do i=0,ns-1
          print *,' ihist =',i,(ihist(i,j),j=0,nx-1)
        enddo
      endif  
*
      if (debug(5)) then
        print *,'       '
        print *,'Filled bins in histo'
        do i=0,ns-1
           do j=0,nx-1
             if (ihist(i,j).gt.0) then
               print *,'nslop, nx0 =',i,j,' hit =',ihist(i,j)
             endif
           enddo
         enddo
       endif  
*
      ntrack = 0 
      do i=0,ns-1
         do j=0,nx-1
            if (ihist(i,j).le.1) goto 555 
            ntrack = ntrack + 1
            nh = 0
         
            sl = slpmin + (i+0.5)*sstep
            x0 = x0min + (j+0.5)*xstep
            nh = ihist(i,j)
            
            if (debug(4)) then
               print *,'    '  
               print *,'Track number ',ntrack
            endif  
            if (debug(5)) then
               print *,'        '
               print *,'central bins: slop, x0 =',i,j
               print *,'hits =',nh
               print *,'Average slope , x0 =',sl,x0
            endif
* 
* ** Find and pick the good hits
*      
          do icell=1,4
            diffmini(icell) = 1.0e6
            xyhist(icell) = 0.
            zhist(icell) = 0.
            ihhist(icell) = 0      
            if (debug(5)) print *,'icell =',icell
            do ih =1,nhit(icell)   
              zhit = zpos(icell,ih)
              hist = x0+sl*zhit
              hpos = pos(icell,ih)
              diff = hpos-hist
              if (debug(5)) print *,ih,' hist, hpos, zhit, diff =',
     >        hist, hpos,zhit,diff 
              if (abs(diff).lt.diffmini(icell)) then  
                diffmini(icell)=abs(diff)
                xyhist(icell) = hpos
                zhist(icell) = zhit
                ihhist(icell) = ih
              endif
            enddo    
          enddo
*
          if (debug(4)) then
            print *,'  '
            print *,'Picked hits'
            do icell=1,4
              if (nhit(icell).ne.0) then
                print *,'layer =',icell,
     >          ' - pos =',xyhist(icell),' - diff =',diffmini(icell),
     >          ' - id hit =',ihhist(icell)
              else
                print *,'layer =',icell,' - No hits'
              endif  
            enddo   
          endif  
*
          nfit = 0
          call vzero(nlh,4)
          do icell=1,4
            if (ihhist(icell).ne.0) then
              nfit = nfit+1
              nlh(icell)=1
            endif
          enddo
*
* ** Fit the track
*
          CALL DRIFTCHFIT(xyhist,zhist,nlh,nfit)
          thetat = atan(a)
          if (debug(4)) then
            print *,'  '
            print *,'Fit - Number of Hits =',nfit
            print *,'a, b, theta =',a,b,thetat*180./acos(-1.) 
            print *,'chi2, chi2red =',chi2,chi2red
          endif
*      
* ** Residuals
*
          if (debug(4)) print *,'  ' 
          do ij=1,4
            if (nlh(ij).eq.1) then
              tfit(ij) = a*zhist(ij)+b 
              res(ij) = tfit(ij)-xyhist(ij)
              if (debug(4)) print *,'fit =',tfit(ij),
     >        ' - hit xy, z =',xyhist(ij),zhist(ij),' - res =',res(ij)
            endif   
          enddo
*

          nfit_track(ntrack)  = nfit
          chi2_track(ntrack)  = chi2red
          theta_track(ntrack) = thetat
          a_track(ntrack) = a
          b_track(ntrack) = b

          if (ixy.eq.1) then
             CALL TRACKONX
             
          elseif (ixy.eq.2) then
             CALL TRACKONY
             
          endif 



 
          do icell = 1,4
            ihhist_track(ntrack,icell) = ihhist(icell)
          enddo
          do ij=1,4
            if (nlh(ij).eq.1) then
              pos_track(ntrack,ij) = xyhist(ij)
              fit_track(ntrack,ij) = tfit(ij)
              res_track(ntrack,ij) = res(ij)            
            endif  
          enddo
*
 555      continue
        enddo
      enddo
*
      RETURN
      END
*
* ********************************************************************
*
      SUBROUTINE DRIFTCHFIT(xpos,zpos,nlh,nfit)
*
* Fit of the muon tracks in the drift chambers
*
* *********************************************************************
*
      INTEGER nlh(4)
      DIMENSION xpos(4),zpos(4)
      REAL chi2,chi2red
      COMMON /fit/ a,b,chi2,chi2red
*
* ------------------------------------------------------
*  
* ** Linear fit: x = a*zx + b (unweighted)
*  
      sz=0             
      sx=0
      st2=0
      a=0
      ss=0
*
      nlay = 4
      do i=1,nlay
        if (nlh(i).eq.1) then
          sz=sz+zpos(i)
          sx=sx+xpos(i)
        endif  
      enddo
      ss=float(nfit)
      szoss=sz/ss
      do i=1,nlay
        if (nlh(i).eq.1) then
          t=zpos(i)-szoss
          st2=st2+t*t
          a=a+t*xpos(i)
        endif  
      enddo
*
* ** determine a and b
*
      a=a/st2
      b=(sx-sz*a)/ss
      sigb=sqrt(1.+sz*sz/(ss*st2)/ss)
      siga=sqrt(1./st2)
*     
* ** determine chi2
*
      chi2=0
      do i=1,nlay
        if (nlh(i).eq.1) chi2=chi2+((xpos(i)-a*zpos(i)-b))**2
      enddo     
      q=1.
      sigdat = sqrt(chi2/(nfit-2))
      siga = siga*sigdat
      sigb = sigb*sigdat
      chi2red = chi2/(nfit-2)
       
**
** ** Fit with CERN library
**
*      
*      call lfit(zpos,xpos,4,1,a,b,var) 
*      print *,'fit with cernlib  '
*      print *,'a, b, chi2 =',a,b,var
*
      RETURN
      END
*
* ***************************************************************
*
      SUBROUTINE DRIFTCHREFIT(ixy,dl,nlayer)           
*
* Refit of the track.
*
* ****************************************************************
*    
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"

      COMMON /hits/ xpos(4,2*nhdcm),zxpos(4,2*nhdcm),xdlen(4,2*nhdcm),
     >              ypos(4,2*nhdcm),zypos(4,2*nhdcm),ydlen(4,2*nhdcm),
     >              nxhit(4),nyhit(4)
      PARAMETER (ntrackm=150)
      COMMON /tracks/ ntrack,nlh(4),pos_track(ntrackm,4),
     >                fit_track(ntrackm,4),res_track(ntrackm,4),
     >                chi2_track(ntrackm),theta_track(ntrackm),
     >                a_track(ntrackm),b_track(ntrackm),
     >                nfit_track(ntrackm),ihhist_track(ntrackm,4),
     >                itrackflag(ntrackm)
      INTEGER nlhx_good,nlhy_good,hitx_good,hity_good
      COMMON /good/ xzpos_good(4),yzpos_good(4), 
     >              nlhx_good(4),nlhy_good(4),hitx_good(4),hity_good(4)
      COMMON /fit/ a,b,chi2,chi2red

      INTEGER hit_ref(4),hit_good(4),nlayer(4) 
      DIMENSION pos_ref(4),zpos_ref(4),dl(4)
      DIMENSION pos_good(4),res_good(4),zpos_good(4)
      DIMENSION zpos(4,2*nhdcm),pos(4,2*nhdcm)
      INTEGER layer_good(4),nlh_good(4),nhit(4)
      INTEGER correctx,correcty
*
* --------------------------------------------------------------------
*
      call vzero(hit_ref,4)
      call vzero(pos_ref,4) 
      call vzero(zpos_ref,4) 
*
      if (ixy.eq.1) then
        a_good = ax_good
        b_good = bx_good
        theta_good = thetax_good
        nfit_good = nfitx_good
        chi2red_good = chi2redx_good
        j=0
        do i=1,4
          nhit(i) = nxhit(i)          
          do ij=1,nhit(i)
            pos(i,ij) = xpos(i,ij)
            zpos(i,ij) = zxpos(i,ij)
          enddo
          nlh_good(i) = nlhx_good(i)
          hit_good(i) = hitx_good(i)
          if (nlh_good(i).eq.1) then
            j=j+1            
            layer_good(j) = layerx_good(j)
            pos_good(j) = xpos_good(j)
            zpos_good(j) = xzpos_good(j)
            res_good(j) = xres_good(j)
          endif   
        enddo
      elseif (ixy.eq.2) then
        a_good = ay_good
        b_good = by_good
        theta_good = thetay_good
        nfit_good = nfity_good
        chi2red_good = chi2redy_good
        j=0
        do i=1,4
          nhit(i) = nyhit(i)
          do ij=1,nhit(i)
            pos(i,ij) = ypos(i,ij)
            zpos(i,ij) = zypos(i,ij) 
          enddo
          nlh_good(i) = nlhy_good(i)
          hit_good(i) = hity_good(i)
          if (nlh_good(i).eq.1) then
            j=j+1           
            layer_good(j) = layery_good(j)
            pos_good(j) = ypos_good(j)
            zpos_good(j) = yzpos_good(j)
            res_good(j) = yres_good(j)
                      
          endif   
        enddo   
      endif      
*
      igoodref = 0  
      nrefit = 0
      icorrect = 0
      do while (igoodref.eq.0)
        nrefit = nrefit + 1
        j=0
        do  i = 1,4
          if (nlh_good(i).eq.1) then
            j=j+1
            if (i.eq.nlayer(nrefit)) then
              hit_ref(i) = hit_good(nlayer(nrefit))-
     >                     (-1)**hit_good(nlayer(nrefit))
              pos_ref(i) = pos(nlayer(nrefit),hit_ref(i))
              zpos_ref(i) = zpos(nlayer(nrefit),hit_ref(i))
            else 
              pos_ref(i) = pos_good(j)
              zpos_ref(i) = zpos_good(j)
              hit_ref(i) = hit_good(i) 
            endif
          endif
*          
cm          print *,'hit_ref =',hit_ref(i),' hit_good = ',hit_good(i),
cm     >         ' hity_good =',hity_good(i)
cm          print *,''
        enddo 
*
cm        print *,'nrefit  ',nrefit
        if (debug(4)) then
          print *,'Picked hits '
          do i = 1,4
            if (nlh_good(i).eq.1) then
              print *,'layer =',i,
     >        ' - pos =',pos_ref(i),
     >        ' - hit =',hit_ref(i)
            else  
              print *,'layer =',i,' - No hits'   
            endif   
          enddo  
        endif  
*
        call driftchfit(pos_ref,zpos_ref,nlh_good
     >                  ,nfit_good)
*
        ntrack = 150
        nfit_track(ntrack)  = nfit
        chi2_track(ntrack)  = chi2red
        theta_track(ntrack) = thetat
        a_track(ntrack) = a
        b_track(ntrack) = b 
        
        if (ixy.eq.1) then
           CALL TRACKONX
        elseif(ixy.eq.2) then   
           CALL TRACKONY
        endif

        if ((chi2red.lt.chi2red_good).and.
     >       (itrackflag(ntrack).eq.1)) then
           icorrect = icorrect + 1
           if (debug(4)) then
              print *,' '
              print *,'Old and New chi2 =',chi2red_good,chi2red,
     >        itrackflag(ntrack)     
              print *,'New track accepted' 
              print *,' '
           endif 
           chi2red_good = chi2red 
           a_good = a
           b_good = b
           theta_good = atan(a) 
           j=0
           do i=1,4
              if (nlh_good(i).eq.1) then
                 j=j+1
                 pos_good(j) = pos_ref(i)
                 hit_good(i) = hit_ref(i)
                 zpos_good(j) = zpos_ref(i)
                 trfit = a_good*zpos_good(j)+b
                 res_good(j) = trfit - pos_good(j) 
              endif
           enddo 
           if (chi2red.gt.chi2cut .and. nrefit.lt.2) then
            if (debug(4)) then
              print *,'Refit second step'         
              print *,'      '
            endif        
            igoodref = 0
          else
            igoodref = 1     
          endif
       elseif ((chi2red.gt.chi2red_good).or.
     >         (itrackflag(ntrack).eq.0)) then
          if (debug(4)) then
            print *,' '
            print *,'Old and New chi2 =',chi2red_good, chi2red,
     >             itrackflag(ntrack)
            print *,'New track NOT accepted' 
            print *,' '
          endif 
          if (nrefit.lt.2) then 
            if (debug(4)) print *,'Refit second step  '   
            igoodref = 0
          else 
            igoodref = 1    
          endif    
        endif
      enddo
*

      if (ixy.eq.1) then
        chi2redx_good = chi2red_good
        ax_good = a_good
        bx_good = b_good
        thetax_good = atan(ax_good) 
        if (icorrect.ge.1) correctx = 1
        j=0
        do i = 1,4
          if (nlhx_good(i) .eq.1) then
            j= j+1
            xpos_good(j) = pos_good(j)
            xres_good(j) = res_good(j)
            xzpos_good(j) = zpos_good(j)
            hitx_good(i) = hit_good(i)  
          endif  
        enddo
      elseif (ixy.eq.2) then
        chi2redy_good = chi2red_good
        ay_good = a_good
        by_good = b_good
        thetay_good = atan(ay_good) 
        if (icorrect.ge.1) correcty = 1
        j=0 
        do i=1,4
          if (nlhy_good(i).eq.1) then
            j=j+1
            ypos_good(j) = pos_good(j) 
            yres_good(j) = res_good(j)
            yzpos_good(j) = zpos_good(j) 
            hity_good(i) = hit_good(i)  
          endif  
        enddo
      endif
*   
      RETURN
      END 
*
* ******************************************************************
*
      SUBROUTINE FINAL
*
*     Determine some final results
*
* ******************************************************************
*
*      print *,'      '
*      print *,' FINAL'
*      print *,'     '
*
* ** Mean and Standard Deviation of the Residuals Distribution
* ** and Single Wire Resolution 
*     
*      mean_resid = hstati(10,1,' ',0)
*      stdev_resid = hstati(10,2,' ',0)
*      resol = sqrt(2.)*stdev_resid
*      print *,'media, sigma resid, resol='
*     >,mean_resid,stdev_resid,resol
*
      RETURN
      END
*
**********************************************************************
*
      SUBROUTINE SCINTMATCH
*
*     Find the matching of the scintillators with the
*     tracks reconstracted in the drift chambers
* 
**********************************************************************
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      INTEGER scintT, scintB
*
* --------------------------------------------------------------
*
      if (debug(4)) then
        print *,'  '
        print *,' SCINTMATCH'
        print *,'   '
      endif
*
      if (ixtrack.eq.0) RETURN 
      sb0z = 0.
      st0z = 3245.
      sdy = 170 
      sy = 250
      d2 = 1000+sdy   
      xt = bx_good+ax_good*st0z
      xb = bx_good+ax_good*sb0z
*      x01 = bx_good+ax_good*29./2.
*      x04 = bx_good+ax_good*(521.+29.+29./2.)
*
      scintT=0
      scintB=0 
      if ((xt.ge.sdy).and.(xt.le.d2)) scintT =int((xt-sdy)/sy)+1
      if ((xb.ge.sdy).and.(xb.le.d2)) scintB =int((xb-sdy)/sy)+1
      iscintexp(1)=scintB
      iscintexp(2)=scintT+4
      nscintexp = 2
*
      if (debug(4)) then
        print *,'scint Top, Bottom =',scintT,scintB
        print *,'scintexp =',iscintexp
      endif  
*      
      RETURN
      END
*
* ***********************************************************
*
      SUBROUTINE RPCMATCH
*
*     Determine the strip hit by the DC tracks
*
* **************************************************************
*
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
      LOGICAL dcl,rpcl,scil,patl,trackl
      COMMON /logic/ dcl,rpcl,scil,patl,trackl
*
      PARAMETER (istripw =30.)
      PARAMETER (offx = 342.)
      PARAMETER (offy = 352.)
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
* 
* ---------------------------------------------------------------      
*
      if (debug(2)) then
        print *,'  '
        print *,' RPCMATCH'
        print *,'  '
      endif
*
* **  Determine if the tracks cross the RPC and find the expected strips
*
      nsxexp = 0
      nsyexp = 0
      call vzero(isxexp,10)
      call vzero(isyexp,10)
*
      rpcw = 16*istripw
      if (ixtrack.eq.1) then
        ximpact = bx_good+ax_good*343.
        ximpact0 = ximpact - offx
        if (ximpact0.ge.0 .and. ximpact0.le.rpcw) then
          nsxexp = nsxexp + 1
          isxexp(nsxexp) = int(ximpact0)/istripw         
        endif    
      endif  
      if (iytrack.eq.1) then
        yimpact = by_good+ay_good*349.
        yimpact0 = yimpact - offy
        if (yimpact0.ge.0 .and. yimpact0.le.rpcw) then      
          nsyexp = nsyexp + 1
          isyexp(nsyexp) = int(yimpact0)/istripw + 16
        endif    
      endif 
*    
      if (debug(3)) then
        print *,'Expected strips'
        print *,'  X layer =',nsxexp
        do i=1,nsxexp
          print *,i,'  strip =',isxexp(i)
        enddo 
        print *,'  Y layer =',nsyexp
        do i=1,nsyexp
          print *,i,'  strip =',isyexp(i)
        enddo  
      endif      
*
      RETURN
      END
*
*
* ***************************************************************
*
      SUBROUTINE MCSIMULATION(nevsim_trigger,nevsim)     
*
*     Monte Carlo simulation of the cosmic rays flux 
*     and cosmic track reconstruction.
*
*                                      Authors:
*                                
*                             Cosmic Ray Flux : D. Piccolo
*                        Track Reconstruction : M. Caprio 
*
********************************************************************
*
      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
*
      COMMON /driftchgeo/  wire(0:95),zwire(0:95)
      COMMON /hits/ xpos(4,2*nhdcm),zxpos(4,2*nhdcm),xdlen(4,2*nhdcm),
     >              ypos(4,2*nhdcm),zypos(4,2*nhdcm),ydlen(4,2*nhdcm),
     >              nxhit(4),nyhit(4)
      COMMON /CONST/sdy,sy,stsbz,sb0z,d0x,d0y,cdx,cdy
      LOGICAL eve, endofdata,good_xtrack, good_ytrack
      COMMON /ev/ eve,endofdata,good_xtrack,good_ytrack,ievent
      REAL gauss_var,gn
      REAL muon,electron
      LOGICAL mom
      PARAMETER (cutoff_etrig = .12)
      PARAMETER (pi = 3.1415926) 
      DIMENSION pos(8,100),zpos(8,100),laymc(8)
      PARAMETER (xc=92.) 
      PARAMETER (zc=29.)
      PARAMETER (zz0=2850.) 
      PARAMETER (wibeam = 1.5)
*      PARAMETER (cdx=88.75)
*      PARAMETER (cdy=87.75)
*
* ------------------------------------------------------------- 
*
      eve = .true.
      good_xtrack = .true.
      good_ytrack = .true.
*
      call vzero(pos_mc,8)
      call vzero(dlength_mc,8)
      call vzero(layer_mc,8)
      call vzero(hitwire_mc,8) 
      call vzero(zhitwire_mc,8) 
      call vzero(ichan_mc,8)
      call vzero(iscintmc,2)
      call vzero(nxhdc_mc,4)
      call vzero(nyhdc_mc,4)
      call vzero(nxhit,4)
      call vzero(nyhit,4)
      call vzero(laymc,8)
      nhdc_mc = 8
      nhdctotmc = 0
      ax_mc=0
      bx_mc=0
      ay_mc=0
      by_mc=0

      cell = xc + wibeam
      hcell = cell/2
      st0z = stsbz + sb0z
      zscint = stsbz
*
      x = 2.
      muon = fluxmu/(fluxmu+fluxe)
      electron = fluxe/(fluxmu+fluxe)
*
      if (rndm(x).le.muon) then
        ipart = 1                  ! Muon
        de = -1
        prob = -1 
        call zenith(theta_mc,phi_mc,ipart)  
        call momentum(e,de,theta_mc,prob)     
      else
        ipart = 2                  ! Electron
        call zenith(theta_mc,phi_mc,ipart) 
        call momentum_e(e,theta_mc) 
      endif 
      energy_mc = e
      deltae = de  
*          
* ** Generate Impact Point on the top (x1,y1) and bottom (x2,y2)
* ** Scintillator Plane 
* ** itrig_mc = 1 General (ns=4), itrig_mc = 2 Central (ns=2)  

 


      xtmin = sdy + d0x + xvt  ! sdy = sdx
      xtmax = xtmin + 1000. 

      xbmin = sdy + d0x + xvb  ! sdy = sdx
      xbmax = xbmin + 1000.
 
      if (debug(8)) print *,'xmin ,xmax t-b ',xtmin,xtmax,xbmin,xbmax
      if (itrig_mc.eq.1.or.itrig_mc.eq.3) then
        ymin = sdy + d0y 
        nus = 4   
        ymax = ymin + nus*250.          
      elseif (itrig_mc.eq.2.or.itrig_mc.eq.4) then
        ymin = sdy + 250. + d0y 
        nus = 2
        ymax = ymin + nus*250.
      endif    
      if (debug(8)) print *,'ymin ,ymax ',ymin,ymax

      x1= xtmin + rndm(1)*(xtmax-xtmin)   
      y1= ymin + rndm(1)*(ymax-ymin) 
*              
      theta_mc =pi*theta_mc/180.
      phi_mc = pi*phi_mc/180.

      x2 = x1 + zscint*tan(theta_mc)*cos(phi_mc)
      y2 = y1 + zscint*tan(theta_mc)*sin(phi_mc)

      if (debug(8)) print *,'x1,y1  x2,y2',x1,y1,'   ',x2,y2      
 
      if (itrig_mc.eq.1.or.itrig_mc.eq.2) then
         itrf = 0 
         if ((x2.ge.xbmin .and. x2.le.xbmax) .and.
     >        (y2.ge.ymin .and. y2.le.ymax)) then 
            itrf = 1         
            iscintmc(1) = int((y1-sdy-d0y)/250)+5
            iscintmc(2) = int((y2-sdy-d0y)/250)+1
            if (debug(8)) print *,' scint top e bot ',iscintmc(1)
     >           ,iscintmc(2)
         endif
      endif

      if (itrig_mc.eq.3.or.itrig_mc.eq.4) then
         itrf = 0
         nscintt = int((y1-sdy-d0y)/250)
         nscintb = int((y2-sdy-d0y)/250)
         if (debug(8)) print *,'nscint t-b ',nscintt,nscintb 
         if ((x2.ge.xbmin .and. x2.le.xbmax).and.(nscintt.eq.
     >        nscintb).and.(y2.ge.ymin.and.y2.le.ymax)) then 
            itrf = 1
            iscintmc(1) = int((y1-sdy-d0y)/250)+5
            iscintmc(2) = int((y2-sdy-d0y)/250)+1
         endif
      endif

       
*  
* ** Energy cut for electrons
*                         
*      if (ipart.eq.2 .and. energy.lt.cutoff_etrig) itrf = 0
*
* ** If trigger fired perform track reconstruction
*



      if (itrf.eq.1) then
         nevsim_trigger = nevsim_trigger + 1 

         if (debug(8))print *,'evento n.  ',nevsim

*
***   X track *** 
*     
         ax_mc = (x1-x2)/(st0z-sb0z) 
         bx_mc = x1-ax_mc*st0z
         thetax_mc = atan(ax_mc)
         if (debug(8)) print *,'ax_mc, bx_mc, thetax_mc =',ax_mc,bx_mc,
     >        thetax_mc*180./3.1415
*
***   Y track ***
*
         ay_mc = (y1-y2)/(st0z-sb0z)
         by_mc = y1-ay_mc*st0z
         thetay_mc = atan(ay_mc)
         if (debug(8)) print *,'ay_mc, by_mc, thetay_mc =',ay_mc,by_mc,
     >        thetay_mc*180./3.1415
*
*     Determino la coordinata z di ognuno degli 8 piani dei fili
*
         
         do i=1,8
            layer_mc(i)= i
            zhitwire_mc(i) = zwire((i-1)*12)
            if (debug(8))   print *,i,' layer =',layer_mc(i),
     >           zhitwire_mc(i)
         enddo

*
*
*     Intersezione delle due rette con i rispettivi quattro piani
*

         do i=1,5,4
            pos_mc(i)=ax_mc*zhitwire_mc(i)+bx_mc !  global position x
            pos_mc(i+1)=ax_mc*zhitwire_mc(i+1)+bx_mc
            pos_mc(i+2)=ay_mc*zhitwire_mc(i+2)+by_mc !  global position y
            pos_mc(i+3)=ay_mc*zhitwire_mc(i+3)+by_mc
         enddo
         
         mom = .true.
         if (mom) then

            minbx = d0x+xvb+cdx
            minby = d0y+cdy
            mintx = d0x+xvt+cdx
            minty = d0y+cdy

            
         
            do i = 1,4
               if (i.le.2) then
                  if ((pos_mc(i).ge.(minbx+(46.75)*(i-1))).and.
     >                 (pos_mc(i).le.
     >                 (minbx+1122.+(46.75)*(i-1)))) nxhit(i)=2
                  if ((pos_mc(i+2).ge.(minby+(46.75)*(i-1))).and.
     >                 (pos_mc(i+2).le.
     >                 (minby+1122.+(46.75)*(i-1)))) nyhit(i)=2
               else
                  if ((pos_mc(i+2).ge.(mintx+(46.75)*(i-3))).and.
     >                 (pos_mc(i+2).le.
     >                 (mintx+1122.+(46.75)*(i-3)))) nxhit(i)=2
                  if ((pos_mc(i+4).ge.(minyt+(46.75)*(i-3))).and.
     >                 (pos_mc(i+4).le.
     >                 (minty+1122.+(46.75)*(i-3)))) nyhit(i)=2    
               endif
            enddo
** lungo y dovrebbe essere inutile poiche' la geometria e' fissata ***
        

            if (debug(8)) print *,'nxhit ',nxhit
            if (debug(8)) print *,'nyhit ',nyhit

            nxlaymc = 0
            nylaymc = 0

            do i=1,4
               if (i.le.2) then
                  if (nxhit(i).eq.2) then
                     nxlaymc = nxlaymc + 1
                     laymc(i) = 1

                  endif
               else
                  if (nxhit(i).eq.2) then
                     nxlaymc = nxlaymc + 1
                     laymc(i+2) = 1

                  endif
               endif
            enddo
        
            do i=1,4
               if (i.le.2) then
                  if (nyhit(i).eq.2) then
                     nylaymc = nylaymc + 1
                     laymc(i+2) = 1

                  endif
               else
                  if (nyhit(i).eq.2) then
                     nylaymc = nylaymc + 1
                     laymc(i+4) = 1

                  endif  
               endif
            enddo

            nhdctotmc = nxlaymc + nylaymc

         endif
        
         if (debug(8)) then
            print *,' nhdctmc ',nhdctotmc
            do i=1,8
               print *,'laymc ',laymc(i)
               if (laymc(i).eq.1) then
                  print *,'pos_mc(',i,')  ',pos_mc(i)
               else
                  print *,'layer ',i,'no mc hits'
               endif
            enddo
         endif

         do i = 3,7,4
            ichan_mc(i)=int((pos_mc(i)-d0y-cdy)/cell)+(i-1)*12 
            ichan_mc(i+1)=int((pos_mc(i+1)-hcell-d0y-cdy)/cell)+i*12
         enddo
     
         ichan_mc(1)=int((pos_mc(1)-d0x-cdx-xvb)/cell) 
         ichan_mc(2)=int((pos_mc(2)-d0x-cdx-xvb-hcell)/cell)+12 
         ichan_mc(5)=int((pos_mc(5)-d0x-cdx-xvt)/cell)+48
         ichan_mc(6)=int((pos_mc(6)-d0x-cdx-xvt-hcell)/cell)+60 
      
         if (debug(8)) then
            do i=1,8
               print *,'ichan_mc(',i,')  ',ichan_mc(i)
            enddo
         endif
      
         do i=1,8
            hitwire_mc(i) = wire(ichan_mc(i))
            dlength_mc(i)=abs(wire(ichan_mc(i))-pos_mc(i))          
         enddo
*
         print *,' '
         if (debug(8)) then
            do i=1,8
               print *,'layer,channel =',  layer_mc(i),ichan_mc(i)
               print *,' position, wire , z , dlength =',
     >              pos_mc(i),hitwire_mc(i),zhitwire_mc(i),dlength_mc(i)
            enddo   
         endif
      
         print *,' '

*
*
******  Smearing
*      
         do i =1,8
            if (laymc(i).eq.1) then
            
 400           continue
               CALL RANMAR(gn,1)
               CALL RNORML(gn,1)
*
               if (debug(8)) print *,'dlength+gn',i,dlength_mc(i)+gn 
               if (((dlength_mc(i)+gn).lt.0.).or.
     >              ((dlength_mc(i)+gn).gt.46.))  then
                  goto 400
               endif
               gauss_var(i) = gn      
               if (debug(8)) print *,'gauss  ',i, gauss_var(i)
            endif
         enddo
         
*      
         do i=1,8
            pos(i,1)=wire(ichan_mc(i))+dlength_mc(i)+gauss_var(i)
            pos(i,2)=wire(ichan_mc(i))-dlength_mc(i)-gauss_var(i)
            zpos(i,1)=zwire(ichan_mc(i))
            zpos(i,2)=zwire(ichan_mc(i))
            if(debug(8)) then
               print *,''
               print *,'hit + -,zp1,zp2'
     >              ,pos(i,1),pos(i,2),zpos(i,1),zpos(i,2)
            endif
         enddo
*      
         do i=1,4
            if (i.eq.1.or.i.eq.2) then
               do j=1,2
                  xpos(i,j) = pos(i,j) 
                  zxpos(i,j) = zpos(i,j)
               enddo
               do j=1,2
                  xpos(i+2,j) = pos(i+4,j)
                  zxpos(i+2,j) = zpos(i+4,j)                           
               enddo      
            elseif (i.eq.3.or.i.eq.4) then
               do j=1,2
                  ypos(i-2,j) = pos(i,j)
                  zypos(i-2,j) = zpos(i,j)
               enddo
               do j=1,2
                  ypos(i,j) = pos(i+4,j)
                  zypos(i,j) = zpos(i+4,j)
               enddo   
            endif
         enddo
*
      
         do i=1,4
            anx = float(nxhit(i))/2.
            any = float(nyhit(i))/2.
            nxhdc_mc(i)=int(anx)
            nyhdc_mc(i)=int(any)        
         enddo  
*     
         if (debug(8)) then
            print *,'    '
            do i=1,4
               print *,'layer, hit x, hit y =',i,nxhit(i),nyhit(i)
               do j=1,nxhit(i)
                  print *,'hit x, xpos , zxpos =',j,xpos(i,j),zxpos(i,j)
               enddo
               do j=1,nyhit(i)
                  print *,'hit y, ypos , zypos =',j,ypos(i,j),zypos(i,j)
               enddo
            enddo
         endif
      
*     
         do i = 1,4
            if (i.le.2) then
               xdlen (i,1) = dlength_mc(i)+gauss_var(i)
               ydlen (i,1) = dlength_mc(i+2)+gauss_var(i+2)
               xdlen (i,2) = dlength_mc(i)+gauss_var(i)
               ydlen (i,2) = dlength_mc(i+2)+gauss_var(i+2)
            else
               xdlen (i,1) = dlength_mc(i+2)+gauss_var(i+2)
               ydlen (i,1) = dlength_mc(i+4)+gauss_var(i+4) 
               xdlen (i,2) = dlength_mc(i+2)+gauss_var(i+2)
               ydlen (i,2) = dlength_mc(i+4)+gauss_var(i+4)
            endif
         enddo
         if (debug(8)) then
            print *,''
            do i = 1,4
               print *,'xdlen,  ydlen  ',xdlen(i,1),ydlen(i,1)
            enddo
         endif
*     
      endif
*
      RETURN
      END     

*
* *********************************************************
*
      SUBROUTINE RPCCLUSMC
*
* Create multiplicity for each spoken strip
*
* *********************************************************
*     

      INCLUDE "ntuple.inc"
      INCLUDE "dc.inc"
      PARAMETER (nxstrip = 16)
      PARAMETER (nystrip = 16)
      COMMON /rpcgeo/ stx(0:nxstrip-1),stzx(0:nxstrip-1),
     >     sty(0:nystrip-1),stzy(0:nystrip-1)     
    
*      DIMENSION nxstripmc(10),nystripmc(10)
      REAL ximprpcmc,yimprpcmc,dx,dy
*
* ---------------------------------------------------------
*

      multip_mc = 0
      nxstripcentralmc = 0
      nystripcentralmc = 0
      call vzero(nxstripmc,10)
      call vzero(nystripmc,10)




***    X plane   ****
*
      ximprpcmc = ax_good*340 + bx_good
      dx = ximprpcmc-342

      if (debug(8)) print*,'ximprpcmc ',ximprpcmc

      if ((dx.ge.0.).and.(dx.le.480)) then
         nxstripcentralmc = dx/30
        if (debug(8)) print*,'central strip',nxstripcentralmc 
      else
         nxstripcentralmc = -1
         if (debug(8)) print*,'central strip',nxstripcentralmc  
         return
      endif
      
      if (ximprpcmc.ge.stx(nxstripcentralmc)) then
         kk = 1
      else
         kk = 0
      endif
      if (debug(8)) print *,'xstrip kk',stx(nxstripcentralmc),kk

      CALL multiplicity(multip_mc)

      if (debug(8)) print *,'multip',multip_mc,kk

      do i = 1,multip_mc
         if (((-1)**multip_mc).lt.0) then   
            nxstripmc(i) = nxstripcentralmc - int(multip_mc/2) + (i-1)
         else
            nxstripmc(i) = nxstripcentralmc - int(multip_mc/2) 
     >       + (i-1) + kk  
         endif
      enddo
      
      if (debug(8)) then
         do i= 1,multip_mc
            print *,'strip non corrette',nxstripmc(i)
         enddo
      endif
      
      if (nxstripmc(1).lt.0) then
         ishift = abs(nxstripmc(1))
         do i = 1,multip_mc
           if (debug(8)) print *,'lstrip  i prima', nxstripmc(i),'  ',i 
            nxstripmc(i) = nxstripmc(i) + ishift
          if (debug(8)) print *,'lstrip  i dopo', nxstripmc(i),'  ',i
         enddo
      endif
      if (nxstripmc(multip_mc).gt.15) then
         diffstrip = nxstripmc(multip_mc)-15
         do i = 1,multip_mc
            if (debug(8)) print *,'gstrip  i prima', nxstripmc(i),'  ',i 
            nxstripmc(i) = nxstripmc(i) - diffstrip
            if (debug(8)) print *,'gstrip  i dopo', nxstripmc(i),'  ',i 
         enddo
      endif

      if (debug(8)) then
         do i= 1,multip_mc
            print *,'strip corrette',nxstripmc(i)
         enddo
      endif



      RETURN
      END




*
* *********************************************************
*
       SUBROUTINE zenith(th,phi,ipart)
*
*  Generate angles theta and phi
* 
* *********************************************************
*
       PARAMETER (pi = 3.1415926)               
       REAL buffer1,buffer2
       REAL th,phi           
*
* ----------------------------------------------------------
* 
       if (ipart.eq.1) then 
         nn = 2.   
       elseif (ipart.eq.2) then
         nn = 3.6
       endif
*
       call generate(buffer1,buffer2,nn)
       phi = 360*buffer2
       th = buffer1*90.
*
       END
*
* ************************************************************
*
       SUBROUTINE generate(buffer1,buffer2,nn)
*                                
       PARAMETER (pi = 3.1415926)
       REAL buffer1,buffer2,x,hitmiss
       REAL teta,cos2teta         
*
* --------------------------------------------------------
*                            
       x = 1
       flag = 0
       buffer2=rndm(x)
*
*      hit or miss 
*                
       do while (flag.eq.0)           
         buffer1 = rndm(x) 
         hitmiss = rndm(x)
         teta = pi*buffer1/2.
         cos2teta = cos(teta)**nn 
         cos2teta = cos2teta*sin(teta)
         if (hitmiss.lt.cos2teta) then
           flag=1
         endif
       enddo       
*
       END                        
*
* **********************************************************
*
       SUBROUTINE momentum(e,de,teta,prob)
*
* ******************************************************
*
       parameter (nev = 10000)
       parameter (pi = 3.1415926)
       parameter (ro = 0.0016)
       parameter (alpha = 0.00206)
       parameter (h0 =630.)                   
*
       parameter (A = 0.199)
       parameter (gamma = 2.63)
       parameter (b1 = 90.)
       parameter (r1 = .76) 
       parameter (rmumass = 0.1)
       parameter (emax = 1000.)
       parameter (emin = 0.22)   
       parameter (rlambdaf = 1030.)
       parameter (rlambda0 = 120.)
       parameter (taumu = 2.197e-6)
       parameter (clight = 3.e10)                           
*
       real e,teta,de,x,random1,random2,random3,prob
       real C,H,f1,f2,gemin,ge,rlambdad,emean,xx
*
* ----------------------------------------------------------
* 
*    calcola energy loss according de/dx = ro * alpha * x (x (cm),  e (GeV))
*       
       de = -1                                                          
       e = -1
       if (cos(teta).ne.0) then      
            xx = (rlambdaf/cos(teta*pi/180.)-rlambda0)/ro
            de = ro * alpha * xx
       endif

       if (de.eq.-1) then
            return
       endif

       x = 1.
       flag = 0 
       flag1 = 0                               
       C = A*r1**(gamma-1.)
       RM = (emax+de)**(1.-gamma) - (emin+de)**(1.-gamma)
       H = (1.-gamma)/(C*RM)
       do while (flag.eq.0)
           random1 = rndm(x)                                                  
           f1 = (1.-gamma)*random1/(H*C)  
           f2 = (emin+de)**(1.-gamma)
           e = (f1 + f2)**(1./(1.-gamma)) - de
           emean = e + de/2.
*
*    first hit or miss  (deacy probability of muon)
*                                                                  
           ter = taumu*emean*ro       
           rlambdad = ter*cos(pi*teta/180.)*clight/rmumass  
           if (cos(pi*teta/180.).eq.0) then
                prob = 1.
           else
                prob = exp(-(rlambdaf-rlambda0)/rlambdad)
           endif
           random3 = rndm(x)
           if (random3.lt.prob) then
              flag1 = 1
           endif
           if (flag1.eq.1) then
*
*     second hit or miss
*               
              gg = C*b1*RM/(1.-gamma)  
              gemin = gg * (1./(emin+de+b1))
              ge = gg * (1./(e+de+b1))
              random2 = rndm(x) * gemin                    
              if (random2.le.ge) then
                   flag = 1 
              endif
           endif
           flag1 = 0   
       enddo
*
       RETURN
       END           
*
* **************************************************************
*
       REAL FUNCTION rintegral(ene)
*      
* **************************************************************
* 
       parameter (A = 0.199)
       parameter (gamma = 2.63)
       parameter (b1 = 90.)   
       parameter (ro = 0.0016)
       parameter (r1 = .76) 
       parameter (demean = 2.4)    
       parameter (rlambdaf = 1030.)
       parameter (rlambda0 = 120.)
       parameter (taumu = 2.197e-6)
       parameter (clight = 3.e10)                           
       parameter (rmumass = 0.1)
       parameter (costeta = .85)
*
* ------------------------------------------------------------
*            
*       const = taumu*ro*clight/rmumass       
*       delambda = rlambdaf - rlambda0         
       rintegral= 
     +  (A*r1**(gamma-1)*(ene+demean)**(-gamma)*b1/(ene+demean+b1))
     +      *exp(-((rlambdaf-rlambda0)/((taumu*ro*clight/rmumass)*
     +       (ene+demean/2.)*costeta)))
*
       END
*
* ****************************************************************
*
       SUBROUTINE  momentum_e(e_e)
*
* *******************************************************************
*
       parameter (cutoff_ene = 0.04)
       parameter (ene_max = 1000.)
       parameter (a = 1.45)       
*
* --------------------------------------------------------------------
*
       d_el = (cutoff_ene**(-a) - ene_max**(-a))
       random = rndm(x)          
       e_e = ((1-random)*d_el + ene_max**(-a))**(-1/a)
*      
       END 
*
* *************************************************************************
* 
       SUBROUTINE flux_electron(summ_e)
*
*   Ccalculate integral flux of electron depending from value of 
*   cut off energy (minimum electron energy accepted).
*                   
* *************************************************************************
*
       parameter (cutoff_ene = 0.04)
*
       summ_e = .000022*cutoff_ene**(-1.45)
*
       END
*
* ************************************************************************
*
*
      SUBROUTINE multiplicity(imultip)
      
      DIMENSION ientries(8)
      DATA ientries
     > /5400,12500,20500,2800,1400,400,1000,400/

*      if (debug(8)) then 
         do i = 1,8
            print *,'entries',ientries(i)
         enddo
*      endif



*      DATA ientries(1) /5400/, ientries(2) /12500/,  ientries(3) /20500/,
*     >     ientries(4) /2800/, ientries(5) /1400/,  ientries(6) /400/,
*     >     ientries(7) /1000/, ientries(8) /400/

      

      flag = 0
      do while (flag.eq.0)
         imultip = rndm(1)*8
         ivalue = rndm(1)*20500
         if (ivalue.lt.ientries(imultip)) then
            flag =1
         endif
      enddo
      print *,'imultip',imultip


      RETURN
      END
*
***************************************************************************
*
      SUBROUTINE TRACKONY
*
*     check if track Y is on fired scintillators 
*
* *************************************************************************
*
      INCLUDE "dc.inc"
      INCLUDE "ntuple.inc"
      PARAMETER (nscintil=8)
      PARAMETER (ntrackm=150)
      COMMON /scintgeo/  yscint(nscintil),zscint(nscintil)
      COMMON /tracks/ ntrack,nlh(4),pos_track(ntrackm,4),
     >                fit_track(ntrackm,4),res_track(ntrackm,4),
     >                chi2_track(ntrackm),theta_track(ntrackm),
     >                a_track(ntrackm),b_track(ntrackm),
     >                nfit_track(ntrackm),ihhist_track(ntrackm,4),      
     >                itrackflag(ntrackm)     
      

      
      timpact = a_track(ntrack)*3245 + b_track(ntrack)
      bimpact = a_track(ntrack)*(0.) + b_track(ntrack)
      

      itracktop = 0
      itrackbot = 0

      do i=1,nscint
         if (debug(4)) print *,'# of scint ',nscint
         if (iscint(i).le.4) then
            if (debug(4)) then
               print *,'bottom layer- scint number ',iscint(i)
               print *,'scintpos ',yscint(iscint(i))
            endif

            if ((bimpact.gt.(yscint(iscint(i))-145)).and.
     >           (bimpact.lt.(yscint(iscint(i))+145))) then
               itrackbot = 1
            endif
         else
            if (debug(4)) then
               print *,'top layer- scint number ',iscint(i)
               print *,'scintpos ',yscint(iscint(i))
            endif

            if ((timpact.gt.(yscint(iscint(i))-145)).and.
     >           (timpact.lt.(yscint(iscint(i))+145))) then
               itracktop = 1
            endif
         endif
         if (debug(4)) print *,'itrack bot-top ',itrackbot,' ',itracktop
      enddo
      if (itrackbot.eq.1.and.itracktop.eq.1) then
         itrackf = 1
      else
         itrackf = 0
      endif
      itrackflag(ntrack) = itrackf
      if (debug(4)) print *,'ntrack timp bimp flag '
     >     ,ntrack,timpact,bimpact,itrackf
      
      RETURN
      END


*
***************************************************************************
*
      SUBROUTINE TRACKONX
*
*     check if track X is on scintillators surface
*
* *************************************************************************
*

      INCLUDE "dc.inc"
      INCLUDE "ntuple.inc"
      PARAMETER (ntrackm=150)
      COMMON /tracks/ ntrack,nlh(4),pos_track(ntrackm,4),
     >                fit_track(ntrackm,4),res_track(ntrackm,4),
     >                chi2_track(ntrackm),theta_track(ntrackm),
     >                a_track(ntrackm),b_track(ntrackm),
     >                nfit_track(ntrackm),ihhist_track(ntrackm,4),      
     >                itrackflag(ntrackm)     
            
      timpact = a_track(ntrack)*3245 + b_track(ntrack)
      bimpact = a_track(ntrack)*(0.) +b_track(ntrack)

      xtopleft = 170. + xvt - 200.     
      xtopright = 1170. + xvt + 200.  
      xbotleft = 170. + xvb - 200.    
      xbotright =1170. + xvb + 200.   


      if (((timpact.gt.xtopleft).and.(timpact.lt.xtopright)).and.
     >   ((bimpact.gt.xbotleft).and.(bimpact.lt.xbotright))) then
         
         itrackf = 1
      else
         itrackf = 0
      endif
      itrackflag(ntrack) = itrackf
      if (debug(4)) print *,'ntrack timp bimp flag '
     >     ,ntrack,timpact,bimpact,itrackf
      
      RETURN
      END

*
* ********************************************************************
*

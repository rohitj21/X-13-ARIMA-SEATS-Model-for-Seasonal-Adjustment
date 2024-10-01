C     Last change:  BCM  25 Nov 1998   12:19 pm
**==x12hdr.f    processed by SPAG 4.03F  at 10:39 on 20 Oct 1994
      SUBROUTINE x12hdr(Nfcst,Srsttl,Nsrscr,Ttlvec,Notc,Lx11,Lmodel,
     &                  Lseats,Lwidpr,Pos1,Nuspad,Nustad,Iqtype,Fcntyp,
     &                  Lam,Ciprob,Dattim,Cnstnt,Isrflw,Lognrm)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     Print header page and title info for X-13ARIMA-SEATS
c-----------------------------------------------------------------------
      LOGICAL F
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO=0D0,F=.false.)
c-----------------------------------------------------------------------
      INCLUDE 'stdio.i'
      INCLUDE 'lex.i'
      INCLUDE 'srslen.prm'
c-----------------------------------------------------------------------
      INCLUDE 'build.prm'
      INCLUDE 'metadata.prm'
      INCLUDE 'model.prm'
      INCLUDE 'notset.prm'
      INCLUDE 'prior.prm'
      INCLUDE 'tbllog.prm'
c-----------------------------------------------------------------------
      INCLUDE 'agr.cmn'
      INCLUDE 'error.cmn'
      INCLUDE 'force.cmn'
      INCLUDE 'hiddn.cmn'
      INCLUDE 'mdldat.cmn'
      INCLUDE 'metadata.cmn'
      INCLUDE 'missng.cmn'
      INCLUDE 'mq3.cmn'
      INCLUDE 'prior.cmn'
      INCLUDE 'rho.cmn'
      INCLUDE 'tbllog.cmn'
      INCLUDE 'title.cmn'
      INCLUDE 'units.cmn'
      INCLUDE 'x11adj.cmn'
      INCLUDE 'x11log.cmn'
      INCLUDE 'x11msc.cmn'
      INCLUDE 'x11opt.cmn'
      INCLUDE 'x11reg.cmn'
c-----------------------------------------------------------------------
      INCLUDE 'mdltbl.i'
      INCLUDE 'cmptbl.i'
      INCLUDE 'x11tbl.i'
      INCLUDE 'spctbl.i'
c-----------------------------------------------------------------------
      LOGICAL Lx11,Lmodel,Lseats,Lwidpr,Lognrm
      DOUBLE PRECISION Lam,Ciprob,Cnstnt
      INTEGER Nfcst,i,ie2,ie3,ie4,ie5,ie6,ie7,Notc,j,nttl,ipos,
     &        Nsrscr,nstr,ie7a,Pos1,lnlen,msp,linsp,nxrg,nastr,
     &        nxrg2,Nuspad,Nustad,Fcntyp,ikey,Iqtype,Isrflw,ival
      CHARACTER qqmm*(9),runs*(26),num*(2),mmqq*(7),malo*(14),avg*(8),
     &          avg2*(4),xhdr*(LINLEN),Ttlvec*(80),cmonth*(3),cqtr*(3),
     &          Dattim*(24),mqcds*(3),r*(11),pcd*(15),sf*(4),cm2*(6),
     &          xb*(LINLEN),str*(10),Srsttl*(*),xrgstr*(50),adjstr*(20),
     &          outstr*(25),thisky*(LINLEN),thisvl*(LINLEN)
      DIMENSION qqmm(3),runs(4),num(4),mmqq(3),malo(5),
     &          avg(7),avg2(8),Ttlvec(10),cmonth(12),cqtr(4),mqcds(3),
     &          r(2),pcd(2),sf(PSP),cm2(12),Pos1(2)
c-----------------------------------------------------------------------
      INTEGER nblank
      LOGICAL dpeq
      EXTERNAL nblank,dpeq
c-----------------------------------------------------------------------
      DATA pcd/'percent change ','differences    '/
      DATA qqmm/'monthly  ','quarterly','         '/
      DATA runs/
     &   'seasonal adjustment       ','summary measures          ',
     &   'trend estimation          ','regARIMA model estimation '/
      DATA num/'st','nd','rd','th'/
      DATA mmqq/'month  ','quarter','period '/
      DATA malo/'multiplicative','additive      ','logarithmic   ',
     &          'pseudo-add    ','auto-mode     '/
      DATA avg/'3x3     ','3x5     ','3x9     ','3x15    ','Stable  ',
     &         'MSR     ','3x1     '/
      DATA avg2/'def ','3x3 ','3x5 ','3x9 ','3x15','Stbl','MSR ','3x1 '/
      DATA cmonth/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &     'Oct','Nov','Dec'/
      DATA cm2/'uary  ','ruary ','ch    ','il    ','      ','e     ',
     &     'y     ','ust   ','tember','ober  ','ember ','ember '/
      DATA cqtr/'1st','2nd','3rd','4th'/
      DATA mqcds/'mcd','qcd','   '/
      DATA r/'ratios     ','differences'/
c-----------------------------------------------------------------------
      IF(Ny.eq.1)THEN
       ie2=0
      ELSE IF(Ny.eq.12)THEN
       ie2=1
      ELSE IF(Ny.eq.4)THEN
       ie2=2
      ELSE
       ie2=3
      END IF
      IF(Lx11)THEN
       ie3=Kfulsm+1
      ELSE IF(Lseats)THEN
       ie3=1
      ELSE
       ie3=4
      END IF
      ie4=0
      ie5=Pos1(MO)
      IF(ie5.gt.4)ie5=4
      ie6=Lstmo
      IF(ie6.gt.4)ie6=4
      ie7=Muladd+1
      IF(Psuadd)ie7=4
      IF(Fcntyp.eq.0)ie7=5
c-----------------------------------------------------------------------
      IF(ie2.gt.0)THEN
       Qm=qqmm(ie2)
       Moqu=mmqq(ie2)
       Mqcd=mqcds(ie2)
      ELSE
       Qm=qqmm(3)
       Moqu='       '
       Mqcd=mqcds(3)
      END IF
      IF(ie7.eq.5)THEN
       Pcdif=pcd(2)
       Rad=r(2)
      ELSE
       ie7a=2-mod(ie7,2)
       IF(Psuadd)ie7a=2
       Pcdif=pcd(ie7a)
       Rad=r(ie7a)
      END IF
      IF(Kfulsm.lt.2.and.Lx11)THEN
       sf(1)=avg2(Lterm+1)
       DO i=2,Ny
        IF(Lterm.gt.0.and.Lter(i).eq.0)THEN
         ie4=ie4+1
        ELSE IF(Lter(i).ne.Lterm)THEN
         ie4=Lter(i)+ie4
        END IF
        sf(i)=avg2(Lter(i)+1)
       END DO
      END IF
c-----------------------------------------------------------------------
      CALL setchr(' ',LINLEN,xb)
      nttl=nblank(Ttlvec(1))
      Ntitle=0
      IF(Nsrscr.gt.0)THEN
       IF(xb(1:Nsrscr).eq.Srsttl(1:Nsrscr))THEN
        Nsrscr=0
       ELSE
        Title=Srsttl(1:Nsrscr)//xb(1:(80-Nsrscr))
        Ntitle=Nsrscr
       END IF
      END IF
      IF(Ntitle.eq.0)THEN
       IF(nttl.gt.0)THEN
        Title=Ttlvec(1)(1:nttl)//xb(1:(80-nttl))
        Ntitle=nttl
       ELSE IF(Nser.gt.0)THEN
        Ntitle=8+Nser+NPRGNM
        IF(Ntitle.gt.79)THEN
         Title=PRGNAM//' run of '//Serno(1:Nser-(Ntitle-80))
         Ntitle=80
        ELSE
         Title=PRGNAM//' run of '//Serno(1:Nser)//xb(1:(80-Ntitle))
        END IF
       ELSE
        Ntitle=NPRGNM
        Title=PRGNAM//xb(1:(80-Ntitle))
       END IF
      END IF
c-----------------------------------------------------------------------
      IF(Savtab(LSRSHD))THEN
       CALL svdttm(Nform,Dattim)
       WRITE(Nform,1790)VERNUM,BUILD
 1790  FORMAT('version: ',a,/,'build: ',a,/,'output: out')
      END IF
      IF(Prttab(LSRSHD))THEN
c-----------------------------------------------------------------------
       lnlen=80
       msp=7
       linsp=8
       IF(Lwidpr)THEN
        lnlen=LINLEN
        msp=msp+24
        linsp=13
       END IF
c-----------------------------------------------------------------------
       i=(lnlen-49)/2
       WRITE(xhdr,1010)xb(1:i)
 1010  FORMAT(a,'U. S. Department of Commerce, U. S. Census Bureau')
       IF(.not.Lcmpaq)WRITE(Mt1,1020)'1'
       WRITE(Mt1,1020)xhdr(1:nblank(xhdr))
 1020  FORMAT(/,a)
c-----------------------------------------------------------------------
       CALL setchr(' ',LINLEN,xhdr)
       IF(Ny.eq.12.or.Ny.eq.4)THEN
        i=(lnlen-(NPRGNM+nblank(Qm)+nblank(runs(ie3))+10))/2
        WRITE(xhdr,1030)xb(1:i),PRGNAM,Qm(1:nblank(Qm)),
     &        runs(ie3)(1:nblank(runs(ie3)))
       ELSE
        i=(lnlen-(NPRGNM+nblank(runs(ie3))+10))/2
        WRITE(xhdr,1029)xb(1:i),PRGNAM,runs(ie3)(1:nblank(runs(ie3)))
       END IF
 1029  FORMAT(a,a,1x,a,' Method,')
 1030  FORMAT(a,a,1x,a,1x,a,' Method,')
       WRITE(Mt1,1020)xhdr(1:nblank(xhdr))
       i=(lnlen-(23+nblank(VERNUM)+nblank(BUILD)))/2
       WRITE(Mt1,1031)xb(1:i),VERNUM,BUILD
 1031  FORMAT(a,'Release Version ',a,' Build ',a)
       IF(.not.Lcmpaq)WRITE(Mt1,'(/)')
       WRITE(Mt1,1040)xb(1:(msp+4)),xb(1:(msp+6)),xb(1:(msp+10)),
     &                xb(1:(msp+15)),xb(1:(msp+5)),xb(1:msp),
     &                xb(1:msp),xb(1:(msp+9)),xb(1:(msp+8)),
     &                xb(1:(msp+3)),xb(1:(msp+2)),PRGNAM,xb(1:(msp+3)),
     &                xb(1:(msp+5)),xb(1:(msp+7)),xb(1:(msp+10)),
     &                xb(1:(msp+15))
 1040  FORMAT(a,'This software application provides an enhanced ',
     &          'version of',/,
     &        a,'Statistics Canada''s X-11-ARIMA extension (Dagum, ',
     &          '1980)',/,
     &        a,'of the X-11 variant of the Census Method II of',/,
     &        a,'Shiskin, Young and Musgrave (1967).',//
     &        a,'It also provides an ARIMA model-based method ',
     &          'following',/,
     &        a,'Hillmer and Tiao (1982) and Burman (1980) that is ',
     &          'very similar',/,
     &        a,'to the update of the method of SEATS (Gomez and ',
     &          'Maravall, 1996)',/,
     &        a,'produced at the Bank of Spain by G. Caporello and',/,
     &        a,'A. Maravall for TSW (Caporello and Maravall, 2004).',/,
     &        a,'The present application includes additional ',
     &          'enhancements.',//,
     &        a,a,' includes an automatic ARIMA model selection ',
     &            'procedure',/,
     &        a,'based largely on the procedure of Gomez and Maravall ',
     &          '(1998)',/,
     &        a,'as implemented in TRAMO (1996) and subsequent ',
     &          'revisions.',//,
     &        a,'Primary Programmers: Brian Monsell, Mark Otto and,',/,
     &        a,'for the ARIMA model-based signal extraction,',/,
     &        a,'Gianluca Caporello and Victor Gomez',//)
c-----------------------------------------------------------------------
       WRITE(Mt1,1050)Title(1:nblank(Title))
       IF(Nser.gt.0)WRITE(Mt1,1060)Serno(1:Nser)
       IF(Dattim(1:2).ne.'  ')WRITE(Mt1,1070)Dattim(2:24)
 1050  FORMAT(5x,'Series Title- ',a)
 1060  FORMAT(5x,'Series Name- ',a)
 1070  FORMAT(5x,a)
c-----------------------------------------------------------------------
       IF(Notc.gt.0)THEN
        DO i=1,Notc
         WRITE(Mt1,1080)Ttlvec(i)
 1080    FORMAT(24x,a)
        END DO
       END IF
c-----------------------------------------------------------------------
       IF(Ny.eq.4.or.Ny.eq.12)THEN
        WRITE(Mt1,1090)xb(1:linsp),Pos1(MO),num(ie5),
     &                 Moqu(1:nblank(Moqu)),Pos1(YR),Lstmo,num(ie6),
     &                 Moqu(1:nblank(Moqu)),Lstyr
 1090   FORMAT(/,a,'-Period covered- ',i2,a2,1x,a,',',i4,' to ',i2,a2,
     &         1x,a,',',i4)
       ELSE IF (Ny.eq.1) THEN
        WRITE(Mt1,1091)xb(1:linsp),Pos1(YR),Lstyr
 1091   FORMAT(/,a,'-Period covered- ',i4,' to ',i4)
       ELSE
        WRITE(Mt1,1092)xb(1:linsp),Pos1(YR),Pos1(MO),Lstyr,Lstmo
 1092   FORMAT(/,a,'-Period covered- ',i4,'.',i2.2,' to ',i4,'.',i2.2)
       END IF
c-----------------------------------------------------------------------
       IF(Lx11)THEN
        WRITE(Mt1,1100)xb(1:linsp),malo(ie7)(1:nblank(malo(ie7))),
     &                 runs(ie3)(1:nblank(runs(ie3)))
 1100   FORMAT(a,'-Type of run - ',a,1x,a)
        IF(.not.Lcmpaq)WRITE(Mt1,'()')
c-----------------------------------------------------------------------
        WRITE(Mt1,1110)xb(1:linsp),Sigml,Sigmu
 1110   FORMAT(a,'-Sigma limits for graduating extreme values are ',
     &         f4.1,' and ',f4.1,' .')
c-----------------------------------------------------------------------
c        IF(Kexopt.eq.1)WRITE(Mt1,1120)xb(1:linsp)
c 1120   FORMAT(a,'-Modify extreme values before computing the B7 ',
c     &         'trend cycle curve.')
c-----------------------------------------------------------------------
        IF(Kfulsm.lt.2)THEN
c-----------------------------------------------------------------------
         IF(ie4.ne.0)THEN
          WRITE(Mt1,1150)xb(1:linsp)
 1150     FORMAT(a,'-The following moving averages were selected for ',
     &             'the seasonal factor curves:')
          linsp=linsp+5
          IF(Ny.eq.4)WRITE(Mt1,1160)xb(1:linsp),(cqtr(i),i=1,4),
     &                              xb(1:linsp),(sf(j),j=1,4)
          IF(Ny.eq.12)WRITE(Mt1,1170)xb(1:linsp),(cmonth(i),i=1,12),
     &                               xb(1:linsp),(sf(j),j=1,12)
          linsp=linsp-5
 1160     FORMAT(a,4(1x,a3,1x),/,a,4(1x,a4))
 1170     FORMAT(a,12(1x,a3,1x),/,a,12(1x,a4))
         ELSE IF(Lterm.eq.0)THEN
          WRITE(Mt1,1180)xb(1:linsp),xb(1:(linsp+1))
 1180     FORMAT(a,'-3x3 moving average used in section 1 of each ',
     &           'iteration,',/,a,'3x5 moving average in section 2 of ',
     &           'each iteration.')
         ELSE IF(Lterm.eq.6)THEN
          WRITE(Mt1,1190)xb(1:linsp),xb(1:(linsp+1)),xb(1:(linsp+1))
 1190     FORMAT(a,'-3x3 moving average used in section 1 of each ',
     &           'iteration, ',/,a,'3x5 moving average in section 2',
     &           ' of iterations B and C,'/,a,
     &           'moving average for final seasonal factors chosen by ',
     &           'Global MSR.')
         ELSE
          WRITE(Mt1,1200)xb(1:linsp),avg(Lterm)(1:nblank(avg(Lterm)))
 1200     FORMAT(a,'-a ',a,' moving average selected for the ',
     &           'seasonal factor curves.')
         END IF
        END IF
c-----------------------------------------------------------------------
        IF(Ktcopt.gt.0)THEN
         WRITE(Mt1,1210)xb(1:linsp),Ktcopt
 1210    FORMAT(a,'-Moving average for the variable trend cycle ',
     &          'routine is a ',i3,'-term Henderson')
        END IF
       END IF
c-----------------------------------------------------------------------
       IF(Kfmt.gt.0)THEN
        CALL setchr(' ',20,adjstr)
        IF(Ny.eq.12.or.Ny.eq.4)THEN
         WRITE(adjstr,1219)Qm(1:nblank(Qm))
 1219    FORMAT(a,' adjustment')
         nastr=18
         IF(Ny.eq.4)nastr=nastr+2
        ELSE
         adjstr(1:10)='adjustment'
         nastr=10
        END IF
        IF(Nustad.gt.0.and.Nuspad.gt.0)THEN
         WRITE(Mt1,1220)xb(1:linsp),adjstr(1:nastr),
     &                  '(temporary and permanent)'
        ELSE IF(Nustad.gt.0)THEN
         WRITE(Mt1,1220)xb(1:linsp),adjstr(1:nastr),'(temporary)'
        ELSE IF(Nuspad.gt.0)THEN
         WRITE(Mt1,1220)xb(1:linsp),adjstr(1:nastr),'(permanent)'
        END IF
 1220   FORMAT(a,'-Prior ',a,' ',a,' factors.')
       END IF
       IF(.not.dpeq(Cnstnt,DNOTST))THEN
        ipos=1
        CALL dtoc(Cnstnt,outstr,ipos)
        WRITE(Mt1,1221)xb(1:linsp),outstr(1:(ipos-1))
 1221   FORMAT(a,'-Constant value added to series (',a,')')
       END IF
c-----------------------------------------------------------------------
       IF(Lx11)THEN
        IF(Kswv.eq.1)THEN
         WRITE(Mt1,1230)xb(1:linsp)
 1230    FORMAT(a,'-Prior trading day adjustment.')
c         WRITE(Mt1,1230)xb(1:linsp),Out(1:(nblank(Out)+1)),
c     &                  mmqq(1)(1:nblank(mmqq(1)))
c 1230    FORMAT(a,'-Prior trading day adjustment with',a,'length of ',
c     &          a,' adjustment.')
        END IF
c-----------------------------------------------------------------------
        IF(Ixreg.gt.0)THEN
c-----------------------------------------------------------------------
         CALL wrtdat(Begxrg,Ny,str,nstr)
         nxrg=0
         IF(Tdgrp.gt.0.and.Holgrp.gt.0)THEN
          xrgstr='Trading day and holiday'
          nxrg=23
         ELSE IF(Tdgrp.gt.0)THEN
          xrgstr='Trading day'
          nxrg=11
         ELSE IF(Stdgrp.gt.0)THEN
          xrgstr='Stock trading day'
          nxrg=17
         ELSE IF(Holgrp.gt.0)THEN
          xrgstr='Holiday'
          nxrg=7
         END IF
         IF(nxrg.gt.0)THEN
          IF(Otlxrg)THEN
           IF(Critxr.gt.0)THEN
            IF(Lwidpr)THEN
             WRITE(Mt1,1240)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr),
     &                      xb(1:linsp),Critxr
 1240        FORMAT(a,'-',a,' irregular regression computed ',
     &              'starting ',a,' with AO outliers identified using',
     &              ' a',/,a,'critical value of ',f5.2,'.')
            ELSE
             WRITE(Mt1,1241)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr),
     &                      xb(1:linsp),Critxr
 1241        FORMAT(a,'-',a,' irregular regression computed ',
     &              'starting ',a,/,a,' with AO outliers identified ',
     &              'using a critical value of ',f5.2,'.')
            END IF
           ELSE
            IF(Lwidpr)THEN
             WRITE(Mt1,1242)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr),
     &                      xb(1:linsp)
 1242        FORMAT(a,'-',a,' irregular regression computed ',
     &              'starting ',a,' with AO outliers identified using',
     &              ' a',/,a,'default critical value.')
            ELSE
             WRITE(Mt1,1243)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr),
     &                      xb(1:linsp)
 1243        FORMAT(a,'-',a,' irregular regression computed ',
     &              'starting ',a,/,a,' with AO outliers identified ',
     &              'using a default critical value.')
            END IF
           END IF
          ELSE IF(Sigxrg.gt.ZERO)THEN
           IF(Lwidpr)THEN
            WRITE(Mt1,1244)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr),Sigxrg
 1244       FORMAT(a,'-',a,' irregular regression computed ',
     &             'starting ',a,' excluding irregular values outside ',
     &             f5.2,'-sigma limits.')
           ELSE
            WRITE(Mt1,1245)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr),
     &                     xb(1:linsp),Sigxrg
 1245       FORMAT(a,'-',a,' irregular regression computed ',
     &             'starting ',a,/,a,' excluding irregular values ',
     &             'outside ',f5.2,'-sigma limits.')
           END IF
          ELSE
           WRITE(Mt1,1246)xb(1:linsp),xrgstr(1:nxrg),str(1:nstr)
 1246      FORMAT(a,'-',a,' irregular regression computed ',
     &             'starting ',a)
          END IF
         END IF
c-----------------------------------------------------------------------
         IF(Tdgrp.gt.0)THEN
          IF(Axrgtd)THEN
           IF(Ixreg.eq.2)THEN
            WRITE(Mt1,1250)xb(1:linsp),'Trading day',' ',
     &                                 ' as prior factors.'
           ELSE
            WRITE(Mt1,1250)xb(1:linsp),'Trading day',' ','.'
           END IF
          ELSE
           WRITE(Mt1,1250)xb(1:linsp),'Trading day',' not','.'
          END IF
         END IF
         IF(Stdgrp.gt.0)THEN
          IF(Axrgtd)THEN
           IF(Ixreg.eq.2)THEN
            WRITE(Mt1,1250)xb(1:linsp),'Stock trading day',' ',
     &                                 ' as prior factors.'
           ELSE
            WRITE(Mt1,1250)xb(1:linsp),'Stock Trading day',' ','.'
           END IF
          ELSE
           WRITE(Mt1,1250)xb(1:linsp),'Stock Trading day',' not','.'
          END IF
         END IF
         IF(Holgrp.gt.0)THEN
          IF(Axrghl)THEN
           IF(Ixreg.eq.2)THEN
            WRITE(Mt1,1250)xb(1:linsp),'Holiday',' ',
     &                                 ' as prior factors.'
           ELSE
            WRITE(Mt1,1250)xb(1:linsp),'Holiday',' ','.'
           END IF
          ELSE
           WRITE(Mt1,1250)xb(1:linsp),'Holiday',' not','.'
          END IF
         END IF
 1250    FORMAT(a,'-',a,' irregular regression estimates',a,'applied',a)
c-----------------------------------------------------------------------
         nxrg=0
         IF(Xeastr.or.Xuser.or.Xtdtst.gt.0)THEN
          IF(Xtdtst.eq.2)THEN
           xrgstr((nxrg+1):(nxrg+17))='stock trading day'
           nxrg=nxrg+17
          ELSE IF(Xtdtst.gt.0)THEN
           xrgstr((nxrg+1):(nxrg+11))='trading day'
           nxrg=nxrg+11
          END IF
          IF(Xeastr)THEN
           IF(nxrg.gt.0)THEN
            xrgstr((nxrg+1):(nxrg+1))=','
            nxrg=nxrg+1
           END IF
           xrgstr((nxrg+1):(nxrg+6))='Easter'
           nxrg=nxrg+6
          END IF
          IF(Xuser)THEN
           IF(nxrg.gt.0)THEN
            xrgstr((nxrg+1):(nxrg+1))=','
            nxrg=nxrg+1
           END IF
           xrgstr((nxrg+1):(nxrg+12))='user-defined'
           nxrg=nxrg+12
          END IF
          WRITE(Mt1,1280)xb(1:linsp),xrgstr(1:nxrg)
 1280     FORMAT(a,'-Irregular regression AIC test performed for ',a,
     &             ' regressors.')
         END IF
        END IF
c-----------------------------------------------------------------------
        IF(Khol.eq.1)THEN
         WRITE(Mt1,1290)xb(1:linsp)
 1290    FORMAT(a,'-Prior holiday adjustment factors estimated ',
     &          'for:')
         IF(Keastr.eq.1)WRITE(Mt1,1300)xb(1:(linsp+3))
 1300    FORMAT(a,'- Easter')
        END IF
c-----------------------------------------------------------------------
        IF(Nuspad.gt.0)WRITE(Mt1,1350)xb(1:linsp)
 1350   FORMAT(a,'-Permanent prior adjustment factors will be applied',
     &         ' directly to the final seasonally adjusted series')
        IF(Finhol.AND.(Adjhol.eq.1.or.Axrghl.or.Khol.eq.1))
     &     WRITE(Mt1,1360)xb(1:linsp)
 1360   FORMAT(a,'-Holiday adjustment factors applied directly to the ',
     &         'final seasonally adjusted series')
c-----------------------------------------------------------------------
        IF(Issap.gt.0)THEN
         WRITE(Mt1,1370)xb(1:linsp)
 1370    FORMAT(a,'-Sliding spans analysis performed')
        END IF
c-----------------------------------------------------------------------
        IF(.not.(Prttab(LSPCS0).or.Prttab(LSPCS1).or.Prttab(LSPCS2).or.
     &     Prttab(LSPS1I).or.Prttab(LSPS2I).or.Prttab(LSPS0C)))THEN
         IF(.not.Lcmpaq)WRITE(Mt1,1380)xb(1:linsp),xb(1:(linsp+1))
 1380    FORMAT(a,'-Spectral estimates of original series, table ',
     &          'D11 and table E3 will be searched for ',/,a,
     &          'signficant seasonal and trading day peaks')
        ELSE IF(Thtapr.gt.0)THEN
         WRITE(Mt1,1390)xb(1:linsp),Thtapr
 1390    FORMAT(a,'-Spectral plots generated with rho for Tukey-',
     &          'Hanning taper = ',f6.3)
        ELSE
         IF(.not.Lcmpaq)WRITE(Mt1,1400)xb(1:linsp)
 1400    FORMAT(a,'-Spectral plots generated for selected series')
        END IF
        IF(Pos1(YR).ne.Bgspec(YR).or.Pos1(MO).ne.Bgspec(MO))THEN
         CALL wrtdat(Bgspec,Ny,str,nstr)
         IF(Lfatal)RETURN
         IF(.not.Lcmpaq)WRITE(Mt1,1410)xb(1:linsp),str(1:nstr)
 1410    FORMAT(a,'-Spectral plots generated for series starting in ',a)
        END IF
c-----------------------------------------------------------------------
        IF(Imad.eq.1)WRITE(Mt1,1420)xb(1:linsp)
 1420   FORMAT(a,'-X-11 outlier detection procedure uses moving median',
     &         ' absolute deviations')
        IF(Imad.eq.2)WRITE(Mt1,1430)xb(1:linsp)
 1430   FORMAT(a,'-X-11 outlier detection procedure uses moving median',
     &         ' absolute deviations of the log data')
        IF(Imad.eq.3)WRITE(Mt1,1440)xb(1:linsp)
 1440   FORMAT(a,'-X-11 outlier detection procedure uses tau ',
     &         'adjustment to moving median absolute deviations')
        IF(Imad.eq.4)WRITE(Mt1,1450)xb(1:linsp)
 1450   FORMAT(a,'-X-11 outlier detection procedure uses tau ',
     &         'adjustment to moving median absolute deviations of ',
     &         'the log data')
c-----------------------------------------------------------------------
        IF(Ishrnk.eq.1)WRITE(Mt1,1710)xb(1:linsp),'global'
        IF(Ishrnk.eq.2)WRITE(Mt1,1710)xb(1:linsp),'local'
 1710   FORMAT(a,'-X-11 seasonal factors adjusted using ',a,
     &         'shrinkage factors from Miller and Williams (2003)')
c-----------------------------------------------------------------------
       ELSE IF(Lmodel)THEN
        IF(.not.Prttab(LSPCS0))THEN
         IF(.not.Lcmpaq)WRITE(Mt1,1460)xb(1:linsp)
 1460    FORMAT(a,'-Spectral estimates of original series will be ',
     &          'searched for signficant trading day peaks')
        ELSE IF(Thtapr.gt.0)THEN
         WRITE(Mt1,1470)xb(1:linsp),Thtapr
 1470    FORMAT(a,'-Spectral plot of the original series will be ',
     &          'generated with rho for Tukey-Hanning taper = ',f6.3)
        ELSE
         IF(.not.Lcmpaq)WRITE(Mt1,1480)xb(1:linsp)
 1480    FORMAT(a,'-Spectral plot of the original series generated')
        END IF
       END IF
c-----------------------------------------------------------------------
       IF(Lseats)THEN
        WRITE(Mt1,1800)xb(1:linsp)
 1800   FORMAT(a,'-SEATS model based seasonal adjustment performed.')
       END IF
c-----------------------------------------------------------------------
       IF(Iyrt.gt.0)THEN
        WRITE(Mt1,1130)xb(1:linsp)
 1130   FORMAT(a,'-Modify the D11. series to make the yearly totals ',
     &           'of the seasonally')
        IF(Iftrgt.eq.0)THEN
         WRITE(Mt1,1131)
 1131    FORMAT(14x,'adjusted series agree with the original series.')
        ELSE IF(Iftrgt.eq.1)THEN
         WRITE(Mt1,1132)
 1132    FORMAT(14x,'adjusted series agree with the calendar adjusted',
     &              ' series.')
        ELSE IF(Iftrgt.eq.2)THEN
         WRITE(Mt1,1133)
 1133    FORMAT(14x,'adjusted series agree with the permanent prior ',
     &              'adjusted',/,14x,'series.')
        ELSE IF(Iftrgt.eq.3)THEN
         WRITE(Mt1,1134)
 1134    FORMAT(14x,'adjusted series agree with the calendar and ',
     &              'permanent prior',/,14x,'adjusted series.')
        END IF
        IF(Iyrt.eq.1)THEN
         WRITE(Mt1,1135)xb(1:linsp)
 1135    FORMAT(a,'-Denton method used.')
        ELSE IF(Iyrt.eq.2)THEN
         WRITE(Mt1,1136)xb(1:linsp),Lamda,Rol
 1136    FORMAT(a,'-Regression method used, with lambda = ',f10.7,
     &           ', rho = ',f10.7,'.')
        END IF
c-----------------------------------------------------------------------
        IF(Begyrt.gt.1)THEN
         IF(Ny.eq.12)WRITE(Mt1,1140)xb(1:linsp),Moqu(1:nblank(Moqu)),
     &                              cmonth(Begyrt),cm2(Begyrt)
         IF(Ny.eq.4)WRITE(Mt1,1140)xb(1:linsp),Moqu,cqtr(Begyrt),
     &                             ' Quarter'
 1140    FORMAT(a,'-First ',a,' of fiscal year set to be ',a,a)
        END IF
       END IF
c-----------------------------------------------------------------------
       IF(Lnoprt)WRITE(Mt1,1520)xb(1:linsp)
 1520  FORMAT(a,'-Printout suppressed.  Only user-specified tables and',
     &        ' plots will be printed out.')
      END IF
c-----------------------------------------------------------------------
      IF(Ixreg.eq.2.or.Khol.eq.1)THEN
       WRITE(Mt1,1780)
 1780  FORMAT(//,'   Tables labeled "First pass" are from an initial',
     &           ' seasonal adjustment used to estimate ',/,
     &           '   irregular regression and/or X-11 Easter effects.')
      END IF
c-----------------------------------------------------------------------
      IF (Divpwr.ne.NOTSET) THEN
       WRITE(Mt1,1770)Divpwr
 1770  FORMAT(//,'   All values of original series divided by 10 ** ',
     &           i2,' in this run.')
      END IF
c-----------------------------------------------------------------------
      IF(Savtab(LSRSHD))THEN
       IF(Nsrscr.gt.0)THEN
        WRITE(Nform,1600)'srstit',Srsttl(1:Nsrscr)
       ELSE
        WRITE(Nform,1600)'srstit',Title(1:Ntitle)
       END IF
       WRITE(Nform,1600)'srsnam',Serno
c-----------------------------------------------------------------------
       WRITE(Nform,1610)'freq',Ny
c-----------------------------------------------------------------------
       IF(Ny.eq.12.or.Ny.eq.4)THEN
        WRITE(Nform,1550)Pos1(MO),num(ie5),Moqu(1:nblank(Moqu)),
     &                   Pos1(YR),Lstmo,num(ie6),Moqu(1:nblank(Moqu)),
     &                   Lstyr
       ELSE IF(Ny.eq.1)THEN
        WRITE(Nform,1551)Pos1(YR),Lstyr
*        WRITE(Nform,1553)Bgspec(YR)
       ELSE
        WRITE(Nform,1550)Pos1(MO),num(ie5),'period',Pos1(YR),
     &                    Lstmo,num(ie6),'period',Lstyr
*        ie5=Bgspec(MO)
*        IF(ie5.gt.4)ie5=4
*        WRITE(Nform,1552)Bgspec(MO),num(ie5),'period',Bgspec(YR)
       END IF       
       WRITE(Nform,1610)'nobs',Nspobs       
c-----------------------------------------------------------------------
       IF(Isrflw.eq.1)THEN
        WRITE(Nform,1600)'datatype','flow'
       ELSE IF(Isrflw.eq.2)THEN
        WRITE(Nform,1600)'datatype','stock'
       END IF
c-----------------------------------------------------------------------
       ipos=1
       IF(dpeq(Cnstnt,DNOTST))THEN
        WRITE(Nform,1630)'constant',ZERO
       ELSE
        WRITE(Nform,1630)'constant',Cnstnt
       END IF
c-----------------------------------------------------------------------
       IF(Lmodel)THEN
        CALL prtnfn(Fcntyp,Lam,1)
        WRITE(Nform,1610)'nfcst',Nfcst
        IF(Nfcst.gt.0)WRITE(Nform,1620)'ciprob',Ciprob
        IF(Lognrm)THEN
         WRITE(Nform,1600)'lognormal','yes'
        ELSE
         WRITE(Nform,1600)'lognormal','no'
        END IF
        WRITE(Nform,1630)'mvval',Mvval
        IF(Iqtype.eq.0)THEN
         WRITE(Nform,1600)'iqtype','ljungbox'
        ELSE
         WRITE(Nform,1600)'iqtype','boxpierce'
        END IF
       END IF
       IF(Lx11)THEN
        WRITE(Nform,1660)'samode',malo(ie7)(1:nblank(malo(ie7))),
     &                   runs(ie3)(1:nblank(runs(ie3)))
c-----------------------------------------------------------------------
        WRITE(Nform,1640)'siglim',Sigml,Sigmu
c-----------------------------------------------------------------------
c        IF(Kexopt.eq.1)WRITE(Nform,1580)
c 1580   FORMAT('strike:')
c-----------------------------------------------------------------------
        IF(Iyrt.gt.0)THEN
         WRITE(Nform,1600)'adjtot','yes'
         WRITE(Nform,1610)'adjtotstart',Begyrt
         IF(Iftrgt.eq.0)THEN
          WRITE(Nform,1600)'adjtottarget','original'
         ELSE IF(Iftrgt.eq.1)THEN
          WRITE(Nform,1600)'adjtottarget','original'
         ELSE IF(Iftrgt.eq.2)THEN
          WRITE(Nform,1600)'adjtottarget','pprioradj'
         ELSE
          WRITE(Nform,1600)'adjtottarget','both'
         END IF
 1592    FORMAT('adjtottarget: ',a)
         IF(Iyrt.gt.1)THEN
          WRITE(Nform,1600)'adjtottype','regression'
          WRITE(Nform,1620)'adjtotlambda',Lamda
          WRITE(Nform,1620)'adjtotrho',Rol
          IF(Mid.eq.0)THEN
           WRITE(Nform,1600)'adjtotmode','ratio'
          ELSE
           WRITE(Nform,1600)'adjtotmode','diff'
          END IF
          IF(Lfctfr)THEN
           WRITE(Nform,1600)'adjtotfct','yes'
          ELSE
           WRITE(Nform,1600)'adjtotfct','no'
          END IF
         ELSE
          WRITE(Nform,1600)'adjtottype','denton'
         END IF
        ELSE
         WRITE(Nform,1600)'adjtot','no'
        END IF
c-----------------------------------------------------------------------
        IF(Lterm.eq.NOTSET)THEN
         WRITE(Nform,1600)'seasonalma','None'
        ELSE
         IF(Ny.eq.4)WRITE(Nform,1650)'seasonalma',(sf(j),j=1,4)
         IF(Ny.eq.12)WRITE(Nform,1650)'seasonalma',(sf(j),j=1,12)
        END IF
c-----------------------------------------------------------------------
        IF(Ktcopt.gt.0)THEN
         WRITE(Nform,1610)'trendma',Ktcopt
        ELSE
         WRITE(Nform,1600)'trendma','default'
        END IF
c-----------------------------------------------------------------------
        IF(Kswv.eq.1)WRITE(Nform,1600)'priortd','yes'
c       IF(Kswv.eq.1)WRITE(Nform,1660)Out(1:nblank(Out))
c 1660  FORMAT('priortd: with',a)
c-----------------------------------------------------------------------
        IF(Ixreg.gt.0)THEN
         WRITE(Nform,1600)'x11regress','yes'
         IF(Otlxrg)THEN
          WRITE(Nform,1600)'x11regressextreme','autoao'
          WRITE(Nform,1620)'x11irrcrtval',Critxr
         ELSE IF(Sigxrg.gt.0)THEN
          WRITE(Nform,1600)'x11regressextreme','sigma'
          WRITE(Nform,1620)'x11irrsiglim',Sigxrg
         ELSE
          WRITE(Nform,1600)'x11regressextreme','none'
         END IF
        ELSE
         WRITE(Nform,1600)'x11regress','no'
        END IF
c-----------------------------------------------------------------------
        IF(Khol.eq.1)WRITE(Nform,1600)'x11easter','yes'
c-----------------------------------------------------------------------
        IF(Imad.eq.0)WRITE(Nform,1600)'x11otl','stderr'
        IF(Imad.eq.1)WRITE(Nform,1600)'x11otl','mad'
        IF(Imad.eq.2)WRITE(Nform,1600)'x11otl','madlog'
        IF(Imad.eq.3)WRITE(Nform,1600)'x11otl','taumad'
        IF(Imad.eq.4)WRITE(Nform,1600)'x11otl','taumadlog'
c-----------------------------------------------------------------------
        IF(Ishrnk.eq.0)WRITE(Nform,1600)'shrink','none'
        IF(Ishrnk.eq.1)WRITE(Nform,1600)'shrink','global'
        IF(Ishrnk.eq.2)WRITE(Nform,1600)'shrink','local'
       ELSE
        IF(Lseats)THEN
         WRITE(Nform,1600)'samode','SEATS seasonal adjustment'
        ELSE
         WRITE(Nform,1600)'samode','none'
        END IF
       END IF
       IF (Divpwr.ne.NOTSET) THEN
        WRITE(Nform,1610)'divpower',Divpwr
       END IF
c-----------------------------------------------------------------------
       IF(Ny.eq.12)THEN
        ie5=Bgspec(MO)
        IF(ie5.gt.4)ie5=4
        WRITE(Nform,1600)'spectrum','yes'
        WRITE(Nform,1552)Bgspec(MO),num(ie5),Moqu(1:nblank(Moqu)),
     &                   Bgspec(YR)
        IF(Spctyp.eq.0)THEN
         WRITE(Nform,1600)'spectype','AR-spectrum'
        ELSE
         WRITE(Nform,1600)'spectype','periodogram'
        END IF
        IF(Ldecbl)THEN
         WRITE(Nform,1600)'decibel','yes'
        ELSE
         WRITE(Nform,1600)'decibel','no'
        END IF
        IF(Lrbstsa)THEN
         WRITE(Nform,1600)'specrobustsa','yes'
        ELSE
         WRITE(Nform,1600)'specrobustsa','no'
        END IF
*       IF(Spcl10)THEN
*        WRITE(Nform,1600)'speclog10','yes'
*       ELSE
*        WRITE(Nform,1600)'speclog10','no'
*       END IF
        IF(Spcsrs.eq.0)THEN
         WRITE(Nform,1600)'specseries','original'
        ELSE IF(Spcsrs.eq.1)THEN
         WRITE(Nform,1600)'specseries','outlieradjoriginal'
        ELSE IF(Spcsrs.eq.2)THEN
         WRITE(Nform,1600)'specseries','adjoriginal'
        ELSE IF(Spcsrs.eq.3)THEN
         WRITE(Nform,1600)'specseries','modoriginal'
        END IF
        IF(Lprsfq)THEN
         WRITE(Nform,1600)'showseasonalfreq','yes'
        ELSE
         WRITE(Nform,1600)'showseasonalfreq','no'
        END IF
        IF(Mxarsp.eq.NOTSET)THEN
         WRITE(Nform,1610)'specmaxar',30
        ELSE
         WRITE(Nform,1610)'specmaxar',Mxarsp
        END IF
        IF(Lfqalt)THEN
         WRITE(Nform,1600)'altfreq','yes'
        ELSE
         WRITE(Nform,1600)'altfreq','no'
        END IF
        IF(Svallf)THEN
         WRITE(Nform,1600)'saveallspecfreq','yes'
        ELSE
         WRITE(Nform,1600)'saveallspecfreq','no'
        END IF
        CALL svfreq(Ny,Svallf)
        WRITE(Nform,1630)'peaklocal',Plocal
        WRITE(Nform,1610)'peakwd',Peakwd
        WRITE(Nform,1620)'siglevel',Spclim
       ELSE
        WRITE(Nform,1600)'spectrum','no'
       END IF
c-----------------------------------------------------------------------
       IF(Iag.ge.0)THEN
        IF(Iag.eq.0)THEN
         WRITE(Nform,1600)'comptype','add'
        ELSE IF(Iag.eq.1) THEN
         WRITE(Nform,1600)'comptype','sub'
        ELSE IF(Iag.eq.2) THEN
         WRITE(Nform,1600)'comptype','mult'
        ELSE IF(Iag.eq.3) THEN
         WRITE(Nform,1600)'comptype','div'
        END IF
        WRITE(Nform,1620)'compwt',W
       END IF 
c-----------------------------------------------------------------------
      END IF
c-----------------------------------------------------------------------
c   Write out user defined metadata.
c-----------------------------------------------------------------------
      IF(Hvmtdt)THEN
       DO i=1,Nval
        CALL getstr(Keystr,Keyptr,Nkey,i,thisky,ikey)
        IF(Lfatal)RETURN
        CALL getstr(Valstr,Valptr,Nval,i,thisvl,ival)
        IF(Lfatal)RETURN
        WRITE(Nform,1600)'metadata.'//thisky(1:ikey),thisvl(1:ival)
       END DO
      END IF
c-----------------------------------------------------------------------
 1550 FORMAT('span: ',i2,a2,1x,a,',',i4,' to ',i2,a2,1x,a,',',i4)
 1551 FORMAT('span: ',i4,' to ',i4)
 1552 FORMAT('startspec: ',i2,a2,1x,a,',',i4)
 1553 FORMAT('startspec: ',i4)
 1600 FORMAT(a,': ',a)
 1610 FORMAT(a,': ',i5)
 1620 FORMAT(a,': ',f12.6)
 1630 FORMAT(a,': ',e20.10)
 1640 FORMAT(a,':',2(1x,f12.6))
 1650 FORMAT(a,': ',11(a4,2x),a4)
 1660 FORMAT(a,': ',a,1x,a)
c-----------------------------------------------------------------------
      RETURN
c-----------------------------------------------------------------------
      END

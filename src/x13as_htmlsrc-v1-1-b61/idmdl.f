C     Last change:  BCM  16 Feb 1999   11:16 am
      SUBROUTINE idmdl(Dflist,Niddf,Nidsdf,Mxidlg,Lgraf)
      IMPLICIT NONE
c     ------------------------------------------------------------------
      LOGICAL F,T
      DOUBLE PRECISION ONE,ZERO
      PARAMETER(ONE=1D0,ZERO=0D0,F=.false.,T=.true.)
c     ------------------------------------------------------------------
      INCLUDE 'stdio.i'
      INCLUDE 'acfptr.prm'
      INCLUDE 'srslen.prm'
      INCLUDE 'model.prm'
      INCLUDE 'model.cmn'
      INCLUDE 'mdldat.cmn'
      INCLUDE 'tbllog.i'
      INCLUDE 'mdltbl.i'
      INCLUDE 'units.cmn'
      INCLUDE 'error.cmn'
      INCLUDE 'tbllog.prm'
      INCLUDE 'tbllog.cmn'
c     ------------------------------------------------------------------
      INTEGER PXY
      PARAMETER(PXY=PLEN*(PB+1))
c     ------------------------------------------------------------------
      INTEGER cnstcl,Dflist,endopr,idf,igrp,isdf,itmp,Mxidlg,mxndf,
     &        mxnsdf,na,ndf,nefobs,Niddf,Nidsdf,nsdf
      DOUBLE PRECISION a,txy
      LOGICAL Lgraf,locok
      DIMENSION a(PLEN+2*PORDER),Dflist(PDFLG,2),txy(PXY)
c-----------------------------------------------------------------------
      INTEGER strinx
      EXTERNAL strinx
c-----------------------------------------------------------------------
      INTEGER Fhacf,Fhpcf,Fhacfg,Fhpcfg
      COMMON /cfhacf/ Fhacf,Fhpcf,Fhacfg,Fhpcfg
c-----------------------------------------------------------------------
c     Do a regression a regression on the series and variables
c differenced by the maximum order of differencing found in the
c difference lists.
c-----------------------------------------------------------------------
      IF(Nb.gt.0)THEN
       CALL maxidx(Dflist,Niddf,itmp,mxndf)
       CALL maxidx(Dflist(1,2),Nidsdf,itmp,mxnsdf)
       CALL copy(Xy,Nspobs*Ncxy,1,txy)
       CALL difflt(Nspobs,Ncxy,mxndf,mxnsdf,Sp,txy,nefobs)
c-----------------------------------------------------------------------
c  check to see if there are enough observations to perform ols
c  regression/acfs (BCM 2-2000)
c-----------------------------------------------------------------------
       IF(nefobs.le.0)THEN
        CALL eWritln(
     &     'Not enough data to perform maximum order of differencing',
     &     STDERR,Mt2,T,F)
        CALL writln('       specified in the diff and sdiff arguments'//
     &              ' of the identify spec.',STDERR,Mt2,F,T)
        CALL abend()
        RETURN
       END IF
c-----------------------------------------------------------------------
       igrp=strinx(F,Grpttl,Grpptr,1,Ngrptl,'Constant')
       IF(igrp.gt.0)THEN
        IF(.not.Lquiet)WRITE(STDERR,1010)
 1010   FORMAT(/,' WARNING: For calculating the ACF''s and PACF''s ',
     &           'requested from the identify',/,
     &           '          spec, a sample mean adjustment has been ',
     &           'used in place of the',/,
     &           '          effect of the constant regressor ',
     &           'specified in the regression spec.',/)
        CALL wWritln('For calculating the <abbr '//
     &               'title="autocorrelation functions">ACF''s</abbr>'//
     &               ' and <abbr title="partial autocorrelation '//
     &               'functions">PACF''s</abbr>',
     &               Mt1,Mt2,T,F)
        CALL writln('requested from the identify spec, a sample mean '//
     &              'adjustment has been used in place of the',
     &               Mt1,Mt2,F,F)
        CALL writln(' effect of the constant regressor specified '//
     &              ' in the regression spec.',Mt1,Mt2,F,T)
        CALL setdp(ONE,nefobs,a)
        cnstcl=Grp(igrp-1)
        CALL copycl(a,nefobs,1,1,Ncxy,cnstcl,txy)
       ELSE
        cnstcl=0
       END IF
c     ------------------------------------------------------------------
       CALL olsreg(txy,nefobs,Ncxy,Ncxy,B,Chlxpx,PXPX,Sngcol)
       IF(Lfatal)RETURN
c     ------------------------------------------------------------------
       IF(Sngcol.gt.0)THEN
        Convrg=F
        Armaer=PISNER
        CALL prterr(nefobs,F)
        IF(Lfatal)RETURN
       END IF
c     ------------------------------------------------------------------
       IF(cnstcl.gt.0)B(cnstcl)=ZERO
      END IF
      itmp=0
      CALL resid(Xy,Nspobs,Ncxy,Ncxy,1,Nb,-ONE,B,a)
      IF(Lfatal)RETURN
      CALL yprmy(a,Nspobs,Var)
      Var=Var/Nspobs
c-----------------------------------------------------------------------
c     Print the regression estimates.  Trick prtmdl in thinking that
c there is no ARMA model.
c-----------------------------------------------------------------------
      IF(Nb.gt.0)THEN
       endopr=Mdl(MA)
       Mdl(MA)=1
       IF(Prttab(LIDRGC))THEN
        CALL prtmdl(F,T,F,F,F,F,F,F,F,itmp,F,F,F)
        IF(Lfatal)RETURN
       END IF
       Mdl(MA)=endopr
      END IF
c-----------------------------------------------------------------------
c     Set up file handles needed to save acfs and pacfs generated by the
c     identify spec
c-----------------------------------------------------------------------
      locok=T
      IF(Savtab(LSPIDN+LACF))
     &   CALL opnfil(.true.,F,LSPIDN+LACF,Fhacf,locok)
      IF(Savtab(LSPIDN+LPCF).and.locok)
     &   CALL opnfil(.true.,F,LSPIDN+LPCF,Fhpcf,locok)
      IF(Lgraf.and.locok)THEN
       CALL opnfil(.true.,Lgraf,LSPIDN+LACF,Fhacfg,locok)
       IF(locok)
     &    CALL opnfil(.true.,Lgraf,LSPIDN+LPCF,Fhpcfg,locok)
      END IF
      IF(.not.locok)THEN
       CALL abend
       RETURN
      END IF
c-----------------------------------------------------------------------
c     Print the acf and pacf for all the orders of differencing
c requested.
c-----------------------------------------------------------------------
      DO idf=1,Niddf
       ndf=Dflist(idf,1)
c     ------------------------------------------------------------------
       DO isdf=1,Nidsdf
        nsdf=Dflist(isdf,2)
c     ------------------------------------------------------------------
*        IF(ndf.eq.0)THEN
*         IF(nsdf.eq.0)THEN
*          WRITE(Mt1,1020)
* 1020     FORMAT(/,' Differencing:  none')
*c     ------------------------------------------------------------------
*         ELSE
*          WRITE(Mt1,1030)nsdf
* 1030     FORMAT(/,' Differencing:  Seasonal Order=',i1)
*         END IF
*c     ------------------------------------------------------------------
*        ELSE IF(nsdf.eq.0)THEN
*         WRITE(Mt1,1040)ndf
* 1040    FORMAT(/,' Differencing:  Nonseasonal Order=',i1)
*c     ------------------------------------------------------------------
*        ELSE
*         WRITE(Mt1,1050)ndf,nsdf
* 1050    FORMAT(/,' Differencing:  Nonseasonal Order=',i1,
*     &          ', Seasonal Order=',i1)
*        END IF
c     ------------------------------------------------------------------
        CALL copy(a,Nspobs,1,txy)
        CALL difflt(Nspobs,1,ndf,nsdf,Sp,txy,na)
        CALL prtacf(LSPIDN,na,txy,na,Mxidlg,Lgraf,F,ndf,nsdf)
        IF(Lfatal)RETURN
       END DO
      END DO
c     ------------------------------------------------------------------
      IF(Savtab(LSPIDN+LACF))CALL fclose(Fhacf)
      IF(Savtab(LSPIDN+LPCF))CALL fclose(Fhpcf)
      IF(Lgraf)THEN
       CALL fclose(Fhacfg)
       CALL fclose(Fhpcfg)
      END IF
c     ------------------------------------------------------------------
      RETURN
      END

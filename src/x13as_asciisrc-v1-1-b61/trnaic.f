C     Last change:  SRD  19 Nov 99    6:37 am
      SUBROUTINE trnaic(Lx11,Lmodel,Lprt,Lprtfm)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     Estimate regARIMA model for the untransformed and log transformed
c     series.  The routine will choose the model with the lowest value 
c     of AICC and print out the resulting model.
c-----------------------------------------------------------------------
      DOUBLE PRECISION ONE,ZERO
      LOGICAL F,T
      INTEGER DIV
      PARAMETER(ONE=1D0,ZERO=0D0,DIV=4,F=.false.,T=.true.)
c-----------------------------------------------------------------------
      INCLUDE 'stdio.i'
      INCLUDE 'srslen.prm'
      INCLUDE 'model.prm'
      INCLUDE 'model.cmn'
      INCLUDE 'mdldat.cmn'
      INCLUDE 'picktd.cmn'
      INCLUDE 'arima.cmn'
      INCLUDE 'lkhd.cmn'
      INCLUDE 'adj.cmn'
      INCLUDE 'prior.prm'
      INCLUDE 'prior.cmn'
      INCLUDE 'priusr.cmn'
      INCLUDE 'extend.cmn'
      INCLUDE 'x11adj.cmn'
      INCLUDE 'x11fac.cmn'
      INCLUDE 'x11opt.cmn'
      INCLUDE 'notset.prm'
      INCLUDE 'hiddn.cmn'
      INCLUDE 'inpt.cmn'
      INCLUDE 'units.cmn'
      INCLUDE 'error.cmn'
      INCLUDE 'title.cmn'
      INCLUDE 'x11ptr.cmn'
      INCLUDE 'x11srs.cmn'
      INCLUDE 'svllog.prm'
      INCLUDE 'svllog.cmn'
      INCLUDE 'mdlsvl.i'
      INCLUDE 'mq3.cmn'
c-----------------------------------------------------------------------
      INTEGER PA
      PARAMETER(PA=PLEN+2*PORDER)
c-----------------------------------------------------------------------
      LOGICAL Lprt,lhide,argok,Lprtfm,Lx11,rok,lam2,lad2,lax2,Lmodel,
     &        begrgm,inptok
      DOUBLE PRECISION a,a2,aicno,aiclog,lomeff,Maxsrs,Minsrs,Temp,
     &                 trnsrs
      INTEGER Frstry,kf2,na,Nefobs,nbeg,nend,endlag,ilag,fhnote
      DIMENSION Temp(PLEN),a(PA),lomeff(PLEN),trnsrs(PLEN),a2(PLEN),
     &          begrgm(PLEN)
c-----------------------------------------------------------------------
      INTEGER nblank,strinx
      EXTERNAL nblank,strinx
c-----------------------------------------------------------------------
      COMMON /maxmin/ Maxsrs,Minsrs
      COMMON /work  / Temp
c-----------------------------------------------------------------------
c     Set up Fcntyp and Lam to do no transformation
c-----------------------------------------------------------------------
      Fcntyp=4
      Lam=ONE
      Ixreg=-Ixreg
      CALL setdp(ONE,PLEN,lomeff)
      CALL copy(Adj,PLEN,1,a2)
      fhnote=STDERR
      IF(Lquiet)fhnote=0
c-----------------------------------------------------------------------
c     Generate leap year effect for possible later use.
c-----------------------------------------------------------------------
      IF(Picktd.and.Lmodel)THEN
       IF(Lrgmtd.AND.(MOD(Tdzero,2).ne.0))THEN
        CALL gtrgpt(Begadj,Tddate,Tdzero,begrgm,Nadj)
       ELSE
        CALL setlg(T,PLEN,begrgm)
       END IF
       IF(Kfulsm.eq.2)THEN
        CALL td7var(Begadj,Sp,Nadj,1,1,T,F,T,lomeff,begrgm)
       ELSE
        CALL td7var(Begadj,Sp,Nadj,1,1,F,F,T,lomeff,begrgm)
       END IF
      END IF
c-----------------------------------------------------------------------
c     Change Begspn and Endspn to match the model span, if necessary.
c-----------------------------------------------------------------------
      CALL dfdate(Begmdl,Begspn,Sp,nbeg)
      CALL dfdate(Endspn,Endmdl,Sp,nend)
      IF(nbeg.gt.0)CALL cpyint(Begmdl,2,1,Begspn)
      IF(nend.gt.0)CALL cpyint(Endmdl,2,1,Endspn)
c-----------------------------------------------------------------------
c     Process the series
c-----------------------------------------------------------------------
      CALL dfdate(Endspn,Begspn,Sp,Nspobs)
      Nspobs=Nspobs+1
      IF(nbeg.gt.0.or.nend.gt.0.or.Issap.eq.2)THEN
       CALL dfdate(Begspn,Begsrs,Sp,Frstsy)
       Frstsy=Frstsy+1
       Nomnfy=Nobs-Frstsy+1
       Nobspf=min(Nspobs+max(Nfcst-Fctdrp,0),Nomnfy)
      END IF
      CALL copy(Orig(Pos1ob+nbeg),Nobspf,-1,trnsrs)
c-----------------------------------------------------------------------
c     If no model specified or automatic model specified, use default
c     model (0 1 1)(0 1 1) or (0 1 1).
c-----------------------------------------------------------------------
      IF((.not.Lmodel).OR.Lautom.or.Lautox)THEN
       inptok=T
       CALL mdlint()
       IF(Lmodel.and.Lseff)THEN
        CALL mdlset(0,1,1,0,0,0,inptok)
       ELSE
        CALL mdlset(0,1,1,0,1,1,inptok)
       END IF
       IF((.not.inptok).or.Lfatal)THEN
        WRITE(STDERR,1070)
        WRITE(Mt2,1070)
        RETURN
       END IF
      END IF
c-----------------------------------------------------------------------
c     Set up the regression matrix
c-----------------------------------------------------------------------
      CALL regvar(Trnsrs,Nobspf,Fctdrp,Nfcst,0,Userx,Bgusrx,Nrusrx,
     &            Priadj,Reglom,Nrxy,Begxy,Frstry,T,Elong)
      IF(Lfatal)RETURN
c-----------------------------------------------------------------------
c     Estimate the regression and ARMA parameters
c-----------------------------------------------------------------------
      lax2=Lautox
      Lautox=F
      lam2=Lautom
      Lautom=F
      lad2=Lautod
      Lautod=F
      argok=Lautom
      CALL rgarma(T,Mxiter,Mxnlit,F,a,na,Nefobs,argok)
      IF(Lfatal)RETURN
c-----------------------------------------------------------------------
c     If an estimation error is found, discontinue the routine.
c-----------------------------------------------------------------------
      IF(Armaer.eq.PMXIER.or.Armaer.eq.PSNGER.or.Armaer.eq.PISNER.or.
     &   Armaer.eq.PNIFER.or.Armaer.eq.PNIMER.or.Armaer.eq.PCNTER.or.
     &   Armaer.eq.POBFN0.or.Armaer.lt.0)THEN
       WRITE(STDERR,1071)'untransformed'
       WRITE(Mt1,1071)'untransformed'
       WRITE(Mt2,1071)'untransformed'
       CALL prterr(nefobs,F)
       CALL abend()
       RETURN
c-----------------------------------------------------------------------
c     If only a warning message would be printed out, reset the error
c     indicator variable to zero.
c-----------------------------------------------------------------------
      ELSE IF(Armaer.ne.0)THEN
       Armaer=0
      END IF
c-----------------------------------------------------------------------
c     Compute the likelihood statistics and AICC for the first model
c-----------------------------------------------------------------------
      IF(.not.Lprt)THEN
       lhide=Lhiddn
       Lhiddn=T
      END IF
      IF(Lprt)WRITE(Mt1,1010)
      CALL prlkhd(Y(Frstsy),Adj(Adj1st),Adjmod,Fcntyp,Lam,F,Lprt,F)
      IF(Lfatal)RETURN
      aicno=Aicc
      IF(Svltab(LSLTRN))WRITE(Ng,1011)Aicc
      IF(Lsumm.gt.0)WRITE(Nform,1012)'nolog',Aicc
c-----------------------------------------------------------------------
c     Restore initial values to model variables
c-----------------------------------------------------------------------
c      CALL restor(T,F,F)
      IF(Nopr.gt.0)THEN
       endlag=Opr(Nopr)-1
       DO ilag=1,endlag
        IF(.not.Arimaf(ilag))Arimap(ilag)=Ap1(ilag)
       END DO
      END IF
c-----------------------------------------------------------------------
c     Perform log transformation on series, and reset indicator 
c     variables accordingly
c-----------------------------------------------------------------------
      Fcntyp=1
      Lam=ZERO
      Adjmod=1 
      CALL setdp(ONE,PLEN,Adj)
c-----------------------------------------------------------------------
c     If model has TD, perform leap year prior adjustment
c-----------------------------------------------------------------------
      IF(Lmodel.and.Picktd)THEN
c-----------------------------------------------------------------------
c     Remove leap year regressor
c-----------------------------------------------------------------------
       CALL rmlnvr(Priadj,Kfulsm,Nobs)
c-----------------------------------------------------------------------
c     Leap Year adjust original series, if necessary.
c-----------------------------------------------------------------------
       kf2=Kfmt
c-----------------------------------------------------------------------
c     Perform leap year adjustment
c-----------------------------------------------------------------------
       CALL addate(Begspn,Sp,-Nbcst,Begadj)
       Nadj=Nspobs+Nbcst+max(Sp,Nfcst-Fctdrp)
       CALL eltfcn(DIV,Y(Frstsy),lomeff,Nspobs,PLEN,trnsrs)
       CALL copy(lomeff,PLEN,1,Adj)
       CALL dfdate(Begspn,Begadj,Sp,Adj1st)
       Adj1st=Adj1st+1
      END IF
      IF(Lmvaft.or.Ln0aft)THEN
       CALL trnfcn(trnsrs,Nspobs,Fcntyp,Lam,trnsrs)
      ELSE
       CALL trnfcn(trnsrs,Nobspf,Fcntyp,Lam,trnsrs)
      END IF
      IF(Lfatal)RETURN
c-----------------------------------------------------------------------
c     estimate model for log transformation
c-----------------------------------------------------------------------
      CALL regvar(Trnsrs,Nobspf,Fctdrp,Nfcst,0,Userx,Bgusrx,Nrusrx,
     &            Priadj,Reglom,Nrxy,Begxy,Frstry,T,F)
      IF(.not.Lfatal)CALL rgarma(T,Mxiter,Mxnlit,F,A,Na,Nefobs,argok)
      IF(Lfatal)RETURN
c-----------------------------------------------------------------------
c     If an estimation error is found, discontinue the routine.
c-----------------------------------------------------------------------
      IF(Armaer.eq.PMXIER.or.Armaer.eq.PSNGER.or.Armaer.eq.PISNER.or.
     &   Armaer.eq.PNIFER.or.Armaer.eq.PNIMER.or.Armaer.eq.PCNTER.or.
     &   Armaer.eq.POBFN0.or.Armaer.lt.0)THEN
       WRITE(STDERR,1071)'log transformed'
       WRITE(Mt1,1071)'log transformed'
       WRITE(Mt2,1071)'log transformed'
       CALL prterr(nefobs,F)
       CALL abend()
       RETURN
      END IF
c-----------------------------------------------------------------------
c     Compute and print out the likelihood statistics and AICC for the
c     model with TD
c-----------------------------------------------------------------------
      IF(Lprt)WRITE(Mt1,1020)
      CALL prlkhd(Y(Frstsy),Adj(Adj1st),Adjmod,Fcntyp,Lam,F,Lprt,Lprtfm)
      IF(Lfatal)RETURN
      aiclog=Aicc
      IF(Svltab(LSLTRN))WRITE(Ng,1021)Aicc
      IF(Lsumm.gt.0)WRITE(Nform,1012)'log',Aicc
      IF(.not.Lprt)Lhiddn=lhide
c-----------------------------------------------------------------------
c     Reset beginning and ending dates for span, if necessary.
c-----------------------------------------------------------------------
      IF(nbeg.gt.0.or.nend.gt.0)
     &   CALL setspn(Sp,nend,nbeg,Begspn,Endspn,Begmdl,Endmdl,Nspobs,
     &               Frstsy,Nobspf,Begsrs,Nobs,Nfcst,Fctdrp,Nomnfy,
     &               Begadj,Adj1st)
c-----------------------------------------------------------------------
      IF((aiclog+Traicd).lt.aicno)THEN
       IF(Lprt)THEN
        WRITE(Mt1,1030)Traicd
        IF(Lx11)WRITE(Mt1,1031)
        WRITE(Mt1,1032)
       END IF
       IF(Svltab(LSLTRN))THEN
        IF(Lmodel)THEN
         WRITE(Ng,1033)'Log Transformation'
        ELSE
         WRITE(Ng,1033)'Multiplicative Seasonal Adjustment'
        END IF
       END IF
c-----------------------------------------------------------------------
c     Set up variables used in X-11 for prior adjustment (if this model
c     has td).
c-----------------------------------------------------------------------
       IF(nbeg.gt.0.or.nend.gt.0)THEN
        Nadj=Nspobs+Nbcst+max(Sp,Nfcst-Fctdrp)
        IF(Picktd)CALL copy(lomeff,PLEN,1,Adj)
       END IF
c-----------------------------------------------------------------------
       CALL copy(Adj,Nadj,-1,Sprior(Setpri))
       IF(Lmodel.and.Picktd)THEN
        CALL ssprep(T,F,F)
        IF(Kfmt.eq.0)THEN
         Kfmt=1
         Prmser='LPY'
         IF(Kfulsm.lt.2)THEN
          Prmser='LPY'
         ELSE
          IF(Sp.eq.12)THEN
           Prmser='LOM'
          ELSE IF(Sp.eq.4)THEN
           Prmser='LOQ'
          END IF
         END IF
        END IF
       END IF
c-----------------------------------------------------------------------
c     redo table formats
c-----------------------------------------------------------------------
       CALL tfmts(Sp,Kdec,Maxsrs,Minsrs,0,Lwdprt,rok)
c-----------------------------------------------------------------------
c     If X-11 seasonal adjustment to be performed, set Muladd
c-----------------------------------------------------------------------
       Muladd=0
       Tmpma=0
       Pcdif='percent change '
       Rad='ratios     '
c-----------------------------------------------------------------------
c     Re-initialize several variable used in program to 1.0 instead of
c     0.0
c-----------------------------------------------------------------------
       CALL setdp(ONE,PLEN,Sts)
       CALL setdp(ONE,PLEN,Stsi)
       CALL setdp(ONE,PLEN,Sti)
       CALL setdp(ONE,PLEN,Stptd)
       CALL setdp(ONE,PLEN,Temp)
       CALL setdp(ONE,PLEN,Factd)
       CALL setdp(ONE,PLEN,Facao)
       CALL setdp(ONE,PLEN,Facls)
       CALL setdp(ONE,PLEN,Factc)
       CALL setdp(ONE,PLEN,Facso)
       CALL setdp(ONE,PLEN,Facsea)
       CALL setdp(ONE,PLEN,Facusr)
       CALL setdp(ONE,PLEN,Fachol)
       CALL setdp(ONE,PLEN,Facxhl)
       CALL setdp(ONE,PLEN,X11hol)
       CALL setdp(ONE,PLEN,Faccal)
c-----------------------------------------------------------------------
       IF(Ln0aft.or.Lmvaft)THEN
        CALL writln('NOTE: At least one value that is either less than o
     &r equal to zero or',fhnote,Mt2,T)
        CALL writln('      equal to the missing value code was found aft
     &er the span of data',fhnote,Mt2,F)
        CALL writln('      to be analyzed, but within the time frame of
     &the forecasts',fhnote,Mt2,F)
        CALL writln('      generated by the regARIMA model.',fhnote,Mt2,
     &              F)
        CALL writln('      In this situation, the forecast output will n
     &ot include a',fhnote,Mt2,T)
        CALL writln('      comparison of the transformed forecasts with 
     &the corresponding',fhnote,Mt2,F)
        CALL writln('      values of the transformed original series.',
     &               fhnote,Mt2,F)
       END IF
c-----------------------------------------------------------------------
      ELSE
       IF(Lprt)THEN
        WRITE(Mt1,1040)Traicd
        IF(Lx11)WRITE(Mt1,1041)
        WRITE(Mt1,1032)
       END IF
       IF(Svltab(LSLTRN))THEN
        IF(Lmodel)THEN
         WRITE(Ng,1033)'No Transformation'
        ELSE
         WRITE(Ng,1033)'Additive Seasonal Adjustment'
        END IF
       END IF
       Fcntyp=4
       Lam=ONE
       Adjmod=2
       Priadj=1
       IF(Ln0aft)Ln0aft=F
       IF(Lmvaft)THEN
        CALL writln('NOTE: At least one value that is either less than o
     &r equal to zero or',fhnote,Mt2,T)
        CALL writln('      equal to the missing value code was found aft
     &er the span of data',fhnote,Mt2,F)
        CALL writln('      to be analyzed, but within the time frame of
     &the forecasts',fhnote,Mt2,F)
        CALL writln('      generated by the regARIMA model.',fhnote,Mt2,
     &              F)
        CALL writln('      In this situation, the forecast output will n
     &ot include a',fhnote,Mt2,T)
        CALL writln('      comparison of the forecasts with the correspo
     &nding values of the',fhnote,Mt2,F)
        CALL writln('      original series.',fhnote,Mt2,F)
       END IF
       CALL setdp(ZERO,PLEN,Adj)
       CALL copy(Orig(Pos1ob+nbeg),Nobspf,-1,trnsrs)
c-----------------------------------------------------------------------
c     If the first model is better, restore original series and 
c     reset LOM variables, if no td in first model
c-----------------------------------------------------------------------
       IF(Lmodel)THEN
        IF(Picktd)THEN
         CALL adrgef(DNOTST,'Leap Year','Leap Year',PRGTLY,F,F)
         IF(Lfatal)RETURN
         CALL copy(a2,PLEN,1,Adj)
         Kfmt=kf2
        END IF
       END IF
      END IF
c-----------------------------------------------------------------------
      Lautox=lax2
      Lautom=lam2
      Lautod=lad2
      IF(Lsumm.gt.0)THEN
       IF(Lmodel)CALL prtnfn(Fcntyp,Lam,2)
       IF(Lx11)THEN
        IF(Muladd.eq.0)THEN
         WRITE(Nform,1050)'multiplicative'
        ELSE
         WRITE(Nform,1050)'additive'
        END IF
       END IF
      END IF
c-----------------------------------------------------------------------
      Ixreg=-Ixreg
c-----------------------------------------------------------------------
 1010 FORMAT(//,
     &  ' Likelihood statistics for model fit to untransformed series.')
 1011 FORMAT(' AICC(no log) : ',f15.4)
 1012 FORMAT('aictest.trans.aicc.',a,': ',e29.15)
 1020 FORMAT(//,
     &' Likelihood statistics for model fit to log transformed series.')
 1021 FORMAT(' AICC(log)    : ',f15.4)
 1030 FORMAT(//,'   *****   AICC (with aicdiff=',F5.2,
     &          ') prefers log transformation   *****')
 1031 FORMAT('   *****   Multiplicative seasonal adjustment will be perf
     &ormed.   ****')
 1032 FORMAT(//)
 1033 FORMAT(/,'  Automatic transformation test : ',a)
 1040 FORMAT(//,'   *****   AICC (with aicdiff=',F5.2,
     &          ') prefers no transformation   *****')
 1041 FORMAT('   *****   Additive seasonal adjustment will be performed.
     &   ****')
 1050 FORMAT('finmode: ',a)
 1070 FORMAT(/,' ERROR: Unable to set up ARIMA model for automatic ',
     &         'transformation selection',/,
     &         '        procedure for the reason(s) given above.')
 1071 FORMAT(/,'  Estimation error found during automatic ',
     &         'transformation selection',
     &       /,'  procedure while fitting regARIMA model to the ',a,
     &         ' series.')
c-----------------------------------------------------------------------
      RETURN
      END

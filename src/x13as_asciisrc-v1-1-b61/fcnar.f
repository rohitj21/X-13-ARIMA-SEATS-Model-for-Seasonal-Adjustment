C     Last change:  BCM  14 May 1998    9:17 am
      SUBROUTINE fcnar(Na,Testpm,Estprm,A,Lauto,Gudrun,Err,Lckinv)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     fcnar.f, Release 1, Subroutine Version 1.10, Modified 14 Feb 1995.
c-----------------------------------------------------------------------
c	This routine works in the nonlinear routine on the
c regression residuals.  Calculates the G'G matrix since
c the parameters have changed.
c-----------------------------------------------------------------------
c     Subroutine to calculate exact MA ARIMA filter residuals.
c Setmdl differences the X:y matrix and changes the model to remove
c the differencing.  The the remaining ARMA model is estimated.
c Model information is in ARIMA.cmn common so the variables are saved
c between calls of the routines fcnar, and arflt. Setmdl also
c constructs a vector of parameters to be estimated in the nonlinear
c routine.  Fcnar calculates ARIMA filter residuals given new estimated
c parameters, estprm, from the nonlinear routine, regression residuals,
c tsrs, from rgcpnt, and the model information that was constructed
c in setmdl.  ARflt filters an extended [X:y] matrix from rgcpnt
c using parameter estimates saved during the last fcnar call.
c-----------------------------------------------------------------------
c Name  Type Description
c-----------------------------------------------------------------------
c a       d  Output na long vector of the deviances.
c err     i  Output error warning to have the nonlinear routine
c             terminate the program (not used).
c estprm  d  Input nestpm long vector of estimated parameters from the
c             nonlinear routine.  Nestpm is found in model.cmn
c estptr  i  Local pointer in either estprm or arimap for the first operator
c             to be expanded.
c i       i  Local do loop parameter
c iflt    i  Local index for the current filter type, DIFF, AR, or MA.
c ilag    i  Local index for the current lag, pointer to the current
c             element in lag,arimap, and arimaf.
c iopr    i  Local index for the current operator, it is the pointer to the
c             current row in the operator specfication matrix, opr.
c lagptr  i  Local pointer to the current coefficient and lag in arimap
c             and arimal
c na      i  Input number of a's or the number of residuals expected by
c             the nonlinear routine, nefobs+order of the MA operator
c nlag    i  Local number of lags in the current operator of a filter
c nopr    i  Local for the number of operators in a DIFF, AR, or MA filter.
c one     d  Local PARAMETER for a double precision 1
c Testpm  i  Input dummy number which should be the same as nestpm, the
c             of parmeters in estprm
c zero    d  Local PARAMETER for double precision 0
c-----------------------------------------------------------------------
      INCLUDE 'stdio.i'
      INCLUDE 'notset.prm'
      INCLUDE 'srslen.prm'
      INCLUDE 'model.prm'
      INCLUDE 'series.cmn'
      INCLUDE 'units.cmn'
      INCLUDE 'model.cmn'
      INCLUDE 'mdldat.cmn'
      INCLUDE 'error.cmn'
c     ------------------------------------------------------------------
      INTEGER PA
      DOUBLE PRECISION TWO
      PARAMETER(TWO=2D0,PA=PLEN+2*PORDER)
c     ------------------------------------------------------------------
      CHARACTER cdot*(1),tmpttl*(POPRCR)
      LOGICAL Gudrun,Lckinv,Lauto
      INTEGER Err,info,Na,ntmpcr,Testpm,fh2
      DOUBLE PRECISION A,Estprm,fac
      DIMENSION A(PA),Estprm(Nestpm)
c-----------------------------------------------------------------------
      cdot='.'
      IF(Lauto)cdot=' '
      fh2=0
      IF(.not.Lauto)THEN
       fh2=Mt1
       IF(.not.Gudrun)fh2=0
      END IF
c-----------------------------------------------------------------------
c     Insert the estimated parameters in the model into the ARIMA
c filtering data structures.
c-----------------------------------------------------------------------
      CALL upespm(Estprm)
c-----------------------------------------------------------------------
c     Call the ARIMA filter
c-----------------------------------------------------------------------
      CALL copy(Tsrs,Nspobs,1,A)
      CALL armafl(Nspobs,1,.true.,Lckinv,A,Na,PA,info)
c-----------------------------------------------------------------------
c     Print the warning messages here because
c-----------------------------------------------------------------------
      IF(info.ne.0)THEN
       IF(Lprier)THEN
        IF(info.eq.PINVER)THEN
         CALL getstr(Oprttl,Oprptr,Noprtl,Prbfac,tmpttl,ntmpcr)
         IF(Lfatal)RETURN
         IF(fh2.gt.0)WRITE(fh2,1010)tmpttl(1:ntmpcr),cdot
         CALL errhdr
         WRITE(Mt2,1010)tmpttl(1:ntmpcr),cdot
c     ------------------------------------------------------------------
        ELSE IF(info.eq.PGPGER)THEN
         IF(fh2.gt.0)WRITE(fh2,1020)PRGNAM,cdot
         CALL errhdr
         WRITE(Mt2,1020)PRGNAM,cdot
c     ------------------------------------------------------------------
        ELSE IF(info.eq.PACFER)THEN
         IF(fh2.gt.0)WRITE(fh2,1030)cdot
         CALL errhdr
         WRITE(Mt2,1030)cdot
c     ------------------------------------------------------------------
        ELSE IF(info.eq.PVWPER)THEN
         IF(fh2.gt.0)WRITE(fh2,1040)cdot
         CALL errhdr
         WRITE(Mt2,1040)cdot
        END IF
c     ------------------------------------------------------------------
        IF(Lckinv)THEN
         CALL errhdr
         IF(Lauto)THEN
          WRITE(Mt2,1049)Mdldsn(1:Nmddcr)
         ELSE
          IF(fh2.gt.0)THEN
           WRITE(fh2,1050)
           CALL prtitr(A,Na,Estprm,Nestpm,' ',NOTSET,NOTSET)
          END IF
          WRITE(Mt2,1050)
         END IF
c     ------------------------------------------------------------------
        ELSE
         CALL errhdr
         IF(fh2.gt.0)WRITE(fh2,1060)
         WRITE(Mt2,1060)
        END IF
       END IF
c-----------------------------------------------------------------------
c     Make the residuals so big a bad jump will be brought back inbounds
c-----------------------------------------------------------------------
c       CALL dcopy(Na,Lrgrsd,0,A,1)
       CALL setdp(Lrgrsd,Na,A)
       info=0
       Err=-info
c     ------------------------------------------------------------------
      ELSE IF(Lextma)THEN
       fac=exp(Lndtcv/TWO/Dnefob)
       CALL scrmlt(fac,Na,A)
c     ------------------------------------------------------------------
      END IF
      RETURN
c     ------------------------------------------------------------------
 1010 FORMAT(/,' WARNING: ',a,' roots inside the unit circle',a)
 1020 FORMAT(/,' WARNING: Problem with MA parameter estimation.  ',a,
     &         ' can''t',
     &       /,'          invert the G''G matrix. Try a simpler ARIMA ',
     &         'model without',
     &       /,'          parameter constraints. Please send us the ',
     &         'data and spec file',
     &       /,'          that produced this message ',
     &         '(x12@census.gov)',a)
 1030 FORMAT(/,' WARNING: Problem calculating the theoretical ARMA ACF',
     &         a)
 1040 FORMAT(/,' WARNING: Problem calculating var(w_p|z)',a)
 1049 FORMAT('          for model ',a,'.  Will',/,
     &       '          attempt to fix the problem, and continue.')
 1050 FORMAT('          Will print out the parameters,',/,
     &       '          attempt to fix the problem, and continue.')
 1060 FORMAT(/)
      END

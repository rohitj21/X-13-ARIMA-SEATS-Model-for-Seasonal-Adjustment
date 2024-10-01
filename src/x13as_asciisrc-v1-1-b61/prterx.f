      SUBROUTINE prterx()
      IMPLICIT NONE
c     -------------------------------------------------------------------
c      If irregular regression matrix is singular, print out error 
c      message.
c     -------------------------------------------------------------------
      LOGICAL F
      PARAMETER(F=.false.)
c     -------------------------------------------------------------------
      INCLUDE 'stdio.i'
      INCLUDE 'srslen.prm'
      INCLUDE 'model.prm'
      INCLUDE 'model.cmn'
      INCLUDE 'mdldat.cmn'
      INCLUDE 'arima.cmn'
      INCLUDE 'units.cmn'
      INCLUDE 'error.cmn'
      INCLUDE 'tbllog.prm'
      INCLUDE 'tbllog.cmn'
      INCLUDE 'xrgtbl.i'
c     -------------------------------------------------------------------
      CHARACTER str*(PMDLCR)
      INTEGER nchr,nfil
c     -------------------------------------------------------------------
      INTEGER nblank
      EXTERNAL nblank
c-----------------------------------------------------------------------
      IF(Sngcol.lt.Ncxy)THEN
       CALL getstr(Colttl,Colptr,Ncoltl,Sngcol,str,nchr)
       IF(Lfatal)RETURN
      ELSE
       nchr=4
       str(1:nchr)='data'
      END IF
      CALL errhdr
      nfil=nblank(Cursrs)
      WRITE(STDERR,1230)Cursrs(1:nfil)
 1230 FORMAT(' Error(s) found while estimating the irregular ',
     &       'regression model.',/,
     &       ' For more details, check the error file (',a,'.err).')
      WRITE(Mt1,1270)str(1:nchr)
      WRITE(Mt2,1270)str(1:nchr)
 1270 FORMAT(/,' ERROR: Irregular regression matrix singular ',
     &         'because of ',a,'.',
     &       /,'        Check irregular regression model.',/)
c-----------------------------------------------------------------------
c     Added printing of regression matrix 
c-----------------------------------------------------------------------
      IF(Prttab(LXRXMX))THEN
       CALL prtshd('Irregular Component Regression Matrix',Begxy,Sp,
     &             Nrxy,F)
       IF(.not.Lfatal)CALL prtmtx(Begxy,Sp,Xy,Nrxy,Ncxy,Colttl,Colptr,
     &                            Ncoltl)
       IF(Lfatal)RETURN
      END IF
      CALL abend()
c-----------------------------------------------------------------------
      RETURN
      END

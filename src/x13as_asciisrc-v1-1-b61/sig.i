C 
C... Variables in Common Block /sig/ ... 
      integer SQG,HAR,FORTR,ILAM,NA,PG,NOADMISS,OUT,ITER,L,NDEC,
     $        BIAS,INOADMISS,RSA,NHTOFIX
      character TTLSET*80
      real*8 SQF,DOF,TIME,ZVAR,SEK,EPSPHI,RMOD,MAXBIAS
      common /sig/ TTLSET,SQF,DOF,TIME,ZVAR,SEK,EPSPHI,RMOD,MAXBIAS,
     $             SQG,HAR,FORTR,ILAM,NA,PG,NOADMISS,OUT,ITER,L,
     $             NDEC,BIAS,INOADMISS,RSA,NHTOFIX

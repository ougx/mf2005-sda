      SUBROUTINE IZREAD(IZONE,NLAY,NROW,NCOL,NZDIM,INZN1,INZN2,IOUT)
C     ******************************************************************
C     ROUTINE TO INPUT 3-D ZONE MATRIX, IZONE
C       INZN1 IS INPUT UNIT
C       IOUT IS OUTPUT UNIT
C     ******************************************************************
C        SPECIFICATIONS:
      DIMENSION IZONE(NCOL,NROW,NLAY)
      CHARACTER*20 FMTIN
      CHARACTER*80 NAME,NAMPRV,LINE
      CHARACTER*10 LOCAT
C     ------------------------------------------------------------------
      NAMPRV=' '
      DO 1000 K=1,NLAY
C
C-----READ ARRAY CONTROL RECORD.
      READ(INZN1,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,0,INZN1)
      LOCAT=LINE(ISTART:ISTOP)
C
C-----USE LOCAT TO SEE WHERE ARRAY VALUES COME FROM.
      IF(LOCAT.NE.'CONSTANT') GO TO 90
C
C-----LOCAT='CONSTANT' -- SET ALL ARRAY VALUES EQUAL TO ICONST.
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ICONST,R,0,INZN1)
      DO 80 I=1,NROW
      DO 80 J=1,NCOL
80    IZONE(J,I,K)=ICONST
      WRITE(IOUT,*)
      WRITE(IOUT,83) ICONST,K
83    FORMAT(13X,'Zone Array =',I4,' for layer',I4)
      IF(ICONST.LT.0) THEN
         WRITE(*,*) ' NEGATIVE ZONE NUMBER IS NOT ALLOWED'
         STOP
      END IF
      GO TO 1000
C
C-----Get FMTIN and IPRN -- there may be an unused value for ICONST
C-----in columns 11-20 if the file is an old file.
90    CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,I,R,0,INZN1)
      IF(LINE(ISTART:ISTART).NE.'(' ) THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,I,R,0,INZN1)
        IF(LINE(ISTART:ISTART).NE.'(' ) THEN
          WRITE(*,91) LINE
91        FORMAT(1X,
     1   'Format for reading zone array does not contain "(":',/1X,A)
          STOP
        END IF
      END IF
      FMTIN=LINE(ISTART:ISTOP)
C-----Blank inside parentheses indicates free format
      NC=ISTOP-ISTART-1
      IF(NC.LE.0) THEN
         FMTIN=' '
      ELSE
        IF(LINE(ISTART+1:ISTOP-1).EQ.' ') FMTIN=' '
      END IF
C-----Get print flag
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPRN,R,0,INZN1)
C
C-----LOCAT SHOULD BE EITHER 'EXTERNAL' OR 'INTERNAL'
C-----IF 'INTERNAL', READ ARRAY FROM SAME FILE
      IF(LOCAT.EQ.'INTERNAL') THEN
         INUNIT=INZN1
         WRITE(IOUT,*)
         WRITE(IOUT,92) K
92       FORMAT(1X,'Zone Array for layer',I4,
     1      ' will be read from the Zone File')
C
C-----IF 'EXTERNAL', OPEN A SEPARATE FILE
      ELSE IF(LOCAT.EQ.'EXTERNAL') THEN
         READ(INZN1,'(A)') NAME
         INUNIT=INZN2
         WRITE(IOUT,*)
         WRITE(IOUT,93) K,NAME
93       FORMAT(1X,'Zone Array for layer',I4,
     1         ' will be read from file:'/1X,A)
         IF(NAME.NE.NAMPRV) THEN
            IF(NAMPRV.NE.' ') CLOSE(UNIT=INUNIT)
            OPEN(UNIT=INUNIT,FILE=NAME,STATUS='OLD')
            WRITE(IOUT,96)
96          FORMAT(1X,'The file was opened successfully.')
            NAMPRV=NAME
         ELSE
            WRITE(IOUT,97)
97          FORMAT(1X,'This file is already open -- will continue readi
     1ng from the current location.')
         END IF
C
C-----LOCAT IS INVALID
      ELSE
         WRITE(*,*) ' INVALID LOCAT IN ARRAY CONTROL RECORD:',LOCAT
         STOP
      END IF
C
C-----LOCAT>0 -- READ RECORDS USING FREE-FORMAT OR FMTIN.
      IF(FMTIN.EQ.' ') THEN
         WRITE(IOUT,98) K
98       FORMAT(1X,'Zone Array for layer',I4,
     1       ' will be read using free format.'/1X,55('-'))
         DO 100 I=1,NROW
         READ(INUNIT,*) (IZONE(J,I,K),J=1,NCOL)
100      CONTINUE
      ELSE
         WRITE(IOUT,104) K,FMTIN
104      FORMAT(1X,'Zone Array for layer',I4,
     1       ' will be read using format: ',A/1X,71('-'))
         DO 110 I=1,NROW
         READ (INUNIT,FMTIN) (IZONE(J,I,K),J=1,NCOL)
110      CONTINUE
      END IF
C
C-----CHECK FOR NEGATIVE IZONE VALUES
320   DO 400 I=1,NROW
      DO 400 J=1,NCOL
      IF(IZONE(J,I,K).LT.0) THEN
         WRITE(*,*) ' NEGATIVE ZONE AT (LAYER,ROW,COLUMN):',K,I,J
         STOP
      END IF
400   CONTINUE
C
C-----IF PRINT CODE (IPRN) =>0 THEN PRINT ARRAY VALUES.
      IF(IPRN.LT.0) GO TO 1000
C
C-----PRINT COLUMN NUMBERS AT THE TOP OF THE PAGE
      WRITE(IOUT,421) (I,I=1,NCOL)
421   FORMAT(/,(5X,25I3))
      WRITE(IOUT,422)
422   FORMAT(1X,79('-'))
C
C-----PRINT EACH ROW IN THE ARRAY.
      DO 430 I=1,NROW
      WRITE(IOUT,423) I,(IZONE(J,I,K),J=1,NCOL)
423   FORMAT(1X,I3,1X,25I3/(5X,25I3))
430   CONTINUE
C
1000  CONTINUE
      IF(NAMPRV.NE.' ') CLOSE(UNIT=INZN2)
C
C-----RETURN
      RETURN
      END


      SUBROUTINE SDA_BCF7BDCH(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW FROM CONSTANT-HEAD CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,CR,CC,CV,
     1                      BOTM,LBOTM,IOUT
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,
     1                      ICHFLG
      USE GWFBCFMODULE,ONLY:IBCFCB,LAYCON
      use sda
C
      CHARACTER*16 TEXT
      DOUBLE PRECISION HD,CHIN,CHOUT,XX1,XX2,XX3,XX4,XX5,XX6,Q
C
      DATA TEXT /'   CONSTANT HEAD'/
C     ------------------------------------------------------------------
      !CALL SGWF2BCF7PNT(IGRID)
C
C1------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
      IBD=0
C
C2------CLEAR BUDGET ACCUMULATORS.
      ZERO=0.
      CHIN=ZERO
      CHOUT=ZERO
      IBDLBL=0
C
C3------CLEAR BUFFER.
      DO 5 K=1,NLAY
      DO 5 I=1,NROW
      DO 5 J=1,NCOL
      BUFF(J,I,K)=ZERO
5     CONTINUE
C
C3A-----IF SAVING CELL-BY-CELL FLOW IN A LIST, COUNT CONSTANT-HEAD
C3A-----CELLS AND WRITE HEADER RECORDS.
      IF(IBD.EQ.2) THEN
         NCH=0
         DO 7 K=1,NLAY
         DO 7 I=1,NROW
         DO 7 J=1,NCOL
         IF(IBOUND(J,I,K).LT.0) NCH=NCH+1
7        CONTINUE
         CALL UBDSV2(KSTP,KPER,TEXT,IBCFCB,NCOL,NROW,NLAY,
     1          NCH,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      END IF
C
C4------LOOP THROUGH EACH CELL AND CALCULATE FLOW INTO MODEL FROM EACH
C4------CONSTANT-HEAD CELL.
      DO 200 K=1,NLAY
      LC=LAYCON(K)
      DO 200 I=1,NROW
      DO 200 J=1,NCOL
C
C5------IF CELL IS NOT CONSTANT HEAD SKIP IT & GO ON TO NEXT CELL.
      IF (IBOUND(J,I,K).GE.0)GO TO 200
C
C6------CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
      X1=ZERO
      X2=ZERO
      X3=ZERO
      X4=ZERO
      X5=ZERO
      X6=ZERO
      CHCH1=ZERO
      CHCH2=ZERO
      CHCH3=ZERO
      CHCH4=ZERO
      CHCH5=ZERO
      CHCH6=ZERO
C
C7------CALCULATE FLOW THROUGH THE LEFT FACE.
C7------COMMENTS A-C APPEAR ONLY IN THE SECTION HEADED BY COMMENT 7,
C7------BUT THEY APPLY IN A SIMILAR MANNER TO SECTIONS 8-12.
C
C7A-----IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C7A-----TO NEXT FACE.  NO FLOW OCCURS AT THE EDGE OF THE GRID, TO AN
C7A-----ADJACENT NO-FLOW CELL, OR TO AN ADJACENT CONSTANT-HEAD CELL.
      IF(J.EQ.1) GO TO 30
      IF(IBOUND(J-1,I,K).EQ.0) GO TO 30
      IF(IBOUND(J-1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 30
C
C7B-----CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
      HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
      CHCH1=HDIFF*CR(J-1,I,K)
      IF(IBOUND(J-1,I,K).LT.0) GO TO 30
      X1=CHCH1
      XX1=X1
C
C7C-----ACCUMULATE POSITIVE AND NEGATIVE FLOW.
      IF (X1.LT.ZERO) THEN
        CHOUT=CHOUT-XX1
      ELSE
        CHIN=CHIN+XX1
      END IF
C
C8------CALCULATE FLOW THROUGH THE RIGHT FACE.
   30 IF(J.EQ.NCOL) GO TO 60
      IF(IBOUND(J+1,I,K).EQ.0) GO TO 60
      IF(IBOUND(J+1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 60
      HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
      CHCH2=HDIFF*CR(J,I,K)
      IF(IBOUND(J+1,I,K).LT.0) GO TO 60
      X2=CHCH2
      XX2=X2
      IF(X2.LT.ZERO) THEN
        CHOUT=CHOUT-XX2
      ELSE
        CHIN=CHIN+XX2
      END IF
C
C9------CALCULATE FLOW THROUGH THE BACK FACE.
   60 IF(I.EQ.1) GO TO 90
      IF(IBOUND(J,I-1,K).EQ.0) GO TO 90
      IF(IBOUND(J,I-1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 90
      HDIFF=HNEW(J,I,K)-HNEW(J,I-1,K)
      CHCH3=HDIFF*CC(J,I-1,K)
      IF(IBOUND(J,I-1,K).LT.0) GO TO 90
      X3=CHCH3
      XX3=X3
      IF(X3.LT.ZERO) THEN
        CHOUT=CHOUT-XX3
      ELSE
        CHIN=CHIN+XX3
      END IF
C
C10-----CALCULATE FLOW THROUGH THE FRONT FACE.
   90 IF(I.EQ.NROW) GO TO 120
      IF(IBOUND(J,I+1,K).EQ.0) GO TO 120
      IF(IBOUND(J,I+1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 120
      HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
      CHCH4=HDIFF*CC(J,I,K)
      IF(IBOUND(J,I+1,K).LT.0) GO TO 120
      X4=CHCH4
      XX4=X4
      IF(X4.LT.ZERO) THEN
        CHOUT=CHOUT-XX4
      ELSE
        CHIN=CHIN+XX4
      END IF
C
C11-----CALCULATE FLOW THROUGH THE UPPER FACE.
  120 IF(K.EQ.1) GO TO 150
      IF(IBOUND(J,I,K-1).EQ.0) GO TO 150
      IF(IBOUND(J,I,K-1).LT.0 .AND. ICHFLG.EQ.0) GO TO 150
      HD=HNEW(J,I,K)
      IF(LC.NE.3 .AND. LC.NE.2) GO TO 122
  122 HDIFF=HD-HNEW(J,I,K-1)
      CHCH5=HDIFF*CV(J,I,K-1)
      IF(IBOUND(J,I,K-1).LT.0) GO TO 150
      X5=CHCH5
      XX5=X5
      IF(X5.LT.ZERO) THEN
        CHOUT=CHOUT-XX5
      ELSE
        CHIN=CHIN+XX5
      END IF
C
C12-----CALCULATE FLOW THROUGH THE LOWER FACE.
  150 IF(K.EQ.NLAY) GO TO 180
      IF(IBOUND(J,I,K+1).EQ.0) GO TO 180
      IF(IBOUND(J,I,K+1).LT.0 .AND. ICHFLG.EQ.0) GO TO 180
      HD=HNEW(J,I,K+1)
      IF(LAYCON(K+1).NE.3 .AND. LAYCON(K+1).NE.2) GO TO 152
  152 HDIFF=HNEW(J,I,K)-HD
      CHCH6=HDIFF*CV(J,I,K)
      IF(IBOUND(J,I,K+1).LT.0) GO TO 180
      X6=CHCH6
      XX6=X6
      IF(X6.LT.ZERO) THEN
        CHOUT=CHOUT-XX6
      ELSE
        CHIN=CHIN+XX6
      END IF
C
C13-----SUM THE FLOWS THROUGH SIX FACES OF CONSTANT HEAD CELL, AND
C13-----STORE SUM IN BUFFER.
 180  RATE=CHCH1+CHCH2+CHCH3+CHCH4+CHCH5+CHCH6
      Q=dble(RATE)
      call SDA_SaveBD(iCBB,iSTEP,J,I,K,iBCHD,
     1 0.0d0,Q,NCOL,NROW,NLAY)
C
C14-----PRINT THE FLOW FOR THE CELL IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,899) TEXT,KPER,KSTP
  899    FORMAT(1X,/1X,A,'   PERIOD',I3,'   STEP',I3)
         WRITE(IOUT,900) K,I,J,RATE
  900    FORMAT(1X,'LAYER',I3,'   ROW',I4,'   COL',I4,
     1       '   RATE',1PG15.6)
         IBDLBL=1
      END IF
C
C15-----IF SAVING CELL-BY-CELL FLOW IN LIST, WRITE FLOW FOR CELL.
      IF(IBD.EQ.2) CALL UBDSVA(IBCFCB,NCOL,NROW,J,I,K,RATE,IBOUND,NLAY)
  200 CONTINUE
C
C16-----IF SAVING CELL-BY-CELL FLOW IN 3-D ARRAY, WRITE THE ARRAY.
      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,
     1                   IBCFCB,BUFF,NCOL,NROW,NLAY,IOUT)
C
C17-----SAVE TOTAL CONSTANT HEAD FLOWS AND VOLUMES IN VBVL TABLE
C17-----FOR INCLUSION IN BUDGET. PUT LABELS IN VBNM TABLE.
C
C18-----RETURN.
      RETURN
      END

      SUBROUTINE SDA_LPF7BDCH(KSTP,KPER)
C     ******************************************************************
C     COMPUTE FLOW FROM CONSTANT-HEAD CELLS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BUFF,CR,CC,CV,
     1                      BOTM,LBOTM,IOUT
      USE GWFBASMODULE,ONLY:MSUM,VBVL,VBNM,DELT,PERTIM,TOTIM,ICBCFL,
     1                      ICHFLG
      USE GWFLPFMODULE,ONLY:ILPFCB,LAYTYP,NOVFC
      use sda
      CHARACTER*16 TEXT
      DOUBLE PRECISION HD,CHIN,CHOUT,XX1,XX2,XX3,XX4,XX5,XX6,Q
C
      DATA TEXT /'   CONSTANT HEAD'/
C     ------------------------------------------------------------------
      !CALL SGWF2LPF7PNT(IGRID)
C
C1------SET IBD TO INDICATE IF CELL-BY-CELL BUDGET VALUES WILL BE SAVED.
      IBD=0
C
C2------CLEAR BUDGET ACCUMULATORS.
      ZERO=0.
      CHIN=ZERO
      CHOUT=ZERO
      IBDLBL=0
C
C3------CLEAR BUFFER.
      DO 5 K=1,NLAY
      DO 5 I=1,NROW
      DO 5 J=1,NCOL
      BUFF(J,I,K)=ZERO
5     CONTINUE
C
C3A-----IF SAVING CELL-BY-CELL FLOW IN A LIST, COUNT CONSTANT-HEAD
C3A-----CELLS AND WRITE HEADER RECORDS.
      IF(IBD.EQ.2) THEN
         NCH=0
         DO 7 K=1,NLAY
         DO 7 I=1,NROW
         DO 7 J=1,NCOL
         IF(IBOUND(J,I,K).LT.0) NCH=NCH+1
7        CONTINUE
         CALL UBDSV2(KSTP,KPER,TEXT,ILPFCB,NCOL,NROW,NLAY,
     1          NCH,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      END IF
C
C4------LOOP THROUGH EACH CELL AND CALCULATE FLOW INTO MODEL FROM EACH
C4------CONSTANT-HEAD CELL.
      DO 200 K=1,NLAY
      DO 200 I=1,NROW
      DO 200 J=1,NCOL
C
C5------IF CELL IS NOT CONSTANT HEAD SKIP IT & GO ON TO NEXT CELL.
      IF (IBOUND(J,I,K).GE.0)GO TO 200
C
C6------CLEAR VALUES FOR FLOW RATE THROUGH EACH FACE OF CELL.
      X1=ZERO
      X2=ZERO
      X3=ZERO
      X4=ZERO
      X5=ZERO
      X6=ZERO
      CHCH1=ZERO
      CHCH2=ZERO
      CHCH3=ZERO
      CHCH4=ZERO
      CHCH5=ZERO
      CHCH6=ZERO
C
C7------CALCULATE FLOW THROUGH THE LEFT FACE.
C7------COMMENTS A-C APPEAR ONLY IN THE SECTION HEADED BY COMMENT 7,
C7------BUT THEY APPLY IN A SIMILAR MANNER TO SECTIONS 8-12.
C
C7A-----IF THERE IS NO FLOW TO CALCULATE THROUGH THIS FACE, THEN GO ON
C7A-----TO NEXT FACE.  NO FLOW OCCURS AT THE EDGE OF THE GRID, TO AN
C7A-----ADJACENT NO-FLOW CELL, OR TO AN ADJACENT CONSTANT-HEAD CELL.
      IF(J.EQ.1) GO TO 30
      IF(IBOUND(J-1,I,K).EQ.0) GO TO 30
      IF(IBOUND(J-1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 30
C
C7B-----CALCULATE FLOW THROUGH THIS FACE INTO THE ADJACENT CELL.
      HDIFF=HNEW(J,I,K)-HNEW(J-1,I,K)
      CHCH1=HDIFF*CR(J-1,I,K)
      IF(IBOUND(J-1,I,K).LT.0) GO TO 30
      X1=CHCH1
      XX1=X1
C
C7C-----ACCUMULATE POSITIVE AND NEGATIVE FLOW.
      IF(X1.LT.ZERO) THEN
        CHOUT=CHOUT-XX1
      ELSE
        CHIN=CHIN+XX1
      END IF
C
C8------CALCULATE FLOW THROUGH THE RIGHT FACE.
   30 IF(J.EQ.NCOL) GO TO 60
      IF(IBOUND(J+1,I,K).EQ.0) GO TO 60
      IF(IBOUND(J+1,I,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 60
      HDIFF=HNEW(J,I,K)-HNEW(J+1,I,K)
      CHCH2=HDIFF*CR(J,I,K)
      IF(IBOUND(J+1,I,K).LT.0) GO TO 60
      X2=CHCH2
      XX2=X2
      IF(X2.LT.ZERO) THEN
        CHOUT=CHOUT-XX2
      ELSE
        CHIN=CHIN+XX2
      END IF
C
C9------CALCULATE FLOW THROUGH THE BACK FACE.
   60 IF(I.EQ.1) GO TO 90
      IF (IBOUND(J,I-1,K).EQ.0) GO TO 90
      IF (IBOUND(J,I-1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 90
      HDIFF=HNEW(J,I,K)-HNEW(J,I-1,K)
      CHCH3=HDIFF*CC(J,I-1,K)
      IF(IBOUND(J,I-1,K).LT.0) GO TO 90
      X3=CHCH3
      XX3=X3
      IF(X3.LT.ZERO) THEN
        CHOUT=CHOUT-XX3
      ELSE
        CHIN=CHIN+XX3
      END IF
C
C10-----CALCULATE FLOW THROUGH THE FRONT FACE.
   90 IF(I.EQ.NROW) GO TO 120
      IF(IBOUND(J,I+1,K).EQ.0) GO TO 120
      IF(IBOUND(J,I+1,K).LT.0 .AND. ICHFLG.EQ.0) GO TO 120
      HDIFF=HNEW(J,I,K)-HNEW(J,I+1,K)
      CHCH4=HDIFF*CC(J,I,K)
      IF(IBOUND(J,I+1,K).LT.0) GO TO 120
      X4=CHCH4
      XX4=X4
      IF(X4.LT.ZERO) THEN
        CHOUT=CHOUT-XX4
      ELSE
        CHIN=CHIN+XX4
      END IF
C
C11-----CALCULATE FLOW THROUGH THE UPPER FACE.
  120 IF(K.EQ.1) GO TO 150
      IF (IBOUND(J,I,K-1).EQ.0) GO TO 150
      IF (IBOUND(J,I,K-1).LT.0 .AND. ICHFLG.EQ.0) GO TO 150
      HD=HNEW(J,I,K)
      IF(NOVFC.NE.0 .OR. LAYTYP(K).EQ.0) GO TO 122
  122 HDIFF=HD-HNEW(J,I,K-1)
      CHCH5=HDIFF*CV(J,I,K-1)
      IF(IBOUND(J,I,K-1).LT.0) GO TO 150
      X5=CHCH5
      XX5=X5
      IF(X5.LT.ZERO) THEN
        CHOUT=CHOUT-XX5
      ELSE
        CHIN=CHIN+XX5
      END IF
C
C12-----CALCULATE FLOW THROUGH THE LOWER FACE.
  150 IF(K.EQ.NLAY) GO TO 180
      IF(IBOUND(J,I,K+1).EQ.0) GO TO 180
      IF(IBOUND(J,I,K+1).LT.0 .AND. ICHFLG.EQ.0) GO TO 180
      HD=HNEW(J,I,K+1)
      IF(NOVFC.NE.0 .OR. LAYTYP(K+1).EQ.0) GO TO 152
  152 HDIFF=HNEW(J,I,K)-HD
      CHCH6=HDIFF*CV(J,I,K)
      IF(IBOUND(J,I,K+1).LT.0) GO TO 180
      X6=CHCH6
      XX6=X6
      IF(X6.LT.ZERO) THEN
        CHOUT=CHOUT-XX6
      ELSE
        CHIN=CHIN+XX6
      END IF
C
C13-----SUM THE FLOWS THROUGH SIX FACES OF CONSTANT HEAD CELL, AND
C13-----STORE SUM IN BUFFER.
 180  RATE=CHCH1+CHCH2+CHCH3+CHCH4+CHCH5+CHCH6
      Q=dble(RATE)
      call SDA_SaveBD(iCBB,iSTEP,J,I,K,iBCHD,
     1 0.0d0,Q,NCOL,NROW,NLAY)
C
C14-----PRINT THE FLOW FOR THE CELL IF REQUESTED.
      IF(IBD.LT.0) THEN
         IF(IBDLBL.EQ.0) WRITE(IOUT,899) TEXT,KPER,KSTP
  899    FORMAT(1X,/1X,A,'   PERIOD ',I4,'   STEP ',I3)
         WRITE(IOUT,900) K,I,J,RATE
  900    FORMAT(1X,'LAYER ',I3,'   ROW ',I5,'   COL ',I5,
     1       '   RATE ',1PG15.6)
         IBDLBL=1
      END IF
C
C15-----IF SAVING CELL-BY-CELL FLOW IN LIST, WRITE FLOW FOR CELL.
      IF(IBD.EQ.2) CALL UBDSVA(ILPFCB,NCOL,NROW,J,I,K,RATE,IBOUND,NLAY)
  200 CONTINUE
C
C16-----IF SAVING CELL-BY-CELL FLOW IN 3-D ARRAY, WRITE THE ARRAY.
      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,
     1                   ILPFCB,BUFF,NCOL,NROW,NLAY,IOUT)
C
C
C18-----RETURN.
      RETURN
      END


      SUBROUTINE UBUDSVI(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,IOUT)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      INTEGER BUFF(NCOL,NROW,NLAY)
      INTEGER  KSTP,KPER,NCOL,NROW,NLAY
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD IDENTIFYING DATA.
      WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBUDSV SAVING "',A16,'" ON UNIT',I3,
     1     ' AT TIME STEP',I3,', STRESS PERIOD ',I4)
      WRITE(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,NLAY
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      WRITE(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END

      SUBROUTINE UBUDRDI(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,IOUT)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT,TEXT1
      INTEGER BUFF(NCOL,NROW,NLAY)
      INTEGER KKSTP,KKPER,KCOL,KROW,KLAY
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD IDENTIFYING DATA.
      WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBUDRD READING "',A16,'" ON UNIT',I3,
     1     ' AT TIME STEP',I3,', STRESS PERIOD ',I4)
      READ(IBDCHN) KKSTP,KKPER,TEXT1,KCOL,KROW,KLAY
      if (TEXT /= TEXT1) call USTOP("WRONG UBUDRDI")
      if (KKSTP /= KSTP) call USTOP("WRONG UBUDRDI")
      if (KKPER /= KPER) call USTOP("WRONG UBUDRDI")
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      READ(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END


      SUBROUTINE UBUDRD(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,IOUT)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT,TEXT1
      DIMENSION BUFF(NCOL,NROW,NLAY)
      INTEGER KKSTP,KKPER,KCOL,KROW,KLAY
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD IDENTIFYING DATA.
      WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBUDRD READING "',A16,'" ON UNIT',I3,
     1     ' AT TIME STEP',I3,', STRESS PERIOD ',I4)
      READ(IBDCHN) KKSTP,KKPER,TEXT1,KCOL,KROW,KLAY
      if (TEXT /= TEXT1) call USTOP("WRONG UBUDRD")
      if (KKSTP /= KSTP) call USTOP("WRONG UBUDRD")
      if (KKPER /= KPER) call USTOP("WRONG UBUDRD")
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      READ(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END

      subroutine SDA_CalcErr(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,RES,
     1   NCOL,NROW,NLAY,NODES)
      DOUBLE PRECISION HNEW, DZERO,HHCOF,RRHS
      DOUBLE PRECISION Z,B,D,E,F,H,S
      DOUBLE PRECISION ZHNEW,BHNEW,DHNEW,FHNEW,HHNEW,SHNEW
C
      DIMENSION HNEW(NODES), IBOUND(NODES), CR(NODES), CC(NODES),
     1  CV(NODES), HCOF(NODES), RHS(NODES), RES(NODES)
      integer ::  K,I,J,KK,JJ,II

      RES = 0.
      NRC = NCOL*NROW
      DO 150 K=1,NLAY
      DO 150 I=1,NROW
      DO 150 J=1,NCOL
C
C6A-----SET UP CURRENT CELL LOCATION INDEXES.  THESE ARE DEPENDENT
C6A-----ON THE DIRECTION OF EQUATION ORDERING.
      II=I
      JJ=J
      KK=K
C
C6B-----CALCULATE 1 DIMENSIONAL SUBSCRIPT OF CURRENT CELL AND
C6B-----SKIP CALCULATIONS IF CELL IS NOFLOW OR CONSTANT HEAD
  122 N=JJ+(II-1)*NCOL+(KK-1)*NRC
      IF(IBOUND(N) < 0) GO TO 150
C
C6C-----CALCULATE 1 DIMENSIONAL SUBSCRIPTS FOR LOCATING THE 6
C6C-----SURROUNDING CELLS
      NRN=N+NCOL
      NRL=N-NCOL
      NCN=N+1
      NCL=N-1
      NLN=N+NRC
      NLL=N-NRC
C
C6D-----CALCULATE 1 DIMENSIONAL SUBSCRIPTS FOR CONDUCTANCE TO THE 6
C6D-----SURROUNDING CELLS.  THESE DEPEND ON ORDERING OF EQUATIONS.
      NCF=N
      NCD=NCL
      NRB=NRL
      NRH=N
      NLS=N
      NLZ=NLL
C
C6E-----ASSIGN VARIABLES IN MATRICES A & U INVOLVING ADJACENT CELLS
C6E1----NEIGHBOR IS 1 ROW BACK
  126 B=DZERO
      BHNEW=DZERO
      IF(I.EQ.1) GO TO 128
      B=CC(NRB)
      BHNEW=B*HNEW(NRL)
C
C6E2----NEIGHBOR IS 1 ROW AHEAD
  128 H=DZERO
      HHNEW=DZERO
      IF(I.EQ.NROW) GO TO 130
      H=CC(NRH)
      HHNEW=H*HNEW(NRN)
C
C6E3----NEIGHBOR IS 1 COLUMN BACK
  130 D=DZERO
      DHNEW=DZERO
      IF(J.EQ.1) GO TO 132
      D=CR(NCD)
      DHNEW=D*HNEW(NCL)
C
C6E4----NEIGHBOR IS 1 COLUMN AHEAD
  132 F=DZERO
      FHNEW=DZERO
      IF(J.EQ.NCOL) GO TO 134
      F=CR(NCF)
      FHNEW=F*HNEW(NCN)
C
C6E5----NEIGHBOR IS 1 LAYER BEHIND
  134 Z=DZERO
      ZHNEW=DZERO
      IF(K.EQ.1) GO TO 136
      Z=CV(NLZ)
      ZHNEW=Z*HNEW(NLL)
C
C6E6----NEIGHBOR IS 1 LAYER AHEAD
  136 S=DZERO
      SHNEW=DZERO
      IF(K.EQ.NLAY) GO TO 138
      S=CV(NLS)
      SHNEW=S*HNEW(NLN)
C
C6E7----CALCULATE THE NEGATIVE SUM OF ALL CONDUCTANCES TO NEIGHBORING
C6E7----CELLS
  138 E=-Z-B-D-F-H-S
C
      HHCOF=HCOF(N)

C
C6G-----CALCULATE THE RESIDUAL
      RRHS=RHS(N)
      RES(N)=RRHS-ZHNEW-BHNEW-DHNEW-E*HNEW(N)-HHCOF*HNEW(N)-FHNEW-HHNEW
     1      -SHNEW
C
C
  150 CONTINUE
      end subroutine

      SUBROUTINE TIMEDIFF(IBDT,IEDT,SEC)
C     ******************************************************************
C     Get the elapsed time since IBDT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: IBDT(8), IEDT(8)
      INTEGER IDPM(12)
      real SEC
      DATA IDPM/31,28,31,30,31,30,31,31,30,31,30,31/ ! Days per month
      DATA NSPD/86400/  ! Seconds per day
C     ------------------------------------------------------------------
C
C     Get current date and time, assign to IEDT, and write.
      !CALL DATE_AND_TIME(VALUES=IEDT)

C
C     Calculate elapsed time in days and seconds
      if(IBDT(1)/=IEDT(1).or.IBDT(2)/=IEDT(2).or.IBDT(3)/=IEDT(3)) then
        NDAYS=0
        LEAP=0
        IF (MOD(IEDT(1),4).EQ.0) LEAP = 1
        IBD = IBDT(3)            ! BEGIN DAY
        IED = IEDT(3)            ! END DAY
C     FIND DAYS
        IF (IBDT(2).NE.IEDT(2)) THEN
C       MONTHS DIFFER
          MB = IBDT(2)             ! BEGIN MONTH
          ME = IEDT(2)             ! END MONTH
          NM = ME-MB+1             ! NUMBER OF MONTHS TO LOOK AT
          IF (MB.GT.ME) NM = NM+12
          MC=MB-1
          DO 10 M=1,NM
            MC=MC+1                ! MC IS CURRENT MONTH
            IF (MC.EQ.13) MC = 1
            IF (MC.EQ.MB) THEN
              NDAYS = NDAYS+IDPM(MC)-IBD
              IF (MC.EQ.2) NDAYS = NDAYS + LEAP
            ELSEIF (MC.EQ.ME) THEN
              NDAYS = NDAYS+IED
            ELSE
              NDAYS = NDAYS+IDPM(MC)
              IF (MC.EQ.2) NDAYS = NDAYS + LEAP
            ENDIF
10          CONTINUE
        ELSEIF (IBD.LT.IED) THEN
C       START AND END IN SAME MONTH, ONLY ACCOUNT FOR DAYS
          NDAYS = IED-IBD
        ENDIF
        ELSEC=NDAYS*NSPD
      else
        ELSEC = 0.0
      endif
C
C     ADD OR SUBTRACT SECONDS
      ELSEC = ELSEC+(IEDT(5)-IBDT(5))*3600.0
      ELSEC = ELSEC+(IEDT(6)-IBDT(6))*60.0
      ELSEC = ELSEC+(IEDT(7)-IBDT(7))
      ELSEC = ELSEC+(IEDT(8)-IBDT(8))*0.001

      SEC = ELSEC

      RETURN
      END

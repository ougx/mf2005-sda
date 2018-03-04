!  Code in this file defines values for OPEN-statement specifiers.  Some
!  of the values are extensions to ANSI Fortran 90 and 95.  One of the
!  specifiers is not included in ANSI FORTRAN 77. The included
!  specifiers are ACCESS, FORM and ACTION.
!
      CHARACTER*10 :: ACCESS='STREAM    '
      CHARACTER*12 :: FORM=  'UNFORMATTED '
      CHARACTER*10 :: ACTION(2)=(/'READ      ','READWRITE '/)

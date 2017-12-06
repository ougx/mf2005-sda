module SDA
  !use GLOBAL, only                  :   DEFAULT_KIND
  logical                           ::  debug = .FALSE.
  integer                           ::  idebug                  !output file number for scenario runs (2D head differebces)
  integer                           ::  IGRID = 1               !fake sub grid number

  doubleprecision                   ::  DTRead                  !Total time used for reading files
  doubleprecision                   ::  DTSolver                !Total time used to solve equation
  doubleprecision                   ::  DTStpIni                !Total time used to solve equation
  doubleprecision                   ::  DTBD                    !Total time used to finalize the scenario run
  doubleprecision                   ::  DTStpEnd                !Total time used for step calculations


  integer                           ::  NSEN = 0                !total number of scenario runs, > 0 no baseline run; < 0 include baseline run
  integer                           ::  QWEL                    !flag if the well package is activated in the scenario runs, 1 is activated, other is none
  integer                           ::  QRCH                    !flag if the recharge package is activated in the scenario runs, 1 is activated, other is none
  integer                           ::  QHED                    !flag for writting head difference, if QHED > 0, write as formatted, QHED < 0, write as unformated, QHED = 0, no write
  integer                           ::  QCCB                    !flag for writting cell by cell budgets


  integer                           ::  mCell                   !number of active cells


  integer                           ::  NCBMAX                  !maximum number of boundary
  integer                           ::  NCB                     !number of boundary each time step

  integer                           ::  iSTEP                   !current time step
  integer                           ::  NSTEP                   !total number of time steps


  integer                           ::  iTypeBC
  integer                           ::  nTypeBC




  real,allocatable    ::  SDACR(:,:,:)            !baseline CR
  real,allocatable    ::  SDACC(:,:,:)            !baseline CC
  real,allocatable    ::  SDACV(:,:,:)            !baseline CV
  real,allocatable       ::  SDAHCOF(:,:,:)          !baseline HCOF
  !real,allocatable,dimension(:,:,:) ::  SDARHS                 !baseline RHS
  real,allocatable,dimension(:,:,:) ::  SDASC0                  !baseline Storage Capacity0
  real,allocatable,dimension(:,:,:) ::  SDASC                   !baseline Storage Capacity
  integer,allocatable,dimension(:,:,:) ::  SDAIBD               !baseline IBOUND




  integer*1                         ::  CBINDX(100)             !the order of CB in outputs
  integer                           ::  CBNTOT(100)             !total number of CB during the simulation
  character(len=16)                 ::  CBTYPE(100)             !CB name

  doubleprecision                   ::  CBVBVL(4,100)           !In&Out, Rate and Cumulative


  integer                           ::  NBZONEMax               !maximum number of boundary ZONEs for all scenarios
  integer                           ::  NBZONE                  !number of boundary ZONEs of a scenario run
  integer,allocatable               ::  BZONE(:,:,:)            !boundary ZONE, defined at NewScen
  doubleprecision,allocatable       ::  RBZONE(:,:,:)           !flow rate of each zone, Calculated each NewStp


  integer,allocatable               ::  iCR(:)                  !output file number for baseline CR, positive is unformatted, negative is formatted
  integer,allocatable               ::  iCC(:)                  !output file number for baseline CC, positive is unformatted, negative is formatted
  integer,allocatable               ::  iCV(:)                  !output file number for baseline CV, positive is unformatted, negative is formatted
  integer,allocatable               ::  iHCOF(:)                !output file number for baseline HCOF, positive is unformatted, negative is formatted
  integer,allocatable               ::  iRHS(:)                 !output file number for baseline RHS, positive is unformatted, negative is formatted
  integer,allocatable               ::  iSC(:)                  !output file number for baseline SC, positive is unformatted, negative is formatted
  integer,allocatable               ::  iCB(:)                  !output file number for baseline boundary condutance
  integer,allocatable               ::  iACTIVE(:)              !output file number for IBOUND
  integer,allocatable               ::  iERR(:)                 !output file number for baseline residual


  integer                           ::  iHED =0                 !output file number for scenario runs (2D head differebces)

  character*80                      ::  ScenName                !Scenario Name
  integer                           ::  iSen                    !index of scenario run

  integer                           ::  iCBB =  0               !output file number for scenario runs

  include '../mf-src/openspec.inc'

  !cell of head-dependent flow boundary condition
  type  ::  SDAB
    integer                           ::  ICOL                  !col
    integer                           ::  IROW                  !row
    integer                           ::  ILAY                  !layer

    real                              ::  CB                    !conductance
    real                              ::  HB                    !boundary head
    integer                           ::  BND                   !boundary type index
    type(SDAB), pointer               ::  next
  end type



  type(SDAB),    allocatable          ::  pSDACB(:)
  !type(avLayer), allocatable          ::  avLAY(:)
end module

subroutine SDA_AddCB(IC,IR,IL,TyB,CB,HB)

  use sda,only  : NSEN, NCB, pSDACB
  integer :: IC,IR,IL,TyB
  real :: CB,HB

  if (NSEN == 0) return
  if (CB == 0.) return
  NCB=NCB+1

!  if (K>14) call USTOP(BTy)
!
!  if (BTy == ' EVT' .or. BTy == ' ETS') then
!    SDACB(NCB)%CB = -CB
!  else
!    SDACB(NCB)%CB = CB
!  endif
  pSDACB(NCB)%CB = CB
  pSDACB(NCB)%ICOL = IC
  pSDACB(NCB)%IROW = IR
  pSDACB(NCB)%ILAY = IL
  pSDACB(NCB)%BND = TyB
  pSDACB(NCB)%HB = HB


end subroutine

subroutine SDA_AddSC(IC,IR,IL,SC)

  !SC is the Specific yield or Ss * Thickness
  use sda,only  : SDASC, NSEN
  integer :: IC,IR,IL
  real :: SC

  if (NSEN>=0) return

!  if (IC/=0 .and. IR/=0) then
    SDASC(IC,IR,IL) = SC
!  elseif (IC = 0 .and. IR = 0) then
!    SDASC(:,:,IL) = SC(:,:,IL)
!  elseif (IC = 0) then
!    SDASC(:,IR,IL) = SC(:,IR,IL)
!  endif

end subroutine

subroutine SDA_Ini(iin)
  !initialize
  use sda
  use GLOBAL,only                         :   NSTP,NPER,NCOL,NROW,NLAY,IUNIT,IOUT,IBOUND
  use GWFBASMODULE,only                   :   MSUM,VBNM
  integer                                 ::  iin                   !input file number of SDA main file
  integer                                 ::  iset, nset            !input file number of SDA main file
  integer                                 ::  iper, iper0           !input file number of SDA main file
  character*400                           ::  sline                 !input file number of SDA main file
  character*100                           ::  fhcof
  integer                                 ::  i1,i2,i3, itmp, iloc  !input file number of SDA main file
  real                                    ::  rtmp                  !input file number of SDA main file

  !include 'openspec.inc'

  IGRID = 1

  allocate(iCR(NPER))                   !output file number for baseline CR, positive is unformatted, negative is formatted
  allocate(iCC(NPER))                   !output file number for baseline CC, positive is unformatted, negative is formatted
  allocate(iCV(NPER))                   !output file number for baseline CV, positive is unformatted, negative is formatted
  allocate(iHCOF(NPER))                 !output file number for baseline HCOF, positive is unformatted, negative is formatted
  !allocate(iRHS(1))                  !output file number for baseline RHS, positive is unformatted, negative is formatted
  allocate(iSC(NPER))                   !output file number for baseline SC, positive is unformatted, negative is formatted
  allocate(iCB(NPER))                   !output file number for baseline boundary condutance
  !allocate(iERR(1))                  !output file number for baseline residual
  allocate(iACTIVE(NPER))               !output file number for baseline IBOUND




  iper0 = 1
  NBZONEMax = 99

  iloc = 1
  CALL URDCOM(iin,IOUT,sline)
  call URWORD(sline,iloc,i1,i2,2,NSEN,rtmp,IOUT,iin)
  call URWORD(sline,iloc,i1,i2,2,NCBMAX,rtmp,IOUT,iin)
  call URWORD(sline,iloc,i1,i2,2,nset,rtmp,IOUT,iin)


  allocate(pSDACB(NCBMAX))

  !write(IOUT,*) " READING UNIT NUMBERS FOR SDA ... "

  do iset = 1, nset

    write(IOUT,*) " READING UNIT SET ", iset
    read(iin, "(A)") sline
    if (iper0>NPER) cycle
    iloc = 1
    call URWORD(sline,iloc,i1,i2,2,iper,rtmp,IOUT,iin)
    if (iper>NPER) then
      iper = NPER
    endif

    call URWORD(sline,iloc,i1,i2,2,iACTIVE(iper),rtmp,IOUT,iin)
    call URWORD(sline,iloc,i1,i2,2,iCC(iper),rtmp,IOUT,iin)
    call URWORD(sline,iloc,i1,i2,2,iCR(iper),rtmp,IOUT,iin)
    call URWORD(sline,iloc,i1,i2,2,iCV(iper),rtmp,IOUT,iin)
    call URWORD(sline,iloc,i1,i2,2,iSC(iper),rtmp,IOUT,iin)
    call URWORD(sline,iloc,i1,i2,2,iHCOF(iper),rtmp,IOUT,iin)
    call URWORD(sline,iloc,i1,i2,2,iCB(iper),rtmp,IOUT,iin)
!    call URWORD(sline,iloc,i1,i2,2,iERR(iper),rtmp,IOUT,iin)
!
!    if (iERR(iper)==iCB(iper) .or. iERR(iper)==iSC(iper) .or. iERR(iper)==iHCOF(iper) .or. &
!        iERR(iper)==iCV(iper) .or. iERR(iper)==iCR(iper) .or. iERR(iper)==iACTIVE(iper)) then
!      write(IOUT,*) "Residuals needed to be saved to a different file."
!      call USTOP('Residuals needed to be saved to a different file.')
!    endif

    iACTIVE(iper0:(iper-1)) = iACTIVE(iper)
    iCR(iper0:(iper-1)) = iCR(iper)
    iCC(iper0:(iper-1)) = iCC(iper)
    iCV(iper0:(iper-1)) = iCV(iper)
    iHCOF(iper0:(iper-1)) = iHCOF(iper)
    iSC(iper0:(iper-1)) = iSC(iper)
    iCB(iper0:(iper-1)) = iCB(iper)
    !iERR(iper0:(iper-1)) = iERR(iper)

    iper0 = iper + 1
  enddo

  if (debug) open(newunit=idebug, file="debug.txt")


  !allocate(NCBStp(NSTEP))
  !open(newunit=iCB(1),file='SDACB.dat',form=form,access=access)

  !store the active cell indices of each layer
  !the purpose is to reduce the size of the coefficient files


  if (NSEN > 0) then
    CALL SDA_Run(iin)
    CALL USTOP(' ')
  else

    allocate(SDAIBD(NCOL,NROW,NLAY))
    allocate(SDACV(NCOL,NROW,NLAY))
    allocate(SDACC(NCOL,NROW,NLAY))
    allocate(SDACR(NCOL,NROW,NLAY))
    allocate(SDAHCOF(NCOL,NROW,NLAY))
    allocate(SDASC0(NCOL,NROW,NLAY))
    allocate(SDASC(NCOL,NROW,NLAY))

    SDAIBD = 0
    SDASC = 0.
    SDACV = 0.
    SDACC = 0.
    SDACR = 0.
    SDAHCOF = 0.
    SDASC0 = 0.

    ! save number of initial active cells in ibound
!    do i1 = 1, NLAY
!      write(iHCOF(i1)) avLAY(i1)%nActive
!    end do
  end if

end subroutine

subroutine SDA_IterIni
  !baseline iteration start
  use sda


  !SDASC = 0.
  NCB = 0
  !pSDACB => pSDACB1

end subroutine


subroutine SDA_IterEnd(KSTP,KPER)
  !baseline time step end

  !call SDA_CalcErr

  call SDA_SaveCB(KSTP,KPER)
  call SDA_SaveHCOF(KSTP,KPER)

end subroutine

SUBROUTINE SDA_SaveCB(KSTP,KPER)
  !save head-dependent boundary for the current time step
  use sda
  use GLOBAL, only : IOUT
  integer                           ::  i
  integer                           ::  KSTP                  !col
  integer                           ::  KPER                  !row

  write(IOUT,'(//x,A,I0,A)') 'SDA -- SAVING ', NCB, ' HEAD DEPENDENT BOUNDARY CONDITIONS '


  WRITE(iCB(KPER)) KSTP,KPER,NCB
  do i=1, NCB
    write(iCB(KPER)) pSDACB(i)%ILAY,pSDACB(i)%IROW,pSDACB(i)%ICOL,pSDACB(i)%BND,pSDACB(i)%CB
  enddo
END

SUBROUTINE SDA_ReadCB(KSTP,KPER)
  !read head-dependent boundary for the current time step
  use sda

  integer                           ::  KSTP,KKSTP                  !step
  integer                           ::  KPER,KKPER                  !period
  integer                           ::  i                     !index

  read(iCB(KPER)) KKSTP,KKPER,NCB

  if (KKSTP /= KSTP .or. KKPER /= KPER) call USTOP('Wrong CB')

  do i=1, NCB
    read(iCB(KPER)) pSDACB(i)%ILAY,pSDACB(i)%IROW,pSDACB(i)%ICOL,pSDACB(i)%BND,pSDACB(i)%CB
  enddo

END


SUBROUTINE SDA_ReadALLCB()
  !read head-dependent boundary for the current time step
  use sda
  use GLOBAL,only                   :   NPER,NSTP,IUNIT
  use GWFBASMODULE,only             :   MSUM,VBNM

  integer                           ::  KSTP,KKSTP                  !step
  integer                           ::  KPER,KKPER                  !period
  integer                           ::  i,IL,IR,IC,IB                     !index
  real                              ::  CB

  integer*2                         ::  iidx(100), max_term


  max_term = 29



  MSUM = 0

  !          1234567890123456   1234567890123456   1234567890123456   1234567890123456   1234567890123456

  CBTYPE =(/'STORAGE         ','CONSTANT HEAD   ','WELLS           ','RECHARGE        ','DRAINS          ', & !5
            'RIVER LEAKAGE   ','STR LEAKAGE     ','SFR LEAKAGE     ','HEAD DEP BOUNDS ','ET              ', & !10
            'FARM WELLS      ','FARM  NET  RECH ','GHOST-NODE FLUX ','DRAINS_DRT      ','ET SEGMENTS     ', & !15
            'IBS STORAGE     ','LAKE  SEEPAGE   ','MNW             ','MNW2            ','RESERV. LEAKAGE ', & !20
            'RIPARIAN_ET     ','INST. IB STORAGE','SWR LEAKAGE     ','SWR GWET        ','SWT STORAGE     ', & !25
            'UPW STORAGE     ','SURFACE LEAKAGE ','UZF RECHARGE    ','GW ET           ','                ', & !30
           ('                ',i1=31,100)/)



  iidx = (/ 66,   66,   2,    8,    3,  & !STO, CHD, WEL, RCH, DRN
            4,    18,   44,   7,    5,  & !RIV, STR, SFR, GHB, EVT
            0,    0,    0,    40,   39, & !FMP, FMP, BFH, DRT, ETS
            19,   22,   52,   50,   17, & !IBS, LAK, MNW1,MNW2,RES
            0,    54,   64,   64,   57, & !RIP, SUB, SWR, SWR, SWT
            0,    55,   55,   55,   0,  & !UPW, UZF
            (-i,i=31,100) &
          /)

  CBNTOT = 0
  !CBNTOT(1:max_term) = IUNIT(iidx(1:max_term))
  where (iidx > 0)
    CBNTOT = IUNIT(iidx)
  end where

!  rewind(iCB(1))
!  do KPER=1, NPER
!    do KSTP=1,NSTP(KPER)
!      read(iCB(1)) KKSTP,KKPER,NCB
!
!      if (KKSTP /= KSTP .or. KKPER /= KPER) call USTOP('Wrong CB')
!      !if (abs(DT-DELT) >= 1.e-7) call USTOP('Wrong delt for CB')
!
!      do i=1, NCB
!        read(iCB(1)) IL,IR,IC,IB,CB
!        CBNTOT(IB) = CBNTOT(IB) + 1
!      enddo
!    end do
!  end do


  i = 0
  CBINDX = 0
  VBNM = ''
  do IB = 1, max_term
    if (CBNTOT(IB) > 0) then
      i = i + 1
      VBNM(i) = CBTYPE(IB)
      CBINDX(IB) = i
    end if
  end do

  MSUM = i

END



subroutine SDA_SaveHCOF(KSTP,KPER)
  !save flow coefficients
  use sda
  use GLOBAL,only                         :   CR,CC,CV,HCOF,IOUT,IBOUND,RHS,NCOL,NROW,NLAY

  CHARACTER*16 TEXT
  integer ::  K,KKPER,KKSTP,KSTP,KPER
  character*1 ::  newco

  !       1234567890123456
  !KKPER = min(1,KPER)
  DO K=1,NLAY

    if ((KPER+KSTP) == 2 .or. any(IBOUND(:,:,K) /= SDAIBD(:,:,K))) then
      newco = 'N'
      SDAIBD(:,:,K) = IBOUND(:,:,K)
    else
      newco = 'O'
    endif
    write(iACTIVE(KPER)) newco
    if (newco == 'N') then
      write(TEXT,"(A6,A4,I3,3x)") 'IBOUND','LAY',K
      write(iACTIVE(KPER)) IBOUND(:,:,K)
      !call UBUDSVI(KSTP,KPER,TEXT,iACTIVE(KPER),IBOUND(:,:,K),NCOL,NROW,1,IOUT)
    endif

    if ((KPER+KSTP) == 2 .or. any(CR(:,:,K) /= SDACR(:,:,K))) then
      newco = 'N'
      SDACR(:,:,K) = CR(:,:,K)
    else
      newco = 'O'
    endif
    write(iCR(KPER)) newco
    if (newco == 'N') then
      !       1234567890123456
      write(TEXT,"(A6,A4,I3,3x)") 'CR','LAY',K
      write(iCR(KPER)) CR(:,:,K)
      !call UBUDSV(KSTP,KPER,TEXT,iCR(KPER),CR(:,:,K),NCOL,NROW,1,IOUT)
    endif

    if ((KPER+KSTP) == 2 .or. any(CC(:,:,K) /= SDACC(:,:,K))) then
      newco = 'N'
      SDACC(:,:,K) = CC(:,:,K)
    else
      newco = 'O'
    endif
    write(iCC(KPER)) newco
    if (newco == 'N') then
      !       1234567890123456
      write(TEXT,"(A6,A4,I3,3x)") 'CC','LAY',K
      write(iCC(KPER)) CC(:,:,K)
      !call UBUDSV(KSTP,KPER,TEXT,iCC(KPER),CC(:,:,K),NCOL,NROW,1,IOUT)
    endif


    if ((KPER+KSTP) == 2 .or. any(CV(:,:,K) /= SDACV(:,:,K))) then
      newco = 'N'
      SDACV(:,:,K) = CV(:,:,K)
    else
      newco = 'O'
    endif
    write(iCV(KPER)) newco
    if (newco == 'N') then
      !       1234567890123456
      write(TEXT,"(A6,A4,I3,3x)") 'CV','LAY',K
      write(iCV(KPER)) CV(:,:,K)
      !call UBUDSV(KSTP,KPER,TEXT,iCV(KPER),CV(:,:,K),NCOL,NROW,1,IOUT)
    endif


    if ((KPER+KSTP) == 2 .or. any(SDASC(:,:,K) /= SDASC0(:,:,K))) then
      newco = 'N'
      SDASC0(:,:,K) = SDASC(:,:,K)
    else
      newco = 'O'
    endif
    write(iSC(KPER)) newco
    if (newco == 'N') then
      !       1234567890123456
      write(TEXT,"(A6,A4,I3,3x)") 'SC','LAY',K
      write(iSC(KPER)) SDASC(:,:,K)
      !call UBUDSV(KSTP,KPER,TEXT,iSC(KPER),SDASC(:,:,K),NCOL,NROW,1,IOUT)
    endif

  enddo
  !check if anything change, if yes, save, otherwise skip

  !if (debug) write(idebug, "(10(1PE15.7))") RHS
endsubroutine

subroutine SDA_ReadHCOF(KPER,KSTP)
  !read flow coefficients
  use sda
  use GLOBAL,only                         :   CR,CC,CV,HCOF,IOUT,IBOUND,NCOL,NROW,NLAY
  use GWFBASMODULE,only                   :   DELT
  integer ::  KSTP,KPER

  CHARACTER*16 :: TEXT,TEXT1,TEXT2,TEXT3,TEXT4,TEXT5,TEXT6

  integer :: K
  character*1 :: newco


  DO K=1,NLAY

    write(TEXT,"(A6,A4,I3,3x)") 'IBOUND','LAY',K
    read(iACTIVE(KPER)) newco
    if (newco == "N") then
      !call UBUDRDI(KSTP,KPER,TEXT,iACTIVE(KPER),IBOUND(:,:,K),NCOL,NROW,1,IOUT)
      read(iACTIVE(KPER)) IBOUND(:,:,K)
    else
      if (newco /= "O") call USTOP(TEXT)
    endif


    write(TEXT,"(A6,A4,I3,3x)") 'CR','LAY',K
    read(iCR(KPER)) newco
    if (newco == "N") then
      !       1234567890123456
      !call UBUDRD(KSTP,KPER,TEXT,iCR(KPER),CR(:,:,K),NCOL,NROW,1,IOUT)
      read(iCR(KPER)) CR(:,:,K)
    else
      if (newco /= "O") call USTOP(TEXT)
    endif


    write(TEXT,"(A6,A4,I3,3x)") 'CC','LAY',K
    read(iCC(KPER)) newco
    if (newco == "N") then
      !       1234567890123456
      !call UBUDRD(KSTP,KPER,TEXT,iCC(KPER),CC(:,:,K),NCOL,NROW,1,IOUT)
      read(iCC(KPER)) CC(:,:,K)
    else
      if (newco /= "O") call USTOP(TEXT)
    endif



    write(TEXT,"(A6,A4,I3,3x)") 'CV','LAY',K
    read(iCV(KPER)) newco
    if (newco == "N") then
      !       1234567890123456
      !call UBUDRD(KSTP,KPER,TEXT,iCV(KPER),CV(:,:,K),NCOL,NROW,1,IOUT)
      read(iCV(KPER)) CV(:,:,K)
    else
      if (newco /= "O") call USTOP(TEXT)
    endif



    write(TEXT,"(A6,A4,I3,3x)") 'SC','LAY',K
    read(iSC(KPER)) newco
    if (newco == "N") then
      !       1234567890123456
      !call UBUDRD(KSTP,KPER,TEXT,iSC(KPER),SDASC0(:,:,K),NCOL,NROW,1,IOUT)
      read(iSC(KPER)) SDASC0(:,:,K)
    else
      if (newco /= "O") call USTOP(TEXT)
    endif


    !call SDA_TSPCoeff
  enddo

  HCOF = -SDASC0 / DELT
  do k = 1, NCB
    HCOF(pSDACB(k)%ICOL, pSDACB(k)%IROW, pSDACB(k)%ILAY) = HCOF(pSDACB(k)%ICOL, pSDACB(k)%IROW, pSDACB(k)%ILAY) - pSDACB(k)%CB
  end do

endsubroutine



subroutine SDA_Run(iin)
  !start SDA
  !baseline run can be skipped. coefficients are read from the the file directly
  use sda
  use GLOBAL
  use GWFBASMODULE

  implicit none

  integer                                 ::  iin               !input file number of SDA main file
  integer                                 ::  k ,i1,i2,iloc     !layer index
  character*400                           ::  line              !layer index


  WRITE(*,*) ' Entering scenario analysis mode'


  if (NSEN<0) then

    ! release the array from baseline run
    deallocate(SDAIBD)
    deallocate(SDACV)
    deallocate(SDACC)
    deallocate(SDACR)
    deallocate(SDAHCOF)

  else

    allocate(SDASC(NCOL,NROW,NLAY))
    allocate(SDASC0(NCOL,NROW,NLAY))

    !close(iACTIVE,iCR,iCC,iCV,iHCOF,iSC,iCB)

  endif


  !free format
  IFREFM = 1

  ! read the flow coefficients
  ! call SDA_ReadHCOF()


  !zones for zone budget
  allocate(BZONE(NCOL,NROW,NLAY))

  !rate of each zone
  allocate(RBZONE(2,0:NBZONEMax,20))



!  ConstantCC = .false.
!  ConstantSC = .false.
  IF(IUNIT(1).GT.0) then
    rewind(IUNIT(1))
    call GWF2BCF7AR(IUNIT(1),0,IGRID)
  endif

  IF(IUNIT(23).GT.0) then
    rewind(IUNIT(23))
    CALL GWF2LPF7AR(IUNIT(23),IGRID)
  endif

  IF(IUNIT(37).GT.0) then
    rewind(IUNIT(37))
    CALL GWF2HUF7AR(IUNIT(37),IUNIT(47),IUNIT(53),0,IGRID)
  endif


  IF(IUNIT(62).GT.0) then
    rewind(IUNIT(62))
    !CALL GWF2UPW1AR(IUNIT(62), IGRID)
  endif


  ! determine the head-dependent flow terms used in the model
  call SDA_ReadALLCB()

  ! enter scenario run
  do isen = 1, abs(NSEN)
    call SDA_NewScen(iin)
  end do

  WRITE(*,*) ' Normal termination of scenario analysis'
endsubroutine


subroutine SDA_NewScen(iin)
  !advance to next scenario run
  use sda
  use GLOBAL
  use GWFBASMODULE, only            :   TOTIM,CDDNFM,IHEDFM,LBDDSV,PERTIM,VBNM,VBVL
  use PCGN
  implicit none
  integer                           ::  iin                   !main SDA file
  integer                           ::  KSTP, KPER
  integer                           ::  N,LLOC,IFLEN
  integer                           ::  ITYP1,ITYP2
  integer                           ::  INAM1,INAM2
  real                              ::  R
  real                              ::  BUDPERC               !mass balance error
  character*80                      ::  fileName              !mass balance error
  character*300                     ::  strline               !mass balance error
  INTEGER                           ::  IBDT(8,10)
  logical                           ::  sol


  DTRead = 0.d0
  DTStpEnd = 0.d0
  DTSolver = 0.d0
  CALL DATE_AND_TIME(VALUES=IBDT(:,1))
  !scenario name
  read(iin, *) ScenName
  ScenName = trim(ScenName)


  read(iin, "(A)") fileName
  !close the output files for the previous scenario run
  close(IOUT)
  !reopen IOUT
  open(unit=IOUT, file = fileName)

  strline = "(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ',I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2,/)"
  write(*,"(1x,A,A)") "Start running Scenario Analysis ", ScenName
  write(*, strline) (IBDT(LLOC,1),LLOC=1,3),(IBDT(LLOC,1),LLOC=5,7)


  write(IOUT,"(A,I4,A,a)") "START RUNNING SCENARIO ANALYSIS", isen, " PROJECT NAME:", ScenName
  strline = "(1X,'RUN START DATE AND TIME (YYYY/MM/DD HH:MM:SS): ',I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2,/)"
  write(IOUT, strline) (IBDT(LLOC,1),LLOC=1,3),(IBDT(LLOC,1),LLOC=5,7)

  !CSV file
  read(iin, "(A)") fileName
  close(iCBB)
  open(newunit=iCBB, file = fileName, action = "write", status = 'replace')
  !write header
  write(iCBB,"(A)") "TIME,ZONE,QIN,QOUT,TYPE"
  write(IOUT,"(A,A)") " HEAD AND FLOW WILL BE SAVED TO ",fileName
  !endif

  !solver
  LLOC=1
  read(iin, "(A)") strline
  CALL URWORD(strline,LLOC,ITYP1,ITYP2,1,N,R,IOUT,iin)

  CALL URWORD(strline,LLOC,INAM1,INAM2,1,N,R,IOUT,iin)


  !find a positive file number for solver
  sol = .true.
  do LLOC= 200, 500
    inquire(LLOC, opened = sol)
    if (.not. sol) exit
  enddo
  N = 0
  if (strline(ITYP1:ITYP2) == 'SIP') then
    close(IUNIT(9))
    IUNIT(9) = LLOC
    open(unit = IUNIT(9), file = strline(INAM1:INAM2),ACTION='READ',STATUS='OLD')
    call SIP7AR(IUNIT(9),MXITER,IGRID)
    N = N + 1
  else
    IUNIT(9) = 0
  endif

  if (strline(ITYP1:ITYP2) == 'DE4') then
    close(IUNIT(10))
    IUNIT(10) = LLOC
    open(unit = IUNIT(10), file = strline(INAM1:INAM2),ACTION='READ',STATUS='OLD')
    call DE47AR(IUNIT(10),MXITER,IGRID)
    N = N + 1
  else
    IUNIT(10) = 0
  endif

  if (strline(ITYP1:ITYP2) == 'PCG') then
    close(IUNIT(13))
    IUNIT(13) = LLOC
    open(unit = IUNIT(13), file = strline(INAM1:INAM2),ACTION='READ',STATUS='OLD')
    CALL PCG7AR(IUNIT(13),MXITER,IGRID)
    N = N + 1
  else
    IUNIT(13) = 0
  endif

  if (strline(ITYP1:ITYP2) == 'GMG') then
    close(IUNIT(42))
    IUNIT(42) = LLOC
    open(unit = IUNIT(42), file = strline(INAM1:INAM2),ACTION='READ',STATUS='OLD')
    CALL GMG7AR(IUNIT(42),MXITER,IGRID)
    N = N + 1
  else
    IUNIT(42) = 0
  endif


  if (strline(ITYP1:ITYP2) == 'PCGN') then
    close(IUNIT(59))
    IUNIT(59) = LLOC
    open(unit = IUNIT(59), file = strline(INAM1:INAM2),ACTION='READ',STATUS='OLD')
    CALL PCGN2AR(IUNIT(59),IFREFM,MXITER,IGRID)
    N = N + 1
  else
    IUNIT(59) = 0
  endif

!  IF(IUNIT(63).GT.0) then
!    IF(IUNIT(21).GT.0) then
!      rewind(IUNIT(21))
!      CALL GWF2HFB7AR(IUNIT(21),IGRID)
!    endif
!    rewind(IUNIT(63))
!    CALL GWF2NWT1AR(IUNIT(63),MXITER,IUNIT(22),0,IGRID)
!  endif
!




  if (N/=1) call USTOP('Something wrong with Solver')





  !zone file
  write(IOUT, "(/,A)") " READING ZONE FILE ... "
  read(iin, "(A)") fileName
  open(newunit = N, file = fileName, ACTION='READ', status = 'OLD')
  call SDA_Zone(N)
  if (NBZONE == 0) call USTOP('ZONE NUMBER IS ZERO')
  write(IOUT, "(A,I3,A)") " THERE ARE ", NBZONE ," ZONES."


  !accretional input flow
  read(iin, *) QHED, QWEL, QRCH
  !write(IOUT, "(/,A,1PE10.2,A)") "  FLOW SAVE THRESHOLD: ", QTor, "; FLOWS LOWER THAN THIS VALUE WILL NOT BE SAVED TO OUTPUT FILE."
  if (QHED /= 0) then
    read(iin, "(A)") fileName
    if (QHED>0) then
      LBDDSV = 1
      CDDNFM = '(10(1X1PE13.5))'
      open(newunit = iHED, file = fileName)
    else
      OPEN(newunit=iHED,FILE=fileName,FORM=FORM,ACCESS=ACCESS)
    end if

  end if

  !close(IUNIT(2))
  !close(IUNIT(8))

  if (QWEL > 0) then
    if (IUNIT(2) == 0) then
      call USTOP('THE WEL PACKAGE MUST BE WITH A FILE UNIT NUMBER IN THE NAME FILE')
    endif
    read(iin, "(A)") fileName
    close(IUNIT(2))
    open(unit = IUNIT(2), file = fileName, action = 'read')
    call GWF2WEL7AR(IUNIT(2),IGRID)
  endif
  if (QRCH > 0) then
    if (IUNIT(8) == 0) then
      call USTOP('THE RCH PACKAGE MUST BE WITH A FILE UNIT NUMBER IN THE NAME FILE')
    endif
    read(iin, "(A)") fileName
    close(IUNIT(8))
    open(unit = IUNIT(8), file = fileName, action = 'read')
    call GWF2RCH7AR(IUNIT(8),IGRID)
  endif



  !start the loop


  CBVBVL = 0.0d0

  TOTIM = 0.0
  PERTIM = 0.0
  HNEW = 0.0d0
  CALL DATE_AND_TIME(VALUES=IBDT(:,2))
  iSTEP = 0

  CR=0.;CC=0.;CV=0.;HCOF=0.;SDASC=0.


  !Go back to the file at the beginning
  do KPER = 1, NPER
    rewind(iACTIVE(KPER))
    rewind(iCC(KPER))
    rewind(iCR(KPER))
    rewind(iCV(KPER))
    rewind(iHCOF(KPER))
    rewind(iSC(KPER))
    rewind(iCB(KPER))
  enddo


  do KPER = 1, NPER

    CALL GWF2BAS7ST(KPER,IGRID)


    IF(QWEL > 0) CALL GWF2WEL7RP(IUNIT(2),IGRID)
    IF(QRCH > 0) CALL GWF2RCH7RP(IUNIT(8),IGRID)

    CALL UMESPR('SOLVING FOR HEAD',' ',IOUT)
    WRITE(*,"(' Scenario ',A,' Solving:  Stress period: ',i5,4x,'Eqn.')") trim(ScenName),KPER

    do KSTP = 1, NSTP(KPER)
      call SDA_NewTSP(KSTP, KPER)
    enddo
  enddo

  CALL DATE_AND_TIME(VALUES=IBDT(:,3))

  IF(QWEL >0)         CALL GWF2WEL7DA(IGRID)
  IF(QRCH >0)         CALL GWF2RCH7DA(IGRID)
  IF(IUNIT(9).GT.0)   CALL SIP7DA(IGRID)
  IF(IUNIT(10).GT.0)  CALL DE47DA(IGRID)
  IF(IUNIT(13).GT.0)  CALL PCG7DA(IGRID)
  IF(IUNIT(42).GT.0)  CALL GMG7DA(IGRID)
  IF(IUNIT(59).GT.0)  CALL PCGN2DA(IGRID)

  CALL GLO1BAS6ET(IOUT,IBDT(:,1),1)

  CALL DATE_AND_TIME(VALUES=IBDT(:,4))
  write (IOUT, "(A)") "CPU Time Summary"
  call TIMEDIFF(IBDT(:,1),IBDT(:,2),R)
  write (IOUT, "(A,F10.3,A)") " INITIALIZING SCENARIO USES ", R," SECONDS IN TOTAL"
  write (IOUT, "(A,F10.3,A)") " READING FILE USES          ", DTRead," SECONDS IN TOTAL"
  write (IOUT, "(A,F10.3,A)") " SOLVING EQ. USES           ", DTSolver," SECONDS IN TOTAL"
  write (IOUT, "(A,F10.3,A)") " WRITING RESULTS USES       ", DTStpEnd," SECONDS IN TOTAL"
  call TIMEDIFF(IBDT(:,1),IBDT(:,4),R)
  write (IOUT, "(A,F10.3,A)") " TOTAL CPU TIME IS          ", R," SECONDS IN TOTAL"

endsubroutine


subroutine SDA_Zone(iin)
  use sda
  use GLOBAL
  implicit none
  integer                           ::  iin
  integer                           ::  NL,NR,NC



  READ(iin,*) NL,NR,NC

  IF(NC.NE.NCOL .OR. NR.NE.NROW .OR. NL.NE.NLAY) THEN
    WRITE(*,*) 'MISMATCH BETWEEN DIMENSIONS OF CELL-BY-CELL DATA AND ZONE DATA:'
    WRITE(*,*) 'LAYERS, ROWS, COLUMNS IN ZONE DATA:',NL,NR,NC
    WRITE(*,*) 'LAYERS, ROWS, COLUMNS IN CELL-BY-CELL FILE:', NLAY,NROW,NCOL
    STOP
  END IF

  CALL IZREAD(BZONE,NLAY,NROW,NCOL,NBZONE,iin,99,IOUT)

  NBZONE = maxval(BZONE)

  close(iin)

end subroutine


subroutine SDA_NewTSP(KSTP, KPER)
  !advance to next time step
  use sda
  USE PCGMODULE
  USE SIPMODULE
  USE DE4MODULE
  USE GMGMODULE
  USE PCGN
  use GLOBAL
  use GWFBASMODULE

  implicit none
  integer                           ::  KPER, KSTP, KITER, KKITER, KKSTP, KKPER
  integer                           ::  i,j,k,kk,IER,ICNVG
  real                              ::  BUDPERC               !mass balance error
  integer                           ::  ITT(8,10)
  real                              ::  dt

  call date_and_time(values = ITT(:,1))

  !iSTEP = iSTEP + 1
  KKSTP=KSTP; KKPER=KPER
  !del time
  CALL GWF2BAS7AD(KKPER,KKSTP,IGRID)


  call SDA_ReadCB(KKSTP,KKPER)

  !read IBOUND, CV,CC,CR,HCOF
  call SDA_ReadHCOF(KSTP, KPER)

  SDASC = SDASC0 / DELT




  call date_and_time(values = ITT(:,2))
  !if (debug) go to 101
  !read addtional well or recharge
  RHS = - HOLD * SDASC

  if (QWEL > 0) call GWF2WEL7FM(IGRID)
  if (QRCH > 0) call GWF2RCH7FM(IGRID)

  !if (debug) write(idebug, "(10(1PE15.7))") RHS


  !initial guess
  !HNEW = 0.0d0

  !7C2----ITERATIVELY FORMULATE AND SOLVE THE FLOW EQUATIONS.
  DO KITER = 1, MXITER
    KKITER = KITER
    !CALL GWF2BAS7FM(IGRID)
    !7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
    IER=0
    IF (IUNIT(9) /= 0) THEN
      CALL SIP7PNT(IGRID)
      CALL SIP7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,EL,FL,GL, &
        V,W,HDCG,LRCH,NPARM,KKITER,HCLOSE,ACCL,ICNVG, &
        KKSTP,KKPER,IPCALC,IPRSIP,MXITER,NSTP(KKPER), &
        NCOL,NROW,NLAY,NODES,IOUT,0,IER)
    END IF
    IF (IUNIT(10) /= 0) THEN
      CALL DE47PNT(IGRID)
      CALL DE47AP(HNEW,IBOUND,AU,AL,IUPPNT,IEQPNT,D4B,MXUP, &
        MXLOW,MXEQ,MXBW,CR,CC,CV,HCOF,RHS,ACCLDE4,KITER, &
        ITMX,MXITER,NITERDE4,HCLOSEDE4,IPRD4,ICNVG,NCOL, &
        NROW,NLAY,IOUT,LRCHDE4,HDCGDE4,IFREQ,KKSTP,KKPER, &
        DELT,NSTP(KKPER),ID4DIR,ID4DIM,MUTD4, &
        DELTL,NBWL,NUPL,NLOWL,NLOW,NEQ,NUP,NBW,IER)
    END IF
    IF (IUNIT(13) /= 0) THEN
      CALL PCG7PNT(IGRID)
      CALL PCG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,VPCG,SS, &
        P,CD,HCHG,LHCH,RCHG,LRCHPCG,KKITER,NITER, &
        HCLOSEPCG,RCLOSEPCG,ICNVG,KKSTP,KKPER,IPRPCG, &
        MXITER,ITER1,NPCOND,NBPOL,NSTP(KKPER),NCOL,NROW, &
        NLAY,NODES,RELAXPCG,IOUT,MUTPCG,IT1,DAMPPCG,BUFF, &
        HCSV,IER,HPCG,DAMPPCGT,ISSFLG(KKPER),HDRY, &
        IHCOFADD)
    END IF
    !            IF (IUNIT(14).GT.0) THEN
    !              CALL LMG7PNT(IGRID)
    !              CALL LMG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,A,IA,JA,U1,
    !     1           FRHS,IG,ISIZ1,ISIZ2,ISIZ3,ISIZ4,KKITER,BCLOSE,DAMPLMG,
    !     2           ICNVG,KKSTP,KKPER,MXITER,MXCYC,NCOL,NROW,NLAY,NODES,
    !     3           HNOFLO,IOUT,IOUTAMG,ICG,IADAMPLMG,DUPLMG,DLOWLMG)
    !            END IF
    IF (IUNIT(42) /= 0) THEN
      CALL GMG7PNT(IGRID)
      CALL GMG7AP(HNEW,RHS,CR,CC,CV,HCOF,HNOFLO,IBOUND, &
        IITER,MXITER,RCLOSEGMG,HCLOSEGMG, &
        KKITER,KKSTP,KKPER,NCOL,NROW,NLAY,ICNVG, &
        SITER,TSITER,DAMPGMG,IADAMPGMG,IOUTGMG, &
        IOUT,GMGID, &
        IUNITMHC,DUP,DLOW,CHGLIMIT, &
        BIGHEADCHG,HNEWLAST)
    ENDIF
    IF (IUNIT(59) /= 0) THEN
      CALL PCGN2AP(HNEW,RHS,CR,CC,CV,HCOF,IBOUND,KKITER,KKSTP,KKPER,ICNVG,HNOFLO,IGRID)
      !call USTOP('Not supporting PCGN')
    ENDIF
    IF(IER.EQ.1) CALL USTOP(' ')
    !
    !7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
    IF (ICNVG.EQ.1) exit
  enddo

  call date_and_time(values = ITT(:,3))
  IF(ICNVG.EQ.0) THEN
    WRITE(IOUT,* ) ' FAILURE TO MEET SOLVER CONVERGENCE CRITERIA'
    call USTOP(' FAILURE TO MEET SOLVER CONVERGENCE CRITERIA')
  END IF

  !outputs
  if (QHED > 0) then
    do k=1, NLAY
      kk=k
      call ULASV2(real(HNEW(:,:,K)),' HEAD DIFFERENCE',KSTP,KPER,PERTIM,TOTIM,NCOL,NROW,KK,iHED,CDDNFM,LBDDSV,IBOUND(:,:,K))
    enddo
  elseif (QHED < 0) then
    do k=1, NLAY
      kk=k
      call ULASAV(real(HNEW(:,:,K)),' HEAD DIFFERENCE',KSTP,KPER,PERTIM,TOTIM,NCOL,NROW,KK,iHED)
    enddo
  endif



  !output
  call SDA_BD(KKSTP,KKPER)

  call TIMEDIFF(ITT(:,1),ITT(:,2),dt)
  DTRead = DTRead + dble(dt)
  call date_and_time(values = ITT(:,4))
  call TIMEDIFF(ITT(:,2),ITT(:,3),dt)
  DTSolver = DTSolver + dble(dt)
  call TIMEDIFF(ITT(:,3),ITT(:,4),dt)
  DTStpEnd = DTStpEnd + dble(dt)

endsubroutine


subroutine SDA_BD(KSTP,KPER)
  !calculate volumetric budgets
  !use numbers
  use sda
  use GLOBAL,only                   :   HNEW,HOLD,NCOL,NROW,NLAY,IBOUND,BUFF,IUNIT,IOUT,ITMUNI
  use GWFBASMODULE,only             :   TOTIM,DELT,VBNM,VBVL,MSUM,PERTIM
  USE GWFWELMODULE, ONLY            :   NWELLS,WELL
  USE GWFRCHMODULE,ONLY             :   NRCHOP,IRCHCB,RECH,IRCH

  implicit none
  integer                           ::  KSTP,KPER
  integer                           ::  i, IR, IC, IL, NHC
  integer,allocatable               ::  IIR(:), IIC(:), IIL(:)
  double precision                  ::  Q               !flow term
  !double precision                  ::  QQSC(NCOL,NROW,NLAY)
  real                              ::  Qnet
  character*4                       ::  ty
  !print the accretion flow first
  CBVBVL(3,:) = 0.d0
  CBVBVL(4,:) = 0.d0
  RBZONE = 0.d0


  iTypeBC = 1
  !Storage
  !where (IBOUND > 0)
    !QQSC = SDASC*(HOLD-HNEW)
  !end where

  write (IOUT,"(/,A,1PE16.9/)") "TOTAL DRAWDOWN IS ", sum(HNEW)
  DO IR=1,NROW
    DO IC=1,NCOL
      DO IL=1,NLAY
        if (IBOUND(IC,IR,IL) > 0) then
          Q = SDASC(IC,IR,IL)*(HOLD(IC,IR,IL)-HNEW(IC,IR,IL))
          call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HNEW(IC,IR,IL),Q)
        endif
      enddo  !IL=1,NLAY
    enddo !IC=1,NCOL
  enddo !IR=1,NROW

  !constant head and storage
  iTypeBC = 2
  IF(IUNIT(1).GT.0) then
    CALL SDA_BCF7BDCH(KSTP,KPER)
  endif

  IF(IUNIT(23).GT.0) then
    CALL SDA_LPF7BDCH(KSTP,KPER)
  endif

  iTypeBC = 3
  if (QWEL > 0) then
    do i=1,NWELLS
      IR=WELL(2,i)
      IC=WELL(3,i)
      IL=WELL(1,i)
      IF(IBOUND(IC,IR,IL) > 0) then
        Q=WELL(4,i)
        call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HNEW(IC,IR,IL),Q)
      endif
    enddo
  endif

  iTypeBC = 4
  if (QRCH > 0) then
    IF(NRCHOP.EQ.1) THEN
      !
      !5------NRCHOP=1, SO RECH GOES INTO LAYER 1. PROCESS EACH HORIZONTAL
      !5------CELL LOCATION.
      IL = 1
      DO IR=1,NROW
        DO IC=1,NCOL
          !
          !5A-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.

          IF(IBOUND(IC,IR,1).GT.0) THEN
            Q=RECH(IC,IR)
            call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HNEW(IC,IR,IL),Q)

          END IF
        enddo !IR=1,NROW
      enddo !IC=1,NCOL
    ELSE IF(NRCHOP.EQ.2) THEN
      !
      !6------NRCHOP=2, RECH IS IN LAYER SPECIFIED IN INDICATOR ARRAY(IRCH).
      !6------PROCESS EACH HORIZONTAL CELL LOCATION.
      DO  IR=1,NROW
        DO  IC=1,NCOL
          !
          !6A-----GET LAYER INDEX FROM INDICATOR ARRAY(IRCH).
          IL=IRCH(IC,IR)
          !
          !6B-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.
          IF(IL.EQ.0) exit
          IF(IBOUND(IC,IR,IL).GT.0) THEN
            Q=RECH(IC,IR)
            call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HNEW(IC,IR,IL),Q)

          END IF
        enddo !IC=1,NCOL
      enddo !IR=1,NROW
    ELSE
      !
      !7------NRCHOP=3; RECHARGE IS INTO HIGHEST CELL IN A VERTICAL COLUMN
      !7------THAT IS NOT NO FLOW.  PROCESS EACH HORIZONTAL CELL LOCATION.
      DO IR=1,NROW
        DO IC=1,NCOL
          !
          !7A-----INITIALIZE IRCH TO 1, AND LOOP THROUGH CELLS IN A VERTICAL
          !7A-----COLUMN TO FIND WHERE TO PLACE RECHARGE.
          IRCH(IC,IR)=1
          DO IL=1,NLAY
            !
            !7C-----IF CELL IS VARIABLE HEAD, THEN DO BUDGET FOR IT.
            IF (IBOUND(IC,IR,IL).GT.0) THEN
              Q=RECH(IC,IR)
              call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HNEW(IC,IR,IL),Q)
              exit
            END IF
          enddo  !IL=1,NLAY
        enddo !IC=1,NCOL
      enddo !IR=1,NROW
      !
    END IF
  endif



  do i = 1, NCB
    IC = pSDACB(i)%ICOL
    IR = pSDACB(i)%IROW
    IL = pSDACB(i)%ILAY
    if (IBOUND(IC,IR,IL) > 0) then
      Q = -pSDACB(i)%CB*HNEW(IC,IR,IL)
      iTypeBC = pSDACB(i)%BND
      call SDA_SaveBD(0,iSTEP,IC,IR,IL,HNEW(IC,IR,IL),Q)
    endif
  enddo

  !print result
  VBVL(1,1:MSUM) = real(CBVBVL(1,1:MSUM))
  VBVL(2,1:MSUM) = real(CBVBVL(2,1:MSUM))
  VBVL(3,1:MSUM) = real(CBVBVL(3,1:MSUM))
  VBVL(4,1:MSUM) = real(CBVBVL(4,1:MSUM))

  CALL SGWF2BAS7V(MSUM+1,VBNM,VBVL,KSTP,KPER,IOUT,Qnet)
  call SGWF2BAS7T(KSTP,KPER,DELT,PERTIM,TOTIM,ITMUNI,IOUT)
  write(IOUT, "(//)")


  !save flow by zone
  do i =1, MSUM
    !if (.not. PrintType(i)) cycle
    do IL = 1, NBZONE
      if (abs(RBZONE(1,IL,i)) <= 1.d-99) RBZONE(1,IL,i) = 0.d0
      if (abs(RBZONE(2,IL,i)) <= 1.d-99) RBZONE(2,IL,i) = 0.d0
      WRITE(iCBB,"(F0.3,',',I0,',',1PE16.9,',',1PE16.9,',',A)") TOTIM,IL,RBZONE(1,IL,i),RBZONE(2,IL,i),trim(VBNM(i))
    enddo
  enddo


endsubroutine


SUBROUTINE SDA_SaveBD(IIO,ISTP,ICOL,IROW,ILAY,DH,DF)
  !calculate and save the budgets
  use GWFBASMODULE, only            :   TOTIM, DELT, VBVL,MSUM
  use sda, only                     :   RBZONE,BZONE,CBVBVL,iTypeBC,CBINDX
  integer                           ::  IIO                   !file number
  integer                           ::  ICOL                  !col
  integer                           ::  IROW                  !row
  integer                           ::  ILAY                  !layer
  integer                           ::  ISTP                  !current time step
  integer                           ::  IDX

  integer                           ::  NCOL                  !col
  integer                           ::  NROW                  !row
  integer                           ::  NLAY                  !layer

  doubleprecision                   ::  DH                    !head change
  doubleprecision                   ::  DF                    !flow change


  integer                           ::  INODE                 !layer

  !not outputing small flows
  IDX = CBINDX(iTypeBC)
  !Qnet = Qnet + DF
  INODE = BZONE(ICOL,IROW,ILAY)
  if (DF>0) then
    CBVBVL(1,IDX) = CBVBVL(1,IDX) + DELT*DF
    CBVBVL(3,IDX) = CBVBVL(3,IDX) + DF
    RBZONE(1,INODE,IDX) = RBZONE(1,INODE,IDX) + DF
  else
    CBVBVL(2,IDX) = CBVBVL(2,IDX) - DELT*DF
    CBVBVL(4,IDX) = CBVBVL(4,IDX) - DF
    RBZONE(2,INODE,IDX) = RBZONE(2,INODE,IDX) - DF
  endif



  !INODE = (ILAY-1)*(NROW*NCOL)+(IROW-1)*NCOL + ICOL
  !WRITE(IIO,"(F12.3,',',I10,',',1PE15.8,',',1PE15.8,',',A4)") TOTIM,INODE,DH,DF,BType(K)
END subroutine


subroutine SDA_CVT_NSEN
  use SDA, only: NSEN
  NSEN = abs(NSEN)
end subroutine

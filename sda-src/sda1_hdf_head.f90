

subroutine SDA_Ini(iin)
  !initialize
  use sda
  use HDF
  use GLOBAL,only                         :   NSTP,NPER,NCOL,NROW,NLAY,IUNIT,IOUT,IBOUND
  use GWFBASMODULE,only                   :   MSUM,VBNM
  integer                                 ::  iin                   !input file number of SDA main file
  integer                                 ::  iset, nset            !input file number of SDA main file
  integer                                 ::  iper, iper0           !input file number of SDA main file
  character*400                           ::  sline                 !input file number of SDA main file
  character*100                           ::  fhcof
  character(len=256)                      ::  h5fname
  integer                                 ::  i1,i2,i3, itmp, iloc  !input file number of SDA main file
  real                                    ::  rtmp                  !input file number of SDA main file
  integer(hid_t)                          ::  strtype
  integer                                 ::  ierr

  ! initialize HDF5 (extremely important)
  call h5open_f(ierr)
  if (ierr /= 0) call USTOP('error when initializing HDF5')

  ! define hdf datatypes to be used
  call h5tcopy_f(H5T_NATIVE_INTEGER, h5t_int, ierr)

  IGRID = 1

  nset = 0
  iper0 = 1
  NBZONEMax = 99

  iloc = 1
  CALL URDCOM(iin,IOUT,sline)
  call URWORD(sline,iloc,i1,i2,2,NSEN,rtmp,IOUT,iin)
  call URWORD(sline,iloc,i1,i2,2,NCBMAX,rtmp,IOUT,iin)
  if (NCBMAX<0) then
    NCBMAX=-NCBMAX
    CALL h5fcreate_f('debug.h5', H5F_ACC_TRUNC_F, ihdf5_db, iERR)
    if (iERR /= 0) call USTOP('Program cannot create HDF5 file: debug.h5')
  end if
  allocate(cbCOL(NCBMAX), cbROW(NCBMAX), cbLAY(NCBMAX), cbTyp(NCBMAX), cbCond(NCBMAX), cbHead(NCBMAX))
  allocate(cbCOL_o(NCBMAX), cbROW_o(NCBMAX), cbLAY_o(NCBMAX), cbTyp_o(NCBMAX), cbCond_o(NCBMAX), cbHead_o(NCBMAX))


  CALL URDCOM(iin,IOUT,sline)

  if (NSEN<0) then
    ! create a new HDF5 file
    CALL h5fcreate_f(trim(sline), H5F_ACC_TRUNC_F, ihdf5_in, iERR)
    if (iERR /= 0) call USTOP('Program cannot create HDF5 file: '//trim(sline))


  elseif (NSEN>0) then
    CALL h5fopen_f(trim(sline), H5F_ACC_RDWR_F, ihdf5_in, iERR)
    if (iERR /= 0) call USTOP('Program cannot open HDF5 file: '//trim(sline))
  end if
  !write(IOUT,*) " READING UNIT NUMBERS FOR SDA ... "


  if (debug) open(newunit=idebug, file="debug.txt")


  NSTEP = sum(NSTP)
  c_dims(1) = NLAY; c_dims(2) = NSTEP
  g_dims(1) = NCOL; g_dims(2) = NROW
  allocate(newCC(NLAY, NSTEP))
  allocate(newCV(NLAY, NSTEP))
  allocate(newCR(NLAY, NSTEP))
  allocate(newSC(NLAY, NSTEP))
  allocate(newIB(NLAY, NSTEP))
  allocate(newCB(NSTEP))


  ! save the active cells
  nActive = NLAY * NCOL * NROW - count(IBOUND == 0)
  allocate(actLay(nActive))
  allocate(actRow(nActive))
  allocate(actCol(nActive))
  nActive = 0
  do i1=1, NLAY
    do i2 = 1, NROW
      do i3 = 1, NCOL
        if (IBOUND(i3,i2,i1) /= 0) then
          nActive = nActive + 1
          actLay(nActive) = i1
          actRow(nActive) = i2
          actCol(nActive) = i3
        end if
      end do
    end do
  end do


  allocate(SC(NCOL,NROW,NLAY))

  if (NSEN < 0) then
    iSTEP = 0
    allocate(SDAIBD(NCOL,NROW,NLAY))
    allocate(SDACV(NCOL,NROW,NLAY))
    allocate(SDACC(NCOL,NROW,NLAY))
    allocate(SDACR(NCOL,NROW,NLAY))
    allocate(SDASC(NCOL,NROW,NLAY))

    SDAIBD = 0
    SDASC = 0.
    SDACV = 0.
    SDACC = 0.
    SDACR = 0.
    
    ! save initial head
    call SDA_SaveHCOF(0, 0)

  else
    CALL SDA_Run(iin)
    CALL USTOP(' ')
  end if



end subroutine

subroutine SDA_AddCB(IC,IR,IL,TyB,CB,HB)

  use sda
  integer :: IC,IR,IL,TyB
  real :: CB,HB

  if (NSEN == 0) return
  if (CB == 0.) return
  NCB=NCB+1

  cbCOL(NCB) = IC
  cbROW(NCB) = IR
  cbLAY(NCB) = IL
  cbCond(NCB)= CB
  cbTyp(NCB)= TyB

end subroutine

subroutine SDA_AddSC(IC,IR,IL,S)

  !SC is the Specific yield or Ss * Thickness
  use sda,only  : SC, NSEN
  integer :: IC,IR,IL
  real :: S

  if (NSEN>=0) return

!  if (IC/=0 .and. IR/=0) then
    SC(IC,IR,IL) = S
!  elseif (IC = 0 .and. IR = 0) then
!    SDASC(:,:,IL) = SC(:,:,IL)
!  elseif (IC = 0) then
!    SDASC(:,IR,IL) = SC(:,IR,IL)
!  endif

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
  use sda
  use HDF
  use GLOBAL, only: NLAY
  integer(hsize_t) :: m_dims(2)

  !call SDA_CalcErr
  iSTEP = iSTEP + 1
  call SDA_SaveCB(KSTP,KPER)
  call SDA_SaveHCOF(KSTP,KPER)

  ! save the flags at the end of the simulation
  m_dims(1) = NLAY; m_dims(2) = NSTEP
  if (iSTEP == NSTEP) then
    call saveHDFint(ihdf5_in, newCB, 'newCB', (/int(NSTEP, hsize_t)/), h5t_int)
    call saveHDFint(ihdf5_in, newCC, 'newCC', m_dims, h5t_int)
    call saveHDFint(ihdf5_in, newCR, 'newCR', m_dims, h5t_int)
    call saveHDFint(ihdf5_in, newCV, 'newCV', m_dims, h5t_int)
    call saveHDFint(ihdf5_in, newIB, 'newIB', m_dims, h5t_int)
    call saveHDFint(ihdf5_in, newSC, 'newSC', m_dims, h5t_int)
    !call saveHDFint(ihdf5_in, newHC, 'newHC', (/NLAY, NSTEP/), h5t_int)
  end if
end subroutine

SUBROUTINE SDA_SaveCB(KSTP,KPER)
  !save head-dependent boundary for the current time step
  use sda
  use HDF
  use GLOBAL, only : IOUT
  integer                           ::  i
  integer                           ::  KSTP                  !col
  integer                           ::  KPER                  !row
  character*4                       ::  txtStep
  integer                           ::  ierr


  write(txtStep, '(I0.4)') iSTEP

  if (all(cbCOL(1:NCB)==cbCOL_o(1:NCB)) .and. &
      all(cbROW(1:NCB)==cbROW_o(1:NCB)) .and. &
      all(cbLAY(1:NCB)==cbLAY_o(1:NCB)) .and. &
      all(cbTyp(1:NCB)==cbTyp_o(1:NCB)) .and. &
      all(cbCond(1:NCB)==cbCond_o(1:NCB))) then
    write(IOUT,'(//x,A,I0,A)') 'SDA -- REPEAT ', NCB, ' HEAD DEPENDENT BOUNDARY CONDITIONS '//txtStep
    newCB(iSTEP) = 0
  else
    write(IOUT,'(//x,A,I0,A)') 'SDA -- SAVING ', NCB, ' HEAD DEPENDENT BOUNDARY CONDITIONS '//txtStep
    newCB(iSTEP) = NCB
    cbCOL_o(1:NCB) = cbCOL(1:NCB)
    cbROW_o(1:NCB) = cbROW(1:NCB)
    cbLAY_o(1:NCB) = cbLAY(1:NCB)
    cbTyp_o(1:NCB) = cbTyp(1:NCB)
    cbCond_o(1:NCB)= cbCond(1:NCB)
    call saveHDFint(ihdf5_in, cbCOL, 'cbCol'//txtStep, (/int(NCB, hsize_t)/), h5t_int)
    call saveHDFint(ihdf5_in, cbROW, 'cbRow'//txtStep, (/int(NCB, hsize_t)/), h5t_int)
    call saveHDFint(ihdf5_in, cbLAY, 'cbLay'//txtStep, (/int(NCB, hsize_t)/), h5t_int)
    call saveHDFint(ihdf5_in, cbTyp, 'cbTyp'//txtStep, (/int(NCB, hsize_t)/), h5t_int)
    call saveHDFreal(ihdf5_in, cbCond, 'cbCond'//txtStep, (/int(NCB, hsize_t)/))
  end if

END

SUBROUTINE SDA_ReadCB(KSTP, KPER)
  !read head-dependent boundary for the current time step
  use sda
  use HDF
  use HDF
  integer                           ::  KSTP,KPER                  !step
  character*4                       ::  txtStep

  if (newCB(iSTEP) > 0) then
    write(txtStep, '(I0.4)') iSTEP
    NCB = newCB(iSTEP)
    call readHDFint(ihdf5_in,  cbCOL(1:NCB), 'cbCol'//txtStep, (/int(NCB,hsize_t)/))
    call readHDFint(ihdf5_in,  cbROW(1:NCB), 'cbRow'//txtStep, (/int(NCB,hsize_t)/))
    call readHDFint(ihdf5_in,  cbLAY(1:NCB), 'cbLay'//txtStep, (/int(NCB,hsize_t)/))
    call readHDFint(ihdf5_in,  cbTyp(1:NCB), 'cbTyp'//txtStep, (/int(NCB,hsize_t)/))
    call readHDFreal(ihdf5_in, cbCond(1:NCB), 'cbCond'//txtStep, (/int(NCB,hsize_t)/))
  end if
END


SUBROUTINE SDA_FindCBType()
  !check how many cb types
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
  use HDF
  use GLOBAL,only  :   IOUT,IBOUND,RHS,NCOL,NROW,NLAY,HNEW,HOLD

  CHARACTER*16     :: TEXT
  CHARACTER*7      :: txtStepLayer
  integer          :: K,KKPER,KKSTP,KSTP,KPER

  !       1234567890123456
  !KKPER = min(1,KPER)
  DO K=1,NLAY
    write(txtStepLayer, '(I0.4, I0.3)') iSTEP, K
    call saveHDFdbl(ihdf5_in, HNEW(:,:,K), 'HN'//txtStepLayer, g_dims)
    if (iSTEP == 0) return ! to save initial head only
    if (iSTEP == 1 .or. any(IBOUND(:,:,K) /= SDAIBD(:,:,K))) then
      newIB(K, iSTEP) = 1
      SDAIBD(:,:,K) = IBOUND(:,:,K)
      call saveHDFint(ihdf5_in, SDAIBD(:,:,K), 'IB'//txtStepLayer, g_dims, h5t_int)
    else
      newIB(K, iSTEP) = 0
    endif

  enddo

endsubroutine

subroutine SDA_ReadHCOF(KSTP, KPER)
  !read flow coefficients
  use sda
  use HDF
  use GLOBAL,only                         :   HCOF,IOUT,IBOUND,NCOL,NROW,NLAY,HNEW,HOLD,RHS,IUNIT,CC,CR,CV
  use GWFBASMODULE,only                   :   DELT
  integer ::  KSTP,KPER

  CHARACTER*7      :: txtStepLayer
  integer :: K, i

  DO K=1,NLAY
    write(txtStepLayer, '(I0.4, I0.3)') iSTEP, K
    call readHDFdbl(ihdf5_in,  HNEW(:,:,K), 'HN'//txtStepLayer, g_dims)        
    if (iSTEP==0) return ! read initial head
    if (newIB(K, iSTEP) == 1) call readHDFint(ihdf5_in, IBOUND(:,:,K), 'IB'//txtStepLayer, g_dims)    
  end do
  
  
  ! calculate A
  HCOF = 0.
  IF(IUNIT(1).GT.0)  CALL GWF2BCF7FM(1, KSTP, KPER, IGRID) ! update CC CR CV HCOF
  IF(IUNIT(23).GT.0) CALL GWF2LPF7FM(1, KSTP, KPER, IGRID) ! update CC CR CV HCOF
  SC = -HCOF * DELT
  RHS = HCHG0 * HCOF
  do i = 1, NCB
    HCOF(cbCOL(i), cbROW(i), cbLAY(i)) = HCOF(cbCOL(i), cbROW(i), cbLAY(i)) - cbCond(i)
  end do
endsubroutine



subroutine SDA_Run(iin)
  !start SDA
  !baseline run can be skipped. coefficients are read from the the file directly
  use sda
  use HDF
  use GLOBAL
  use GWFBASMODULE

  implicit none

  integer                                 ::  iin               !input file number of SDA main file
  integer                                 ::  k ,i1,i2,iloc     !layer index
  character*400                           ::  line              !layer index


  WRITE(*,*) ' Entering scenario analysis mode'


  if (NSEN<0) then

    ! release the array from baseline run
    NSEN = abs(NSEN)
  endif


  !free format
  IFREFM = 1

  ! read the flow coefficients
  ! call SDA_ReadHCOF()


  !zones for zone budget
  allocate(BZONE(NCOL,NROW,NLAY))

  !rate of each zone
  allocate(RBZONE(2,0:NBZONEMax,20))

  ! head change
  allocate(HCHG0(NCOL,NROW,NLAY),HCHG1(NCOL,NROW,NLAY))

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
  call SDA_FindCBType()
  
  write(iout, *) ' read data flags'
  call readHDFint(ihdf5_in, newIB, 'newIB', c_dims)
  call readHDFint(ihdf5_in, newCB, 'newCB', (/int(NSTEP, hsize_t)/))
  
  ! enter scenario run
  do isen = 1, abs(NSEN)
    call SDA_NewScen(iin)
  end do
  call h5fclose_f(ihdf5_in, iin)
  WRITE(*,*) ' Normal termination of scenario analysis'
endsubroutine


subroutine SDA_NewScen(iin)
  !advance to next scenario run
  use sda
  use HDF
  use GLOBAL
  use GWFBASMODULE, only            :   TOTIM,CDDNFM,IHEDFM,LBDDSV,PERTIM,VBNM,VBVL
  use PCGN
  implicit none
  integer                           ::  iin                   !main SDA file
  integer                           ::  KSTP, KPER
  integer                           ::  N,LLOC,IFLEN,K
  integer                           ::  ITYP1,ITYP2
  integer                           ::  INAM1,INAM2
  real                              ::  R
  real                              ::  BUDPERC               !mass balance error
  character*80                      ::  fileName              !mass balance error
  character*300                     ::  strline               !mass balance error
  INTEGER                           ::  IBDT(8,10)
  logical                           ::  sol
  integer                           ::  ierr
  character*7                       ::  txtStepLayer

  DTRead = 0.d0
  DTStpEnd = 0.d0
  DTSolver = 0.d0
  CALL DATE_AND_TIME(VALUES=IBDT(:,1))
  !scenario name
  read(iin, *) ScenName
  ScenName = trim(ScenName)

  ! list file name
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
      call USTOP('THE WEL PACKAGE MUST BE ALSO ACTIVATED IN THE BASELINE RUN')
    endif
    read(iin, "(A)") fileName
    close(IUNIT(2))
    open(unit = IUNIT(2), file = fileName, action = 'read')
    call GWF2WEL7AR(IUNIT(2),IGRID)
  endif
  if (QRCH > 0) then
    if (IUNIT(8) == 0) then
      call USTOP('THE RCH PACKAGE MUST BE ALSO ACTIVATED IN THE BASELINE RUN')
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
  HCHG1 = 0.0d0
  HCHG0 = 0.
  CALL DATE_AND_TIME(VALUES=IBDT(:,2))
  
  ! check out initial head
  iSTEP = 0
  call SDA_ReadHCOF(0, 0)
  
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
  use HDF
  USE PCGMODULE
  USE SIPMODULE
  USE DE4MODULE
  USE GMGMODULE
  USE PCGN, only: PCGN2AP
  use GLOBAL
  use GWFBASMODULE

  implicit none
  integer                           ::  KPER, KSTP, KITER, KKITER, KKSTP, KKPER
  integer                           ::  i,j,k,kk,IER,ICNVG
  real                              ::  BUDPERC               !mass balance error
  integer                           ::  ITT(8,10)
  real                              ::  dt
  CHARACTER*7      :: txtStepLayer

  call date_and_time(values = ITT(:,1))

  iSTEP = iSTEP + 1
  KKSTP=KSTP; KKPER=KPER
  !del time
  CALL GWF2BAS7AD(KKPER,KKSTP,IGRID)
  !save head change
  HCHG0 = HCHG1


  call SDA_ReadCB(KSTP, KPER)
  !read IBOUND, CV,CC,CR,HCOF
  call SDA_ReadHCOF(KSTP, KPER)
  call date_and_time(values = ITT(:,2))

  if (QWEL > 0) call GWF2WEL7FM(IGRID)
  if (QRCH > 0) call GWF2RCH7FM(IGRID)

  do K=1,NLAY
    write(txtStepLayer, '(I0.4, I0.3)') iSTEP, K
    call saveHDFreal(ihdf5_db, CC(:,:,K), 'CC'//txtStepLayer, g_dims)
    call saveHDFreal(ihdf5_db, CR(:,:,K), 'CR'//txtStepLayer, g_dims)
    call saveHDFreal(ihdf5_db, RHS(:,:,K), 'RHS'//txtStepLayer, g_dims)
    call saveHDFreal(ihdf5_db, HCOF(:,:,K), 'HCOF'//txtStepLayer, g_dims)
    call saveHDFreal(ihdf5_db, real(HCHG1(:,:,K)), 'DH'//txtStepLayer, g_dims)
  end do


  !7C2----ITERATIVELY FORMULATE AND SOLVE THE FLOW EQUATIONS.
  DO KITER = 1, MXITER
    KKITER = KITER
    !CALL GWF2BAS7FM(IGRID)
    !7C2B---MAKE ONE CUT AT AN APPROXIMATE SOLUTION.
    IER=0
    IF (IUNIT(9) /= 0) THEN
      CALL SIP7PNT(IGRID)
      CALL SIP7AP(HCHG1,IBOUND,CR,CC,CV,HCOF,RHS,EL,FL,GL, &
        V,W,HDCG,LRCH,NPARM,KKITER,HCLOSE,ACCL,ICNVG, &
        KKSTP,KKPER,IPCALC,IPRSIP,MXITER,NSTP(KKPER), &
        NCOL,NROW,NLAY,NODES,IOUT,0,IER)
    END IF
    IF (IUNIT(10) /= 0) THEN
      CALL DE47PNT(IGRID)
      CALL DE47AP(HCHG1,IBOUND,AU,AL,IUPPNT,IEQPNT,D4B,MXUP, &
        MXLOW,MXEQ,MXBW,CR,CC,CV,HCOF,RHS,ACCLDE4,KITER, &
        ITMX,MXITER,NITERDE4,HCLOSEDE4,IPRD4,ICNVG,NCOL, &
        NROW,NLAY,IOUT,LRCHDE4,HDCGDE4,IFREQ,KKSTP,KKPER, &
        DELT,NSTP(KKPER),ID4DIR,ID4DIM,MUTD4, &
        DELTL,NBWL,NUPL,NLOWL,NLOW,NEQ,NUP,NBW,IER)
    END IF
    IF (IUNIT(13) /= 0) THEN
      CALL PCG7PNT(IGRID)
      CALL PCG7AP(HCHG1,IBOUND,CR,CC,CV,HCOF,RHS,VPCG,SS, &
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
      CALL GMG7AP(HCHG1,RHS,CR,CC,CV,HCOF,HNOFLO,IBOUND, &
        IITER,MXITER,RCLOSEGMG,HCLOSEGMG, &
        KKITER,KKSTP,KKPER,NCOL,NROW,NLAY,ICNVG, &
        SITER,TSITER,DAMPGMG,IADAMPGMG,IOUTGMG, &
        IOUT,GMGID, &
        IUNITMHC,DUP,DLOW,CHGLIMIT, &
        BIGHEADCHG,HNEWLAST)
    ENDIF
    IF (IUNIT(59) /= 0) THEN
      CALL PCGN2AP(HCHG1,RHS,CR,CC,CV,HCOF,IBOUND,KKITER,KKSTP,KKPER,ICNVG,HNOFLO,IGRID)
      !call USTOP('Not supporting PCGN')
    ENDIF
    IF(IER.EQ.1) CALL USTOP(' ')
    !
    !7C2C---IF CONVERGENCE CRITERION HAS BEEN MET STOP ITERATING.
    IF (ICNVG.EQ.1) exit
  enddo

  call date_and_time(values = ITT(:,3))
  
  

  !outputs
  if (QHED > 0) then
    do k=1, NLAY
      kk=k
      call ULASV2(real(HCHG1(:,:,K)),' HEAD DIFFERENCE',KSTP,KPER,PERTIM,TOTIM,NCOL,NROW,KK,iHED,CDDNFM,LBDDSV,IBOUND(:,:,K))
    enddo
  elseif (QHED < 0) then
    do k=1, NLAY
      kk=k
      call ULASAV(real(HCHG1(:,:,K)),' HEAD DIFFERENCE',KSTP,KPER,PERTIM,TOTIM,NCOL,NROW,KK,iHED)
    enddo
  endif



  !output
  call SDA_BD(KKSTP,KKPER,BUDPERC)
  IF(ICNVG.EQ.0) THEN
    WRITE(IOUT,* ) ' FAILURE TO MEET SOLVER CONVERGENCE CRITERIA'
    IF(ABS(BUDPERC)>STOPER) THEN
      write(IOUT, *) 'PERCENT DISCREPANCY: ', BUDPERC
      call USTOP(' FAILURE TO MEET SOLVER CONVERGENCE CRITERIA')
    end if
  END IF
  
  call TIMEDIFF(ITT(:,1),ITT(:,2),dt)
  DTRead = DTRead + dble(dt)
  call date_and_time(values = ITT(:,4))
  call TIMEDIFF(ITT(:,2),ITT(:,3),dt)
  DTSolver = DTSolver + dble(dt)
  call TIMEDIFF(ITT(:,3),ITT(:,4),dt)
  DTStpEnd = DTStpEnd + dble(dt)

endsubroutine


subroutine SDA_BD(KSTP,KPER,BUDPERC)
  !calculate volumetric budgets
  !use numbers
  use sda
  use GLOBAL,only                   :   HNEW,HOLD,NCOL,NROW,NLAY,IBOUND,BUFF,IUNIT,IOUT,ITMUNI
  use GWFBASMODULE,only             :   TOTIM,DELT,VBNM,VBVL,MSUM,PERTIM
  USE GWFWELMODULE, ONLY            :   NWELLS,WELL
  USE GWFRCHMODULE,ONLY             :   NRCHOP,IRCHCB,RECH,IRCH

  implicit none
  integer                           ::  KSTP,KPER
  real                              ::  BUDPERC ! budget discrepancy
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
  

  write (IOUT,"(/,A,1PE16.9/)") "TOTAL DRAWDOWN IS ", sum(HCHG1)
  DO IR=1,NROW
    DO IC=1,NCOL
      DO IL=1,NLAY
        if (IBOUND(IC,IR,IL) > 0) then
          Q = SC(IC,IR,IL)*(HCHG0(IC,IR,IL)-HCHG1(IC,IR,IL)) / DELT
          call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HCHG1(IC,IR,IL),Q)
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
        call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HCHG1(IC,IR,IL),Q)
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
            call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HCHG1(IC,IR,IL),Q)

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
            call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HCHG1(IC,IR,IL),Q)

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
              call SDA_SaveBD(iCBB,iSTEP,IC,IR,IL,HCHG1(IC,IR,IL),Q)
              exit
            END IF
          enddo  !IL=1,NLAY
        enddo !IC=1,NCOL
      enddo !IR=1,NROW
      !
    END IF
  endif



  do i = 1, NCB
    IC = cbCOL(i)
    IR = cbROW(i)
    IL = cbLAY(i)
    if (IBOUND(IC,IR,IL) > 0) then
      Q = -cbCond(i)*HCHG1(IC,IR,IL)
      iTypeBC = cbTyp(i)
      ! debug. 
      if (abs(Q) > 1000.) then
        write(IOUT, '(/,A5,3I5)') 'CRL= ', IC, IR, IL
        write(IOUT, *) 'IB = ', IBOUND(IC,IR,IL)
        write(IOUT, *) 'dH = ', HCHG1(IC,IR,IL)
        write(IOUT, *) 'dH0= ', HCHG0(IC,IR,IL)
        write(IOUT, *) 'CB = ', cbCond(i)
        write(IOUT, *) 'Tp = ', CBTYPE(iTypeBC)
        write(IOUT, *) 'Q  = ', Q
      end if
      call SDA_SaveBD(0,iSTEP,IC,IR,IL,HCHG1(IC,IR,IL),Q)
    endif
  enddo

  !print result
  VBVL(1,1:MSUM) = real(CBVBVL(1,1:MSUM))
  VBVL(2,1:MSUM) = real(CBVBVL(2,1:MSUM))
  VBVL(3,1:MSUM) = real(CBVBVL(3,1:MSUM))
  VBVL(4,1:MSUM) = real(CBVBVL(4,1:MSUM))

  CALL SGWF2BAS7V(MSUM+1,VBNM,VBVL,KSTP,KPER,IOUT,BUDPERC)
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

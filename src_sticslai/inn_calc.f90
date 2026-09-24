MODULE inn_calculation

  USE constantes
  USE constantes_soil
  USE pft_parameters

  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
! inn_calc  — STICS-like nitrogen status (ORCSTICS version)
!
! INPUTS:
!   n_leaf : tN/ha
!   c_leaf : tDM/ha
!   ivm    : optional PFT index. If not provided, defaults to PFT12.
!
! OUTPUTS:
!   innlai    : N stress index on leaf expansion  [0..1]
!   innsenes  : N stress index on leaf senescence [0..1]
!   inns      : N stress index used for RUE/growth [0..1]
!   inn       : Nitrogen Nutrition Index
!   nc        : Critical N concentration (% DM)
!   cnplante  : Actual N concentration (% DM)
!------------------------------------------------------------------------------

SUBROUTINE inn_calc(n_leaf, c_leaf, innlai, innsenes, inns, inn, nc, cnplante, ivm)

  ! Inputs
  REAL, INTENT(IN)  :: n_leaf   ! tN/ha
  REAL, INTENT(IN)  :: c_leaf   ! tDM/ha

  ! Optional PFT index
  INTEGER, INTENT(IN), OPTIONAL :: ivm

  ! Outputs
  REAL, INTENT(OUT) :: innlai, innsenes, inns
  REAL, INTENT(OUT) :: inn, nc, cnplante

  ! Local variables
  INTEGER :: P_codeplisoleN
  INTEGER :: P_codeINN
  INTEGER :: ipft

  REAL :: P_adilmax, P_bdilmax, P_masecNmax, P_wmax, P_masecNmin
  REAL :: adilI, bdilI, adilmaxI, bdilmaxI

  REAL :: P_INNmin, P_INNimin, P_innturgmin, P_innsen, P_QNpltminINN

  REAL :: VabsN, deltabso
  REAL :: tDM, QNplante, masecpartiel, masecdil, masecabso, W_eff
  REAL :: NCmax, dNdWc, dQNc, inni
  REAL :: magrain, absodrp

  REAL :: cnleaf_min, cnleaf_max

  ! ---------------------------------------------------------------------------
  ! Minimal PFT adaptation
  ! Default = PFT12, for backward compatibility.
  ! If caller passes ivm=13, use maize parameters.
  ! ---------------------------------------------------------------------------
  ipft = 12
  IF (PRESENT(ivm)) ipft = ivm

  IF (ipft < 1 .OR. ipft > nvm) THEN
     ipft = 12
  ENDIF

  cnleaf_min = cn_leaf_min_mtc(ipft)
  cnleaf_max = cn_leaf_max_mtc(ipft)

  ! ===== Local defaults / parameters =====
  P_codeplisoleN = 1
  P_codeINN      = 1

  ! Critical dilution parameters
  P_adilmax   = 8.5
  P_bdilmax   = 0.44
  P_masecNmax = 0.20
  P_masecNmin = 4.0

  ! Alternative set for isolated-plant option, not used for now
  adilI    = 4.40
  bdilI    = 0.37
  adilmaxI = 5.00
  bdilmaxI = 0.25

  ! Stress mapping and guards
  P_INNmin       = 0.36
  P_INNimin      = 0.30
  P_innturgmin   = 0.20
  P_innsen       = 0.17
  P_QNpltminINN  = 1.0

  ! Instantaneous vars, unused when P_codeINN=1
  VabsN    = 0.0
  deltabso = 0.0

  ! ===== Computation =====

  ! 1) Working state variables
  tDM          = MAX(c_leaf, 1.0e-12)       ! tDM/ha
  QNplante     = n_leaf * 1000.0            ! kgN/ha
  masecpartiel = tDM
  masecdil     = tDM

  ! Non-grain proxy at this level
  magrain   = 0.0
  absodrp   = 1.0
  masecabso = masecdil - (magrain/100.0) + (absodrp*magrain/100.0)

  ! 2) Actual plant %N
  ! QNplante is kgN/ha, tDM is tDM/ha.
  ! Divide by 10 to convert to % DM.
  cnplante = QNplante / (tDM * 10.0)

  ! 3) Critical N concentration curve with biomass floor
  W_eff = MAX(tDM, P_masecNmax)

  IF (P_codeplisoleN == 2) THEN

     nc    = adilI    * W_eff**(-bdilI)
     NCmax = adilmaxI * W_eff**(-bdilmaxI)
     dNdWc = 10.0 * adilI * (1.0 - bdilI) * W_eff**(-bdilI)

  ELSE

     nc    = adil(ipft) * W_eff**(-bdil(ipft))
     NCmax = P_adilmax     * W_eff**(-bdil(ipft))

  ENDIF

  ! 4) Nitrogen Nutrition Index
  IF (nc > 0.0) THEN
     inn = cnplante / nc
     inn = MIN(inn, NCmax / nc)
  ELSE
     inn = 1.0
  ENDIF

  ! 5) Cumulative INN path
  inns = MIN(1.0, inn)

  ! Floor/clamp as in STICS
  inns = MAX(P_INNmin, inns)
  inns = MIN(1.0, MAX(P_INNmin, inn))

  ! 6) Map inns -> innlai
  innlai = (P_innturgmin - 1.0) / (P_INNmin - 1.0) * inns + &
           (1.0 - (P_innturgmin - 1.0) / (P_INNmin - 1.0))

  innlai = MIN(1.0, MAX(P_INNmin, innlai))

  ! 7) Map inns -> innsenes
  innsenes = (P_innsen - 1.0) / (P_INNmin - 1.0) * inns + &
             (1.0 - (P_innsen - 1.0) / (P_INNmin - 1.0))

  innsenes = MIN(1.0, MAX(P_INNmin, innsenes))

  ! 8) Startup guard: if plant N is tiny, disable N stress
  IF (QNplante < P_QNpltminINN) THEN
     inn      = 1.0
     inns     = 1.0
     innlai   = 1.0
     innsenes = 1.0
  ENDIF

END SUBROUTINE inn_calc

END MODULE inn_calculation

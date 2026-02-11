MODULE weathergen
  ! Wrapper utilities to generate hourly meteorology from daily inputs
  ! using existing routines in getmet.f90. Includes a deterministic PRNG
  ! to ensure cross-platform reproducibility for rainfall timing.

  USE maestcom
  USE, INTRINSIC :: iso_c_binding
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: wg_generate_day_internal
  PUBLIC :: wg_seed_internal
  PUBLIC :: wg_seed
  PUBLIC :: wg_generate_day

  INTEGER(KIND=8) :: rng_state = 88172645463393265_8

CONTAINS

  SUBROUTINE wg_seed_internal(seed)
    ! Set deterministic RNG seed
    INTEGER(KIND=8), INTENT(IN) :: seed
    rng_state = MERGE(seed, 88172645463393265_8, seed /= 0_8)
  END SUBROUTINE wg_seed_internal

  REAL FUNCTION rng_uniform()
    ! xorshift64* generator -> [0,1)
    INTEGER(KIND=8) :: x
    x = rng_state
    x = IEOR(x, SHIFTL(x,13))
    x = IEOR(x, SHIFTR(x,7))
    x = IEOR(x, SHIFTL(x,17))
    rng_state = x
    ! Map to (0,1): use unsigned-like scaling by dividing by 2^63 and taking abs
    rng_uniform = REAL(ABS(INT(x,8)), KIND=KIND(1.0)) / REAL(HUGE(1_8), KIND=KIND(1.0))
  END FUNCTION rng_uniform

  SUBROUTINE wg_assignrain_det(total_mm, ppt)
    ! Deterministic version of ASSIGNRAIN using the local PRNG.
    REAL, INTENT(IN) :: total_mm
    REAL, INTENT(OUT) :: ppt(MAXHRS)  ! mm per hour
    INTEGER :: ihr, irain, ihrswithrain, i
    REAL :: rain, rate, r

    DO ihr = 1, KHRS
      ppt(ihr) = 0.0
    END DO

    IF (total_mm .LE. 2.0) THEN
      r = rng_uniform()
      irain = INT(r * KHRS) + 1
      IF (irain < 1) irain = 1
      IF (irain > KHRS) irain = KHRS
      ppt(irain) = total_mm
    ELSE IF (total_mm .GT. 46.0) THEN
      rain = total_mm / REAL(KHRS)
      DO ihr = 1, KHRS
        ppt(ihr) = rain
      END DO
    ELSE
      ihrswithrain = MIN( INT( (total_mm/2.0) * KHRS / 24.0 ), KHRS )
      ihrswithrain = MAX(ihrswithrain, 1)
      rate = total_mm / REAL(ihrswithrain)
      DO i = 1, ihrswithrain
        r = rng_uniform()
        irain = INT(r * KHRS) + 1
        IF (irain < 1) irain = 1
        IF (irain > KHRS) irain = KHRS
        ppt(irain) = ppt(irain) + rate
      END DO
    END IF
  END SUBROUTINE wg_assignrain_det

  SUBROUTINE wg_generate_day_internal(idate, alat, &
                              tmin, tmax, sw_mean_wm2, precip_mm, &
                              press_pa, &
                              tair, tsoil, rh, vpd, vmfd, &
                              radabv, fbeam, ppt, press)
    ! Generate hourly series for a single day from daily inputs.
    ! Inputs
    INTEGER, INTENT(IN) :: idate          ! days-since-1950
    REAL,    INTENT(IN) :: alat           ! radians
    REAL,    INTENT(IN) :: tmin, tmax     ! deg C
    REAL,    INTENT(IN) :: sw_mean_wm2    ! W m-2 daily mean shortwave
    REAL,    INTENT(IN) :: precip_mm      ! mm/day
    REAL,    INTENT(IN) :: press_pa       ! Pa (daily mean)
    ! Outputs (length KHRS)
    REAL,    INTENT(OUT) :: tair(MAXHRS), tsoil(MAXHRS)
    REAL,    INTENT(OUT) :: rh(MAXHRS), vpd(MAXHRS), vmfd(MAXHRS)
    REAL,    INTENT(OUT) :: radabv(MAXHRS,3), fbeam(MAXHRS,3)
    REAL,    INTENT(OUT) :: ppt(MAXHRS)
    REAL,    INTENT(OUT) :: press(MAXHRS)

    ! Externals from getmet.f90 and radn.f90
    EXTERNAL :: CALCTHRLY, CALCRH, RHTOVPD, VPDTOMFD
    EXTERNAL :: BRISTO, CALCFBMD, CALCPARHRLY, CALCNIR, CALCFSUN, THERMAL, CALCTSOIL
    EXTERNAL :: ZENAZ, SUN
    INTEGER, EXTERNAL :: JDATE

    INTEGER :: ihr, idoy
    REAL :: fbm, radbm, raddf
    REAL :: zen(MAXHRS), az(MAXHRS), fsun(MAXHRS)
    REAL :: par_wm2_to_mj
    REAL :: dec, eqntim_loc, dayl, sunset

    ! Compute solar geometry for this day at latitude alat
    idoy = JDATE(idate)
    CALL SUN(idoy, alat, 0.0, dec, eqntim_loc, dayl, sunset)
    CALL ZENAZ(alat, 0.0, 0.0, dec, eqntim_loc, zen, az)

    ! Air pressure constant across hours
    DO ihr = 1, KHRS
      press(ihr) = press_pa
    END DO

    ! Hourly air temperatures
    CALL CALCTHRLY(tmax, tmin, dayl, tair)

    ! Soil temperature as mean daily air temperature
    CALL CALCTSOIL(tair, tsoil)

    ! Humidity: RH from Tmin-as-dewpoint, then VPD, then VMFD
    CALL CALCRH(tmin, tair, rh)
    CALL RHTOVPD(rh, tair, vpd)
    CALL VPDTOMFD(vpd, press, vmfd)

    ! Radiation: convert daily mean SW (W m-2) to daily total MJ m-2 d-1
    par_wm2_to_mj = (sw_mean_wm2 * 86400.0) / 1.0e6
    ! Estimate beam fraction and distribute hourly
    CALL CALCFBMD(idate, zen, par_wm2_to_mj*FPAR, fbm)
    radbm = par_wm2_to_mj*FPAR * fbm
    raddf = par_wm2_to_mj*FPAR * (1.0 - fbm)
    CALL CALCPARHRLY(radbm, raddf, zen, radabv, fbeam)
    CALL CALCNIR(radabv, fbeam)

    ! Sunlit fraction and thermal radiation
    CALL CALCFSUN(fbeam, fsun)
    CALL THERMAL(tair, vpd, fsun, radabv)

    ! Rainfall allocation
    CALL wg_assignrain_det(precip_mm, ppt)

  END SUBROUTINE wg_generate_day_internal

  SUBROUTINE wg_init() BIND(C, NAME='wg_init')
    ! C API: Initialize the module. Must be called before any other functions.
    ! Sets global parameters and RNG state.
    KHRS = 24
    SPERHR = 3600 * 24.0 / KHRS
    HHRS = (KHRS) / 2.0
  END SUBROUTINE wg_init

  ! C API: Seed the deterministic RNG used for rainfall timing and any stochastic components.
  ! Parameters
  ! - seed [int64]: deterministic seed. Use 0 to keep the existing default seed.
  ! Returns
  ! - 0 on success.
  INTEGER(c_int) FUNCTION wg_seed(seed) BIND(C, NAME='wg_seed')
    INTEGER(c_int64_t), VALUE :: seed
    CALL wg_seed_internal(seed)
    wg_seed = 0_c_int
  END FUNCTION wg_seed

  ! C API: Generate hourly meteorology for a single day from daily inputs.
  ! Purpose
  !   Wraps wg_generate_day for FFI. All scalars are passed by value, arrays as
  !   flat contiguous buffers. Output arrays must be preallocated by the caller.
  ! Parameters (inputs)
  ! - idate [int]: day index (days-since-1950; used by SUN/JDATE).
  ! - alat [float]: latitude in radians.
  ! - tmin, tmax [float]: daily Tmin/Tmax in deg C.
  ! - sw_mean_wm2 [float]: daily mean shortwave radiation [W m^-2].
  ! - precip_mm [float]: daily precipitation total [mm/day].
  ! - press_pa [float]: daily mean air pressure [Pa].
  ! - nhrs [int]: length of the hourly series (must equal KHRS compiled into the lib).
  ! Parameters (outputs)
  ! - tair[nhrs], tsoil[nhrs] [float]: air and soil temperature [deg C].
  ! - rh[nhrs] [float]: relative humidity [%].
  ! - vpd[nhrs] [float]: vapour pressure deficit.
  ! - vmfd[nhrs] [float]: vapour mole fraction deficit.
  ! - radabv[nhrs*3] [float]: above-canopy radiation components, column-major layout
  !   with 3 columns; element (i,j) is at index (j-1)*nhrs + i. Columns are
  !   PAR, NIR, thermal radiation.
  ! - fbeam[nhrs*3] [float]: beam fraction components, same layout as radabv.
  ! - ppt[nhrs] [float]: hourly precipitation [mm/hr], deterministic allocation.
  ! - press[nhrs] [float]: pressure [Pa].
  ! Returns
  ! - 0 on success; 1 if nhrs != KHRS.
  INTEGER(c_int) FUNCTION wg_generate_day(idate, alat, &
                               tmin, tmax, sw_mean_wm2, precip_mm, &
                               press_pa, &
                               nhrs, &
                               tair, tsoil, rh, vpd, vmfd, &
                               radabv, fbeam, ppt, press) &
                               BIND(C, NAME='wg_generate_day')
    INTEGER(c_int), VALUE :: idate
    REAL(c_float), VALUE :: alat
    REAL(c_float), VALUE :: tmin, tmax, sw_mean_wm2, precip_mm
    REAL(c_float), VALUE :: press_pa
    INTEGER(c_int), VALUE :: nhrs
    REAL(c_float), INTENT(OUT) :: tair(*), tsoil(*), rh(*), vpd(*), vmfd(*)
    REAL(c_float), INTENT(OUT) :: radabv(*), fbeam(*), ppt(*), press(*)
    REAL(c_float) :: radabv2(MAXHRS,3), fbeam2(MAXHRS,3)
    CALL wg_generate_day_internal(idate, alat, &
                         tmin, tmax, sw_mean_wm2, precip_mm, &
                         press_pa, &
                         tair, tsoil, rh, vpd, vmfd, &
                         radabv2, fbeam2, &
                         ppt, press)
    ! Copy 2D outputs back into flat column-major buffers length nhrs*3
    BLOCK
      INTEGER :: i, j, idx
      DO j = 1, 3
        DO i = 1, nhrs
          idx = (j-1)*nhrs + i
          radabv(idx) = radabv2(i,j)
          fbeam(idx)  = fbeam2(i,j)
        END DO
      END DO
    END BLOCK
    wg_generate_day = 0_c_int
  END FUNCTION wg_generate_day

END MODULE weathergen

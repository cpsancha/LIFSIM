/*
 * mainV03_56_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "mainV03_56".
 *
 * Model version              : 1.477
 * Simulink Coder version : 8.11 (R2016b) 25-Aug-2016
 * C source code generated on : Mon May 08 13:32:21 2017
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_mainV03_56_private_h_
#define RTW_HEADER_mainV03_56_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#include "mainV03_56.h"
#ifndef ATMOS_TYPEDEF

typedef enum { COESA = 1, MILHDBK310, MILSTD210C } AtmosTypeIdx;

typedef enum { PROFILE = 1, ENVELOPE } ModelIdx;

typedef enum { HIGHTEMP = 1, LOWTEMP, HIGHDENSITY,
  LOWDENSITY, HIGHPRESSURE, LOWPRESSURE } VarIdx;

typedef enum { PP1 = 1, PP10 } PPercentIdx;

typedef enum { K5 = 1, K10, K20, K30, K40 } PAltIdx;

typedef enum { EXTREME = 1, P1, P5, P10, P20 } EPercentIdx;

#define ATMOS_TYPEDEF
#endif                                 /* ATMOS_TYPEDEF */

#ifndef ATMOS_DEFINE
#define PRESSURE0                      101325.0                  /*  N/m^2                  */
#define TEMPERATURE0                   288.15                    /*  K                      */
#define GRAV_CONST                     9.80665                   /*  m/s^2                  */
#define MOL_WT                         28.9644                   /*  kg/kgmol (air)         */
#define R_HAT                          8314.32                   /*  J/kgmol.K (gas const.) */
#define GAMMA                          1.4                       /*  (specific heat ratio) */
#define GMR                            ( GRAV_CONST * MOL_WT / R_HAT )
#define ATMOS_DEFINE
#endif                                 /* ATMOS_DEFINE */

#ifndef COESA76_DEFINE_DATA

/* 1976 COESA atmosphere model */
#define NUM1976PTS                     8

static real_T altitude76[NUM1976PTS] = {/* in meters (m) */
  0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0 };

static real_T tempGradient76[NUM1976PTS] = {/* in K/m  */
  (-0.0065), 0.0, 0.0010, 0.0028, 0.0, -0.0028, -0.0020, -0.0020 };

#define COESA76_DEFINE_DATA
#endif                                 /* COESA76_DEFINE_DATA */

#ifndef GRAVITY2_TYPEDEF

typedef enum { WGS84TAYLORSERIES = 1, WGS84CLOSEAPPROX,
  WGS84EXACT } GravityTypeIdx;

typedef enum { METRIC = 1, ENGLISH } UnitIdx;

typedef enum { JANUARY = 1, FEBRUARY, MARCH, APRIL, MAY, JUNE, JULY,
  AUGUST, SEPTEMBER, OCTOBER, NOVEMBER, DECEMBER } MonthIdx;

#define GRAVITY2_TYPEDEF
#endif                                 /* GRAVITY2_TYPEDEF */

#ifndef WGS84_DEFINE
#define WGS84_A                        6378137.0                 /* Semi-major Axis (m) */
#define WGS84_INV_F                    298.257223563             /* Reciprocal of Flattening */
#define WGS84_W_DEF                    7292115.0e-11             /* WGS 84 Angular Velocity of Earth (rad/sec)*/
#define WGS84_GM_DEF                   3986004.418e+8            /* Earth's Gravitational Const. (m^3/sec^2) */
#define WGS84_GM_PRM                   3986000.9e+8              /* Earth's Grav. Const. (m^3/sec^2) [no atmos]*/
#define WGS84_W_PRM                    7292115.1467e-11          /* IAU Angular Velocity of Earth (rad/sec) */
#define WGS84_G_E                      9.7803253359              /* Theoretical (Normal) Gravity at the Equator
                                                                    (on the Ellipsoid) (m/s^2) */
#define WGS84_K                        0.00193185265241          /* Theoretical (Normal) Gravity Formula Const.*/
#define WGS84_E_2                      6.69437999014e-3          /* First Eccentricity Squared */
#define WGS84_EL                       5.2185400842339e+5        /* Linear Eccentricity */
#define WGS84_B                        6356752.3142              /* Semi-minor Axis (m) */
#define WGS84_B_A                      0.996647189335            /* Axis Ratio */
#define M2FT                           1.0/0.3048
#define AERO_PI                        3.14159265358979323846
#define DEG2RAD                        AERO_PI/180.0
#define YEAR2000                       2000
#define WGS84_DEFINE
#endif                                 /* WGS84_DEFINE */

/* Used by FromWorkspace Block: '<S463>/FromWs' */
#ifndef rtInterpolate
# define rtInterpolate(v1,v2,f1,f2)    (((v1)==(v2))?((double)(v1)): (((f1)*((double)(v1)))+((f2)*((double)(v2)))))
#endif

#ifndef rtRound
# define rtRound(v)                    ( ((v) >= 0) ? floor((v) + 0.5) : ceil((v) - 0.5) )
#endif

extern real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_roundd_snf(real_T u);
extern real_T rt_modd_snf(real_T u0, real_T u1);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern void rt_mrdivide_U1d1x3_U2d3x3_Yd1x3_snf(const real_T u0[3], const real_T
  u1[9], real_T y[3]);
void InitCalcAtmosCOESA(real_T *temperature76, real_T *pressureRatio76);
void CalcAtmosCOESA(const real_T *altitude, real_T *temp, real_T *pressure,
                    real_T *density, real_T *speedofsound, real_T *temperature76,
                    real_T *pressureRatio76, int_T numPoints);
real_T calc_Julian_date(real_T uDataMonth, real_T uDataDay, real_T uDataYear);
int_T wgs84_calc_shared_vars(real_T uDataMonth, real_T uDataDay, real_T
  uDataYear, int_T uDataPrecessing, int_T uDataCentrifugal, real_T h, real_T E2,
  real_T cosphi, real_T sinphi, real_T sin2phi, real_T coslambda, real_T
  sinlambda, real_T GM, real_T *gamma_u_ptr, real_T *gamma_beta_ptr, real_T
  *cosbeta_ptr, real_T *sinbeta_ptr, real_T *u_ptr, real_T *u2E2_ptr, real_T
  *w_ptr);
void wgs84_exact(real_T *h, real_T *phi, real_T *lambda, real_T uDataMonth,
                 real_T uDataDay, real_T uDataYear, int_T uDataPrecessing, int_T
                 uDataCentrifugal, real_T E2, real_T GM, real_T opt_m2ft, real_T
                 *y, real_T *gamma_h, real_T *gamma_phi,int_T k);
extern int32_T plook_s32dd_bincp(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, int32_T *prevIndex);
extern uint32_T plook_bincpa(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, uint32_T *prevIndex);
extern real_T intrp2d_la_pw(const uint32_T bpIndex[], const real_T frac[], const
  real_T table[], uint32_T stride, const uint32_T maxIndex[]);
extern int32_T plook_s32dd_bincpa(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, int32_T *prevIndex);
extern real_T intrp8d_s32dla_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[], const uint32_T maxIndex[]);
extern real_T intrp7d_s32dla_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[], const uint32_T maxIndex[]);
extern real_T intrp6d_s32dla_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[], const uint32_T maxIndex[]);
extern real_T intrp5d_s32dla_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[], const uint32_T maxIndex[]);
extern real_T intrp8d_s32dl_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[]);
extern real_T intrp7d_s32dl_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[]);
extern real_T intrp6d_s32dl_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[]);
extern real_T intrp5d_s32dl_pw(const int32_T bpIndex[], const real_T frac[],
  const real_T table[], const uint32_T stride[]);
extern uint32_T plook_evenca(real_T u, real_T bp0, real_T bpSpace, uint32_T
  maxIndex, real_T *fraction);
extern real_T intrp3d_la_pw(const uint32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride[], const uint32_T maxIndex[]);
extern uint32_T plook_binca(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction);
extern int32_T binsearch_s32d_prevIdx(real_T u, const real_T bp[], uint32_T
  startIndex, uint32_T maxIndex);
extern uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T
  startIndex, uint32_T maxIndex);
extern uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex);
extern void joyinput(SimStruct *rts);
extern void mainV03_56_Distanceintogusty_Init(B_Distanceintogusty_mainV03_56_T
  *localB, P_Distanceintogusty_mainV03_56_T *localP,
  X_Distanceintogusty_mainV03_56_T *localX);
extern void mainV03_56_Distanceintogusty_Reset(P_Distanceintogusty_mainV03_56_T *
  localP, X_Distanceintogusty_mainV03_56_T *localX);
extern void mainV03_56_Distanceintogusty_Deriv(real_T rtu_V,
  DW_Distanceintogusty_mainV03_56_T *localDW, P_Distanceintogusty_mainV03_56_T
  *localP, X_Distanceintogusty_mainV03_56_T *localX,
  XDot_Distanceintogusty_mainV03_56_T *localXdot, real_T rtp_d_m);
extern void mainV03_56_Distanceintogusty_Disable
  (DW_Distanceintogusty_mainV03_56_T *localDW);
extern void mainV03_56_Distanceintogusty(RT_MODEL_mainV03_56_T * const
  mainV03_56_M, boolean_T rtu_Enable, B_Distanceintogusty_mainV03_56_T *localB,
  DW_Distanceintogusty_mainV03_56_T *localDW, P_Distanceintogusty_mainV03_56_T
  *localP, X_Distanceintogusty_mainV03_56_T *localX, real_T rtp_d_m);

/* private model entry point functions */
extern void mainV03_56_derivatives(void);

#endif                                 /* RTW_HEADER_mainV03_56_private_h_ */

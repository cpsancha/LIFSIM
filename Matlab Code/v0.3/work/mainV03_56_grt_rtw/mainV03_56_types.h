/*
 * mainV03_56_types.h
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

#ifndef RTW_HEADER_mainV03_56_types_h_
#define RTW_HEADER_mainV03_56_types_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "zero_crossing_types.h"
#ifndef DEFINED_TYPEDEF_FOR_ActuatorsBus_
#define DEFINED_TYPEDEF_FOR_ActuatorsBus_

typedef struct {
  real_T deltae;
  real_T deltar;
  real_T deltafr;
  real_T deltafl;
} ActuatorsBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_LLABus_
#define DEFINED_TYPEDEF_FOR_LLABus_

typedef struct {
  real_T Latitude_deg;
  real_T Longitude_deg;
  real_T Altitude_m;
} LLABus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_EulerBus_
#define DEFINED_TYPEDEF_FOR_EulerBus_

typedef struct {
  real_T phi;
  real_T theta;
  real_T psi;
} EulerBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_OmegaBus_
#define DEFINED_TYPEDEF_FOR_OmegaBus_

typedef struct {
  real_T p;
  real_T q;
  real_T r;
} OmegaBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_plantDataBus_
#define DEFINED_TYPEDEF_FOR_plantDataBus_

typedef struct {
  LLABus LLA;
  real_T V_ned[3];
  real_T X_ned[3];
  EulerBus Euler;
  real_T DCM_body_earth[9];
  real_T V_body[3];
  OmegaBus Omega_body;
  real_T dOmega_body[3];
  real_T Accel_body[3];
  real_T alpha;
  real_T beta;
  real_T alpha_dot;
  real_T beta_dot;
  real_T CG[3];
  real_T remainingCapacity;
} plantDataBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_envDataBus_
#define DEFINED_TYPEDEF_FOR_envDataBus_

typedef struct {
  real_T h_ground;
  real_T Fground[3];
  real_T Mground[3];
  real_T windVelocity[3];
  real_T windOmega[3];
  real_T T;
  real_T a;
  real_T P;
  real_T AirDensity;
  real_T Fgravity[3];
  real_T Gravity_ned[3];
  real_T MagneticField_ned[3];
} envDataBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_LLAMeasBus_
#define DEFINED_TYPEDEF_FOR_LLAMeasBus_

typedef struct {
  real_T LatitudeMeas_deg;
  real_T LongitudeMeas_deg;
  real_T AltitudeMeas_m;
} LLAMeasBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_OmegaMeasBus_
#define DEFINED_TYPEDEF_FOR_OmegaMeasBus_

typedef struct {
  real_T pMeas;
  real_T qMeas;
  real_T rMeas;
} OmegaMeasBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_SensorsBus_
#define DEFINED_TYPEDEF_FOR_SensorsBus_

typedef struct {
  LLAMeasBus LLAMeas;
  real_T X_nedMeas[3];
  real_T VearthMeas[3];
  real_T AccelMeas_body[3];
  OmegaMeasBus OmegaMeas_body;
  real_T DCMMeas_body_earth[9];
  real_T EulerMeas[3];
  real_T VaeroMeas_body[3];
  real_T alphaMeas;
  real_T betaMeas;
  real_T CGMeas[3];
  real_T remainingCapacityMeas;
} SensorsBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_CommandBus_
#define DEFINED_TYPEDEF_FOR_CommandBus_

typedef struct {
  real_T roll_cmd;
  real_T pitch_cmd;
  real_T yaw_cmd;
  real_T altitude_cmd;
  real_T Throttle1_cmd;
  real_T flightCondition;
  real_T Position_cmd[3];
} CommandBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_ThrottleBus_
#define DEFINED_TYPEDEF_FOR_ThrottleBus_

typedef struct {
  real_T Throttle1;
  real_T Throttle2;
  real_T Throttle3;
  real_T Throttle4;
  real_T Throttle5;
} ThrottleBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_betaBus_
#define DEFINED_TYPEDEF_FOR_betaBus_

typedef struct {
  int32_T idxbeta;
  real_T fbeta;
} betaBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_alphaBus_
#define DEFINED_TYPEDEF_FOR_alphaBus_

typedef struct {
  int32_T idxalpha;
  real_T falpha;
} alphaBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_altBus_
#define DEFINED_TYPEDEF_FOR_altBus_

typedef struct {
  int32_T idxalt;
  real_T falt;
} altBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_deltaeBus_
#define DEFINED_TYPEDEF_FOR_deltaeBus_

typedef struct {
  int32_T idxdeltae;
  real_T fdeltae;
} deltaeBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_deltarBus_
#define DEFINED_TYPEDEF_FOR_deltarBus_

typedef struct {
  int32_T idxdeltar;
  real_T fdeltar;
} deltarBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_deltafrBus_
#define DEFINED_TYPEDEF_FOR_deltafrBus_

typedef struct {
  int32_T idxdeltafr;
  real_T fdeltafr;
} deltafrBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_deltaflBus_
#define DEFINED_TYPEDEF_FOR_deltaflBus_

typedef struct {
  int32_T idxdeltafl;
  real_T fdeltafl;
} deltaflBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_xcgBus_
#define DEFINED_TYPEDEF_FOR_xcgBus_

typedef struct {
  int32_T idxxcg;
  real_T fxcg;
} xcgBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_analysisCasesBus_
#define DEFINED_TYPEDEF_FOR_analysisCasesBus_

typedef struct {
  alphaBus alpha;
  betaBus beta;
  altBus alt;
  xcgBus xcg;
  deltaeBus deltae;
  deltarBus deltar;
  deltafrBus deltafr;
  deltaflBus deltafl;
} analysisCasesBus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_PPAnPhCDcFDL0ZPEazAckG_
#define DEFINED_TYPEDEF_FOR_struct_PPAnPhCDcFDL0ZPEazAckG_

typedef struct {
  real_T K;
  real_T D;
  real_T Mu;
} struct_PPAnPhCDcFDL0ZPEazAckG;

#endif

/* Parameters for system: '<S47>/Distance into gust (y)' */
typedef struct P_Distanceintogusty_mainV03_56_T_
  P_Distanceintogusty_mainV03_56_T;

/* Parameters for system: '<S465>/J-K Flip-Flop' */
typedef struct P_JKFlipFlop_mainV03_56_T_ P_JKFlipFlop_mainV03_56_T;

/* Parameters (auto storage) */
typedef struct P_mainV03_56_T_ P_mainV03_56_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_mainV03_56_T RT_MODEL_mainV03_56_T;

#endif                                 /* RTW_HEADER_mainV03_56_types_h_ */
